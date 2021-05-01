# -*- coding: utf-8 -*-
import subprocess
import json
import os
import io
from multiprocessing import Pool
import multiprocessing
import multiprocessing.pool
from operator import itemgetter
import random
import string
import pickle
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pysam
import mip_classes as mod
import pandas as pd
from pandas.errors import MergeError
import gzip
from primer3 import calcHeterodimerTm
import primer3
import traceback
from msa_to_vcf import msa_to_vcf as msa_to_vcf
import itertools
import sys
import allel
from Bio import SeqIO

print("functions reloading")
# backbone dictionary
mip_backbones = {
    "hybrid_bb": "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAGNNNNNNNNNNNN",
    "hybrid_split": "NNNNAGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAGNNNNNNNNNN",
    "hybrid_split_hp": "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAGNNNNNNNNNN",
    "gc_bb": "GCAGATCGGAAGAGCACACCTCGCCAAGCTTTCGGCNNNNNNNNNNNN",
    "slx_bb": "CTTCAGCTTCCCGATCCGACGGTAGTGTNNNNNNNNNNNN"
}

"""
# Below class allows processors from a pool from multiprocessing module to
create processor pools of their own.
# http://mindcache.io/2015/08/09/python-multiprocessing-module-daemonic-processes-are-not-allowed-to-have-children.html
class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)
# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
"""


# above code was broken when switching to python 3. Below is taken from:
# https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic/8963618#8963618
class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NoDaemonProcessPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(NoDaemonProcessPool, self).__init__(*args, **kwargs)


# Exception wrapper for multiprocessing taken from
# https://stackoverflow.com/questions/6126007/python-getting-a-traceback-from-a-multiprocessing-process/26096355#26096355
class ExceptionWrapper(object):

    def __init__(self, ee, exc):
        self.ee = ee
        self.exc = exc
        __,  __, self.tb = sys.exc_info()

    def re_raise(self):
        print(self.exc)
        raise self.ee.with_traceback(self.tb)


###############################################################
# Region prep related functions
###############################################################


def coordinate_to_target(coordinates, snp_locations, capture_size):
    """ Create MIP targets starting from a snp file that is produced offline,
    usually from Annovar. This is a tab separated file with the following
    chr1	2595307	2595307	A	G	rs3748816.
    This can be generalized to any target with coordinates.
    """
    # create target regions to cover all snps
    # start by getting snps on same chromosome together
    snp_chroms = {}
    reference_snp_locations = rsl = coordinates
    for r in rsl:
        chrom = rsl[r]["chrom"]
        try:
            snp_chroms[chrom].append([rsl[r]["begin"],
                                      rsl[r]["end"]])
        except KeyError:
            snp_chroms[chrom] = [[rsl[r]["begin"],
                                  rsl[r]["end"]]]
    # merge snps that are too close to get separate regions
    # the length should be twice the capture size
    merged_snp_chroms = {}
    for c in snp_chroms:
        merged_snp_chroms[c] = merge_overlap(snp_chroms[c], 2 * capture_size)
    # create regions for alignment
    for c in merged_snp_chroms:
        regions = merged_snp_chroms[c]
        for r in regions:
            snps_in_region = []
            for s in reference_snp_locations:
                if ((reference_snp_locations[s]["chrom"] == c)
                    and (r[0] <= reference_snp_locations[s]["begin"]
                         <= reference_snp_locations[s]["end"] <= r[1])):
                    snps_in_region.append(s)
            r.append(snps_in_region)
        for reg in regions:
            snps = reg[2]
            reg_begin = reg[0]
            reg_end = reg[1]
            reg_locations = []
            for s in snps:
                s_locations = []
                locations = snp_locations[s]
                ref_location = reference_snp_locations[s]
                ref_begin = ref_location["begin"]
                ref_end = ref_location["end"]
                left_flank_buffer = ref_begin - reg_begin + capture_size
                right_flank_buffer = reg_end - ref_end + capture_size
                for l in locations:
                    snp_chrom = l["chrom"]
                    snp_begin = l["begin"]
                    snp_end = l["end"]
                    tar_begin = snp_begin - left_flank_buffer
                    tar_end = snp_end + right_flank_buffer
                    s_locations.append([snp_chrom, tar_begin, tar_end])
                reg_locations.append(s_locations)
            reg.append(reg_locations)
    # create target coordinate for each region
    target_coordinates = {}
    for c in merged_snp_chroms:
        regions = merged_snp_chroms[c]
        for reg in regions:
            region_name = "-".join(reg[2])
            region_targets = reg[3][0]
            for i in range(len(region_targets)):
                reg_name = region_name + "-" + str(i)
                if reg_name in target_coordinates:
                    print((reg_name, " is already in targets!"))
                else:
                    target_coordinates[reg_name] = region_targets[i]
    return target_coordinates


def rsid_to_target(resource_dir, snp_file):
    """ Create MIP targets starting from a snp file that is produced offline,
    usually from Annovar. This is a tab separated file with the following
    content: chr1	2595307	2595307	A	G	rs3748816.
    This can be generalized to any target with coordinates.
    """
    # one snp can have multiple locations on the reference genome,
    # this can happen with snps in regions where there are multiple different
    # assemblies (HLA locus, for example). So first step is to get each of
    # these locations in the genome.
    snp_locations = {}
    capture_types = {}
    with io.open(os.path.join(resource_dir, snp_file),
                 encoding="utf-8") as infile:
        for line in infile:
            newline = line.strip().split("\t")
            rsid = newline[5]
            try:
                # update the location dictionary if the rsid is already present
                temp_dic = {"chrom": newline[0],
                            "begin": int(newline[1]),
                            "end": int(newline[2]),
                            "ref_base": newline[3],
                            "alt_bases": [newline[4]]}
                # check if this location is already in the dict
                # append the new alternative base to the dict
                for snp in snp_locations[rsid]:
                    if ((snp["begin"] == temp_dic["begin"])
                        and (snp["end"] == temp_dic["end"])
                        and (snp["chrom"] == temp_dic["chrom"])
                            and (snp["ref_base"] == temp_dic["ref_base"])):
                        snp["alt_bases"].append(temp_dic["alt_bases"][0])
                        break
                else:
                    # add the snp dict if the location is different than what
                    # is present in the location dict.
                    snp_locations[rsid].append(temp_dic)
            except KeyError:
                # add the new rsid to location dict if it is not already there
                snp_locations[rsid] = [temp_dic]
                capture_types[rsid] = newline[6]
    # one reference location for each snp is required
    # alternative assambly chromosomes have an underscore in their names,
    # so that will be utilized to get the location in the orignal assembly,
    # i.e. the chromosome that does not have the underscore
    # (chr7 and not chr7_alt08)
    reference_snp_locations = {}
    problem_snps = []
    for s in snp_locations:
        if len(snp_locations[s]) == 1:
            reference_snp_locations[s] = snp_locations[s][0]
        else:
            for i in range(len(snp_locations[s])):
                if len(snp_locations[s][i]["chrom"].split("_")) == 1:
                    reference_snp_locations[s] = snp_locations[s][i]
                    break
            else:
                print(("Short chromosome name not found! "
                       "Please check the output list."))
                problem_snps.append(s)
        reference_snp_locations[s]["capture_type"] = capture_types[s]
    return reference_snp_locations, snp_locations


def gene_to_target(gene_list, species):
    target_coordinates = {}
    for gene in gene_list:
        e = get_exons(get_gene(gene,
                               get_file_locations()[species]["refgene"],
                               alternative_chr=1))
        try:
            target_coordinates[gene] = {"chrom": e["chrom"],
                                        "begin": e["begin"],
                                        "end": e["end"]}
        except KeyError:
            target_coordinates[gene] = {"chrom": np.nan,
                                        "begin": np.nan,
                                        "end": np.nan}
    return target_coordinates


def gene_to_target_exons(gene_list, species, exon_list):
    target_coordinates = {}
    for i in range(len(gene_list)):
        gene = gene_list[i]
        exons_wanted = exon_list[i]
        gene_exons = get_exons(get_gene(gene,
                               get_file_locations()[species]["refgene"],
                               alternative_chr=1))
        exons = gene_exons["exons"]
        if gene_exons["orientation"] == "-":
            exons.reverse()
        if exons_wanted == "all":
            for j in range(len(exons)):
                e = exons[j]
                tar_name = "-".join([gene, "exon", str(j)])
                target_coordinates[tar_name] = {"chrom": gene_exons["chrom"],
                                                "begin": e[0],
                                                "end": e[1]}
        else:
            for j in exons_wanted:
                try:
                    e = exons[j]
                    tar_name = "-".join(gene, "exon", str(j))
                    target_coordinates[tar_name] = {
                        "chrom": gene_exons["chrom"],
                        "begin": e[0],
                        "end": e[1]}
                except IndexError:
                    print(("Exon ", j, " does not exist for gene ", gene))
    return target_coordinates


def parse_alignment(reg_file):
    """ Create a rinfo dictionary from a rinfo file."""
    reg_dic = {}
    with open(reg_file, "r") as infile:
        for line in infile:
            if line.startswith("REGION"):
                newline = line.strip().split("\t")
                key1 = newline[1].split(":")[0]
                key2 = newline[1].split(":")[1]
                if key1 not in reg_dic:
                    reg_dic[key1] = {key2: {"copyname": newline[2],
                                            "chr": int(newline[3][3:]),
                                            "begin": int(newline[4]),
                                            "end": int(newline[5]),
                                            "ori": (newline[6] == "F")}}
                else:
                    reg_dic[key1][key2] = {"copyname": newline[2],
                                           "chr": int(newline[3][3:]),
                                           "begin": int(newline[4]),
                                           "end": int(newline[5]),
                                           "ori": (newline[6] == "F")}
    return reg_dic


def update_rinfo_file(rinfo_file, update_file, output_file):
    """Update a rinfo file with the lines provided in the update_file.

    This function will read all lines from a rinfo file and an update file.
    First two columns of rinfo files describe the parameters while the
    rest assign values. All the lines in the update file which share the
    first column with a line in the original file will replace that line
    in the original file. All other lines in the original file will remain.
    """
    # read the update file
    update_dict = {}
    with open(update_file) as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                update_dict[(newline[0], newline[1])] = line
    # read the rinfo file and update as appropriate
    with open(rinfo_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                line_key = (newline[0], newline[1])
                try:
                    outfile.write(update_dict[line_key])
                except KeyError:
                    outfile.write(line)
            else:
                outfile.write(line)


def get_target_coordinates(res_dir, species, capture_size,
                           coordinates_file=None, snps_file=None,
                           genes_file=None):
    """Extract MIP target coordinates from provided files."""
    capture_types = {}
    # Get target coordinates specified as genomic coordinates
    if coordinates_file is None:
        region_coordinates = {}
        coord_names = []
    else:
        coordinates_file = os.path.join(res_dir, coordinates_file)
        try:
            coord_df = pd.read_table(coordinates_file, index_col=False)
            coord_names = coord_df["Name"].tolist()
            coord_df.rename(columns={"Name": "name", "Chrom": "chrom",
                            "Start": "begin", "End": "end"}, inplace=True)
            region_coordinates = coord_df.set_index("name").to_dict(
                orient="index")
            # update capture types of targets
            for g in region_coordinates:
                if g not in capture_types:
                    capture_types[g] = region_coordinates[g]["CaptureType"]
        except IOError:
            print(("Target coordinates file {} could not be found.").format(
                (coordinates_file)))
            region_coordinates = {}
            coord_names = []

    # Get Gene target coordinates
    if genes_file is None:
        gene_coordinates = {}
        gene_names = []
    else:
        # get the alias file (gene name to gene id mapping) if available
        try:
            with open(get_file_locations()[species]["alias"]) as infile:
                alias = json.load(infile)
        except (KeyError, IOError):
            pass
        try:
            genes_file = os.path.join(res_dir, genes_file)
            genes_df = pd.read_table(genes_file, index_col=False)
            gene_names = genes_df["Gene"].tolist()
            genes = genes_df.set_index("Gene").to_dict(orient="index")
            gene_id_to_gene = {}
            gene_ids = []
            gene_coordinates = {}
            for g in genes:
                try:
                    if np.isnan(genes[g]["GeneID"]):
                        try:
                            gene_id = alias[g]
                            genes[g]["GeneID"] = gene_id
                        except KeyError:
                            print("""Alias for gene %s is not found.
                                Either provide a gene ID or use an alias
                                which is present in refgene file.""" % g)
                            continue
                        except NameError:
                            print(""" Gene ID is not provided for %s.
                                If gene name will be used to extract gene
                                ID an alias dictionary must be specified.
                                """ % g)
                            continue
                except TypeError:
                    pass
                gene_ids.append(genes[g]["GeneID"])
                gene_id_to_gene[genes[g]["GeneID"]] = g
                capture_types[g] = genes[g]["CaptureType"]
            gene_id_coordinates = gene_to_target(gene_ids, species)
            for gid in gene_id_coordinates:
                gene_coordinates[gene_id_to_gene[gid]] = gene_id_coordinates[
                    gid]
        except IOError:
            print(("Target genes file {} could not be found.").format(
                (genes_file)))
            gene_coordinates = {}
            gene_names = []

    if snps_file is None:
        snp_coordinates = {}
    else:
        # Get SNP target coordinates
        try:
            snps_file = os.path.join(res_dir, snps_file)
            snp_df = pd.read_table(snps_file, index_col=False,
                                   dtype={"Start": int, "End": int})
            snp_df.rename(columns={"Name": "name", "Chrom": "chrom",
                                   "Start": "begin", "End": "end"},
                          inplace=True)
            snp_coordinates = snp_df.set_index("name").to_dict(orient="index")
            for g in snp_coordinates:
                if g not in capture_types:
                    capture_types[g] = "targets"
        except IOError:
            print(("Target SNPs file {} could not be found.").format(
                (snps_file)))
            snp_coordinates = {}

    # merge coordinates dictionaries
    all_coordinates = {}
    all_coordinates.update(snp_coordinates)
    all_coordinates.update(gene_coordinates)
    all_coordinates.update(region_coordinates)

    # Fix names that has unwanted characters
    for c in list(all_coordinates.keys()):
        clist = []
        for ch in c:
            if ch.isalnum():
                clist.append(ch)
            else:
                clist.append("-")
        newc = "".join(clist)
        if newc != c:
            print("%s is replaced with %s" % (c, newc))
            all_coordinates[newc] = all_coordinates.pop(c)
            capture_types[newc] = capture_types.pop(c)
    target_regions, target_names = merge_coordinates(all_coordinates,
                                                     capture_size)
    # prioritize gene names ond  coordinate names  over snp or other names
    for t in list(target_names.keys()):
        for n in target_names[t]:
            if n in gene_names:
                target_names[n] = target_names.pop(t)
                target_regions[n] = target_regions.pop(t)
                break
            elif n in coord_names:
                target_names[n] = target_names.pop(t)
                target_regions[n] = target_regions.pop(t)
                break
    out_dict = {"target_regions": target_regions,
                "target_names": target_names,
                "capture_types": capture_types,
                "gene_names": gene_names,
                "snp_coordinates": snp_coordinates,
                "gene_coordinates": gene_coordinates,
                "region_coordinates": region_coordinates}

    return out_dict


def merge_coordinates(coordinates, capture_size):
    """Merge overlapping coordinates for MIP targets.

    Parameters
    ----------
    coordinates: python dictionary
        Coordinates to be merged in the form {target-name: {chrom: chrx,
        begin: start-coordinate, end: end-coordinate}, ..}
    capture_size: int
        Anticipated MIP capture size. If two regions are as close as 2 times
        this value, they will be merged.

    Returns
    -------
    target_coordinates: python dictionary
        merged coordinates dictionary
    target_names: python dictionary
        names of included targets in each merged region.
    """
    # create target regions to cover all snps
    # start by getting snps on same chromosome together
    chroms = {}
    for c in coordinates:
        chrom = coordinates[c]["chrom"]
        try:
            chroms[chrom].append([coordinates[c]["begin"],
                                  coordinates[c]["end"]])
        except KeyError:
            chroms[chrom] = [[coordinates[c]["begin"],
                              coordinates[c]["end"]]]
    # merge snps that are too close to get separate regions
    # the length should be twice the capture size
    merged_chroms = {}
    for c in chroms:
        merged_chroms[c] = merge_overlap(chroms[c], 2 * capture_size)
    # create regions for alignment
    # create target coordinate for each region
    target_coordinates = {}
    target_names = {}
    for c in merged_chroms:
        regions = merged_chroms[c]
        for reg in regions:
            targets_in_region = []
            for co in coordinates:
                if (coordinates[co]["chrom"] == c
                    and reg[0] <= coordinates[co]["begin"]
                        <= coordinates[co]["end"] <= reg[1]):
                    targets_in_region.append(co)
            region_name = targets_in_region[0]
            target_names[region_name] = targets_in_region
            r_start = reg[0]
            r_end = reg[1]
            target_coordinates[region_name] = [c, r_start, r_end]
    return target_coordinates, target_names


def create_target_fastas(res_dir, targets, species, flank):
    """ Create fasta files for a list of region coordinates provided as a dict
    in the form {target1: [chrx, start, end], target2: [chrx, start, end], ..},
    flank on both sides with the specified length. If beginning  coordinate is
    less than zero, reset the beginning coordinate to zero..
    """
    for t in list(targets.keys()):
        chrom = targets[t][0]
        begin = targets[t][1] - flank + 1
        if begin < 0:
            begin = 0
        end = targets[t][2] + flank
        rk = chrom + ":" + str(begin) + "-" + str(end)
        try:
            with open(os.path.join(res_dir, t + ".fa"), "w") as outfile:
                outfile.write(get_fasta(rk, species, header=t))
        except Exception as e:
            print(("Fasta file for {} could not be created, "
                   "due to error {}. It will be removed"
                   " from the target list.").format(t, e))
            targets.pop(t)
    return


def add_fasta_targets(res_dir, fasta_files, fasta_capture_type):
    fasta_sequences = {}
    capture_types = {}
    for f in fasta_files:
        f_file = os.path.join(res_dir, f)
        try:
            fasta_sequences.update(fasta_parser(f_file))
        except IOError:
            print(("Fasta file {} could not be found.").format(f_file))
    for f in list(fasta_sequences.keys()):
        flist = []
        for fch in f:
            if fch.isalnum():
                flist.append(fch)
            else:
                flist.append("-")
        newf = "".join(flist)
        if f != newf:
            print("%s is changed to %s." % (f, newf))
            fasta_sequences[newf] = fasta_sequences.pop(f)
        if newf not in capture_types:
            capture_types[newf] = fasta_capture_type
        with open(os.path.join(res_dir, newf + ".fa"), "w") as outfile:
            outfile.write(">" + newf + "\n" + fasta_sequences[newf] + "\n")
    return {"fasta_sequences": fasta_sequences, "capture_types": capture_types}


def set_genomic_target_alignment_options(target_regions, fasta_sequences,
                                         identity, coverage, flank):
    alignment_list = []
    fasta_list = list(fasta_sequences.keys()) + list(target_regions.keys())
    for t in fasta_list:
        temp_dict = {"gene_name": t, "identity": identity}
        try:
            target_size = target_regions[t][2] - target_regions[t][1]
            fasta_size = target_size + 2 * flank
        except KeyError:
            fasta_size = len(fasta_sequences[t])
        cover = round(coverage * 100 / fasta_size, 1)
        temp_dict["options"] = []
        if cover > 100:
            cover = 100
        temp_dict["coverage"] = cover
        if fasta_size < 100:
            temp_dict["options"].extend(["--notransition", "--step=10",
                                         "--ambiguous=iupac"])
        elif fasta_size < 1000:
            temp_dict["options"].extend(["--notransition", "--step=10",
                                         "--ambiguous=iupac"])
        elif fasta_size < 5000:
            temp_dict["options"].extend(["--notransition",
                                         "--step=" + str(int(fasta_size/10)),
                                         "--ambiguous=iupac"])
        else:
            temp_dict["options"].extend(["--notransition",
                                         "--step=" + str(int(fasta_size/10)),
                                         "--ambiguous=iupac"])
        alignment_list.append(temp_dict)
    return alignment_list


def align_region_multi(alignment_list, pro):
    """Parallelize a list of lastz alignments."""
    p = Pool(pro)
    p.map_async(align_region_worker, alignment_list)
    p.close()
    p.join()
    return


def align_region_worker(l):
    """Worker function for align_region_multi.

    Aligns a single fasta query file to a target fasta file. Both query
    and target fasta files  can be multi sequence files.
    """
    # get parameters from the input list
    # first item is the fasta file name, including file extension
    region_key = l[0]
    # second item holds the run directory for lastz
    resource_dir = l[1]
    # output file is the target name + ".al" where the alignment output
    # will be saved.
    output_file = l[2]
    # target fasta file is usually the reference genome
    target_fasta = l[3]
    # each action item will be appended to the target or query argument
    # within brackets. [unmask] and [multiple] are important target actions
    # unmask: allows starting alignments in masked(lowercase) parts of the
    # target multiple: indicates there are multiple sequences in the target
    # file (e.g. chromosomes, contigs)
    target_actions = l[4]
    # query file is always treated as a multiple sequence file
    # so there is no need for the multiple action
    query_actions = l[5]
    # percent cutoff value for identity/coverage of query to target. This only
    # affects reporting and not the alignment process itself.
    identity_cutoff = l[6]
    coverage_cutoff = l[7]
    # format of the output, follows --format: argument in lastz
    # if format is general, it should be followed by a comma separated list of
    # fields to output, e.g. general:name1,text1,name2,text2,diff,score would
    # seq of target, output the name of the query, sequence of the query, name
    # of the target, a string showing the alignment and the alignment score
    output_format = l[8]
    # additional options to pass to lastz
    options = l[9]
    query_fasta = os.path.join(resource_dir, region_key)
    # create target actions text
    if len(target_actions) > 0:
        target_act = "[" + ",".join(target_actions) + "]"
    else:
        target_act = ""
    # create query actions text
    if len(query_actions) > 0:
        query_act = "[" + ",".join(query_actions) + "]"
    else:
        query_act = ""
    # create the command list to pass to the processor
    comm = ["lastz_32",
            target_fasta + target_act,
            query_fasta + query_act,
            "--output=" + os.path.join(resource_dir, output_file),
            "--format=" + output_format,
            "--filter=identity:" + str(identity_cutoff),
            "--filter=coverage:" + str(coverage_cutoff)]
    # add any extra options to the end of the command
    comm.extend(options)
    # run the command using subprocess module
    subprocess.check_output(comm)
    return


def align_genes_for_design(fasta_list, res_dir,
                           alignment_types=["differences", "general"],
                           species="hs", num_processor=30):
    """Perform specified alignments given in an alignment dict.

    This functions is called from align_targets function for the initial
    target alignment to the reference genome.
    It align sequences given in an alignment dict which contains alignment
    specifics. Each entry in this dict must have a corresponding fasta file in
    the res_dir specified. The alignment is performed against the reference
    genome. This function merely prepares a list of commands to pass to
    align_genes_for_design_worker function to carry out alignments in
    parallel where multiple processors are available. Two types of alignment
    outputs will be generated; one "general" informative about the alignment
    such as where the alignment starts and ends, what is the percent identity,
    coverage etc. The second output is the differences between the aligned
    sequences, showing at which positions there are nucleotide changes and
    what the changes are.

    Parameters
    ----------
    fasta_list: list
        A list of dictionaries each of which contains specifics
        for a single alignment, such as the name of the fasta file, coverage
        and identity cut offs and any additional alignment parameters that are
        passed to LastZ.
    res_dir: str
        Path to working directory where input and output files are located.
    alignment_types: list
        List of alignment types to be performed. Only "general" and/or
        "differences" options are allowed.
    species: str
        Species whose reference genome will be used for alignment.
    num_processor: int
        Number of processors available for parallel processing.
    """
    region_list = []
    for gene_dict in fasta_list:
        gene_name = gene_dict["gene_name"]
        # percent cutoff value for identity/coverage of query to target.
        # This only affects reporting and not the alignment process itself.
        identity = gene_dict["identity"]
        coverage = gene_dict["coverage"]
        options = gene_dict["options"]
        # alignment target is the reference genome of the specified species.
        target = get_file_locations()[species]["fasta_genome"]
        # alignment output should have the following fields.
        # These are the bare minimum to be able to parse the alignment later.
        out_fields = ["name1", "strand1", "zstart1", "end1", "length1",
                      "name2", "strand2", "zstart2", "end2", "zstart2+",
                      "end2+", "length2", "identity", "coverage"]
        out_fields = ",".join(out_fields)
        gen_out = "general:" + out_fields
        # output fields for "differences" is fixed; it outputs the differences
        # between the aligned sequence and the target.
        dif_out = "differences"
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)
        # prepare a list of commands to feed to lastz for both alignment types
        # i.e.  "general" and "differences". Some of the additional parameters
        # we are supplying here are the target and query actions.
        # each action item will be appended to the target or query argument
        # within brackets. [unmask] and [multiple] are important target actions
        # unmask: allows starting alignments in masked(lowercase) parts of the
        # target multiple: indicates there are multiple sequences in the target
        # file (e.g. chromosomes, contigs)
        if "general" in alignment_types:
            al = [gene_name + ".fa", res_dir, gene_name + ".al", target,
                  ["multiple", "unmask", "nameparse=darkspace"],
                  ["unmask", "nameparse=darkspace"],
                  identity, coverage, gen_out, options]
            region_list.append(al)
        if "differences" in alignment_types:
            al = [gene_name + ".fa", res_dir, gene_name + ".differences",
                  target, ["multiple", "unmask", "nameparse=darkspace"],
                  ["unmask", "nameparse=darkspace"],
                  identity,  coverage, dif_out, options]
            region_list.append(al)
    align_region_multi(region_list, num_processor)
    return


def merge_alignments(resource_dir, fasta_list, output_prefix="merged"):
    """ Merge the results of "general" type lastZ alignments into a
    single file. This is used to process the alignment results from the
    align_genes_for_design function where target sequences are aligned
    against the reference genome.

    Parameters
    ----------
    resource_dir: str
        Path to working directory where the alignment outputs are.
    fasta_list: list
        A list of dictionaries each of which has the specifics for a single
        sequence alignment. It is used only to get alignment file names here.
    output_prefix: str
        Name for the output file. This will be appended by ".al" extension.
    """
    # create a list for each alignment type (general and differences)
    als_out = []
    with open(os.path.join(
            resource_dir, output_prefix + ".al"), "w") as alignment_file:
        for f in fasta_list:
            fnum = 0
            with open(os.path.join(resource_dir, f + ".al")) as alignment:
                linenum = 0
                for line in alignment:
                    if linenum > 0:
                        als_out.append(line.strip())
                    elif fnum == 0:
                        als_out.append(line.strip())
                        linenum += 1
                    else:
                        linenum += 1
            fnum += 0
        alignment_file.write("\n".join(als_out))
    return


def merge_alignment_diffs(resource_dir, fasta_list, output_prefix="merged"):
    """ Merge the results of "differences" type lastZ alignments into a
    single file. This is used to process the alignment results from the
    align_genes_for_design function where target sequences are aligned
    against the reference genome.

    Parameters
    ----------
    resource_dir: str
        Path to working directory where the alignment outputs are.
    fasta_list: list
        A list of dictionaries each of which has the specifics for a single
        sequence alignment. It is used only to get alignment file names here.
    output_prefix: str
        Name for the output file. This will be appended by ".al" extension.
    """
    # create a list for each alignment type (general and differences)
    diffs_out = []
    with open(os.path.join(
            resource_dir, output_prefix + ".differences"), "w") as diff_file:
        for f in fasta_list:
            fnum = 0
            with open(os.path.join(resource_dir, f + ".differences")) as diffs:
                for d in diffs:
                    diffs_out.append(d.strip())
            fnum += 0
        diff_file.write("\n".join(diffs_out))
    return


def alignment_parser(wdir, name, spacer=0, gene_names=[]):
    """ Parse merged genome alignment results file which is generated by
    align_genes_for_design function to align design targets to reference
    genomes. One query (target region) may have multiple alignments to the
    genome.

    Parameters
    ----------
    wdir: str
        Path to working directory
    name: str
        File name for the merged alignment file
    spacer: int
        Spacer length to use when merging overlapping regions. If two regions
        are not overlapping but the distance between them is smaller than the
        spacer, they will be merged.

    Returns
    -------
    A list of dictionaries:
    target_regions: merged genomic coordinates for grouped targets.
        This dictionary is used as the final target regions.
        For example: {r1: [[chr1, 100, 200], [chr3, 30, 300]],
                      r3: [chr4, 0, 300]]}
    region_names: names for each region.
        For example: {r1: [r1, r2], r3: [r3]}
    imperfect_aligners: names of the target regions for which a perfect
        alignment to the reference genome has not been found.
    """
    alignment_dict = {}
    # open alignment files
    with open(os.path.join(wdir, name + ".al")) as infile:
        # each line in the file is a separate alignment for which we'll
        # prepare a dictionary.
        for line in infile:
            newline = line.strip().split("\t")
            # first line has column names
            if line.startswith("#"):
                colnames = [newline[0][1:]]
                colnames.extend(newline[1:])
            else:
                temp_dict = {}
                for i in range(len(colnames)):
                    col = colnames[i]
                    value = newline[i]
                    temp_dict[col] = value
                query_name = temp_dict["name2"]
                try:
                    alignment_dict[query_name].append(temp_dict)
                except KeyError:
                    alignment_dict[query_name] = [temp_dict]
    # go through each target sequence and each alignment for that
    # target to where in the genome it was aligned to.
    aligned_regions = {}
    for query in alignment_dict:
        aligned_regions[query] = []
        for a in alignment_dict[query]:
            chrom = a["name1"]
            begin = int(a["zstart1"])
            end = int(a["end1"])
            aligned_regions[query].append([chrom, begin, end])
    # check for overlapping alignments. These can be the same target aligning
    # to overlapping regions in the genome (internal duplications) or
    # different targets aligning to the same (or overlapping) regions in the
    # genome (paralogus sequences).
    # overlapping regions will be grouped together to form the final target
    # regions for probe design.
    overlaps = {}
    for q1 in aligned_regions:
        # each target will have itself as overlapping
        overlaps[q1] = [q1]
        # get the genomic regions q1 was aligned to
        reg1 = aligned_regions[q1]
        # go through each region
        for r1 in reg1:
            # check overlap with other target regions
            for q2 in aligned_regions:
                if q1 == q2:
                    continue
                reg2 = aligned_regions[q2]
                for r2 in reg2:
                    if check_overlap(r1, r2, spacer):
                        overlaps[q1].append(q2)
                        break
    # go through the overlaps and remove the overlapping overlaps
    # e.g. if a overlaps b, b overlaps a also. We'll have {a: [a,b], b: [b, a]}
    # in the overlaps dict. We want only one of these, so reduce to {a:[a, b]}
    overlap_found = True
    # place a failsafe counter to avoid unforseen infinite loops
    exit_counter = 0
    while (overlap_found and (exit_counter < 10000)):
        overlap_found = False
        for o in list(overlaps.keys()):
            # check if o is still in the overlaps and has not been removed
            if o in overlaps:
                val = overlaps[o]
                # get the overlapping regions for "val" and add them
                # to overlapping regions for "o", then remove "val"
                for v in val:
                    if (v in overlaps) and (o in overlaps) and (o != v):
                        overlaps[o].extend(overlaps[v])
                        overlaps.pop(v)
                        overlap_found = True
    if exit_counter > 9999:
        print("Overlap removal while loop limit is reached.")
    # clean up overlapping region lists by removing duplicates.
    for o in list(overlaps.keys()):
        overlaps[o] = sorted(list(set(overlaps[o])))
    ##########################################################################
    # create a new dictionary for target regions.
    # for each target group in overlaps, we'll have genomic coordinates
    # that will be used as final targets.
    ##########################################################################
    # group regions according to their chromosomes
    separated_regions = {}
    for o in overlaps:
        sep = separated_regions[o] = {}
        for g in overlaps[o]:
            regs = aligned_regions[g]
            for r in regs:
                try:
                    sep[r[0]].append(r[1:])
                except KeyError:
                    sep[r[0]] = [r[1:]]
    # merge each overlapping region
    separated_merged_regions = {}
    for s in separated_regions:
        merged_sep = separated_merged_regions[s] = {}
        for chrom in separated_regions[s]:
            merged_region = merge_overlap(separated_regions[s][chrom])
            merged_sep[chrom] = merged_region
    ###########################################
    # organize target regions, assign region names based on the original
    # target names. Assign a reference target.
    ###########################################
    # sort target regions based on the length of
    # chromosome name and the length of region. Sort is based on the region
    # size and chromosome name is used as a tie-breaker
    # to distinguish alternate contigs and not use them as reference, but
    # it is not absolutely necessary and it would not behave as expected
    # when chromosome names do not follow that convention, i.e, chr6 and
    # chr6_altXYZ.
    for ar in list(aligned_regions.keys()):
        regs = aligned_regions[ar]
        for r in regs:
            r.append(0 - len(r[0]))
            r.append(r[2] - r[1] + 1)
        aligned_regions[ar] = sorted(regs, key=itemgetter(4, 3),
                                     reverse=True)
    target_regions = {}
    region_names = {}
    regions = separated_merged_regions
    for r in regions:
        target_regions[r] = []
        for chrom in regions[r]:
            for l in regions[r][chrom]:
                temp_region = [chrom]
                temp_region.extend(l)
                temp_region.append(-len(chrom))
                temp_region.append(l[1] - l[0])
                target_regions[r].append(temp_region)
        # sort target regions per target group based on the length of
        # chromosome name and the length of region. Chromosome name is used
        # to distinguish alternate contigs and not use them as reference, but
        # it is not absolutely necessary and it would not behave as expected
        # when chromosome names do not follow that convention, i.e, chr6 and
        # chr6_altXYZ
        target_regions[r] = sorted(target_regions[r], key=itemgetter(4, 3),
                                   reverse=True)
        # assign names to grouped targets
        reg_names = []
        # for each region we go back to individual region alignments and see
        # if the individual alignment overlaps with this region. If it does
        # we use the individual regions name for this region within the group.
        for i in range(len(target_regions[r])):
            reg = target_regions[r][i]
            reg_chrom = reg[0]
            reg_begin = reg[1]
            reg_end = reg[2]
            for c in aligned_regions:
                main_region = aligned_regions[c][0]
                if (reg_chrom == main_region[0]
                        and reg_begin <= main_region[1]
                        and reg_end >= main_region[2]):
                    reg_names.append(c)
                    break
            else:
                reg_names.append("na")
        # assign a reference region for each group based on gene names provided
        # this is mainly to used to have better names for regions. For example,
        # if a gene is a target as well as a snp, we would like the gene name
        # to be the name of the group as opposed to the SNP's name.
        ref_found = False
        for g in gene_names:
            if g in reg_names:
                ref_found = True
                ref_index = reg_names.index(g)
                ref_name = g
                break
        if not ref_found:
            ref_name = r
            ref_index = 0
        ref_region = target_regions[r].pop(ref_index)
        reg_names.pop(ref_index)
        target_regions[r] = [ref_region] + target_regions[r]
        reg_names = [ref_name] + reg_names
        region_names[ref_name] = reg_names
        target_regions[reg_names[0]] = target_regions.pop(r)
        overlaps[reg_names[0]] = overlaps.pop(r)
    # after the alignments are done, some regions will not have proper names
    # and some will have "na". We'll change those to avoid repeating
    # names.
    for r in list(region_names.keys()):
        rnames = region_names[r]
        nnames = []
        rn_counts = {}
        for rn in rnames:
            rnc = rnames.count(rn)
            rn_counts[rn] = {"total_count": rnc,
                             "used_count": 0}
        for rn in rnames:
            if rn_counts[rn]["total_count"] > 1:
                nnames.append(rn + "-" + str(rn_counts[rn]["used_count"]))
                rn_counts[rn]["used_count"] += 1
            else:
                nnames.append(rn)
        region_names[r] = nnames
    # find target regions that could not be perfectly aligned to the genome
    # these are usually extragenomic sequences supplied in fasa files, such as
    # certain TCR haplotypes.
    imperfect_aligners = []
    for r in alignment_dict:
        best_score = 0
        alignments = alignment_dict[r]
        for a in alignments:
            cov = int(a["covPct"].split(".")[0])
            idt = int(a["idPct"].split(".")[0])
            score = cov * idt
            if score > best_score:
                best_score = score
        if best_score != 10000:
            imperfect_aligners.append(r)
    return [target_regions, region_names, imperfect_aligners, aligned_regions,
            overlaps]


def set_intra_alignment_options(target_regions, identity, coverage,
                                max_allowed_indel_size):
    """Set lastZ alignment options for intraparalog_aligner function."""
    alignment_options_dict = {}
    for t in target_regions:
        temp_dict = {"gene_name": t, "identity": identity}
        reference_len = target_regions[t][0][-1]
        small_target = 0
        for r in target_regions[t]:
            if r[-1] < coverage:
                small_target += 1
                try:
                    smallest_target = min([smallest_target, r[-1]])
                except NameError:
                    smallest_target = int(r[-1])
        if small_target > 0:
            print(("{} targets within {} are smaller than intra_coverage"
                   " value. This means that those targets will not be aligned."
                   " Smallest target's length was {}. Set intra_coverage"
                   " to a value smaller than this value to align all regions."
                   ).format(small_target, t, smallest_target))
        cover = round(coverage * 100 / reference_len, 1)
        gap_open_penalty = 400
        gap_extend_penalty = 30
        ydrop = max_allowed_indel_size * gap_extend_penalty + gap_open_penalty
        alignment_opts = ["--ydrop=" + str(ydrop), "--notransition",
                          "--ambiguous=iupac", "--noytrim"]
        temp_dict["options"] = alignment_opts
        if cover > 100:
            cover = 100
        temp_dict["coverage"] = cover
        alignment_options_dict[t] = temp_dict
    return alignment_options_dict


def intraparalog_aligner(resource_dir,
                         target_regions,
                         region_names,
                         imperfect_aligners,
                         fasta_sequences,
                         species,
                         num_process,
                         alignment_options_dict={}):
    """Align all regions within a target group.

    Align all regions within a target group to the region selected
    as the reference region.

    Returns
    -------
    Returns nothing. It creates query.fa target.fa and .aligned files for each
    target region group. These alignment have no genomic coordinates, so
    all coordinates are relative to the given sequence. Also, the region names
    are indicated as the reference gene name + copy name as this is originally
    intended for use in paralog genes.
    """
    alignment_commands = []
    out_fields = "name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,"
    out_fields = out_fields + "end2,zstart2+,end2+,length2,identity,coverage"
    gen_out = "general:" + out_fields
    diff_out = "differences"
    for t in target_regions:
        alignment_options = alignment_options_dict[t]["options"]
        identity = alignment_options_dict[t]["identity"]
        coverage = alignment_options_dict[t]["coverage"]
        tar_regs = target_regions[t]
        # create a fasta file for the reference copy (or reference region)
        target_keys = [tr[0] + ":" + str(tr[1] + 1)
                       + "-" + str(tr[2]) for tr in tar_regs]
        query_key = target_keys[0]
        with open(os.path.join(resource_dir, t + ".query.fa"), "w") as outfile:
            outfile.write(">" + t + "_ref\n")
            outfile.write(get_sequence(query_key, species))
        # create a fasta file that includes all target regions within a group.
        with open(os.path.join(
                resource_dir, t + ".targets.fa"), "w") as outfile:
            outfile_list = []
            for i in range(len(target_keys)):
                k = target_keys[i]
                cname = "_C" + str(i)
                outfile_list.append(">" + t + cname)
                outfile_list.append(get_sequence(k, species))
            # add extragenomic (i.e. imperfect_aligners)
            ols = region_names[t]
            o_count = 0
            for o in ols:
                if o in imperfect_aligners:
                    outfile_list.append(">" + t + "_X" + str(o_count))
                    outfile_list.append(fasta_sequences[o])
                    o_count += 1
            outfile.write("\n".join(outfile_list))
        comm = [t + ".query.fa", resource_dir, t + ".aligned",
                os.path.join(resource_dir, t + ".targets.fa"),
                ["multiple", "unmask", "nameparse=darkspace"],
                ["unmask", "nameparse=darkspace"],
                identity, coverage, gen_out,
                alignment_options, species]
        alignment_commands.append(comm)
        comm = [t + ".query.fa", resource_dir,
                t + ".differences",
                os.path.join(resource_dir, t + ".targets.fa"),
                ["multiple", "unmask", "nameparse=darkspace"],
                ["unmask", "nameparse=darkspace"],
                identity, coverage,
                diff_out, alignment_options, species]
        alignment_commands.append(comm)
    return align_region_multi(alignment_commands, num_process)


def intra_alignment_checker(family_name, res_dir, target_regions,
                            region_names):
    """
    Parse intraparalog_aligner results.

    Following a within group alignment, check if any individual region
    within the group has multiple aligned parts. If found, split that region
    into multiple regions to be re-aligned by intraparalog_aligner.
    """
    alignment_file = family_name + ".aligned"
    new_regions = {}
    with open(os.path.join(res_dir, alignment_file), "r") as alignment:
        for line in alignment:
            # extract the column names from the first line
            if line.startswith("#"):
                newline = line.strip().split("\t")
                newline[0] = newline[0][1:]
                colnames = list(newline)
            # assign values of each column for each alignment
            else:
                newline = line.strip().split("\t")
                temp_dict = {}
                for i in range(len(colnames)):
                    temp_dict[colnames[i]] = newline[i]
                alignment_id = temp_dict["name1"]
                ci = alignment_id.split("_")[-1]
                ct = ci[0]
                if ct == "C":
                    cn = int(ci[1:])
                    tr = target_regions[cn]
                    start = tr[1] + int(temp_dict["zstart1"])
                    end = tr[1] + int(temp_dict["end1"])
                    size = end - start + 1
                    try:
                        new_regions[cn].append([tr[0], start, end,
                                               0 - len(tr[0]), size])
                    except KeyError:
                        new_regions[cn] = [[tr[0], start, end,
                                           0 - len(tr[0]), size]]
    # check if any paralog is missing after aligning to the reference copy
    targeted_copies = list(range(len(target_regions)))
    missing_copies = set(targeted_copies).difference(new_regions.keys())
    if len(missing_copies) > 0:
        print(("Paralog copies {} were not successfully aligned to "
               "the reference copy for the target {}. You may consider "
               "relaxing the alignment filters '--local-coverage' "
               "and '--local-identity'").format(
                   ", ".join(map(str, sorted(missing_copies))), family_name))
    ret_regions = []
    rnames = []
    for ci in sorted(new_regions):
        ret_regions.extend(sorted(new_regions[ci]))
        if len(new_regions[ci]) > 1:
            print(("Paralog copy {} for target region {} was aligned "
                   "to the reference copy multiple times. This copy will "
                   "be treated as multiple independent paralog copies and "
                   "realigned to the reference copy as separate "
                   "targets.").format(ci, family_name))
            for i in range(len(new_regions[ci])):
                rnames.append(region_names[ci] + "-" + str(i))
        else:
            rnames.append(region_names[ci])
    return [ret_regions, rnames]


def align_paralogs(res_dir, target_regions, region_names, imperfect_aligners,
                   fasta_sequences, species, identity, coverage,
                   max_allowed_indel_size, num_process):
    alignment_options = set_intra_alignment_options(
        target_regions, identity, coverage, max_allowed_indel_size)
    intraparalog_aligner(res_dir, target_regions, region_names,
                         imperfect_aligners, fasta_sequences, species,
                         num_process, alignment_options)
    for r in target_regions.keys():
        ntr = intra_alignment_checker(r, res_dir, target_regions[r],
                                      region_names[r])
        target_regions[r] = ntr[0]
        region_names[r] = ntr[1]
    alignment_options = set_intra_alignment_options(
        target_regions, identity, coverage, max_allowed_indel_size)
    intraparalog_aligner(res_dir, target_regions, region_names,
                         imperfect_aligners, fasta_sequences, species,
                         num_process, alignment_options)


def get_missed_targets(original_target_regions, target_regions,
                       aligned_regions, min_target_size, flank, capture_types):
    org_chroms = {}
    new_chroms = {}
    for o in original_target_regions:
        org_regs = original_target_regions[o]
        for org in org_regs:
            try:
                org_chroms[org[0]].append(org[1:3])
            except KeyError:
                org_chroms[org[0]] = [org[1:3]]
        new_regs = target_regions[o]
        for nrg in new_regs:
            try:
                new_chroms[nrg[0]].append(nrg[1:3])
            except KeyError:
                new_chroms[nrg[0]] = [nrg[1:3]]
    uncovered_chroms = {}
    for chrom in org_chroms:
        try:
            uncov = subtract_overlap(org_chroms[chrom], new_chroms[chrom])
            if len(uncov) > 0:
                uncovered_chroms[chrom] = uncov
        except KeyError:
            uncovered_chroms[chrom] = org_chroms[chrom]
    not_aligned_coordinates = {}
    for ar in aligned_regions:
        main_region = aligned_regions[ar][0]
        extra_count = 0
        for uc in uncovered_chroms:
            unc_regs = uncovered_chroms[uc]
            for ur in unc_regs:
                if len(overlap(main_region[1:3], ur)) > 0:
                    not_aligned_coordinates[
                        ar + "-extra-" + str(extra_count)
                    ] = {"chrom": uc,
                         "begin": ur[0],
                         "end": ur[1]}
    missed_target_regions, missed_target_names = merge_coordinates(
        not_aligned_coordinates, flank)
    for t in list(missed_target_regions.keys()):
        target_size = (missed_target_regions[t][-1]
                       - missed_target_regions[t][-2] + 1)
        if target_size < min_target_size:
            missed_target_regions.pop(t)
            missed_target_names.pop(t)
    missed_capt_types = {}
    for t in missed_target_names:
        try:
            missed_capt_types[t] = capture_types[t.split("extra")[0][:-1]]
        except KeyError:
            print(("Capture type not found for {}."
                   " Setting capture type to 'whole'").format(t))
            missed_capt_types[t] = "whole"
    return [missed_target_regions, missed_target_names, missed_capt_types]


def align_targets(res_dir, target_regions, species, flank, fasta_files,
                  fasta_capture_type, genome_identity, genome_coverage,
                  num_process, gene_names, max_allowed_indel_size,
                  intra_identity, intra_coverage, capture_types,
                  min_target_size, merge_distance, savefile):
    # create fasta files for each target coordinate
    create_target_fastas(res_dir, target_regions, species, flank)

    if fasta_files is None:
        fasta_sequences = fasta_capture_types = {}
    else:
        # add target sequences provided by fasta files
        fasta_targets = add_fasta_targets(
            res_dir, fasta_files, fasta_capture_type=fasta_capture_type)
        fasta_sequences = fasta_targets["fasta_sequences"]
        fasta_capture_types = fasta_targets["capture_types"]
    capture_types.update(fasta_capture_types)

    # create a list of target names from all sources
    targets_list = (list(target_regions.keys())
                    + list(fasta_sequences.keys()))

    # align target sequences to reference genome
    # create alignment options
    genomic_alignment_list = set_genomic_target_alignment_options(
        target_regions, fasta_sequences, genome_identity, genome_coverage,
        flank)
    # perform genome alignment
    align_genes_for_design(genomic_alignment_list, res_dir,
                           alignment_types="general", species=species,
                           num_processor=num_process)

    # merge all alignment files
    merge_alignments(res_dir, targets_list, output_prefix="merged")

    # parse genome alignment file
    # negative merge_distance values keep the target regions separate
    # even if they overlap. Positive values lead to merging targets.
    # However, the alignments are already carried out with flanking
    # sequence so increasing that merge distance is avoided by setting the
    # merge_distance 0 here for positive values.
    if merge_distance > 0:
        merge_distance = 0
    genome_alignment = alignment_parser(res_dir, "merged",
                                        spacer=merge_distance,
                                        gene_names=gene_names)
    target_regions = copy.deepcopy(genome_alignment[0])
    region_names = copy.deepcopy(genome_alignment[1])
    imperfect_aligners = copy.deepcopy(genome_alignment[2])
    aligned_regions = copy.deepcopy(genome_alignment[3])
    overlaps = copy.deepcopy(genome_alignment[4])

    # align sequences within target groups (paralog sequences)
    align_paralogs(res_dir, target_regions, region_names, imperfect_aligners,
                   fasta_sequences, species, intra_identity, intra_coverage,
                   max_allowed_indel_size, num_process)

    # compare original target_regions to the final target regions
    # to determine if any region is missing due to alignments performed
    original_target_regions = genome_alignment[0]
    missed_target_regions, missed_target_names, missed_capture_types = (
        get_missed_targets(original_target_regions, target_regions,
                           aligned_regions, min_target_size, flank,
                           capture_types))

    out_dict = {"original_target_regions": genome_alignment[0],
                "original_region_names": genome_alignment[1],
                "original_imperfect_aligners": genome_alignment[2],
                "original_aligned_regions": genome_alignment[3],
                "original_overlaps": genome_alignment[4],
                "target_regions": target_regions,
                "region_names": region_names,
                "aligned_regions": aligned_regions,
                "capture_types": capture_types,
                "imperfect_aligners": imperfect_aligners,
                "overlaps": overlaps,
                "missed_target_regions": missed_target_regions,
                "missed_target_names": missed_target_names,
                "missed_capture_types": missed_capture_types}
    with open(os.path.join(res_dir, savefile), "w") as outfile:
        json.dump(out_dict, outfile, indent=1)
    return out_dict


def alignment_mapper(family_name, res_dir):
    """Create a coordinate map of within group alignments."""
    alignment_file = family_name + ".aligned"
    difference_file = family_name + ".differences"
    with open(os.path.join(res_dir, alignment_file), "r") as alignment, open(
            os.path.join(res_dir, difference_file), "r") as difference:
        # create an alignment dictionary for each region that a query
        # aligns to these correspond to each line in the alignment file
        # and thus, are relative coordinates.
        alignment_dic = {}
        for line in alignment:
            # extract the column names from the first line
            if line.startswith("#"):
                newline = line.strip().split("\t")
                newline[0] = newline[0][1:]
                colnames = list(newline)
            # assign values of each column for each alignment
            else:
                newline = line.strip().split("\t")
                temp_dict = {"differences": []}
                for i in range(len(colnames)):
                    temp_dict[colnames[i]] = newline[i]
                alignment_id = temp_dict["name1"]
                if alignment_id in alignment_dic:
                    print(("{} aligned to the reference copy multiple times. "
                           "Only the first alignment will be used for "
                           "coordinate mapping.").format(alignment_id))
                    continue
                alignment_dic[alignment_id] = temp_dict
                cov = float(alignment_dic[alignment_id]["covPct"][:-1])
                idt = float(alignment_dic[alignment_id]["idPct"][:-1])
                alignment_dic[alignment_id]["score"] = np.mean([idt, cov])
        # differences file is a continuous file for all alignments
        # extract differences for each alignment
        for line in difference:
            newline = line.strip().split("\t")
            dname = newline[0]
            alignment_dic[dname]["differences"].append(newline[:-2])
        # map each position in each alignment to the query
        for a in alignment_dic:
            snps = alignment_dic[a]["snps"] = {}
            co = alignment_dic[a]["coordinates"] = {}
            rev_co = alignment_dic[a]["reverse_coordinates"] = {}
            # if alignment on reverse strand
            if alignment_dic[a]["strand2"] == "-":
                # genomic coordinate of target start
                # this position is zstart2+ away from query end
                # (when it is a - alignment)
                al_start = int(alignment_dic[a]["zstart1"])
                query_plus_end = int(alignment_dic[a]["end2+"])
                # assign start to the first key of the coord dictionary
                first_key = query_plus_end - 1
                co[first_key] = al_start
                rev_co[al_start] = first_key
                last_key = first_key
                inserted = 0
                for d in alignment_dic[a]["differences"]:
                    # start/end coordinates of diff relative to the query
                    diff_start = int(d[6])
                    diff_end = int(d[7])
                    query_length = int(d[9])
                    # for each diff, fill in the coordinates
                    # between the last_key in the coord dic and
                    # start_key - diff start
                    for j in range(last_key - 1, query_length
                                   - diff_start - 1, -1):
                        # j decreases by one, starting from the last
                        # available key the value will be 1 more than the
                        # previous key (j+1)
                        if j == last_key - 1:
                            co[j] = round(co[j + 1] - 0.1) + 1 + inserted
                        else:
                            co[j] = round(co[j + 1] - 0.1) + 1
                        rev_co[co[j]] = j
                    # current last key is now first_key - diff_start
                    last_key = query_length - diff_start - 1
                    query_diff_end = last_key + 1
                    # genomic coordinate of target at diff start
                    tar_start = int(d[1])
                    # genomic coordinate of target at diff end
                    tar_end = int(d[2])
                    # if end and start are the same, there is a deletion
                    # in target compared to query
                    # all nucleotides from diff start to diff end will have
                    # the same coordinate
                    if tar_start == tar_end:
                        inserted = 0
                        for i in range(diff_end - diff_start):
                            co[last_key - i] = tar_start - 0.5
                        last_key -= diff_end - diff_start - 1
                    # in cases of deletion in query, only rev_co will be
                    # updated
                    elif diff_start == diff_end:
                        inserted = 0
                        for i in range(tar_end - tar_start):
                            rev_co[co[last_key + 1] + i + 1] = (
                                last_key + 0.5)
                            inserted += 1
                        last_key += 1
                    # last_key will be mapped to target start
                    # if there is only a SNP and no indel
                    else:
                        inserted = 0
                        co[last_key] = tar_start
                        rev_co[tar_start] = last_key
                    query_diff_start = last_key
                    diff_key = str(query_diff_start) + "-" + str(
                        query_diff_end)
                    snps[diff_key] = {"chrom": d[0],
                                      "target_begin": int(d[1]),
                                      "target_end": int(d[2]),
                                      "target_orientation": d[3],
                                      "query_start": diff_start,
                                      "query_end": diff_end,
                                      "query_orientation": d[8],
                                      "target_base": d[10],
                                      "query_base": d[11]}
                # fill in the coordinates between last diff
                # and the alignment end
                query_plus_start = int(alignment_dic[a]["zstart2+"])
                for k in range(last_key - 1, query_plus_start - 1, -1):
                    co[k] = round(co[k+1] - 0.1) + 1
                    rev_co[co[k]] = k
            # when the alignment is on the forward strand
            else:
                # where on target sequence the alignment starts
                tar_start = int(alignment_dic[a]["zstart1"])
                # where in the query sequence the alinment starts
                q_start = int(alignment_dic[a]["zstart2"])
                co[q_start] = tar_start
                rev_co[tar_start] = q_start
                # last key used is q_start, last key is updated each time
                # something is added to the coordinate dict.
                last_key = first_key = q_start
                inserted = 0
                for d in alignment_dic[a]["differences"]:
                    # where on query sequence the difference starts and
                    # ends
                    diff_start = int(d[6])
                    diff_end = int(d[7])
                    diff_key = d[6] + "-" + d[7]
                    query_length = d[9]
                    snps[diff_key] = {"chrom": d[0],
                                      "target_begin": int(d[1]),
                                      "target_end": int(d[2]),
                                      "target_orientation": d[3],
                                      "query_start": diff_start,
                                      "query_end": diff_end,
                                      "query_orientation": d[8],
                                      "target_base": d[10],
                                      "query_base": d[11]}
                    # from the last key to the diff start the query and
                    # target sequences are the same in length and co dict
                    # is filled so
                    for i in range(last_key + 1, diff_start):
                        if i == last_key + 1:
                            co[i] = round(co[i-1] - 0.1) + 1 + inserted
                            inserted = 0
                        else:
                            co[i] = round(co[i-1] - 0.1) + 1
                        rev_co[co[i]] = i
                    # update last used key in co dict
                    last_key = diff_start
                    # genomic coordinate of target at diff start
                    tar_start = int(d[1])
                    # genomic coordinate of target at diff end
                    tar_end = int(d[2])
                    # if end and start are the same, there is a deletion
                    # in target compared to query
                    # all nucleotides from diff start to diff end will have
                    # the same coordinate
                    if tar_start == tar_end:
                        inserted = 0
                        for i in range(diff_end - diff_start):
                            co[last_key + i] = tar_start - 0.5
                        last_key += diff_end - diff_start - 1
                    # in cases of deletion in query (insertion in target)
                    # position will be mapped to the target end coordinate
                    elif diff_start == diff_end:
                        inserted = 0
                        for i in range(tar_end - tar_start):
                            rev_co[co[last_key - 1] + 1 + i] = (
                                last_key - 0.5)
                            inserted += 1
                        last_key -= 1
                    # if there is no indel
                    # last_key will be mapped to target start
                    else:
                        inserted = 0
                        co[last_key] = tar_start
                        rev_co[tar_start] = last_key
                # fill in the coordinates between last diff
                # and the alignment end
                q_end = int(alignment_dic[a]["end2"])
                for k in range(last_key + 1, q_end):
                    co[k] = round(co[k-1] - 0.1) + 1
                    rev_co[co[k]] = k
    return alignment_dic


###############################################################
# Design related functions
###############################################################


def order_mips(mip_info, design_name, res_dir):
    mip_sequences = []
    for g in sorted(mip_info):
        for m in sorted(mip_info[g]["mips"]):
            minfo = mip_info[g]["mips"][m]["mip_dic"]["mip_information"]
            for c in minfo:
                s = minfo[c]["SEQUENCE"]
                n = m + "_" + c
                num = int(m.split("_")[-1][3:])
                mip_sequences.append([n, s, g, num, m, c])
        if len(mip_info[g]["mips"]) == 0:
            mip_info.pop(g)
    mip_sequences = sorted(mip_sequences, key=itemgetter(2, 3))
    print("%d probes will be ordered." % len(mip_sequences))
    # Check for probes that have the same sequence
    sequence_only = [i[1].upper() for i in mip_sequences]
    for s in sequence_only:
        if sequence_only.count(s) > 1:
            print("At least two probes share the sequence %s" % s)
    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    columns = list(range(1, 13))
    for i in range(len(mip_sequences)):
        m = mip_sequences[i]
        plate = i/96
        pl_pos = i % 96
        col = columns[pl_pos % 12]
        row = rows[pl_pos/12]
        m.extend([row, col, plate])
    for i in range(len(mip_sequences)):
        m = mip_sequences[i]
        s = list(m[1])
        N_found = False
        for j in s:
            if s[j] == "N":
                if N_found:
                    s[j] == "(N)"
                else:
                    N_found = True
                    s[j] == "(N:25252525)"
        m.append("".join(s))
    order_dict = {}
    for i in range(len(mip_sequences)):
        m = mip_sequences[i]
        pl = m[-2]
        pl_name = design_name + "_" + str(pl)
        try:
            order_dict[pl_name].append(m)
        except KeyError:
            order_dict[pl_name] = [m]
    for o in order_dict:
        with open(os.path.join(res_dir, o), "w") as outfile:
            outfile_list = ["\t".join(["WellPosition", "Name", "Sequence"])]
            plate_mips = order_dict[o]
            for m in plate_mips:
                wp = m[-4] + str(m[-3])
                outfile_list.append("\t".join([wp, m[0], m[-1]]))
            outfile.write("\n".join(outfile_list))
    return


def create_dirs(dir_name):
    """ create subdirectory names for a given dir,
    to be used by os.makedirs, Return a list of
    subdirectory names."""
    primer3_input_DIR = dir_name + "/primer3_input_files/"
    primer3_output_DIR = dir_name + "/primer3_output_files/"
    bowtie2_input_DIR = dir_name + "/bowtie2_input/"
    bowtie2_output_DIR = dir_name + "/bowtie2_output/"
    mfold_input_DIR = dir_name + "/mfold_input/"
    mfold_output_DIR = dir_name + "/mfold_output/"
    return [primer3_input_DIR, primer3_output_DIR, bowtie2_input_DIR,
            bowtie2_output_DIR, mfold_input_DIR, mfold_output_DIR]


def get_snps(region, snp_file):
    """ Take a region string and a  tabix'ed snp file,
    return a list of snps which are lists of
    tab delimited information from the snp file. """
    # extract snps using tabix, in tab separated lines
    snp_temp = subprocess.check_output(["tabix", snp_file, region]).decode(
        "UTF-8"
    )
    # split the lines (each SNP)
    snps_split = snp_temp.split("\n")
    # add each snp in the region to a list
    # as lists of
    snps = []
    for line in snps_split:
        snp = line.split('\t')
        snps.append(snp)
    # remove last item which is coming from the new line at the end
    del snps[-1]
    return snps


def get_vcf_snps(region, snp_file):
    """ Take a region string and a tabix'ed snp file,
    return a list of snps which are lists of
    tab delimited information from the snp file. """
    # extract snps using tabix, in tab separated lines
    snp_temp = subprocess.check_output(["bcftools", "view", "-H", "-G", "-r",
                                        region, snp_file]).decode("UTF-8")
    # split the lines (each SNP)
    snps_split = snp_temp.split("\n")[:-1]
    # add each snp in the region to a list
    # as lists of
    snps = []
    for line in snps_split:
        snp = line.split('\t')[:8]
        snps.append(snp)
    return snps


def get_exons(gene_list):
    """ Take a list of transcript information in refgene format and return a
    list of exons in the region as [[e1_start, e1_end], [e2_start], [e2_end],
    ..]. The transcripts must belong to the same gene (i.e. have the same gene
    name).Merge overlapping exons.
    """
    # get start and end coordinates of exons in gene list
    starts = []
    ends = []
    gene_names = []
    gene_ids = []
    chrom_list = []
    for gene in gene_list:
        chrom_list.append(gene[2])
    chrom_set = list(set(chrom_list))
    if len(chrom_set) == 0:
        return {}
    chrom_set = [c for c in chrom_set if len(c) < 6]
    if len(chrom_set) > 1:
        print(("More than one chromosomes, ",
               chrom_set,
               ", has specified gene ",
               gene[12]))
        return {}
    chrom = chrom_set[0]
    for gene in gene_list:
        if gene[2] == chrom:
            starts.extend(list(map(int, gene[9].split(",")[:-1])))
            ends.extend(list(map(int, gene[10].split(",")[:-1])))
            gene_names.append(gene[12])
            gene_ids.append(gene[1])
            ori = gene[3]
    # pair exon starts and ends
    exons = []
    for i in range(len(starts)):
        exons.append([starts[i], ends[i]])
    # check for overlapping exons and merge if any
    overlapping = 1
    while overlapping:
        overlapping = 0
        for i in range(len(exons)):
            e = exons[i]
            for j in range(len(exons)):
                x = exons[j]
                if (i != j) and ((e[0] <= x[0] <= e[1])
                                 or (e[0] <= x[1] <= e[1])
                                 or (x[0] <= e[0] <= x[1])):
                    # merge exons and add to the exon list
                    exons.append([min(e[0], x[0]), max(e[1], x[1])])
                    # remove the exons e and x
                    exons.remove(e)
                    exons.remove(x)
                    # change overlapping to 1 so we can stop the outer for loop
                    overlapping = 1
                    # once an overlapping exon is found, break the for loop
                    break
            if overlapping:
                # if an overlapping exon is found, stop this for loop and
                # continue with the while loop with the updated exon list
                break
    # get the gene start and end coordinates
    if (len(starts) >= 1) and (len(ends) >= 1):
        start = min(starts)
        end = max(ends)
    else:
        print(("No exons found for ",  gene_list[0][1]))
        return {}
    # create an output dict
    out = {}
    out["chrom"] = chrom
    out["begin"] = start + 1
    out["end"] = end
    out["exons"] = [[e[0] + 1, e[1]] for e in sorted(exons, key=itemgetter(0))]
    out["names"] = gene_names
    out["ids"] = gene_ids
    out["orientation"] = ori
    return out


def get_gene_name(region, species):
    """ Return the gene(s) in a region. """
    gene_names = []
    try:
        genes = get_snps(region, get_file_locations()[species][
            "refgene_tabix"])
        for g in genes:
            gene_names.append(g[12])
    except KeyError:
        pass
    return gene_names


def get_gene(gene_name, refgene_file, chrom=None, alternative_chr=1):
    """ Return genomic coordinates of a gene extracted from the refseq genes file.
    Refgene fields are as follows:
    0:bin, 1:name, 2:chrom, 3:strand, 4:txStart, 5:txEnd, 6:cdsStart, 7:cdsEnd,
    8:exonCount, 9:exonStarts, 10:exonEnds, 11:score, 12:name2,
    13:cdsStartStat, 14:cdsEndStat, 15:exonFrames.
    Field 12 will be used for name search."""
    # all chromosomes must be included if chromosome of the gene is not
    # provided therefore, chrom cannot be None when alternative_chr is set to 0
    if not (chrom or alternative_chr):
        print(("Chromosome of the gene %s must be specified "
               "or all chromosomes must be searched."))
        print(("Specify a chromosome or set alternative chromosome to 1."
               % gene_name))
        return 1
    with open(refgene_file, 'r') as infile:
        coord = []
        for line in infile:
            if not line.startswith('#'):
                newline = line.strip().split('\t')
                if newline[12] == gene_name:
                    coord.append(newline)
    if len(coord) < 1:
        print(("No gene found with the name ", gene_name))
        return []
    alter = []
    if chrom:
        # add each gene to alter dict, in the corresponding chromosome key
        for c in coord:
            if c[2] == chrom:
                alter.append(c)
    # find genes on alternate chromosomes if requested
    elif alternative_chr:
        for c in coord:
            alter.append(c)
    return alter


def create_gene_fasta(gene_name_list, wdir, species="hs", flank=150,
                      multi_file=False):
    """ Get a list of genes, extract exonic sequence + flanking sequence.
    Create fasta files in corresponding directory for each gene if multi_file
    is True, create a single fasta file if False.
    """
    region_list = []
    for gene_name in gene_name_list:
        if gene_name.startswith("chr"):
            coord = get_coordinates(gene_name)
            query = make_region(coord[0], coord[1] - flank, coord[2] + flank)
        else:
            e = get_exons(
                get_gene(gene_name, get_file_locations()[species]["refgene"],
                         alternative_chr=1)
                )
            query = e["chrom"] + ":" + str(e["begin"] - flank) + "-" + str(
               e["end"] + flank)
        region_list.append(query)
    regions = get_fasta_list(region_list, species)
    fasta_dict = {}
    for i in range(len(region_list)):
        r = region_list[i]
        gene_name = gene_name_list[i]
        fasta_dict[gene_name] = regions[r]
    if multi_file:
        for gene_name in fasta_dict:
            save_dict = {gene_name: fasta_dict[gene_name]}
            filename = os.path.join(wdir, gene_name + ".fa")
            save_fasta_dict(save_dict, filename)
    else:
        save_fasta_dict(fasta_dict, os.path.join(wdir, "multi.fa"))


def get_region_exons(region, species):
    try:
        genes = get_snps(region, get_file_locations()[species][
            "refgene_tabix"])
    except KeyError:
        genes = []
    return get_exons(genes)


def get_cds(gene_name, species):
    gene_list = get_gene(gene_name,
                         get_file_locations()[species]["refgene"],
                         alternative_chr=1)
    if len(gene_list) > 1:
        print(("More than one refgene entry was found for the gene ",
               gene_name))
        print(("Exons from alternative transcripts will be merged "
               "and CDS will be generated from that."))
        print("This may lead to unreliable CDS sequence information.")
    if len(gene_list) == 0:
        return {}
    g = gene_list[0]
    cds = {"chrom": g[2],
           "orientation": g[3],
           "begin": int(g[6]) + 1,
           "end": int(g[7])}
    exons = get_exons(gene_list)["exons"]
    exons_nuc = []
    for i in range(len(exons)):
        e = exons[i]
        if not e[0] <= cds["begin"] <= e[1]:
            exons[i] == "remove"
        else:
            e[0] = cds["begin"]
            break
    exons = [i for i in exons if i != "remove"]
    for i in range(-1, -1 * len(exons), -1):
        e = exons[i]
        if not e[0] <= cds["end"] <= e[1]:
            exons[i] = "remove"
        else:
            e[1] = cds["end"]
            break
    exons = [i for i in exons if i != "remove"]
    sequences = []
    for e in exons:
        exons_nuc.extend(list(range(e[0], e[1] + 1)))
        sequences.append(fasta_to_sequence(
            get_fasta(cds["chrom"]
                      + ":" + str(e[0]) + "-"
                      + str(e[1]), species)))
    coord = {}
    if cds["orientation"] == "+":
        cds["sequence"] = "".join(sequences)
        for i in range(len(exons_nuc)):
            coord[i] = exons_nuc[i]
    else:
        cds["sequence"] = reverse_complement("".join(sequences))
        rev_exons = list(reversed(exons_nuc))
        for i in range(len(exons_nuc)):
            coord[i] = rev_exons[i]
    cds["coordinates"] = coord
    cds
    return cds


def make_boulder(fasta, primer3_input_DIR, exclude_list=[],
                 output_file_name="", sequence_targets=[]):
    """ Create a boulder record file in primer3_input_DIR from a given fasta
    STRING. SEQUENCE_ID is the fasta header, usually the genomic region
    (chrX:m-n) exclude_list is [coordinate,length] of any regions primers
    cannot overlap.
    """
    # parse fasta string, get header and remove remaining nextlines.
    fasta_list = fasta.split("\n")
    fasta_head = fasta_list[0][1:]
    seq_template = "".join(fasta_list[1:])
    # convert exclude list to strings
    exclude_string_list = []
    exclude_region = ""
    for i in exclude_list:
        exclude_string_list.append(str(i[0])+","+str(i[1]))
        exclude_region = " ".join(exclude_string_list)
    # create the boulder record
    if len(sequence_targets) == 0:
        sequence_target_string = ""
    else:
        sequence_target_string = " ".join([",".join(map(str, s))
                                           for s in sequence_targets])
    boulder = ("SEQUENCE_ID=" + fasta_head + "\n" +
               "SEQUENCE_TEMPLATE=" + seq_template + "\n" +
               "SEQUENCE_TARGET=" + sequence_target_string + "\n" +
               "SEQUENCE_EXCLUDED_REGION=" + exclude_region + "\n" + "=")
    if output_file_name == "":
        outname = fasta_head
    else:
        outname = output_file_name
    with open(os.path.join(primer3_input_DIR, outname), 'w') as outfile:
        outfile.write(boulder)
    return boulder


def make_primers_worker(l):
    """
    Worker function to make_primers_multi.

    A worker function to make primers for multiple regions using separate
    processors. Read boulder record in given input directory and creates primer
    output files in output directory
    """
    # function arguments should be given as a list due to single
    # iterable limitation of map_async function of multiprocessor.Pool
    # input boulder record name
    input_file = l[0]
    # primer settings used
    settings = l[1]
    # output file name
    output_file = l[2]
    # locations of input/output dirs
    primer3_input_DIR = l[3]
    primer3_output_DIR = l[4]
    primer3_settings_DIR = l[5]
    subregion_name = l[6]
    paralog_name = l[7]
    primer_type = l[8]
    input_file = os.path.join(primer3_input_DIR, input_file)
    output_file = os.path.join(primer3_output_DIR, output_file)
    settings = os.path.join(primer3_settings_DIR, settings)
    # call primer3 program using the input and settings file
    res = subprocess.run(["primer3_core",
                          "-p3_settings_file=" + settings, input_file],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if res.returncode != 0:
        print(("Primer design for the gene {} subregion {} {} arm failed "
               "with error {}").format(paralog_name, subregion_name,
                                       primer_type, res.stderr))
        return
    else:
        primer3_output = res.stdout
        # write boulder record to file.
        with open(output_file, 'w') as outfile:
            outfile.write(primer3_output.decode("UTF-8"))
        return


def make_primers_multi(ext_list, lig_list, pro):
    """Design primers in parallel using the make_primers_worker function."""
    # create a pool of twice the number of targets (for extension and ligation)
    # p = Pool(2*pro)
    p = Pool(pro)
    # make extension primers using extension arm primer settings
    p.map_async(make_primers_worker, ext_list)
    # make ligation primers using ligation arm primer settings
    p.map_async(make_primers_worker, lig_list)
    # close pool
    p.close()
    # wait for processes to finish
    p.join()
    return


def primer_parser3(input_file, primer3_output_DIR, bowtie2_input_DIR,
                   parse_out, fasta=1, outp=1):
    """
    Parse a primer3 output file and generate a primer fasta file.

    The fasta file for the primers that only contains primer names and
    sequences will be placed in the bowtie input directory  to be
    used as bowtie2 input.
    Return a dictionary {sequence_information:{}, primer_information{}}
    first dict has tag:value pairs for input sequence while second dict
    has as many dicts as the primer number returned with primer name keys
    and dicts as values {"SEQUENCE": "AGC..", "TM":"58"...}. Also write
    this dictionary to a json file in primer3_output_DIR.
    """
    primer_dic = {}
    # all target sequence related information will be placed in
    # sequence_information dictionary.
    primer_dic["sequence_information"] = {}
    # primer information will be kept in primer_information dicts.
    primer_dic["primer_information"] = {}
    # load the whole input file into a list.
    infile = open(primer3_output_DIR + input_file, 'r')
    lines = []
    for line in infile:
        # if a line starts with "=" that line is a record separator
        if not line.startswith("="):
            # boulder record tag-value pairs separated by "="
            inline = line.strip('\n').split('=')
            lines.append(inline)
    infile.close()
    # find sequence related information and add it to appropriate dic.
    for pair in lines:
        tag = pair[0]
        value = pair[1]
        if tag.startswith("SEQUENCE"):
            if tag == "SEQUENCE_ID":
                new_value = value.split(",")[-1].replace("CHR", "chr")
                primer_dic["sequence_information"][tag] = new_value
            else:
                primer_dic["sequence_information"][tag] = value
    # find how many left primers returned and create empty dictionary
    # for each primer in primer_information dict.
    for pair in lines:
        tag = pair[0]
        value = pair[1]
        if tag == "PRIMER_LEFT_NUM_RETURNED":
            # Add this to sequence information dic because it is sequence
            # specific information
            primer_dic["sequence_information"][
                "SEQUENCE_LEFT_NUM_RETURNED"] = value
            # create empty dictionaries with primer name keys
            for i in range(int(value)):
                primer_key = "PRIMER_LEFT_" + str(i)
                primer_dic["primer_information"][primer_key] = {}
    # do the same for right primers found
    for pair in lines:
        tag = pair[0]
        value = pair[1]
        if tag == "PRIMER_RIGHT_NUM_RETURNED":
            primer_dic["sequence_information"][
                "SEQUENCE_RIGHT_NUM_RETURNED"] = value
            for i in range(int(value)):
                primer_key = "PRIMER_RIGHT_" + str(i)
                primer_dic["primer_information"][primer_key] = {}
    # get sequence coordinate information to determine genomic coordinates of
    # primers because primer information is relative to template sequence
    sequence_coordinates = get_coordinates(primer_dic[
        "sequence_information"]["SEQUENCE_ID"])
    seq_chr = sequence_coordinates[0]
    seq_start = int(sequence_coordinates[1])
    # get primer information from input file and add to primer dictionary
    for pair in lines:
        tag = pair[0]
        value = pair[1]
        if ((tag.startswith("PRIMER_LEFT_")
             or tag.startswith("PRIMER_RIGHT_"))
            and (tag != "PRIMER_LEFT_NUM_RETURNED")
                and (tag != "PRIMER_RIGHT_NUM_RETURNED")):
            attributes = tag.split('_')
            # primer coordinates tag does not include an attribute value
            # it is only primer name = coordinates, so:
            if len(attributes) > 3:
                # then this attribute is not coordinates and should have an
                # attribute value such as TM or HAIRPIN etc.
                primer_name = '_'.join(attributes[0:3])
                attribute_value = '_'.join(attributes[3:])
                primer_dic["primer_information"][primer_name][
                    attribute_value] = value
            else:
                # then this attribute is coordinates and has no attribute value
                # give it an attribute valute "COORDINATES"
                primer_name = '_'.join(attributes[0:3])
                primer_dic["primer_information"][primer_name][
                    'COORDINATES'] = value
                # the coordinates are relative to sequence template
                # find the genomic coordinates
                coordinate_values = value.split(",")
                if tag.startswith("PRIMER_LEFT"):
                    # sequence start is added to primer start to get genomic
                    # primer start
                    genomic_start = seq_start + int(coordinate_values[0])
                    # primer len is added "to genomic start because it is a
                    # left primer
                    genomic_end = genomic_start + int(coordinate_values[1]) - 1
                    primer_dic["primer_information"][primer_name][
                        'GENOMIC_START'] = genomic_start
                    primer_dic["primer_information"][primer_name][
                        'GENOMIC_END'] = genomic_end
                    primer_dic["primer_information"][primer_name][
                        'CHR'] = seq_chr
                    primer_dic["primer_information"][primer_name][
                        'ORI'] = "forward"
                else:
                    # sequence start is added to primer start to get genomic
                    # primer start
                    genomic_start = seq_start + int(coordinate_values[0])
                    # primer len is subtracted from genomic start because it is
                    # a right primer
                    genomic_end = genomic_start - int(coordinate_values[1]) + 1
                    primer_dic["primer_information"][primer_name][
                        'GENOMIC_START'] = genomic_start
                    primer_dic["primer_information"][primer_name][
                        'GENOMIC_END'] = genomic_end
                    primer_dic["primer_information"][primer_name][
                        'CHR'] = seq_chr
                    primer_dic["primer_information"][primer_name][
                        'ORI'] = "reverse"
            # add NAME as a key to primer information dictionary
            primer_dic["primer_information"][primer_name]['NAME'] = primer_name
    # if some primers were eliminated from initial primer3 output, remove from
    # dictionary
    for primer in list(primer_dic["primer_information"].keys()):
        if primer_dic["primer_information"][primer] == {}:
            primer_dic["primer_information"].pop(primer)
    # dump the dictionary to json file in primer3_output_DIR if outp parameter
    # is true
    if outp:
        dict_file = open(os.path.join(primer3_output_DIR, parse_out), 'w')
        json.dump(primer_dic, dict_file, indent=1)
        dict_file.close()
    # generate a simple fasta file with primer names
    if fasta:
        outfile = open(bowtie2_input_DIR+parse_out, 'w')
        for primer in primer_dic["primer_information"]:
            # primer name is fasta header and sequence is fasta sequence
            fasta_head = primer
            fasta_line = primer_dic["primer_information"][primer]["SEQUENCE"]
            outfile.write(">" + fasta_head + "\n" + fasta_line + "\n")
        outfile.close()
    return primer_dic


def paralog_primers(primer_dict, copies, coordinate_converter, settings,
                    primer3_output_DIR, outname, species, outp=0):
    """
    Process primers generated for paralogs.

    Take a primer dictionary file and add genomic start and end coordinates
    of all its paralogs.
    """
    # uncomment for using json object instead of dic
    # load the primers dictionary from file
    # with open(primer_file, "r") as infile:
    #     primer_dic = json.load(infile)
    # primer dict consists of 2 parts, sequence_information dict
    # and primer information dict. We wont'change the sequence_info part
    primers = primer_dict["primer_information"]
    primer_keys = set()
    for primer in list(primers.keys()):
        p_name = primer
        p_dic = primers[primer]
        p_coord = coordinate_converter
        p_copies = copies
        chroms = p_coord["C0"]["chromosomes"]
        start = p_dic["GENOMIC_START"]
        end = p_dic["GENOMIC_END"]
        ref_coord = p_dic["COORDINATES"]
        primer_ori = p_dic["ORI"]
        p_dic["PARALOG_COORDINATES"] = {}
        primer_seq = p_dic["SEQUENCE"]
        # add reference copy as paralog
        p_dic["PARALOG_COORDINATES"]["C0"] = {"SEQUENCE": primer_seq,
                                              "ORI": primer_ori,
                                              "CHR": chroms["C0"],
                                              "NAME": p_name,
                                              "GENOMIC_START": start,
                                              "GENOMIC_END": end,
                                              "COORDINATES": ref_coord}
        for c in p_copies:
            if c != "C0":
                # check if both ends of the primer has aligned with reference
                try:
                    para_start = p_coord["C0"][c][start]
                    para_end = p_coord["C0"][c][end]
                except KeyError:
                    # do not add that copy if it is not aligned
                    continue
                para_primer_ori = para_start < para_end
                if para_primer_ori:
                    para_primer_key = (chroms[c] + ":" + str(para_start) + "-"
                                       + str(para_end))
                    p_dic["PARALOG_COORDINATES"][c] = {
                        "ORI": "forward", "CHR": chroms[c], "NAME": p_name,
                        "GENOMIC_START": para_start, "GENOMIC_END": para_end,
                        "COORDINATES": ref_coord, "KEY": para_primer_key}
                    primer_keys.add(para_primer_key)
                else:
                    para_primer_key = chroms[c] + ":" + str(
                        para_end) + "-" + str(para_start)
                    p_dic["PARALOG_COORDINATES"][c] = {
                        "ORI": "reverse", "CHR": chroms[c], "NAME": p_name,
                        "GENOMIC_START": para_start, "GENOMIC_END": para_end,
                        "COORDINATES": ref_coord, "KEY": para_primer_key}
                    primer_keys.add(para_primer_key)
    if len(primer_keys) > 0:
        primer_sequences = get_fasta_list(primer_keys, species)
        for p in primers:
            para = primers[p]["PARALOG_COORDINATES"]
            for c in para:
                if c != "C0":
                    copy_dict = para[c]
                    p_ori = copy_dict["ORI"]
                    p_key = copy_dict["KEY"]
                    p_seq = primer_sequences[p_key]
                    if p_ori == "reverse":
                        p_seq = reverse_complement(p_seq)
                    copy_dict["SEQUENCE"] = primer_sequences[p_key]
    if outp:
        with open(os.path.join(primer3_output_DIR, outname), "w") as outf:
            json.dump(primer_dict, outf, indent=1)
    return primer_dict


def bowtie2_run(fasta_file, output_file, bowtie2_input_DIR,
                bowtie2_output_DIR, species, process_num=4,
                seed_MM=1, mode="-a", seed_len=18, gbar=1, local=0):
    """Align primers from a fasta file to specified species genome."""
    file_locations = get_file_locations()
    # check if entered species is supported
    genome = file_locations[species]["bowtie2_genome"]
    # determine what type of alignment is wanted
    # local or end-to-end
    if local:
        check_local = "--local"
    else:
        check_local = "--end-to-end"
    res = subprocess.Popen(["bowtie2", "-p", str(process_num),  "-D", "20",
                            "-R", "3", "-N", str(seed_MM), "-L",
                            str(seed_len), "-i", "S,1,0.5", "--gbar",
                            str(gbar), mode, check_local, "-x", genome, "-f",
                            os.path.join(bowtie2_input_DIR, fasta_file), "-S",
                            os.path.join(bowtie2_output_DIR, output_file)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    log_file = os.path.join(
        bowtie2_output_DIR, "log_" + species + "_" + id_generator(6))
    with open(log_file, "wb") as outfile:
        outfile.write(res.communicate()[1])

    return 0


def bowtie(fasta_file, output_file, bowtie2_input_DIR, bowtie2_output_DIR,
           options, species, process_num=4, mode="-a", local=0, fastq=0):
    """Align a fasta or fastq file to a genome using bowtie2."""
    file_locations = get_file_locations()
    # check if entered species is supported
    genome = file_locations[species]["bowtie2_genome"]
    # determine what type of alignment is wanted
    # local or end-to-end
    if local:
        check_local = "--local"
    else:
        check_local = "--end-to-end"
    com = ["bowtie2", "-p " + str(process_num)]
    com.extend(options)
    com.append(mode)
    com.append(check_local)
    com.append("-x " + genome)
    if fastq:
        com.append("-q " + os.path.join(bowtie2_input_DIR, fasta_file))
    else:
        com.append("-f " + os.path.join(bowtie2_input_DIR, fasta_file))
    com.append("-S " + os.path.join(bowtie2_output_DIR, output_file))
    subprocess.check_output(com)
    return 0


def bwa(fastq_file, output_file, output_type, input_dir,
        output_dir, options, species, base_name="None"):
    """
    Align a fastq file to species genome using bwa.

    Options should be a list that starts with the command (e.g. mem, aln etc).
    Additional options should be appended as strings of "option value",
    for example, "-t 30" to use 30 threads. Output type can be sam or bam.
    Recommended options ["-t30", "-L500", "-T100"]. Here L500 penalizes
    clipping severely so the alignment becomes end-to-end and T100 stops
    reporting secondary alignments, assuming their score is below 100.
    """
    genome_file = get_file_locations()[species]["bwa_genome"]
    read_group = ("@RG\\tID:" + base_name + "\\tSM:" + base_name + "\\tLB:"
                  + base_name + "\\tPL:ILLUMINA")
    options = copy.deepcopy(options)
    options.append("-R" + read_group)
    if output_type == "sam":
        com = ["bwa"]
        com.extend(options)
        com.append(genome_file)
        com.append(os.path.join(input_dir, fastq_file))
        with open(os.path.join(output_dir, output_file), "w") as outfile:
            subprocess.check_call(com, stdout=outfile)
    else:
        com = ["bwa"]
        com.extend(options)
        com.append(genome_file)
        com.append(os.path.join(input_dir, fastq_file))
        sam = subprocess.Popen(com, stdout=subprocess.PIPE)
        bam_com = ["samtools", "view", "-b"]
        bam = subprocess.Popen(bam_com, stdin=sam.stdout,
                               stdout=subprocess.PIPE)
        bam_file = os.path.join(output_dir, output_file)
        sort_com = ["samtools", "sort", "-T", "/tmp/", "-o", bam_file]
        subprocess.run(sort_com, stdin=bam.stdout)
        subprocess.run(["samtools", "index", bam_file], check=True,
                       stderr=subprocess.PIPE)


def bwa_multi(fastq_files, output_type, fastq_dir, bam_dir, options, species,
              processor_number, parallel_processes):
    """Align fastq files to species genome using bwa in parallel."""
    if len(fastq_files) == 0:
        fastq_files = [f.name for f in os.scandir(fastq_dir)]
    if output_type == "sam":
        extension = ".sam"
    elif output_type == "bam":
        extension = ".srt.bam"
    else:
        print(("Output type must be bam or sam, {} was given").format(
            output_type))
        return
    if not os.path.exists(bam_dir):
        os.makedirs(bam_dir)
    if parallel_processes == 1:
        for f in fastq_files:
            # get base file name
            base_name = f.split(".")[0]
            bam_name = base_name + extension
            options.extend("-t" + str(processor_number))
            bwa(f, bam_name, output_type, fastq_dir, bam_dir, options, species,
                base_name)
    else:
        processor_per_process = processor_number // parallel_processes
        p = NoDaemonProcessPool(parallel_processes)
        options = options + ["-t " + str(processor_per_process)]
        results = []
        errors = []
        for f in fastq_files:
            base_name = f.split(".")[0]
            bam_name = base_name + extension
            p.apply_async(bwa, (f, bam_name, output_type, fastq_dir, bam_dir,
                                options, species, base_name),
                          callback=results.append,
                          error_callback=errors.append)
        p.close()
        p.join()
        if len(errors) > 0:
            for e in errors:
                print("Error in bwa_multi function", e.stderr)


def parse_cigar(cigar):
    """
    Parse a CIGAR string.

    CIGAR string is made up of numbers followed
    by key letters that represent a sequence alignment; return a dictionary
    with alignment keys and number of bases with that alignment key as values.
    Below is some more information about cigar strings.

    2S20M1I2M5D,for, example would mean that the 2 bases are "S"oft clipped
    from 5' end of the sequence(read) aligned and it is not part of the
    alignment; following that 2 bases, 20 bases of the read aligns or "M"atches
    to the reference sequence, match here does not mean the bases are
    identical, just that there is 1 base of reference for each base of the read
    and there are enough similarity between the two sequences that they
    aligned. 1 base following the 20M is an insertion, that is, it exists in
    the read but not in the reference; 5 bases at the end are "D"eletions,
    they are in the reference but not in the read.
    """
    cig = {}
    values = []
    for c in cigar:
        try:
            values.append(str(int(c)))
        except ValueError:
            if c in list(cig.keys()):
                cig[c] += int("".join(values))
            else:
                cig[c] = int("".join(values))
            values = []
    return cig


def get_cigar_length(cigar):
    """Get the length of the reference sequence from CIGAR string."""
    try:
        # parse cigar string and find out how many insertions are in the
        # alignment
        insertions = parse_cigar(cigar)["I"]
    except KeyError:
        # the key "I" will not be present in the cigar string if there is no
        # insertion
        insertions = 0
    # all the values in the cigar dictionary represent a base in the reference
    # seq,
    # except the insertions, so they should be subtracted
    return sum(parse_cigar(cigar).values()) - insertions


def parse_bowtie(primer_dict, bt_file, primer_out, primer3_output_DIR,
                 bowtie2_output_DIR, species, settings, outp=1):
    """
    Take a bowtie output (sam) file and filter top N hits per primer.

    When a primer has more than "upper_hit_limit" bowtie hits,
    remove that primer.
    Add the bowtie hit information, including hit sequence to
    the primers dictionary.
    """
    # extract how many bowtie hits should be added
    # to the primer information for further TM analysis
    N = int(settings["hit_limit"])
    # how many total bowtie hits gets a primer fired
    M = int(settings["upper_hit_limit"])
    # read in bowtie file
    infile = open(os.path.join(bowtie2_output_DIR, bt_file), 'r')
    primers = copy.deepcopy(primer_dict)
    # create a temp dic to count hits/primer
    counter_dic = {}
    # create a bowtie key that will be used when adding
    # bowtie information to primers
    bowtie_key = "bowtie_information_" + species
    # all bowtie hits that will be used further for TM analysis
    # will need to have sequence information with them
    # region keys for hits (in chrx:begin-end format) will be
    # kept in a list for mass fasta extraction later.
    keys = set()
    #
    # read bowtie hits
    for line in infile:
        try:
            if not line.startswith("@"):
                record = line.strip('\n').split('\t')
                primer_name = record[0]
                # increment hit counter for primer
                try:
                    counter_dic[primer_name] += 1
                except KeyError:
                    counter_dic[primer_name] = 1
                # check how many hits have been analyzed for this primer
                # if upper hit limit has been reached, mark primer for removal
                if counter_dic[primer_name] >= M:
                    primers['primer_information'][primer_name]["remove"] = True
                    continue
                # move on to the next hit if primer hit limit has been reached.
                # no further hits will be added for those primers
                if counter_dic[primer_name] >= N:
                    continue
                flag = record[1]
                # a flag value of 4 means there was no hit, so pass those lines
                if flag == "4":
                    continue
                # chromosome of the bowtie hit
                chrom = record[2]
                # genomic position of bowtie hit
                pos = int(record[3])
                # get cigar string of alignment
                cigar = record[5]
                # extract which strand is the bowtie hit on
                # true if forward
                strand = ((int(record[1]) % 256) == 0)
                # get hit coordinates
                hit_start = pos
                # bowtie gives us the start position of the hit
                # end position is calculated using the cigar string
                # of the hit
                hit_end = pos + get_cigar_length(cigar) - 1
                # create region keys required for sequence retrieval
                # we want 3 nt extra on the 5' of the primer
                # because when alternative primers for paralogs
                # are considered we check +/- 3 nt from 5' end
                # to balance TM.
                if strand:
                    # Primer's 5' is the hit start when the hit is on forward
                    # strand so the nucleotides are added at start position
                    bt_start = hit_start
                    bt_end = hit_end
                    hit_str = "forward"
                    hit_region_key = (chrom + ":" + str(hit_start)
                                      + "-" + str(hit_end))
                else:
                    bt_start = hit_end
                    bt_end = hit_start
                    hit_str = "reverse"
                    hit_region_key = (chrom + ":" + str(hit_start)
                                      + "-" + str(hit_end))
                # add region key to keys list for fasta retrieval later
                keys.add(hit_region_key)
                # add all hit information to primer dictionary
                try:
                    primers["primer_information"][primer_name][bowtie_key][
                        str(counter_dic[primer_name])
                    ] = {"chrom": chrom, "begin": bt_start, "end": bt_end,
                         "key": hit_region_key, "strand": hit_str}
                except KeyError:
                    primers["primer_information"][primer_name][bowtie_key] = {
                         str(counter_dic[primer_name]): {"chrom": chrom,
                                                         "begin": bt_start,
                                                         "end": bt_end,
                                                         "key": hit_region_key,
                                                         "strand": hit_str}
                    }
        except KeyError:
            # in earlier versions of this function the primers with
            # excessive hits were removed during iteration and that lead
            # to keyerrors. Now there should be no key error.
            continue
    # get the fasta sequences of all hits
    sequence_dic = get_fasta_list(keys, species)
    # remove primers with too many hits and add bowtie information for others.
    for p in list(primers["primer_information"].keys()):
        try:
            if primers["primer_information"][p]["remove"]:
                primers["primer_information"].pop(p)
                continue
        except KeyError:
            pass
        # add hit sequences to primer dictionary
        # forward strand hits are added directly
        # reverse strand hits are reversed-complemented
        # so the hit is always in the primer orientation and
        # and similar in sequence"
        try:
            for h in primers["primer_information"][p][bowtie_key]:
                if (primers["primer_information"][p]
                        [bowtie_key][h]["strand"] == "forward"):
                    primers["primer_information"][p][bowtie_key][h][
                        "sequence"
                    ] = sequence_dic[primers["primer_information"][p][
                        bowtie_key][h]["key"]
                    ]
                else:
                    primers["primer_information"][p][bowtie_key][h][
                        "sequence"
                    ] = reverse_complement(
                        sequence_dic[primers["primer_information"]
                                     [p][bowtie_key][h]["key"]]
                    )
        except KeyError:
            # if there is no bowtie hit for this primer (happens for host
            # species):
            primers["primer_information"][p][bowtie_key] = {}
    # save the updated primers file
    if outp:
        with open(os.path.join(
                primer3_output_DIR, primer_out), 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers


def process_bowtie(primers, primer_out, primer3_output_DIR,
                   bowtie2_output_DIR, species, settings, host=False, outp=1):
    """
    Process a primer dict with bowtie information added.

    Look at bowtie hits for each primer, determine if they
    are on intended targets or nonspecific. In cases of paralogus
    regions, check all paralogs and determine if the primer
    will bind to any paralog. Create alternative primers if necessary
    and allowed. Get melting temperatures of all hits and add
    all these information to the primer dictionary.
    """
    # get Na, Mg and oligo concentrations these are specified in M but primer3
    # uses mM for ions and nM for oligos, so those will be adjusted.
    Na = float(settings["Na"]) * 1000
    Mg = float(settings["Mg"]) * 1000
    conc = float(settings["oligo_conc"]) * pow(10, 9)
    # are alternative mip arms allowed/desired
    alt_arm = int(settings["alternative_arms"])
    bowtie_key = "bowtie_information_" + species
    alt_keys = set([])
    # get reference chromosome lengths
    genome_file = get_file_locations()[species]["fasta_genome"]
    reference_lengths = {}
    genome_sam = pysam.FastaFile(genome_file)
    for r in genome_sam.references:
        reference_lengths[r] = genome_sam.get_reference_length(r)
    # read bowtie hits
    for primer_name in primers['primer_information']:
        try:
            primer_seq = primers['primer_information'][primer_name]["SEQUENCE"]
            if not host:
                para = (primers['primer_information'][primer_name]
                        ["PARALOG_COORDINATES"])
                if ("BOWTIE_BINDS" not in
                        primers['primer_information'][primer_name]):
                    primers[
                        'primer_information'][primer_name]["BOWTIE_BINDS"] = []
                if ("ALT_BINDS" not in
                        primers['primer_information'][primer_name]):
                    primers[
                        'primer_information'][primer_name]["ALT_BINDS"] = []
            for bt_hit_name in list(primers['primer_information']
                                    [primer_name][bowtie_key].keys()):
                bt_hit = (primers['primer_information'][primer_name]
                          [bowtie_key][bt_hit_name])
                bt_chrom = bt_hit["chrom"]
                bt_begin = bt_hit["begin"]
                bt_end = bt_hit["end"]
                bt_ori = bt_hit["strand"]
                bt_seq = bt_hit["sequence"]
                if host:
                    bt_hit["TM"] = calcHeterodimerTm(
                        primer_seq,
                        reverse_complement(bt_seq),
                        mv_conc=Na,
                        dv_conc=Mg,
                        dntp_conc=0,
                        dna_conc=conc
                    )
                    continue
                intended = 0
                # para is a dict like:
                # {C0:{"CHR": "chr4", "GENOMIC_START" ..}, C1:{..
                # for non-CNV regions, bowtie mapping should be exactly the
                # same as genomic coordinates, so even if there is 1 bp
                # difference, we'll count this as off target. For CNV regions,
                # a more generous 20 bp padding will be allowed to account for
                # differences in our mapping and bowtie mapping. Bowtie mapping
                # will be accepted as the accurate mapping and paralog
                # coordinates will be changed accordingly.
                map_padding = 1
                if len(para) > 1:
                    map_padding = 20
                for k in para:
                    para_ori = para[k]["ORI"]
                    para_chr = para[k]["CHR"]
                    para_begin = para[k]["GENOMIC_START"]
                    para_end = para[k]["GENOMIC_END"]
                    if ((para_ori == bt_ori) and (para_chr == bt_chrom)
                            and (abs(para_begin - bt_begin) < map_padding)
                            and (abs(para_end - bt_end) < map_padding)):
                        intended = 1
                        # Get bowtie determined coordinates and sequences
                        # for the paralog copy. These will have priority
                        # over GENOMIC_ values calculated internally.
                        para[k]["BOWTIE_END"] = bt_end
                        para[k]["BOWTIE_START"] = bt_begin
                        para[k]["BOWTIE_SEQUENCE"] = bt_seq
                    if intended:
                        # if the paralog sequence is the same as the reference
                        # this primer should bind to the paralog copy as well.
                        if bt_seq.upper() == primer_seq.upper():
                            para[k]["BOWTIE_BOUND"] = True
                            primers['primer_information'][
                                primer_name]["BOWTIE_BINDS"].append(k)
                        else:
                            # if the sequences are not exactly the same
                            # we'll assume the primer does not bind to the
                            # paralog and attempt to generate an alternative
                            # primer for this paralog.
                            para[k]["BOWTIE_BOUND"] = False
                            # Do this only if alternative MIP arms are allowed
                            # specified by alt_arm setting.
                            if alt_arm:
                                # get chromosome length to avoid setting
                                # alt arms beyon chromosome ends
                                para_chr_length = reference_lengths[para_chr]
                                al = {}
                                al["ref"] = {"ALT_SEQUENCE": primer_seq}
                                al["ref"]["ALT_TM"] = calcHeterodimerTm(
                                    primer_seq,
                                    reverse_complement(primer_seq),
                                    mv_conc=Na,
                                    dv_conc=Mg,
                                    dntp_conc=0,
                                    dna_conc=conc
                                )
                                for j in range(-3, 4):
                                    if j == 0:
                                        continue
                                    alt_start = bt_begin + j
                                    alt_end = bt_end
                                    if ((alt_start < 0) or (alt_end < 0)
                                            or (alt_start > para_chr_length)
                                            or (alt_end > para_chr_length)):
                                        continue
                                    if para_ori == "forward":
                                        alt_primer_key = create_region(
                                            bt_chrom,
                                            alt_start,
                                            alt_end
                                        )
                                    else:
                                        alt_primer_key = create_region(
                                            bt_chrom,
                                            alt_end,
                                            alt_start
                                        )
                                    al[j] = {}
                                    al[j]["ALT_START"] = alt_start
                                    al[j]["ALT_END"] = alt_end
                                    al[j]["ALT_ORI"] = para_ori
                                    al[j]["ALT_KEY"] = alt_primer_key
                                    alt_keys.add(alt_primer_key)
                                para[k]["ALTERNATIVES"] = al
                            else:
                                para[k]["ALTERNATIVES"] = {}
                                para[k]["ALT_TM"] = 0
                                para[k]["ALT_TM_DIFF"] = 100
                                para[k]["ALT_BOUND"] = False
                        # remove bowtie hit for intended target
                        primers['primer_information'][
                            primer_name][bowtie_key].pop(bt_hit_name)
                        break
                # add TM value for unindended target
                if not intended:
                    bt_hit["TM"] = calcHeterodimerTm(
                        primer_seq,
                        reverse_complement(bt_seq),
                        mv_conc=Na,
                        dv_conc=Mg,
                        dntp_conc=0,
                        dna_conc=conc
                    )
            # Design alternative primers (if allowed) for paralogs
            # when there is no bowtie hit for that paralog.
            if not host:
                for k in para:
                    try:
                        para[k]["BOWTIE_END"]
                    except KeyError:
                        para_ori = para[k]["ORI"]
                        para_chr = para[k]["CHR"]
                        para_begin = para[k]["GENOMIC_START"]
                        para_end = para[k]["GENOMIC_END"]
                        para[k]["BOWTIE_BOUND"] = False
                        if alt_arm:
                            # get chromosome length to avoid setting
                            # alt arms beyon chromosome ends
                            para_chr_length = reference_lengths[para_chr]
                            al = {}
                            al["ref"] = {"ALT_SEQUENCE": primer_seq}
                            al["ref"]["ALT_TM"] = calcHeterodimerTm(
                                primer_seq,
                                reverse_complement(primer_seq),
                                mv_conc=Na,
                                dv_conc=Mg,
                                dntp_conc=0,
                                dna_conc=conc
                            )
                            for j in range(-3, 4):
                                if j == 0:
                                    continue
                                alt_start = para_begin + j
                                alt_end = para_end
                                if ((alt_start < 0) or (alt_end < 0)
                                        or (alt_start > para_chr_length)
                                        or (alt_end > para_chr_length)):
                                    continue
                                if para_ori == "forward":
                                    alt_primer_key = create_region(
                                        para_chr,
                                        alt_start,
                                        alt_end
                                    )
                                else:
                                    alt_primer_key = create_region(
                                        para_chr,
                                        alt_end,
                                        alt_start
                                    )
                                al[j] = {}
                                al[j]["ALT_START"] = alt_start
                                al[j]["ALT_END"] = alt_end
                                al[j]["ALT_ORI"] = para_ori
                                al[j]["ALT_KEY"] = alt_primer_key
                                alt_keys.add(alt_primer_key)
                            para[k]["ALTERNATIVES"] = al
                        else:
                            para[k]["ALTERNATIVES"] = {}
                            para[k]["ALT_TM"] = 0
                            para[k]["ALT_TM_DIFF"] = 100
                            para[k]["ALT_BOUND"] = False

        except KeyError:
            continue
    if len(alt_keys) > 0:
        alt_sequences = get_fasta_list(alt_keys, species)
        for primer_name in primers['primer_information']:
            para = (primers['primer_information'][primer_name]
                    ["PARALOG_COORDINATES"])
            for k in para:
                try:
                    alt_candidates = para[k]["ALTERNATIVES"]
                except KeyError:
                    continue
                for c in list(alt_candidates.keys()):
                    try:
                        alt_candidates[c]["ALT_TM"]
                    except KeyError:
                        alt_ori = alt_candidates[c]["ALT_ORI"]
                        alt_key = alt_candidates[c]["ALT_KEY"]
                        alt_seq = alt_sequences[alt_key]
                        if alt_ori == "reverse":
                            alt_seq = reverse_complement(alt_seq)
                        if alt_seq != "":
                            alt_tm = calcHeterodimerTm(
                                alt_seq,
                                reverse_complement(alt_seq),
                                mv_conc=Na,
                                dv_conc=Mg,
                                dntp_conc=0,
                                dna_conc=conc
                            )
                            alt_candidates[c]["ALT_TM"] = alt_tm
                            alt_candidates[c]["ALT_SEQUENCE"] = alt_seq
                        else:
                            alt_candidates.pop(c)
    if outp:
        with open(os.path.join(
                primer3_output_DIR, primer_out), 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers


def filter_bowtie(primers, output_file, primer3_output_DIR, species, TM=46,
                  hit_threshold=0, lower_tm=46, lower_hit_threshold=3, outp=1):
    """
    Check TMs of bowtie hits of given primers, on a given genome.

    Filter the primers with too many nonspecific hits.
    """
    for primer in list(primers["primer_information"].keys()):
        # create a hit count parameter for hits with significant tm
        # there are two parameters specified in the rinfo file
        # high temp limit and low temp limit. The idea is to allow
        # a very small (if any) number of nonspecific targets with high TM
        # values but allow some low TM off targets.
        hc = 0
        lhc = 0
        # check if bowtie information exists in dic
        try:
            bt_key = "bowtie_information_" + species
            bowtie = primers["primer_information"][primer][bt_key]
            for h in bowtie:
                hit = bowtie[h]
                try:
                    # if TM information is included in bowtie, compare with
                    # high and low TM, increment hc, lc if necessary and
                    # discard primers passing specified off target tresholds.
                    if float(hit["TM"]) >= TM:
                        hc += 1
                        if hc > hit_threshold:
                            primers["primer_information"].pop(primer)
                            break
                    elif float(hit["TM"]) >= lower_tm:
                        lhc += 1
                        if lhc > lower_hit_threshold:
                            primers["primer_information"].pop(primer)
                            break
                except KeyError:
                    continue
            # remove bowtie information once we use it.
            primers["primer_information"][primer].pop(bt_key)
        except KeyError:
            continue
    if outp:
        # write dictionary to file in primer3_output_DIR
        outfile = open(os.path.join(primer3_output_DIR, output_file), 'w')
        json.dump(primers, outfile, indent=1)
        outfile.close()
    return primers


def alternative(primer_dic, output_file,
                primer3_output_DIR, tm_diff, outp=1):
    """
    Pick the best alternative arm for primers that do not bind all paralogs.

    This is done by picking the alternative primer with melting temperature
    that is closest to the original primer.
    """
    primers = primer_dic["primer_information"]
    try:
        for primer_name in primers:
            primer = primers[primer_name]
            para = primer["PARALOG_COORDINATES"]
            for c in para:
                try:
                    alts = para[c]["ALTERNATIVES"]
                    # get the original primer TM
                    ref_tm = alts["ref"].pop("ALT_TM")
                    alts.pop("ref")
                    # sort alt primers by their TM difference from the ref
                    sorted_alts = sorted(
                        alts, key=lambda a: abs(alts[a]["ALT_TM"] - ref_tm)
                    )
                    # use the primer only if the TM difference is within
                    # specified limit.
                    if abs(alts[sorted_alts[0]]["ALT_TM"] - ref_tm) <= tm_diff:
                        primer["ALT_BINDS"].append(c)
                        para[c].update(alts[sorted_alts[0]])
                    para[c].pop("ALTERNATIVES")
                except KeyError:
                    try:
                        para[c].pop("ALTERNATIVES")
                    except KeyError:
                        pass
                except IndexError:
                    try:
                        para[c].pop("ALTERNATIVES")
                    except KeyError:
                        pass
    except KeyError:
        pass
    if outp:
        with open(os.path.join(
                primer3_output_DIR, output_file), "w") as outfile:
            json.dump(primer_dic, outfile, indent=1)
    return primer_dic


def score_paralog_primers(primer_dict, output_file, primer3_output_DIR,
                          ext, mask_penalty, species, backbone, outp=1):
    """
    Score primers in a dictionary according to a scoring matrix.

    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last.
    """
    primers = primer_dict["primer_information"]
    extension = (ext == "extension")
    # primer scoring coefficients were calculated based on
    # linear models of various parameters and provided as a dict
    with open("/opt/resources/mip_scores.dict", "rb") as infile:
        linear_coefs = pickle.load(infile)
    # the model was developed using specific reaction conditions as below.
    # actual conditions may be different from these but we'll use these
    # for the model.
    na = 25  # Sodium concentration
    mg = 10  # magnesium concentration
    conc = 0.04  # oligo concentration

    # get extension arm sequence
    if extension:
        for p in primers:
            extension_arm = primers[p]["SEQUENCE"]
            # calculate gc content of extension arm
            extension_gc = calculate_gc(extension_arm)
            # count lowercase masked nucleotides. These would likely be masked
            # for variation underneath.
            extension_lowercase = sum([c.islower() for c in extension_arm])
            # calculate TM with the model parameters for TM
            ext_TM = primer3.calcTm(extension_arm, mv_conc=na, dv_conc=mg,
                                    dna_conc=conc, dntp_conc=0)
            # create a mip parameter dict
            score_features = {"extension_gc": extension_gc,
                              "extension_lowercase": extension_lowercase,
                              "ext_TM": ext_TM}
            # calculate primer score using the linear model provided
            tech_score = 0
            for feature in score_features:
                degree = linear_coefs[feature]["degree"]
                primer_feature = score_features[feature]
                poly_feat = [pow(primer_feature, i) for i in range(degree + 1)]
                tech_score += sum(linear_coefs[feature]["coef"] * poly_feat)
                tech_score += linear_coefs[feature]["intercept"]
            primers[p]["SCORE"] = tech_score

    # get ligation arm parameters
    else:
        for p in primers:
            ligation_arm = primers[p]["SEQUENCE"]
            # calculate gc content of extension arm
            ligation_gc = calculate_gc(ligation_arm)
            # only the 3' end of the ligation arm was important in terms of
            # lowercase masking.
            ligation_lowercase_end = sum([c.islower()
                                          for c in ligation_arm[-5:]])
            # calculate TM of ligation sequence (actual ligation probe arm)
            # agains probe backbone.
            ligation_bb_TM = primer3.calcHeterodimerTm(
                reverse_complement(ligation_arm), backbone,
                mv_conc=na, dv_conc=mg, dna_conc=conc, dntp_conc=0)
            # create a mip parameter dict
            score_features = {"ligation_gc": ligation_gc,
                              "ligation_lowercase_end": ligation_lowercase_end,
                              "ligation_bb_TM": ligation_bb_TM}
            # calculate primer score using the linear model provided
            tech_score = 0
            for feature in score_features:
                degree = linear_coefs[feature]["degree"]
                primer_feature = score_features[feature]
                poly_feat = [pow(primer_feature, i) for i in range(degree + 1)]
                tech_score += sum(linear_coefs[feature]["coef"] * poly_feat)
                tech_score += linear_coefs[feature]["intercept"]
            primers[p]["SCORE"] = tech_score

    if outp:
        # write dictionary to json file
        outfile = open(os.path.join(primer3_output_DIR, output_file), "w")
        json.dump(primer_dict, outfile, indent=1)
        outfile.close()
    return primer_dict


def filter_primers(primer_dict, output_file,
                   primer3_output_DIR, n, bin_size, outp=1):
    """
    Filter primers so that only top n scoring primers remain for each bin.

    Primers are divided into bins of the given size based on the 3' end of
    the primer. Only top performing n primers ending in the same bin will
    remain after filtering.
    For example, bin_size=3 and n=1 would chose the best scoring primer
    among primers that end within 3 bps of each other.
    """
    # load extension and ligation primers from file
    template_seq = primer_dict["sequence_information"]["SEQUENCE_TEMPLATE"]
    template_len = len(template_seq)
    forward_bins = {}
    reverse_bins = {}
    for i in range(template_len//bin_size + 1):
        forward_bins[i] = []
        reverse_bins[i] = []
    for primer in list(primer_dict["primer_information"].keys()):
        # get primer orientation
        ori = primer_dict["primer_information"][primer]["ORI"]
        # get primer start coordinate
        start = int(primer_dict["primer_information"][primer]
                    ["COORDINATES"].split(",")[0])
        primer_len = int(primer_dict["primer_information"][primer]
                         ["COORDINATES"].split(",")[1])
        if ori == "forward":
            end = start + primer_len - 1
        elif ori == "reverse":
            end = start - primer_len + 1
        # which bin the start coordinate falls into
        end_bin = end//bin_size
        # get primer score
        score = primer_dict["primer_information"][primer]["SCORE"]
        # append the primer name/score to appropriate bin dic
        if ori == "forward":
            forward_bins[end_bin].append([primer, score])
        elif ori == "reverse":
            reverse_bins[end_bin].append([primer, score])
    best_primer_dict = {}
    best_primer_dict["sequence_information"] = primer_dict[
        "sequence_information"]
    best_primer_dict["primer_information"] = {}
    # find best scoring mips in each forward bin
    for key in forward_bins:
        # sort primers for score
        primer_set = sorted(forward_bins[key], key=itemgetter(1))
        # get best scoring primers (all primers if there are less than n)
        if len(primer_set) < n:
            best_primers = primer_set
        else:
            best_primers = primer_set[-n:]
        # add best primers do dictionary
        for primers in best_primers:
            primer_name = primers[0]
            best_primer_dict["primer_information"][primer_name] = primer_dict[
                "primer_information"][primer_name]
    # find best scoring mips in each reverse bin
    for key in reverse_bins:
        # sort primers for score
        primer_set = sorted(reverse_bins[key], key=itemgetter(1))
        # get best scoring primers (all primers if there are less than n)
        if len(primer_set) < n:
            best_primers = primer_set
        else:
            best_primers = primer_set[-n:]
        # add best primers do dictionary
        for primers in best_primers:
            primer_name = primers[0]
            best_primer_dict["primer_information"][primer_name] = primer_dict[
                "primer_information"][primer_name]
    # write new dic to file
    if outp:
        with open(os.path.join(
                primer3_output_DIR, output_file), "w") as outfile:
            json.dump(best_primer_dict, outfile, indent=1)
    return best_primer_dict


def pick_paralog_primer_pairs(extension, ligation, output_file,
                              primer3_output_DIR, min_size, max_size,
                              alternative_arms, region_insertions,
                              subregion_name, outp=1):
    """Pick primer pairs satisfying a given size range."""
    # assign primer information dictionaries to a shorter name
    ext = extension["primer_information"]
    lig = ligation["primer_information"]
    # check if extension and ligation dictionaries have primers
    if len(ext) == 0:
        return 1
    if len(lig) == 0:
        return 1
    # create a primer pairs dic. This dictionary is similar to primer dic
    primer_pairs = {}
    # has the same sequence_information key:value pairs
    primer_pairs["sequence_information"] = {}
    # has pair information key instead of primer_information
    primer_pairs["pair_information"] = {}
    # populate sequence information (same as extension or ligation)
    primer_pairs["sequence_information"]['SEQUENCE_TEMPLATE'] = extension[
        "sequence_information"]['SEQUENCE_TEMPLATE']
    primer_pairs["sequence_information"]['SEQUENCE_EXCLUDED_REGION'] = (
        extension["sequence_information"]['SEQUENCE_EXCLUDED_REGION']
    )
    primer_pairs["sequence_information"]['SEQUENCE_TARGET'] = extension[
        "sequence_information"]['SEQUENCE_TARGET']
    primer_pairs["sequence_information"]['SEQUENCE_ID'] = extension[
        "sequence_information"]['SEQUENCE_ID']
    # pick primer pairs
    for e in ext.keys():
        # extension primer information for this mip will be e_info
        e_info = ext[e]
        # get primer coordinates
        ext_start = e_info["GENOMIC_START"]
        ext_end = e_info["GENOMIC_END"]
        # get primer orientation
        ext_ori = ext_end > ext_start
        # if end is greater than start then it is a left(fw) primer,
        # and ext_ori is True.
        # get coordinates of this primer in paralog copies.
        ep_info = e_info["PARALOG_COORDINATES"]
        # the paralogs bound by primer according to bowtie mapping
        e_binds = e_info["BOWTIE_BINDS"]
        # paralogs that were not bound by the primer and alt primers were
        # designed.
        e_alt_binds = e_info["ALT_BINDS"]
        # find a ligation primer
        for l in list(lig.keys()):
            l_info = lig[l]
            # get primer coordinates
            lig_start = l_info["GENOMIC_START"]
            lig_end = l_info["GENOMIC_END"]
            # get orientation of primer
            lig_ori = lig_end < lig_start
            # if end is less than start, it is a right primer
            # create a list for start and end coordinates
            coord = []
            # continue only if the two orientations have the same value
            if lig_ori == ext_ori:
                # check if relative positions of primers are correct
                if ext_ori:
                    # ligation end should be greater than extension end
                    # for forward pairs
                    position = lig_end > ext_end
                else:
                    # extension end should be greater than ligation end
                    # for reverse pairs
                    position = ext_end > lig_end
                # get pair information if relative positions of primers are
                # correct
                if position:
                    coord = [ext_start, ext_end, lig_start, lig_end]
                    coord.sort()
                    prod_size = coord[-1] - coord[0] + 1
                    pairs = {}
                    # get paralogus coordinates
                    lp_info = l_info["PARALOG_COORDINATES"]
                    l_binds = l_info["BOWTIE_BINDS"]
                    l_alt_binds = l_info["ALT_BINDS"]
                    # find the paralogs that are hybridized by both primers
                    # start with paralog copies that are bound by the
                    # original primers (not alts).
                    paralogs = list(set(l_binds).intersection(e_binds))
                    for p in paralogs:
                        try:
                            p_coord = []
                            ep_start = ep_info[p]["BOWTIE_START"]
                            ep_end = ep_info[p]["BOWTIE_END"]
                            ep_ori = ep_end > ep_start
                            lp_start = lp_info[p]["BOWTIE_START"]
                            lp_end = lp_info[p]["BOWTIE_END"]
                            lp_ori = lp_end < lp_start
                            lp_chrom = lp_info[p]["CHR"]
                            if lp_ori == ep_ori:
                                if lp_ori:
                                    p_position = lp_end > ep_end
                                    pair_ori = "forward"
                                else:
                                    p_position = lp_end < ep_end
                                    pair_ori = "reverse"
                                if p_position:
                                    p_coord = [ep_start, ep_end,
                                               lp_start, lp_end]
                                    p_coord.sort()
                                    prod_size = p_coord[-1] - p_coord[0] + 1
                                    pairs[p] = {
                                        "capture_size": prod_size,
                                        "extension_start": ep_start,
                                        "extension_end": ep_end,
                                        "ligation_start": lp_start,
                                        "ligation_end": lp_end,
                                        "mip_start": p_coord[0],
                                        "mip_end": p_coord[3],
                                        "capture_start": p_coord[1] + 1,
                                        "capture_end": p_coord[2] - 1,
                                        "chrom": lp_chrom,
                                        "orientation": pair_ori
                                    }
                        except KeyError:
                            continue
                    # check if any pairs' product is within size limits
                    # taking into account reported insertions within
                    # the target region. If there are insertions, we reduce
                    # the max size to accomodate those insertions.
                    # Deletions are handled differently because their impact
                    # on the captures will be different. Any deletion that
                    # is small enough to be captured will still be captured
                    # without any alterations. However the capture size will
                    # become smaller, which is not detrimental.
                    pair_found = 0
                    captured_copies = []
                    for p in list(pairs.keys()):
                        if not region_insertions.empty:
                            max_insertion_size = region_insertions.loc[
                                (region_insertions["copy_chrom"]
                                 == pairs[p]["chrom"])
                                & (region_insertions["copy_begin"]
                                   > pairs[p]["capture_start"])
                                & (region_insertions["copy_end"]
                                   < pairs[p]["capture_end"]),
                                "max_size"].sum()
                        else:
                            max_insertion_size = 0
                        adjusted_max_size = max_size - max_insertion_size
                        if adjusted_max_size < (min_size/2):
                            continue
                        # we do not have to adsjust min_size unless the max
                        # size get too close to min_size, in which case
                        # we leave a 30 bp distance between min an max so
                        # that we're not very limited in primer  pair choices.
                        adjusted_min_size = min(adjusted_max_size - 30,
                                                min_size)
                        if (adjusted_max_size
                                >= pairs[p]["capture_size"]
                                >= adjusted_min_size):
                            captured_copies.append(p)
                            pair_found = 1
                    if pair_found:
                        # if a pair is found for any copy
                        # remove minimum size restriction for other copies
                        for p in list(pairs.keys()):
                            if p in captured_copies:
                                continue
                            if not region_insertions.empty:
                                max_insertion_size = region_insertions.loc[
                                    (region_insertions["copy_chrom"]
                                     == pairs[p]["chrom"])
                                    & (region_insertions["copy_begin"]
                                       > pairs[p]["capture_start"])
                                    & (region_insertions["copy_end"]
                                       < pairs[p]["capture_end"]),
                                    "max_size"].sum()
                            else:
                                max_insertion_size = 0
                            adjusted_max_size = max_size - max_insertion_size
                            if adjusted_max_size < (min_size/2):
                                continue
                            if (adjusted_max_size
                                    >= pairs[p]["capture_size"] >= 0):
                                captured_copies.append(p)
                        # C0 must be in the captured copies because the
                        # reference copy is used for picking mip sets
                        if "C0" not in captured_copies:
                            continue
                        # create a pair name as
                        # PAIR_extension primer number_ligation primer number
                        ext_name = e.split('_')[2]
                        lig_name = l.split('_')[2]
                        pair_name = ("PAIR_" + subregion_name + "_" + ext_name
                                     + "_" + lig_name)
                        if ext_ori:
                            orientation = "forward"
                            pair_name = pair_name + "_F"
                        else:
                            orientation = "reverse"
                            pair_name = pair_name + "_R"
                        primer_pairs["pair_information"][pair_name] = {
                            "pairs": pairs,
                            "extension_primer_information": ext[e],
                            "ligation_primer_information": lig[l],
                            "orientation": orientation,
                            "captured_copies": captured_copies
                        }
                        # Check if there are any paralog copies that require
                        # alt primers to be used. If so, create those pairs.
                        alt_paralogs = list((set(l_alt_binds).union(
                                            e_alt_binds)).difference(paralogs))
                        alts = {}
                        for a in alt_paralogs:
                            try:
                                alt_arms = []
                                p_coord = []
                                # check if the extension primer is the
                                # original or alt.
                                if ep_info[a]["BOWTIE_BOUND"]:
                                    ep_start = ep_info[a]["BOWTIE_START"]
                                    ep_end = ep_info[a]["BOWTIE_END"]
                                else:
                                    try:
                                        ep_start = ep_info[a]["ALT_START"]
                                        ep_end = ep_info[a]["ALT_END"]
                                        alt_arms.append("extension")
                                    except KeyError:
                                        continue
                                ep_ori = ep_end > ep_start
                                # check if ligation primer is the original
                                # or alternative designed.
                                if lp_info[a]["BOWTIE_BOUND"]:
                                    lp_start = lp_info[a]["BOWTIE_START"]
                                    lp_end = lp_info[a]["BOWTIE_END"]
                                else:
                                    try:
                                        lp_start = lp_info[a]["ALT_START"]
                                        lp_end = lp_info[a]["ALT_END"]
                                        alt_arms.append("ligation")
                                    except KeyError:
                                        continue
                                lp_ori = lp_end < lp_start
                                lp_chrom = lp_info[a]["CHR"]
                                if lp_ori == ep_ori:
                                    if lp_ori:
                                        p_position = lp_end > ep_end
                                        pair_ori = "forward"
                                    else:
                                        p_position = lp_end < ep_end
                                        pair_ori = "reverse"
                                    if p_position:
                                        p_coord = [ep_start, ep_end,
                                                   lp_start, lp_end]
                                        p_coord.sort()
                                        prod_size = (p_coord[-1]
                                                     - p_coord[0] + 1)
                                        alts[a] = {
                                            "capture_size": prod_size,
                                            "extension_start": ep_start,
                                            "extension_end": ep_end,
                                            "ligation_start": lp_start,
                                            "ligation_end": lp_end,
                                            "mip_start": p_coord[0],
                                            "mip_end": p_coord[3],
                                            "capture_start": p_coord[1] + 1,
                                            "capture_end": p_coord[2] - 1,
                                            "chrom": lp_chrom,
                                            "orientation": pair_ori,
                                            "alternative_arms": alt_arms
                                        }
                            except KeyError:
                                # if extension or ligation primer coordinates
                                # are not available for the paralog copy
                                # for any reason, e.g. the copy does not align
                                # to the ref for this primer, there will be
                                # a key error and it should be caught in this
                                # block.
                                continue
                        # check if any pairs' product is within size limits
                        captured_copies = []
                        for a in list(alts.keys()):
                            # does it satisfy arm setting?
                            good_alt = 0
                            # "any" means both ligation and extension arms
                            # are allowed to have alt sequences.
                            if alternative_arms == "any":
                                good_alt = 1
                            # if only one arm is allowed to have alt sequence,
                            # it could be specified as "one" or the specific
                            # arm (extension or ligation).
                            elif ((len(alts[a]["alternative_arms"]) == 1)
                                  and ((alternative_arms
                                        == alts[a]["alternative_arms"][0])
                                       or (alternative_arms == "one"))):
                                good_alt = 1
                            # if the alt capture is valid, check the capture
                            # size and determined if it is likely to be
                            # captured.
                            if good_alt:
                                if not region_insertions.empty:
                                    max_insertion_size = region_insertions.loc[
                                        (region_insertions["copy_chrom"]
                                         == alts[a]["chrom"])
                                        & (region_insertions["copy_begin"]
                                           > alts[a]["capture_start"])
                                        & (region_insertions["copy_end"]
                                           < alts[a]["capture_end"]),
                                        "max_size"].sum()
                                else:
                                    max_insertion_size = 0
                                adjusted_max_size = (max_size
                                                     - max_insertion_size)
                                if adjusted_max_size < (min_size/2):
                                    continue
                                if (adjusted_max_size
                                        >= alts[a]["capture_size"] >= 0):
                                    captured_copies.append(a)
                                    primer_pairs["pair_information"][
                                        pair_name]["pairs"][a] = alts[a]
                        primer_pairs["pair_information"][pair_name][
                            "alt_copies"] = captured_copies
    # return if no pairs found
    if len(primer_pairs["pair_information"]) == 0:
        # No primer pairs found.
        return 1
    # write dict to file in primer_output_DIR
    if outp:
        with open(os.path.join(
               primer3_output_DIR, output_file), 'w') as outfile:
            json.dump(primer_pairs, outfile, indent=1)
    return primer_pairs


def add_capture_sequence(primer_pairs, output_file, primer3_output_DIR,
                         species, outp=1):
    """
    Extract the sequence between primers.

    Get captured sequence using the primer coordinates.
    """
    capture_keys = set()
    for p_pair in primer_pairs["pair_information"]:
        pairs = primer_pairs["pair_information"][p_pair]["pairs"]
        for p in pairs:
            paralog_key = pairs[p]["chrom"] + ":" + str(pairs[p][
                "capture_start"]) + "-" + str(pairs[p]["capture_end"])
            pairs[p]["capture_key"] = paralog_key
            capture_keys.add(paralog_key)
    capture_sequence_dic = get_fasta_list(capture_keys, species)
    for p_pair in primer_pairs["pair_information"]:
        pairs = primer_pairs["pair_information"][p_pair]["pairs"]
        for p in pairs:
            if pairs[p]["orientation"] == "forward":
                pairs[p]["capture_sequence"] = capture_sequence_dic[pairs[p][
                    "capture_key"]]
            else:
                pairs[p]["capture_sequence"] = reverse_complement(
                    capture_sequence_dic[pairs[p]["capture_key"]]
                )
    if outp:
        with open(os.path.join(
                primer3_output_DIR, output_file), "w") as outfile:
            json.dump(primer_pairs, outfile, indent=1)
    return primer_pairs


def make_mips(pairs, output_file, primer3_output_DIR, mfold_input_DIR,
              backbone, outp=1):
    """
    Make mips from primer pairs.

    Take the reverse complement of ligation primer sequence, add the backbone
    sequence and the extension primer. Standard backbone is used if none
    specified.
    Add a new key to each primer pair:
    "mip_information" with a dictionary that has SEQUENCE key
    and mip sequence as value.
    """
    # check if the primer dictionary is empty
    if len(pairs["pair_information"]) == 0:
        return 1
    # get primer sequences for each primer pair
    for primers in pairs["pair_information"]:
        extension_sequence = pairs["pair_information"][primers][
            "extension_primer_information"]["SEQUENCE"]
        ligation_sequence = pairs["pair_information"][primers][
            "ligation_primer_information"]["SEQUENCE"]
        # reverse complement ligation primer
        ligation_rc = reverse_complement(ligation_sequence)
        # add sequences to make the mip
        mip_sequence = ligation_rc + backbone + extension_sequence
        # create a dictionary to hold mip information
        mip_dic = {"ref": {"SEQUENCE": mip_sequence,
                           "captures": copy.deepcopy(
                               pairs["pair_information"][primers]
                               ["captured_copies"]
                           )}}
        # create alternative mips where necessary
        if "alt_copies" in list(pairs["pair_information"][primers].keys()):
            alt_sequences = {}
            alt_counter = 0
            alt = pairs["pair_information"][primers]["alt_copies"]
            p_para = pairs["pair_information"][primers]["pairs"]
            e_para = pairs["pair_information"][primers][
                "extension_primer_information"]["PARALOG_COORDINATES"]
            l_para = pairs["pair_information"][primers][
                "ligation_primer_information"]["PARALOG_COORDINATES"]
            # since alt primers are created for each copy, it is possible
            # that some copies have the same primer pair. Pick just one
            # such pair and remove the others.
            for a in alt:
                if "extension" in p_para[a]["alternative_arms"]:
                    extension_sequence = e_para[a]["ALT_SEQUENCE"].upper()
                if "ligation" in p_para[a]["alternative_arms"]:
                    ligation_sequence = l_para[a]["ALT_SEQUENCE"].upper()
                value_found = 0
                # search through already created alt pairs to see if this one
                # is already there.
                for key, value in list(alt_sequences.items()):
                    if ([extension_sequence, ligation_sequence]
                            == value["sequences"]):
                        value_found = 1
                        # add the copy name to the dict and not create
                        # a new key for this copy.
                        value["copies"].append(a)
                        break
                # create new entry if this alt pair is new
                if not value_found:
                    alt_sequences[alt_counter] = {
                        "sequences": [extension_sequence, ligation_sequence],
                        "copies": [a]
                    }
                    alt_counter += 1
            # create mip sequence and dict for the alt pairs
            for alt_pair in alt_sequences:
                seq_dic = alt_sequences[alt_pair]["sequences"]
                alt_copies = alt_sequences[alt_pair]["copies"]
                # reverse complement ligation primer
                ligation_rc = reverse_complement(seq_dic[1])
                # add sequences to make the mip
                mip = ligation_rc + backbone + seq_dic[0]
                mip_dic["alt" + str(alt_pair)] = {"SEQUENCE": mip,
                                                  "captures": alt_copies}
        pairs["pair_information"][primers]["mip_information"] = mip_dic

    # write mip sequences to a fasta file in mfold_input_DIR
    # to check hairpin formation
    with open(os.path.join(mfold_input_DIR, output_file), "w") as outfile:
        for primers in pairs["pair_information"]:
            outline = (">" + primers + "\n" + pairs["pair_information"]
                       [primers]["mip_information"]["ref"]['SEQUENCE'] + "\n")
            outfile.write(outline)
    # write mip dictionary to file in primer3_output_DIR
    if outp:
        outfile = open(os.path.join(primer3_output_DIR, output_file), 'w')
        json.dump(pairs, outfile, indent=1)
        outfile.close()
    return pairs


def check_hairpin(pairs, output_file, settings, output_dir, outp=1):
    """Check possible hairpin formation in MIP probe.

    Calculate possible hiybridization between the MIP arms or between the MIP
    arms and the probe backbone. Remove MIPs with likely hairpins.
    """
    pairs = copy.deepcopy(pairs)
    # get Na, Mg and oligo concentrations these are specified in M but primer3
    # uses mM for ions and nM for oligos, so those will be adjusted.
    Na = float(settings["mip"]["Na"]) * 1000
    Mg = float(settings["mip"]["Mg"]) * 1000
    conc = float(settings["mip"]["oligo_conc"]) * pow(10, 9)
    # number of mips will be used to determine the bacbone concentration
    mip_count = int(settings["mip"]["mipset_size"])
    # get TM thresholds for hairpins, arm tms should be the same
    # otherwise we'll use the lower of the two
    ext_arm_tm = float(settings["extension"]["hairpin_tm"])
    lig_arm_tm = float(settings["ligation"]["hairpin_tm"])
    arm_tm = min([ext_arm_tm, lig_arm_tm])
    # backbone tm will be used for interactions between arms and
    # all the backbones (from other mips as well). This will cause a higher
    # tm since the backbones will be more in concentration, so it could
    # make sense to keep this threshold high. On the other hand, eliminating
    # even low likelyhood interactions could be useful.
    backbone_tm = float(settings["mip"]["hairpin_tm"])
    backbone_name = settings["mip"]["backbone"]
    backbone = mip_backbones[backbone_name]
    # go through mips and calculate hairpins
    # we will calculate hairpins by looking at TMs between arm sequences
    # and backbone sequences since the whole MIP sequence is too long
    # for nearest neighbor calculations (at least for primer3 implementation).
    for p in list(pairs["pair_information"].keys()):
        pair_dict = pairs["pair_information"][p]
        mip_dict = pair_dict["mip_information"]
        # for each primer pair we can have a number of mips due to paralog
        # copies having alternative mips. We'll go through each mip.
        for m in list(mip_dict.keys()):
            mip_seq = mip_dict[m]["SEQUENCE"]
            # extract arm and backbone sequences from the mip sequence
            lig = mip_seq[:mip_seq.index(backbone)]
            ext = mip_seq[mip_seq.index(backbone) + len(backbone):]
            bb = backbone.replace("N", "")
            # calculate dimer TMs between sequence combinations
            ext_lig = calcHeterodimerTm(ext, lig, mv_conc=Na, dv_conc=Mg,
                                        dntp_conc=0, dna_conc=conc)
            bb_ext_arm = calcHeterodimerTm(ext, bb, mv_conc=Na, dv_conc=Mg,
                                           dntp_conc=0, dna_conc=conc)
            bb_lig_arm = calcHeterodimerTm(lig, bb, mv_conc=Na, dv_conc=Mg,
                                           dntp_conc=0, dna_conc=conc)
            # take the maximum TM for hairpin threshold comparison
            arms = max([ext_lig, bb_ext_arm, bb_lig_arm])
            # calculate TM between arms and the whole reaction backbones
            # backbone concentration will be more for this calculation.
            bb_ext = calcHeterodimerTm(ext, bb, mv_conc=Na, dv_conc=Mg,
                                       dntp_conc=0, dna_conc=conc * mip_count)
            bb_lig = calcHeterodimerTm(lig, bb, mv_conc=Na, dv_conc=Mg,
                                       dntp_conc=0, dna_conc=conc * mip_count)
            bb_temp = max([bb_ext, bb_lig])
            # if either hairpin tms is higher than the limit, remove the mip
            # and remove the paralog copy that is supposed to be captured
            # by this specific mip from the pair dictionary.
            if (arms > arm_tm) or (bb_temp > backbone_tm):
                lost_captures = mip_dict[m]["captures"]
                mip_copies = pair_dict["captured_copies"]
                mip_copies = list(set(mip_copies).difference(lost_captures))
                pair_dict["captured_copies"] = mip_copies
                alt_copies = pair_dict["alt_copies"]
                alt_copies = list(set(alt_copies).difference(lost_captures))
                pair_dict["alt_copies"] = alt_copies
                mip_dict.pop(m)
            else:
                mip_dict[m]["Melting Temps"] = {"arms_hp": ext_lig,
                                                "ext_hp": bb_ext_arm,
                                                "lig_hp": bb_lig_arm,
                                                "ext_backbone": bb_ext,
                                                "lig_backbone": bb_lig}
        if len(mip_dict) == 0:
            pairs["pair_information"].pop(p)
    for p in pairs["pair_information"].keys():
        pair_dict = pairs["pair_information"][p]
        hp_dict = pair_dict["hairpin"] = {}
        mip_dict = pair_dict["mip_information"]
        for m in mip_dict:
            hp_dict[m] = mip_dict[m]["Melting Temps"]
    if outp:
        output_file = os.path.join(output_dir, output_file)
        with open(output_file, "w") as outfile:
            json.dump(pairs, outfile)
    return pairs


def filter_mips(mip_dic, bin_size, mip_limit):
    """
    Filter MIPs covering similar regions.

    Filter MIPs so that only top scoring mip ending within the "bin_size"
    nucleotides on the same strand remain.
    """
    # load extension and ligation primers from file
    shuffled = list(mip_dic.keys())
    random.shuffle(shuffled)
    for m in shuffled:
        if len(mip_dic) <= mip_limit:
            return
        try:
            m_start = mip_dic[m].mip["C0"]["capture_start"]
            m_end = mip_dic[m].mip["C0"]["capture_end"]
            m_func = mip_dic[m].func_score
            m_tech = mip_dic[m].tech_score
            m_ori = mip_dic[m].mip["C0"]["orientation"]
            for n in shuffled:
                if len(mip_dic) <= mip_limit:
                    return
                try:
                    if mip_dic[m].name != mip_dic[n].name:
                        n_start = mip_dic[n].mip["C0"]["capture_start"]
                        n_end = mip_dic[n].mip["C0"]["capture_end"]
                        n_func = mip_dic[n].func_score
                        n_tech = mip_dic[n].tech_score
                        n_ori = mip_dic[n].mip["C0"]["orientation"]
                        if (((abs(n_start - m_start) <= bin_size)
                             and (abs(n_end - m_end) <= bin_size))
                                and (m_ori == n_ori)):
                            if (m_tech + m_func) >= (n_tech + n_func):
                                mip_dic.pop(n)
                            else:
                                mip_dic.pop(m)
                                break
                except KeyError:
                    continue
        except KeyError:
            continue
    return


def compatible_mip_check(m1, m2, overlap_same, overlap_opposite):
    d = m1.mip_dic
    # get m1 coordinates
    ext_start = d["extension_primer_information"]["GENOMIC_START"]
    ext_end = d["extension_primer_information"]["GENOMIC_END"]
    lig_start = d["ligation_primer_information"]["GENOMIC_START"]
    lig_end = d["ligation_primer_information"]["GENOMIC_END"]
    # get mip1 orientation
    ori = d["orientation"]

    # get m2 coordinates
    m = m2.mip_dic
    next_ext_start = m["extension_primer_information"]["GENOMIC_START"]
    next_ext_end = m["extension_primer_information"]["GENOMIC_END"]
    next_lig_start = m["ligation_primer_information"]["GENOMIC_START"]
    next_lig_end = m["ligation_primer_information"]["GENOMIC_END"]
    # get mip2 orientation
    next_ori = m["orientation"]

    if ori == next_ori:
        m1_start = min([ext_start, ext_end, lig_start, lig_end])
        m1_end = max([ext_start, ext_end, lig_start, lig_end])
        m2_start = min([next_ext_start, next_ext_end, next_lig_start,
                       next_lig_end])
        m2_end = max([next_ext_start, next_ext_end, next_lig_start,
                      next_lig_end])
        ol = overlap([m1_start, m1_end], [m2_start, m2_end])
        if len(ol) == 0:
            return True
        else:
            return (ol[1] - ol[0] + 1) <= overlap_same
    else:
        m1_set = set(list(range(min([ext_start, ext_end]),
                                max([ext_start, ext_end]) + 1))
                     + list(range(min([lig_start, lig_end]),
                                  max([lig_start, lig_end]) + 1)))

        m2_set = set(list(range(min([next_ext_start, next_ext_end]),
                                max([next_ext_start, next_ext_end]) + 1))
                     + list(range(min([next_lig_start, next_lig_end]),
                                  max([next_lig_start, next_lig_end]) + 1)))
        ol = len(m1_set.intersection(m2_set))
        return ol <= overlap_opposite


def compatible_chains(primer_file, mip_dict, primer3_output_DIR,
                      primer_out, output_file, must_bonus, set_copy_bonus,
                      overlap_same, overlap_opposite, outp, bin_size,
                      trim_increment, trim_limit, set_size, chain_mips,
                      intervals):
    try:
        with open(os.path.join(
                primer3_output_DIR, primer_file), "r") as infile:
            scored_mips = json.load(infile)
    except IOError:
        print("Primer file does not exist.")
        return 1
    else:
        # make a copy of the original mip dict to use in filtering
        temp_dict = copy.deepcopy(mip_dict)
        # create small subregions for binning MIPs and creating compatible
        # mip sets for smaller regions
        begin = intervals[0]
        end = intervals[1]
        bins = list(range(begin, end, bin_size))
        # if a single nucleotide is the target, the interval will be the
        # position of that nucleotide as [pos, pos] and the range will return
        # an empty list. In this case we'll crease a [pos, pos] list instead.
        if begin == end:
            bins = [begin, end]
        if bins[-1] != end:
            bins.append(end)
        num_bins = len(bins) - 1
        # group MIPs into bins. Bins can share MIPs.
        binned = {}
        for i in range(num_bins):
            binned[i] = {}
            bin_start = bins[i]
            bin_end = bins[i + 1]
            for k in temp_dict:
                cp = temp_dict[k].mip["C0"]
                cs = cp["capture_start"]
                ce = cp["capture_end"]
                if len(overlap([cs, ce], [bin_start, bin_end])) > 0:
                    binned[i][k] = temp_dict[k]
        # remove MIPs covering similar regions until we have only
        # "set_size" number of MIPs per bin.
        for i in binned:
            trim_size = 1
            while (trim_size <= trim_limit) and (len(binned[i]) > set_size):
                filter_mips(binned[i], trim_size, set_size)
                trim_size += trim_increment
        # create (in)compatibility lists for each MIP
        for k in list(scored_mips["pair_information"].keys()):
            # get coordinates of mip arms
            d = scored_mips["pair_information"][k]
            # extension arm start position
            es = d["extension_primer_information"]["GENOMIC_START"]
            # extension arm end position
            ee = d["extension_primer_information"]["GENOMIC_END"]
            # ligation arm start position
            ls = d["ligation_primer_information"]["GENOMIC_START"]
            # ligation arm end position
            le = d["ligation_primer_information"]["GENOMIC_END"]
            # get mip orientation
            ori = d["orientation"]
            # create an in/compatibility list
            incompatible = set()
            compatible = set()
            # loop through all mips to populate compatibility lists
            for mk in list(scored_mips["pair_information"].keys()):
                m = scored_mips["pair_information"][mk]
                # next MIP's extension arm start position
                nes = m["extension_primer_information"]["GENOMIC_START"]
                # next MIP's extension arm end position
                nee = m["extension_primer_information"]["GENOMIC_END"]
                # next MIP's ligation arm start position
                nls = m["ligation_primer_information"]["GENOMIC_START"]
                # next MIP's ligation arm end position
                nle = m["ligation_primer_information"]["GENOMIC_END"]
                # get mip orientation
                next_ori = m["orientation"]
                compat = 0
                next_compat = 0
                # check if the two mips are compatible in terms of
                # orientation and coordinates
                if ori == next_ori == "forward":
                    if (((ls < nls) and (ls < nes + overlap_same))
                            or ((ls > nls) and (es + overlap_same > nls))):
                        compat = 1
                elif ori == next_ori == "reverse":
                    if (((ls < nls) and (es < nls + overlap_same))
                            or ((ls > nls) and (ls + overlap_same > nes))):
                        compat = 1
                elif (ori == "forward") and (next_ori == "reverse"):
                    if ((ls < nls + overlap_opposite)
                            or (es + overlap_opposite > nes)):
                        compat = 1
                    elif ((es < nls) and (ee < nls + overlap_opposite)
                          and (le + overlap_opposite > nle)
                          and (ls < nee + overlap_opposite)):
                        compat = 1
                        next_compat = 1
                    elif ((es > nls) and (es + overlap_opposite > nle)
                          and (ee < nee + overlap_opposite)
                          and (le + overlap_opposite > nes)):
                        compat = 1
                elif (ori == "reverse") and (next_ori == "forward"):
                    if ((ls + overlap_opposite > nls)
                            or (es < nes + overlap_opposite)):
                        compat = 1
                    elif ((ls > nes) and (ls + overlap_opposite > nee)
                          and (le < nle + overlap_opposite)
                          and (ee + overlap_opposite > nls)):
                        compat = 1
                    elif ((ls < nes) and (le < nes + overlap_opposite)
                          and (ee + overlap_opposite > nee)
                          and (es < nle + overlap_opposite)):
                        compat = 1
                        next_compat = 1
                if not compat:
                    incompatible.add(mk)
                if next_compat:
                    compatible.add(mk)
            d["incompatible"] = incompatible
            d["compatible"] = compatible

        def compatible_recurse(l):
            """
            Take a list, l,  of numbers that represent a mip set with
            their corresponding "place" in the mip dictionary, and index
            number, i. Find the subset of mips in the rest of the list
            that are compatible with the mip at index i, using compatibility
            dictionary d. For each mip in the subset, find compatible mips
            in the rest of the list. Recurse until the subset does not have
            any mips. Append each compatible subset to a final result list, f.
            """
            # create a set of mips that are incompatible with any mip in
            # the starting list.
            incomp = set(l)
            for il in l:
                incomp.update(scored_mips["pair_information"][il][
                    "incompatible"])
            # create a set of mips that can be the "next" mip that can be
            # added to the mip list
            comp = scored_mips["pair_information"][l[-1]][
                "compatible"].difference(incomp).intersection(subset)
            # if there are mips that can be added, call compatible_recurse
            # function for each of those mips
            if len(comp) > 0:
                for n in comp:
                    compatible_recurse(l + [n])
            # stop recursing when the mip chain cannot be elongated
            else:
                mip_sets.append((l))

        keys = sorted(scored_mips["pair_information"],
                      key=lambda a: scored_mips["pair_information"][a]
                      ["pairs"]["C0"]["capture_start"])
        ms_dict = {}
        for i in binned:
            subset = binned[i]
            mip_sets = []
            for k in keys:
                if k in subset:
                    comp_list = scored_mips["pair_information"][k][
                        "compatible"].intersection(subset)
                    if len(comp_list) > 0:
                        # for each of the mips in the compatibility list,
                        for m in comp_list:
                            # check if these two mips are present in other sets
                            # if they are, then no need to pursue this branch
                            # anymore as the same branch will be in the other
                            # mip set as well
                            test_set = frozenset([k, m])
                            for p_set in mip_sets:
                                if test_set.issubset(set(p_set)):
                                    break
                            else:
                                # create an initial result list to be used by
                                # the compatible_recurse function
                                compatible_recurse([k, m])
                    else:
                        mip_sets.append(([k]))
            ms_dict[i] = mip_sets

        # define a funtcion for getting the mipset score and coverage
        def score_mipset(mip_set):
            # create a dic for diffs captured cumulatively by all
            # mips in the set
            merged_caps = []
            # create a list for mip scores based on mip sequence and
            # not the captured diffs
            mip_scores = []
            # create a list for what is captured by the set (only must
            # captures)
            must_captured = []
            # create a list for other targets captured
            targets_captured = []
            # a list for mip coordinates
            capture_coordinates = []
            for mip_key in mip_set:
                # extract the mip name
                # extract the captured diffs from the mip_dic and
                # append to capture list
                mip_obj = mip_dict[mip_key]
                uniq = mip_obj.capture_info["unique_captures"]
                merged_caps.extend(uniq)
                must_captured.extend(mip_obj.captures)
                targets_captured.extend(mip_obj.captured_targets)
                if ((mip_obj.tech_score > 0)
                        and (mip_obj.func_score > 0)):
                    mip_scores.append(
                        float(mip_obj.tech_score * mip_obj.func_score)
                        / 1000
                    )
                else:
                    mip_scores.append(
                        float(mip_obj.tech_score + mip_obj.func_score)
                        / 1000)
                mcoord = sorted(
                    [mip_obj.extension["C0"]["GENOMIC_START"],
                     mip_obj.ligation["C0"]["GENOMIC_START"],
                     mip_obj.extension["C0"]["GENOMIC_END"],
                     mip_obj.ligation["C0"]["GENOMIC_END"]]
                )
                capture_coordinates.append([mcoord[1] + 1,
                                            mcoord[2] - 1])
            merged_capture_coordinates = merge_overlap(
                capture_coordinates, 50)
            scp = len(set(merged_caps)) * set_copy_bonus
            must_set = list(set(must_captured))
            mb = len(must_set) * must_bonus
            total_score = mb + scp + sum(mip_scores)
            return total_score, merged_capture_coordinates

        # create a dictionary to hold mip sets and their scores
        mip_set_dict = {}
        for i in ms_dict:
            mip_set_dict[i] = {}
            bin_co = bins[i: i + 2]
            bin_size = bin_co[1] - bin_co[0] + 1
            for j in range(len(ms_dict[i])):
                ms = ms_dict[i][j]
                sc = score_mipset(ms)
                coverage = overlap(sc[1][0], bin_co)
                coverage = (coverage[1] - coverage[0] + 1) / bin_size
                mip_set_dict[i][j] = {"mip_list": ms, "score": sc[0],
                                      "coordinates": sc[1][0],
                                      "coverage": coverage}
        for i in mip_set_dict:
            iter_keys = list(mip_set_dict[i].keys())
            for j in iter_keys:
                try:
                    s1 = mip_set_dict[i][j]["mip_list"]
                    sc1 = mip_set_dict[i][j]["score"]
                    crd1 = mip_set_dict[i][j]["coordinates"]
                    cov1 = mip_set_dict[i][j]["coverage"]
                    for k in iter_keys:
                        if k == j:
                            continue
                        try:
                            s2 = mip_set_dict[i][k]["mip_list"]
                            sc2 = mip_set_dict[i][k]["score"]
                            crd2 = mip_set_dict[i][k]["coordinates"]
                            cov2 = mip_set_dict[i][k]["coverage"]
                            if check_redundant_region(crd1, crd2, spacer=0):
                                # if one set is to be removed pick the one
                                # with full coverage of the target region
                                # in case there is one
                                if chain_mips:
                                    if (cov1 == 1) and (cov2 < 1):
                                        mip_set_dict[i].pop(k)
                                    elif (cov2 == 1) and (cov1 < 1):
                                        mip_set_dict[i].pop(j)
                                        break
                                    # if both are covering the target
                                    # or if both are failing to cover
                                    # then pick the set with better score
                                    elif sc2 > sc1:
                                        mip_set_dict[i].pop(j)
                                        break
                                    else:
                                        mip_set_dict[i].pop(k)
                                # if chaining mip is not required
                                # pick the better scoring set
                                elif sc2 > sc1:
                                    mip_set_dict[i].pop(j)
                                    break
                                else:
                                    mip_set_dict[i].pop(k)
                        except KeyError:
                            continue
                except KeyError:
                    continue

        # merge compatible chains within each bin (to some extent)
        merged_sets = {}
        for i in mip_set_dict:
            mip_sets = set()
            for j in mip_set_dict[i]:
                mip_sets.add(frozenset(mip_set_dict[i][j]["mip_list"]))

            # these mip sets only contain mip chains. We can expand each
            # such set by merging with other sets after removing incompatible
            # mips from the second set.
            counter = 0
            for counter in range(5):
                new_mip_sets = set()
                for s1 in mip_sets:
                    inc = set()
                    for m in s1:
                        inc.update(scored_mips["pair_information"][m][
                            "incompatible"])
                    new_set = set(s1)
                    for s2 in mip_sets:
                        counter += 1
                        s3 = s2.difference(inc).difference(new_set)
                        if len(s3) > 0:
                            new_set.update(s3)
                            for m in new_set:
                                inc.update(scored_mips["pair_information"][m][
                                    "incompatible"])
                    new_mip_sets.add(frozenset(new_set))
                mip_sets = new_mip_sets
            if len(mip_sets) > 0:
                merged_sets[i] = mip_sets

        # combine mip sets in different bins
        # first, calculate how many combinations there will be
        combo_length = 1
        for i in merged_sets:
            combo_length *= len(merged_sets[i])
        # if too many combinations, reduce by picking the top 5 scoring
        # sets for each bin
        if combo_length > pow(10, 7):
            for i in list(merged_sets.keys()):
                top_sets = set(sorted(merged_sets[i],
                                      key=lambda a: score_mipset(a)[0],
                                      reverse=True)[:5])
                merged_sets[i] = top_sets
            combo_length = 1
            for i in merged_sets:
                combo_length *= len(merged_sets[i])
            # if still too many combinations, take the top set for each bin
            if combo_length > pow(10, 7):
                for i in list(merged_sets.keys()):
                    top_sets = set(sorted(merged_sets[i],
                                          key=lambda a: score_mipset(a)[0],
                                          reverse=True)[:1])
                    merged_sets[i] = top_sets

        # combine mip sets in different bins
        combined_sets = set()
        combo_list = list(itertools.product(
            *[merged_sets[i] for i in sorted(merged_sets)]))
        for l in combo_list:
            if len(l) == 1:
                m_set = set(l[0])
            else:
                m_set = set()
                for i in range(len(l) - 1):
                    s1 = l[i]
                    s2 = l[i + 1]
                    inc = set()
                    for m in s1:
                        inc.update(scored_mips["pair_information"][m][
                            "incompatible"])
                    s3 = s2.difference(inc)
                    m_set.update(s1.union(s3))
            combined_sets.add(frozenset(m_set))

        if outp:
            with open(os.path.join(
                    primer3_output_DIR, output_file), "w") as outfile:
                outfile.write("\n".join([",".join(s) for s in combined_sets])
                              + "\n")
        with open(os.path.join(
                primer3_output_DIR, primer_out), "wb") as outfile:
            pickle.dump(scored_mips, outfile)
    return combined_sets


def design_mips(design_dir, g):
    print(("Designing MIPs for ", g))
    try:
        Par = mod.Paralog(os.path.join(design_dir, g, "resources",
                                       g + ".rinfo"))
        Par.run_paralog()
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ",
                   Par.paralog_name))
        if Par.chain_mips:
            if Par.chained_mips:
                print(("All MIPs are chained for paralog ", Par.paralog_name))
            else:
                print(("MIPs are NOT chained for paralog ", Par.paralog_name))
    except Exception as e:
        print((g, str(e), " FAILED!!!"))
    return


def design_mips_worker(design_list):
    design_dir, g = design_list
    print(("Designing MIPs for ", g))
    try:
        rinfo_file = os.path.join(design_dir, g, "resources", g + ".rinfo")
        Par = mod.Paralog(rinfo_file)
        Par.run_paralog()
        if len(Par.mips) == 0:
            return
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ",
                   Par.paralog_name))
        if Par.chain_mips:
            if Par.chained_mips:
                print(("All MIPs are chained for paralog ", Par.paralog_name))
            else:
                print(("MIPs are NOT chained for paralog ", Par.paralog_name))
    except Exception as e:
        print((g, str(e), " FAILED!!!"))
        traceback.print_exc()
    return 0


def design_mips_multi(design_dir, g_list, num_processor):
    chore_list = [[design_dir, g] for g in g_list]
    res = []
    try:
        p = NoDaemonProcessPool(num_processor)
        p.map_async(design_mips_worker, chore_list, callback=res.append)
        p.close()
        p.join()
    except Exception as e:
        res.append(str(e))
    return res


def parasight(resource_dir,
              design_info_file,
              designed_gene_list=None,
              extra_extension=".extra",
              use_json=False):
    if not use_json:
        with open(design_info_file, "rb") as infile:
            design_info = pickle.load(infile)
    else:
        with open(design_info_file) as infile:
            design_info = json.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    pdf_dir = os.path.join(resource_dir, "pdfs")
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash", "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                        + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                        "-sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = os.path.join(design_info[t]["design_dir"], t,  t)
        backup_name = basename + ".extra"
        filtered_name = basename + "_filtered.pse"
        backup_list.append("scp " + backup_name + " " + backup_name + ".bak")
        backup_list.append("mv " + filtered_name + " " + backup_name)
        psname = basename + ".01.01.ps"
        pdfname = basename + ".pdf"
        gs_command = ("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                      + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                      "-sOutputFile=" + pdfname + " " + psname)
        if designed_gene_list is not None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".pdf")
        else:
            pdf_convert_list.append(t + ".pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".pdf "
                        + os.path.join(pdf_dir, t + ".pdf"))
        outlist = ["parasight76.pl",
                   "-showseq", basename + ".show",
                   "-extra", basename + extra_extension,
                   "-template", "/opt/resources/nolabel.pst",
                   "-precode file:" + basename + ".precode",
                   "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode", "w") as outfile:
            outfile.write("$opt{'filename'}='" + t
                          + "';&fitlongestline; &print_all (0,'"
                          + basename + "')")
    with open(os.path.join(resource_dir, "backup_commands"), "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(
            os.path.join(resource_dir, "parasight_commands"), "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(os.path.join(resource_dir, "gs_commands"), "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(os.path.join(resource_dir, "copy_commands"), "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(os.path.join(resource_dir, "convert_commands"), "w") as outfile:
        outfile.write("\n".join(pdf_merge_list))
    visualization_list = ["#!/usr/bin/env bash"]
    visualization_list.append("chmod +x backup_commands")
    visualization_list.append("./backup_commands")
    visualization_list.append("chmod +x parasight_commands")
    visualization_list.append("./parasight_commands")
    visualization_list.append("chmod +x gs_commands")
    visualization_list.append("./gs_commands")
    visualization_list.append("chmod +x copy_commands")
    visualization_list.append("./copy_commands")
    visualization_list.append("chmod +x convert_commands")
    visualization_list.append("./convert_commands")
    with open(os.path.join(resource_dir, "visualize.sh"), "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return


def parasight_print(resource_dir, design_dir, design_info_file,
                    designed_gene_list=None, extra_extension=".extra",
                    use_json=False, print_out=False):
    if not use_json:
        with open(design_info_file, "rb") as infile:
            design_info = pickle.load(infile)
    else:
        with open(design_info_file) as infile:
            design_info = json.load(infile)
    output_file = os.path.join(resource_dir, "parasight_print.txt")
    with open(output_file, "w") as outfile:
        for g in design_info:
            if (designed_gene_list is None) or (g in designed_gene_list):
                show_file = os.path.join(design_dir, g, g + ".show")
                extras_file = os.path.join(design_dir, g, g + extra_extension)
                line = ["parasight76.pl", "-showseq", show_file,
                        "-extra ", extras_file]
                if print_out:
                    print(" ".join(line))
                outfile.write(" ".join(line) + "\n")


###############################################################
# Data analysis related functions
###############################################################


def get_analysis_settings(settings_file):
    """Convert analysis settings file to dictionary."""
    settings = {}
    with open(settings_file) as infile:
        for line in infile:
            try:
                if not line.startswith("#"):
                    newline = line.strip().split("\t")
                    value = newline[1].split(",")
                    if len(value) == 1:
                        settings[newline[0]] = value[0]
                    else:
                        settings[newline[0]] = [v for v in value if v != ""]
            except Exception as e:
                print(("Formatting error in settings file, line {}"
                       "causing error '{}''").format(line, e))
                print(newline)
                return
    return settings


def write_analysis_settings(settings, settings_file):
    """Create a settings file from a settings dictionary."""
    outfile_list = [["# Setting Name", "Setting Value"]]
    for k, v in settings.items():
        if isinstance(v, list):
            val = ",".join(map(str, v))
        else:
            val = str(v)
        outfile_list.append([k, val])
    with open(settings_file, "w") as outfile:
        outfile.write("\n".join(["\t".join(o) for o in outfile_list]) + "\n")
    return


###############################################################################
# New contig based analysis for vcf generation
###############################################################################

def map_haplotypes(settings):
    """Bwa-map haplotypes from MIPWrangler output to the reference genome.

    Extract each unique haplotype sequence from the MIPWrangler output and
    map to reference genome. MIPWrangler maps the sequencing data to the MIPs
    used for an experiment based on the probe arms. We compare here whether
    the best genomic loci for a given haplotype matches to the MIPWrangler
    assignment. If not, we consider those off target and remove.
    """
    wdir = settings["workingDir"]
    haplotypes_fq_file = os.path.join(wdir, settings["haplotypesFastqFile"])
    haplotypes_sam_file = os.path.join(wdir, settings["haplotypesSamFile"])
    bwa_options = settings["bwaOptions"]
    call_info_file = settings["callInfoDictionary"]
    species = settings["species"]
    try:
        tol = int(settings["alignmentTolerance"])
    except KeyError:
        tol = 200
    # DATA EXTRACTION ###
    raw_results = pd.read_table(os.path.join(wdir,
                                             settings["mipsterFile"]))
    ##########################################################
    # Add the statistics for each haplotype to the data
    # such as how many samples had a given haplotype
    # and how many barcodes supported a given haplotype
    # Filter the haplotypes for those criteria to
    # remove possible noise and infrequent haplotypes
    ##########################################################
    # Haplotype Filters from the settings file
    haplotype_min_barcode_filter = int(settings["minHaplotypeBarcodes"])
    haplotype_min_sample_filter = int(settings["minHaplotypeSamples"])
    haplotype_min_sample_fraction_filter = float(
        settings["minHaplotypeSampleFraction"]
    )
    # Gather per haplotype data across samples
    hap_counts = raw_results.groupby(
        "haplotype_ID"
    )["barcode_count"].sum().reset_index().rename(
        columns={"barcode_count": "Haplotype Barcodes"})
    hap_sample_counts = raw_results.groupby("haplotype_ID")[
        "sample_name"].apply(lambda a: len(set(a))).reset_index().rename(
        columns={"sample_name": "Haplotype Samples"})
    num_samples = float(raw_results["sample_name"].unique().size)
    hap_sample_counts["Haplotype Sample Fraction"] = (
        hap_sample_counts["Haplotype Samples"] / num_samples
    )
    hap_counts = hap_counts.merge(hap_sample_counts)
    initial_hap_count = len(hap_counts)
    hap_counts = hap_counts.loc[(hap_counts["Haplotype Samples"]
                                 >= haplotype_min_sample_filter)
                                & (hap_counts["Haplotype Sample Fraction"]
                                   >= haplotype_min_sample_fraction_filter)
                                & (hap_counts["Haplotype Barcodes"]
                                   >= haplotype_min_barcode_filter)]
    print(("Out of {} initial haplotypes, {} were filtered using {}, {}, and "
           "{} as minimum total UMI count; number and fraction of samples "
           " the haplotype was observed in, respectively.").format(
               initial_hap_count, initial_hap_count - len(hap_counts),
               haplotype_min_barcode_filter, haplotype_min_sample_filter,
               haplotype_min_sample_fraction_filter))

    hap_df = raw_results.loc[raw_results["haplotype_ID"].isin(
        hap_counts["haplotype_ID"])].groupby(
        ["gene_name", "mip_name", "haplotype_ID"])[
        "haplotype_sequence"].first().reset_index()
    # fill in fake sequence quality scores for each haplotype. These scores
    # will be used for mapping only and the real scores for each haplotype
    # for each sample will be added later.This step is probably unnecessary
    # as the bwa mem algorithm does not seem to use the quality scores.
    hap_df["quality"] = hap_df["haplotype_sequence"].apply(
        lambda a: "H" * len(a))
    haps = hap_df.set_index("haplotype_ID").to_dict(orient="index")
    # BWA alignment
    # create a fastq file for bwa input
    with open(haplotypes_fq_file, "w") as outfile:
        for h in haps:
            outfile.write("@" + h + "\n")
            outfile.write(haps[h]["haplotype_sequence"] + "\n" + "+" + "\n")
            outfile.write(haps[h]["quality"] + "\n")
    # run bwa
    bwa(haplotypes_fq_file, haplotypes_sam_file, "sam", "", "", bwa_options,
        species)
    # process alignment output sam file
    header = ["haplotype_ID", "FLAG", "CHROM", "POS", "MAPQ", "CIGAR", "RNEXT",
              "PNEXT", "TLEN", "SEQ", "QUAL"]
    sam_list = []
    with open(haplotypes_sam_file) as infile:
        for line in infile:
            if not line.startswith("@"):
                newline = line.strip().split()
                samline = newline[:11]
                for item in newline[11:]:
                    value = item.split(":")
                    if value[0] == "AS":
                        samline.append(int(value[-1]))
                        break
                else:
                    samline.append(-5000)
                sam_list.append(samline)
    sam = pd.DataFrame(sam_list, columns=header + ["alignment_score"])
    # find alignment with the highest alignment score. We will consider these
    # the primary alignments and the source of the sequence.
    sam["best_alignment"] = (sam["alignment_score"] == sam.groupby(
        "haplotype_ID")["alignment_score"].transform("max"))
    # add MIP column to alignment results
    sam["MIP"] = sam["haplotype_ID"].apply(lambda a: a.split(".")[0])
    # create call_info data frame for all used probes in the experiment
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    call_df_list = []
    for g in call_info:
        for m in call_info[g]:
            if m in used_probes:
                mip_number = int(m.split("_")[-1][3:])
                sub_number = int(m.split("_")[-2][3:])
                for c in call_info[g][m]["copies"]:
                    call_dict = call_info[g][m]["copies"][c]
                    try:
                        call_dict.pop("genes")
                    except KeyError:
                        pass
                    try:
                        call_dict.pop("variants")
                    except KeyError:
                        pass
                    call_dict["gene"] = g
                    call_dict["MIP"] = m
                    call_dict["copy"] = c
                    call_dict["mip_number"] = mip_number
                    call_dict["sub_number"] = sub_number
                    call_df_list.append(pd.DataFrame(call_dict, index=[0]))
    call_df = pd.concat(call_df_list, ignore_index=True, sort=True)
    # combine alignment information with design information (call_info)
    haplotype_maps = call_df.merge(
        sam[["MIP", "haplotype_ID", "CHROM", "POS", "best_alignment",
             "alignment_score"]])
    haplotype_maps["POS"] = haplotype_maps["POS"].astype(int)
    haplotype_maps = haplotype_maps.merge(
        hap_df[["haplotype_ID", "haplotype_sequence"]])
    # determine which haplotype/mapping combinations are for intended targets
    # first, compare mapping coordinate to the MIP coordinate to see  if
    # a MIP copy matches with the alignment.
    haplotype_maps["aligned_copy"] = (
        (haplotype_maps["CHROM"] == haplotype_maps["chrom"])
        & (abs(haplotype_maps["POS"] - haplotype_maps["capture_start"]) <= tol)
    )
    # aligned_copy means the alignment is on the intended MIP target
    # this is not necessarily the best target, though. For a haplotype sequence
    # to be matched to a MIP target, it also needs to be the best alignment.
    haplotype_maps["mapped_copy"] = (haplotype_maps["aligned_copy"]
                                     & haplotype_maps["best_alignment"])
    # rename some fields to be compatible with previous code
    haplotype_maps.rename(columns={"gene": "Gene", "copy": "Copy",
                                   "chrom": "Chrom"}, inplace=True)
    # any haplotype that does was not best mapped to at least one target
    # will be considered an off target haplotype.
    haplotype_maps["off_target"] = ~haplotype_maps.groupby(
        "haplotype_ID")["mapped_copy"].transform("any")
    off_target_haplotypes = haplotype_maps.loc[haplotype_maps["off_target"]]
    # filter off targets and targets that do not align to haplotypes
    haplotypes = haplotype_maps.loc[(~haplotype_maps["off_target"])
                                    & haplotype_maps["aligned_copy"]]
    # each MIP copy/haplotype_ID combination must have a single alignment
    # if there are multiple, the best one will be chosen

    def get_best_alignment(group):
        return group.sort_values("alignment_score", ascending=False).iloc[0]

    haplotypes = haplotypes.groupby(["MIP", "Copy", "haplotype_ID"],
                                    as_index=False).apply(get_best_alignment)
    haplotypes.index = (range(len(haplotypes)))
    # filter to best mapping copy/haplotype pairs
    mapped_haplotypes = haplotypes.loc[haplotypes["mapped_copy"]]
    mapped_haplotypes["mapped_copy_number"] = mapped_haplotypes.groupby(
        ["haplotype_ID"])["haplotype_ID"].transform(len)

    mapped_haplotypes.to_csv(os.path.join(
        wdir, "mapped_haplotypes.csv"), index=False)
    off_target_haplotypes.to_csv(os.path.join(
        wdir, "offtarget_haplotypes.csv"), index=False)
    haplotypes.to_csv(os.path.join(
        wdir, "aligned_haplotypes.csv"), index=False)
    haplotype_maps.to_csv(os.path.join(
        wdir, "all_haplotypes.csv"), index=False)
    num_hap = len(set(haplotype_maps["haplotype_ID"]))
    num_off = len(set(off_target_haplotypes["haplotype_ID"]))
    print(("{} of {} haplotypes were off-target, either not mapping to "
           "the reference genome, or best mapping to a region which was "
           "not targeted.").format(num_off, num_hap))
    return


def get_vcf_haplotypes(settings):
    """
    Reverse compatibile map_haplotypes function.

    This is the old name for map_haplotypes function. Some notebooks might
    use the old name. So this will just run the map_haplotypes when called
    by the old name.
    """
    map_haplotypes(settings)


def get_haplotype_counts(settings):
    """Get UMI and read counts for each on target haplotype for each sample.

    MIPWrangler output has the UMI and read counts per haplotype but some of
    those are off target and some are mapping to multiple loci by design.
    The decision on whether a haplotype sequence is on or off target and where
    it maps best or if it maps to multiple loci are made by the map_haplotypes
    function. This function distributes the UMI and read counts in the
    MIPWrangler output using the mapped haplotypes data for each sample.
    If a haplotype sequence is uniquely mapping to a targeted locus, we
    allocate all reads for that sample and haplotype sequence to that locus.
    If it is mapping to multiple places, we determine the ratios of those
    'paralogous copies' for that sample based on the average mapping around
    each locus and allocate the reads for that sample and that haplotype
    sequence proportionally to the mapped loci. If a haplotype sequence is
    mapping best to an unintended locus, we remove those.
    """
    wdir = settings["workingDir"]
    ##########################################################
    ##########################################################
    # Process 1: use sample sheet to determine which data points from the
    # mipster file should be used, print relevant statistics.
    ##########################################################
    ##########################################################
    # process sample sheets
    run_meta = pd.read_table(os.path.join(wdir, "samples.tsv"))
    # create a unique sample ID for each sample using sample name,
    # sample set and replicate fields from the sample list file.
    run_meta["sample_name"] = (
            run_meta["sample_name"].astype(str)
        )
    run_meta["Sample Name"] = run_meta["sample_name"]
    run_meta["Sample ID"] = run_meta[
        ["sample_name", "sample_set", "replicate"]
    ].apply(lambda a: "-".join(map(str, a)), axis=1)
    # Sample Set key is reserved for meta data
    # but sometimes erroneously included in the
    # sample sheet. It should be removed.
    try:
        run_meta.drop("Sample Set", inplace=True, axis=1)
    except (ValueError, KeyError):
        pass
    # a change to the formatting of sample sheets uses library_prep
    # instead of Library Prep, so the below line is for backwards compatibility
    run_meta.rename(columns={"library_prep": "Library Prep"}, inplace=True)
    # drop duplicate values originating from
    # multiple sequencing runs of the same libraries
    run_meta = run_meta.drop_duplicates()
    run_meta = run_meta.groupby(
        ["Sample ID", "Library Prep"]
    ).first().reset_index()
    run_meta.to_csv(os.path.join(wdir, "run_meta.csv"))
    # get used sample ids
    sample_ids = run_meta["Sample ID"].unique().tolist()
    ##########################################################
    ##########################################################
    # Process 2: extract all observed variants from observed
    # haplotypes and create a variation data frame that will
    # be able to map haplotype IDs to variation.
    ##########################################################
    ##########################################################
    # get the haplotype dataframe for all mapped haplotypes
    mapped_haplotype_df = pd.read_csv(
        os.path.join(wdir, "mapped_haplotypes.csv"))
    ##########################################################
    ##########################################################
    # Process 3: load the MIPWrangler output which has
    # per sample per haplotype information, such as
    # haplotype sequence quality, barcode counts etc.
    # Create a suitable dataframe that can be merged
    # with variant data to get the same information for each
    # variant (variant barcode count, variant quality, etc.)
    ##########################################################
    ##########################################################
    # get the MIPWrangler Output
    raw_results = pd.read_table(os.path.join(wdir, settings["mipsterFile"]))
    # limit the results to the samples intended for this analysis
    raw_results = raw_results.loc[
        raw_results["sample_name"].isin(sample_ids)
    ]
    # rename some columns for better visualization in tables
    raw_results.rename(
        columns={"sample_name": "Sample ID",
                 "mip_name": "MIP",
                 "gene_name": "Gene",
                 "barcode_count": "Barcode Count",
                 "read_count": "Read Count"},
        inplace=True
    )
    # use only the data corresponding to mapped haplotypes
    # filtering the off target haplotypes.
    mapped_results = raw_results.merge(mapped_haplotype_df, how="inner")
    # Try to estimate the distribution of data that is mapping
    # to multiple places in the genome.
    # This is done in 4 steps.
    # 1) Get uniquely mapping haplotypes and barcode counts
    unique_df = mapped_results.loc[mapped_results["mapped_copy_number"] == 1]
    unique_table = pd.pivot_table(unique_df,
                                  index="Sample ID",
                                  columns=["Gene", "MIP", "Copy", "Chrom"],
                                  values=["Barcode Count"],
                                  aggfunc=np.sum)
    # 2) Estimate the copy number of each paralog gene
    # for each sample from the uniquely mapping data
    # Two values from the settings are used to determine the copy number
    # in a given gene. Average copy count is the ploidy of the organism
    # and the normalization percentile is what percentile is used for
    # normalizing data. For example, for human genes ACC is 2 and
    # if the percentiles are given as 0.4, 0.6: we would calculate the
    # take the 40th and 60th percentile of them barcode counts for each probe
    # across the samples and assume that the average of 40th and 60 pctl values
    # to represent the average copy count of 2. Then caluculate this value
    # for each probe and each sample.
    try:
        average_copy_count = float(settings["averageCopyCount"])
        norm_percentiles = list(map(float,
                                settings["normalizationPercentiles"]))
    except KeyError:
        average_copy_count = 2
        norm_percentiles = [0.4, 0.6]
    unique_df.loc[:, "Copy Average"] = average_copy_count
    # Adjusted barcode count will represent the estimated barcode count
    # for multimapping haplotypes. For example, if hap1 is mapping to 2
    # places in the genome and its barcode count for a sample containing this
    # haplotype is 100. If we determined the copy numbers of the two mapping
    # regions to be 1 and 1, the adjusted barcode count for each region
    # would be 50. We'll set this value for uniquely mapping haplotypes
    # to the Barcode Count, as they are not multi mapping.
    unique_df.loc[:, "Adjusted Barcode Count"] = unique_df["Barcode Count"]
    unique_df.loc[:, "Adjusted Read Count"] = unique_df["Read Count"]
    unique_table.fillna(0, inplace=True)
    # calculate the copy counts using the get_copy_counts function.
    # this function normalizes data for each probe across samples
    # and estimates copy counts using the percentile values as mentioned.
    copy_counts = get_copy_counts(unique_table,
                                  average_copy_count,
                                  norm_percentiles)
    # 3) Estimate the copy number of each "Gene"
    # from the average copy count of uniquely mapping
    # data for all MIPs within the gene.
    cc = copy_counts.groupby(level=["Gene", "Copy"], axis=1).sum()
    gc = copy_counts.groupby(level=["Gene"], axis=1).sum()
    ac = cc.div(gc, level="Gene")
    # 4) Distribute multi mapping data proportional to
    # Paralog's copy number determined from the
    # uniquely mapping data
    multi_df = mapped_results.loc[mapped_results["mapped_copy_number"] > 1]
    if not multi_df.empty:
        # get the average copy count for the gene the haplotype belongs to
        mca = multi_df.apply(lambda r: get_copy_average(r, ac), axis=1)
        multi_df.loc[mca.index, "Copy Average"] = mca
        multi_df["copy_sum"] = multi_df.groupby(
            ["Sample ID", "haplotype_ID"])["Copy Average"].transform("sum")
        multi_df["copy_len"] = multi_df.groupby(
            ["Sample ID", "haplotype_ID"])["Copy Average"].transform("size")
        null_index = multi_df["copy_sum"] == 0
        multi_df.loc[null_index, "Copy Average"] = (
            average_copy_count / multi_df.loc[null_index, "copy_len"])
        multi_df.loc[null_index, "copy_sum"] = average_copy_count
        multi_df["Copy Average"].fillna(0, inplace=True)
        multi_df["Adjusted Barcode Count"] = (multi_df["Barcode Count"]
                                              * multi_df["Copy Average"]
                                              / multi_df["copy_sum"])
        multi_df["Adjusted Read Count"] = (multi_df["Read Count"]
                                           * multi_df["Copy Average"]
                                           / multi_df["copy_sum"])

    # Combine unique and multimapping data
    combined_df = pd.concat([unique_df, multi_df], ignore_index=True,
                            sort=True)
    combined_df.rename(
        columns={
            "Barcode Count": "Raw Barcode Count",
            "Adjusted Barcode Count": "Barcode Count",
            "Read Count": "Raw Read Count",
            "Adjusted Read Count": "Read Count"
        },
        inplace=True
    )
    # print total read and barcode counts
    print(
        (
         "Total number of reads and barcodes were {0[0]} and {0[1]}."
         " On target number of reads and barcodes were {1[0]} and {1[1]}."
        ).format(
            raw_results[["Read Count", "Barcode Count"]].sum(),
            combined_df[["Read Count", "Barcode Count"]].sum().astype(int)
        )
    )
    combined_df.to_csv(os.path.join(wdir, "haplotype_counts.csv"), index=False)
    # So far the count data only includes MIPs that has at least one read
    # in at least one sample. We would like to include MIPs with no reads
    # as well. So we'll create a dataframe that has all the intended MIPs
    # and merge with the count data.
    # create call_info data frame for all used probes in the experiment
    call_info_file = settings["callInfoDictionary"]
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    call_df_list = []
    for g in call_info:
        for m in call_info[g]:
            if m in used_probes:
                for c in call_info[g][m]["copies"]:
                    call_dict = {"MIP": m, "Copy": c}
                    call_df_list.append(pd.DataFrame(call_dict, index=[0]))
    call_df = pd.concat(call_df_list, ignore_index=True, sort=True)

    # merge the count data with probe data. Fill missing values with 0.
    combined_df = call_df.merge(combined_df, how="left").fillna(0)
    # Create pivot table of combined barcode counts
    # This is a per MIP per sample barcode count table
    # of the samples with sequencing data
    barcode_counts = pd.pivot_table(combined_df,
                                    index="Sample ID",
                                    columns=["MIP",
                                             "Copy"],
                                    values=["Barcode Count"],
                                    aggfunc=np.sum)
    # Sample name for probes without data would be NA and replaced to 0
    # remove that if it exists
    try:
        barcode_counts.drop(0, inplace=True)
    except KeyError:
        pass
    print("There are {} samples with sequence data".format(
        barcode_counts.shape[0]
    ))
    # After pivot table is created, the column names have an extra
    # row with the name "Barcode Count". Remove that from column names.
    bc_cols = barcode_counts.columns
    bc_cols = [bc[1:] for bc in bc_cols]
    # barcode count data is only available for samples with data
    # so if a sample has not produced any data, it will be missing
    # these samples should be added with 0 values for each probe
    all_barcode_counts = pd.merge(
        run_meta[["Sample ID", "replicate"]].set_index("Sample ID"),
        barcode_counts, left_index=True, right_index=True, how="left")
    all_barcode_counts.drop("replicate", axis=1, inplace=True)
    # fix column names
    all_barcode_counts.columns = pd.MultiIndex.from_tuples(
        bc_cols, names=["MIP", "Copy"]
    )
    all_barcode_counts.fillna(0, inplace=True)
    print("There are {} total samples.".format(all_barcode_counts.shape[0]))
    all_barcode_counts.to_csv(os.path.join(wdir, "barcode_counts.csv"))
    # Create an overview statistics file for samples including
    # total read count, barcode count, and how well they cover each MIP.
    sample_counts = combined_df.groupby("Sample ID")[["Read Count",
                                                      "Barcode Count"]].sum()
    # Find samples without any data and print the number
    no_data = run_meta.loc[
        ~run_meta["Sample ID"].isin(sample_counts.index)
    ]
    print(("{} out of {} samples had no data and they will be excluded from "
           "the variant calls.").format(no_data.shape[0], run_meta.shape[0]))

    # add samples with no data
    sample_counts = pd.merge(
        run_meta[["Sample ID", "replicate"]].set_index("Sample ID"),
        sample_counts, left_index=True, right_index=True, how="left")
    sample_counts.drop("replicate", axis=1, inplace=True)
    target_cov = pd.concat(
        [(all_barcode_counts >= 1).sum(axis=1),
         (all_barcode_counts >= 5).sum(axis=1),
         (all_barcode_counts >= 10).sum(axis=1)],
        axis=1,
    ).rename(
        columns={
            0: "targets_with_1_barcodes",
            1: "targets_with_5_barcodes",
            2: "targets_with_10_barcodes"
        }
    )
    sample_counts = sample_counts.merge(target_cov,
                                        how="outer",
                                        left_index=True,
                                        right_index=True).fillna(0)
    target_cov_file = os.path.join(wdir, "sample_summary.csv")
    sample_counts.to_csv(target_cov_file)

    return


def freebayes_call(bam_dir="/opt/analysis/padded_bams",
                   fastq_dir="/opt/analysis/padded_fastqs",
                   options=[],
                   vcf_file="/opt/analysis/variants.vcf.gz",
                   targets_file=None, make_fastq=True,
                   align=True, settings=None, settings_file=None,
                   bam_files=None, bam_list=None, verbose=True,
                   fastq_padding=20, min_base_quality=1,
                   errors_file="/opt/analysis/freebayes_errors.txt",
                   warnings_file="/opt/analysis/freebayes_warnings.txt",
                   merge_distance=1000, contig_padding=500):
    """Call variants for MIP data using freebayes.

    A mapped haplotype file must be present in the working directory. This
    is generated during haplotype processing. Per sample fastqs and bams
    will be created if align=True. Fastqs are generated with a default 20 bp
    padding on each side of the haplotype. This assumes that there were no
    errors where the MIP arms bind to the DNA. It may cause some false negative
    calls where there was imperfect binding, but it is crucial for determining
    variants close to the MIP arms.

    Parameters
    ----------
    bam_dir: str/path, /opt/analysis/padded_bams
        path to the directory where per sample bam files are or where they
        will be created if align=True.
    fastq_dir: str/path, /opt/analysis/padded_fastqs
        path to the directory where per sample fastq files are or  where they
        will be created if align=True.
    vcf_file: str/path, /opt/analysis/variants.vcf.gz
        Output vcf file path.
    options: list, []
        options to pass to freebayes directly, such as --min-coverage
        the list must have each parameter and value as separate items.
        For example, ["--min-alternate-count", "2"] and not
        ["--min-alternate-count 2"]
    align: bool, True
        Set to false if fastq and bam files have already been created.
    settings: dict, None
        Analysis settings dictionary. Either this or settings_file must
        be provided.
    settings_file: str/path, None
        Path to the analysis settings file. Either this or the settings dict
        must be provided.
    targets_file: str/path, None
        Path to targets file to force calls on certain locations even if
        those variants do not satisfy filter criteria. It must be a tab
        separated text file with minimum columns CHROM, POS, REF, ALT.
    bam_files: list, None
        list of bam files within the bam_dir to pass to freebayes. If None (
        default), all bam files in the bam_dir will be used.
    verbose: bool, True
        if set to True, print errors and warnings in addition to saving to
        errors and warnings files.
    errors_file: str/path, /opt/analysis/freebayes_errors.txt
        file to save freebayes errors.
    warnings_file: str/path, /opt/analysis/freebayes_warnings
        file to save freebayes warnings
    merge_distance: int, 200
        When creating contigs from MIP target regions, merge targets closer
        to each other than this distance.
    contig_padding: int, 50
        Add this much padding to the contigs when calling freebayes.
    """
    # get the analysis settings
    # check if both settings and the settings file are None:
    if (settings is None) and (settings_file is None):
        print("settings or settings file must be provided for freebayes_call.")
        return
    else:
        if settings is None:
            settings = get_analysis_settings(settings_file)
        else:
            settings = copy.deepcopy(settings)
    # get the working directory from settings
    wdir = settings["workingDir"]
    # load mapped haplotypes file. This file has the genomic locations
    # of the haplotypes in mip data
    mapped_haplotypes_file = os.path.join(wdir, "mapped_haplotypes.csv")
    # get the mip data file location. This file has per sample haplotype
    # information including counts.
    mipster_file = os.path.join(wdir, settings["mipsterFile"])
    if make_fastq:
        # create fastq files from MIP data. One read per UMI will be created.
        generate_mapped_fastqs(fastq_dir, mipster_file,
                               mapped_haplotypes_file, settings["species"],
                               pro=int(settings["processorNumber"]),
                               pad_size=fastq_padding)
    if align:
        # map per sample fastqs to the reference genome, creating bam files.
        # bam files will have sample groups added, which is required for
        # calling variants across the samples.
        bwa_multi([], "bam", fastq_dir, bam_dir,
                  settings["bwaOptions"], settings["species"],
                  int(settings["processorNumber"]),
                  int(settings["processorNumber"]))

    # divide data into contigs to make parallelization more efficient
    # we'll create contigs from overlapping MIPs.
    # load the call info dictionary which contains per MIP information
    call_file = settings["callInfoDictionary"]
    with open(call_file) as infile:
        call_dict = json.load(infile)
    # create a dataframe that has the genomic coordinates of each MIP
    call_df = []
    for g in call_dict:
        for m in call_dict[g]:
            for c in call_dict[g][m]["copies"]:
                cdict = call_dict[g][m]["copies"][c]
                call_df.append([cdict["chrom"], cdict["capture_start"],
                                cdict["capture_end"]])
    call_df = pd.DataFrame(call_df, columns=["chrom", "capture_start",
                                             "capture_end"])

    # create a function that generates contigs of MIPs which overlap
    # with 1 kb padding on both sides.
    def get_contig(g):
        intervals = zip(g["capture_start"], g["capture_end"])
        return pd.DataFrame(merge_overlap(
            [list(i) for i in intervals], spacer=merge_distance))

    # create contigs per chromosome
    contigs = call_df.groupby("chrom").apply(get_contig)
    contigs = contigs.reset_index()
    contigs.rename(columns={"level_1": "contig", 0: "contig_capture_start",
                            1: "contig_capture_end"}, inplace=True)
    contigs["contig_name"] = contigs["chrom"] + "_" + contigs["contig"].astype(
        str)
    # we'll call freebayes on each contig by providing a region string in the
    # form chrx:begin-end. Create those strings for each contig with some
    # padding. It is important to check that we don't end up with a start
    # position of <1 or end position longer than chom length.
    # Begin by adding chromosome length to contig info.
    # get reference chromosome lengths
    genome_file = get_file_locations()[settings["species"]]["fasta_genome"]
    reference_lengths = {}
    genome_sam = pysam.FastaFile(genome_file)
    for r in genome_sam.references:
        reference_lengths[r] = genome_sam.get_reference_length(r)
    contigs["chromosome_length"] = contigs["chrom"].map(reference_lengths)
    contigs["region_start"] = contigs["contig_capture_start"] - contig_padding
    contigs.loc[contigs["region_start"] < 1, "region_start"] = 1
    contigs["region_end"] = contigs["contig_capture_end"] + contig_padding
    contigs["region_end"] = contigs[
        ["region_end", "chromosome_length"]].min(axis=1).values
    contigs["region"] = contigs["chrom"] + ":" + (
        contigs["region_start"]).astype(str) + "-" + (
        contigs["region_end"]).astype(str)
    # we'll force calls on targeted variants if so specified
    if targets_file is not None:
        # each contig must include at least one of the targets, otherwise
        # freebayes throws an error. So we'll load the targets and add the
        # targets option to only those contigs that contain targets
        targets = pd.read_table(targets_file)
        # merge targets and contigs dataframes to determine which contigs
        # contain targets. chrom will be used as the common column name
        targets["chrom"] = targets["CHROM"]
        targets = targets.merge(contigs)
        # remove rows where chrom is shared but target position is outside
        # of contig boundries.
        targets = targets.loc[
            (targets["contig_capture_start"] <= targets["POS"])
            & (targets["POS"] <= targets["contig_capture_end"])]
        targets["contains_targets"] = True
        # merge only two columns of the targets df to contigs so that
        # the only shared column is contig_name. More than one target can
        # be in a single contig, so we need to drop duplicates from targets.
        contigs = contigs.merge(targets[
            ["contig_name", "contains_targets"]].drop_duplicates(), how="left")
        contigs["contains_targets"].fillna(False, inplace=True)
        # create a targets.vcf file for freebayes
        targets_vcf = os.path.join(wdir, "targets.vcf")
        with open(targets_vcf, "w") as outfile:
            outfile.write('##fileformat=VCFv4.2\n')
            outfile.write(
                '##FILTER=<ID=PASS,Description="All filters passed">\n')
            outfile.write('##INFO=<ID=TR,Number=.,Type=String,Description'
                          '="Targeted variant.">\n')
            vcf_fields = ["ID", "QUAL", "FILTER"]
            for vf in vcf_fields:
                targets[vf] = "."
            targets["INFO"] = "TR"
            vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                          "FILTER", "INFO"]
            targets = targets.rename(columns={"CHROM": "#CHROM"})[vcf_fields]
            targets.sort_values(["#CHROM", "POS"]).to_csv(
                outfile, sep="\t", index=False)
        # bgzip and index
        res = subprocess.run(["bgzip", "-f", targets_vcf],
                             stderr=subprocess.PIPE)
        if res.returncode != 0:
            print("Error in compressing targets.vcf file", res.stderr)
        targets_vcf = targets_vcf + ".gz"
        res = subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-f",
                              targets_vcf], stderr=subprocess.PIPE)
        if res.returncode != 0:
            print("Error in indexing targets.vcf.gz file ", res.stderr)
    else:
        contigs["contains_targets"] = False

    # create a contig dictionary from the contigs dataframe
    # this dict will be passed to the worker function for parallelization
    chrom_dict = {}
    gb = contigs.groupby("chrom")
    for g in gb.groups.keys():
        gr = gb.get_group(g)
        chrom_dict[g] = gr[["contig_name", "region",
                            "contains_targets"]].set_index(
            "contig_name").to_dict(orient="index")

    # populate the contigs dictionary for freebayes parameters
    # start with options to be added for each contig
    # get fasta genome location
    genome_fasta = get_file_locations()[settings["species"]]["fasta_genome"]
    # specify fasta genome file
    options.extend(["-f", genome_fasta])
    # add if bam files are specified. Nothing should be added to options
    # after the bam files.
    if bam_files is not None:
        options.extend(bam_files)
    if bam_list is not None:
        options.extend(["-L", bam_list])
    # create a file list in the bam_dir that has full path to all bam files
    # if all bam files are to be used
    else:
        bam_list = os.path.join(bam_dir, "bamlist.txt")
        with open(bam_list, "w") as outfile:
            for f in os.scandir(bam_dir):
                if os.path.splitext(f.name)[1] == ".bam":
                    outfile.write(f.path + "\n")
        options.extend(["-L", bam_list])
    # add minimum base quality parameter to options if not already provided
    if ("--min-base-quality" not in options) and ("-q" not in options):
        options.extend(["-q", str(min_base_quality)])
    # create a list for keeping all contig vcf file paths to concatanate
    # them at the end.
    contig_vcf_paths = []
    # create a similar list for zipped vcf files
    contig_vcf_gz_paths = []
    # create a list of per contig dictionary to feed to multiprocessing
    # function apply_async
    contig_dict_list = []
    # create the contigs vcf directory
    cvcfs_dir = os.path.join(wdir, "contig_vcfs")
    if not os.path.exists(cvcfs_dir):
        os.makedirs(cvcfs_dir)
    # update contig_dict with contig specific options
    for chrom in chrom_dict:
        for contig_name in chrom_dict[chrom]:
            contig_dict = chrom_dict[chrom][contig_name]
            ################################################################
            # create contig specific options and
            # add contigs region string (chrx:begin-end)
            region = contig_dict["region"]
            contig_options = ["-r", region]
            # add contigs vcf file name
            contig_vcf = os.path.join(wdir, "contig_vcfs",
                                      contig_name + ".vcf")
            contig_dict["vcf_path"] = contig_vcf
            # add output file to the freebayes options
            contig_options.extend(["-v", contig_vcf])
            # add contig vcf path to the list
            contig_vcf_paths.append(contig_vcf)
            # add contigs vcf.gz file name
            contig_vcf_gz = os.path.join(wdir, "contig_vcfs",
                                         contig_name + ".vcf.gz")
            contig_vcf_gz_paths.append(contig_vcf_gz)
            contig_dict["vcf_gz_path"] = contig_vcf_gz
            # if contig includes targets, we'll force calls on those
            if contig_dict["contains_targets"]:
                contig_options.extend(["-@", targets_vcf])
            # we'll add the contig specific options to the beginning of
            # the options list in case bam files were added to the options
            # and they must stay at the end because they are positional args.
            contig_dict["options"] = contig_options + options
            # add the contig dict to contig dict list
            contig_dict_list.append(contig_dict)

    # create a processor pool for parallel processing
    pool = Pool(int(settings["processorNumber"]))
    # create a results container for the return values from the worker function
    results = []
    errors = []
    # run the freebayes worker program in parallel
    pool.map_async(freebayes_worker, contig_dict_list, callback=results.extend,
                   error_callback=errors.extend)
    # join and close the processor pool.
    pool.close()
    pool.join()

    # compare the length of the results object and the number of contigs
    # print an error message if they are not the same
    if len(contig_dict_list) != (len(results) + len(errors)):
        print(("Number of contigs, {}, is not the same as number of results "
               "from the variant caller, {}, plus number of errors, {}. "
               "This means some calls have failed silently. "
               "Results and errors should be inspected.").format(
               len(contig_dict_list), len(results), len(errors)))
    # check each contig's variant call results for errors and warnings
    # open files to save errors and warnings
    with open(errors_file, "w") as ef, open(warnings_file, "wb") as wf:
        # keep a count of warnings an errors
        error_count = 0
        warning_count = 0
        for res in results:
            for r in res:
                try:
                    r.check_returncode()
                except subprocess.CalledProcessError as e:
                    error_count += 1
                    ef.write(str(e) + "\n")
                    if verbose:
                        print("Error in freebayes calls: ", e)
                # print if any warnings were issued
                if len(r.stderr) > 0:
                    warning_count += 1
                    wf.write(r.stderr + b"\n")
                    if verbose:
                        print("Warning in freebayes calls: ", r.stderr)
        # if errors are not printed but present, print an message to indicate
        # the presence of errors/warnings
        if not verbose:
            if error_count > 0:
                print(("Errors were encountered in freebayes calls."
                       " Please inspect {} for errors.").format(errors_file))
            if warning_count > 0:
                print(("There were warnings in freebayes calls."
                       " Please inspect {} for warnings.").format(
                           warnings_file))

    if len(errors) > 0:
        print(("There were {} calls that failed").format(len(errors)))

    # concatanate contig vcfs. The number of contigs may be high, so we'll
    # write the vcf paths to a file and bcftools will read from that file
    cvcf_paths_file = os.path.join(wdir, "contig_vcfs", "vcf_file_list.txt")
    with open(cvcf_paths_file, "w") as outfile:
        outfile.write("\n".join(contig_vcf_gz_paths) + "\n")
    subprocess.run(["bcftools", "concat", "-f", cvcf_paths_file, "-Oz",
                    "-o", vcf_file], check=True)
    subprocess.run(["bcftools", "index", "-f", vcf_file], check=True)

    # fix vcf header if --gvcf option has been used
    if "--gvcf" in options:
        temp_vcf_path = os.path.join(wdir, "temp.vcf.gz")
        vcf_reheader(os.path.basename(vcf_file), temp_vcf_path, wdir=wdir)
        old_vcf_path = os.path.join(wdir, "unfixed.vcf.gz")
        subprocess.run(["mv", vcf_file, old_vcf_path])
        subprocess.run(["mv", temp_vcf_path, vcf_file])
        subprocess.run(["bcftools", "index", "-f", vcf_file], check=True)
    return (contig_dict_list, results, errors)


def freebayes_worker(contig_dict):
    """Run freebayes program with the specified options.

    Run freebayes program with the specified options and return a
    subprocess.CompletedProcess object.
    """
    options = contig_dict["options"]
    command = ["freebayes"]
    command.extend(options)
    # run freebayes command piping the output
    fres = subprocess.run(command, stderr=subprocess.PIPE)
    # check the return code of the freebayes run. if succesfull continue
    if fres.returncode == 0:
        # bgzip the vcf output, using the freebayes output as bgzip input
        vcf_path = contig_dict["vcf_path"]
        gres = subprocess.run(["bgzip", "-f", vcf_path],
                              stderr=subprocess.PIPE)
        # make sure bugzip process completed successfully
        if gres.returncode == 0:
            # index the vcf.gz file
            vcf_gz_path = contig_dict["vcf_gz_path"]
            ires = subprocess.run(["bcftools", "index", "-f", vcf_gz_path],
                                  stderr=subprocess.PIPE)
            # return the CompletedProcess objects
            return (fres, gres, ires)
        else:
            return (fres, gres)
    # if freebayes call failed, return the completed process object
    # instead of attempting to zip the vcf file which does not exist if
    # freebayes failed.
    else:
        return (fres, )


def vcf_reheader(vcf_file, fixed_vcf_file, wdir="/opt/analysis/"):
    """Fix vcf header QA/QR fields.

    When --gvcf option is used in freebayes variant calling pipeline,
    the header of the vcf file comes out incorrect for QA/QR fields number
    type, Integer instead of Float. This function fixes those lines from
    the header and creates a new vcf file with the correct header.
    """
    # get the current header
    vcf_path = os.path.join(wdir, vcf_file)
    header = subprocess.Popen(["bcftools", "view", "-h", vcf_path],
                              stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    com = header.communicate()
    if header.returncode != 0:
        print("Failed to extract vcf header. Header will not be fixed.")
        return

    # convert the header byte string to text and creale a list of lines
    head = com[0].decode("utf-8").split("\n")
    # create a new header with fixed fields
    new_head = []
    for line in head:
        if ("ID=QA" in line) or ("ID=QR" in line):
            new_head.append(line.replace("Type=Integer", "Type=Float"))
        else:
            new_head.append(line)

    new_header_path = os.path.join(wdir, "new_vcf_header.txt")
    with open(new_header_path, "w") as outfile:
        outfile.write("\n".join(new_head) + "\n")

    fixed_vcf_path = os.path.join(wdir, fixed_vcf_file)
    subprocess.run(["bcftools", "reheader", "-h", new_header_path,
                    vcf_path,  "-o", fixed_vcf_path], check=True)
    return


def gatk(options):
    """GATK wrapper function.

    Run gatk program with the given options. Return the subprocess result.
    """
    return subprocess.run(["gatk", *options], stderr=subprocess.PIPE)


def gatk_file_prep(bam_dir="/opt/analysis/padded_bams",
                   fastq_dir="/opt/analysis/padded_fastqs",
                   targets_file=None,
                   settings=None, settings_file=None,
                   errors_file="/opt/analysis/gatk_file_prep_output.txt"):
    """Prepare files for calling variants for MIP data using gatk.

    A mapped haplotype file must be present in the working directory. This
    is generated during haplotype processing. Per sample fastqs and bams
    will be created. Fastqs are generated with a default 20 bp
    padding on each side of the haplotype. This assumes that there were no
    errors where the MIP arms bind to the DNA. It may cause some false negative
    calls where there was imperfect binding, but it is crucial for determining
    variants close to the MIP arms.

    Parameters
    ----------
    bam_dir: str/path, /opt/analysis/padded_bams
        path to the directory where per sample bam files are or where they
        will be created if align=True.
    fastq_dir: str/path, /opt/analysis/padded_fastqs
        path to the directory where per sample fastq files are or  where they
        will be created if align=True.
    settings: dict, None
        Analysis settings dictionary. Either this or settings_file must
        be provided.
    settings_file: str/path, None
        Path to the analysis settings file. Either this or the settings dict
        must be provided.
    targets_file: str/path, None
        Path to targets file to force calls on certain locations even if
        those variants do not satisfy filter criteria. It must be a tab
        separated text file with minimum columns CHROM, POS, REF, ALT.
    errors_file: str/path, /opt/analysis/gatk_file_prep_errors.txt
        file to save freebayes errors.
    """
    # get the analysis settings
    # check if both settings and the settings file are None:
    if (settings is None) and (settings_file is None):
        print("settings or settings file must be provided for freebayes_call.")
        return
    else:
        if settings is None:
            settings = get_analysis_settings(settings_file)
        else:
            settings = copy.deepcopy(settings)
    # get the working directory from settings
    wdir = settings["workingDir"]
    # load mapped haplotypes file. This file has the genomic locations
    # of the haplotypes in mip data
    mapped_haplotypes_file = os.path.join(wdir, "mapped_haplotypes.csv")
    # get the mip data file location. This file has per sample haplotype
    # information including counts.
    mipster_file = os.path.join(wdir, settings["mipsterFile"])
    # get the mip data file location. This file has per sample haplotype
    # information including counts.
    mipster_file = os.path.join(wdir, settings["mipsterFile"])
    # create fastq files from MIP data. One read per UMI will be created.
    generate_mapped_fastqs(fastq_dir, mipster_file,
                           mapped_haplotypes_file, settings["species"],
                           pro=int(settings["processorNumber"]))
    # if there is a targets file provided, we'll create a hypothetical
    # sample that has all of the targeted variants. This way, a variant site
    # for each target will be created in the final vcf file even if a
    # variant was not observed in the data.
    if targets_file is not None:
        # load the targets as dataframe converting field names to
        # field names in a haplotypes file.
        targets = pd.read_table(targets_file).rename(
            columns={"CHROM": "Chrom", "POS": "capture_start",
                     "ALT": "haplotype_sequence",
                     "mutation_name": "haplotype_ID"})
        # fill in orientation and copy number information for all targets.
        targets["orientation"] = "forward"
        targets["mapped_copy_number"] = 1
        targets["capture_end"] = (targets["capture_start"]
                                  + targets["REF"].apply(len) - 1)
        # create a haplotype file for the targeted mutations
        haplotype_fields = ['capture_end', 'capture_start', 'Chrom',
                            'orientation', 'haplotype_ID',
                            'haplotype_sequence', 'mapped_copy_number']
        mutant_haplotypes = "/opt/analysis/mutant_haplotypes.csv"
        targets[haplotype_fields].to_csv(mutant_haplotypes, index=False)

        # create a hypothetical sample that has all mutations and a
        # corresponding mip data file that shows a UMI count of 20
        # for each observation
        targets["sample_name"] = "control_mutant"
        targets["sequence_quality"] = targets["haplotype_sequence"].apply(
            lambda a: "".join(["H" for i in range(len(a))]))
        targets["barcode_count"] = 20
        data_fields = ["sample_name", 'haplotype_ID', "haplotype_sequence",
                       'sequence_quality', 'barcode_count']
        mutant_data_file = "/opt/analysis/mutant_data.tsv"
        targets[data_fields].to_csv(mutant_data_file, index=False, sep="\t")
        # create a fastq file for the "control_mutant" sample
        padding = 100
        generate_mapped_fastqs(fastq_dir, mutant_data_file,
                               mutant_haplotypes, settings["species"],
                               pro=int(settings["processorNumber"]),
                               pad_size=padding)
    # map per sample fastqs to the reference genome, creating bam files.
    # bam files will have sample groups added, which is required for
    # calling variants across the samples.
    bwa_multi([], "bam", fastq_dir, bam_dir,
              settings["bwaOptions"], settings["species"],
              int(settings["processorNumber"]),
              int(settings["processorNumber"]))

    # create an  intervals file to be used in gatk call
    intervals_bed = "/opt/analysis/intervals.bed"
    call_file = settings["callInfoDictionary"]
    with open(call_file) as infile:
        call_dict = json.load(infile)
    # create a dataframe that has the genomic coordinates of each MIP
    probe_info = []
    for g in call_dict:
        for m in call_dict[g]:
            for c in call_dict[g][m]["copies"]:
                cdict = call_dict[g][m]["copies"][c]
                probe_info.append([cdict["chrom"], cdict["capture_start"],
                                   cdict["capture_end"]])
    probe_info = pd.DataFrame(probe_info, columns=["chrom", "capture_start",
                                                   "capture_end"])
    probe_info["bed_start"] = probe_info["capture_start"] - 200
    probe_info["bed_end"] = probe_info["capture_end"] + 200
    probe_info[["chrom", "bed_start", "bed_end"]].to_csv(
        intervals_bed, index=False, header=(None), sep="\t")
    intervals_list = "/opt/analysis/intervals.list"
    genome_dict = get_file_locations()[settings["species"]]["genome_dict"]
    interval_call = gatk(["BedToIntervalList", "-I", intervals_bed,
                          "-O", intervals_list, "-SD", genome_dict])
    # check the return code and if not 0 print warning
    if interval_call.returncode != 0:
        print(("An error ocurred when creating the intervals list. "
               "Please see the {} for details.").format(errors_file))
    # save command output
    with open(errors_file, "ab") as outfile:
        outfile.write(interval_call.stderr)


def gatk_haplotype_caller(
        options, bam_dir, settings,
        errors_file="/opt/analysis/gatk_haplotype_caller_output.txt"):
    genome_fasta = get_file_locations()[settings["species"]]["fasta_genome"]
    intervals_list = "/opt/analysis/intervals.list"
    haplotype_caller_opts = ["HaplotypeCaller", "-R", genome_fasta,
                             "--native-pair-hmm-threads", "1",
                             "-L", intervals_list] + options
    # scan the bam directory and get file paths. Assign an output name
    # for each file (gvcf output)
    bam_files = []
    for f in os.scandir(bam_dir):
        if os.path.splitext(f.name)[1] == ".bam":
            base_name = os.path.splitext(f.name)[0]
            gvcf = os.path.join(bam_dir, base_name + ".g.vcf.gz")
            bam_files.append([f.path, gvcf])
    pool = NoDaemonProcessPool(int(settings["processorNumber"]))
    results = []
    errors = []
    for bam in bam_files:
        io_options = ["-I", bam[0], "-O", bam[1]]
        pool.apply_async(gatk, (haplotype_caller_opts + io_options, ),
                         callback=results.append, error_callback=errors.append)
    pool.close()
    pool.join()
    if len(errors) > 0:
        print(("An error ocurred during haplotype calling . "
               "Please see the {} for details.").format(errors_file))
        # save command output
        with open(errors_file, "ab") as outfile:
            for e in errors:
                outfile.write(str(e))
    for r in results:
        if r.returncode != 0:
            print(("An error ocurred when creating the intervals list. "
                   "Please see the {} for details.").format(errors_file))
        # save command output
        with open(errors_file, "ab") as outfile:
            outfile.write(r.stderr)
    return


def genotype_gvcfs(settings, bam_dir, options, gdb, vcf_file,
                   sample_map=None, keep_control_mutant=False,
                   errors_file="/opt/analysis/gatk_genotype_gvcfs_output.txt"):
    if sample_map is None:
        # scan the bam directory and get file paths. Assign an output name
        # for each file (gvcf output)
        bam_files = []
        for f in os.scandir(bam_dir):
            if os.path.splitext(f.name)[1] == ".bam":
                base_name = os.path.splitext(f.name)[0]
                gvcf = os.path.join(bam_dir, base_name + ".g.vcf.gz")
                bam_files.append([f.path, gvcf])
        sample_map = os.path.join(settings["workingDir"], "sample_map.txt")
        with open(sample_map, "w") as outfile:
            for f in bam_files:
                sample_name = ".".join(os.path.basename(f[0]).split(".")[:-2])
                outfile.write(sample_name + "\t" + f[1] + "\n")
    intervals_list = "/opt/analysis/intervals.list"
    gdb_path = os.path.join("/opt/analysis/", gdb)
    gdb_import = ["--java-options", "-Xmx32G", "GenomicsDBImport",
                  "--genomicsdb-workspace-path", gdb_path,
                  "--sample-name-map", sample_map,
                  "-L", intervals_list,
                  "--max-num-intervals-to-import-in-parallel",
                  settings["processorNumber"]]
    gdb_result = gatk(gdb_import)
    if gdb_result.returncode != 0:
        print(("An error ocurred when during genomics DB import. "
               "Please see the {} for details.").format(errors_file))
    # save command output
    with open(errors_file, "ab") as outfile:
        outfile.write(gdb_result.stderr)

    # genotype gvcfs
    genome_fasta = get_file_locations()[settings["species"]][
        "fasta_genome"]
    gdb = "gendb://" + gdb
    if keep_control_mutant:
        temp_vcf_file = vcf_file
    else:
        temp_vcf_file = "/opt/analysis/temp.vcf.gz"
    genotype_gvcfs = ["GenotypeGVCFs", "-R", genome_fasta,
                      "-V", gdb, "-O", temp_vcf_file, "-L", intervals_list]
    genotypes = gatk(genotype_gvcfs + options)
    if genotypes.returncode != 0:
        print(("An error ocurred during genotyping GVCFs. "
               "Please see the {} for details.").format(errors_file))
    # save command output
    with open(errors_file, "ab") as outfile:
        outfile.write(genotypes.stderr)

    # remove control mutant sample if requested
    if not keep_control_mutant:
        res = subprocess.run(["bcftools", "view", "-s^control_mutant",
                              "-Oz", "-o", vcf_file, temp_vcf_file,
                              "--force-samples"],
                             stderr=subprocess.PIPE)
        if res.returncode != 0:
            print(("An error ocurred while removing control mutant. "
                   "Please see the {} for details.").format(errors_file))
        # save command output
        with open(errors_file, "ab") as outfile:
            outfile.write(res.stderr)
        # index the final vcf file
        res = subprocess.run(["bcftools", "index", "-f", vcf_file],
                             stderr=subprocess.PIPE)
        if res.returncode != 0:
            print(("An error ocurred while indexing the final vcf file. "
                   "Please see the {} for details.").format(errors_file))
        # save command output
        with open(errors_file, "ab") as outfile:
            outfile.write(res.stderr)


def vcf_to_tables_fb(vcf_file, settings=None, settings_file=None,
                     annotate=True, geneid_to_genename=None,
                     target_aa_annotation=None, aggregate_aminoacids=False,
                     target_nt_annotation=None, aggregate_nucleotides=False,
                     decompose_options=[], annotated_vcf=False,
                     aggregate_none=False, min_site_qual=-1,
                     min_target_site_qual=-1, min_genotype_qual=-1,
                     min_alt_qual=-1, min_ref_qual=-1, min_mean_alt_qual=-1,
                     min_mean_ref_qual=-1, output_prefix=""):
    """Create various tables from a vcf file.

    Create various tables from a vcf file generated by the freebayes
    program. There are 3 different types of count output for each variant:
    variant count, reference count and coverage. The vcf file will be split
    into biallelic variants. Table versions of the input vcf will be created
    but the info fields will be limited to the mandatory vcf fields and some
    annotation data if avaliable.

    In addition to the original vcf table, aa change tables can be generated.
    These will be generated by filtering the vcf to missense variants only,
    decomposing block substitutions (haplotypes) and combining the counts for
    the same aminoacid changes. This operation is specifically intended for
    generating data for targeted missense mutations and only reports that. All
    other variants, even those complex variants including targeted variants
    will not be reported. Finally, one specific mutation (dhps-437) will have
    reference counts instead of variant counts if present. This is because this
    drug resistance variant is encoded by the 3d7 reference sequence.

    Parameters
    ----------
    settings: dict, None
        Analysis settings dictionary. Either this or settings_file must
        be provided.
    settings_file: str/path, None
        Path to the analysis settings file. Either this or the settings dict
        must be provided.
    annotate: bool, True
        Annotate variant file. This is required for protein level analysis.
    vcf_file: str/path
        Starting vcf file.
    geneid2genename: str/path, None.
        Path to a tab separated tex file that maps gene ids to gene names.
        Column names must be gene_id and gene_name. Gene IDs
        will populate the Gene field if this file is not provided.
    target_aa_annotation: str/path, None.
        Path to a tab separated text file with targeted variant information to
        annotate and label targeted amino acid changes.
        It must have gene_name, aminoacid_change, and mutation_name columns.
        Amino acid changes should be represented as refAAPosAltAA. refAA and
        AltAA must be three letter amino acid codes.
        This file is required for targeted protein variant labeling.
    target_nt_annotation: str/path, None.
        Path to a tab separated text file with targeted variant information to
        annotate and label targeted nucleotide changes.
        It must have CHROM, POS, REF, ALT, NAME columns.
        This file is required for targeted nucleotide variant labeling.
    aggregate_aminoacids: bool, False
        whether counts for same amino acids should be aggregated. This involves
        decomposing multi amino acid changes for missense variants. If amino
        acid based targets will be annotated, based on a provided annotation
        dictionary, aggregation step must be completed. Targeted mutations
        that are part of complex events (indels, stop loss/gain etc.) will not
        be labeled as targeted.
    aggregate_nucleotides: bool, False
        whether the counts for nucleotide changes should be aggregated. This
        involves decomposing all variants to the smallest units possible,
        breaking all haplotype data. The level of decomposition should be
        specified with the decompose_options parameter.
    aggregate_none: bool, False.
        Do no aggregation on counts, save the original (annotated if requested)
        vcf file as 3  count tables. Three aggregation options are compatible
        with each other and can be used all at once.
    decompose_options: list, []
        if aggregate nucleotides option is selected, these options will be
        passed to vt program. "-a" for decomposing variants containing indels,
        for example. "-p" for keeping phase information. Any option to vt
        decompose_blocksub would be valid. By default indels will not be
        decomposed.
    annotated_vcf: bool, False
        is the provided vcf file annotated using snpEff. These annotations
        will be used if no count aggregation is to be done and annotate option
        is False.
    min_site_qual: float, -1
        Filter variants with QUAL values less than this value if the site is
        not a targeted site. If targeted, the site will be kept regardless of
        the qual value for the site. freebayes manual indicates that
        simulations showed a value between 1-30 would be good. So a minimum
        value of 1 here would clean up most junk sites.
    min_target_site_qual: float, -1
        If a variant site is targeted but the site qual is lower than this,
        reset the alternate observation counts to 0. It may be best to leave
        this at the default value since there is usually additional evidence
        that a targeted variant exists in a samples compared to a de novo
        variant.
    """
    # get the analysis settings
    # check if both settings and the settings file are None:
    if (settings is None) and (settings_file is None):
        print("settings or settings file must be provided for freebayes_call.")
        return
    else:
        if settings is None:
            settings = get_analysis_settings(settings_file)
        else:
            settings = copy.deepcopy(settings)
    # get the working directory from settings
    wdir = settings["workingDir"]
    # All postprocessing steps require biallelic variant representation.
    # so we'll use bcftools to split multiallelics to their own lines.
    genome_fasta = get_file_locations()[settings["species"]]["fasta_genome"]
    vcf_path = os.path.join(wdir, vcf_file)
    split_vcf_path = os.path.join(wdir, output_prefix + "split." + vcf_file)
    subprocess.run(["bcftools", "norm", "-f", genome_fasta, "-m-both",
                    vcf_path, "-Oz", "-o", split_vcf_path], check=True,
                   stderr=subprocess.PIPE)
    subprocess.run(["bcftools", "index", "-f", split_vcf_path], check=True,
                   stderr=subprocess.PIPE)

    # Will protein level aggregation be performed on the variants?
    # This will only be done for simple missense variants but it is important
    # to annotate the vcf file before breaking down the haplotypes.
    if annotate:
        annotated_vcf_path = os.path.join(wdir, output_prefix + "split.ann."
                                          + vcf_file)
        res = annotate_vcf_file(settings, split_vcf_path, annotated_vcf_path)
        if res != 0:
            print("Annotating the vcf file failed.")
            return
    else:
        annotated_vcf_path = split_vcf_path
    if aggregate_aminoacids:
        if not (annotate or annotated_vcf):
            print("annotate option must be set to true or an annotadet vcf "
                  "file must be provided and annotated_vcf option must be "
                  "set to true for amino acid level aggregation. \n"
                  "Exiting!")
            return
        # check if a target annotation dict is provided.
        target_annotation_dict = {}
        if target_aa_annotation is not None:
            taa = pd.read_table(target_aa_annotation).set_index(
                ["gene_name", "aminoacid_change"]).to_dict(orient="index")
            for k in taa.keys():
                target_annotation_dict[k] = taa[k]["mutation_name"]
        # check if a gene id to gene name file is provided
        gene_ids = {}
        if geneid_to_genename is not None:
            gids = pd.read_table(geneid_to_genename).set_index("gene_id")
            gids = gids.to_dict(orient="index")
            for g in gids:
                gene_ids[g] = gids[g]["gene_name"]

        # load annotated vcf file
        variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                  alt_number=1,
                                  transformers=allel.ANNTransformer())
        # allel import provides a variants dictionary with keys such as
        # variants/AD, variants/POS for variant level information
        # the values are arrays with each element corresponding to one variant.
        # similarly, calldata/GT type keys hold the genotype level data.
        #############################################################
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        # find missense variant locations in the data. We are going to split
        # multi amino acid changes for missense variants only for target
        # annotation and count aggregation.
        missense = ["missense_variant" == variant for variant
                    in variants["variants/ANN_Annotation"]]
        # spcecify fields of interest from the INFO fields
        variant_fields = ["ANN_Gene_ID", "ANN_HGVS_p", "ANN_Annotation",
                          "QUAL"]
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO', 'calldata/DP',
                            'calldata/GT', 'calldata/GQ', 'calldata/QA',
                            'calldata/QR']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # so now we have a list of length equal to variant number.
        # each item is a tuple such as ('PF3D7_0104300', 'Gln107Leu') or
        # ('PF3D7_0104300', 'AspGluAsp144HisGlnTyr'). We'll split these
        # compound SNVs later.
        # get count data for missense variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # split the compound mutations
        split_variants = []
        split_calls = []
        for i in range(len(missense)):
            mv = variant_data[i][:3]
            # get the aa change such as AspGluAsp144HisGlnTyr
            aa_change = mv[1]
            # if no aa change, skip
            if aa_change == "":
                continue
            try:
                # if a mapping dict is present, add the gene name
                # this would get Pfubp1 from PF3D7_0104300, for example
                gene_name = gene_ids[mv[0]]
            except KeyError:
                gene_name = mv[0]

            # get site quality, remove those not satisfying min_site_qual
            # unless they are targeted mutations
            site_qual = float(variant_data[i][3])
            if missense[i]:
                # get the position of the change (144 above)
                aa_pos = int("".join([c for c in aa_change if c.isdigit()]))
                # split the aa change to reference aminoacid sequence and
                # alt amino acid sequence.
                aa_split = aa_change.split(str(aa_pos))
                reference = aa_split[0]
                alternate = aa_split[1]
                # aa changes are in 3 letter format. Loop through each aa and
                # split to single aa changes.
                for j in range(0, len(reference), 3):
                    new_pos = int(aa_pos + j/3)
                    # convert single amino acid names to 1 letter code.
                    new_reference = reference[j:j+3]
                    new_alternate = alternate[j:j+3]
                    new_change = new_reference + str(new_pos) + new_alternate
                    try:
                        # if this variant is in the targets, annotate it so.
                        mut_name = target_annotation_dict[
                            (gene_name, new_change)]
                        targeted_mutation = "Yes"
                        # reset alt observation counts to 0 if quality is low
                        if site_qual < min_target_site_qual:
                            call_data[i][0][:] = 0
                    except KeyError:
                        # remove low quality non-target alleles as well as
                        # synonymous changes
                        if ((site_qual < min_site_qual)
                                or (new_reference == new_alternate)):
                            continue
                        mut_name = gene_name + "-" + new_change
                        targeted_mutation = "No"
                    # add the split variant information split variants list
                    split_variants.append(mv + (new_change, gene_name,
                                                mut_name, targeted_mutation))
                    # add the individual level data to split calls list.
                    split_calls.append(call_data[i])
            else:
                try:
                    # if this variant is in the targets, annotate it as such.
                    mut_name = target_annotation_dict[
                        (gene_name, aa_change)]
                    targeted_mutation = "Yes"
                    if site_qual < min_target_site_qual:
                        call_data[i][0][:] = 0
                except KeyError:
                    # remove low qual or synonymous changes
                    if ((site_qual < min_site_qual)
                            or (mv[2] == "synonymous_variant")):
                        continue
                    mut_name = gene_name + "-" + aa_change
                    targeted_mutation = "No"
                # add compound variant data to split variant data
                split_variants.append(mv + (aa_change, gene_name,
                                            mut_name, targeted_mutation))
                # add the individual level data to split calls list.
                split_calls.append(call_data[i])

            # get individual level data
            genotype_quals = call_data[i][4]
            ao_count = call_data[i][0]
            alt_quals = call_data[i][5]
            average_alt_quals = alt_quals / ao_count
            ro_count = call_data[i][1]
            ref_quals = call_data[i][6]
            average_ref_quals = ref_quals / ro_count
            gq_mask = genotype_quals < min_genotype_qual
            qa_mask = alt_quals < min_alt_qual
            qr_mask = ref_quals < min_ref_qual
            av_qa_mask = average_alt_quals < min_mean_alt_qual
            av_qr_mask = average_ref_quals < min_mean_ref_qual
            # replace count data for individuals failing quality thresholds
            # alt allele count AO
            call_data[i][0][qa_mask] = 0
            call_data[i][0][av_qa_mask] = 0
            # ref allele count RO
            call_data[i][1][qr_mask] = 0
            call_data[i][1][av_qr_mask] = 0
            # reset coverage for gq failure
            call_data[i][2][gq_mask] = 0
            # reset genotypes for gq failure
            call_data[i][3][gq_mask] = -2

        # create a multiindex for the variant df that we'll create next
        index = pd.MultiIndex.from_tuples(
            split_variants, names=["Gene ID", "Compound Change", "ExonicFunc",
                                   "AA Change", "Gene", "Mutation Name",
                                   "Targeted"])
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # combine counts for same changes
        grouping_keys = ["Gene ID", "Gene", "Mutation Name", "ExonicFunc",
                         "AA Change", "Targeted"]
        # replace -1 (allel assigned NA values) values with 0
        # sum alt counts
        mutation_counts = variant_counts.groupby(grouping_keys).sum()
        # take the max of ref counts
        mutation_refs = reference_counts.groupby(grouping_keys).min()
        # take the max of coverage counts
        mutation_coverage = coverage.groupby(grouping_keys).max()
        # due to aggregating aa changes, ref counts can be overcounted even
        # if the minimum ref count is taken for the aggregate. The reason for
        # this is that each nucleotide variant's reference observation count
        # may include the alternate alleles for another nucleotide variant
        # that codes for the same aa change. So we'll set the ref counts
        # to coverage - alt count where ref count exceeds this value
        diff_count = mutation_coverage - mutation_counts
        ref_difference = (mutation_refs > diff_count).sum()
        # get the variant indices where ref count exceeds coverage - alt count
        exceed_index = ref_difference.loc[ref_difference > 0].index
        mutation_refs.loc[:, exceed_index] = diff_count.loc[:, exceed_index]
        # get genotypes as called by the variant caller
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index)

        def combine_gt(g):
            if 1 in g.values:
                return 1
            elif 0 in g.values:
                if 2 in g.values:
                    return 1
                else:
                    return 0
            elif 2 in g.values:
                return 2
            else:
                return -1

        gt_calls = gt_calls.groupby(grouping_keys).agg(combine_gt)

        # for one pf mutation alt count will be replaced with ref count
        # because reference allele is drug resistant
        dhps_key = ("PF3D7_0810800", "dhps", "dhps-Gly437Ala",
                    "missense_variant", "Gly437Ala", "Yes")
        dhps_new_key = ("PF3D7_0810800", "dhps", "dhps-Ala437Gly",
                        "missense_variant", "Ala437Gly", "Yes")
        try:
            mutation_counts.loc[dhps_new_key, :] = mutation_refs.loc[
                dhps_key, :]
            mutation_refs.loc[dhps_new_key, :] = mutation_counts.loc[
                dhps_key, :]
            mutation_coverage.loc[dhps_new_key, :] = mutation_coverage.loc[
                dhps_key, :]
            gt_calls.loc[dhps_new_key, :] = gt_calls.loc[
                dhps_key, :].replace({2: 0, 0: 2})
            gt_calls.drop(dhps_key, inplace=True)
            mutation_counts.drop(dhps_key, inplace=True)
            mutation_refs.drop(dhps_key, inplace=True)
            mutation_coverage.drop(dhps_key, inplace=True)
            mutation_counts = mutation_counts.sort_index()
            mutation_refs = mutation_refs.sort_index()
            mutation_coverage = mutation_coverage.sort_index()
            gt_calls = gt_calls.sort_index()
        except KeyError:
            pass

        # save count tables
        mutation_counts.T.to_csv(os.path.join(wdir, output_prefix
                                              + "alternate_AA_table.csv"))
        mutation_refs.T.to_csv(os.path.join(wdir, output_prefix
                                            + "reference_AA_table.csv"))
        mutation_coverage.T.to_csv(os.path.join(wdir, output_prefix
                                                + "coverage_AA_table.csv"))
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_AA_table.csv"))

    if aggregate_nucleotides:
        # aggregating counts of nucleotides requires decomposing block
        # substitutions, at a minimum. If desired, complex variants involving
        # indels can be decomposed as well.
        decomposed_vcf = os.path.join(wdir, output_prefix
                                      + "decomposed." + vcf_file)
        # prepare vt decompose command
        comm = ["vt", "decompose_blocksub"] + decompose_options
        comm.append(split_vcf_path)
        comm.extend(["-o", decomposed_vcf])
        # run decompose
        subprocess.run(comm, check=True)
        subprocess.run(["bcftools", "index", "-f", decomposed_vcf], check=True)
        # load decomposed vcf file
        variants = allel.read_vcf(decomposed_vcf, fields=["*"], alt_number=1)
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        # spcecify fields of interest from the INFO fields
        variant_fields = ["CHROM", "POS", "REF", "ALT", "QUAL"]
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO', 'calldata/DP',
                            'calldata/GT', 'calldata/GQ', 'calldata/QA',
                            'calldata/QR']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # get count data for the variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        # check if a target annotation dict is provided.
        target_annotation_dict = {}
        if target_nt_annotation is not None:
            taa = pd.read_table(target_nt_annotation).set_index(
                ["CHROM", "POS", "REF", "ALT"]).to_dict(orient="index")
            for k in taa.keys():
                target_annotation_dict[k] = taa[k]["mutation_name"]
        grouping_keys = ["CHROM", "POS", "REF", "ALT", "Mutation Name",
                         "Targeted"]
        split_variants = []
        split_calls = []
        for i in range(len(variant_data)):
            vd = variant_data[i][:4]
            site_qual = float(variant_data[i][4])
            try:
                t_anno = target_annotation_dict[vd]
                targeted_mutation = "Yes"
                if site_qual < min_target_site_qual:
                    call_data[i][0][:] = 0
            except KeyError:
                # remove low qual and nonvariant sites
                if ((site_qual < min_site_qual) or (vd[2] == vd[3])):
                    continue
                t_anno = ":".join(map(str, vd))
                targeted_mutation = "No"
            split_variants.append(vd + (t_anno, targeted_mutation))
            split_calls.append(call_data[i])

            # get individual level data
            genotype_quals = call_data[i][4]
            ao_count = call_data[i][0]
            alt_quals = call_data[i][5]
            average_alt_quals = alt_quals / ao_count
            ro_count = call_data[i][1]
            ref_quals = call_data[i][6]
            average_ref_quals = ref_quals / ro_count
            gq_mask = genotype_quals < min_genotype_qual
            qa_mask = alt_quals < min_alt_qual
            qr_mask = ref_quals < min_ref_qual
            av_qa_mask = average_alt_quals < min_mean_alt_qual
            av_qr_mask = average_ref_quals < min_mean_ref_qual
            # replace count data for individuals failing quality thresholds
            # alt allele count AO
            call_data[i][0][qa_mask] = 0
            call_data[i][0][av_qa_mask] = 0
            # ref allele count RO
            call_data[i][1][qr_mask] = 0
            call_data[i][1][av_qr_mask] = 0
            # reset coverage for gq failure
            call_data[i][2][gq_mask] = 0
            # reset genotypes for gq failure
            call_data[i][3][gq_mask] = -2

        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # create a multiindex for the variant df that we'll create next
        index = pd.MultiIndex.from_tuples(
            split_variants, names=grouping_keys)
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # combine counts for same changes
        # sum alt counts
        mutation_counts = variant_counts.groupby(grouping_keys).sum()
        # take the max of ref counts
        mutation_refs = reference_counts.groupby(grouping_keys).min()
        # take the max of coverage counts
        mutation_coverage = coverage.groupby(grouping_keys).max()
        # save count tables
        mutation_counts.T.to_csv(os.path.join(wdir, output_prefix
                                              + "alternate_AN_table.csv"))
        mutation_refs.T.to_csv(os.path.join(wdir, output_prefix
                                            + "reference_AN_table.csv"))
        mutation_coverage.T.to_csv(os.path.join(wdir, output_prefix
                                                + "coverage_AN_table.csv"))
        # get genotypes
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index)

        def combine_gt(g):
            if 1 in g.values:
                return 1
            elif 0 in g.values:
                if 2 in g.values:
                    return 1
                else:
                    return 0
            elif 2 in g.values:
                return 2
            else:
                return -1

        gt_calls = gt_calls.groupby(grouping_keys).agg(combine_gt)
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_AN_table.csv"))

    if aggregate_none:
        # if no aggregation will be done, load the vcf file
        if annotate or annotated_vcf:
            # if annotation was requested use the annotated vcf path
            variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                      alt_number=1,
                                      transformers=allel.ANNTransformer())
        else:
            # if the file is not annotated, don't try to parse ANN field.
            variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                      alt_number=1)
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        variant_fields = ["CHROM", "POS", "REF", "ALT", "QUAL"]
        if annotate or annotated_vcf:
            variant_fields.extend(["ANN_Gene_ID", "ANN_HGVS_p"])
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO', 'calldata/DP',
                            'calldata/GT', 'calldata/GQ', 'calldata/QA',
                            'calldata/QR']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # get count data for the variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        split_variants = []
        split_calls = []
        for i in range(len(variant_data)):
            vd = variant_data[i][:4]
            site_qual = float(variant_data[i][4])
            if site_qual < min_site_qual:
                continue
            if annotate or annotated_vcf:
                g_ann = variant_data[i][5]
                p_ann = variant_data[i][6]
                if p_ann == "":
                    p_ann = "."
                if g_ann == "":
                    g_ann = "."
            else:
                p_ann = "."
                g_ann = "."
            vd = vd + (g_ann, p_ann)
            split_variants.append(vd)
            split_calls.append(call_data[i])

            # get individual level data
            genotype_quals = call_data[i][4]
            ao_count = call_data[i][0]
            alt_quals = call_data[i][5]
            average_alt_quals = alt_quals / ao_count
            ro_count = call_data[i][1]
            ref_quals = call_data[i][6]
            average_ref_quals = ref_quals / ro_count
            gq_mask = genotype_quals < min_genotype_qual
            qa_mask = alt_quals < min_alt_qual
            qr_mask = ref_quals < min_ref_qual
            av_qa_mask = average_alt_quals < min_mean_alt_qual
            av_qr_mask = average_ref_quals < min_mean_ref_qual
            # replace count data for individuals failing quality thresholds
            # alt allele count AO
            call_data[i][0][qa_mask] = 0
            call_data[i][0][av_qa_mask] = 0
            # ref allele count RO
            call_data[i][1][qr_mask] = 0
            call_data[i][1][av_qr_mask] = 0
            # reset coverage for gq failure
            call_data[i][2][gq_mask] = 0
            # reset genotypes for gq failure
            call_data[i][3][gq_mask] = -2
        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # create a multiindex for the variant df that we'll create next
        variant_fields = variant_fields[:4] + [
            "variants/Gene ID", "variants/AA Change"]
        index = pd.MultiIndex.from_tuples(split_variants,
                                          names=[v.split("variants/")[1]
                                                 for v in variant_fields])
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # save count tables
        variant_counts.T.to_csv(os.path.join(wdir, output_prefix
                                             + "alternate_table.csv"))
        reference_counts.T.to_csv(os.path.join(wdir, output_prefix
                                               + "reference_table.csv"))
        coverage.T.to_csv(os.path.join(wdir, output_prefix
                                       + "coverage_table.csv"))

        # get genotypes
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index).replace(-2, -1)
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_table.csv"))


def vcf_to_tables(vcf_file, settings=None, settings_file=None, annotate=True,
                  geneid_to_genename=None, target_aa_annotation=None,
                  aggregate_aminoacids=False, target_nt_annotation=None,
                  aggregate_nucleotides=False, decompose_options=[],
                  annotated_vcf=False, aggregate_none=False, min_site_qual=-1,
                  min_target_site_qual=-1, min_genotype_qual=None,
                  output_prefix=""):
    """Create various tables from a vcf file.

    Create various tables from a vcf file generated by the freebayes
    program. There are 3 different types of count output for each variant:
    variant count, reference count and coverage. The vcf file will be split
    into biallelic variants. Table versions of the input vcf will be created
    but the info fields will be limited to the mandatory vcf fields and some
    annotation data if avaliable.

    In addition to the original vcf table, aa change tables can be generated.
    These will be generated by filtering the vcf to missense variants only,
    decomposing block substitutions (haplotypes) and combining the counts for
    the same aminoacid changes. This operation is specifically intended for
    generating data for targeted missense mutations and only reports that. All
    other variants, even those complex variants including targeted variants
    will not be reported. Finally, one specific mutation (dhps-437) will have
    reference counts instead of variant counts if present. This is because this
    drug resistance variant is encoded by the 3d7 reference sequence.

    Parameters
    ----------
    settings: dict, None
        Analysis settings dictionary. Either this or settings_file must
        be provided.
    settings_file: str/path, None
        Path to the analysis settings file. Either this or the settings dict
        must be provided.
    annotate: bool, True
        Annotate variant file. This is required for protein level analysis.
    vcf_file: str/path
        Starting vcf file.
    geneid2genename: str/path, None.
        Path to a tab separated tex file that maps gene ids to gene names.
        Column names must be gene_id and gene_name. Gene IDs
        will populate the Gene field if this file is not provided.
    target_aa_annotation: str/path, None.
        Path to a tab separated text file with targeted variant information to
        annotate and label targeted amino acid changes.
        It must have gene_name, aminoacid_change, and mutation_name columns.
        Amino acid changes should be represented as refAAPosAltAA. refAA and
        AltAA must be three letter amino acid codes.
        This file is required for targeted protein variant labeling.
    target_nt_annotation: str/path, None.
        Path to a tab separated text file with targeted variant information to
        annotate and label targeted nucleotide changes.
        It must have CHROM, POS, REF, ALT, NAME columns.
        This file is required for targeted nucleotide variant labeling.
    aggregate_aminoacids: bool, False
        whether counts for same amino acids should be aggregated. This involves
        decomposing multi amino acid changes for missense variants. If amino
        acid based targets will be annotated, based on a provided annotation
        dictionary, aggregation step must be completed. Targeted mutations
        that are part of complex events (indels, stop loss/gain etc.) will not
        be labeled as targeted.
    aggregate_nucleotides: bool, False
        whether the counts for nucleotide changes should be aggregated. This
        involves decomposing all variants to the smallest units possible,
        breaking all haplotype data. The level of decomposition should be
        specified with the decompose_options parameter.
    aggregate_none: bool, False.
        Do no aggregation on counts, save the original (annotated if requested)
        vcf file as 3  count tables. Three aggregation options are compatible
        with each other and can be used all at once.
    decompose_options: list, []
        if aggregate nucleotides option is selected, these options will be
        passed to vt program. "-a" for decomposing variants containing indels,
        for example. "-p" for keeping phase information. Any option to vt
        decompose_blocksub would be valid. By default indels will not be
        decomposed.
    annotated_vcf: bool, False
        is the provided vcf file annotated using snpEff. These annotations
        will be used if no count aggregation is to be done and annotate option
        is False.
    min_site_qual: float, -1
        Filter variants with QUAL values less than this value if the site is
        not a targeted site. If targeted, the site will be kept regardless of
        the qual value for the site. freebayes manual indicates that
        simulations showed a value between 1-30 would be good. So a minimum
        value of 1 here would clean up most junk sites.
    min_target_site_qual: float, -1
        If a variant site is targeted but the site qual is lower than this,
        reset the alternate observation counts to 0. It may be best to leave
        this at the default value since there is usually additional evidence
        that a targeted variant exists in a samples compared to a de novo
        variant.
    """
    # get the analysis settings
    # check if both settings and the settings file are None:
    if (settings is None) and (settings_file is None):
        print("settings or settings file must be provided for freebayes_call.")
        return
    else:
        if settings is None:
            settings = get_analysis_settings(settings_file)
        else:
            settings = copy.deepcopy(settings)
    # get the working directory from settings
    wdir = settings["workingDir"]
    # All postprocessing steps require biallelic variant representation.
    # so we'll use bcftools to split multiallelics to their own lines.
    genome_fasta = get_file_locations()[settings["species"]]["fasta_genome"]
    vcf_path = os.path.join(wdir, vcf_file)
    # filter genotype for quality if specified
    if min_genotype_qual is not None:
        if vcf_file.endswith(".gz"):
            vtype = "--gzvcf"
        else:
            vtype = "--vcf"
        filt_res = subprocess.Popen(["vcftools", vtype, vcf_path,
                                     "--minGQ", str(min_genotype_qual),
                                     "--recode", "--recode-INFO-all",
                                     "--stdout"],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        filt_vcf_path = os.path.join(
            wdir, output_prefix + "variants.GQ."
            + str(min_genotype_qual) + ".vcf.gz")
        with open(filt_vcf_path, "wb") as outfile:
            zip_res = subprocess.run(["bgzip", "-f"], stdin=filt_res.stdout,
                                     stdout=outfile,
                                     stderr=subprocess.PIPE)
        index_res = subprocess.run(
            ["bcftools", "index", "-f", filt_vcf_path],
            stderr=subprocess.PIPE)
        if zip_res.returncode != 0:
            print(("Compression of GQ filtered vcf failed due to "
                   "error: {}. \n Genotypes will not be "
                   "filtered.").format(zip_res.stderr))
        elif index_res.returncode != 0:
            print(("Indexing GQ filtered vcf file failed "
                   "due to error: {}. \n Genotypes will not "
                   "be filtered.").format(index_res.stderr))
        else:
            vcf_path = filt_vcf_path

    split_vcf_path = os.path.join(wdir, output_prefix + "split." + vcf_file)
    subprocess.run(["bcftools", "norm", "-f", genome_fasta, "-m-both",
                    vcf_path, "-Oz", "-o", split_vcf_path], check=True,
                   stderr=subprocess.PIPE)
    subprocess.run(["bcftools", "index", "-f", split_vcf_path], check=True,
                   stderr=subprocess.PIPE)

    # Will protein level aggregation be performed on the variants?
    # This will only be done for simple missense variants but it is important
    # to annotate the vcf file before breaking down the haplotypes.
    if annotate:
        annotated_vcf_path = os.path.join(wdir, output_prefix + "split.ann."
                                          + vcf_file)
        res = annotate_vcf_file(settings, split_vcf_path, annotated_vcf_path)
        if res != 0:
            print("Annotating the vcf file failed.")
            return
    else:
        annotated_vcf_path = split_vcf_path
    if aggregate_aminoacids:
        if not (annotate or annotated_vcf):
            print("annotate option must be set to true or an annotadet vcf "
                  "file must be provided and annotated_vcf option must be "
                  "set to true for amino acid level aggregation. \n"
                  "Exiting!")
            return
        # check if a target annotation dict is provided.
        target_annotation_dict = {}
        if target_aa_annotation is not None:
            taa = pd.read_table(target_aa_annotation).set_index(
                ["gene_name", "aminoacid_change"]).to_dict(orient="index")
            for k in taa.keys():
                target_annotation_dict[k] = taa[k]["mutation_name"]
        # check if a gene id to gene name file is provided
        gene_ids = {}
        if geneid_to_genename is not None:
            gids = pd.read_table(geneid_to_genename).set_index("gene_id")
            gids = gids.to_dict(orient="index")
            for g in gids:
                gene_ids[g] = gids[g]["gene_name"]

        # load annotated vcf file
        variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                  alt_number=1,
                                  transformers=allel.ANNTransformer())
        # allel import provides a variants dictionary with keys such as
        # variants/AD, variants/POS for variant level information
        # the values are arrays with each element corresponding to one variant.
        # similarly, calldata/GT type keys hold the genotype level data.
        #############################################################
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        # find missense variant locations in the data. We are going to split
        # multi amino acid changes for missense variants only for target
        # annotation and count aggregation.
        missense = ["missense_variant" == variant for variant
                    in variants["variants/ANN_Annotation"]]
        # spcecify fields of interest from the INFO fields
        variant_fields = ["ANN_Gene_ID", "ANN_HGVS_p", "ANN_Annotation",
                          "QUAL"]
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO',
                            'calldata/DP', 'calldata/GT']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # so now we have a list of length equal to variant number.
        # each item is a tuple such as ('PF3D7_0104300', 'Gln107Leu') or
        # ('PF3D7_0104300', 'AspGluAsp144HisGlnTyr'). We'll split these
        # compound SNVs later.
        # get count data for missense variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # split the compound mutations
        split_variants = []
        split_calls = []
        for i in range(len(missense)):
            mv = variant_data[i][:3]
            # get the aa change such as AspGluAsp144HisGlnTyr
            aa_change = mv[1]
            # if no aa change, skip
            if aa_change == "":
                continue
            try:
                # if a mapping dict is present, add the gene name
                # this would get Pfubp1 from PF3D7_0104300, for example
                gene_name = gene_ids[mv[0]]
            except KeyError:
                gene_name = mv[0]

            # get site quality, remove those not satisfying min_site_qual
            # unless they are targeted mutations
            site_qual = float(variant_data[i][3])
            if missense[i]:
                # get the position of the change (144 above)
                aa_pos = int("".join([c for c in aa_change if c.isdigit()]))
                # split the aa change to reference aminoacid sequence and
                # alt amino acid sequence.
                aa_split = aa_change.split(str(aa_pos))
                reference = aa_split[0]
                alternate = aa_split[1]
                # aa changes are in 3 letter format. Loop through each aa and
                # split to single aa changes.
                for j in range(0, len(reference), 3):
                    new_pos = int(aa_pos + j/3)
                    # convert single amino acid names to 1 letter code.
                    new_reference = reference[j:j+3]
                    new_alternate = alternate[j:j+3]
                    new_change = new_reference + str(new_pos) + new_alternate
                    try:
                        # if this variant is in the targets, annotate it so.
                        mut_name = target_annotation_dict[
                            (gene_name, new_change)]
                        targeted_mutation = "Yes"
                        # reset alt observation counts to 0 if quality is low
                        if site_qual < min_target_site_qual:
                            call_data[i][0][:] = 0
                    except KeyError:
                        # remove low quality non-target alleles as well as
                        # synonymous changes
                        if ((site_qual < min_site_qual)
                                or (new_reference == new_alternate)):
                            continue
                        mut_name = gene_name + "-" + new_change
                        targeted_mutation = "No"
                    # add the split variant information split variants list
                    split_variants.append(mv + (new_change, gene_name,
                                                mut_name, targeted_mutation))
                    # add the individual level data to split calls list.
                    split_calls.append(call_data[i])
            else:
                try:
                    # if this variant is in the targets, annotate it as such.
                    mut_name = target_annotation_dict[
                        (gene_name, aa_change)]
                    targeted_mutation = "Yes"
                    if site_qual < min_target_site_qual:
                        call_data[i][0][:] = 0
                except KeyError:
                    # remove low qual or synonymous changes
                    if ((site_qual < min_site_qual)
                            or (mv[2] == "synonymous_variant")):
                        continue
                    mut_name = gene_name + "-" + aa_change
                    targeted_mutation = "No"
                # add compound variant data to split variant data
                split_variants.append(mv + (aa_change, gene_name,
                                            mut_name, targeted_mutation))
                # add the individual level data to split calls list.
                split_calls.append(call_data[i])
        # create a multiindex for the variant df that we'll create next
        index = pd.MultiIndex.from_tuples(
            split_variants, names=["Gene ID", "Compound Change", "ExonicFunc",
                                   "AA Change", "Gene", "Mutation Name",
                                   "Targeted"])
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # combine counts for same changes
        grouping_keys = ["Gene ID", "Gene", "Mutation Name", "ExonicFunc",
                         "AA Change", "Targeted"]
        # replace -1 (allel assigned NA values) values with 0
        # sum alt counts
        mutation_counts = variant_counts.groupby(grouping_keys).sum()
        # take the max of ref counts
        mutation_refs = reference_counts.groupby(grouping_keys).min()
        # take the max of coverage counts
        mutation_coverage = coverage.groupby(grouping_keys).max()
        # due to aggregating aa changes, ref counts can be overcounted even
        # if the minimum ref count is taken for the aggregate. The reason for
        # this is that each nucleotide variant's reference observation count
        # may include the alternate alleles for another nucleotide variant
        # that codes for the same aa change. So we'll set the ref counts
        # to coverage - alt count where ref count exceeds this value
        diff_count = mutation_coverage - mutation_counts
        ref_difference = (mutation_refs > diff_count).sum()
        # get the variant indices where ref count exceeds coverage - alt count
        exceed_index = ref_difference.loc[ref_difference > 0].index
        mutation_refs.loc[:, exceed_index] = diff_count.loc[:, exceed_index]
        # get genotypes as called by the variant caller
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index)

        def combine_gt(g):
            if 1 in g.values:
                return 1
            elif 0 in g.values:
                if 2 in g.values:
                    return 1
                else:
                    return 0
            elif 2 in g.values:
                return 2
            else:
                return -1

        gt_calls = gt_calls.groupby(grouping_keys).agg(combine_gt)

        # for one pf mutation alt count will be replaced with ref count
        # because reference allele is drug resistant
        dhps_key = ("PF3D7_0810800", "dhps", "dhps-Gly437Ala",
                    "missense_variant", "Gly437Ala", "Yes")
        dhps_new_key = ("PF3D7_0810800", "dhps", "dhps-Ala437Gly",
                        "missense_variant", "Ala437Gly", "Yes")
        try:
            mutation_counts.loc[dhps_new_key, :] = mutation_refs.loc[
                dhps_key, :]
            mutation_refs.loc[dhps_new_key, :] = mutation_counts.loc[
                dhps_key, :]
            mutation_coverage.loc[dhps_new_key, :] = mutation_coverage.loc[
                dhps_key, :]
            gt_calls.loc[dhps_new_key, :] = gt_calls.loc[
                dhps_key, :].replace({2: 0, 0: 2})
            gt_calls.drop(dhps_key, inplace=True)
            mutation_counts.drop(dhps_key, inplace=True)
            mutation_refs.drop(dhps_key, inplace=True)
            mutation_coverage.drop(dhps_key, inplace=True)
            mutation_counts = mutation_counts.sort_index()
            mutation_refs = mutation_refs.sort_index()
            mutation_coverage = mutation_coverage.sort_index()
            gt_calls = gt_calls.sort_index()
        except KeyError:
            pass

        # save count tables
        mutation_counts.T.to_csv(os.path.join(wdir, output_prefix
                                              + "alternate_AA_table.csv"))
        mutation_refs.T.to_csv(os.path.join(wdir, output_prefix
                                            + "reference_AA_table.csv"))
        mutation_coverage.T.to_csv(os.path.join(wdir, output_prefix
                                                + "coverage_AA_table.csv"))
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_AA_table.csv"))

    if aggregate_nucleotides:
        # aggregating counts of nucleotides requires decomposing block
        # substitutions, at a minimum. If desired, complex variants involving
        # indels can be decomposed as well.
        decomposed_vcf = os.path.join(wdir, output_prefix
                                      + "decomposed." + vcf_file)
        # prepare vt decompose command
        comm = ["vt", "decompose_blocksub"] + decompose_options
        comm.append(split_vcf_path)
        comm.extend(["-o", decomposed_vcf])
        # run decompose
        subprocess.run(comm, check=True)
        subprocess.run(["bcftools", "index", "-f", decomposed_vcf], check=True)
        # load decomposed vcf file
        variants = allel.read_vcf(decomposed_vcf, fields=["*"], alt_number=1)
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        # spcecify fields of interest from the INFO fields
        variant_fields = ["CHROM", "POS", "REF", "ALT", "QUAL"]
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO',
                            'calldata/DP', 'calldata/GT']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # get count data for the variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        # check if a target annotation dict is provided.
        target_annotation_dict = {}
        if target_nt_annotation is not None:
            taa = pd.read_table(target_nt_annotation).set_index(
                ["CHROM", "POS", "REF", "ALT"]).to_dict(orient="index")
            for k in taa.keys():
                target_annotation_dict[k] = taa[k]["mutation_name"]
        grouping_keys = ["CHROM", "POS", "REF", "ALT", "Mutation Name",
                         "Targeted"]
        split_variants = []
        split_calls = []
        for i in range(len(variant_data)):
            vd = variant_data[i][:4]
            site_qual = float(variant_data[i][4])
            try:
                t_anno = target_annotation_dict[vd]
                targeted_mutation = "Yes"
                if site_qual < min_target_site_qual:
                    call_data[i][0][:] = 0
            except KeyError:
                # remove low qual and nonvariant sites
                if ((site_qual < min_site_qual) or (vd[2] == vd[3])):
                    continue
                t_anno = ":".join(map(str, vd))
                targeted_mutation = "No"
            split_variants.append(vd + (t_anno, targeted_mutation))
            split_calls.append(call_data[i])
        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # create a multiindex for the variant df that we'll create next
        index = pd.MultiIndex.from_tuples(
            split_variants, names=grouping_keys)
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # combine counts for same changes
        # sum alt counts
        mutation_counts = variant_counts.groupby(grouping_keys).sum()
        # take the max of ref counts
        mutation_refs = reference_counts.groupby(grouping_keys).min()
        # take the max of coverage counts
        mutation_coverage = coverage.groupby(grouping_keys).max()
        # save count tables
        mutation_counts.T.to_csv(os.path.join(wdir, output_prefix
                                              + "alternate_AN_table.csv"))
        mutation_refs.T.to_csv(os.path.join(wdir, output_prefix
                                            + "reference_AN_table.csv"))
        mutation_coverage.T.to_csv(os.path.join(wdir, output_prefix
                                                + "coverage_AN_table.csv"))
        # get genotypes
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index)

        def combine_gt(g):
            if 1 in g.values:
                return 1
            elif 0 in g.values:
                if 2 in g.values:
                    return 1
                else:
                    return 0
            elif 2 in g.values:
                return 2
            else:
                return -1

        gt_calls = gt_calls.groupby(grouping_keys).agg(combine_gt)
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_AN_table.csv"))

    if aggregate_none:
        # if no aggregation will be done, load the vcf file
        if annotate or annotated_vcf:
            # if annotation was requested use the annotated vcf path
            variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                      alt_number=1,
                                      transformers=allel.ANNTransformer())
        else:
            # if the file is not annotated, don't try to parse ANN field.
            variants = allel.read_vcf(annotated_vcf_path, fields=["*"],
                                      alt_number=1)
        # Freebayes vcfs have AO and RO counts for alt and ref allele depths
        # but GATK has a combined AD depth. Create AO and RO from AD if
        # needed
        try:
            variants["calldata/AO"]
        except KeyError:
            variants["calldata/RO"] = variants["calldata/AD"][:, :, 0]
            variants["calldata/AO"] = variants["calldata/AD"][:, :, 1]
        variant_fields = ["CHROM", "POS", "REF", "ALT", "QUAL"]
        if annotate or annotated_vcf:
            variant_fields.extend(["ANN_Gene_ID", "ANN_HGVS_p"])
        variant_fields = ["variants/" + v for v in variant_fields]
        # specify fields of interest from individual level data
        # that is basically the count data for tables. AO: alt allele count,
        # RO ref count, DP: coverage.
        call_data_fields = ['calldata/AO', 'calldata/RO',
                            'calldata/DP', 'calldata/GT']
        variants["calldata/GT"] = variants["calldata/GT"].sum(axis=2)
        # zip variant level  information together, so we have a single value
        # for each variant
        variant_data = list(zip(*[variants[v] for v in variant_fields]))
        # get count data for the variants
        call_data = list(zip(*[variants[c] for c in call_data_fields]))
        split_variants = []
        split_calls = []
        for i in range(len(variant_data)):
            vd = variant_data[i][:4]
            site_qual = float(variant_data[i][4])
            if site_qual < min_site_qual:
                continue
            if annotate or annotated_vcf:
                g_ann = variant_data[i][5]
                p_ann = variant_data[i][6]
                if p_ann == "":
                    p_ann = "."
                if g_ann == "":
                    g_ann = "."
            else:
                p_ann = "."
                g_ann = "."
            vd = vd + (g_ann, p_ann)
            split_variants.append(vd)
            split_calls.append(call_data[i])
        # first item of the above list is alt counts, then ref counts and
        # coverage.
        #############################
        # create a multiindex for the variant df that we'll create next
        variant_fields = variant_fields[:4] + [
            "variants/Gene ID", "variants/AA Change"]
        index = pd.MultiIndex.from_tuples(split_variants,
                                          names=[v.split("variants/")[1]
                                                 for v in variant_fields])
        # get alt counts
        variant_counts = pd.DataFrame(np.array(split_calls)[:, 0],
                                      columns=variants["samples"],
                                      index=index).replace(-1, 0)
        # get reference counts
        reference_counts = pd.DataFrame(np.array(split_calls)[:, 1],
                                        columns=variants["samples"],
                                        index=index).replace(-1, 0)
        # get coverage depth
        coverage = pd.DataFrame(np.array(split_calls)[:, 2],
                                columns=variants["samples"],
                                index=index).replace(-1, 0)
        # save count tables
        variant_counts.T.to_csv(os.path.join(wdir, output_prefix
                                             + "alternate_table.csv"))
        reference_counts.T.to_csv(os.path.join(wdir, output_prefix
                                               + "reference_table.csv"))
        coverage.T.to_csv(os.path.join(wdir, output_prefix
                                       + "coverage_table.csv"))

        # get genotypes
        gt_calls = pd.DataFrame((np.array(split_calls)[:, 3]),
                                columns=variants["samples"],
                                index=index).replace(-2, -1)
        gt_calls.T.to_csv(os.path.join(wdir, output_prefix
                                       + "genotypes_table.csv"))



def get_mutation_position(change):
    digits = []
    found = False
    for i in change:
        if i.isdigit():
            digits.append(i)
            found = True
        elif found:
            break
    return int("".join(digits))


def merge_contigs(settings, contig_info_dict, results):
    # merge contig vcfs for each chromosome
    wdir = settings["workingDir"]
    species = settings["species"]
    genome_fasta = get_file_locations()[species]["fasta_genome"]
    vcfdir = os.path.join(wdir, "msa_vcfs")
    vcf_file_list = []
    if not os.path.exists(vcfdir):
        os.makedirs(vcfdir)
    for chrom in contig_info_dict:
        chrom_vcf_list = os.path.join(wdir, chrom + "_vcf_files.txt")
        chrom_vcf_file = os.path.join(wdir, chrom + ".vcf.gz")
        with open(chrom_vcf_list, "w") as outf:
            for contig in contig_info_dict[chrom]:
                contig_name = contig_info_dict[chrom][contig]["contig_name"]
                if contig_name in results:
                    contigs_dir = contig_info_dict[chrom][contig][
                        "contigs_dir"]
                    contig_vcf_file = os.path.join(contigs_dir,
                                                   contig_name + ".vcf")
                    subprocess.call(["bgzip", "-f", contig_vcf_file],
                                    cwd=contigs_dir)
                    subprocess.call(["bcftools", "index", "-f",
                                     contig_vcf_file + ".gz"],
                                    cwd=contigs_dir)
                    outf.write(contig_vcf_file + ".gz" + "\n")
        subprocess.call(["bcftools", "concat", "-f", chrom_vcf_list, "-Oz",
                         "-o", chrom_vcf_file])
        vcf_file_list.append(chrom_vcf_file)

        split_vcf_file = os.path.join(wdir, chrom + ".split.vcf.gz")
        subprocess.call(["bcftools", "norm", "-m-both", "-N", "-Oz",
                         chrom_vcf_file, "-o", split_vcf_file])
        vcf_file_list.append(split_vcf_file)

        filt_vcf_file = os.path.join(wdir, chrom + ".split.filt.vcf.gz")

        minVariantBarcodes = settings["minVariantBarcodes"]
        minVariantSamples = settings["minVariantSamples"]
        minVariantSampleFraction = settings["minVariantSampleFraction"]
        minVariantSampleTotal = settings["minVariantSampleTotal"]
        minVariantMeanQuality = settings["minVariantMeanQuality"]
        minVariantMeanWsaf = settings["minVariantMeanWsaf"]
        minMipCountFraction = settings["minMipCountFraction"]

        filter_expressions = [
            "((INFO/AD[1] >= " + minVariantBarcodes + ")",
            "(INFO/SC[1] >= " + minVariantSamples + ")",
            "(INFO/SF[1] >= " + minVariantSampleFraction + ")",
            "(INFO/NS >= " + minVariantSampleTotal + ")",
            "(INFO/QS[1] >= " + minVariantMeanQuality + ")",
            "(INFO/WSAF[1] >= " + minVariantMeanWsaf + ")",
            "(INFO/MCF[1] >= " + minMipCountFraction + ")"]

        filter_expressions = " & ".join(filter_expressions)
        filter_expressions = filter_expressions + ') | (OT !=".")'

        subprocess.call(["bcftools", "view", "-i", filter_expressions, "-Oz",
                         split_vcf_file, "-o", filt_vcf_file])
        vcf_file_list.append(filt_vcf_file)

        merged_vcf_file = os.path.join(wdir, chrom + ".merged.filt.vcf.gz")
        subprocess.call(["bcftools", "norm", "-m+any", "-N", "-Oz",
                         filt_vcf_file, "-o", merged_vcf_file])
        vcf_file_list.append(merged_vcf_file)

        norm_vcf_file = os.path.join(wdir, chrom + ".norm.vcf.gz")
        subprocess.call(["bcftools", "norm", "-m-both", "-f", genome_fasta,
                         "-cs", "-Oz", merged_vcf_file, "-o", norm_vcf_file])
        vcf_file_list.append(norm_vcf_file)

        # annotate with snpEff
        try:
            ann_db_dir = get_file_locations()[species]["snpeff_dir"]
            ann_db = get_file_locations()[species]["snpeff_db"]
            annotate = True
        except KeyError:
            annotate = False
        if annotate:
            ann = subprocess.Popen(["java", "-Xmx10g", "-jar",
                                    os.path.join(ann_db_dir, "snpEff.jar"),
                                    ann_db, norm_vcf_file],
                                   stdout=subprocess.PIPE)
            annotated_vcf_file = os.path.join(wdir, chrom + ".norm.ann.vcf.gz")
            with open(annotated_vcf_file, "wb") as avf:
                subprocess.call(["bgzip"], stdin=ann.stdout,
                                stdout=avf)

            vcf_file_list.append(annotated_vcf_file)
        subprocess.call(["mv"] + vcf_file_list + [vcfdir])


def annotate_vcf_file(settings, vcf_file, annotated_vcf_file, options=[]):
    """Annotate a vcf file using snpEff, bgzip and index the output file."""
    # get the species information from settings
    species = settings["species"]
    try:
        # find where snpEff files are located and which database should be used
        ann_db_dir = get_file_locations()[species]["snpeff_dir"]
        ann_db = get_file_locations()[species]["snpeff_db"]
    except KeyError:
        print("snpeff_dir and snpeff_db must be specified in the settings "
              "to carry out snpeff annotations.")
        return
    # run snpeff program on the vcf file. Snpeff outputs to stdout so we'll
    # redirect it to the annotated vcf file. If output file name provided
    # ends with .gz, we will remove it here because bgzip will add that in
    # the next step
    if annotated_vcf_file.endswith(".gz"):
        annotated_vcf_file = annotated_vcf_file[:-3]
    with open(annotated_vcf_file, "wb") as avf:
        comm = ["java", "-Xmx10g", "-jar",
                os.path.join(ann_db_dir, "snpEff.jar"), ann_db, vcf_file]
        comm.extend(options)
        res = subprocess.run(comm, stdout=avf, stderr=subprocess.PIPE)
        if res.returncode != 0:
            print("Error  in snpEff call ", res.stderr)
            return res.returncode
    # most vcf operations require a bgzipped indexed file, so do those
    res = subprocess.run(["bgzip", "-f", annotated_vcf_file],
                         stderr=subprocess.PIPE)
    if res.returncode != 0:
        print("Error in compressing the annotated vcf file, ", res.stderr)
        return res.returncode
    res = subprocess.run(["bcftools", "index", "-f",
                          annotated_vcf_file + ".gz"], stderr=subprocess.PIPE)
    if res.returncode != 0:
        print("Error in indexing the annotated vcf file, ", res.stderr)
        return res.returncode
    return 0


def process_contig(contig_dict):
    try:
        chrom = contig_dict["chrom"]
        contig_start = contig_dict["contig_start"]
        contig_end = contig_dict["contig_end"]
        species = contig_dict["species"]
        contig_ref_seq = get_sequence(create_region(
            chrom, contig_start, contig_end), species)
        contig_haplotypes_file = contig_dict["contig_haplotypes_file"]
        contig_haps = pd.read_csv(contig_haplotypes_file)
        nastring = ".:.:.:.:.:.:."
        # Create a contig sequence for each haplotype.
        # This will be done by gettig the forward strand sequence for each
        # haplotype and padding it on both flanks with the reference sequence
        # up to the contig start/end.
        #
        # get forward strand sequence for all haplotypes
        contig_haps["forward_sequence"] = contig_haps["haplotype_sequence"]
        reverse_index = contig_haps["orientation"] == "reverse"
        contig_haps.loc[reverse_index, "forward_sequence"] = (
            contig_haps.loc[reverse_index, "forward_sequence"].apply(
                reverse_complement))

        def get_padded_sequence(row):
            chrom = row["Chrom"]
            contig_start = int(row["contig_start"])
            contig_end = int(row["contig_end"])
            capture_start = int(row["capture_start"])
            capture_end = int(row["capture_end"])
            left_key = create_region(chrom, contig_start, capture_start - 1)
            right_key = create_region(chrom, capture_end + 1, contig_end)
            left_pad = get_sequence(left_key, species)
            right_pad = get_sequence(right_key, species)
            return left_pad + str(row["forward_sequence"]) + right_pad

        contig_haps["padded_sequence"] = contig_haps.apply(
            get_padded_sequence, axis=1)
        g_dict = contig_haps.set_index(
            ["MIP", "Copy", "haplotype_ID"]).to_dict(orient="index")
        sequences = {"ref": contig_ref_seq}
        contig_targets = contig_dict["contig_targets"]
        if contig_targets is not None:
            contig_targets["padded_sequence"] = contig_targets.apply(
                get_padded_sequence, axis=1)
            target_pos = contig_targets[
                ["Pos", "End", "Mutation Name"]].to_dict(orient="records")
            targets_dict = contig_targets.to_dict(orient="index")
            for t in targets_dict:
                sequences[t] = targets_dict[t]["padded_sequence"]
        else:
            targets_dict = {}
            target_pos = []
        for k in g_dict.keys():
            sequences[":".join(k)] = g_dict[k]["padded_sequence"]
        wdir = contig_dict["contigs_dir"]
        contig_name = contig_dict["contig_name"]
        fasta_file = os.path.join(wdir, contig_name + ".fa")
        alignment_file = os.path.join(wdir, contig_name + ".aln")
        save_fasta_dict(sequences, fasta_file)
        if contig_dict["aligner"] == "muscle":
            mh = contig_dict["max_hours"]
            subprocess.call(["muscle", "-in", fasta_file, "-out",
                             alignment_file, "-maxhours", mh])
        elif contig_dict["aligner"] == "decipher":
            subprocess.call(["Rscript", "/opt/src/align.R", fasta_file,
                             alignment_file])
        alignments = fasta_parser(alignment_file)
        ref_seq = alignments["ref"]
        alignment_to_genomic = {0: contig_start - 1}
        insertion_count = 0
        for i in range(len(ref_seq)):
            if ref_seq[i] != "-":
                alignment_to_genomic[i+1] = i + contig_start - insertion_count
            else:
                insertion_count += 1
        genomic_to_alignment = {}
        for alignment_position in alignment_to_genomic:
            genomic_to_alignment[alignment_to_genomic[
                alignment_position]] = alignment_position

        def get_hap_start_index(row):
            hid = row["haplotype_ID"]
            cop = row["Copy"]
            hap_start = row["capture_start"] - 1
            hap_start_index = genomic_to_alignment[hap_start]
            hap_mip = row["MIP"]
            alignment_header = ":".join([hap_mip, cop, hid])
            hap_al = alignments[alignment_header][:hap_start_index]
            ref_al = alignments["ref"][:hap_start_index]
            diff = ref_al.count("-") - hap_al.count("-")
            return hap_start_index - diff

        contig_haps["haplotype_start_index"] = contig_haps.apply(
            get_hap_start_index, axis=1)

        raw_vcf_file = os.path.join(wdir, contig_name + ".raw.vcf")
        if contig_dict["msa_to_vcf"] == "miptools":
            msa_to_vcf(alignment_file, raw_vcf_file, ref="ref",
                       snp_only=contig_dict["snp_only"])
        else:
            subprocess.call(
                ["java", "-jar", "/opt/programs/jvarkit/dist/msa2vcf.jar",
                 "-m", "-c", "ref", "-o", raw_vcf_file, alignment_file])
        contig_dict["raw_vcf_file"] = raw_vcf_file
        # find  comment line number
        with open(raw_vcf_file) as infile:
            line_count = 0
            for line in infile:
                if line.startswith("##"):
                    line_count += 1
                else:
                    break
        vcf = pd.read_table(raw_vcf_file, skiprows=line_count)
        if vcf.empty:
            return contig_name + "_empty"
        vcf = vcf.drop(["ID", "QUAL", "FILTER", "INFO", "FORMAT"],
                       axis=1).set_index(["#CHROM", "POS", "REF", "ALT"])
        vcf = vcf.applymap(lambda a: 0 if a == "." else int(a.split(":")[0]))
        vcf = vcf.reset_index()
        vcf["alignment_position"] = vcf["POS"]
        vcf["POS"] = vcf["alignment_position"].map(alignment_to_genomic)
        vcf["CHROM"] = chrom
        vcf.drop("#CHROM",  inplace=True, axis=1)
        vcf = vcf.set_index(["CHROM", "POS", "REF", "ALT",
                             "alignment_position"])
        drop_seqs = ["ref"] + list(map(str, targets_dict.keys()))
        vcf.drop(drop_seqs, axis=1, inplace=True)
        vcf_stack = pd.DataFrame(vcf.stack()).reset_index()
        vcf_stack.rename(
            columns={"level_5": "alignment_header", 0: "genotype"},
            inplace=True)
        vcf_stack[["MIP", "Copy", "haplotype_ID"]] = vcf_stack[
            "alignment_header"].apply(lambda a: pd.Series(a.split(":")))
        vcf_merge = vcf_stack.merge(
            contig_haps[["MIP", "Copy", "haplotype_ID",
                         "capture_start", "capture_end",
                         "haplotype_start_index"]])
        vcf_merge["END"] = vcf_merge["REF"].apply(len) + vcf_merge["POS"] - 1
        vcf_merge["covered"] = (
            (vcf_merge["capture_start"] - 30 <= vcf_merge["END"])
            & (vcf_merge["capture_end"] + 30 >= vcf_merge["POS"]))
        vcf_merge.loc[~vcf_merge["covered"], "genotype"] = np.nan
        vcf_clean = vcf_merge.loc[~vcf_merge["genotype"].isnull()]
        if vcf_clean.empty:
            return contig_name + "_empty"
        contig_seq = pd.DataFrame(contig_haps.groupby("haplotype_ID")[
            "forward_sequence"].first()).to_dict(orient="index")

        def get_variant_index(row):
            pos_index = row["alignment_position"]
            hap_start_index = row["haplotype_start_index"]
            hap_copy = row["Copy"]
            hid = row["haplotype_ID"]
            hap_mip = row["MIP"]
            alignment_header = ":".join([hap_mip, hap_copy, hid])
            hap_al = alignments[alignment_header]
            hap_al = hap_al[hap_start_index:pos_index]
            variant_index = len(hap_al) - hap_al.count("-") - 1
            alts = [row["REF"]]
            alts.extend(row["ALT"].split(","))
            gen = int(row["genotype"])
            alt = alts[gen]
            variant_end_index = variant_index + len(alt)
            if variant_index < 0:
                variant_index = 0
            if variant_end_index < 1:
                variant_end_index = 1
            seq = contig_seq[hid]["forward_sequence"]
            var_seq = seq[variant_index:variant_end_index]
            return pd.Series([variant_index, variant_end_index, alt, var_seq])

        vcf_clean[
            ["variant_index", "variant_end_index", "allele", "variant"]
        ] = vcf_clean.apply(get_variant_index, axis=1)

        contig_counts_file = contig_dict["contig_counts_file"]
        contig_counts = pd.read_csv(contig_counts_file)
        contig_counts["forward_sequence_quality"] = contig_counts[
            "sequence_quality"]
        reverse_index = contig_counts["orientation"] == "reverse"
        contig_counts.loc[reverse_index, "forward_sequence_quality"] = (
            contig_counts.loc[reverse_index, "forward_sequence_quality"].apply(
                lambda a: a[::-1]))
        combined_vcf = vcf_clean[
            ["CHROM", "POS", "REF", "ALT", "genotype",
             "MIP", "Copy", "haplotype_ID", "variant_index",
             "variant_end_index"]].merge(contig_counts[
                 ["Sample ID", "haplotype_ID", "MIP", "Copy",
                  "Barcode Count", "forward_sequence_quality"]])

        def get_variant_quality(row):
            start_index = row["variant_index"]
            end_index = row["variant_end_index"]
            qual = row["forward_sequence_quality"]
            if end_index > len(qual) - 1:
                end_index = len(qual) - 1
            qual_scores = [ord(qual[i]) - 33 for i in
                           range(start_index, end_index)]
            return np.mean(qual_scores)

        combined_vcf["variant_quality"] = combined_vcf.apply(
            get_variant_quality, axis=1)

        min_count = contig_dict["min_count"]
        if min_count < 1:
            min_count = 1
        min_depth = contig_dict["min_coverage"]
        if min_depth < 1:
            min_depth = 1
        min_wsaf = contig_dict["min_wsaf"]
        if min_wsaf == 0:
            min_wsaf = 0.0001

        def collapse_vcf(group):
            key = group.iloc[0][["CHROM", "POS", "REF", "ALT"]].values
            alts = key[3].split(",")
            allele_count = len(alts) + 1
            allele_depths = []
            for i in range(allele_count):
                allele_depths.append(group.loc[group["genotype"] == i,
                                               "Barcode Count"].sum().round(0))
            total_depth = int(round(np.sum(allele_depths), 0))
            wsaf = np.array(allele_depths)/total_depth
            if total_depth < min_depth:
                return nastring
            genotypes = []
            for i in range(allele_count):
                if (allele_depths[i] >= min_count) and (wsaf[i] >= min_wsaf):
                    genotypes.append(i)
            if len(genotypes) == 0:
                return nastring
            else:
                alleles = list(range(allele_count))
                geno = sorted(zip(alleles, allele_depths),
                              key=itemgetter(1, 0), reverse=True)[:2]
                if len(genotypes) == 1:
                    gt = str(geno[0][0])
                    gt = gt + "/" + gt
                else:
                    gt1 = geno[0][0]
                    gt2 = geno[1][0]
                    gt = sorted(map(str, [gt1, gt2]))
                    gt = "/".join(gt)
            allele_depths = [str(int(a)) for a in allele_depths]
            variant_quals = []
            for i in range(allele_count):
                variant_quals.append(group.loc[group["genotype"] == i,
                                               "variant_quality"].max())
            variant_quals = ["." if np.isnan(v) else str(int(round(v, 0)))
                             for v in variant_quals]
            mip_count = []
            for i in range(allele_count):
                mip_count.append(len(set(group.loc[group["genotype"] == i,
                                                   "MIP"])))
            hap_count = []
            for i in range(allele_count):
                hap_count.append(len(set(group.loc[group["genotype"] == i,
                                                   "haplotype_ID"])))
            return ":".join([gt, ",".join(allele_depths),
                             str(total_depth),
                             ",".join(variant_quals),
                             ",".join(map(str, mip_count)),
                             ",".join(map(str, hap_count)),
                             ",".join(map(str, wsaf.round(3))),
                             ])

        collapsed_vcf = pd.DataFrame(combined_vcf.groupby(
            ["CHROM", "POS", "REF", "ALT", "Sample ID"]).apply(collapse_vcf)
        ).reset_index()
        vcf_table = collapsed_vcf.pivot_table(
            index=["CHROM", "POS", "REF", "ALT"],
            columns="Sample ID", aggfunc="first")
        vcf_table.fillna(nastring, inplace=True)

        def get_var_summary(row):
            val = row.values
            ad = []
            quals = []
            wsafs = []
            mip_counts = []
            hap_counts = []
            genotypes = []
            for v in val:
                if v != nastring:
                    genotypes.append(v.split(":")[0])
                    ad.append(list(map(int, v.split(":")[1].split(","))))
                    quals.append(v.split(":")[3].split(","))
                    mip_counts.append(list(map(
                        int, v.split(":")[4].split(","))))
                    hap_counts.append(list(map(
                        int, v.split(":")[5].split(","))))
                    wsafs.append(list(map(float, v.split(":")[6].split(","))))
            if len(ad) == 0:
                return "."
            geno_dict = {}
            an_count = 0
            for geno in genotypes:
                try:
                    geno_list = list(map(int, geno.split("/")))
                    for gt in geno_list:
                        try:
                            geno_dict[gt] += 1
                        except KeyError:
                            geno_dict[gt] = 1
                        an_count += 1
                except ValueError:
                    continue
            number_of_alleles = len(ad[0])
            ac_list = []
            for i in range(number_of_alleles):
                try:
                    ac_list.append(geno_dict[i])
                except KeyError:
                    ac_list.append(0)

            quality = []
            for q in quals:
                nq = []
                for q_val in q:
                    if q_val == ".":
                        nq.append(np.nan)
                    else:
                        nq.append(int(q_val))
                quality.append(nq)
            quals = np.nanmean(quality, axis=0)
            quality = []
            for q in quals:
                if np.isnan(q):
                    quality.append(".")
                else:
                    quality.append(str(round(q, 1)))

            wsafs = pd.DataFrame(wsafs)
            wsafs = wsafs.applymap(
                lambda a: a if a >= min_wsaf else np.nan).mean().round(4)
            wsafs = wsafs.fillna(0).astype(str)

            mip_counts = pd.DataFrame(mip_counts)
            mip_counts = mip_counts.applymap(
                lambda a: a if a > 0 else np.nan).mean().round(2)
            mip_frac = (mip_counts / (mip_counts.max())).round(2)
            mip_frac = mip_frac.fillna(0).astype(str)
            mip_counts = mip_counts.fillna(0).astype(str)

            hap_counts = pd.DataFrame(hap_counts)
            hap_counts = hap_counts.applymap(
                lambda a: a if a > 0 else np.nan).mean().round(2)
            hap_counts = hap_counts.fillna(0).astype(str)

            info_cols = [
                "DP=" + str(np.sum(ad)),
                "AD=" + ",".join(map(str, np.sum(ad, axis=0))),
                "AC=" + ",".join(map(str, ac_list[1:])),
                "AN=" + str(an_count),
                "AF=" + ",".join(map(str, (
                    np.array(ac_list)/an_count)[1:].round(4))),
                "RAF=" + ",".join(map(str, (
                    np.array(ac_list)/an_count).round(4))),
                "RAC=" + ",".join(map(str, ac_list)),
                "NS=" + str(len(ad)),
                "SC=" + ",".join(map(str, (np.array(ad) >= min_count).sum(
                    axis=0))),
                "SF=" + ",".join(map(str, ((np.array(ad) >= min_count).sum(
                    axis=0)/len(ad)).round(5))),
                "QS=" + ",".join(quality),
                "WSAF=" + ",".join(wsafs),
                "MC=" + ",".join(mip_counts),
                "MCF=" + ",".join(mip_frac),
                "HC=" + ",".join(hap_counts)]

            variant_pos = row.name[1]
            ref_len = len(row.name[2])
            variant_end = variant_pos + ref_len - 1
            overlapping_targets = set()
            for p in target_pos:
                ol = overlap([variant_pos, variant_end],
                             [p["Pos"], p["End"]])
                if len(ol) > 0:
                    overlapping_targets.add(p["Mutation Name"])
            if len(overlapping_targets) > 0:
                ot_field = ",".join(sorted(overlapping_targets))
                info_cols.append("OT=" + ot_field)

            return ";".join(info_cols)

        var_summary = pd.DataFrame(vcf_table.apply(
            get_var_summary, axis=1)).rename(columns={0: "INFO"})
        var_summary["FORMAT"] = "GT:AD:DP:QS:MC:HC:WSAF"
        var_summary["ID"] = "."
        var_summary["QUAL"] = "."
        var_summary["FILTER"] = "."
        samples = vcf_table.columns.droplevel(0).tolist()
        vcf_table.columns = samples
        samples = contig_dict["sample_ids"]
        vcf_table = vcf_table.loc[:, samples].fillna(nastring)
        vcf_table = vcf_table.merge(var_summary, left_index=True,
                                    right_index=True)
        vcf_table = vcf_table.reset_index()[
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
             "FORMAT"] + samples]
        vcf_table.rename(columns={"CHROM": "#CHROM"}, inplace=True)
        vcf_table = vcf_table.sort_values("POS")
        vcf_header = [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description='
            '"Total coverage for locus, across samples.">',
            "##INFO=<ID=AD,Number=R,Type=Integer,Description="
            '"Total coverage per allele, across samples.">',
            "##INFO=<ID=AC,Number=A,Type=Integer,Description="
            '"Total number of alternate alleles in called genotypes.">',
            "##INFO=<ID=AN,Number=1,Type=Integer,Description="
            '"Total number of alleles in called genotypes.">',
            "##INFO=<ID=AF,Number=A,Type=Float,Description="
            '"Allele frequency (AC/AN) for alternate alleles.">',
            "##INFO=<ID=RAF,Number=R,Type=Float,Description="
            '"Allele frequency (AC/AN) for all alleles.">',
            "##INFO=<ID=RAC,Number=R,Type=Integer,Description="
            '"Total number of each allele in called genotypes.">',
            "##INFO=<ID=NS,Number=1,Type=Integer,Description="
            '"Number of samples with genotype calls.">',
            "##INFO=<ID=SC,Number=R,Type=Integer,Description="
            '"Number of samples carrying the allele.">',
            "##INFO=<ID=SF,Number=R,Type=Float,Description="
            '"Frequency of samples carrying the allele.">',
            "##INFO=<ID=QS,Number=R,Type=Float,Description="
            '"Average sequence quality per allele.">',
            "##INFO=<ID=WSAF,Number=R,Type=Float,Description="
            '"Average nonzero WithinSampleAlleleFrequency.">',
            "##INFO=<ID=MC,Number=R,Type=Float,Description="
            '"Average number of MIPs supporting the allele (when called).">',
            "##INFO=<ID=MCF,Number=R,Type=Float,Description="
            '"MC expressed as the fraction of MAX MC.">',
            "##INFO=<ID=HC,Number=R,Type=Float,Description="
            '"Average number of haplotypes supporting the allele'
            ' (when called).">',
            "##INFO=<ID=OT,Number=.,Type=String,Description="
            '"Variant position overlaps with a target.">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
            '"Number of observation for each allele.">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
            '"Total read depth (coverage) at this position">',
            '##FORMAT=<ID=QS,Number=R,Type=Integer,Description='
            '"Sequence quality per allele.">',
            '##FORMAT=<ID=MC,Number=R,Type=Integer,Description='
            '"Number of MIPs supporting the allele.">',
            '##FORMAT=<ID=HC,Number=R,Type=Integer,Description='
            '"Number of haplotypes supporting the allele.">',
            '##FORMAT=<ID=WSAF,Number=R,Type=Float,Description='
            '"Within sample allele frequency.">']
        # save vcf file
        contig_vcf_file = os.path.join(wdir, contig_name + ".vcf")
        with open(contig_vcf_file, "w") as outfile:
            outfile.write("\n".join(vcf_header) + "\n")
            vcf_table.to_csv(outfile, index=False, sep="\t")
        contig_variants_file = os.path.join(wdir,
                                            contig_name + "_variants.csv")
        combined_vcf.to_csv(contig_variants_file)
        collapsed_variants_file = os.path.join(wdir, contig_name
                                               + "_collapsed_variants.csv")
        collapsed_vcf.to_csv(collapsed_variants_file)
        contig_haps.to_csv(contig_haplotypes_file)
        contig_counts.to_csv(contig_counts_file)
        return contig_name
    except Exception as e:
        return ExceptionWrapper(e)


###############################################################################
# general use functions.
###############################################################################

def parse_alignment_positions(alignment_file, contig_start, ref_key="ref"):
    """ Parse a multiple sequence alignment file given in fasta format.
    Using the genomic start position of the reference contig, create a
    genome to alignment and alignment to genome position maps.
    """
    alignments = fasta_parser(alignment_file)
    ref_seq = alignments[ref_key]
    alignment_to_genomic = {0: contig_start - 1}
    insertion_count = 0
    for i in range(len(ref_seq)):
        if ref_seq[i] != "-":
            alignment_to_genomic[i+1] = i + contig_start - insertion_count
        else:
            insertion_count += 1
    genomic_to_alignment = {}
    for alignment_position in alignment_to_genomic:
        genomic_to_alignment[alignment_to_genomic[
            alignment_position]] = alignment_position
    return {"a2g": alignment_to_genomic, "g2a": genomic_to_alignment}


def check_overlap(r1, r2, padding=0):
    """ Check if two regions overlap. Regions are given as lists of chrom (str),
    begin (int), end (int)."""
    # check chromosome equivalency
    o1 = r1[0] == r2[0]
    # check interval overlap
    merged = merge_overlap([r1[1:], r2[1:]], padding)
    o2 = len(merged) == 1
    return o1 & o2


def make_region(chromosome, begin, end):
    """ Create region string from coordinates.
    takes 2 (1 for human 1-9) digit chromosome,
    begin and end positions (1 indexed)"""
    region = "chr" + str(chromosome) + ":" + str(begin) + "-" + str(end)
    return region


def create_region(chromosome, begin, end):
    """ Create region string from coordinates.
    chromosome string,
    begin and end positions (1 indexed)"""
    region = chromosome + ":" + str(begin) + "-" + str(end)
    return region


def get_coordinates(region):
    """ Define coordinates chr, start pos and end positions
    from region string chrX:start-end. Return coordinate list.
    """
    chromosome = region.split(":")[0]
    coord = region.split(":")[1]
    coord_list = coord.split("-")
    begin = int(coord_list[0])
    end = int(coord_list[1])
    return [chromosome, begin, end]


def get_fasta(region, species="pf", offset=1, header="na"):
    """ Take a region string (chrX:begin-end (1 indexed)),
    and species (human=hs, plasmodium= pf),Return fasta record.
    """
    if offset == 0:
        region_coordinates = get_coordinates(region)
        region = (region_coordinates[0] + ":" + str(region_coordinates[1] + 1)
                  + "-" + str(region_coordinates[2]))
    region = region.encode("utf-8")
    file_locations = get_file_locations()
    genome_fasta = file_locations[species]["fasta_genome"].encode("utf-8")
    fasta = pysam.faidx(genome_fasta, region)
    if header != "na":
        fasta_seq = "\n".join(fasta.split("\n")[1:])
        fasta = ">" + header + "\n" + fasta_seq
    return fasta


def get_fasta_list(regions, species):
    """ Take a list of regions and return fasta sequences."""
    if len(regions) == 0:
        return {}
    file_locations = get_file_locations()
    genome_fasta = file_locations[species]["fasta_genome"]
    region_file = "/tmp/region s_" + id_generator(10) + ".txt"
    with open(region_file, "w") as outfile:
        for r in regions:
            outfile.write(r + "\n")
    fasta_dic = {}
    command = ["samtools",  "faidx", "-r", region_file, genome_fasta]
    out = subprocess.check_output(command).decode("UTF-8")
    fasta_list = out.split(">")[1:]
    for f in fasta_list:
        fl = f.strip().split("\n")
        fhead = fl[0]
        fseq = "".join(fl[1:])
        fasta_dic[fhead] = fseq
    return fasta_dic


def create_fasta_file(region, species, output_file):
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    with open(output_file, "a") as outfile:
        outfile.write(get_fasta(region, species))


def merge_overlap(intervals, spacer=0):
    """Merge overlapping intervals.

    Take a list of lists of 2 elements, [start, stop],
    check if any [start, stop] pairs overlap and merge if any.
    Return the merged [start, stop] list.
    """
    # reuse a piece of code from get_exons:
    #######################################
    exons = copy.deepcopy(intervals)
    exons = [e for e in exons if len(e) == 2]
    for e in exons:
        e.sort()
    exons.sort()
    if len(exons) < 2:
        return exons
    overlapping = 1
    while overlapping:
        overlapping = 0
        for i in range(len(exons)):
            e = exons[i]
            for j in range(len(exons)):
                x = exons[j]
                if i == j:
                    continue
                else:
                    if e[1] >= x[1]:
                        if (e[0] - x[1]) <= spacer:
                            overlapping = 1
                    elif x[1] >= e[1]:
                        if (x[0] - e[1]) <= spacer:
                            overlapping = 1
                    if overlapping:
                        # merge exons and add to the exon list
                        exons.append([min(e[0], x[0]), max(e[1], x[1])])
                        # remove the exons e and x
                        exons.remove(e)
                        exons.remove(x)
                        # break once an overlapping exon is found
                        break
            if overlapping:
                # if an overlapping exon is found,
                # stop this for loop and continue with the
                # while loop with the updated exon list
                break
    exons.sort()
    return exons


def overlap(reg1, reg2):
    """
    Return overlap between two regions.
    e.g. [10, 30], [20, 40] returns [20, 30]
    """
    try:
        intersect = set(range(reg1[0], reg1[1] + 1)).intersection(
            set(range(reg2[0], reg2[1] + 1)))
        intersect = sorted(intersect)
        return [intersect[0]] + [intersect[-1]]
    except IndexError:
        return []


def remove_overlap(reg1, reg2, spacer=0):
    """
    Remove overlap between two regions.
    e.g. [10, 30], [20, 40] returns [10, 20], [30, 40]
    """
    regions = sorted([sorted(reg1), sorted(reg2)])
    try:
        if regions[0][1] - regions[1][0] >= spacer:
            coords = sorted(reg1 + reg2)
            return[[coords[0], coords[1] - 1],
                   [coords[2] + 1, coords[3]]]
        else:
            return regions
    except IndexError:
        return []


def complete_overlap(reg1, reg2):
    """
    Return whether one of the two given regions contain the other.
    e.g. [10, 40], [20, 30] returns True.
    """
    regions = sorted([sorted(reg1), sorted(reg2)])
    try:
        return (((regions[0][0] == regions[1][0])
                 and (regions[0][1] <= regions[1][1]))
                or ((regions[0][0] < regions[1][0])
                    and (regions[0][1] >= regions[1][1])))
    except IndexError:
        return False


def check_redundant_region(reg1, reg2, spacer=0):
    """
    Return whether one of the two given regions is redundant.
    i.e. one contains the other or there is less than 'spacer'
    non-overlap between them.
    """
    regions = sorted([sorted(reg1), sorted(reg2)])
    try:
        if complete_overlap(*regions):
            return True
        else:
            non_overlap = remove_overlap(*regions)
            if len(non_overlap) == 0:
                return False
            extra = sum([r[1] - r[0] + 1 for r in non_overlap])
            if extra <= spacer:
                return True
            else:
                return False
    except IndexError:
        return False


def subtract_overlap(uncovered_regions, covered_regions, spacer=0):
    """
    Given two sets of regions in the form [[start, end], [start, end]],
    return a set of regions that is the second set subtracted from the first.
    """
    uncovered_set = set()
    for r in uncovered_regions:
        try:
            uncovered_set.update(list(range(r[0], r[1] + 1)))
        except IndexError:
            pass
    covered_set = set()
    for r in covered_regions:
        try:
            covered_set.update(list(range(r[0], r[1] + 1)))
        except IndexError:
            pass
    uncovered_remaining = sorted(uncovered_set.difference(covered_set))
    if len(uncovered_remaining) > 0:
        uncovered = [[uncovered_remaining[i-1], uncovered_remaining[i]]
                     for i in range(1, len(uncovered_remaining))
                     if uncovered_remaining[i] - uncovered_remaining[i-1] > 1]
        unc = [uncovered_remaining[0]]
        for u in uncovered:
            unc.extend(u)
        unc.append(uncovered_remaining[-1])
        return [[unc[i], unc[i+1]]for i in range(0, len(unc), 2)
                if unc[i+1] - unc[i] > spacer]
    else:
        return []


def trim_overlap(region_list, low=0.1, high=0.9, spacer=0):
    """
    Given a set of regions in the form [[start, end], [start, end]],
    return a set of regions with any overlapping parts trimmed
    when overlap size / smaller region size ratio is lower than "low";
    or flanking region outside of overlap is trimmed when the ratio
    is higher than "high".
    """
    do_trim = True
    while do_trim:
        do_trim = False
        break_for = False
        region_list = [r for r in region_list if r != "remove"]
        for i in range(len(region_list)):
            if break_for:
                break
            else:
                for j in range(len(region_list)):
                    if i != j:
                        reg_i = region_list[i]
                        reg_j = region_list[j]
                        if reg_i == reg_j:
                            region_list[i] = "remove"
                            break_for = True
                            do_trim = True
                            break
                        else:
                            overlapping_region = overlap(reg_i, reg_j)
                            if len(overlapping_region) > 0:
                                reg_sizes = sorted([reg_i[1] - reg_i[0] + 1,
                                                    reg_j[1] - reg_j[0] + 1])
                                overlap_size = float(overlapping_region[1]
                                                     - overlapping_region[0])
                                overlap_ratio = overlap_size/reg_sizes[0]
                                if overlap_ratio <= low:
                                    region_list[i] = "remove"
                                    region_list[j] = "remove"
                                    region_list.extend(remove_overlap(
                                        reg_i, reg_j, spacer))
                                    break_for = True
                                    do_trim = True
                                    break
                                elif overlap_ratio >= high:
                                    region_list[i] = "remove"
                                    region_list[j] = "remove"
                                    region_list.append(overlapping_region)
                                    break_for = True
                                    do_trim = True
                                    break
                                else:
                                    print(overlap_ratio,
                                          "is outside trim range for ",
                                          reg_i, reg_j)

    return region_list


def fasta_parser(fasta_file, use_description=False):
    """Convert a fasta file to python dict.

    Convert a fasta file with multiple sequences to a dictionary with fasta
    id as keys and sequences as values. The fasta id is the text in the fasta
    header before the first space character. If the entire header line is
    to be used, use_description=True should be passed.
    """
    fasta_dic = {}
    records = SeqIO.parse(fasta_file, format="fasta")
    for rec in records:
        if use_description:
            header = rec.description
        else:
            header = rec.id
        if header in fasta_dic:
            print(("%s occurs multiple times in fasta file" % header))
        fasta_dic[header] = str(rec.seq)
    return fasta_dic


def fasta_parser_verbatim(fasta):
    """Convert a fasta file with multiple sequences to a dictionary.

    Convert a fasta file with multiple sequences to a dictionary with fasta
    headers as keys and sequences as values. Spaces are allowed in keys.
    """
    fasta_dic = {}
    with open(fasta) as infile:
        for line in infile:
            # find the headers
            if line.startswith(">"):
                header = line[1:-1]
                if header in fasta_dic:
                    print(("%s occurs multiple times in fasta file" % header))
                fasta_dic[header] = ""
                continue
            try:
                fasta_dic[header] = fasta_dic[header] + line.strip()
            except KeyError:
                fasta_dic[header] = line.strip()
    return fasta_dic


def fasta_to_sequence(fasta):
    """ Convert a multiline fasta sequence to one line sequence"""
    f = fasta.strip().split("\n")
    if len(f) > 0:
        return "".join(f[1:])
    else:
        return ""


def get_sequence(region, species):
    return fasta_to_sequence(get_fasta(region, species))


def unmask_fasta(masked_fasta, unmasked_fasta):
    """ Unmask lowercased masked fasta file, save """
    with open(masked_fasta) as infile, open(unmasked_fasta, "w") as outfile:
        for line in infile:
            if not line.startswith((">", "#")):
                outfile.write(line.upper())
            else:
                outfile.write(line)
    return


def fasta_to_fastq(fasta_file, fastq_file):
    """ Create a fastq file from fasta file with dummy quality scores."""
    fasta = fasta_parser(fasta_file)
    fastq_list = []
    for f in fasta:
        fastq_list.append("@" + f)
        fastq_list.append(fasta[f])
        fastq_list.append("+")
        fastq_list.append("H" * len(fasta[f]))
    with open(fastq_file, "w") as outfile:
        outfile.write("\n".join(fastq_list))
    return


def combine_sample_data(gr):
    """Combine data from multiple sequencing runs for the same sample.

    Take a pandas groupby object representing multiple data points
    corresponding the same sequence and sample, from multiple sequence runs.
    Sum the barcode and read counts for the combined result. Use the sequencing
    quality values for the record with most supporting barcodes.

    Return a single combined record in pd.Series object so that all results can
    be combined into a new pd.DataFrame for all samples.
    """
    result = {}
    result["barcode_count"] = gr["barcode_count"].sum()
    result["read_count"] = gr["read_count"].sum()
    result["sequence_quality"] = gr.sort_values(
        "barcode_count",
        ascending=False
    )["sequence_quality"].iloc[0]
    result["mip_name"] = gr["mip_name"].iloc[0]
    result["gene_name"] = gr["gene_name"].iloc[0]
    return pd.Series(result)


def combine_info_files(wdir,
                       settings_file,
                       info_files,
                       sample_sheets,
                       combined_file,
                       sample_sets=None):
    """Combine MIPWrangler outputs from multiple runs."""
    settings = get_analysis_settings(os.path.join(wdir, settings_file))
    colnames = dict(list(zip(settings["colNames"],
                         settings["givenNames"])))
    c_keys = list(colnames.keys())
    c_vals = [colnames[k] for k in c_keys]
    data = []
    run_meta = []
    for i in range(len(sample_sheets)):
        current_run_meta = pd.read_table(sample_sheets[i])
        for k in ["sample_name", "sample_set", "replicate"]:
            current_run_meta[k] = current_run_meta[k].astype(str)
        current_run_meta["sheet_order"] = i
        current_run_meta["Original SID"] = current_run_meta[
            ["sample_name", "sample_set", "replicate"]
        ].apply(lambda a: "-".join(a), axis=1)
        run_meta.append(current_run_meta)
    run_meta = pd.concat(run_meta, ignore_index=True)
    if sample_sets is not None:
        sps = pd.DataFrame(sample_sets, columns=["sample_set",
                                                 "probe_set"])
    else:
        sps = run_meta.groupby(
            ["sample_set", "probe_set"]
        ).first().reset_index()[["sample_set", "probe_set"]]
    run_meta = run_meta.merge(sps, how="inner")
    run_meta.rename(columns={"library_prep": "Library Prep"}, inplace=True)
    run_meta_collapsed = run_meta.groupby(
        ["sample_name", "sample_set", "replicate", "Library Prep"]
    ).first().reset_index()[["sample_name", "sample_set",
                             "replicate", "Library Prep"]]
    # check if there are repeating sample_name, sample_set, replicate
    # combinations; which make up the sample ID. If there are, replicate
    # numbers will need to be re-assigned so that each library has a unique
    # ID. If no overlap, they should be left as they are.
    repeat_found = False
    # check if replicates are to be ignored, i.e. merge all libraries from
    # the same DNA source.
    merge_replicates = False
    try:
        if int(settings["mergeReplicates"]):
            merge_replicates = True
    except KeyError:
        pass

    def assign_replicate(replicates):
        replicates = list(map(int, replicates))
        group_size = len(replicates)
        reps_available = set(range(1, group_size + 1))
        reps_used = set(replicates)
        reps_available = reps_available.difference(reps_used)
        reps_available = sorted(reps_available, reverse=True)
        reps_used = set()
        for i in range(group_size):
            rep = replicates[i]
            if np.isnan(rep) or (rep in reps_used):
                rep = int(reps_available.pop())
                replicates[i] = rep
                reps_used.add(rep)
                # print a warning unless the replicates will eventually
                # be merged
                if not (repeat_found or merge_replicates):
                    repeat_found
                    print("Sample ID will change for a sample because there "
                          "was another sample with the same ID. Please check "
                          "the samples.tsv file to compare SID and Original "
                          "SID fields.")
            else:
                replicates[i] = int(rep)
                reps_used.add(rep)
        return pd.Series(replicates)

    run_meta_collapsed["new_replicate"] = run_meta_collapsed.groupby(
        ["sample_name", "sample_set"])["replicate"].transform(
        assign_replicate).astype(str)
    run_meta = run_meta.merge(run_meta_collapsed)
    run_meta["Sample ID"] = run_meta[["sample_name",
                                      "sample_set",
                                      "new_replicate"]].apply(
        lambda a: "-".join(a), axis=1
    )
    # load the probe set dictionary to extract the
    # probes that we're interested in
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    for i in range(len(info_files)):
        i_file = info_files[i]
        current_run_meta = run_meta.loc[run_meta["sheet_order"] == i]
        current_run_dict = current_run_meta.set_index(
            "Original SID"
        ).to_dict(orient="index")
        line_number = 0
        try:
            gzip.open(i_file, "rb").readline()
            inf_file = gzip.open(i_file, "rb")
        except IOError:
            inf_file = open(i_file, "rb")
        with inf_file as infile:
            for line in infile:
                newline = line.decode("utf-8").strip().split("\t")
                line_number += 1
                if line_number == 1:
                    col_indexes = [
                        newline.index(ck)
                        for ck in c_keys
                    ]
                    for ci in col_indexes:
                        if colnames[newline[ci]] == "sample_name":
                            si_index = ci
                        elif colnames[newline[ci]] == "mip_name":
                            mip_name_index = ci
                else:
                    ori_sample_id = newline[si_index]
                    mip_fam_name = newline[mip_name_index]
                    if mip_fam_name in used_probes:
                        try:
                            library = current_run_dict[
                                ori_sample_id
                            ]["Library Prep"]
                            sample_id = current_run_dict[
                                ori_sample_id
                            ]["Sample ID"]
                            d = ([newline[ci] if ci != si_index else sample_id
                                  for ci in col_indexes] + [library])
                            data.append(d)
                        except KeyError:
                            continue
    info = pd.DataFrame(data, columns=c_vals + ["Library Prep"])
    info["barcode_count"] = info["barcode_count"].astype(int)
    info["read_count"] = info["read_count"].astype(int)
    # check if replicates are to be ignored, i.e. merge all libraries from
    # the same DNA source.
    if merge_replicates:
        info["original_sample_name"] = info["sample_name"]
        info["sample_name"] = info["sample_name"].apply(
            lambda a: "-".join(a.split("-")[:-1]) + "-1")
        info["Library Prep"] = "merged"
    info = info.groupby(
        ["sample_name", "haplotype_sequence", "Library Prep"]
    ).apply(combine_sample_data).reset_index()
    m_groups = info.groupby("mip_name")
    h_list = []
    for m, g in m_groups:
        md = pd.DataFrame(g.groupby(["mip_name",
                                    "haplotype_sequence"]).size().sort_values(
            ascending=False
        ).reset_index()).reset_index()
        md["index"] = md["index"].astype(str)
        md["haplotype_ID"] = md["mip_name"] + "." + md["index"]
        h_list.append(md[["haplotype_sequence", "haplotype_ID"]])
    hap_ids = pd.concat(h_list, ignore_index=True)
    info = info.merge(hap_ids)
    info.to_csv(os.path.join(wdir, combined_file), index=False, sep="\t")
    info.groupby(["gene_name", "mip_name", "haplotype_ID"])[
        "haplotype_sequence"].first().reset_index().to_csv(
            os.path.join(wdir, "unique_haplotypes.csv"), index=False)
    run_meta = run_meta.groupby("Sample ID").first().reset_index()
    run_meta = run_meta.drop(["Sample ID",
                              "sheet_order",
                              "replicate"],
                             axis=1).rename(
        columns={"new_replicate": "replicate"}
    )
    if merge_replicates:
        run_meta["replicate"] = 1
    run_meta.to_csv(os.path.join(wdir, "samples.tsv"), sep="\t", index=False)


def process_info_file(wdir,
                      settings_file,
                      info_files,
                      sample_sheets,
                      combined_file,
                      sample_sets=None):
    """
    Process MIPWrangler output file.

    This function extracts the relevant fields from a given MIPWrangler
    output file, renames the columns to be used in downstream analysis and
    merges the provided meta data.
    """
    settings = get_analysis_settings(os.path.join(wdir, settings_file))
    colnames = dict(list(zip(settings["colNames"],
                         settings["givenNames"])))
    c_keys = list(colnames.keys())
    c_vals = [colnames[k] for k in c_keys]
    data = []
    current_run_meta = pd.read_table(sample_sheets[0])
    for k in ["sample_name", "sample_set", "replicate"]:
        current_run_meta[k] = current_run_meta[k].astype(str)
    current_run_meta["sheet_order"] = 0
    current_run_meta["Original SID"] = current_run_meta[
        ["sample_name", "sample_set", "replicate"]
    ].apply(lambda a: "-".join(a), axis=1)
    run_meta = current_run_meta
    run_meta.rename(columns={"library_prep": "Library Prep"}, inplace=True)
    if sample_sets is not None:
        sps = pd.DataFrame(sample_sets, columns=["sample_set",
                                                 "probe_set"])
    else:
        sps = run_meta.groupby(
            ["sample_set", "probe_set"]
        ).first().reset_index()[["sample_set", "probe_set"]]
    run_meta = run_meta.merge(sps, how="inner")
    run_meta["Sample ID"] = run_meta["Original SID"]
    # load the probe set dictionary to extract the
    # probes that we're interested in
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    i_file = info_files[0]
    current_run_meta = run_meta
    current_run_dict = current_run_meta.set_index(
        "Original SID"
    ).to_dict(orient="index")
    line_number = 0
    try:
        gzip.open(i_file, "rb").readline()
        inf_file = gzip.open(i_file, "rb")
    except IOError:
        inf_file = open(i_file, "rb")
    with inf_file as infile:
        for line in infile:
            newline = line.decode("utf-8").strip().split("\t")
            line_number += 1
            if line_number == 1:
                col_indexes = [
                    newline.index(ck)
                    for ck in c_keys
                ]
                for ci in col_indexes:
                    if colnames[newline[ci]] == "sample_name":
                        si_index = ci
                    elif colnames[newline[ci]] == "mip_name":
                        mip_name_index = ci
            else:
                ori_sample_id = newline[si_index]
                mip_fam_name = newline[mip_name_index]
                if mip_fam_name in used_probes:
                    try:
                        library = current_run_dict[
                            ori_sample_id
                        ]["Library Prep"]
                        sample_id = current_run_dict[
                            ori_sample_id
                        ]["Sample ID"]
                        d = ([newline[ci] if ci != si_index else sample_id
                              for ci in col_indexes] + [library])
                        data.append(d)
                    except KeyError:
                        continue
    info = pd.DataFrame(data, columns=c_vals + ["Library Prep"])
    info["barcode_count"] = info["barcode_count"].astype(int)
    info["read_count"] = info["read_count"].astype(int)
    info.to_csv(os.path.join(wdir, combined_file), index=False, sep="\t")
    info.groupby(["gene_name", "mip_name", "haplotype_ID"])[
        "haplotype_sequence"].first().reset_index().to_csv(
            os.path.join(wdir, "unique_haplotypes.csv"), index=False)
    run_meta = run_meta.groupby("Sample ID").first().reset_index()
    run_meta = run_meta.drop("Sample ID", axis=1)
    run_meta.to_csv(os.path.join(wdir, "samples.tsv"), sep="\t", index=False)


def update_probe_sets(
        mipset_table="/opt/project_resources/mip_ids/mipsets.csv",
        mipset_json="/opt/project_resources/mip_ids/probe_sets.json"):
    mipsets = pd.read_csv(mipset_table)
    mipset_list = mipsets.to_dict(orient="list")
    mipset_dict = {}
    for mipset in mipset_list:
        mlist = mipset_list[mipset]
        mipset_dict[mipset] = [m for m in mlist if not pd.isnull(m)]
    with open(mipset_json, "w") as outfile:
        json.dump(mipset_dict, outfile, indent=1)
    return


def generate_fastqs(wdir, mipster_files, min_bc_count, min_bc_frac):
    """
    Generate fastq files for each sample in raw MIPWrangler output file(s).

    These files will have stitched and barcode corrected reads.
    """
    fastq_dir = os.path.join(wdir, "fastq")
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    mipster_dfs = pd.concat([pd.read_table(os.path.join(wdir, mfile),
                                           usecols=[
                                              "s_Sample",
                                              'h_popUID',
                                              "h_seq",
                                              'c_qual',
                                              'c_barcodeCnt',
                                              "c_barcodeFrac"
                                           ])
                             for mfile in mipster_files],
                            axis=0,
                            ignore_index=True)
    mipster = mipster_dfs.loc[
        (mipster_dfs["c_barcodeCnt"] >= min_bc_count)
        & (mipster_dfs["c_barcorac"] >= min_bc_frac)
    ].groupby("s_Sample").apply(lambda x: pd.DataFrame.to_dict(
        x, orient="index"
    )).to_dict()
    for sample in mipster:
        fastq_file = os.path.join(fastq_dir, sample + ".fq.gz")
        with gzip.open(fastq_file, "wb") as outfile:
            outfile_list = []
            for ind in mipster[sample]:
                row = mipster[sample][ind]
                bc = int(row["c_barcodeCnt"])
                hid = row["h_popUID"]
                qual = row["c_qual"]
                seq = row["h_seq"]
                sample = row["s_Sample"]
                for i in range(bc):
                    read_name = "_".join(["@", sample, hid, str(ind), str(i)])
                    outfile_list.extend([read_name, seq, "+", qual])
            outfile.write(("\n".join(outfile_list) + "\n").encode("UTF-8"))
    return


def generate_processed_fastqs_worker(fastq_file, sample_mipster):
    """Worker function for generate_processed_fastqs."""
    with gzip.open(fastq_file, "wb") as outfile:
        outfile_list = []
        for ind in sample_mipster:
            row = sample_mipster[ind]
            bc = int(row["barcode_count"])
            hid = row["haplotype_ID"]
            qual = row["sequence_quality"]
            seq = row["haplotype_sequence"]
            sample = row["sample_name"]
            for i in range(bc):
                read_name = "_".join(["@", sample, hid, str(ind), str(i)])
                outfile_list.extend([read_name, seq, "+", qual])
        outfile.write(("\n".join(outfile_list) + "\n").encode("UTF-8"))


def generate_processed_fastqs(fastq_dir, mipster_file,
                              min_bc_count=1,
                              pro=8):
    """
    Generate fastq files for each sample in processed MIPWrangler output file.

    The resulting fastq files will have stitched and barcode corrected reads.
    """
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    mipster = pd.read_table(mipster_file,
                            usecols=[
                              "sample_name",
                              'haplotype_ID',
                              "haplotype_sequence",
                              'sequence_quality',
                              'barcode_count'
                            ])
    mipster = mipster.loc[mipster["barcode_count"] >= min_bc_count].groupby(
        "sample_name"
    ).apply(lambda x: pd.DataFrame.to_dict(x, orient="index")).to_dict()
    p = Pool(pro)
    for sample in mipster:
        fastq_file = os.path.join(fastq_dir, sample + ".fq.gz")
        sample_mipster = mipster[sample]
        p.apply_async(generate_processed_fastqs_worker, (fastq_file,
                                                         sample_mipster))
    p.close()
    p.join()
    return


def generate_mapped_fastqs(fastq_dir, mipster_file,
                           mapped_haplotypes_file, species, min_bc_count=1,
                           pro=8, pad_size=20, save=False):
    """
    Generate fastq files for each sample in a mapped MIPWrangler output file.

    The input is file is free from off target reads etc. The resulting fastq
    files will have stitched and barcode corrected on target reads.
    """
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    mipster = pd.read_table(mipster_file,
                            usecols=[
                              "sample_name",
                              'haplotype_ID',
                              "haplotype_sequence",
                              'sequence_quality',
                              'barcode_count'
                            ])
    mapped_haplotypes = pd.read_csv(mapped_haplotypes_file)

    def get_padded_haplotype_sequence(row, pad_size, species):
        """Pad haplotype sequence with flanking sequence from the reference."""
        chrom = row["Chrom"]
        capture_start = int(row["capture_start"])
        capture_end = int(row["capture_end"])
        pad_start = capture_start - pad_size
        pad_end = capture_end + pad_size
        left_key = create_region(chrom, pad_start, capture_start - 1)
        right_key = create_region(chrom, capture_end + 1, pad_end)
        left_pad = get_sequence(left_key, species)
        right_pad = get_sequence(right_key, species)
        h_seq = row["haplotype_sequence"]
        ori = row["orientation"]
        if ori == "reverse":
            hs = left_pad + reverse_complement(h_seq) + right_pad
            return reverse_complement(hs)
        else:
            return left_pad + h_seq + right_pad

    if pad_size > 0:
        mapped_haplotypes["padded_haplotype_sequence"] = (
            mapped_haplotypes.apply(get_padded_haplotype_sequence,
                                    args=(pad_size, species), axis=1))
    else:
        mapped_haplotypes["padded_haplotype_sequence"] = (
            mapped_haplotypes["haplotype_sequence"])
    mipster = mipster.merge(mapped_haplotypes[["haplotype_ID",
                                               "padded_haplotype_sequence",
                                               "mapped_copy_number"]])
    mipster["raw_haplotype_sequence"] = mipster["haplotype_sequence"]
    mipster["haplotype_sequence"] = mipster["padded_haplotype_sequence"]
    mipster["raw_sequence_quality"] = mipster["sequence_quality"]
    mipster["sequence_quality"] = (
        pad_size * "!" + mipster["sequence_quality"] + pad_size * "!")
    mipster["adjusted_barcode_count"] = (
        mipster["barcode_count"] / mipster["mapped_copy_number"]).astype(int)
    mipster["raw_barcode_count"] = mipster["barcode_count"]
    mipster["barcode_count"] = mipster["adjusted_barcode_count"]
    mipster_dict = mipster.loc[mipster["barcode_count"]
                               >= min_bc_count].groupby(
        "sample_name").apply(lambda x: pd.DataFrame.to_dict(
                             x, orient="index")).to_dict()
    if save:
        mipster.to_csv(mipster_file + ".padded")
    p = Pool(pro)
    for sample in mipster_dict:
        fastq_file = os.path.join(fastq_dir, sample + ".fq.gz")
        sample_mipster = mipster_dict[sample]
        p.apply_async(generate_processed_fastqs_worker, (fastq_file,
                                                         sample_mipster))
    p.close()
    p.join()
    return


def generate_unprocessed_fastqs(fastq_dir, mipster_file, min_bc_count=1,
                                pro=8):
    """
    Generate fastq files for each sample. These files will have stitched and
    barcode corrected reads.
    """
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    mipster = pd.read_table(mipster_file, usecols=["s_Sample", 'h_popUID',
                                                   "h_seq", 'c_qual',
                                                   'c_barcodeCnt'])
    mipster = mipster.rename(columns={"s_Sample": "sample_name",
                                      "h_popUID": "haplotype_ID",
                                      "h_seq": "haplotype_sequence",
                                      "c_qual": "sequence_quality",
                                      "c_barcodeCnt": "barcode_count"})
    mipster = mipster.loc[mipster["barcode_count"] >= min_bc_count].groupby(
        "sample_name"
    ).apply(lambda x: pd.DataFrame.to_dict(x, orient="index")).to_dict()
    p = Pool(pro)
    for sample in mipster:
        fastq_file = os.path.join(fastq_dir, sample + ".fq.gz")
        sample_mipster = mipster[sample]
        p.apply_async(generate_processed_fastqs_worker, (fastq_file,
                                                         sample_mipster))
    p.close()
    p.join()
    return


def convert_to_int(n):
    """Convert values to integers.

    This is to be used when a pandas dataframe converts integers to floats due
    to the presence of NA values and integer values are
    preferred over floats, e.g. string conversion/comparison.
    """
    try:
        return int(n)
    except ValueError:
        return np.nan


def get_ternary_genotype(gen):
    """Convert a 0/0, 0/1, 1/1 type genotype string to 0, 1, 2."""
    try:
        g = sum(map(int, gen.split(":")[0].split("/")))
    except ValueError:
        g = np.nan
    return g


def variation_to_geno(settings, var_file, output_prefix):
    """Create PLINK files from variation table file."""
    wdir = settings["workingDir"]
    case = {}
    with open(os.path.join(wdir, "case_file")) as infile:
        for line in infile:
            newline = line.strip().split("\t")
            case[newline[0]] = newline[1]
    with open(os.path.join(wdir, var_file)) as infile:
        linenum = 0
        map_snp_ids = []
        geno_snp_ids = []
        genes = {}
        ordered_genes = []
        all_geno_bases = []
        all_geno_numbers = []
        for line in infile:
            newline = line.strip().split("\t")
            if linenum == 0:
                linenum += 1
                header = newline
                sample_start_index = header.index("FORMAT") + 1
                sample_ids = header[sample_start_index:]
                sample_names = ["-".join(s.split("-")[:-2])
                                for s in sample_ids]
                samples_used = []
                ped_sample_info = []
                geno_sample_info = []
                for i in range(len(sample_ids)):
                    s = sample_ids[i]
                    sam_name = sample_names[i]
                    try:
                        affected = case[s]
                        if affected == "case":
                            affected = "2"
                        elif affected == "control":
                            affected = "1"
                        ped_sample_info.append(["0", sam_name, "0", "0", "0",
                                                affected])
                        geno_sample_info.append(["0", sam_name, sam_name,
                                                 affected])
                        samples_used.append(s)
                    except KeyError:
                        continue
                used_sample_mask = np.array([s in samples_used
                                             for s in sample_ids])
            else:
                chrom = newline[0]
                pos = newline[1]
                rsid = newline[2]
                ref = newline[3]
                alt = newline[4]
                aa_change = newline[5]
                gene_name = newline[7]
                if rsid == ".":
                    rsid = chrom[3:] + "-" + pos + "-" + ref + "-" + alt
                map_snp_ids.append([chrom, rsid, "0", pos])
                geno_snp_ids.append([chrom, rsid, pos, ref,
                                     alt, gene_name, aa_change])
                try:
                    genes[gene_name].append(rsid)
                except KeyError:
                    genes[gene_name] = [rsid]
                    ordered_genes.append(gene_name)
                genotypes = np.array(newline[sample_start_index:])[
                    used_sample_mask]
                geno_bases_1 = []
                geno_bases_2 = []
                geno_numbers = []
                for g in genotypes:
                    temp_g = g.split(":")[0]
                    if temp_g == "0/0":
                        geno_bases_1.append(ref)
                        geno_bases_2.append(ref)
                        geno_numbers.append("0")
                    elif temp_g == "0/1":
                        geno_bases_1.append(ref)
                        geno_bases_2.append(alt)
                        geno_numbers.append("1")
                    elif temp_g == "1/1":
                        geno_bases_1.append(alt)
                        geno_bases_2.append(alt)
                        geno_numbers.append("2")
                    else:
                        geno_bases_1.append("0")
                        geno_bases_2.append("0")
                        geno_numbers.append(".")
                all_geno_bases.extend([geno_bases_1,
                                       geno_bases_2])
                all_geno_numbers.append(geno_numbers)
        all_geno_bases = list(zip(*all_geno_bases))
        all_geno_numbers = list(zip(*all_geno_numbers))
        for i in range(len(ped_sample_info)):
            ped_sample_info[i].extend(all_geno_bases[i])
        for i in range(len(geno_sample_info)):
            geno_sample_info[i].extend(all_geno_numbers[i])
    write_list(ped_sample_info,
               os.path.join(wdir, output_prefix + ".ped"))
    write_list(map_snp_ids,
               os.path.join(wdir, output_prefix + ".map"))
    header = ["FAMILY_ID", "INDIVIDUAL_ID", "SAMPLE_ID", "AFFECTION"]
    header.extend([s[1] for s in geno_snp_ids])
    geno_sample_info = [header] + geno_sample_info
    write_list(geno_sample_info, os.path.join(wdir, output_prefix + ".geno"))
    write_list([["**", o] + genes[o] for o in ordered_genes],
               os.path.join(wdir, output_prefix + ".hlist"))
    return


def absence_presence(col, min_val=1):
    """
    Given a numerical dataframe column, convert to binary values for a minimum
    threshold. This should be used by pandas transform or apply.
    """
    return pd.Series([0 if (c < min_val or np.isnan(c))
                      else 1 for c in col.tolist()])


def plot_performance(barcode_counts,
                     tick_label_size=8,
                     cbar_label_size=5,
                     dpi=300,
                     barcode_threshold=1,
                     absent_color="black",
                     present_color="green",
                     save=False,
                     wdir=None,
                     ytick_freq=None,
                     xtick_freq=None,
                     xtick_rotation=90,
                     tick_genes=False,
                     gene_name_index=None):
    """Plot presence/absence plot for a mip run."""
    if xtick_freq is None:
        xtick_freq = barcode_counts.shape[1]//30
        if xtick_freq == 0:
            xtick_freq = 1
    if ytick_freq is None:
        ytick_freq = barcode_counts.shape[0]//30
        if ytick_freq == 0:
            ytick_freq = 1
    fig, ax = plt.subplots()
    cmap = colors.ListedColormap([absent_color, present_color])
    boundaries = [-0.5, 0.5, 1.5]
    norm = colors.BoundaryNorm(boundaries, cmap.N)
    heat = ax.pcolormesh(
        barcode_counts.applymap(
            lambda a: np.nan if np.isnan(a)
            else 0 if a < barcode_threshold
            else 1
        ), cmap=cmap, norm=norm)
    sample_ids = list(barcode_counts.index)
    sample_locs = np.arange(1, len(sample_ids) + 1, ytick_freq) - 0.5
    ylabs = sample_ids[::ytick_freq]
    plt.yticks(sample_locs, ylabs)
    if tick_genes:
        bc_cols = barcode_counts.columns.tolist()
        bc_cols = [c[gene_name_index] for c in bc_cols]
        xlabs = bc_cols[::xtick_freq]
        gene_locs = np.arange(1, len(bc_cols) + 1, xtick_freq) - 0.5
        plt.xticks(gene_locs, xlabs, rotation=xtick_rotation, ha="right")
    for ticklabel in ax.get_xticklabels():
        ticklabel.set_fontsize(tick_label_size)
    for ticklabel in ax.get_yticklabels():
        ticklabel.set_fontsize(tick_label_size)
    ax.set_ylabel("Samples")
    ax.set_xlabel("Probes")
    fig.suptitle("Performance",
                 verticalalignment="bottom")
    fig.tight_layout()
    cbar = fig.colorbar(heat, ticks=[0, 1], shrink=0.2)
    cbar.ax.tick_params(labelsize=cbar_label_size)
    cbar.ax.set_yticklabels(["Absent", "Present"])
    fig.set_dpi(dpi)
    fig.tight_layout()
    if save:
        fig.savefig(os.path.join(wdir, "performance.png"), dpi=dpi,
                    bbox_inches='tight')
        plt.close("all")
    else:
        return fig, ax
    return


def plot_coverage(barcode_counts,
                  tick_label_size=8,
                  cbar_label_size=5,
                  dpi=300,
                  log=None,
                  log_constant=1,
                  linthresh=0.0001,
                  save=False,
                  wdir=None,
                  ytick_freq=None,
                  xtick_freq=None,
                  xtick_rotation=90,
                  tick_genes=False,
                  gene_name_index=None,
                  figure_title="Coverage",
                  title_fontdict=None,
                  ylabel="Samples",
                  xlabel="Probes",
                  fig_size=None,
                  cbar_title=None):
    """Plot UMI coverage per MIP per sample for a mip run."""
    if xtick_freq is None:
        xtick_freq = barcode_counts.shape[1]//30
        if xtick_freq == 0:
            xtick_freq = 1
    if ytick_freq is None:
        ytick_freq = barcode_counts.shape[0]//30
        if ytick_freq == 0:
            ytick_freq = 1
    fig, ax = plt.subplots()
    if log is None:
        heat = ax.pcolormesh(barcode_counts)
        if cbar_title is None:
            cbar_title = ""
    elif log == 2:
        if cbar_title is None:
            cbar_title = "log2"
        heat = ax.pcolormesh(np.log2(barcode_counts + log_constant))
    elif log == 10:
        if cbar_title is None:
            cbar_title = "log10"
        heat = ax.pcolormesh(np.log10(barcode_counts + log_constant))
    elif log == "ln":
        if cbar_title is None:
            cbar_title = "log"
        heat = ax.pcolormesh(np.log(barcode_counts + log_constant))
    elif log == "symlog":
        if cbar_title is None:
            cbar_title = ""
        heat = ax.pcolormesh(barcode_counts,
                             norm=colors.SymLogNorm(linthresh=linthresh))
    else:
        print("log can only be None, 2, 10, 'log', {} provided.".format(log))
    sample_ids = list(barcode_counts.index)
    sample_locs = np.arange(1, len(sample_ids) + 1, ytick_freq) - 0.5
    ylabs = sample_ids[::ytick_freq]
    plt.yticks(sample_locs, ylabs)
    if tick_genes:
        bc_cols = barcode_counts.columns.tolist()
        bc_cols = [c[gene_name_index] for c in bc_cols]
        xlabs = bc_cols[::xtick_freq]
        gene_locs = np.arange(1, len(bc_cols) + 1, xtick_freq) - 0.5
        plt.xticks(gene_locs, xlabs,
                   rotation=xtick_rotation,
                   ha="right")
    for ticklabel in ax.get_xticklabels():
        ticklabel.set_fontsize(tick_label_size)
    for ticklabel in ax.get_yticklabels():
        ticklabel.set_fontsize(tick_label_size)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(figure_title, fontdict=title_fontdict)
    cbar = fig.colorbar(heat, shrink=0.5)
    cbar.ax.tick_params(labelsize=cbar_label_size)
    cbar.ax.set_ylabel(cbar_title,
                       fontsize=cbar_label_size,
                       rotation=90)
    fig.set_dpi(dpi)
    if fig_size is not None:
        fig.set_size_inches(fig_size)
    fig.tight_layout()
    if save:
        fig.savefig(os.path.join(wdir, "coverage.png"),
                    dpi=dpi,
                    bbox_inches='tight')
        plt.close("all")
    else:
        return fig, ax
    return


def get_mutation_counts(col):
    return col.apply(lambda gen: int(gen.split(":")[1]))


def get_totals(col):
    return col.apply(lambda gen: int(gen.split(":")[2]))


def get_coverage(row, sorted_counts):
    chrom = row["Chrom"]
    start = row["Start"]
    end = row["End"]
    sid = row["Sample ID"]
    idx = pd.IndexSlice
    return sorted_counts.loc[sid, idx[:, :, :, :, chrom, :start, end:]].sum()


def add_known(group, used_targets):
    group = group.merge(used_targets, how="outer")
    group["Sample ID"].fillna(method="ffill", inplace=True)
    group["Sample ID"].fillna(method="bfill", inplace=True)
    group["Chrom"].fillna(group["CHROM"], inplace=True)
    group["Start"].fillna(group["POS"], inplace=True)
    group["End"].fillna(group["POS"], inplace=True)
    group["CHROM"].fillna(group["Chrom"], inplace=True)
    group["POS"].fillna(group["Start"], inplace=True)
    group["Barcode Count"].fillna(0, inplace=True)
    return group


def find_ref_total(group):
    nr = group.loc[~group["ExonicFunc"].isin(
        ["synonymous_variant", "."]), "Barcode Count"].sum()
    cov = group["POS Coverage"].max()
    return cov - nr


def get_genotype(row, cutoff):
    if row["Coverage"] > 0:
        if row["Filtered Barcode Fraction"] >= cutoff:
            if row["Filtered Barcode Fraction"] > (1 - cutoff):
                return "MUT"
            else:
                return "MIX"
        else:
            return "WT"
    else:
        return np.nan


def get_aminotype(row, cutoff):
    if row["Coverage"] > 0:
        if row["Mutation Fraction"] >= cutoff:
            if row["Mutation Fraction"] > (1 - cutoff):
                return "MUT"
            else:
                return "MIX"
        else:
            return "WT"
    else:
        return np.nan


def rename_noncoding(row):
    return "-".join([row["Gene"],
                     row["VKEY"]])


def call_microsats(settings, sim=None, freq_cutoff=0.005, min_bc_cutoff=0,
                   use_filtered_mips=True, ref_genome="Pf3d7"):
    wdir = settings["workingDir"]
    with open(os.path.join(wdir, settings["perSampleResults"])) as infile:
        sample_results = json.load(infile)
    with open(os.path.join(wdir, settings["haplotypeDictionary"])) as infile:
        hap_dict = json.load(infile)
    if sim is None:
        sim = pd.read_csv("resources/pf_MS/simulation.tsv", sep="\t")
    ref_sim = sim.loc[sim["genome"] == ref_genome]
    strain_freqs = {}
    for sample_name in sample_results:
        sam_res = sample_results[sample_name]
        sam_freq = {}
        for g in sam_res:
            for m in sam_res[g]:
                if use_filtered_mips and (ref_sim.loc[
                        ref_sim["MIP"] == m,
                        "Clean MS MIP"].values[0] is False):
                    continue
                for c in sam_res[g][m]:
                    total_bcs = float(sam_res[g][m][c]["cumulative_data"][
                        "barcode_count"])
                    if total_bcs >= min_bc_cutoff:
                        filtered_data = sam_res[g][m][c]["filtered_data"]
                        ms_types = {}
                        for hd in filtered_data:
                            bcc = hd["barcode_count"]
                            h = hd["haplotype_ID"]
                            h_len = len(hap_dict[m][h]["sequence"])
                            ms_len = int(h_len - ref_sim.loc[ref_sim["MIP"]
                                         == m, "MS size adjustment"])
                            try:
                                ms_types[ms_len] += bcc
                            except KeyError:
                                ms_types[ms_len] = bcc
                        for ml in list(ms_types.keys()):
                            if (ms_types[ml]/total_bcs) < freq_cutoff:
                                ms_types.pop(ml)
                        try:
                            sam_freq[g][m][c] = ms_types
                        except KeyError:
                            try:
                                sam_freq[g][m] = {c: ms_types}
                            except KeyError:
                                sam_freq[g] = {m: {c: ms_types}}
        strain_freqs[sample_name] = sam_freq
    ms_calls = []
    for s in strain_freqs:
        for g in strain_freqs[s]:
            for m in strain_freqs[s][g]:
                for c in strain_freqs[s][g][m]:
                    for l in strain_freqs[s][g][m][c]:
                        ms_calls.append([s, g, m, c, l, g + "-" + str(int(l)),
                                        strain_freqs[s][g][m][c][l]])
    ms_call_df = pd.DataFrame(
        ms_calls, columns=["Sample ID", "region", "MIP", "Copy", "length",
                           "haplotype name", "count"]).drop("Copy", axis=1)
    merged_calls = pd.DataFrame(ms_call_df.groupby(["Sample ID", "region",
                                                    "haplotype name", "length",
                                                    ])["count"].sum())
    merged_calls.reset_index(inplace=True)
    merged_calls["frequency"] = merged_calls.groupby(
        ["Sample ID", "region"])["count"].transform(lambda a: a/a.sum())
    merged_calls.rename(columns={"length": "MS Length"}, inplace=True)
    merged_calls = merged_calls.merge(sim.groupby(
        ["region", "MS Length", "Unique Strain"], as_index=False).first()[
        ["region", "MS Length", "Unique Strain"]], how="left")
    return {"ms_calls": merged_calls,
            "strain_freqs": strain_freqs}


def get_copy_counts(count_table, average_copy_count=2,
                    norm_percentiles=[0.4, 0.6]):
    """
    Given a table of barcode counts with samples on rows
    and probes on columns, transform the table to return
    estimated copy count.

    Parameters
    ----------
    count_table : numpy array/pandas dataframe
        Table of barcode counts with samples on rows
        and probes on columns
    average_copy_count : float, 2
        Most common copy number in population. 2 for humans
        This number is used to assign the copy number for
        median/average normalized barcode count.
    norm_percentiles : length 2 list of floats between 0 and 1, [0.4, 0.6]
        Percentiles used for calculating average. [0.5, 0.5] would be median.
    """
    # Normalize samples (across columns)
    s_norm = count_table.transform(
        lambda a: a/a.sum(), axis=1)
    # Normalize across samples. This achieves estimating
    # the copy number, assuming the average normalized
    # barcode value (at specified percentile) is the value
    # provided by averageCopyCount setting. This should
    # default to something like median barcode count
    # corresponds to copy number 2.
    p_norm = s_norm.transform(
        lambda a: average_copy_count * a/(a.quantile(norm_percentiles).mean()))
    return p_norm


def get_copy_average(r, ac):
    try:
        return ac.loc[r["Sample ID"], (r["Gene"], r["Copy"])]
    except KeyError:
        return np.nan


def normalize_copies(a):
    if a.isnull().all():
        a = a.fillna(1)
        return a/a.sum()
    else:
        return a.fillna(0)


def repool(wdir,
           data_summary,
           high_barcode_threshold,
           target_coverage_count=None,
           target_coverage_fraction=0.95,
           target_coverage_key="targets_with_10_barcodes",
           barcode_coverage_threshold=10,
           barcode_count_threshold=100,
           low_coverage_action="Repool",
           assesment_key="targets_with_1_barcodes",
           good_coverage_quantile=0.25,
           output_file="repool.csv"):
    """
    Analyze run statistics and determine repooling/recapturing
    strategy for following runs.

    Parameters
    ----------
    wdir : str
        Path to working directory, used only for saving the results.
    data_summary : Pandas DataFrame
        Dataframe containing all count information per sample per target.
    high_barcode_threshold: int/ other number
        Targeted barcode number to determine how much more of a
        sample should be repooled. Should be set to the number where
        majority of samples show good target coverage.
    barcode_coverage_threshold : int / other number, 10
        Average reads per barcode per sample to consider the sample
        saturated and remove from pooling. If sample is not deemed
        complete it will be set to be recaptured.
    barcode_count_threshold : int / other number
        Minimum number of barcodes per sample to determine if a sample
        has very low coverage. Those samples' status will be set
        to the action (recapture or repool) defined by
        low_coverage_action parameter.
    target_coverage_count : int / other number / None, None
        Minimum number of targets (MIPs) that are sequenced/covered
        to the given criteria to consider a sample complete. Defaults
        to None, in which case target_coverage_fraction * total number
        of possible targets will be used.
    target_coverage_fraction : float, 0.95
        See target_coverage_count.
    target_coverage_key : str, "targets_with_10_barcodes"
        Dataframe column name to use for assessing target coverage.
        By default a target that is covered with >10 barcodes will
        be considered covered.
    assesment_key : str, "targets_with_1_barcodes"
        Dataframe key to use for determining uneven coverage across targets
        which happens when barcode number per sample is high but number of
        targets covered is low. By default any target with sequence is
        considered covered.
    good_coverage_quantile : float, 0.25
        Quantile of barcodes for "completed samples". This is used to determine
        if a sample has good enough barcode numbers, then test if it has enough
        targets covered, or the majority of barcodes cover only a small number
        of targets (uneven coverage).
    output_file: str, repool.csv
    """
    if target_coverage_count is None:
        target_coverage_count = (data_summary[target_coverage_key].max()
                                 * target_coverage_fraction)
    # make a copy of data_summary so the original df stays the same
    data_summary = copy.deepcopy(data_summary)
    try:
        data_summary["total_barcode_count"]
    except KeyError:
        data_summary["total_barcode_count"] = data_summary["Barcode Count"]
        data_summary["total_read_count"] = data_summary["Read Count"]
    # mark samples that reached the desired outcome
    data_summary.loc[
        data_summary[target_coverage_key] >= target_coverage_count,
        "Status"
    ] = "Complete"
    # mark samples with low coverage
    data_summary.loc[
        (data_summary["Status"].isnull())
        & (data_summary["total_barcode_count"] < barcode_count_threshold),
        "Status"] = low_coverage_action
    # mark samples with too high barcode coverage
    # these samples will have been sequenced to a high depth but
    # low barcode numbers, so sequencing these more would not make sense.
    # They will be re-captured if more data is needed.
    try:
        data_summary["Barcode Coverage"]
    except KeyError:
        data_summary["Barcode Coverage"] = (
            data_summary["total_read_count"]
            / data_summary["total_barcode_count"]).fillna(0)
    data_summary.loc[
        (data_summary["Status"].isnull())
        & (data_summary["Barcode Coverage"] >= barcode_coverage_threshold),
        "Status"] = "Recapture"
    # Zero barcode coverage is presumably due to poor sequencing
    # So low coverage action should be taken.
    data_summary.loc[
        (data_summary["Status"].isnull())
        & (data_summary["Barcode Coverage"] == 0),
        "Status"
    ] = low_coverage_action
    # All remaining samples will be repooled
    data_summary.loc[
        (data_summary["Status"].isnull()),
        "Status"
    ] = "Repool"
    data_summary["Library to Completion"] = (
        (high_barcode_threshold - data_summary["total_barcode_count"])
        / data_summary["total_barcode_count"])
    # replace inf values with max
    lc_max = data_summary.loc[
        data_summary["Library to Completion"] < np.inf,
        "Library to Completion"].max()
    data_summary.loc[
        data_summary["Library to Completion"] == np.inf,
        "Library to Completion"] = lc_max
    # Determine samples with good barcode counts but poor target coverage
    # These should be investigated to decide what is the reason behind it
    # and how to proceed.
    ##########################################
    # Determine the average barcode count per target covered
    # for all samples where there is targets covered
    data_summary.loc[data_summary[target_coverage_key] > 0,
                     "Barcodes Per Target Covered"] = (
        data_summary.loc[data_summary[target_coverage_key] > 0,
                         "total_barcode_count"]
        / data_summary.loc[data_summary[target_coverage_key] > 0,
                           target_coverage_key]
    )
    # Get the lower quartile of barcodes per target for good data
    # This number will be used to determine if poor coverage samples
    # have high enough barcode coverage despite having poor target coverage.
    good_coverage_threshold = data_summary.loc[
        data_summary["Status"] == "Complete",
        "Barcodes Per Target Covered"].quantile(good_coverage_quantile)
    # Determine samples where barcode coverage is high but target coverage
    # is low
    data_summary.loc[
        (data_summary["Barcodes Per Target Covered"] > good_coverage_threshold)
        & (data_summary[assesment_key] < target_coverage_count),
        "Uneven Coverage"] = True
    data_summary.loc[data_summary["Uneven Coverage"].isnull(),
                     "Uneven Coverage"] = False
    try:
        data_summary.to_csv(os.path.join(wdir, output_file), index=False)
    except TypeError:
        # in an older version of this function, settings dict
        # was passed instead of wdir, for backwards compatibility
        # we'll catch that error and use wdir from the settings dict
        data_summary.to_csv(wdir["workingDir"] + output_file, index=False)
    print(("Out of %d samples %d are completed,"
           " %d will be recaptured and %d repooled" % (
               data_summary.shape[0],
               data_summary.loc[data_summary["Status"] == "Complete"].shape[0],
               data_summary.loc[data_summary["Status"]
                                == "Recapture"].shape[0],
               data_summary.loc[data_summary["Status"] == "Repool"].shape[0])))
    print(("%d samples showed uneven coverage, %d complete,"
           " %d to be recaptured, %d repooled" % (
                data_summary.loc[data_summary["Uneven Coverage"]].shape[0],
                data_summary.loc[data_summary["Uneven Coverage"] &
                                 (data_summary["Status"]
                                  == "Complete")].shape[0],
                data_summary.loc[data_summary["Uneven Coverage"]
                                 & (data_summary["Status"]
                                 == "Recapture")].shape[0],
                data_summary.loc[data_summary["Uneven Coverage"]
                                 & (data_summary["Status"]
                                 == "Repool")].shape[0])))
    return


def aa_to_coordinate(gene, species, aa_start, aa_end=None, alias=False):
    """
    Given a gene name and its amino acid location,
    return the genomic coordinates of the aa.
    This will work with most Plasmodium genes but will
    be problematic when genes have multiple isoforms,
    such as most human genes.
    """
    if alias:
        with open(get_file_locations()[species]["alias"]) as infile:
            alias_dic = json.load(infile)
        try:
            gene = alias_dic[gene]
        except KeyError:
            pass
    cds = get_cds(gene, species)
    if len(cds) == 0:
        return [np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan]
    ori = cds["orientation"]
    coord = cds["coordinates"]
    chrom = cds["chrom"]
    if aa_end is None:
        aa_end = aa_start
    if ori == "+":
        c_end = aa_end * 3 - 1
        c_start = aa_start * 3 - 3
    else:
        c_start = aa_end * 3 - 1
        c_end = aa_start * 3 - 3
    cds_start = coord[c_start]
    cds_end = coord[c_end]
    if ori == "+":
        codon = get_sequence(
            create_region(
                chrom, cds_start, cds_end
            ), species
        )
    else:
        codon = reverse_complement(
            get_sequence(
                create_region(
                    chrom, cds_start, cds_end
                ), species
            )
        )
    aa = translate(codon)
    return [chrom, cds_start, cds_end, ori, codon, aa]


def merge_snps(settings):
    """
    When more than one SNP affects the same codon of a gene,
    merge the two SNPs and create a merged protein change
    annotation.
    """
    wdir = settings["workingDir"]
    species = settings["species"]
    # load haplotype dictionary, save a backup.
    unique_haplotype_file = os.path.join(wdir, settings["haplotypeDictionary"])
    with open(unique_haplotype_file) as infile:
        haplotypes = json.load(infile)
    # create output information list to report all changes
    outlist = [["HaplotypeID", "Copy", "New AA", "Reference AA",
                "ReplacedCDNAchanges", "New Codon", "Reference Codon",
                "ReplacedAaChanges"]]
    # go through each annotated haplotype and merge SNPs
    for m in haplotypes:
        for h in haplotypes[m]:
            if haplotypes[m][h]["mapped"]:
                for cp in haplotypes[m][h]["mapped_copies"]:
                    # get sequence of the haplotype
                    hap_seq = haplotypes[m][h]["sequence"]
                    # get SNPs present in the haplotype
                    diffs = haplotypes[m][h]["mapped_copies"][cp][
                        "differences"]
                    aa_changes = {}
                    multi_indels = []
                    for i in range(len(diffs)):
                        d = diffs[i]
                        # get protein change information for the SNP
                        aa = d["annotation"]["AA Change"]
                        try:
                            aa_pos = int(d["annotation"]["AA Change Position"])
                            try:
                                # add the aa change position to the changes
                                # dict.
                                aa_changes[aa_pos].append(i)
                            except KeyError:
                                aa_changes[aa_pos] = [i]
                        except ValueError:
                            continue
                    # after going through all diffs, look for mutliple diffs
                    # affecting single aminoacid
                    all_merges = []
                    for c in aa_changes:
                        if (len(aa_changes[c]) > 1) and (
                                 c not in multi_indels):
                            # break out of loop if indels found
                            mindel = False
                            indexes = aa_changes[c]
                            # keep positions relative to cDNA in a list
                            c_positions = []
                            # keep positions relative to haplotype in a list
                            h_indexes = []
                            # keep genomic positions of changes in a list
                            g_positions = []
                            # keep the difference between the cDNA and genomic
                            # positions in a list. This will be used to
                            # determine the gene's orientation on the genome.
                            c_offsets = []
                            changes_to_cdna = []
                            changes_to_aa = []
                            for i in indexes:
                                d = diffs[i]
                                # for each diff get the annotation
                                # e.g. 'HBB:HBB:exon1:c.G673A:p.V225I'
                                ano = d["annotation"]
                                if ano["SNV"] != "SNV":
                                    mindel = True
                                    multi_indels.append(c)
                                    break
                                aa = ano["AA Change"]
                                changes_to_aa.append(aa)
                                # get the aa of reference genome (V)
                                aa_ref = aa[0]
                                # get cdna change, e.g. G673A
                                cdna = ano["CDS Change"]
                                changes_to_cdna.append(cdna)
                                # get the mutant base (A)
                                cdna_change = cdna[-1]
                                # compare the sequence on the cDNA with
                                # to the sequence of the haplotype,
                                # to determine the MIP/haplotype's orientation
                                # relative to the cDNA
                                ori = cdna_change == d["hap_base"]
                                cdna_pos = int(ano["CDS Position"])
                                # get genomic position of the change
                                diff_start = int(d["annotation"]["Start"])
                                # get the difference between the genomic and
                                # cDNA position of the change, to be used
                                # in determining the gene's orientation
                                c_positions.append(cdna_pos)
                                h_indexes.extend(d["hap_index"])
                                g_positions.append(diff_start)
                                c_offsets.append(diff_start - cdna_pos)
                            if mindel:
                                break
                            c_positions = sorted(c_positions)
                            h_indexes = sorted(set(h_indexes))
                            g_positions = sorted(g_positions)
                            # if the offset between the cDNA and genomic
                            # positions of the changes are always the same,
                            # the gene is on the plus strand, else it is
                            # reverse.
                            if len(set(c_offsets)) > 1:
                                gene_ori = False
                            else:
                                gene_ori = True
                            # get the position of the first base of the codon
                            codon_offset = (c_positions[0] % 3) - 1
                            codon_pos = c_positions[0] - codon_offset
                            # get the codon's sequence from the haplotype
                            if ori:
                                # if the haplotype is in the same orientation
                                # as the cDNA
                                h_start_index = h_indexes[0] - codon_offset
                                h_end_index = h_start_index + 3
                                hap_codon = hap_seq[h_start_index:h_end_index]
                                codon = hap_codon
                            else:
                                # if the haplotype is in the opposite
                                # orientation as the cDNA
                                h_end_index = h_indexes[-1] + codon_offset + 1
                                h_start_index = h_end_index - 3
                                hap_codon = hap_seq[h_start_index:h_end_index]
                                codon = reverse_complement(hap_codon)
                            # get genomic position and sequence of the codon
                            if gene_ori:
                                g_start = g_positions[0] - codon_offset
                                g_end = g_start + 2
                            else:
                                g_end = g_positions[-1] + codon_offset
                                g_start = g_end - 2
                            # extract the reference codon sequence
                            Ref = get_sequence(create_region(
                                d["chrom"], g_start, g_end), species)
                            if gene_ori:
                                g_codon = Ref
                                Alt = codon
                            else:
                                g_codon = reverse_complement(Ref)
                                Alt = reverse_complement(codon)
                            # calculate merged codon's amino acid
                            merged_aa = translate(codon)
                            # recreate the annotation string for the merge
                            protein_change = aa_ref + str(c) + merged_aa
                            coding_change = g_codon + str(codon_pos) + codon
                            # determine if the merged change is synonymous
                            if aa_ref == merged_aa:
                                ExonicFunc = "synonymous_variant"
                            else:
                                ExonicFunc = "nonsynonymous_variant"
                            merged_dict = {'annotation': {
                                'AA Change': protein_change,
                                "AA Change Position": str(c),
                                "CDS Change": coding_change,
                                "CDS Position": str(codon_pos),
                                'Alt': Alt,
                                'Chr': d["chrom"],
                                'End': g_end,
                                'ExonicFunc': ExonicFunc,
                                'GeneID': d["annotation"]['GeneID'],
                                'Ref': Ref,
                                'Start': g_start,
                                "SNV": "MNV"},
                                           'begin': g_start,
                                           'chrom': d["chrom"],
                                           'end': g_end,
                                           'hap_base': hap_codon,
                                           'hap_index': [h_start_index,
                                                         h_end_index - 1],
                                           'ref_base': Ref,
                                           'type': 'snp',
                                           'vcf_normalized': ":".join(
                                                [d["chrom"], str(g_start), ".",
                                                 Ref, Alt]),
                                           'vcf_raw': ":".join(
                                              [d["chrom"], str(g_start), ".",
                                               Ref, Alt]
                                            ),
                                           "gene_ori": gene_ori,
                                           "ori": ori}
                            all_merges.append(merged_dict)
                            outlist.append([h, c, merged_aa, aa_ref,
                                           ",".join(changes_to_cdna),
                                            codon, g_codon,
                                            ",".join(changes_to_aa)])
                    # Remove SNPs that were merged, add the merged SNP
                    for c in aa_changes:
                        if ((len(aa_changes[c]) > 1)
                                and (c not in multi_indels)):
                            indexes = aa_changes[c]
                            for i in indexes:
                                diffs[i] = "remove"
                    diffs.extend(all_merges)
                    diffs = [d for d in diffs if d != "remove"]
                    haplotypes[m][h]["mapped_copies"][cp][
                        "differences"] = diffs
    with open(unique_haplotype_file, "w") as outfile:
        json.dump(haplotypes, outfile, indent=1)
    # save the report
    write_list(outlist, os.path.join(wdir, "merge_snps_output.txt"))
    return outlist


def load_processed_data(settings):
    """
    Load the data after initial processing.
    Data included will be barcode counts, meta data
    and data summary.
    """
    wdir = settings["workingDir"]
    # load the barcode count data and save a transposed
    # version of it, only to load the transposed version
    # back again and re-transpose. The reason for this is
    # that read_csv method of pandas dataframe does not
    # interpret multi index column names as integers but
    # it does so for multi index index names. begin and
    # end multi index columns should have names as int.
    bc = pd.read_csv(os.path.join(wdir, "barcode_counts.csv"),
                     header=[0, 1, 2, 3, 4, 5, 6], index_col=0)
    bc.T.to_csv(os.path.join(wdir, "barcode_counts.T.csv"))
    bc = pd.read_csv(os.path.join(wdir, "barcode_counts.T.csv"),
                     index_col=[0, 1, 2, 3, 4, 5, 6]).T
    data_summary = pd.read_csv(os.path.join(wdir, "data_summary.csv"),
                               index_col=None)
    merged_meta = pd.read_csv(os.path.join(wdir, "meta_data.csv"),
                              index_col=None)
    bc.index.name = "Sample ID"
    return {"Barcode Counts": bc,
            "Data Summary": data_summary,
            "Meta Data": merged_meta}


def vcf_to_df(vcf_file):
    """
    Convert a possibly compressed (.gz) vcf file to a Pandas DataFrame.

    Parameters:
    ----------
    vcf_file : Path to the vcf file. The file must have the 8 columns
        specified in vcf specifications: CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO
        It can be compressed. Individual genotypes are not used, so they
        can be present or absent. Each indel must have their own line
        in the file, i.e. bcftools norm -m -indels

    Returns : Pandas dataframe with each row corresponding to a position
        and not an individual variant.
    """
    # Check if the file is compressed
    try:
        op = gzip.open(vcf_file).readline()
        op = gzip.open(vcf_file)
    except IOError:
        op = open(vcf_file).readline()
        op = open(vcf_file)
    # create a list of variants that behave unexpectedly
    # this will be used at the end of the function to make sure
    # everything went as expected.
    problem_alts = []
    # keep each variant information in a list to be corverted to Df
    outfile_list = []
    # keep INFO field headers in a list
    info_cols = []
    with op as infile:
        for line in infile:
            try:
                line = line.decode("utf-8")
            except AttributeError:
                pass
            if line.startswith("#"):
                # extract INFO field headers
                if line.startswith("##INFO="):
                    info_cols.append(line.split(",")[0].split("=")[2])
            else:
                newline = line.strip().split("\t")
                chrom = newline[0]
                pos = int(newline[1])
                var_id = newline[2]
                ref = newline[3]
                alt = newline[4]
                qual = newline[5]
                filt = newline[6]
                info_raw = newline[7].split(";")
                info_dict = {}
                for ir in info_raw:
                    split_info = ir.split("=")
                    # the variant info is coded as "field=info"
                    # when field has a vale, and just "field"
                    # when field is a flag.
                    # e.g. "AC=4" shows allele count is 4
                    # e.g. "STR" shows the variant is short tandem repeat
                    try:
                        info_dict[split_info[0]] = split_info[1].split(",")
                    # if field is a flag
                    except IndexError:
                        info_dict[split_info[0]] = [True]
                alt_bases = alt.split(",")
                alt_len = len(alt_bases[0])
                # since all the indels have their own line on this vcf file
                # all alt bases must have the same length (1 for SNPs and
                # indel size for the indels
                for a in alt_bases:
                    if alt_len != len(a):
                        problem_alts.append(newline)
                        break
                # check if the var is SNP
                if len(ref) == alt_len:
                    # SNPs must be length 1
                    if alt_len != 1:
                        problem_alts.append(newline)
                        break
                    for i in range(len(alt_bases)):
                        outlist = [chrom, pos, var_id, ref, alt_bases[i],
                                   qual, filt]
                        var_info = []
                        for col in info_cols:
                            try:
                                var_info.append(info_dict[col][i])
                            except KeyError:
                                var_info.append(np.nan)
                            except IndexError:
                                var_info.append(info_dict[col][0])
                        outlist = outlist + var_info
                        outfile_list.append(outlist)
                # if not a SNP, must be indel
                # indels must have their own line, hence only 1 indel in alt
                # bases
                elif len(alt_bases) > 1:
                    problem_alts.append(newline)
                    break
                # if conforming indel:
                else:
                    alt_base = alt_bases[0]
                    # vcf files have the  indels together with the preceding
                    # base such as REF: TA, ALT: T
                    # the same information is encoded as REF: A, ALT:- in table
                    if ref[0] != alt_base[0]:
                        problem_alts.append(newline)
                        break
                    else:
                        # remove preceding base
                        ref_base = ref[1:]
                        alt_bases = alt_base[1:]
                        # check if insertion
                        if ref_base == "":
                            ref_base = "-"
                            for i in range(2):
                                outlist = [chrom, pos + i, var_id,
                                           "-", alt_bases, qual, filt]
                                var_info = []
                                for col in info_cols:
                                    try:
                                        var_info.append(info_dict[col][0])
                                    except KeyError:
                                        var_info.append(np.nan)
                                outlist = outlist + var_info
                                outfile_list.append(outlist)
                        # if deletion
                        else:
                            # increment position because pos is not the
                            # effected position (ATT-> A, A is not affected)
                            pos += 1
                            for i in range(len(ref_base)):
                                outlist = [chrom, pos + i, var_id,
                                           ref_base[i], "-", qual, filt]
                                var_info = []
                                for col in info_cols:
                                    try:
                                        var_info.append(info_dict[col][0])
                                    except KeyError:
                                        var_info.append(np.nan)
                                outlist = outlist + var_info
                                outfile_list.append(outlist)
        var_df = pd.DataFrame(
            outfile_list, columns=["CHROM", "POS", "ID", "REF", "ALT",
                                   "QUAL", "FILTER"] + info_cols)
    var_df = var_df.astype({"AN": int, "AC": int})
    if len(problem_alts) > 0:
        print(("There are %d problematic alleles, see the output list for "
               "details" % len(problem_alts)))
    return var_df, problem_alts


def collapse_vcf_df(filt_df):
    """
    Take a vcf which has been converted to a Pandas data frame,
    groupby genomic position and add up the allele counts.
    """
    columns = filt_df.columns
    agg = {}
    for col in columns:
        if col not in ["AC", "AN", "CHROM", "POS"]:
            agg[col] = "first"
        elif col == "AC":
            agg[col] = np.sum
        elif col == "AN":
            agg[col] = np.max
    collapsed = filt_df.groupby(["CHROM", "POS"]).agg(agg).reset_index()
    return collapsed


def vcf_to_ucsc_table(collapsed, output_file):
    """
    Take a "per position" vcf dataframe, convert to UCSC genome browser
    style variant table.
    """
    # Alleles and allele counts are represented as A,T, and 2,20, in
    # the UCSC table. We'll create those strings with the following
    # functions

    def get_allele_strings(row):
        return row["REF"] + "," + row["ALT"] + ","

    def get_count_strings(row):
        ref_count = row["AN"] - row["AC"]
        return str(ref_count) + "," + str(row["AC"]) + ","
    collapsed["AS"] = collapsed.apply(get_allele_strings, axis=1)
    collapsed["CS"] = collapsed.apply(get_count_strings, axis=1)
    collapsed["0-offset"] = collapsed["POS"] - 1
    collapsed["BIN"] = 0
    table = collapsed[["BIN", "CHROM", "0-offset", "POS", "AS", "CS"]]
    table = table.sort_values(["CHROM", "0-offset"])
    tc = list(table.columns)
    tc_1 = tc[:4]
    tc_2 = tc[4:]
    for i in range(18):
        table[i] = ""
        tc_1.append(i)
    table = table[tc_1 + tc_2]
    table.to_csv(output_file, sep="\t", index=False, header=False)
    subprocess.call(["bgzip", "-c", output_file], stdout=open(
        output_file + ".gz", "w"))
    subprocess.call(["tabix", "-0", "-s 2", "-b 3", "-e 3",
                     output_file + ".gz"])
    return table


def header_to_primer(bc_dict,
                     header_string,
                     platform):
    """
    Convert a demultiplexed fastq header to forward and
    reverse primer numbers.
    """
    # Create sequence to primer dictionary from primer to sequence dict
    # bc_dict maps primer number to the sample barcode sequence
    # such as 1: AAATGCCC. We would like to get to a dict like AAATGCCC: 1
    seq_to_bc_dict = {v["sequence"]: int(k) for k, v in bc_dict.items()}
    split_string = header_string.split("+")
    if platform == "miseq":
        try:
            fw = seq_to_bc_dict[split_string[1]]
        except KeyError:
            fw = 999
        try:
            rev = seq_to_bc_dict[reverse_complement(split_string[0])]
        except KeyError:
            rev = 999
    elif platform == "nextseq":
        try:
            fw = seq_to_bc_dict[reverse_complement(split_string[1])]
        except KeyError:
            fw = 999
        try:
            rev = seq_to_bc_dict[reverse_complement(split_string[0])]
        except KeyError:
            rev = 999
    return fw, rev


def primer_to_header(bc_dict, primers, platform):
    """
    Convert fw, rev primer numbers to demultiplexed fastq header.
    """
    fw_seq = bc_dict[primers[0]]["sequence"]
    rev_seq = bc_dict[primers[1]]["sequence"]
    if platform == "nextseq":
        return reverse_complement(rev_seq) + "+" + reverse_complement(fw_seq)
    elif platform == "miseq":
        return reverse_complement(rev_seq) + "+" + fw_seq


def check_stitching(stitch_file):
    """
    Take a stitch log file from MIPWrangler output, return summary datframe.
    """
    with open(stitch_file) as infile:
        stitch = []
        for line in infile:
            newline = line.strip()
            stitch.append(newline)
    sti_sum = []
    for l in stitch:
        if '\t"stdOut_" : "[FLASH] Starting FLASH v' in l:
            nl = l.split("\\n")
            for i in range(len(nl)):
                il = nl[i]
                if "Input files" in il:
                    nil = nl[i+1].split(" ")
                    nil = [t for t in nil if t != ""]
                    sid = nil[-1].split("/")[1]
                    sti_sum.append(sid)
                elif (("Total pairs" in il) or
                      ("Combined pairs" in il)):
                    nil = il.split(" ")
                    nil = [t for t in nil if t != ""]
                    sti_sum.append(int(nil[-1]))
    sti = []
    for i in range(0, len(sti_sum), 3):
        sti.append(sti_sum[i:i+3])
    sti = pd.DataFrame(sti, columns=["Sample ID",
                                     "Total Reads",
                                     "Combined Reads"])
    return sti


def filter_vcf(in_vcf, out_vcf, filters_to_remove):
    """
    Filter a vcf (possibly gzipped) for given filters
    such that all variants containing any of the filters
    will be removed.
    """
    filt = set(filters_to_remove)
    # Check if the file is compressed
    try:
        input_vcf = gzip.open(in_vcf).readline()
        input_vcf = gzip.open(in_vcf)
        output_vcf = gzip.open(out_vcf, "w")
    except IOError:
        input_vcf = open(in_vcf)
        output_vcf = open(out_vcf, "w")
    with input_vcf as infile, output_vcf as outfile:
        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
            elif line.startswith("#"):
                cols = line.split("\t")
                for i in range(len(cols)):
                    if cols[i] == "FILTER":
                        filter_index = i
                        break
            else:
                newline = line.split("\t")
                var_filters = newline[filter_index].split(";")
                if len(filt.intersection(var_filters)) == 0:
                    outfile.write(line)
    return


def iupac_converter(iupac_code):
    """
    Return a list of all possible bases corresponding to a given iupac
    nucleotide code.
    """
    iupac_dict = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
                  "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
                  "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    try:
        return list(iupac_dict[iupac_code.upper()])
    except KeyError:
        print(("Non-IUPAC nucleotide code {}. Code must be one of {}").format(
            iupac_code, "".join(list(iupac_dict.keys()))
            ))
        return []


def make_degenerate(base_set):
    """
    Return IUPAC code of degenerate nucleotide corresponding to given base set.
    """
    iupac_dict = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
                  "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT",
                  "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    reverse_iupac = {frozenset(list(v)): k for k, v in iupac_dict.items()}
    try:
        base_set = frozenset(map(str.upper, frozenset(base_set)))
        return reverse_iupac[base_set]
    except (ValueError, KeyError, TypeError):
        return np.nan


def iupac_fasta_converter(header, sequence, max_ns=10):
    """Convert IUPAC degenerate nucleotides to ATGC.

    Given a sequence (header and sequence itself) containing iupac characters,
    return a dictionary with all possible sequences converted to ATCG.
    """
    iupac_dict = {"R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
                  "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
                  "N": "ACGT"}
    iupac_dict = {k: list(iupac_dict[k])
                  for k in list(iupac_dict.keys())}
    if sequence.upper().count("N") >= max_ns:
        return {header: sequence}
    sequence = list(sequence.upper())
    result_list = []

    def iupac_recurse(seq):
        for i in range(len(seq)):
            if seq[i] in list(iupac_dict.keys()):
                iup = iupac_dict[seq[i]]
                for i_seq in iup:
                    new_seq = copy.deepcopy(seq)
                    new_seq[i] = i_seq
                    iupac_recurse(new_seq)
                break
        else:
            result_list.append("".join(seq))
    iupac_recurse(sequence)
    if len(result_list) == 1:
        return {header: result_list[0]}
    else:
        return {header + "-" + str(i): result_list[i]
                for i in range(len(result_list))}


def save_fasta_dict(fasta_dict, fasta_file, linewidth=60):
    """Save a fasta dictionary to file."""
    with open(fasta_file, "w") as outfile:
        for header in fasta_dict:
            outfile.write(">" + str(header) + "\n")
            fasta_seq = fasta_dict[header]
            for i in range(0, len(fasta_seq), linewidth):
                outfile.write(fasta_seq[i: i + linewidth] + "\n")


def generate_sample_sheet(sample_list_file,
                          barcode_dict_file,
                          sample_sheet_template,
                          platform,
                          output_dir,
                          warnings=False):
    """Create a sample sheet to be used by bcl2fasq file from sample list."""
    with open(barcode_dict_file, "rb") as in1:
        barcode_dic = pickle.load(in1)
    # read in sample information
    sample_names = []
    sample_info = {}
    with open(sample_list_file) as infile:
        linenum = 0
        for line in infile:
            newline = line.strip().split("\t")
            # first line is the header with column names
            if linenum == 0:
                colnames = newline
                linenum += 1
            else:
                sample_dict = {colname: colvalue for colname, colvalue
                               in zip(colnames, newline)}
                sample_set = sample_dict["sample_set"]
                sample_name = sample_dict["sample_name"]
                replicate_number = sample_dict["replicate"]
                forward_index = sample_dict["fw"]
                reverse_index = sample_dict["rev"]
                sample_id = "-".join([sample_name,
                                      sample_set,
                                      replicate_number])
                if sample_id in sample_info:
                    print("Repeating sample name ", sample_id)
                if not sample_id.replace("-", "").isalnum():
                    print(("Sample IDs can only contain "
                           "alphanumeric characters and '-'. "
                           "{} has invalid characters.").format(sample_id))
                    continue
                # nextseq and miseq barcodes are handled differently
                if platform == "nextseq":
                    sample_dict.update(
                        {"i7": barcode_dic[reverse_index]["index_sequence"],
                         "i5": barcode_dic[forward_index]["index_sequence"]})
                elif platform == "miseq":
                    sample_dict.update(
                        {"i7": barcode_dic[reverse_index]["index_sequence"],
                         "i5": barcode_dic[forward_index]["sequence"]})
                sample_dict["sample_index"] = linenum
                linenum += 1
                sample_info[sample_id] = sample_dict
                sample_names.append(sample_id)
    # Check for samples sharing one or both barcodes. One barcode sharing is
    # allowed but a warning can be printed if desired by setting the warning
    #  to True. If both barcodes are shared among two samples, those samples
    # will be ignored and a message will be broadcast.
    samples_sharing = []
    for s1 in sample_info:
        for s2 in sample_info:
            if s1 != s2:
                if ((sample_info[s1]["fw"] == sample_info[s2]["fw"])
                   and (sample_info[s1]["rev"] == sample_info[s2]["rev"])):
                    samples_sharing.append([s1, s2])
                elif warnings and (
                    (sample_info[s1]["fw"] == sample_info[s2]["fw"])
                    or (sample_info[s1]["rev"] == sample_info[s2]["rev"])
                ):
                    print("Samples %s and %s share a barcode" % (s1, s2))
    samples_sharing_set = []
    if len(samples_sharing) > 0:
        for s in samples_sharing:
            samples_sharing_set.extend(s)
        samples_sharing_set = set(samples_sharing_set)
        print("There are %d samples sharing the same barcode pair"
              % len(samples_sharing_set))
        pd.DataFrame(samples_sharing).to_csv(
            os.path.join(output_dir, "samples_sharing_barcodes.tsv"),
            sep="\t"
        )
    # create sample sheet
    sample_sheet = os.path.join(output_dir, "SampleSheet.csv")
    with open(sample_sheet_template) as infile, \
            open(sample_sheet, "w") as outfile:
        outfile_list = infile.readlines()
        outfile_list = [o.strip() for o in outfile_list]
        for sample_id in sample_names:
            if sample_id in samples_sharing_set:
                continue
            reverse_index = sample_info[sample_id]["rev"]
            forward_index = sample_info[sample_id]["fw"]
            sample_index = str(sample_info[sample_id]["sample_index"])
            outlist = [sample_index, sample_id, "", "",
                       "S" + reverse_index,
                       sample_info[sample_id]["i7"],
                       "N" + forward_index,
                       sample_info[sample_id]["i5"], "", ""]
            outfile_list.append(",".join(outlist))
        outfile.write("\n".join(outfile_list))


def chromosome_converter(chrom, from_malariagen):
    """ Convert plasmodium chromosome names from standard (chr1, etc) to
    malariagen names (Pf3d7...) and vice versa.
    """
    standard_names = ["chr" + str(i) for i in range(1, 15)]
    standard_names.extend(["chrM", "chrP"])
    malariagen_names = ["Pf3D7_0" + str(i) + "_v3" for i in range(1, 10)]
    malariagen_names = malariagen_names + [
        "Pf3D7_" + str(i) + "_v3" for i in range(10, 15)]
    malariagen_names.extend(["Pf_M76611", "Pf3D7_API_v3"])
    if from_malariagen:
        return dict(zip(malariagen_names, standard_names))[chrom]
    else:
        return dict(zip(standard_names, malariagen_names))[chrom]


def write_list(alist, outfile_name):
    """ Convert values of a list to strings and save to file."""
    with open(outfile_name, "w") as outfile:
        outfile.write("\n".join(["\t".join(map(str, l))
                                for l in alist]) + "\n")
    return


##########################################################
# Core/shared functions
##########################################################


def strip_fasta(sequence):
    seq_list = sequence.split('\n')[1:]
    seq_join = "".join(seq_list)
    return seq_join


def calculate_gc(sequence, fasta=0):
    if fasta:
        seq_list = sequence.split('\n')[1:]
        seq_join = "".join(seq_list)
        seq = seq_join.lower()

    else:
        seq = sequence.lower()
    gc_count = seq.count('g') + seq.count('c')
    at_count = seq.count('a') + seq.count('t')
    percent = int(gc_count * 100 / (gc_count + at_count))
    return percent


def translate(sequence, three_letter=False):
    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
    gencode3 = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
                'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
                'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
                'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
                '*': '*'}
    seq = sequence.upper()
    # Return the translated protein from 'sequence' assuming +1 reading frame
    if not three_letter:
        return ''.join([gencode.get(seq[3*i:3*i+3], 'X')
                        for i in range(len(sequence)//3)])
    else:
        return ''.join([gencode3.get(gencode.get(seq[3*i:3*i+3], 'X'), "X")
                        for i in range(len(sequence)//3)])


def aa_converter(aa_name):
    """
    Output 3 letter and 1 letter amino acid codes for a given
    3 letter or 1 letter amino acid code.
    """
    gencode3 = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
                'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
                'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
                'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'}
    for a in list(gencode3.keys()):
        gencode3[gencode3[a]] = a
    return gencode3[aa_name.capitalize()]


def ntthal(s1, s2, Na=25, Mg=10, conc=0.4, print_command=False,
           td_path="/opt/resources/primer3_settings/primer3_config/"):
    """ Return the melting temperature of two oligos at given conditions,
    using ntthal from primer3 software.

    Parameters
    -----------
    s1 : str, sequence of first oligo.
    s2 : str, sequence of second oligo
    Na : int, Sodium (or other monovalent cation) concentration in mM
    Mg : int, Magnesium (or other divalent cation) concentration in mM
    conc : float, concentration of the more concentrated oligo in nM
    td_path : str, path to thermodynamic alignment parameters.
    """
    cmnd = ["ntthal", "-path", td_path, "-mv", str(Na), "-dv", str(Mg),
            "-d", str(conc), "-s1", s1, "-s2", s2, "-r"]
    if print_command:
        return(" ".join(cmnd))
    else:
        ntt_res = subprocess.check_output(cmnd)
        return float(ntt_res.decode("UTF-8").strip())


def oligoTM(s, Na=25, Mg=10, conc=0.4,
            thermodynamic_parameters=1, salt_correction=2):
    """ Return the melting temperature an oligo at given conditions,
    using oligotm from primer3 software.

    Parameters
    -----------
    s : str, sequence of the oligo.
    Na : int, Sodium (or other monovalent cation) concentration in mM
    Mg : int, Magnesium (or other divalent cation) concentration in mM
    conc : float, concentration of the more concentrated oligo in nM
    tp : [0|1], Specifies the table of thermodynamic parameters and
                the method of melting temperature calculation:
                 0  Breslauer et al., 1986 and Rychlik et al., 1990
                    (used by primer3 up to and including release 1.1.0).
                    This is the default, but _not_ the recommended value.
                 1  Use nearest neighbor parameters from SantaLucia 1998
                    *THIS IS THE RECOMMENDED VALUE*
    sc : [0..2], Specifies salt correction formula for the melting
                 temperature calculation
                  0  Schildkraut and Lifson 1965, used by primer3 up to
                     and including release 1.1.0.
                     This is the default but _not_ the recommended value.
                  1  SantaLucia 1998
                     *THIS IS THE RECOMMENDED VAULE*
                  2  Owczarzy et al., 2004

    """
    ntt_res = subprocess.check_output(
        ["oligotm", "-mv", str(Na), "-dv", str(Mg),
         "-d", str(conc), "-tp", str(thermodynamic_parameters),
         "-sc", str(salt_correction), s])
    return float(ntt_res.decode("UTF-8").strip())


def tm_calculator(sequence, conc, Na, Mg, dNTP_conc=0):
    from math import log
    from math import sqrt
    monovalent_conc = Na/1000
    divalent_conc = Mg/1000
    oligo_conc = conc * pow(10, -9)
    parameters = {}
    parameters['AA'] = (-7900, -22.2, -1.0)
    parameters['AT'] = (-7200, -20.4, -0.88)
    parameters['AC'] = (-8400, -22.4, -1.44)
    parameters['AG'] = (-7800, -21.0, -1.28)

    parameters['TA'] = (-7200, -21.3, -0.58)
    parameters['TT'] = (-7900, -22.2, -1.0)
    parameters['TC'] = (-8200, -22.2, -1.3)
    parameters['TG'] = (-8500, -22.7, -1.45)

    parameters['CA'] = (-8500, -22.7, -1.45)
    parameters['CT'] = (-7800, -21.0, -1.28)
    parameters['CC'] = (-8000, -19.9, -1.84)
    parameters['CG'] = (-10600, -27.2, -2.17)

    parameters['GA'] = (-8200, -22.2, -1.3)
    parameters['GT'] = (-8400, -22.4, -1.44)
    parameters['GC'] = (-9800, -24.4, -2.24)
    parameters['GG'] = (-8000, -19.9, -1.84)
    params = parameters
    # Normalize divalent_conc (Mg) for dNTP_conc
    K_a = 30000
    D = ((K_a * dNTP_conc - K_a * divalent_conc + 1) ** 2
         + 4 * K_a * divalent_conc)
    divalent_conc = (- (K_a * dNTP_conc - K_a * divalent_conc + 1)
                     + sqrt(D)) / (2 * K_a)

    # Define a, d, g coefficients used in salt adjustment
    a_con = 3.92 * (
        0.843 - 0.352 * sqrt(monovalent_conc) * log(monovalent_conc)
    )
    d_con = 1.42 * (
        1.279 - 4.03 * pow(10, -3) * log(monovalent_conc)
        - 8.03 * pow(10, -3) * ((log(monovalent_conc))**2)
    )
    g_con = 8.31 * (
        0.486 - 0.258 * log(monovalent_conc)
        + 5.25 * pow(10, -3) * ((log(monovalent_conc))**3)
    )
    dHsum = 0
    dSsum = 0
    sequence = sequence.upper()
    # define duplex initiation values for T and G terminal nucleotides
    if sequence[-1] == 'G' or sequence[-1] == 'C':
        dHiTer = 100
        dSiTer = -2.8
    elif sequence[-1] == 'A' or sequence[-1] == 'T':
        dHiTer = 2300
        dSiTer = 4.1
    if sequence[0] == 'G' or sequence[0] == 'C':
        dHiIn = 100
        dSiIn = -2.8
    elif sequence[0] == 'A' or sequence[0] == 'T':
        dHiIn = 2300
        dSiIn = 4.1
    dHi = dHiTer + dHiIn
    dSi = dSiTer + dSiIn

    R = 1.987  # ideal gas constant
    for i in range(len(sequence)-1):
        dinuc = sequence[i:(i+2)]
        dinuc_params = params[dinuc]
        dH = dinuc_params[0]
        dS = dinuc_params[1]
        dHsum += dH
        dSsum += dS
    # Tm w/o salt adjustment
    Tm = (dHsum + dHi)/float(dSsum + dSi + (R*log(oligo_conc)))

    # Salt adjustment
    GC_frac = calculate_gc(sequence)/100
    seq_length = len(sequence)
    if sqrt(divalent_conc)/monovalent_conc < 0.22:
        Tm = (Tm /
              (pow(10, -5) * Tm * ((4.29 * GC_frac - 3.95)
                                   * log(monovalent_conc)
                                   + 0.94 * (log(monovalent_conc)**2))
               + 1)
              )

    elif sqrt(divalent_conc)/monovalent_conc <= 6:
        Tm = (Tm /
              (Tm * (a_con
                     - 0.911 * log(divalent_conc)
                     + GC_frac * (6.26 + d_con * log(divalent_conc))
                     + (1 / float(2 * (seq_length - 1))) *
                     (-48.2 + 52.5 * log(divalent_conc) + g_con *
                      (log(divalent_conc)) ** 2))
               * pow(10, -5) + 1))
    elif sqrt(divalent_conc)/monovalent_conc > 6:
        a_con = 3.92
        d_con = 1.42
        g_con = 8.31
        Tm = (Tm /
              (Tm * (a_con
                     - 0.911 * log(divalent_conc)
                     + GC_frac * (6.26 + d_con * log(divalent_conc))
                     + (1 / (2 * float(seq_length - 1))) *
                     (-48.2 + 52.5 * log(divalent_conc) + g_con *
                      (log(divalent_conc)) ** 2))
               * pow(10, -5) + 1))
    return Tm - 273.15


def reverse_complement(sequence):
    """ Return reverse complement of a sequence. """
    complement_bases = {
        'g': 'c', 'c': 'g', 'a': 't', 't': 'a', 'n': 'n',
        'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N', "-": "-",
        "R": "Y", "Y": "R", "S": "W", "W": "S", "K": "M", "M": "K",
        "B": "V", "V": "B", "D":  "H", "H":  "D",
        "r": "y", "y": "r", "s": "w", "w": "s", "k": "m", "m": "k",
        "b": "v", "v": "b", "d":  "h", "h":  "d"
    }

    bases = list(sequence)
    bases.reverse()
    revcomp = []
    for base in bases:
        try:
            revcomp.append(complement_bases[base])
        except KeyError:
            print("Unexpected base encountered: ", base, " returned as X!!!")
            revcomp.append("X")
    return "".join(revcomp)


def get_file_locations():
    """ All static files such as fasta genomes, snp files, etc. must be listed
    in a file in the working directory. File name is file_locations.
    It is a tab separated text file. First tab has 2 letter species name, or
    "all" for general files used for all species. Second tab is the file name
    and third is the location of the file, either relative to script working
    directory, or the absolute path."""
    file_locations = {}
    with open("/opt/species_resources/file_locations.tsv", "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                if newline[0] not in list(file_locations.keys()):
                    file_locations[newline[0]] = {newline[1]: newline[2]}
                else:
                    file_locations[newline[0]][newline[1]] = newline[2]
    return file_locations


def id_generator(N):
    """ Generate a random string of length N consisting of uppercase letters
    and digits. Used for generating names for temporary files, etc.
    """
    return ''.join(random.SystemRandom().choice(
        string.ascii_uppercase + string.digits) for _ in range(N))


def alphanumerize(text, allowed_chars=["-"], replacement_char="-",
                  verbose=False):
    """Replace special characters, spaces etc in a text.

    Replace characters which are not alphanumeric or in the allowed characters
    list with the provided replacement character.
    """
    clist = []
    for ch in text:
        if ch.isalnum() or (ch in allowed_chars):
            clist.append(ch)
        else:
            clist.append(replacement_char)
    newtext = "".join(clist)
    if verbose and (newtext != text):
        print(("{} is replaced with {}.").format(text, newtext))
    return newtext

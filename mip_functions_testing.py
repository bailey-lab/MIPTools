import subprocess
import json
import os
import io
from multiprocessing import Pool, cpu_count
import multiprocessing
import multiprocessing.pool
from operator import itemgetter
import shutil
import random
import string
from ast import literal_eval as evaluate
import pickle
import copy
from Bio import SeqIO
import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth, DBSCAN
import matplotlib.pyplot as plt
from matplotlib import colors
from sklearn.manifold import TSNE
from scipy.stats import chi2_contingency, fisher_exact
from itertools import cycle
from functools import partial
import pysam
import mip_mod_testing as mod
import pandas as pd
import gzip
print("functions reloading")


# > Below class allows processors from a pool from multiprocessing module to create processor pools of their own.
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
def get_file_locations():
    """ All static files such as fasta genomes, snp files, etc. must be listed
    in a file in the working directory. File name is file_locations.
    It is a tab separated text file. First tab has 2 letter species name, or "all"
    for general files used for all species. Second tab is the file name and third is
    the location of the file, either relative to script working directory, or
    the absolute path."""
    file_locations = {}
    with open("resources/file_locations", "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                if newline[0] not in list(file_locations.keys()):
                    file_locations[newline[0]] = {newline[1]:newline[2]}
                else:
                    file_locations[newline[0]][newline[1]] = newline[2]
    return file_locations
def rinfo_converter(rinfo_file, output_file, flank, species="hs", host_species="none",
                    capture_type="targets", target_snps="functional", must_snps=[],
                    target_diffs="filtered", target_dic=None):
    """ Convert 0 offset rinfo file to 1 offset. In addition, make some
    formatting changes such as changing start/stop position names from
    seqb to begin and from seqe to end etc. More importantly, add all capture information
    that will be used by the mip module.
    capture_type options : "targets", "exons", "whole",
    target_snps options: "functional" (human), "functional_pf", "none" (both human and pf)
    target_diffs options: "must"(must only), "exonic", (must + exonic), "filtered" (must+exonic+filtered)
    must_snps: a list of dictionaries with following keys;
        "name", "snp_id", "segment", "chrom", "begin", "end", "weight"
    """
    target_snp_options = {"functional":
                      {"names": "nonsense,missense,stop-loss,frameshift,cds-indel,splice-3,splice-5",
                       "scores": "50,50,50,50,50,50,50"},
                   "none":{"names": "none", "scores":"none"},
                   "all": {"names": "all", "scores": "50"},
                   "functional_pf":{"names":"nonsynonymous", "scores":"50"},
                   "all_pf":{"names":"genic,non-genic,exon,intron,nonsynonymous,synonymous",
                             "scores":"100,100,100,100,100,100"}
                   }
    target_diffs_options_segments = {"filtered":{"names": "must,exonic_diffs,filtered_diffs",
                                       "scores": "20000,200,50"},
                           "exonic":{"names": "must,exonic_diffs",
                                     "scores": "20000,200"},
                           "must":{"names":"must", "scores":"200"},
                           "none":{"names": "none", "scores":"none"}
                   }
    target_diffs_options_genes = {"filtered":{"names": "must,exonic_diffs",
                                       "scores": "20000,20"},
                           "exonic":{"names": "must,exonic_diffs",
                                     "scores": "20000,20"},
                           "must":{"names":"must", "scores":"2000"},
                           "none":{"names": "none", "scores":"none"}
                   }
    settings = {}
    settings_file = "resources/settings_" + species
    with open(settings_file) as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                if newline[0] not in settings:
                    settings[newline[0]] = {}
                if newline[2] == "none":
                    settings[newline[0]][newline[1]] = newline[2]
                else:
                    try:
                        settings[newline[0]][newline[1]] = evaluate(newline[2])
                    except Exception as e:
                        settings[newline[0]][newline[1]] = newline[2]
    """
    segment_capture = {
                    "capture_type": capture_type,
                    "maf_for_arm_design": 0.01,
                    "mta_for_arm_design":100,
                    "maf_for_targeting":0.01,
                    "mta_for_targeting":100,
                    "target_snp_functions":target_snp_options[target_snps]["names"],
                    "score_snp_functions": target_snp_options[target_snps]["scores"],
                    "target_diffs":target_diffs_options[target_diffs]["names"],
                    "score_target_diffs":target_diffs_options[target_diffs]["scores"],
                    "min_mips":0,
                    "max_mips":50,
                    "flank":75,
                    "mask_diffs_lig":1,
                    "mask_diffs_ext":1,
                    "mask_snps_lig":0,
                    "mask_snps_ext":0}
    gene_capture = {
                    "capture_type": capture_type,
                    "maf_for_arm_design": 0.01,
                    "mta_for_arm_design":100,
                    "maf_for_targeting":0.01,
                    "mta_for_targeting":100,
                    "target_snp_functions":target_snp_options[target_snps]["names"],
                    "score_snp_functions": target_snp_options[target_snps]["scores"],
                    "target_diffs":target_diffs_options[target_diffs]["names"],
                    "score_target_diffs":target_diffs_options[target_diffs]["scores"],
                    "min_mips":0,
                    "max_mips":50,
                    "flank":75,
                    "mask_diffs_lig":0,
                    "mask_diffs_ext":0,
                    "mask_snps_lig":0,
                    "mask_snps_ext":0}
    selection_settings = {
                         "type": "compatibility",
                         "low": 1000,
                         "high": 3000,
                         "mip_limit": 50,
                         "trim_size": 10,
                         "trim_increment": 10,
                         "trim_limit": 100}

    """
    with open(rinfo_file, "r") as filein:
        infile = filein.readlines()
        segments = []
        segment_dic = {}
        for line in infile:
            if line.startswith("BASE"):
                newline = line.strip().split('\t')
                if newline[1] != "0":
                    print("File is not zero offset, exiting.")
                    return
        with open(output_file, "w") as outfile:
            for line in infile:
                if line.startswith("BASE"):
                    newline = line.strip().split('\t')
                    newline[1] = "1"
                    outline = "\t".join(newline) + '\n'
                    outfile.write(outline)
                elif line.startswith("COL:REGION"):
                    outlist = ["COL:REGION", "region_name", "copyname", "chrom", "begin",                               "end", "orient", "length"]
                    outline = "\t".join(outlist) + '\n'
                    outfile.write(outline)
                elif line.startswith("REGION"):
                    newline = line.strip().split("\t")
                    outlist = []
                    outlist.append(newline[0].split(":")[0])
                    outlist.append(":".join(newline[0].split(":")[1:]))
                    segmentname = newline[0].split(":")[1]
                    copynumber = newline[0].split(":")[2]
                    if not segmentname in segment_dic:
                        segment_dic[segmentname] = [copynumber]
                    else:
                        segment_dic[segmentname].append(copynumber)
                    segments.append(newline[0].split(":")[1])
                    outlist.extend(newline[1:3])
                    outlist.append(str(int(newline[3]) + 1))
                    outlist.extend(newline[4:])
                    outline = "\t".join(outlist) + '\n'
                    outfile.write(outline)
                else:
                    outfile.write(line)
            sep = "#####################################################################\n"
            comment = "# comment here\n"
            outfile.write(sep)
            outfile.write("SPECIES\t" + species + "\n")
            outfile.write("HOST_SPECIES\t" + host_species + "\n")
            outfile.write(sep)
            must_list = ["COL:MUST", "name", "snp_id", "segment",
                         "chrom", "begin", "end", "weight"]
            outfile.write("\t".join(must_list) + "\n")
            outfile.write(comment)
            outfile.write(sep)
            for m in must_snps:
                must_line = ["MUST"]
                for i in must_list[1:]:
                    must_line.append(m[i])
                outfile.write("\t".join(map(str,must_line)) + "\n")
                outfile.write(comment)

            selection_list = ["COL:SELECTION", "type", "low", "high", "mip_limit",
                              "trim_size", "trim_increment", "trim_limit"]
            outfile.write("\t".join(selection_list) + "\n")
            outfile.write(comment)
            select_line = ["SELECTION:"]
            for i in selection_list[1:]:
                if i in settings["selection"]:
                    select_line.append(settings["selection"][i])
                elif i in settings["all"]:
                    select_line.append(settings["all"][i])
                else:
                    select_line.append(none)
            outfile.write("\t".join(map(str, select_line)) + "\n")
            outfile.write(sep)
            capture_list = ["COL:CAPTURE", "segment", "flank","capture_type",
                            "target_snp_functions", "score_snp_functions",
                            "target_diffs", "score_target_diffs",
                            "maf_for_arm_design", "mta_for_arm_design",
	                        "maf_for_targeting", "mta_for_targeting",
                            "min_mips", "max_mips",  "mask_diffs_lig",
                            "mask_diffs_ext", "mask_snps_lig", "mask_snps_ext"]
            outfile.write("\t".join(capture_list) + "\n")
            segment_list = list(segment_dic.keys())
            segment_list.sort()
            for seg in segment_list:
                if len(segment_dic[seg])> 1:
                    cap_line = ["CAPTURE:", seg, flank, capture_type,
                            target_snp_options[target_snps]["names"],
                            target_snp_options[target_snps]["scores"],
                            target_diffs_options_segments[target_diffs]["names"],
                            target_diffs_options_segments[target_diffs]["scores"]]
                    for i in capture_list[8:]:
                            if i in settings["segment_capture"]:
                                cap_line.append(settings["segment_capture"][i])
                            elif i in settings["all"]:
                                cap_line.append(settings["all"][i])
                            else:
                                cap_line.append("none")
                else:
                    cap_line = ["CAPTURE:", seg, flank, capture_type,
                            target_snp_options[target_snps]["names"],
                            target_snp_options[target_snps]["scores"],
                            target_diffs_options_genes[target_diffs]["names"],
                            target_diffs_options_genes[target_diffs]["scores"]]
                    for i in capture_list[8:]:
                        if i in settings["gene_capture"]:
                            cap_line.append(settings["gene_capture"][i])
                        elif i in settings["all"]:
                            cap_line.append(settings["all"][i])
                        else:
                            cap_line.append("none")

                outfile.write("\t".join(map(str, cap_line)) + "\n")
            outfile.write(comment)
            outfile.write(sep)
            settings_list = ["COL:SETTINGS", "settings_for", "settings_file",
                             "Na", "Mg", "oligo_conc", "tm", "hit_threshold",
                             "lower_tm", "lower_hit_threshold", "tm_diff",
                             "3p_identity", "primer3_output_name", "update_primers",
                             "filtered", "bin_size", "pick_size", "hairpin_tm",
                             "dg", "size_min", "size_max", "backbone", "alternative_arms",
                             "seed_len", "bowtie_mode", "hit_limit", "upper_hit_limit",
                             "local", "processors"]
            outfile.write("\t".join(settings_list) + "\n")
            outfile.write(comment)
            for i, j in zip(["extension", "ligation", "mip"],                     [settings["extension"], settings["ligation"], settings["mip"]]):
                set_line = ["SETTINGS", i]
                for k in settings_list[2:]:
                    if target_dic:
                        try:
                            set_line.append(target_dic[k])
                            continue
                        except KeyError as e:
                            pass
                    if k in j:
                        set_line.append(j[k])
                    elif k in settings["all"]:
                        set_line.append(settings["all"][k])
                    else:
                        set_line.append("none")
                outfile.write("\t".join(map(str, set_line)) + "\n")
            outfile.write(comment)
            outfile.write(sep)
    return
def merge_coordinates(coordinates, capture_size, min_region_size = 0):
    """ Create MIP targets starting from a snp file that is produced offline,
    usually from Annovar. This is a tab separated file with the following content:
    chr1	2595307	2595307	A	G	rs3748816.
    This can be generalized to any target with coordinates.
    """
    # create target regions to cover all snps
    # start by getting snps on same chromosome together
    chroms = {}
    for c in coordinates:
        chrom = coordinates[c]["chrom"]
        try:
            chroms[chrom].append([coordinates[c]["begin"],
                                  coordinates[c]["end"]])
        except KeyError as e:
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
            #region_name = "-".join(targets_in_region)
            region_name = targets_in_region[0]
            target_names[region_name] = targets_in_region
            r_start = reg[0]
            r_end = reg[1]
            r_len = r_end - r_start + 1
            if r_len < min_region_size:
                r_start -= int(min_region_size - r_len/2)
                r_end += int(min_region_size - r_len/2)
            target_coordinates[region_name] = [c, r_start, r_end]
    return target_coordinates, target_names
def create_target_fastas(res_dir, targets, species, flank):
    for t in targets:
        rk = targets[t][0] + ":" + str(targets[t][1] - flank + 1)           + "-" + str(targets[t][2] + flank)
        with open(res_dir + t + ".fa", "w") as outfile:
            outfile.write(get_fasta(rk, species, header = t))
    return
def coordinate_to_target(coordinates, snp_locations, capture_size):
    """ Create MIP targets starting from a snp file that is produced offline,
    usually from Annovar. This is a tab separated file with the following content:
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
        except KeyError as e:
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
                if (reference_snp_locations[s]["chrom"] == c) and                     (r[0] <= reference_snp_locations[s]["begin"] <=                     reference_snp_locations[s]["end"] <= r[1]):
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
                ref_chrom = ref_location["chrom"]
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
    target_names = {}
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
    """
    for t in target_coordinates:
        reg_key = create_region(*target_coordinates[t])
        fasta = get_fasta(reg_key, species, header = t)
        with open(resource_dir + t + ".fa", "w") as outfile:
            outfile.write(fasta)
    """
    return target_coordinates
def rsid_to_target(resource_dir, snp_file):
    """ Create MIP targets starting from a snp file that is produced offline,
    usually from Annovar. This is a tab separated file with the following content:
    chr1	2595307	2595307	A	G	rs3748816.
    This can be generalized to any target with coordinates.
    """
    # one snp can have multiple locations on the reference genome,
    # this can happen with snps in regions where there are multiple different
    # assemblies (HLA locus, for example). So first step is to get each of these
    # locations in the genome.
    snp_locations = {}
    capture_types = {}
    with io.open(resource_dir + snp_file, encoding="utf-8") as infile:
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
                    if (snp["begin"] == temp_dic["begin"]) and                         (snp["end"] == temp_dic["end"]) and                         (snp["chrom"] == temp_dic["chrom"]) and                         (snp["ref_base"] == temp_dic["ref_base"]):
                        snp["alt_bases"].append(temp_dic["alt_bases"][0])
                        break
                else:
                    # add the snp dict if the location is different than what is present
                    # in the location dict.
                    snp_locations[rsid].append(temp_dic)
            except KeyError as e:
                # add the new rsid to location dict if it is not already present
                snp_locations[rsid] = [temp_dic]
                capture_types[rsid] = newline[6]
    # one reference location for each snp is required
    # alternative assambly chromosomes have an underscore in their names,
    # so that will be utilized to get the location in the orignal assembly,
    # i.e. the chromosome that does not have the underscore (chr7 and not chr7_alt08)
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
                print("Short chromosome name not found! Please check the output list.")
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
            target_coordinates[gene] = {"chrom":e["chrom"],
                                        "begin": e["begin"],
                                        "end": e["end"]}
        except KeyError as e:
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
                target_coordinates[tar_name] = {"chrom":gene_exons["chrom"],
                                                "begin": e[0],
                                                "end": e[1]}
        else:
            for j in exons_wanted:
                try:
                    e = exons[j]
                    tar_name = "-".join(gene, "exon", str(j))
                    target_coordinates[tar_name] = {"chrom":gene_exons["chrom"],
                                                    "begin": e[0],
                                                    "end": e[1]}
                except IndexError as e:
                    print(("Exon ", j, " does not exist for gene ", gene))
    return target_coordinates
def parse_alignment(reg_file):
    """ Create a rinfo dictionary from a rinfo file."""
    reg_dic = {}
    with open(reg_file, "r") as infile:
        for line in infile:
            if line.startswith("REGION"):
                newline = line.strip().split("\t")
                #print newline
                key1 = newline[1].split(":")[0]
                key2 = newline[1].split(":")[1]
                if key1 not in reg_dic:
                    reg_dic[key1] = {key2:{"copyname":newline[2],
                                            "chr":int(newline[3][3:]),
                                            "begin":int(newline[4]),
                                             "end":int(newline[5]),
                                             "ori":(newline[6]=="F")}
                                     }
                else:
                    reg_dic[key1][key2] = {"copyname":newline[2],
                                            "chr":int(newline[3][3:]),
                                            "begin":int(newline[4]),
                                             "end":int(newline[5]),
                                             "ori":(newline[6]=="F")}
    return reg_dic
def id_generator(N):
    """ Generate a random string of length N consisting of uppercase letters and digits.
    Used for generating names for temporary files, etc."""
    return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(N))
def hybrid_TM (temp_dir, s1, s2, Na=0.025, Mg=0.01, conc=(0.4*pow(10,-9)), tm=65):
    """ Return the melting temperature of two oligos at given conditions, using melt.pl
    of UNAfold software suite."""
    # generate random file names for the sequences
    # because melt.pl requires sequences to be read from files
    f1 = id_generator(8)
    f2 = id_generator(8)
    with open(temp_dir + f1, "w") as out1:
        with open(temp_dir + f2, "w") as out2:
            out1.write(s1)
            out2.write(s2)
    t = subprocess.check_output(["melt.pl", "-n",  "DNA", "-t" , str(tm),                      "-N", str(Na), "-M", str(Mg), "-C", str(conc),                     f1, f2], cwd=temp_dir)
    return float(t.strip().split("\t")[-1])
def hybrid_TM_files (temp_dir, f1, f2, fo, Na=0.025, Mg=0.01, conc=(0.4*pow(10,-9)), tm=65):
    """ Return the melting temperature of two oligos at given conditions, using melt.pl
    of UNAfold software suite."""
    # generate random file names for the sequences
    # because melt.pl requires sequences to be read from files
    t = subprocess.check_output(["melt.pl", "-n",  "DNA", "-t" , str(tm),                      "-N", str(Na), "-M", str(Mg), "-C", str(conc),                     f1, f2], cwd=temp_dir)
    with open(temp_dir + fo, "w") as of:
        of.write(t)
    return t
def align_region_worker_for_design(l):
    try:
        # get parameters from the input list
        #####
        # first item is the query sequence, if a fasta file, should be file name excluding path
        # it can either be a genomic location (chr1:10-100) or a fasta file.
        # If the item does not start with chr, it is assumed to be a fasta file
        region_key = l[0]
        # second item holds the run directory for lastz, should contain
        resource_dir = l[1]
        output_file = l[2]
        # target file, if a fasta file, should be pathname, not file name
        # if starts with chr, a fasta file will be created in the resource directory
        # if it is self, each sequence in query fasta will be aligned to each sequence
        # in the query fasta
        target_fasta = l[3]
        # each action item will be appended to the target or query argument
        # within brackets. [unmask] and [multiple] are important target actions
        # unmask: allows starting alignments in masked(lowercase) parts of the target
        # multiple: indicates there are multiple sequences in the target file (e.g. chromosomes)
        target_actions = l[4]
        # query file is treated as multiple sequence file without the multiple action
        query_actions = l[5]
        # percent cutoff value for identity/coverage of query to target. This only affects
        # reporting and not the alignment process itself.
        identity_cutoff = l[6]
        coverage_cutoff = l[7]
        # format of the output, follows --format: argument in lastz
        # if format is general, it should be followed by a comma separated list of
        # fields to output, e.g. general:name1,text1,name2,text2,diff,score would output
        # the name of the query, sequence of the query, name of the target, seq of target,
        # a string showing the alignment and the alignment score
        output_format = l[8]
        # additional options to pass to lastz
        options = l[9]
        species = l[10]
        # create a query fasta file from region key, if a key is provided instead of fasta
        query_fasta = resource_dir + region_key + ".fa"
        if target_fasta == "self":
            target_fasta = query_fasta
        """
        if not os.path.exists(query_fasta):
            e = get_exons(get_gene(region_key,
                                   get_file_locations()[species]["refgene"],
                                   alternative_chr=1))
            rk = e["chrom"] + ":" + str(e["begin"]) + "-" +\
                str(e["end"])
            with open(query_fasta, "w") as outfile:
                outfile.write(get_fasta(rk, species, header = region_key))
        if target_fasta == "self":
            target_fasta = query_fasta
        # create a target fasta file from region key, if a key is provided instead of fasta
        elif target_fasta.startswith("chr"):
            target_key = str(target_fasta)
            target_fasta = resource_dir + target_fasta
            if not os.path.exists(target_fasta):
                with open(target_fasta, "w") as outfile:
                    outfile.write(get_fasta(target_key, species))
        """
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
                 query_fasta + query_act ,
                 "--output=" + resource_dir + output_file,
                 "--format=" + output_format,
                 "--filter=identity:" + str(identity_cutoff),
                 "--filter=coverage:" + str(coverage_cutoff)]
        # add any extra options to the end of the command
        comm.extend(options)
        #print " ".join(comm)
        # run the command using subprocess module
        subprocess.check_output(comm)
        return 0
    except Exception as e:
        return str(e)
def align_region_multi_for_design(alignment_list, pro):
    res = []
    try:
        p = Pool(pro)
        p.map_async(align_region_worker_for_design, alignment_list, callback=res.append)
        p.close()
        p.join()
    except Exception as e:
        res.append(str(e))
    return res
def align_region_worker(l):
    try:
        # get parameters from the input list
        #####
        # first item is the query sequence, if a fasta file, should be file name excluding path
        # it can either be a genomic location (chr1:10-100) or a fasta file.
        # If the item does not start with chr, it is assumed to be a fasta file
        region_key = l[0]
        # second item holds the run directory for lastz, should contain
        resource_dir = l[1]
        output_file = l[2]
        # target file, if a fasta file, should be pathname, not file name
        # if starts with chr, a fasta file will be created in the resource directory
        # if it is self, each sequence in query fasta will be aligned to each sequence
        # in the query fasta
        target_fasta = l[3]
        # each action item will be appended to the target or query argument
        # within brackets. [unmask] and [multiple] are important target actions
        # unmask: allows starting alignments in masked(lowercase) parts of the target
        # multiple: indicates there are multiple sequences in the target file (e.g. chromosomes)
        target_actions = l[4]
        # query file is treated as multiple sequence file without the multiple action
        query_actions = l[5]
        # percent cutoff value for identity/coverage of query to target. This only affects
        # reporting and not the alignment process itself.
        identity_cutoff = l[6]
        coverage_cutoff = l[7]
        # format of the output, follows --format: argument in lastz
        # if format is general, it should be followed by a comma separated list of
        # fields to output, e.g. general:name1,text1,name2,text2,diff,score would output
        # the name of the query, sequence of the query, name of the target, seq of target,
        # a string showing the alignment and the alignment score
        output_format = l[8]
        # additional options to pass to lastz
        options = l[9]
        species = l[10]
        # create a query fasta file from region key, if a key is provided instead of fasta
        if region_key.startswith("chr"):
            query_fasta = resource_dir + region_key + ".fa"
            if not os.path.exists(query_fasta):
                with open(query_fasta, "w") as infile:
                    infile.write(get_fasta(region_key, species))
        else:
            query_fasta = resource_dir + region_key
        if target_fasta == "self":
            target_fasta = query_fasta
        # create a target fasta file from region key, if a key is provided instead of fasta
        elif target_fasta.startswith("chr"):
            target_key = str(target_fasta)
            target_fasta = resource_dir + target_fasta
            if not os.path.exists(target_fasta):
                with open(target_fasta, "w") as infile:
                    infile.write(get_fasta(target_key, species))
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
                 query_fasta + query_act ,
                 "--output=" + resource_dir + output_file,
                 "--format=" + output_format,
                 "--filter=identity:" + str(identity_cutoff),
                 "--filter=coverage:" + str(coverage_cutoff)]
        # add any extra options to the end of the command
        comm.extend(options)
        #print comm
        # run the command using subprocess module
        subprocess.check_output(comm)
        return 0
    except Exception as e:
        return str(e)
def align_region_multi(alignment_list, pro):
    res = []
    try:
        p = Pool(pro)
        p.map_async(align_region_worker, alignment_list, callback=res.append)
        p.close()
        p.join()
    except Exception as e:
        res.append(str(e))
    return res
def align_region(l):
    try:
        region_key = l[0]
        resource_dir = l[1]
        output_file = l[2]
        target_fasta = l[3]
        target_actions = l[4]
        query_actions = l[5]
        identity_cutoff = l[6]
        coverage_cutoff = l[7]
        output_format = l[8]
        options = l[9]
        query_fasta = resource_dir + region_key + ".fa"
        if not os.path.exists(query_fasta):
            with open(query_fasta, "w") as infile:
                infile.write(get_fasta(region_key, "hs"))
        if target_fasta == "self":
            target_fasta = query_fasta
        elif target_fasta.startswith("chr"):
            target_key = str(target_fasta)
            target_fasta = resource_dir + target_fasta
            if not os.path.exists(target_fasta):
                with open(target_fasta, "w") as infile:
                    infile.write(get_fasta(target_key, "hs"))

        if len(target_actions) > 0:
            target_act = "[" + ",".join(target_actions) + "]"
        else:
            target_act = ""
        if len(query_actions) > 0:
            query_act = "[" + ",".join(query_actions) + "]"
        else:
            query_act = ""
        comm = ["lastz_32",
                                 target_fasta + target_act,
                                 query_fasta + query_act ,
                                 "--output=" + resource_dir + output_file,
                                 "--format=" + output_format,
                                 "--filter=identity:" + str(identity_cutoff),
                                 "--filter=coverage:" + str(coverage_cutoff)]
        comm.extend(options)
        subprocess.check_output(comm)
        return 0
    except Exception as e:
        return str(e)
def merge_alignments (resource_dir, fasta_list, output_prefix = "merged"):
    als_out = []
    diffs_out = []
    with open(resource_dir + output_prefix + ".al", "w") as alignment_file,         open(resource_dir + output_prefix + ".differences", "w") as diff_file:
        for f in fasta_list:
            fnum = 0
            with open(resource_dir + f + ".al") as alignment,             open(resource_dir + f + ".differences") as diffs:
                linenum = 0
                for line in alignment:
                    if linenum > 0:
                        als_out.append(line.strip())
                    elif fnum == 0:
                        als_out.append(line.strip())
                        linenum += 1
                    else:
                        linenum += 1
                for d in diffs:
                    diffs_out.append(d.strip())
            fnum += 0
        alignment_file.write("\n".join(als_out))
        diff_file.write("\n".join(diffs_out))
    return
def align_genes_for_design(gene_name_list, res_dir,
                           alignment_types = ["differences", "general"],
                species = "hs", num_processor=30 ):
    """ Align genes given in a list to the reference genome of given species.
    Extract gene sequence of the gene from refgene file of the species, starting
    from the first exon and ending with the last exon. Any alternatively spliced
    exons will be added. Exonic sequence will be flanked by specified number of
    bases. Each alignment type and gene will be passed to the align_region_multi
    function to be processed by a separate processor. Alignment results will be
    written to files in wdir + gene_name/resources. All alignments will be checked
    and if there are genes whose alignments overlap (e.g. duplications) they will
    be returned in a list.
    """
    region_list = []
    for gene_dict in gene_name_list:
        gene_name = gene_dict["gene_name"]
        identity = gene_dict["identity"]
        coverage = gene_dict["coverage"]
        options = gene_dict["options"]
        target = get_file_locations()[species]["fasta_genome"]
        out_fields = "name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,zstart2+,end2+,length2,identity,coverage"
        gen_out = "general:"+ out_fields
        dif_out = "differences"
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)
        if "general" in alignment_types:
            al = [gene_name,
              res_dir,
              gene_name + ".al",
              target,
              ["multiple", "unmask", "nameparse=darkspace"],
              ["unmask", "nameparse=darkspace"],
              identity,
              coverage,
              gen_out,
              options,
              species]
            region_list.append(al)
        if "differences" in alignment_types:
            al = [gene_name,
                  res_dir,
                  gene_name + ".differences",
                  target,
                  ["multiple", "unmask", "nameparse=darkspace"],
                  ["unmask", "nameparse=darkspace"],
                  identity,
                  coverage,
                  dif_out,
                  options,
                  species]
            region_list.append(al)
    align_region_multi_for_design(region_list, num_processor)
    #check = alignment_check(region_list)
    return
def align_genes(gene_name_list, wdir, alignment_types = ["differences", "general"],
                species = "hs", num_processor=30, flank=150, identity=90, coverage=5,
                options = ["--notransition", "--step=100"]):
    """ Align genes given in a list to the reference genome of given species.
    Extract gene sequence of the gene from refgene file of the species, starting
    from the first exon and ending with the last exon. Any alternatively spliced
    exons will be added. Exonic sequence will be flanked by specified number of
    bases. Each alignment type and gene will be passed to the align_region_multi
    function to be processed by a separate processor. Alignment results will be
    written to files in wdir + gene_name/resources. All alignments will be checked
    and if there are genes whose alignments overlap (e.g. duplications) they will
    be returned in a list.
    """
    region_list = []
    for gene_name in gene_name_list:
        if gene_name.startswith("chr"):
            query = gene_name
            res_dir = wdir
        elif gene_name.endswith(".fasta"):
            query = gene_name
            res_dir = wdir
        else:
            e = get_exons(
                get_gene(gene_name, get_file_locations()[species]["refgene"], alternative_chr=1)
                )
            query = e["chrom"] + ":" + str(e["begin"] - flank) + "-" + str(e["end"] + flank)
            res_dir = wdir + gene_name + "/resources/"
        target = get_file_locations()[species]["fasta_genome"]
        out_fields = "name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,zstart2+,end2+,length2,identity,coverage"
        gen_out = "general:"+ out_fields
        dif_out = "differences"
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)
        if "general" in alignment_types:
            al = [query,
              res_dir,
              gene_name + ".al",
              target,
              ["multiple", "unmask"],
              ["unmask"],
              identity,
              coverage,
              gen_out,
              options,
              species]
            region_list.append(al)
        if "differences" in alignment_types:
            al = [query,
                  res_dir,
                  gene_name + ".differences",
                  target,
                  ["multiple", "unmask"],
                  ["unmask"],
                  identity,
                  coverage,
                  dif_out,
                  options,
                  species]
            region_list.append(al)
    align_region_multi(region_list, num_processor)
    #check = alignment_check(region_list)
    return
def align_haplotypes(settings,
                     target_actions = ["unmask", "multiple"],
                     query_actions = ["unmask"],
                     output_format="general:name1,text1,name2,text2,diff,score",
                     alignment_options = ["--noytrim"], identity=75, coverage=75):
    """ Get a haplotypes dict and a call_info dict, align each haplotype to reference
    sequences from the call_info dict."""
    wdir = settings["workingDir"]
    haplotypes_file = wdir + settings["tempHaplotypesFile"]
    with open(haplotypes_file) as infile:
        haplotypes = json.load(infile)
    species  = settings["species"]
    alignment_dir = wdir + settings["alignmentDir"]
    num_processor = int(settings["processorNumber"])
    command_list = []
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    # create alignment dir if it does not exist
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    for m in haplotypes:
        # create a fasta file for each mip that contains all haplotype sequences
        # for that mip
        haplotype_fasta = alignment_dir + m + ".haps"
        with open(haplotype_fasta, "w") as outfile:
            outfile_list = []
            for h in haplotypes[m]:
                outlist = [">", h, "\n", haplotypes[m][h]["sequence"]]
                outfile_list.append("".join(outlist))
            outfile.write("\n".join(outfile_list))
        haplotype_fasta = m + ".haps"
        # create a reference file for each mip that contains reference sequences
        # for each paralog copy for that mip
        reference_fasta = alignment_dir + m + ".refs"
        with open(reference_fasta, "w") as outfile:
            outfile_list = []
            gene_name = m.split("_")[0]
            for c in call_info[gene_name][m]["copies"]:
                c_ori = call_info[gene_name][m]["copies"][c]["orientation"]
                c_seq = call_info[gene_name][m]["copies"][c]["capture_sequence"]
                if c_ori == "reverse":
                    c_seq = reverse_complement(c_seq)
                outlist = [">", m + "_" + c, "\n", c_seq]
                outfile_list.append("".join(outlist))
            outfile.write("\n".join(outfile_list))
        # name of the alignment output file for the mip
        output_file = m + ".aligned"
        # create the list to be passed to the alignment worker
        command = [haplotype_fasta, alignment_dir, output_file, reference_fasta,
                   target_actions, query_actions, identity, coverage, output_format,
                   alignment_options, species]
        # add the command to the list that will be passed to the multi-aligner
        command_list.append(command)
    # run the alignment
    alignment = align_region_multi(command_list, num_processor)
    alignment_out_file = wdir + settings["tempAlignmentStdOut"]
    with open(alignment_out_file, "w") as outfile:
        json.dump(alignment, outfile)
    return
def parse_aligned_haplotypes(settings):
    wdir = settings["workingDir"]
    species = settings["species"]
    alignment_dir = wdir + settings["alignmentDir"]
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    temp_haplotypes_file = wdir + settings["tempHaplotypesFile"]
    with open(temp_haplotypes_file) as infile:
        haplotypes = json.load(infile)
    alignments = {}
    inverted_alignments = []
    problem_alignments = []
    problem_snps = []
    for m in haplotypes:
        # each mip has all its haplotypes and reference sequences aligned
        # in a separate file
        with open(alignment_dir + m + ".aligned") as al_file:
            for line in al_file:
                problem_al = False
                if not line.startswith("#"):
                    # each line of the alignment file includes an alignment
                    # between a reference copy sequence of a mip
                    newline = line.strip().split("\t")
                    gene_name = newline[0].split("_")[0]
                    m_name = "_".join(newline[0].split("_")[:-1])
                    ref_copy = newline[0].split("_")[-1]
                    rf_ori = call_info[gene_name][m_name]["copies"][ref_copy]["orientation"]
                    # aligned part of the reference sequence with gaps
                    ref_al = newline[1].upper()
                    if rf_ori == "reverse":
                        ref_al = reverse_complement(ref_al)
                    # aligned part of the reference without gaps
                    ref_used = ref_al.translate(str.maketrans({"-" : None})).upper()
                    hap_name = newline[2]
                    # aligned part of the haplotype with gaps
                    hap_al = newline[3].upper()
                    if rf_ori == "reverse":
                        hap_al = reverse_complement(hap_al)
                    # aligned part of the haplotype without gaps
                    hap_used = hap_al.translate(str.maketrans({"-" : None})).upper()
                    # alignment (.for match, : and X mismatch, - gap)
                    diff = newline[4]
                    if rf_ori == "reverse":
                        diff = diff[::-1]
                    score = int(newline[5])
                    # full haplotype sequence
                    hap_seq = haplotypes[m][hap_name]["sequence"].upper()
                    # full reference sequence
                    ref_seq = call_info[gene_name][m]["copies"][ref_copy]["capture_sequence"].upper()
                    # index of where in full reference the alignment begins
                    ref_align_begin = ref_seq.find(ref_used)
                    # index of where in full reference the alignment ends
                    ref_align_end = ref_align_begin + len(ref_used)
                    # index of where in full haplotype sequence the alignment begins
                    hap_align_begin = hap_seq.find(hap_used)
                    # if the alignment is inverted, i.e. there is a reverse complement
                    # alignment with a significant score, the find method will not find
                    # the haplotype sequence in query or the target sequence in reference
                    # and return -1. These alignments have been happening when one copy
                    # differs so much from another, an inverted alignment scores better.
                    # These should be ignored as the real copy the haplotype comes from
                    # will have a better score. However, there can theoretically an inversion
                    # within a capture region that produces a legitimate inverted alignment
                    # therefore these alignments should be inspected afterwards.
                    if min([hap_align_begin, ref_align_begin]) < 0:
                        al_dict = {"gene_name": gene_name,
                               "mip_name": m,
                               "copy": ref_copy,
                               "score": score,
                               "aligned_hap": hap_al,
                               "aligned_ref": ref_al,
                               "diff": diff,
                               "haplotype_ID": hap_name}
                        inverted_alignments.append(al_dict)
                        continue
                    # index of where in full haplotype sequence the alignment ends
                    hap_align_end = hap_align_begin + len(hap_used)
                    # deal with any existing flanking deletions/insertions
                    # is there any unaligned sequence on the left of alignment
                    left_pad_len = max([hap_align_begin, ref_align_begin])
                    left_pad_diff = abs(hap_align_begin - ref_align_begin)
                    left_pad_ref = ""
                    left_pad_hap = ""
                    left_pad_ref_count = 0
                    left_pad_hap_count = 0
                    # where there are insertions on left, fill the other pad with gaps
                    for i in range(hap_align_begin - ref_align_begin):
                        # only when ref_align_begin is smaller, we need to pad left_pad_ref
                        left_pad_ref = "-" + left_pad_ref
                        left_pad_hap = hap_seq[i] + left_pad_hap
                        # counting how many bases from hap_seq is used for padding
                        left_pad_hap_count += 1
                    # do the same for haplotype sequence
                    for i in range(ref_align_begin - hap_align_begin):
                        # only when ref_align_begin is smaller, we need to pad left_pad_ref
                        left_pad_hap += "-"
                        left_pad_ref += ref_seq[i]
                        # counting how many bases from ref_seq is used for padding
                        left_pad_ref_count += 1
                    # add to left_pads the sequences which are there but did not align
                    for i in range(left_pad_len - left_pad_diff):
                        left_pad_ref += ref_seq[i + left_pad_ref_count]
                        left_pad_hap += hap_seq[i + left_pad_hap_count]
                    # add the left padding info to the alignment
                    for i in range(0,len(left_pad_hap))[::-1]:
                        if left_pad_ref[i] == "-" or left_pad_hap[i] == "-":
                            diff = "-" + diff
                        elif left_pad_ref[i] != left_pad_hap[i]:
                            diff = "X" + diff
                        else:
                            diff = "." + diff
                            problem_al = True
                    # repeat the padding for the right side of alignment
                    right_pad_ref_len = len(ref_seq) - ref_align_end
                    right_pad_hap_len = len(hap_seq) - hap_align_end
                    right_pad_len = max([right_pad_hap_len, right_pad_ref_len])
                    right_pad_diff = abs(right_pad_hap_len - right_pad_ref_len)
                    right_pad_ref = ""
                    right_pad_hap = ""
                    right_pad_ref_count = 0
                    right_pad_hap_count = 0
                    for i in range(right_pad_hap_len - right_pad_ref_len):
                        right_pad_ref = "-" + right_pad_ref
                        right_pad_hap = hap_seq[-i - 1] + right_pad_hap
                        # counting how many bases from hap_seq is used for padding
                        right_pad_hap_count += 1
                    # do the same for haplotype sequence
                    for i in range(right_pad_ref_len - right_pad_hap_len):
                        right_pad_hap = "-" + right_pad_hap
                        right_pad_ref = ref_seq[-i - 1] + right_pad_ref
                        right_pad_ref_count += 1
                    # add to right the sequences which are there but did not align
                    for i in range(right_pad_len - right_pad_diff):
                        right_pad_ref = ref_seq[-i - right_pad_ref_count -1] + right_pad_ref
                        right_pad_hap = hap_seq[-i - right_pad_hap_count -1] + right_pad_hap
                    # add the right padding info to the alignment
                    for i in range(len(right_pad_hap)):
                        if right_pad_ref[i] == "-" or right_pad_hap[i] == "-":
                            diff += "-"
                        elif right_pad_ref[i] != right_pad_hap[i]:
                            diff += "X"
                        else:
                            diff += "."
                            problem_al = True
                    hap_al = left_pad_hap + hap_al + right_pad_hap
                    ref_al = left_pad_ref + ref_al + right_pad_ref
                    # we have padded the alignment so now all the ref and
                    # hap sequence is accounted for and not just the aligned part
                    # ref_name, ref_copy, ref_seq, ref_al
                    # hap_name, hap_seq, hap_al, diff, score have information we'll use
                    c_name = ref_copy
                    h_name = hap_name
                    exact_match = (ref_seq == hap_seq)
                    try:
                        snp_dict = call_info[gene_name][m]["snps"]
                    except KeyError as e:
                        snp_dict = {}
                    copy_dict = call_info[gene_name][m]["copies"][c_name]
                    copy_ori = copy_dict["orientation"]
                    copy_chrom = copy_dict["chrom"]
                    copy_begin = int(copy_dict["capture_start"])
                    copy_end = int(copy_dict["capture_end"])
                    # if copy orientation is reverse, we'll reverse the alignment
                    # so that coordinate conversions are easier and indels are always
                    # left aligned on forward genomic strand
                    if copy_ori == "reverse":
                        ref_al = reverse_complement(ref_al)
                        hap_al = reverse_complement(hap_al)
                        diff = diff[::-1]
                    genomic_pos = copy_begin
                    differences = []
                    indel_count = 0
                    indels = []
                    indel_types = []
                    # keep track of the index of haplotype sequence
                    # to use for checking sequence quality later
                    hap_index = 0
                    for i in range(len(diff)):
                        d = diff[i]
                        # each difference between the hap and ref can be an indel ("-")
                        # or a snp (":" or "x") or the same as the reference (".")
                        # When dealing with indels, it is best to call
                        # consecutive indels as a cumulative indel rather than individual
                        # i.e. AAA/--- instead of A/-, A/-, A/- because if we are looking for
                        # a frameshift insertion A/-, having AAA/--- means we don't observe the
                        # frameshift. But if it is kept as three A/-'s then it looks like
                        # the frameshift mutation is there.
                        if d == "-":
                            # if an indel is encountered, we'll keep track of it until
                            # the end of the indel. That is, when d != "-"
                            indel_count += 1
                            if hap_al[i] == "-":
                                # if a deletion, hap sequence should have "-"
                                indel_types.append("del")
                                indels.append(ref_al[i])
                                # in cases of deletions, we increment the genomic pos
                                # because the reference has a nucleotide in this position
                                genomic_pos += 1
                            elif ref_al[i] == "-":
                                indel_types.append("ins")
                                indels.append(hap_al[i])
                                hap_index += 1
                                # in cases of insertions, we don't increment the genomic pos
                                # because the reference has no nucleotide in this position
                                # insAAA would have same start and end positions
                            else:
                                # if neither hap nor ref has "-" at this position there is
                                # a disagreement between the alignment and the sequences
                                print("Alignment shows indel but sequences do not", h_name)
                                break
                        else:
                            # if the current diff is not an indel,
                            # check if there is preceeding indel
                            if len(indels) > 0:
                                # there should only be a del or ins preceding this base
                                if len(set(indel_types)) != 1:
                                    # Consecutive insertions and deletions
                                    problem_al = True
                                    break
                                else:
                                    indel_type = list(set(indel_types))[0]
                                    indel_length = len(indels)
                                    # genomic_pos is the current position
                                    # since this position is not an indel, indel has ended
                                    # 1 nucleotide prior to this position.
                                    indel_end = genomic_pos - 1
                                    indel_seq = "".join(indels)
                                    buffer_seq = "".join(["-" for j in range(indel_length)])
                                    if indel_type == "del":
                                        indel_begin = genomic_pos - indel_length
                                        ref_base = indel_seq
                                        hap_base = buffer_seq
                                        h_index = [hap_index, hap_index - 1]
                                    else:
                                        # if the preceding indel was an insertion
                                        # the start and end position is the same
                                        indel_begin = genomic_pos - 1
                                        ref_base = buffer_seq
                                        hap_base = indel_seq
                                        h_index = [hap_index - indel_length, hap_index - 1]
                                # create an indel dict and add to differences list
                                differences.append({"begin": indel_begin,
                                                    "end": indel_end,
                                                    "type": indel_type,
                                                    "ref_base": ref_base,
                                                    "hap_base": hap_base,
                                                    "hap_index": h_index,
                                                    "chrom": copy_chrom})
                            # clean up the indel variables
                            indel_count = 0
                            indels = []
                            indel_types = []
                            # continue with the current snp
                            if d == ".":
                                # "." denotes hap and ref has the same sequence
                                pass
                            else:
                                # create a snp dict and add to differences list
                                ref_base = ref_al[i]
                                hap_base = hap_al[i]
                                h_index = [hap_index, hap_index]
                                differences.append({"begin": genomic_pos,
                                                    "end": genomic_pos,
                                                    "type": "snp",
                                                    "ref_base": ref_base,
                                                    "hap_base": hap_base,
                                                    "hap_index": h_index,
                                                    "chrom": copy_chrom})
                            hap_index += 1
                            genomic_pos += 1

                    # since indel dicts are not created until a non-indel character
                    # is encountered, we need to check if there was an indel at the
                    # end of the alignment. If there is, indels list would not have been reset.
                    # check if there is preceeding indel
                    if len(indels) > 0:
                        if len(set(indel_types)) != 1:
                            # Consecutive insertions and deletions
                            problem_al = True
                        else:
                            indel_type = list(set(indel_types))[0]
                            indel_length = len(indels)
                            indel_end = genomic_pos - 1
                            indel_seq = "".join(indels)
                            buffer_seq = "".join(["-" for idl in range(indel_length)])
                            if indel_type == "del":
                                indel_begin = genomic_pos - indel_length
                                ref_base = indel_seq
                                hap_base = buffer_seq
                                h_index = [hap_index, hap_index - 1]
                            else:
                                indel_begin = genomic_pos - 1
                                ref_base = buffer_seq
                                hap_base = indel_seq
                                h_index = [hap_index - indel_length, hap_index - 1]
                        differences.append({"begin": indel_begin,
                                            "end": indel_end,
                                            "type": indel_type,
                                            "ref_base": ref_base,
                                            "hap_base": hap_base,
                                            "hap_index": h_index,
                                            "chrom": copy_chrom})
                        # clean up the indel variables
                        indel_count = 0
                        indels = []
                        indel_types = []
                    # fix the positioning of homopolymer indels
                    if copy_ori == "reverse":
                        ref_seq = reverse_complement(ref_seq)
                    for d in differences:
                        d_chrom = d["chrom"]
                        d_pos = int(d["begin"])
                        ref_base = d["ref_base"].upper()
                        hap_base = d["hap_base"].upper()
                        d_type = d["type"]
                        if d_type in ["ins", "del"]:
                            if d_type == "del":
                                d_pos -= 1
                            d_prior_index = d_pos - copy_begin
                            if d_prior_index >= 0:
                                prior_base = ref_seq[d_prior_index]
                            else:
                                prior_base = get_sequence(d_chrom + ":" + str(d_pos)                                                        + "-" + str(d_pos), species).upper()
                            vcf_ref = prior_base + ref_base
                            vcf_hap = prior_base + hap_base
                            vcf_ref = "".join([b for b in vcf_ref if b != "-"])
                            vcf_hap = "".join([b for b in vcf_hap if b != "-"])
                        else:
                            vcf_ref = ref_base
                            vcf_hap = hap_base
                        vcf_key = d_chrom + ":" + str(d_pos) + ":.:" + vcf_ref + ":" + vcf_hap
                        d["vcf_raw"] = vcf_key
                    # all differences in between the ref and hap has been documented
                    # loop through differences and assign psv/clinical values
                    for d in differences:
                        diff_begin = d["begin"]
                        diff_end = d["end"]
                        if copy_ori == "reverse":
                            # revert bases to their original strand (-)
                            d["ref_base"] = reverse_complement(d["ref_base"])
                            d["hap_base"] = reverse_complement(d["hap_base"])
                            d["hap_index"] = [-1 * d["hap_index"][0] - 1,
                                              -1 * d["hap_index"][1] - 1]
                        diff_ref_base = d["ref_base"]
                        diff_hap_base = d["hap_base"]
                        diff_is_indel = False
                        if d["type"] in ["del", "ins"]:
                            diff_is_indel = True
                        pdiff = False
                        # compare the diff to each snp for that this mip covers
                        for s in snp_dict:
                            # first check for paralogus differences
                            pdiff = False
                            try:
                                if snp_dict[s]["psv"]:
                                    snp_copy = snp_dict[s]["copies"][ref_copy]
                                    if snp_copy["begin"] == diff_begin or                                          snp_copy["end"] == diff_end:
                                        pdiff = True
                                        break
                            except KeyError as e:
                                continue
                        # update diff dictionary
                        d["clinical"] = False
                        d["clinical_id"] = "none"
                        d["base_match"] = False
                        d["psv"] = pdiff

                    # create a dictionary that holds all the alignment information
                    # for the mip and haplotype
                    al_dict = {"gene_name": gene_name,
                                   "mip_name": m,
                                   "haplotype_ID": hap_name,
                                   "score": score,
                                   "exact_match": exact_match,
                                   "differences": differences,
                                   "aligned_hap": hap_al,
                                   "aligned_ref": ref_al,
                                   "diff": diff}
                    # also report alignments that had any problems
                    if problem_al:
                        problem_alignments.append(al_dict)
                    try:
                        alignments[hap_name][ref_copy].append(al_dict)
                    except KeyError as e:
                        try:
                            alignments[hap_name][ref_copy] = [al_dict]
                        except KeyError as e:
                            alignments[hap_name] = {ref_copy: [al_dict]}
    # pick top alignment by the score for rare cases where a capture sequence
    # can align to a reference in multiple ways
    cleaned_alignments = {}
    for h in alignments:
        for c in alignments[h]:
            copy_als = alignments[h][c]
            if len(copy_als) == 1:
                best_al = copy_als[0]
            else:
                best_al_score = 0
                for i in range(len(copy_als)):
                    sc = copy_als[i]["score"]
                    if sc > best_al_score:
                        best_al_score = sc
                        best_al = copy_als[i]
            try:
                cleaned_alignments[h][c] = best_al
            except KeyError as e:
                cleaned_alignments[h] = {c: best_al}
    alignments = cleaned_alignments
    # check if problem alignments and inverted alignments have a better alignment in
    # alignment dictionary. Remove from list if better is found elsewhere.
    problem_dicts = [problem_alignments, inverted_alignments]
    for i in range(len(problem_dicts)):
        probs = problem_dicts[i]
        for j in range(len(probs)):
            a = probs[j]
            hap_name = a["haplotype_ID"]
            al_score = a["score"]
            try:
                for copyname in alignments[hap_name]:
                    other_score = alignments[hap_name][copyname]["score"]
                    if other_score > al_score:
                        # replace alignment in the list with string "remove"
                        probs[j] = "remove"
                        break
            except KeyError as e:
                continue
        # replace the problem dictionary with the updated version
        temp_dict = {}
        for a in probs:
            if a != "remove":
                hap_name = a["haplotype_ID"]
                try:
                    temp_dict[hap_name].append(a)
                except KeyError as e:
                    temp_dict[hap_name] = [a]
        problem_dicts[i] = temp_dict
        if len(temp_dict) > 0:
            print(("%d alignments may have problems, please check %s"
                    %(len(temp_dict),
                      settings["tempAlignmentsFile"])
                  ))
    if len(problem_snps) > 0:
        print(("%d SNPs may have problems, please check please check %s"
            %(len(problem_snps),
              settings["tempAlignmentsFile"])
              ))

    result =  {"alignments": alignments,
            "inverted_alignments": problem_dicts[1],
            "problem_alignments": problem_dicts[0],
            "problem_snps": problem_snps}
    alignment_file = wdir + settings["tempAlignmentsFile"]
    with open(alignment_file, "w") as outfile:
        json.dump(result, outfile)
    return
def update_aligned_haplotypes (settings):
    """
    Update haplotypes with information from the alignment results.
    Find which paralog copy the haplotype best maps to by scoring each alignment
    according to given psv_priority parameter. If 0, psv's will not be taken into account
    and best alignment score will be considered the correct copy. If 1, alignment score
    will be increased by number of psv's supporting the copy multiplied by the psv_multiplier
    parameter. If 2, number of supporting psv's will be considered for best mapping copy and
    the alignment score will only be used as tie breaker.
    """
    psv_priority = int(settings["psvPriority"])
    psv_multiplier = int(settings["psvMultiplier"])
    wdir = settings["workingDir"]
    temp_haplotype_file = wdir + settings["tempHaplotypesFile"]
    with open(temp_haplotype_file) as infile:
        haplotypes = json.load(infile)
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    temp_alignment_file = wdir + settings["tempAlignmentsFile"]
    with open(temp_alignment_file) as infile:
        parsed_alignments = json.load(infile)
    # all alignments from the parsed alignments dict
    alignments = parsed_alignments["alignments"]
    # alignments with problems
    inverted_alignments = parsed_alignments["inverted_alignments"]
    problem_alignments = parsed_alignments["problem_alignments"]
    # update each haplotype with alignment information
    for m in haplotypes:
        gene_name = m.split("_")[0]
        # find how many paralog specific variant the mip covers
        try:
            mip_psv_count = call_info[gene_name][m]["psv_count"]
        except KeyError as e:
            mip_psv_count = 0
        for h in haplotypes[m]:
            # create a copy dict for each haplotype for each possible
            # paralog gene copy that haplotype may belong to
            copies = haplotypes[m][h]["copies"] = {}
            # get the alignment for this haplotype from alignment dict
            try:
                align = alignments[h]
            except KeyError as e:
                haplotypes[m][h]["mapped"] = False
                continue
            # update copies dict with alignment information
            for c in align:
                # count how many psv's were observed for this copy
                copy_psv_count = 0
                # alignment score shows how well the haplotype was aligned to
                # reference sequence of this copy
                score = align[c]["score"]
                # all differences between the reference and hap is kept in "differences"
                differences = align[c]["differences"]
                for d in differences:
                    if d["psv"]:
                        copy_psv_count += 1
                # update alignment with the psv count
                align[c]["psv_count"] = copy_psv_count
                # update haplotype with alignment information
                haplotypes[m][h]["copies"][c] = {"score": align[c]["score"],
                                              "copy_psv_count": align[c]["psv_count"],
                                              "mip_psv_count": mip_psv_count}
            # Now we know how well this haplotype maps to each copy
            # there are two sources of information used for mapping a haptolype
            # to a copy. Alignment score gives us how well both sequences align
            # It is around 18000 for perfectly aligning human haplotypes (~250 bp)
            # deletions are more costly than substitutions. A 20 bp deletion will
            # cost several thousands of alignment score. The second source is the
            # psv counts. Copy psv count shows hom many psvs were NOT as expected,
            # i.e., points to a different copy than the copy at hand. mip psv count
            # is the total number of psvs covered by mip. The difference between them
            # shows how many psvs support the haplotype mapping to the copy at hand.
            # perfect mapping would show "0" copy psvs and many mip psvs.
            ##########################################################################
            # sort copies considering alignment scores and psv scores
            copy_sort = []
            for c in copies:
                mip_psv_count = copies[c]["mip_psv_count"]
                copy_psv_count = copies[c]["copy_psv_count"]
                supporting_psv_count = mip_psv_count - copy_psv_count
                copy_alignment_score = copies[c]["score"]
                copy_psv_score = psv_multiplier * supporting_psv_count
                copy_sort.append([c,
                                  mip_psv_count,
                                  copy_psv_count,
                                  supporting_psv_count,
                                  copy_alignment_score,
                                  copy_psv_score + copy_alignment_score])
            # sort copies
            # priority indexes (f)irst and (s)econd will determine which score should be
            # considered for sorting the copies for mapping purposes. Indexes are as follows:
            # 0: copyname, 1: mip_pst_count, 2: copy_psv_count, 3: supporting_psv_count
            # 4: alignment_score, 5: copy_psv_score + copy_alignment_score
            # psv priority 0: only consider alignment scores
            # 1: supplement alignment scores with psv scores
            # 2: use psv score, break ties with alignment scores
            try:
                f,s  = [[4,1], [5,1], [3,4]][psv_priority]
            except IndexError as e:
                continue
            # sort copies for the first and second priority keys
            copy_keys_sorted = sorted(copy_sort, key = itemgetter(f, s))
            copy_keys_sorted_temp = copy.deepcopy(copy_keys_sorted)
            # remove alt contigs
            for cop_ind in range(len(copy_keys_sorted)):
                cop = copy_keys_sorted[cop_ind]
                if "alt" in call_info[gene_name][m]["copies"][cop[0]]["chrom"]:
                    copy_keys_sorted[cop_ind] = "remove"
            copy_keys_sorted = [cop_key for cop_key in copy_keys_sorted                                if cop_key != "remove"]
            if len(copy_keys_sorted) == 0:
                copy_keys_sorted = copy_keys_sorted_temp
            # pick best scoring copies
            # last item in copy_keys_sorted is the best
            best_copy = copy_keys_sorted[-1]
            # create a list of copies that has the best score
            best_copies = [cop for cop in copy_keys_sorted                                if (cop [f] == best_copy[f]) and (cop[s] == best_copy[s])]
            # create a map dict to be added to haplotype information
            # extract copy keys of best copies
            temp_keys = [i[0] for i in best_copies]
            temp_dic = {}
            for c in temp_keys:
                temp_dic[c] = {"copy_name": call_info[gene_name][m]["copies"][c]["copyname"],
                               "differences": alignments[h][c]["differences"]}
            haplotypes[m][h]["mapped_copies"] = temp_dic
            # create a single copy name for the haplotype such as HBA1_C0
            # if there are multiple copies that the haplotype matchedequally well
            # name will be the compound name of all best mapping copies, e.g. HBA1_C0_HBA2_C1
            mapped_copy_names = []
            for k in sorted(temp_dic.keys()):
                mapped_copy_names.append("_".join([temp_dic[k]["copy_name"], k]))
            haplotypes[m][h]["copy_name"] = "_".join(mapped_copy_names)
    temp_mapped_haps_file = wdir + settings["tempMappedHaplotypesFile"]
    with open(temp_mapped_haps_file, "w") as outfile:
        json.dump(haplotypes, outfile, indent=1)
    return
def update_unique_haplotypes (settings):
    """
    Add new on and off target haplotypes to the unique haplotypes.
    Update sequence_to_haplotype dict with the new haplotypes.
    """
    wdir = settings["workingDir"]
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    sequence_to_haplotype_file = wdir + settings["sequenceToHaplotypeDictionary"]
    try:
        with open(unique_haplotype_file) as infile:
            unique_haplotypes = json.load(infile)
    except IOError as e:
        unique_haplotypes = {}
    try:
        with open(sequence_to_haplotype_file) as infile:
            sequence_to_haplotype = json.load(infile)
    except IOError as e:
        sequence_to_haplotype = {}
    temp_mapped_hap_file = wdir + settings["tempMappedHaplotypesFile"]
    with open(temp_mapped_hap_file) as infile:
        haplotypes = json.load(infile)
    temp_off_file = wdir + settings["tempOffTargetsFile"]
    with open(temp_off_file) as infile:
        off_target_haplotypes = json.load(infile)
    # update unique_haplotypes with on target haplotypes
    for m in haplotypes:
        for h in haplotypes[m]:
            uniq_id = h + "-0"
            try:
                if uniq_id in unique_haplotypes[m]:
                    counter = 0
                    while uniq_id in unique_haplotypes[m]:
                        counter += 1
                        uniq_id = h + "-" + str(counter)
                    unique_haplotypes[m][uniq_id] = haplotypes[m][h]
                else:
                    unique_haplotypes[m][uniq_id] = haplotypes[m][h]
            except KeyError as e:
                unique_haplotypes[m] = {uniq_id: haplotypes[m][h]}

    # update unique_haplotypes with off target haplotypes
    for h in off_target_haplotypes:
        m = h.split(".")[0]
        uniq_id = h + "-0"
        try:
            if uniq_id in unique_haplotypes[m]:
                counter = 0
                while uniq_id in unique_haplotypes[m]:
                    counter += 1
                    uniq_id = h + "-" + str(counter)
                unique_haplotypes[m][uniq_id] = off_target_haplotypes[h]
            else:
                unique_haplotypes[m][uniq_id] = off_target_haplotypes[h]
        except KeyError as e:
            unique_haplotypes[m] = {uniq_id: off_target_haplotypes[h]}
    # update sequence_to_haplotype with new haplotypes
    for u in unique_haplotypes:
        for h in unique_haplotypes[u]:
            if not unique_haplotypes[u][h]["sequence"] in sequence_to_haplotype:
                sequence_to_haplotype[unique_haplotypes[u][h]["sequence"]] = h
    with open(unique_haplotype_file, "w") as outfile:
        json.dump(unique_haplotypes, outfile, indent=1)
    with open(sequence_to_haplotype_file, "w") as outfile:
        json.dump(sequence_to_haplotype, outfile, indent=1)
    return
def update_variation(settings):
    """
    Add new on and off target haplotypes to the unique haplotypes.
    Update sequence_to_haplotype dict with the new haplotypes.
    """
    wdir = settings["workingDir"]
    species = settings["species"]
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    variation_file = wdir + settings["variationDictionary"]
    var_key_to_uniq_file = wdir + settings["variationKeyToUniqueKey"]
    with open(unique_haplotype_file) as infile:
        haplotypes = json.load(infile)
    try:
        with open(variation_file) as infile:
            variation = json.load(infile)
    except IOError as e:
        variation = {}
    var_key_to_uniq_file = wdir +  settings["variationKeyToUniqueKey"]
    try:
        with open(var_key_to_uniq_file) as infile,             open(var_key_to_uniq_file + id_generator(6), "w") as outfile:
            var_key_to_uniq = json.load(infile)
    except IOError as e:
        var_key_to_uniq = {}
    #
    outfile_list = ["##fileformat=VCFv4.1"]
    outfile_list.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT"]))
    temp_variations = []
    for m in haplotypes:
        for h in haplotypes[m]:
            if haplotypes[m][h]["mapped"]:
                try:
                    haplotypes[m][h]["left_normalized"]
                except KeyError as e:
                    left_normalized = True
                    for c in haplotypes[m][h]["mapped_copies"]:
                        differences = haplotypes[m][h]["mapped_copies"][c]["differences"]
                        for d in differences:
                            var_key = d["vcf_raw"]
                            try:
                                uniq_var_key = var_key_to_uniq[var_key]
                                d["annotation"] = variation[uniq_var_key]
                                d["vcf_normalized"] = uniq_var_key
                            except KeyError as e:
                                left_normalized = False
                                temp_variations.append(var_key)
                    if left_normalized:
                         haplotypes[m][h]["left_normalized"] = True
    temp_variations = [temp_var.split(":") for temp_var in set(temp_variations)]
    temp_variations = [[v[0], int(v[1])] + v[2:] for v in temp_variations]
    temp_variations = sorted(temp_variations, key = itemgetter(0, 1))
    temp_variations_lines = ["\t".join(map(str, v)) for v in temp_variations]
    temp_variation_keys = [":".join(map(str, v)) for v in temp_variations]
    outfile_list.extend(temp_variations_lines)
    raw_vcf_file = settings["rawVcfFile"]
    zipped_vcf = raw_vcf_file + ".gz"
    norm_vcf_file = settings["normalizedVcfFile"]
    with open(wdir + raw_vcf_file, "w") as outfile:
        outfile.write("\n".join(outfile_list))
    with open(wdir + zipped_vcf, "w") as outfile:
        dump = subprocess.call(["bgzip", "-c", "-f", raw_vcf_file],
                               cwd = wdir, stdout=outfile)
    dump = subprocess.call(["bcftools", "index", "-f", raw_vcf_file + ".gz"], cwd = wdir)
    unmasked_genome = get_file_locations()[species]["unmasked_fasta_genome"]
    dump = subprocess.call(["bcftools", "norm", "-f", unmasked_genome,
                            "-cw", "-w", "0",
                           "-o", norm_vcf_file, raw_vcf_file + ".gz"], cwd = wdir)
    ann_db_dir = get_file_locations()[species]["annotation_db_dir"]
    ann_build = settings["annotationBuildVersion"]
    ann_protocol = settings["annotationProtocol"].replace(";", ",")
    ann_operation = settings["annotationOperation"].replace(";", ",")
    ann_nastring = settings["annotationNaString"]
    ann_out = settings["annotationOutput"]
    try:
        ann_script = settings["annotationScript"]
    except KeyError as e:
        ann_script = "table_annovar.pl"
    ann_command = [ann_script,
                    norm_vcf_file,
                    ann_db_dir,
                    "-buildver", ann_build,
                    "-vcfinput",
                    "-protocol",  ann_protocol,
                    "-operation", ann_operation,
                    "-nastring", ann_nastring,
                    "-out", ann_out]
    #print " ".join(ann_command)
    dump = subprocess.check_call(ann_command,
                           cwd=wdir)
    normalized_variation_keys = []
    with open(wdir + norm_vcf_file) as infile:
        line_num = 0
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                var_key = temp_variation_keys[line_num]
                normalized_key = ":".join(newline[:5])
                var_key_to_uniq[var_key] = normalized_key
                normalized_variation_keys.append(normalized_key)
                line_num += 1
    #
    annotation_table_file = ann_out + "." + ann_build + "_multianno.txt"
    with open(wdir + annotation_table_file) as infile:
        line_num = -1
        for line in infile:
            newline = line.strip().split("\t")
            if line_num == -1:
                line_num += 1
                colnames = newline
            else:
                normalized_key = normalized_variation_keys[line_num]
                if not normalized_key in variation:
                    variation[normalized_key] = {colnames[i] : newline[i]                                                  for i in range(len(colnames))}
                line_num += 1
    #print temp_variation_keys[:10]
    #print var_key_to_uniq.keys()[:10]
    if line_num != len(temp_variation_keys):
        print("There are more variation keys then annotated variants.")

    for m in haplotypes:
        for h in haplotypes[m]:
            if haplotypes[m][h]["mapped"]:
                try:
                    haplotypes[m][h]["left_normalized"]
                except KeyError as e:
                    for c in haplotypes[m][h]["mapped_copies"]:
                        differences = haplotypes[m][h]["mapped_copies"][c]["differences"]
                        for d in differences:
                            var_key = d["vcf_raw"]
                            uniq_var_key = var_key_to_uniq[var_key]
                            d["annotation"] = variation[uniq_var_key]
                            d["vcf_normalized"] = uniq_var_key
                    haplotypes[m][h]["left_normalized"] = True

    with open(unique_haplotype_file, "w") as outfile:
        json.dump(haplotypes, outfile)
    with open(variation_file, "w") as outfile:
        json.dump(variation, outfile)
    with open(var_key_to_uniq_file, "w") as outfile:
        json.dump(var_key_to_uniq, outfile)
    try:
        if int(settings["mergeSNPs"]):
            dump = merge_snps(settings)
    except KeyError:
        pass
    return
def get_raw_data (settings):
    """ Extract raw data from filtered_data file. If there is data from a previous run,
    new data can be added to old data. If this sample set is new, or being analyzed
    separately, than existing data should be "na".
    Return a list of data point dictionaries. One dict for each haplotype/sample.
    Write this list to disk."""
    wdir = settings["workingDir"]
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    sequence_to_haplotype_file = wdir + settings["sequenceToHaplotypeDictionary"]
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    with open(sequence_to_haplotype_file) as infile:
        sequence_to_haplotype = json.load(infile)
    problem_data = []
    existing_data_file = wdir + settings["existingData"]
    try:
        with open(existing_data_file) as infile:
            raw_data = json.load(infile)
    except IOError:
        raw_data = []
    mipster_file = wdir + settings["mipsterFile"]
    colnames = dict(list(zip(settings["colNames"],
                        settings["givenNames"])))
    with open(mipster_file) as infile:
        ## filteredData only contains relevant fields for this analysis
        # the field names are given in the settings dict and
        # kept in given_names list
        run_ID = settings["runID"]
        line_number = 0
        for line in infile:
            newline = line.strip().split("\t")
            if not line_number == 0:
                line_number += 1
                col_indexes = {}
                try:
                    for ck in list(colnames.keys()):
                        col_indexes[
                            newline.index(ck)
                        ] = {"name": colnames[ck]}
                        if ck == 'c_barcodeCnt':
                            bc_index = newline.index(ck)
                        elif ck == 'c_readCnt':
                            rc_index = newline.index(ck)
                except ValueError:
                    for ck in list(colnames.values()):
                        col_indexes[
                        newline.index(ck)
                        ] = {"name" : ck}
                        if ck == 'barcode_count':
                            bc_index = newline.index(ck)
                        elif ck == 'read_count':
                            rc_index = newline.index(ck)
                # each data point will be a dict
                data_dic = {}
                for i in list(col_indexes.keys()):
                    # check data type and convert data appropriately
                    if i in [bc_index, rc_index]:
                        data_point = int(newline[i])
                    else:
                        data_point = newline[i]
                    data_dic[col_indexes[i]["name"]] = data_point
                data_dic["run_ID"] = run_ID
                # once data dict is created, check the haplotype sequence
                # in the sequence to haplotype dict and update the data dict
                # with values in the haplotype dict
                seq = data_dic.pop("haplotype_sequence")
                try:
                    uniq_id = sequence_to_haplotype[seq]
                    data_dic.pop("haplotype_quality_scores")
                    data_dic["haplotype_ID"] = uniq_id
                    if unique_haplotypes[data_dic["mip_name"]][uniq_id]["mapped"]:
                        raw_data.append(data_dic)
                except KeyError:
                    problem_data.append(data_dic)
                    continue
    # dump the raw_data list to a json file
    raw_data_file = wdir + settings["rawDataFile"]
    with open(raw_data_file, "w") as outfile:
        json.dump(raw_data, outfile)
    problem_data_file = wdir + settings["rawProblemData"]
    with open(problem_data_file, "w") as outfile:
        json.dump(problem_data, outfile, indent=1)
    return
def filter_data (data_file, filter_field, comparison, criteria, output_file):
    """
    Data fields to filter: sample_name, gene_name, mip_name, barcode_count etc.
    Comparison: for numeric values: gt (greater than), lt, gte (greater than or equal to), lte
    criteria must be numeric as well.
    for string values: in, nin (not in), criteria must be a list to include or exclude.
    for string equality assessments, a list with single value can be used.
    """
    filtered_data = []
    with open(data_file) as infile:
        data = json.load(infile)
    for d in data:
        field_value = d[filter_field]
        if comparison == "gt":
            if field_value > criteria:
                filtered_data.append(d)
        elif comparison == "gte":
            if field_value >= criteria:
                filtered_data.append(d)
        elif comparison == "lt":
            if field_value < criteria:
                filtered_data.append(d)
        elif comparison == "lte":
            if field_value <= criteria:
                filtered_data.append(d)
        elif comparison == "in":
            if field_value in criteria:
                filtered_data.append(d)
        elif comparison == "nin":
            if field_value not in criteria:
                filtered_data.append(d)
    with open(output_file, "w") as outfile:
        json.dump(filtered_data, outfile, indent = 1)
def group_samples (settings):
    wdir = settings["workingDir"]
    raw_data_file = wdir  + settings["rawDataFile"]
    with open(raw_data_file) as infile:
        raw_data = json.load(infile)
    samples = {}
    problem_samples = {}
    mip_names = {}
    #diploid_mipnames = []
    #all_mipnames = []
    copy_stable_genes = settings["copyStableGenes"]
    for c in raw_data:
        gene_name = c["gene_name"]
        sample_name = c["sample_name"]
        mip_name = c["mip_name"]
        if (copy_stable_genes == "na") or (gene_name in copy_stable_genes):
            try:
                samples[sample_name]["diploid_barcode_count"] += c["barcode_count"]
                samples[sample_name]["diploid_read_count"] += c["read_count"]
                samples[sample_name]["diploid_mip_names"].append(mip_name)
                #diploid_mipnames.append(mip_name)
                """
                except KeyError:
                    try:
                        samples[sample_name]["diploid_barcode_count"] = c["barcode_count"]
                        samples[sample_name]["diploid_read_count"] = c["read_count"]
                        samples[sample_name]["diploid_mip_names"] = [mip_name]
                        diploid_mipnames.append(mip_name)
                """
            except KeyError:
                samples[sample_name] = {"diploid_barcode_count": c["barcode_count"],
                                        "diploid_read_count": c["read_count"],
                                        "diploid_mip_names": [mip_name],
                                        "total_barcode_count": c["barcode_count"],
                                        "total_read_count": c["read_count"],
                                        "all_mip_names": [mip_name]
                                        }
                #diploid_mipnames.append(mip_name)
        try:
            samples[sample_name]["total_barcode_count"] += c["barcode_count"]
            samples[sample_name]["total_read_count"] += c["read_count"]
            samples[sample_name]["all_mip_names"].append(mip_name)
            #all_mipnames.append(mip_name)
            """
            except KeyError:
                try:
                    samples[sample_name]["total_barcode_count"] = c["barcode_count"]
                    samples[sample_name]["total_read_count"] = c["read_count"]
                    samples[sample_name]["all_mip_names"] = [mip_name]
                    all_mipnames.append(mip_name)
            """
        except KeyError:
            samples[sample_name] = {"diploid_barcode_count": 0,
                                        "diploid_read_count": 0,
                                        "diploid_mip_names": [],
                                        "total_barcode_count": c["barcode_count"],
                                        "total_read_count": c["read_count"],
                                        "all_mip_names": [mip_name]
                                        }
            #all_mipnames.append(mip_name)


    for s in list(samples.keys()):
        try:
            samples[s]["diploid_mip_names"] = list(set((samples[s]["diploid_mip_names"])))
            samples[s]["diploid_mip_number"] = len(samples[s]["diploid_mip_names"])
            samples[s]["all_mip_names"] = list(set((samples[s]["all_mip_names"])))
            samples[s]["total_mip_number"] = len(samples[s]["all_mip_names"])
            mip_names[s] = {}
            mip_names[s]["diploid_mip_names"] = samples[s].pop("diploid_mip_names")
            mip_names[s]["all_mip_names"] = samples[s].pop("all_mip_names")
            try:
                samples[s]["average_diploid_barcode_count"] = round(
                samples[s]['diploid_barcode_count']
                 /samples[s]['diploid_mip_number'], 2)
                samples[s]["sample_normalizer"] = round(
                100/(samples[s]["average_diploid_barcode_count"]), 2)
            except ZeroDivisionError:
                problem_samples[s] = samples.pop(s)
        except KeyError:
            problem_samples[s] = samples.pop(s)

    results =  {"samples": samples,
            "problem_samples": problem_samples,
            "mip_names": mip_names}
    sample_info_file = wdir + settings["sampleInfoFile"]
    with open(sample_info_file, "w") as outfile:
        json.dump(results,outfile, indent=1)
    return
def update_raw_data (settings):
    wdir = settings["workingDir"]
    raw_data_file = wdir  + settings["rawDataFile"]
    with open(raw_data_file) as infile:
        raw_data = json.load(infile)
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    sample_info_file = wdir + settings["sampleInfoFile"]
    with open(sample_info_file) as infile:
        samples = json.load(infile)["samples"]
    normalized_data = []
    problem_data = []
    for r in raw_data:
        try:
            sample_name = r["sample_name"]
            r["sample_normalizer"] = samples[sample_name]["sample_normalizer"]
            r["sample_normalized_barcode_count"] = r["sample_normalizer"] * r["barcode_count"]
            hid = r["haplotype_ID"]
            mid = r["mip_name"]
            r["copy_name"] = unique_haplotypes[mid][hid]["copy_name"]
            normalized_data.append(r)
        except KeyError:
            problem_data.append(r)
    normalized_data_file = wdir + settings["normalizedDataFile"]
    with open(normalized_data_file, "w") as outfile:
        json.dump(normalized_data, outfile)
    problem_data_file = wdir + settings["normalizedProblemData"]
    with open(problem_data_file, "w") as outfile:
        json.dump(problem_data, outfile, indent=1)
    return
def get_counts (settings):
    wdir = settings["workingDir"]
    data_file = wdir + settings["normalizedDataFile"]
    with open(data_file) as infile:
        data = json.load(infile)
    min_barcode_count = int(settings["minBarcodeCount"])
    min_barcode_fraction = float(settings["minBarcodeFraction"])
    counts = {}
    samples = {}
    # get sample data across probes
    for t in data:
        g = t["gene_name"]
        m = t["mip_name"]
        c = t["copy_name"]
        s = t["sample_name"]
        try:
            samples[s][g][m][c]["raw_data"].append(t)
        except KeyError:
            try:
                samples[s][g][m][c] = {"raw_data": [t]}
            except KeyError:
                try:
                    samples[s][g][m] = {c: {"raw_data": [t]}}
                except KeyError:
                    try:
                        samples[s][g] = {m: {c: {"raw_data": [t]}}}
                    except KeyError:
                        samples[s] = {g: {m: {c: {"raw_data": [t]}}}}
    # filter/merge individual data points for each sample/paralog copy
    # data for the same haplotype, possibly from separate runs will be merged
    # data that does not pass a low barcode threshold will be filtered
    merge_keys = ["sample_normalized_barcode_count",
                  "barcode_count",
                  "read_count"]
    for s in samples:
        for g in samples[s]:
            for m in samples[s][g]:
                for c in samples[s][g][m]:
                    data_points = samples[s][g][m][c]["raw_data"]
                    grouped_data = samples[s][g][m][c]["grouped_data"] = []
                    merged_data = samples[s][g][m][c]["merged_data"] = []
                    filtered_data = samples[s][g][m][c]["filtered_data"] = []
                    cumulative_data = samples[s][g][m][c]["cumulative_data"] = {                            "sample_normalized_barcode_count": 0,
                            "barcode_count": 0,
                            "read_count": 0,
                            "haplotypes": []}
                    temp_data_points = copy.deepcopy(data_points)
                    # group data points with same haplotype_ID together
                    for i in range(len(temp_data_points)):
                        di = temp_data_points[i]
                        if di != "remove":
                            same_haps = [di]
                            for j in range(len(temp_data_points)):
                                if i != j:
                                    dj = temp_data_points[j]
                                    if dj != "remove" and (di["haplotype_ID"] ==                                                            dj["haplotype_ID"]):
                                        same_haps.append(dj)
                                        temp_data_points[j] = "remove"
                            grouped_data.append(same_haps)
                    for dl in grouped_data:
                        # find data point with most barcodes
                        temp_barcode_count = 0
                        best_data_index = 0
                        for i in range(len(dl)):
                            b = dl[i]["barcode_count"]
                            if b > temp_barcode_count:
                                temp_barcode_count = b
                                best_data_index = i
                        # use data point with more barcodes as base
                        dm = copy.deepcopy(dl[best_data_index])
                        for i in range(len(dl)):
                            if i != best_data_index:
                                dt = dl[i]
                                for k in merge_keys:
                                    dm[k] += dt[k]
                        merged_data.append(dm)
                    # filter haplotypes with insufficient count/frequency
                    temp_barcode_counts = []
                    for d in merged_data:
                        temp_barcode_counts.append(d["barcode_count"])
                    frac_threshold = min_barcode_fraction * sum(temp_barcode_counts)
                    for d in merged_data:
                        if d["barcode_count"] >= frac_threshold and                           d["barcode_count"] >= min_barcode_count:
                            filtered_data.append(d)
                    for d in filtered_data:
                        for k in merge_keys:
                            cumulative_data[k] += d[k]
                        cumulative_data["haplotypes"].append(d["haplotype_ID"])
    # get probe information across samples
    template_dict = {"sample_normalized_barcode_count": [],
                     "barcode_count": [],
                     "read_count" : [],
                     "haplotypes": [],
                     "sample_names": []}
    for s in samples:
        for g in samples[s]:
            for m in samples[s][g]:
                for c in samples[s][g][m]:
                    filtered_data = samples[s][g][m][c]["filtered_data"]
                    cumulative_data = samples[s][g][m][c]["cumulative_data"]
                    try:
                        norm_counts = counts[g][m][c]["sample_normalized_barcode_count"]
                    except KeyError:
                        try:
                            counts[g][m][c] = copy.deepcopy(template_dict)
                        except KeyError:
                            try:
                                counts[g][m] = {c: copy.deepcopy(template_dict)}
                            except KeyError:
                                counts[g] = {m: {c: copy.deepcopy(template_dict)}}
                    for k in merge_keys:
                        counts[g][m][c][k].append(cumulative_data[k])
                    counts[g][m][c]["sample_names"].append(s)
                    counts[g][m][c]["haplotypes"].append(filtered_data)
    sample_results_file = wdir + settings["perSampleResults"]
    with open(sample_results_file, "w") as outfile:
        json.dump(samples, outfile, indent=1)
    probe_results_file = wdir + settings["perProbeResults"]
    with open(probe_results_file, "w") as outfile:
        json.dump(counts, outfile, indent=1)
    return
def get_unique_probes(settings):
    wdir = settings["workingDir"]
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    # keep which copies a mip can be mapped to in copy_names dictionary
    copy_names = {}
    problem_copy_names = []
    gene_names = []
    for m in unique_haplotypes:
        g = m.split("_")[0]
        gene_names.append(g)
        for h in unique_haplotypes[m]:
            if unique_haplotypes[m][h]["mapped"]:
                try:
                    cn = unique_haplotypes[m][h]["copy_name"]
                except KeyError:
                    problem_copy_names.append(h)
                try:
                    copy_names[g][m].append(cn)
                except KeyError:
                    try:
                        copy_names[g][m] = [cn]
                    except KeyError:
                        copy_names[g] = {m: [cn]}
    for g in copy_names:
        for m in copy_names[g]:
            copy_names[g][m] = sorted(list(set(copy_names[g][m])))
    # keep all copy names that have been observed for a gene in all_copies dict
    all_copy_names = {}
    all_copy_names_list = []
    for g in copy_names:
        cn = []
        for m in copy_names[g]:
            cn.extend(copy_names[g][m])
        cn = sorted(list(set(cn)))
        all_copy_names[g] = cn
        all_copy_names_list.extend(cn)

    # keep probe names that always maps to specific paralog copies in unique_keys
    # uniq_keys[g][m1] = ["HBA1_C0", "HBA2_C1"]
    uniq_keys = {}
    # keep probe names that never maps to specific paralog copies in non_uniq_keys
    # non_uniq_keys[g][m2] = ["HBA1_C0_HBA2_C1"]
    non_uniq_keys = {}
    # keep probe names that maps to specific paralog copies for some copies only
    # in semi_uniq_keys [g][m3] = ["HBA1_C0"] (only maps C0 specifically)
    #semi_uniq_keys = {}
    # if a mip does not map to any target copy (always off target or never works)
    no_copy_keys = {}
    for g in all_copy_names:
        copies = all_copy_names[g]
        keys = sorted(list(call_info[g].keys()),
                      key= lambda a: call_info[g][a]["copies"]["C0"]["capture_start"])
        for k in keys:
            try:
                key_copies = copy_names[g][k]
                nucopies = []
                toss_copies = []
                try:
                    for c in key_copies:
                        split_copies = c.split("_")
                        nucopies.append([split_copies[i+1]                                              for i in range(0, len(split_copies), 2)])
                except IndexError:
                    print(split_copies)
                    break
                for i in range(len(nucopies)):
                    query_list = nucopies[i]
                    for j in range(len(nucopies)):
                        target_list = nucopies[j]
                        if i != j:
                            for c in query_list:
                                if c in target_list:
                                    toss_copies.append(i)
                toss_copies = list(set(toss_copies))
                uniq_copies = [key_copies[i]                                for i in range(len(key_copies)) if i not in toss_copies]
                non_uniq_copies = [key_copies[i]                                for i in range(len(key_copies)) if i in toss_copies]
                if len(uniq_copies) > 0:
                    try:
                        uniq_keys[g][k] = uniq_copies
                    except KeyError:
                        uniq_keys[g] = {k: uniq_copies}
                if len(non_uniq_copies) > 0:
                    try:
                        non_uniq_keys[g][k] = non_uniq_copies
                    except KeyError:
                        non_uniq_keys[g] = {k: non_uniq_copies}
            except KeyError:
                try:
                    no_copy_keys[g].append(k)
                except KeyError:
                    no_copy_keys[g] = [k]
    gene_names = sorted(gene_names)
    result =  {"copy_names": copy_names,
            "all_copy_names": all_copy_names,
            "problem_copy_names": problem_copy_names,
            "gene_names": gene_names,
            "unique_probes": uniq_keys,
            "non_unique_probes": non_uniq_keys,
            "no_copy_probes": no_copy_keys}
    uniq_file = wdir + settings["uniqueProbeFile"]
    with open(uniq_file, "w") as outfile:
        json.dump(result, outfile, indent = 1)
    return
def create_data_table (settings):
    wdir = settings["workingDir"]
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    uniq_file = wdir + settings["uniqueProbeFile"]
    with open(uniq_file) as infile:
        uniq_dict = json.load(infile)
    sample_results_file = wdir + settings["perSampleResults"]
    with open(sample_results_file) as infile:
        counts = json.load(infile)
    sample_names = list(counts.keys())
    big_table = []
    all_tables = {}
    for gene in call_info:
        gene_info = call_info[gene]
        try:
            all_probes = uniq_dict["unique_probes"][gene]
        except KeyError:
            all_probes = {}
        table = []
        for p in all_probes:
            for c in all_probes[p]:
                split_copies = c.split("_")
                if len(split_copies) == 2:
                    c_id = c.split("_")[-1]
                    probe_info = gene_info[p]["copies"][c_id]
                    try:
                        chrom = int(probe_info["chrom"][3:])
                    except ValueError:
                        chrom = 23
                    start = probe_info["capture_start"]
                    end = probe_info["capture_end"]
                    ori = probe_info["orientation"]
                    gc_frac = calculate_gc(probe_info["capture_sequence"])
                elif len(split_copies) > 2:
                    gc_list = []
                    for i in range(0, len(split_copies), 2):
                        c_id = split_copies[i+1]
                        probe_info = gene_info[p]["copies"][c_id]
                        gc_list.append(calculate_gc(probe_info["capture_sequence"]))
                    gc_frac = np.mean(gc_list)
                    start = "na"
                    end = "na"
                    ori = "na"
                    chrom = 99
                s_counts = []
                for s in sample_names:
                    try:
                        bc = counts[s][gene][p][c]["cumulative_data"]                             ["sample_normalized_barcode_count"]
                    except KeyError:
                        bc = 0

                    s_counts.append(bc)
                p_counts = [gene, p, c, chrom, start, end, gc_frac, ori]
                p_counts.extend(s_counts)
                table.append(p_counts)
                big_table.append(p_counts)
        table = sorted(table, key = itemgetter(3, 4, 5))
        all_tables[gene] = table
    big_table = sorted(big_table, key = itemgetter(3, 4, 5))
    result = {"tables": all_tables,
            "big_table" : big_table,
            "sample_names": sample_names}
    tables_file = wdir + settings["tablesFile"]
    with open(tables_file, "w") as outfile:
        json.dump(result, outfile)
    return
def filter_tables (settings):
    wdir = settings["workingDir"]
    tables_file = wdir + settings["tablesFile"]
    with open(tables_file) as infile:
        tables_dic = json.load(infile)
    tables = copy.deepcopy(tables_dic["tables"])
    sample_info_file = wdir + settings["sampleInfoFile"]
    with open(sample_info_file) as infile:
        sample_info = json.load(infile)["samples"]
    min_probe_median = int(settings["minimumProbeMedian"])
    min_sample_median = int(settings["minimumSampleMedian"])
    filtered_tables = {}
    #normalized_tables = {}
    #barcode_tables = {}
    sample_names = copy.deepcopy(tables_dic["sample_names"])
    sample_normalizers = np.asarray([sample_info[s]["sample_normalizer"] for s in sample_names])
    sample_array = np.asarray(sample_names)
    #filtered_samples = sample_array[sample_filter]
    for g in tables:
        t_array = np.asarray(tables[g])
        split_array = np.hsplit(t_array, [8])
        t_info = split_array[0]
        t_data = np.array(split_array[1], float)
        #s_filtered = t_data[:, sample_filter]
        barcode_table = np.transpose(t_data)/sample_normalizers[:, np.newaxis]
        try:
            probe_filter = np.median(t_data, axis = 1) >= min_probe_median
        except IndexError:
            continue
        sample_filter = np.median(barcode_table, axis = 1) >= min_sample_median
        info_filtered = t_info[probe_filter, :]
        data_filtered = t_data[:, sample_filter]
        data_filtered = data_filtered[probe_filter, :]
        median_normalized = 2*(data_filtered / np.median(data_filtered, axis = 1)[:,np.newaxis])
        barcode_filtered = np.transpose(barcode_table)[:, sample_filter]
        barcode_filtered = barcode_filtered[probe_filter, :]
        #filtered_data = np.hstack((info_filtered, data_filtered))
        #normalized_data = np.hstack((info_filtered, median_normalized))
        #barcode_data = np.hstack((info_filtered, barcode_filtered))
        filtered_tables[g] = {"probe_information":info_filtered,
                              "barcode_counts": barcode_filtered,
                              "sample_normalized_barcode_counts": data_filtered,
                              "median_normalized_copy_counts": median_normalized,
                              "sample_names": sample_array[sample_filter]}
        #normalized_tables[g] = normalized_data
        #barcode_tables [g] = barcode_data
    filtered_tables_file = wdir + settings["filteredTablesFile"]
    with open(filtered_tables_file, "wb") as outfile:
        pickle.dump(filtered_tables, outfile)
    return
def generate_clusters(settings):
    cluster_output = {}
    problem_clusters = []
    wdir = settings["workingDir"]
    filtered_tables_file = wdir + settings["filteredTablesFile"]
    with open(filtered_tables_file, 'rb') as infile:
        tables = pickle.load(infile)
    case_file = wdir + settings["caseFile"]
    case_status = {}
    try:
        with open(case_file) as infile:
            for line in infile:
                newline = line.strip().split("\t")
                case_status[newline[0]] = newline[1]
    except IOError:
        pass
    for g in tables:
        try:
            probes = list(tables[g]["probe_information"][:, 1])
            uniq_probes = []
            for p in probes:
                if p not in uniq_probes:
                    uniq_probes.append(p)
            sample_names = tables[g]["sample_names"]
            labels = []
            for s in sample_names:
                try:
                    labels.append(case_status[s])
                except KeyError:
                    labels.append("na")
            copy_counts = tables[g]["median_normalized_copy_counts"]
            barcode_counts = np.transpose(tables[g]["sample_normalized_barcode_counts"])
            collapsed_counts = np.zeros((len(uniq_probes), len(sample_names)))
            for i in range(len(uniq_probes)):
                u = uniq_probes[i]
                for j in range(len(probes)):
                    p = probes[j]
                    if u == p:
                        collapsed_counts[i] += copy_counts[j]
            repeat_counts = []
            for u in uniq_probes:
                repeat_counts.append(probes.count(u))
            upper_limit = np.asarray(repeat_counts)*2 + 0.5
            lower_limit = np.asarray(repeat_counts)*2 - 0.5
            copy_counts = np.transpose(copy_counts)
            collapsed_counts = np.transpose(collapsed_counts)
            ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=100)
            #ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=30)
            #Y  = ts.fit_transform(copy_counts)
            Y  = ts.fit_transform(barcode_counts)
            ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=50)
            V  = ts.fit_transform(Y)
            U  = ts.fit_transform(V)
            T  = ts.fit_transform(U)
            tsne_keys = {"Y": Y, "V": V, "U": U, "T": T}
            ms = MeanShift(bandwidth = 5, cluster_all = True)
            ms.fit(tsne_keys[settings["tsneKey"]])
            cluster_labels = ms.labels_
            cluster_centers = ms.cluster_centers_
            labels_unique = np.unique(cluster_labels)
            n_clusters_ = len(labels_unique)
            cluster_dict = {"cluster_labels" : cluster_labels,
                            "unique_probes": uniq_probes,
                           "all_probes": probes,
                           "plotting": {"upper_limit": upper_limit,
                                        "lower_limit": lower_limit},
                           "mean_shift": ms,
                           "tsne": {"median": copy_counts,
                                    "collapsed_median": collapsed_counts,
                                    "T": T,
                                    "U": U,
                                    "V": V,
                                    "Y": Y},
                           "clusters": {}}
            for k in range(n_clusters_):
                my_members = cluster_labels == k
                cluster_center = cluster_centers[k]
                cluster_case_count = list(np.asarray(labels)[my_members]).count("case")
                cluster_control_count = list(np.asarray(labels)[my_members]).count("control")
                other_case_count = labels.count("case") - cluster_case_count
                other_control_count = labels.count("control") - cluster_control_count
                cluster_table = [[cluster_case_count, cluster_control_count],
                                 [other_case_count, other_control_count]]
                try:
                    cluster_fish = fisher_exact(cluster_table)
                except:
                    cluster_fish = ["na", "na"]
                try:
                    cluster_chi = chi2_contingency(cluster_table)
                except:
                    cluster_chi = ["na", "na", "na", "na"]
                cluster_dict["clusters"][k] = {"cluster_table": cluster_table,
                           "cluster_stats": {"fisher": {"OR": cluster_fish[0],
                                                        "pvalue": cluster_fish[1]},
                                             "chi": {"pvalue": cluster_chi[1],
                                                     "expected": cluster_chi[3]}},
                           "cluster_data": copy_counts[my_members],
                           "cluster_medians": np.median(copy_counts[my_members], axis=0),
                           "cluster_samples": sample_names[my_members],
                           "cluster_members": my_members,
                           "collapsed_data": collapsed_counts[my_members],
                           "collapsed_medians": np.median(collapsed_counts[my_members], axis=0)
                           }
            cluster_output[g] = cluster_dict
        except Exception as e:
            problem_clusters.append([g, e])
    cluster_output_file = wdir + settings["clusterOutputFile"]
    with open(cluster_output_file, "wb") as outfile:
        pickle.dump(cluster_output, outfile)
    return problem_clusters
def dbscan_clusters(settings):
    cluster_output = {}
    problem_clusters = []
    wdir = settings["workingDir"]
    cluster_output_file = wdir + settings["clusterOutputFile"]
    with open(cluster_output_file, "rb") as infile:
        cnv_calls = pickle.load(infile)
    filtered_tables_file = wdir + settings["filteredTablesFile"]
    with open(filtered_tables_file, "rb") as infile:
        tables = pickle.load(infile)
    case_file = wdir + settings["caseFile"]
    case_status = {}
    try:
        with open(case_file) as infile:
            for line in infile:
                newline = line.strip().split("\t")
                case_status[newline[0]] = newline[1]
    except IOError:
        pass
    tsne_key = settings["tsneKey"]
    tsne_plot_key = settings["tsnePlotKey"]
    min_samples = int(settings["minClusterSamples"])
    max_unclustered_frac = float(settings["maxUnclusteredFrac"])
    max_cluster_count = int(settings["maxClusterCount"])
    #tsne_keys = {"Y": Y, "V": V, "U": U, "T": T}
    for g in tables:
        try:
            probes = list(tables[g]["probe_information"][:, 1])
            uniq_probes = []
            for p in probes:
                if p not in uniq_probes:
                    uniq_probes.append(p)
            sample_names = tables[g]["sample_names"]
            labels = []
            for s in sample_names:
                try:
                    labels.append(case_status[s])
                except KeyError:
                    labels.append("na")
            copy_counts = tables[g]["median_normalized_copy_counts"]
            collapsed_counts = np.zeros((len(uniq_probes), len(sample_names)))
            for i in range(len(uniq_probes)):
                u = uniq_probes[i]
                for j in range(len(probes)):
                    p = probes[j]
                    if u == p:
                        collapsed_counts[i] += copy_counts[j]
            repeat_counts = []
            for u in uniq_probes:
                repeat_counts.append(probes.count(u))
            upper_limit = np.asarray(repeat_counts)*2 + 0.5
            lower_limit = np.asarray(repeat_counts)*2 - 0.5
            copy_counts = np.transpose(copy_counts)
            collapsed_counts = np.transpose(collapsed_counts)
            Y = cnv_calls[g]["tsne"][tsne_key]
            P = cnv_calls[g]["tsne"][tsne_plot_key]
            cluster_count = []
            eps_range = np.arange(2., 11.)/2
            for eps in eps_range:
                db = DBSCAN( eps = eps, min_samples = min_samples).fit(Y)
                cluster_labels = db.labels_
                unclustered = float(list(cluster_labels).count(-1))/len(cluster_labels)
                # Number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
                cluster_count.append([n_clusters_, unclustered])
            # find the best clustering
            indexes_PF = []
            for i in range(len(cluster_count)):
                cc = cluster_count[i][0]
                un = cluster_count[i][1]
                if cc <= max_cluster_count and un <= max_unclustered_frac :
                    indexes_PF.append(i)
            if len(indexes_PF) == 0:
                problem_clusters.append([g, "no clusters found!"])
                continue
            else:
                if len(indexes_PF) == 1:
                    best_index = indexes_PF[0]
                else:
                    best_index = indexes_PF[0]
                    for i in indexes_PF[1:]:
                        cc = cluster_count[i][0]
                        best_cc = cluster_count[best_index][0]
                        if best_cc == cc:
                            best_index = i
                db = DBSCAN( eps = eps_range[best_index], min_samples = min_samples).fit(Y)
                cluster_labels = db.labels_
                unclustered = float(list(cluster_labels).count(-1))/len(cluster_labels)
                # Number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
            #cluster_labels = ms.labels_
            #cluster_centers = ms.cluster_centers_
            labels_unique = np.unique(cluster_labels)
            cluster_dict = {"cluster_labels" : cluster_labels,
                            "unique_probes": uniq_probes,
                           "all_probes": probes,
                           "plotting": {"upper_limit": upper_limit,
                                        "lower_limit": lower_limit},
                           "mean_shift": cnv_calls[g]["mean_shift"],
                           "dbscan": db,
                           "tsne": cnv_calls[g]["tsne"],
                           "clusters": {}}
            for k in range(n_clusters_):
                my_members = cluster_labels == k
                #cluster_center = cluster_centers[k]
                cluster_case_count = list(np.asarray(labels)[my_members]).count("case")
                cluster_control_count = list(np.asarray(labels)[my_members]).count("control")
                other_case_count = labels.count("case") - cluster_case_count
                other_control_count = labels.count("control") - cluster_control_count
                cluster_table = [[cluster_case_count, cluster_control_count],
                                 [other_case_count, other_control_count]]
                try:
                    cluster_fish = fisher_exact(cluster_table)
                except:
                    cluster_fish = ["na", "na"]
                try:
                    cluster_chi = chi2_contingency(cluster_table)
                except:
                    cluster_chi = ["na", "na", "na", "na"]
                cluster_dict["clusters"][k] = {"cluster_table": cluster_table,
                           "cluster_stats": {"fisher": {"OR": cluster_fish[0],
                                                        "pvalue": cluster_fish[1]},
                                             "chi": {"pvalue": cluster_chi[1],
                                                     "expected": cluster_chi[3]}},
                           "cluster_data": copy_counts[my_members],
                           "cluster_medians": np.median(copy_counts[my_members], axis=0),
                           "cluster_samples": sample_names[my_members],
                           "cluster_members": my_members,
                           "collapsed_data": collapsed_counts[my_members],
                           "collapsed_medians": np.median(collapsed_counts[my_members], axis=0)
                           }
                cluster_output[g] = cluster_dict
        except Exception as e:
            problem_clusters.append([g, e])
    db_output_file = wdir + settings["dbScanOutputFile"]
    with open(db_output_file, "wb") as outfile:
        pickle.dump(cluster_output, outfile)
    return problem_clusters
def plot_clusters(settings, cluster_method="dbscan"):
    wdir = settings["workingDir"]
    if cluster_method == "dbscan":
        cluster_output_file = wdir + settings["dbScanOutputFile"]
    else :
        cluster_output_file = wdir + settings["clusterOutputFile"]
    with open(cluster_output_file, "rb") as infile:
        cnv_calls = pickle.load(infile)
    filtered_tables_file = wdir + settings["filteredTablesFile"]
    with open(filtered_tables_file, "rb") as infile:
        all_tables = pickle.load(infile)
    figsize = tuple(map(int, settings["figsize"]))
    ymax = int(settings["ymax"])
    if cluster_method == "dbscan":
        image_dir = wdir + "dbscan_images/"
    else:
        image_dir = wdir + "cluster_images/"
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    for gene_name in cnv_calls:
        clusters = cnv_calls[gene_name]["clusters"]
        tables = all_tables[gene_name]
        probe_table = tables["probe_information"]
        copies = [probe_table[0,2]]
        locations = [0]
        for i in range(len(probe_table) -1):
            current_copy = probe_table[i+1, 2]
            if current_copy != copies[-1]:
                copies.append(current_copy)
                locations.extend([i, i+1])
        locations.append(len(probe_table))
        locations = [[locations[i], locations[i+1]] for i in range(0, len(locations), 2)]
        upper_limit = cnv_calls[gene_name]["plotting"]["upper_limit"]
        lower_limit = cnv_calls[gene_name]["plotting"]["lower_limit"]
        fig1 = plt.figure(figsize=figsize)
        fig2 = plt.figure(figsize=figsize)
        fig3 = plt.figure(figsize=figsize)
        for c in range(len(clusters)):
            ax1 = fig1.add_subplot(len(clusters), 1, c)
            ax2 = fig2.add_subplot(len(clusters), 1, c)
            ax3 = fig3.add_subplot(len(clusters), 1, c)
            #ax1 = axarray[c]
            #fig, ax1 = plt.subplots()
            plot_data = clusters[c]["cluster_data"]
            median_data = clusters[c]["cluster_medians"]
            collapsed_median = clusters[c]["collapsed_medians"]
            odds = clusters[c]["cluster_stats"]["fisher"]["OR"]
            fish_p = clusters[c]["cluster_stats"]["fisher"]["pvalue"]
            chi_p = clusters[c]["cluster_stats"]["chi"]["pvalue"]
            cluster_size = len(clusters[c]["cluster_samples"])
            cluster_summary_list = ["cluster size " + str(cluster_size),
                                    "odds ratio " + str(odds),
                                    "fisher's exact pvalue " + str(fish_p),
                                    "chi squared pvalue " + str(chi_p)]
            cluster_summary = "\n".join(cluster_summary_list)
            for d in plot_data:
                ax1.plot(d, lw = 0.2)
            ax1.axhline(1.5, lw = 2, c = "r")
            ax1.axhline(2.5, lw = 2, c = "g")
            ax2.plot(median_data, lw = 1, c = "k")
            ax2.axhline(1.5, lw = 2, c = "r")
            ax2.axhline(2.5, lw = 2, c = "g")
            ax3.plot(collapsed_median, lw = 1, c = "k")
            ax3.plot(upper_limit, lw = 2, c = "g")
            ax3.plot(lower_limit, lw = 2, c = "r")
            colors = ["y", "k"]
            for i in range(len(copies)):
                copyname = copies[i]
                loc = locations[i]
                ax1.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 8, c = colors[i%2])
                ax1.text(x = np.mean(loc), y = - 0.5, s = copyname, fontsize = 10,
                         horizontalalignment = "right", rotation = "vertical")
                ax2.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 2, c = colors[i%2])
                ax2.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 8, c = colors[i%2])
                ax2.text(x = np.mean(loc), y = - 0.5, s = copyname, fontsize = 10,
                         horizontalalignment = "right", rotation = "vertical")
            ax1.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            ax1.set_title("Gene " + gene_name + " cluster " + str(c))
            ax1.set_ylim([0,ymax])
            ax2.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            #ax1.annotate(cluster_summary, xy=(-12, -12), xycoords='axes points',
            #    size=10, ha='right', va='top',bbox=dict(boxstyle='round', fc='w'))
            #ax2.annotate(cluster_summary, xy=(-12, -12), xycoords='axes points',
            #    size=8, ha='right', va='top',bbox=dict(boxstyle='round', fc='w'))
            ax2.set_title("Gene " + gene_name + " cluster " + str(c))
            ax2.set_ylim([0,ymax])
            ax3.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            ax3.set_title("Gene " + gene_name + " cluster " + str(c))
            ax3.set_ylim([0,2*ymax])
        # plot clusters
        if cluster_method == "dbscan":
            ms = cnv_calls[gene_name]["dbscan"]
            Y_plot = cnv_calls[gene_name]["tsne"]["Y"]
            T_plot = cnv_calls[gene_name]["tsne"]["V"]
        else:
            ms = cnv_calls[gene_name]["mean_shift"]
            Y_plot = cnv_calls[gene_name]["tsne"]["Y"]
            T_plot = cnv_calls[gene_name]["tsne"]["T"]
        cluster_labels = ms.labels_
        labels_unique = np.unique(cluster_labels)
        colors = plt.cm.Spectral(np.linspace(0, 1, len(labels_unique)))
        n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        fig4 = plt.figure()
        if cluster_method != "dbscan":
            ax4 = fig4.add_subplot(2,1,2)
            ax5 = fig4.add_subplot(2,1,1)
            cluster_centers = ms.cluster_centers_
            for k, col in zip(labels_unique, colors):
                if k == -1:
                    # Black used for noise.
                    col = 'k'
                my_members = cluster_labels == k
                cluster_center = cluster_centers[k]
                ax5.plot(T_plot[my_members, 0], T_plot[my_members, 1],
                         'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=5)
                ax5.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor="w",
                         markeredgecolor='k', markersize=14)
                ax5.annotate(str(k), xy = cluster_center, ha="center",
                             va = "center", fontsize = 12)
        else:
            ax4 = fig4.add_subplot(1,1,1)
        for k, col in zip(labels_unique, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
            my_members = cluster_labels == k
            ax4.plot(Y_plot[my_members, 0], Y_plot[my_members, 1],
                     'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=5)
            ax4.plot(np.mean(Y_plot[my_members, 0]), np.mean(Y_plot[my_members, 1]),
                     'o', markerfacecolor="w",
                     markeredgecolor='k', markersize=14)
            ax4.annotate(str(k), xy = (np.mean(Y_plot[my_members, 0]),
                                       np.mean(Y_plot[my_members, 1])),
                         ha="center",
                         va = "center", fontsize = 12)
        if cluster_method == "dbscan":
            ax4.set_title(gene_name + " DBSCAN clusters")
        else:
            ax4.set_title(gene_name + " meanshift clusters")
            ax5.set_title(gene_name + " tSNE clusters")
        # save figures
        fig1.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig2.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig3.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig1.savefig(image_dir + gene_name, dpi = 96)
        fig2.savefig(image_dir + gene_name + "_medians", dpi = 96)
        fig3.savefig(image_dir + gene_name + "_collapsed", dpi = 96)
        fig4.savefig(image_dir + gene_name + "_meanshift", dpi = 96)
        plt.close("all")
    return
def type_clinical(settings):
    wdir = settings["workingDir"]
    call_info_file = settings["callInfoDictionary"]
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    sample_results_file = wdir + settings["perSampleResults"]
    with open(sample_results_file) as infile:
        samples = json.load(infile)
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    snp_results = {}
    problem_snp_mips = []
    temp = []
    for g in call_info:
        snp_results[g] = {}
        for m in call_info[g]:
            msnps = call_info[g][m]["snps"]
            cop_dic = call_info[g][m]["copies"]
            for s in msnps:
                if msnps[s]["clinical"]:
                    try:
                        snp_results[g][s][m] = {}
                    except KeyError:
                        snp_results[g][s] = {m: {}}
                    for cop in cop_dic:
                        cname = cop_dic[cop]["copyname"] + "_" + cop
                        snp_results[g][s][m][cname] = {}
                    for sam in samples:
                        try:
                            bc_dic = samples[sam][g][m]
                        except KeyError:
                            problem_snp_mips.append([sam, m])
                            continue
                        for c in snp_results[g][s][m]:
                            split_copies = c.split("_")
                            if len(split_copies) == 2:
                                copy_id = split_copies[-1]
                                if copy_id in msnps[s]["copies"]:
                                    try:
                                        haps = bc_dic[c]["filtered_data"]
                                    except KeyError:
                                        continue
                                    for h in haps:
                                        hap_id = h['haplotype_ID']
                                        rbc = h["barcode_count"]
                                        nbc = h["sample_normalized_barcode_count"]
                                        #rbc = round(nbc / samples[sam]["barcode_normalizer"])
                                        # find SNPs in haplotype if any
                                        try:
                                            diffs = unique_haplotypes[m][hap_id]["mapped_copies"]                                                 [copy_id]["differences"]
                                        except KeyError:
                                            continue
                                        snp_found = False
                                        for d in diffs:
                                            if d["clinical_id"] == s:
                                                snp_found = True
                                                if d["base_match"]:
                                                    base_correct = 1
                                                else:
                                                    base_correct = 0
                                                try:
                                                    snp_results[g][s][m][c][sam].append(
                                                        [nbc, rbc, d["ref_base"], d["hap_base"],
                                                         base_correct])
                                                except KeyError:
                                                    snp_results[g][s][m][c][sam] =                                                        [[nbc, rbc, d["ref_base"],
                                                         d["hap_base"], base_correct]]
                                        if not snp_found:
                                            try:
                                                snp_results[g][s][m][c][sam].append(
                                                    [nbc, rbc, msnps[s]["copies"][copy_id]["ref"],
                                                    msnps[s]["copies"][copy_id]["ref"], 0])
                                            except KeyError:
                                                try:
                                                    snp_results[g][s][m][c][sam] =                                                         [[nbc, rbc, msnps[s]["copies"][copy_id]["ref"],
                                                        msnps[s]["copies"][copy_id]["ref"], 0]]
                                                except KeyError:
                                                    try:
                                                        snp_results[g][s][m][c][sam] =                                                         [[nbc, rbc, msnps[s]["copies"][copy_id]["copy_base"],
                                                        msnps[s]["copies"][copy_id]["copy_base"], 0]]
                                                    except KeyError as e:
                                                        temp.append([g,s,m,c,sam, copy_id, e])

    for g in snp_results:
        for s in snp_results[g]:
            for m in snp_results[g][s]:
                for c in snp_results[g][s][m]:
                    sams = snp_results[g][s][m][c]
                    for sam in sams:
                        l = copy.deepcopy(sams[sam])
                        sams[sam] = {"haplotypes": l}
                        if len(l) == 1:
                            if l[0][-1]:
                                hom = 2
                            else:
                                hom = 0
                            genotype = l[0][-2] + "/" + l[0][-2]
                            bc_count = [l[0][1]]
                        elif len(l) == 2:
                            if l[0][-1] == l[1][-1] == 1:
                                hom = 2
                                genotype = l[0][-2] + "/"  + l[1][-2]
                                bc_count = [l[0][1], l[1][1]]
                            elif l[0][-1] == l[1][-1] == 0:
                                hom = 0
                                genotype = l[0][-2] + "/"  + l[1][-2]
                                bc_count = [l[0][1], l[1][1]]
                            else:
                                hom = 1
                                if l[0][-1] == 0:
                                    genotype = l[0][-2] + "/"  + l[1][-2]
                                    bc_count = [l[0][1], l[1][1]]
                                else:
                                    genotype = l[1][-2] + "/"  + l[0][-2]
                                    bc_count = [l[1][1], l[0][1]]


                        else:
                            i = list(map(list, list(zip(*l))))
                            genotype = "/".join(i[3])
                            bc_count = i[1]
                            hom_set = sorted(list(set(i[-1])))
                            if hom_set == [0]:
                                hom = 0
                            elif hom_set == [1]:
                                hom = 2
                            elif hom_set == [0,1]:
                                hom = 1
                            else:
                                hom = 3


                        sams[sam]["homozygosity"] = hom
                        sams[sam]["genotype"] = genotype
                        sams[sam]["barcode_depth"] = bc_count
    snp_counts = {}
    for g in snp_results:
        for s in snp_results[g]:
            snp_counts[s] = {}
            for m in snp_results[g][s]:
                comp = ""
                wt = ""
                het = ""
                mut = ""
                snp_counts[s][m] = {"homozygosity": []}
                snp_counts[s][m]["sample_names"] = []

                for c in snp_results[g][s][m]:
                    sams = snp_results[g][s][m][c]
                    for sam in sams:
                        snp_counts[s][m]["sample_names"].append(sam)
                        snp_counts[s][m]["homozygosity"].append(sams[sam]["homozygosity"])
                        if sams[sam]["homozygosity"] == 0:
                            wt = sams[sam]["genotype"]
                        elif sams[sam]["homozygosity"] == 2:
                            mut = sams[sam]["genotype"]
                        elif sams[sam]["homozygosity"] == 1:
                            het = sams[sam]["genotype"]
                        else:
                            comp = sams[sam]["genotype"]
                snp_counts[s][m]["genotypes"] = [wt, het, mut, comp]
                homs = snp_counts[s][m]["homozygosity"]
                try:
                    snp_counts[s][m]["freq"] = [float(homs.count(i))/len(homs) for i in range(4)]
                except ZeroDivisionError:
                    temp.append([s,m])


    snp_sum = {}
    for g in snp_results:
        snp_sum[g] = {}
        for s in snp_results[g]:
            snp_sum[g][s] = snp_counts[s]
    snp_file = wdir + settings["snpResultsFile"]
    result = {"snp_results": snp_results,
              "snp_summary": snp_sum,
              "snp_counts": snp_counts,
              "problem_snps": problem_snp_mips,
              "errors": temp}
    with open(snp_file, "w") as outfile:
        json.dump(result, outfile, indent=1)
    return result
def plot_clinical(settings):
    wdir = settings["workingDir"]
    snp_file = wdir + settings["snpResultsFile"]
    with open(snp_file) as infile:
        snps = json.load(infile)
    snp_dir = wdir + "clinical_snp_images/"
    if not os.path.exists(snp_dir):
        os.makedirs(snp_dir)
    sample_dir = wdir + "sample_snp_images/"
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    snp_sum = snps["snp_summary"]
    for g in snp_sum:
        sorted_snps = sorted(snp_sum[g])
        if len(sorted_snps) > 0:
            plt.figure()
            snp_names = []
            wt_freqs = []
            het_freqs = []
            hom_freqs = []
            comp_freqs = []
            wt_genotypes = []
            het_genotypes = []
            hom_genotypes = []
            wt_samples = []
            het_samples = []
            hom_samples = []
            for s in sorted_snps:
                for m in snp_sum[g][s]:
                    if len(snp_sum[g][s][m]["homozygosity"]) > 0:
                        snp_names.append(s)
                        wt_freqs.append(snp_sum[g][s][m]["freq"][0])
                        het_freqs.append(snp_sum[g][s][m]["freq"][1])
                        hom_freqs.append(snp_sum[g][s][m]["freq"][2])
                        comp_freqs.append(snp_sum[g][s][m]["freq"][3])
                        wt_genotypes.append(snp_sum[g][s][m]["genotypes"][0])
                        het_genotypes.append(snp_sum[g][s][m]["genotypes"][1])
                        hom_genotypes.append(snp_sum[g][s][m]["genotypes"][2])
                        wt_samples.append(snp_sum[g][s][m]["homozygosity"].count(0))
                        het_samples.append(snp_sum[g][s][m]["homozygosity"].count(1))
                        hom_samples.append(snp_sum[g][s][m]["homozygosity"].count(2))
            ind = np.arange(len(snp_names))
            wt_freqs = np.asarray(wt_freqs)
            het_freqs = np.asarray(het_freqs)
            hom_freqs = np.asarray(hom_freqs)
            comp_freqs = np.asarray(comp_freqs)
            width = 0.75
            p1 = plt.bar(ind, wt_freqs, width, color = "g")
            p2 = plt.bar(ind, het_freqs, width, color = "y", bottom = wt_freqs)
            p3 = plt.bar(ind, hom_freqs, width, color = "r", bottom = het_freqs +wt_freqs)
            p4 = plt.bar(ind, comp_freqs, width, color = "b", bottom = het_freqs + wt_freqs + hom_freqs)
            plt.ylabel("Allele frequency", fontsize = "20" )
            plt.title("Clinically Relevant SNPs in " + g, fontsize = "20")
            freq_list = [wt_freqs, het_freqs, hom_freqs]
            gen_list = [wt_genotypes, het_genotypes, hom_genotypes]
            count_list = [wt_samples, het_samples,hom_samples]
            for j in range(3):
                if j == 0:
                    vertical_pos = freq_list[0]/2
                else:
                    vertical_pos = freq_list[j]/2 + np.sum(freq_list[:j], axis = 0)
                gen = gen_list[j]
                sam_counts = count_list[j]
                for i in range(len(snp_names)):
                    xpos = i + width/2
                    vp = vertical_pos[i]
                    geno = gen[i]
                    gen_ct = sam_counts[i]
                    plt.text(xpos, vp, geno, horizontalalignment = "center")
                    plt.text(xpos, vp - 0.03, gen_ct, horizontalalignment = "center")

            plt.xticks(ind + width/2 , snp_names, fontsize = "15", rotation= "45")

            plt.legend((p1[0], p2[0], p3[0], p4[0]), ("WT", "Het", "Hom", "MultiAllelic"),
                       loc = "lower right")
            plt.savefig(snp_dir + g + "_clinical", dpi = 96)
            plt.close("all")
    return
def clean_alignment(alignment_file, output_file, directory):
    """ Alignment output in .differences format has sequence information
    that is unnecessary. This function removes that information for easier handling."""
    with open (directory + alignment_file, "r") as infile, open(directory + output_file, "w") as outfile:
        counter = 0
        for line in infile:
            newline = line.strip().split('\t')[:-2]
            outline = "\t".join(newline) + '\n'
            outfile.write(outline)
    return
def alignment_parser(wdir, name, spacer = 0):
    alignment_dict = {}
    # open alignment files
    # open(wdir + name + ".differences") as infile:
    with open(wdir + name + ".al") as infile:
        for line in infile:
            newline = line.strip().split("\t")
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
    aligned_regions = {}
    for query in alignment_dict:
        aligned_regions[query] = []
        for a in alignment_dict[query]:
            chrom = a["name1"]
            begin = int(a["zstart1"])
            end = int(a["end1"])
            aligned_regions[query].append([chrom, begin, end])
    overlaps = {}
    for q1 in aligned_regions:
        overlaps[q1] = [q1]
        reg1 = aligned_regions[q1]
        for r1 in reg1:
            for q2 in aligned_regions:
                if q1 == q2:
                    continue
                reg2 = aligned_regions[q2]
                for r2 in reg2:
                    if check_overlap(r1,r2,spacer):
                        overlaps[q1].append(q2)
                        overlap_found = True
                        break
    overlap_found = True
    while overlap_found:
        overlap_found = False
        for o in list(overlaps.keys()):
            if o in overlaps:
                val = overlaps[o]
                for v in val:
                    if (v in overlaps) and (o in overlaps) and (o!=v):
                        overlaps[o].extend(overlaps[v])
                        overlaps.pop(v)
                        overlap_found = True
    for o in overlaps:
        overlaps[o] = sorted(list(set(overlaps[o])))
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
    # return the initial alignment dictionary,
    # each genomic region for each query (aligned_regions)
    # which queries have overlapping alignments with each other (overlaps)
    # for each query in overlaps, non overlapping genomic regions (separated_merged
    return [alignment_dict, aligned_regions, overlaps, separated_merged_regions]
def intraparalog_aligner(resource_dir,
                         target_regions,
                         region_names,
                         overlaps,
                         imperfect_aligners,
                         fasta_sequences,
                         species,
                         identity,
                         coverage,
                         num_process,
                        alignment_options_dict = {}):
    # align all the target regions to each other
    alignment_commands = []
    out_fields = "name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,zstart2+,end2+,length2,identity,coverage"
    gen_out = "general:"+ out_fields
    diff_out = "differences"
    for t in target_regions:
        if len(alignment_options_dict) == 0:
            alignment_options = []
        else:
            alignment_options = alignment_options_dict[t]["options"]
            identity = alignment_options_dict[t]["identity"]
            coverage = alignment_options_dict[t]["coverage"]
        alignment_options.append("--noytrim")
        commands = []
        #if len(target_regions[t]) > 1:
        query_region = target_regions[t][0]
        tar_regs = target_regions[t]
        target_keys = [tr[0] + ":" + str(tr[1] + 1) + "-" + str(tr[2]) for tr in tar_regs]
        query_key = target_keys[0]
        with open(resource_dir + t + ".query.fa", "w") as outfile:
            outfile.write(">" + t + "_ref\n")
            outfile.write(get_sequence(query_key, species))
        with open(resource_dir + t + ".targets.fa", "w") as outfile:
            outfile_list = []
            for i in range(len(target_keys)):
                k = target_keys[i]
                cname = "_C" + str(i)
                outfile_list.append(">" + t + cname)
                outfile_list.append(get_sequence(k, species))
            # add imgt fastas to targets
            ols = overlaps[t]
            o_count = 0
            for o in ols:
                if o in imperfect_aligners:
                    outfile_list.append(">" + t + "_X" + str(o_count))
                    outfile_list.append(fasta_sequences[o])
                    o_count += 1
            outfile.write("\n".join(outfile_list))
        comm = [t + ".query.fa", resource_dir, t + ".aligned",
                resource_dir + t + ".targets.fa",
                ["multiple", "unmask", "nameparse=darkspace"],
                ["unmask", "nameparse=darkspace"],
                identity, coverage,gen_out,
                alignment_options, species]
        alignment_commands.append(comm)
        comm = [t + ".query.fa", resource_dir,
                t + ".differences",
                resource_dir + t + ".targets.fa",
                ["multiple", "unmask", "nameparse=darkspace"],
                ["unmask", "nameparse=darkspace"],
                identity, coverage,
                diff_out, alignment_options, species]
        alignment_commands.append(comm)
    #print alignment_commands
    return align_region_multi(alignment_commands, num_process)
def check_overlap (r1, r2, padding = 0):
    """ Check if two regions overlap. Regions are given as lists of chrom (str),
    begin (int), end (int)."""
    # check chromosome equivalency
    o1 = r1[0] == r2[0]
    # check interval overlap
    merged = merge_overlap([r1[1:], r2[1:]], padding)
    o2 = len(merged) == 1
    return o1 & o2
def alignment_check(region_list):
    """ Check if any region in the given region list produces overlapping alignments.
    Return offending regions. Region list must have two alignment options for each region:
    first one with general output and .al file extension, which will be used
    to check the alignment."""
    alignments = {}
    for i in range(0, len(region_list), 2):
        alignment_dir = region_list[i][1]
        alignment_file = region_list[i][2]
        alignment_dic = alignments[alignment_file.split(".")[0]] = {"overlap": []}
        alignment_list = alignment_dic["alignments"] = []
        full_alignments = alignment_dic["full_alignments"] = []
        with open(alignment_dir + alignment_file) as infile:
            for line in infile:
                if not line.startswith("#"):
                    newline = line.strip().split("\t")
                    alignment_chrom = newline[0]
                    alignment_start = int(newline[1])
                    alignment_end = int(newline[2])
                    alignment_list.append([alignment_chrom, alignment_start, alignment_end])
                    newline[1] = alignment_start
                    newline[2] = alignment_end
                    full_alignments.append(newline)
    for a1 in alignments:
        for a2 in alignments:
            if a1 != a2:
                a_list1 = alignments[a1]["alignments"]
                a_list2 = alignments[a2]["alignments"]
                overlap_found = False
                for i in a_list1:
                    if not overlap_found:
                        for j in a_list2:
                            if check_overlap(i,j):
                                alignments[a1]["overlap"].append(a2)
                                overlap_found = True
                                break
    return alignments
def intra_alignment_checker(family_name,
                            res_dir,
                           target_regions,
                            region_names
                           ):
    alignment_file = family_name + ".aligned"
    new_regions = {}
    with open (res_dir + alignment_file, "r") as alignment:
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
                        new_regions[cn].append([tr[0],
                                                 start,
                                                 end,
                                                0 - len(tr[0]),
                                                 size])
                    except KeyError:
                        new_regions[cn] = [[tr[0],
                                             start,
                                             end,
                                            0 - len(tr[0]),
                                             size]]
    ret_regions =  []
    rnames = []
    for ci in sorted(new_regions):
        ret_regions.extend(sorted(new_regions[ci]))
        if len(new_regions[ci]) > 1:
            for i in range(len(new_regions[ci])):
                rnames.append(
                    region_names[ci] + "-" + str(i)
                )
        else:
            rnames.append(
            region_names[ci]
            )
    return [ret_regions, rnames]
def alignment_mapper(family_name, res_dir):
    """ Once alignment files (.al and .differences) have been created,
    MODIFY .al file to include copyname, copynumber and query columns.
    If copynumbers are "none", copies will be numbered in ascending size.
    1) Create segments from alignments
    2) create rinfo_0 file
    3) return alignment and segment dictionaries, in case they are needed.
    """
    try:
        alignment_file = family_name + ".aligned"
        difference_file = family_name + ".differences"
        with open (res_dir + alignment_file, "r") as alignment,              open(res_dir + difference_file, "r") as difference:
            # create an alignment dictionary for each region that a query aligns to
            # these correspond to each line in the alignment file
            alignment_dic = {}
            alignment_number = 0
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
                    alignment_id = temp_dict["name1"] + ":"                     + str(alignment_number)
                    alignment_dic[alignment_id] = temp_dict
                    cov = float(alignment_dic[alignment_id]["covPct"][:-1])
                    idt = float(alignment_dic[alignment_id]["idPct"][:-1])
                    alignment_dic[alignment_id]["score"] = np.mean([idt, cov])
                    alignment_number += 1
            # differences file is a continuous file for all alignments
            # extract differences for each alignment
            for line in difference:
                newline = line.strip().split("\t")
                dname = newline[0]
                for aname in alignment_dic:
                    if aname.startswith(dname):
                        alignment_dic[aname]["differences"].append(newline[:-2])
            # map each position in each alignment to the query
            for a in alignment_dic:
                snps = alignment_dic[a]["snps"] = {}
                co = alignment_dic[a]["coordinates"] = {}
                rev_co = alignment_dic[a]["reverse_coordinates"] = {}
                alignment_length = int(alignment_dic[a]["length2"])
                # if alignment on reverse strand
                if alignment_dic[a]["strand2"] == "-":
                    # genomic coordinate of target start
                    # this position is zstart2+ away from query end (when it is a - alignment)
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
                        # between the last_key in the coord dic and start_key - diff start
                        for j in range(last_key - 1, query_length - diff_start -1, -1):
                            # j decreases by one, starting from the last available key
                            # the value will be 1 more than the previous key (j+1)
                            if j == last_key -1:
                                co[j] = round(co[j + 1] - 0.1) + 1 + inserted
                            else:
                                co[j] = round(co[j + 1] - 0.1) + 1
                            rev_co[co[j]] = j
                        # current last key is now first_key - diff_start
                        last_key = query_length - diff_start -1
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
                        # in cases of deletion in query, only rev_co will be updated
                        elif diff_start == diff_end:
                            inserted = 0
                            for i in range(tar_end - tar_start):
                                rev_co[co[last_key + 1] + i + 1] = last_key + 0.5
                                inserted += 1
                            last_key += 1
                        # last_key will be mapped to target start
                        # if there is only a SNP and no indel
                        else:
                            inserted = 0
                            co[last_key] = tar_start
                            rev_co[tar_start] = last_key
                        query_diff_start = last_key
                        diff_key = str(query_diff_start) + "-" + str(query_diff_end)
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
                    for k in range(last_key -1, query_plus_start - 1, -1):
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
                        # where on query sequence the difference starts and ends
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
                        # target sequences are the same in length and co dict is filled so
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
                                rev_co[co[last_key - 1] + 1 + i] = last_key - 0.5
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

        return_dic = {}
        for aname in alignment_dic:
            qname = aname.split(":")[0]
            snps = alignment_dic[aname]["snps"]
            co = alignment_dic[aname]["coordinates"]
            rev_co = alignment_dic[aname]["reverse_coordinates"]
            if qname not in return_dic:
                return_dic[qname] = alignment_dic[aname]
                return_dic[qname]["covered"] = [list(map(int, [alignment_dic[aname]["zstart1"],
                                                 alignment_dic[aname]["end1"]]))]
            else:
                qname_cov = [alignment_dic[aname]["zstart1"],
                            alignment_dic[aname]["end1"]]
                return_dic[qname]["covered"].append(list(map(int, qname_cov)))
        return alignment_dic, return_dic
    except Exception as e:
        return [family_name, e]
def make_region (chromosome, begin, end):
    """ Create region string from coordinates.
    takes 2 (1 for human 1-9) digit chromosome,
    begin and end positions (1 indexed)"""
    region="chr" + str(chromosome) + ":" + str(begin) + "-" + str(end)
    return region
def create_region (chromosome, begin, end):
    """ Create region string from coordinates.
    chromosome string,
    begin and end positions (1 indexed)"""
    region=chromosome + ":" + str(begin) + "-" + str(end)
    return region
def create_dirs (dir_name):
    """ create subdirectory names for a given dir,
    to be used by os.makedirs, Return a list of
    subdirectory names."""
    primer3_input_DIR = dir_name + "/primer3_input_files/"
    primer3_output_DIR = dir_name + "/primer3_output_files/"
    bowtie2_input_DIR = dir_name + "/bowtie2_input/"
    bowtie2_output_DIR = dir_name + "/bowtie2_output/"
    mfold_input_DIR = dir_name + "/mfold_input/"
    mfold_output_DIR = dir_name + "/mfold_output/"
    return [primer3_input_DIR, primer3_output_DIR, bowtie2_input_DIR,            bowtie2_output_DIR, mfold_input_DIR, mfold_output_DIR]
def get_coordinates (region):
    """ Define coordinates chr, start pos and end positions
    from region string chrX:start-end. Return coordinate list.
    """
    chromosome = region.split(":")[0]
    coord = region.split(":")[1]
    coord_list = coord.split("-")
    begin = int(coord_list [0])
    end = int(coord_list [1])
    return [chromosome, begin, end]
def get_fasta (region, species="pf", offset=1, header="na"):
    """ Take a region string (chrX:begin-end (1 indexed)),
    and species (human=hs, plasmodium= pf),Return fasta record.
    """
    if offset == 0:
        region_coordinates = get_coordinates(region)
        region = region_coordinates[0] + ":" + str(region_coordinates[1] + 1) +                "-" + str(region_coordinates[2])
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
    file_locations = get_file_locations()
    genome_fasta = file_locations[species]["fasta_genome"]
    region_subs = make_chunks(regions, 25000)
    fasta_dic = {}
    for region_list in region_subs:
        command = ["samtools",  "faidx", genome_fasta]
        command.extend(region_list)
        #print len(region_list)
        out=subprocess.check_output(command)
        fasta_list = out.split(">")[1:]
        for f in fasta_list:
            fl = f.strip().split("\n")
            fhead = fl[0]
            fseq = "".join(fl[1:])
            fasta_dic[fhead] = fseq
    return fasta_dic
def create_fasta_file (region, species, output_file):
    if not os.path.exists(output_file):
        os.makedirs(output_file)
    with open (output_file, "w") as outfile:
        outfile.write(get_fasta(region, species))
def get_snps (region, snp_file):
    """ Take a region string and a  tabix'ed snp file,
    return a list of snps which are lists of
    tab delimited information from the snp file. """
    # extract snps using tabix, in tab separated lines
    snp_temp = subprocess.check_output(["tabix", snp_file, region])
    # split the lines (each SNP)
    snps_split= snp_temp.split("\n")
    # add each snp in the region to a list
    # as lists of
    snps = []
    for line in snps_split:
        snp = line.split('\t')
        snps.append(snp)
    del snps[-1] #remove last item which is coming from the new line at the end
    return snps
def get_vcf_snps(region, snp_file):
    """ Take a region string and a tabix'ed snp file,
    return a list of snps which are lists of
    tab delimited information from the snp file. """
    # extract snps using tabix, in tab separated lines
    snp_temp = subprocess.check_output(["bcftools",
                                        "view",
                                        "-H", "-G", "-r",
                                        region, snp_file])
    # split the lines (each SNP)
    snps_split= snp_temp.split("\n")[:-1]
    # add each snp in the region to a list
    # as lists of
    snps = []
    for line in snps_split:
        snp = line.split('\t')[:8]
        snps.append(snp)
    return snps
def snp_filter_pf (snp_list, col=23, min_allele_freq=0.001, min_total_allele=50):
    """ filter a SNP list using given criteria, print summary,
    return filtered SNP list. col should be a column in the
    original snp file with number of alleles observed as m,n,
    the column ends with a comma as in UCSC schema. It is 23 in
    snp138.
    """
    filtered_snps = []
    for snp in snp_list:
        # number of alleles for each genotype in the form (m,n,..)
        if snp[col] != "":
            acs = [a for a in snp[col].split(",") if a != ""]
            allele_count = list(map(int, list(map(float, acs))))
            # allele with max count
            max_all = max(allele_count)
            # total alleles
            tot_all = sum(allele_count)
            # minor allele freq
            maf = (tot_all - max_all)/float(tot_all)
            if (maf >= min_allele_freq) and (tot_all >= min_total_allele):
                # snp satisfying the given criteria is filtered
                filtered_snps.append(snp)
    print((" ".join(["Found", str(len(filtered_snps)), "SNPs with >=", str(min_allele_freq),
                    "frequency and >=", str(min_total_allele), "reported alleles!"])))
    return filtered_snps
def snp_filter_hs (snp_list, col=23, min_allele_freq=0.001, min_total_allele=50):
    """ filter a SNP list using given criteria, print summary,
    return filtered SNP list. col should be a column in the
    original snp file with number of alleles observed as m,n,
    the column ends with a comma as in UCSC schema. It is 23 in
    snp138.
    """
    filtered_snps = []
    for snp in snp_list:
        # number of alleles for each genotype in the form (m,n,..)
        if snp[col] != "":
            allele_count_list = snp[col].split(",")
            if allele_count_list[-1] == "":
                allele_count_list.pop(-1)
            allele_count = list(map(int, list(map(float, allele_count_list))))
            # allele with max count
            max_all = max(allele_count)
            # total alleles
            tot_all = sum(allele_count)
            # minor allele freq
            maf = (tot_all - max_all)/float(tot_all)
            if (maf >= min_allele_freq) and (tot_all >= min_total_allele):
                # snp satisfying the given criteria is filtered
                filtered_snps.append(snp)
    print((" ".join(["Found", str(len(filtered_snps)), "SNPs with >=", str(min_allele_freq),
                    "frequency and >=", str(min_total_allele), "reported alleles!"])))
    return filtered_snps
def snp_filter_hla (snp_list, col=23, min_allele_freq=0.001, min_total_allele=50, submitter='SI_MHC_SNP'):
    """ filter a SNP list using given criteria, print summary,
    return filtered SNP list. col should be a column in the
    original snp file with number of alleles observed as m,n,
    the column ends with a comma as in UCSC schema. It is 23 in
    snp138.
    """
    filtered_snps = []
    for snp in snp_list:
        submitter_list = snp[20].strip().split(',')[:-1]
        # number of alleles for each genotype in the form (m,n,..)
        if snp[col] != "":
            allele_count = list(map(int, list(map(float, snp[col].split(",")[:-1]))))
            # allele with max count
            max_all = max(allele_count)
            # total alleles
            tot_all = sum(allele_count)
            # minor allele freq
            maf = (tot_all - max_all)/float(tot_all)
            if (maf >= min_allele_freq) and (tot_all >= min_total_allele):
                # snp satisfying the given criteria is filtered
                filtered_snps.append(snp)
            elif submitter in submitter_list:
                filtered_snps.append(snp)
        elif submitter in submitter_list:
            filtered_snps.append(snp)

    print((" ".join(["Found", str(len(filtered_snps)), "SNPs with >=", str(min_allele_freq),
                    "frequency and >=", str(min_total_allele), "reported alleles!"])))
    return filtered_snps
def snp_function_filter_hs (snp_list, func_list=["nonsense", "missense","stop-loss", "frameshift",                            "cds-indel", "splice-3","splice-5"], col=15):
    """ Take a snp list and filter for functional consequence of that snp.
    Functionality of a snp is located in col (index) 15 of UCSC table and 14 of include file.
    It can be any type from the list ['unknown', 'coding-synon', 'intron', 'near-gene-3',
   'near-gene-5', 'ncRNA', 'nonsense', 'missense', 'stop-loss', 'frameshift', 'cds-indel',
    'untranslated-3', 'untranslated-5', 'splice-3', 'splice-5'].
    """
    # define available functional changes
    snp_func_list = ['unknown', 'coding-synon', 'intron', 'near-gene-3', 'near-gene-5',                  'ncRNA', 'nonsense', 'missense', 'stop-loss', 'frameshift', 'cds-indel',                 'untranslated-3', 'untranslated-5', 'splice-3', 'splice-5']
    # check if the filter criteria is in the available changes
    for f in func_list:
        if f not in snp_func_list:
            return "Filter criteria", f, " is not available in snp list."
    filtered_snps = []
    for snp in snp_list:
        if snp[col] in func_list:
            filtered_snps.append(snp)
    print(("Found ", len(filtered_snps), " SNPs with functional significance"))
    return filtered_snps
def targets (must_file, diff_list):
    """ Take a file with snps or regions that must be captured by mips and a list
    of other variations of interest (created in ucsc table format) and return a target
    dictionary"""
    print(("Extracting target snps from file " + must_file))
    # create a dict of positions that should be targeted
    targets = {}
    # top priority will be given to a must have list of rs numbers
    # or genomic coordinates. Read from a tab separated file that has
    # id, chr, beg, end (None if not specified). rs numbers are not used at this time.
    # so, genomic coordinates must be provided
    with open(must_file, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                snp_name = newline[0]
                snp_chr = newline[1]
                snp_begin = newline[2]
                snp_end = newline[3]
                # from the coordinates provided, create the dict key
                key = snp_chr + ":" + snp_begin + "-" + snp_end
                # add snp to dictionary
                targets[key] = {"chrom":snp_chr, "begin":snp_begin, "end":snp_end,                                 "name": snp_name, "diff": "na", "source": "must"}
    # add snps from the diff_list
    for diff in diff_list:
        snp_name = diff[4]
        snp_chr = diff[1]
        snp_begin = diff[2]
        snp_end = diff[3]
        src = diff[15]
        # from the coordinates, create the dict key
        key = snp_chr + ":" + snp_begin + "-" + snp_end
        # check if the key is already in the dict
        if not key in list(targets.keys()):
            # add snp to dictionary
            targets[key] = {"chrom":snp_chr, "begin":snp_begin, "end":snp_end,                             "name": snp_name, "diff": "na", "source": src}
    # return targets dictionary
    return targets
def merge_overlap (intervals, spacer=0):
    """ Merge overlapping intervals. Take a list of lists of 2 elements, [start, stop],
    check if any [start, stop] pairs overlap and merge if any. Return the merged [start, stop]
    list."""
    exons = copy.deepcopy(intervals)
    exons = [e for e in exons if len(e) == 2]
    for e in exons:
        e.sort()
    exons.sort()
    if len(exons) < 2:
        return exons
    # reuse a piece of code from get_exons:
    #######################################
    overlapping = 1
    while overlapping:
        overlapping = 0
        for i in range(len(exons)):
            e = exons[i]
            for j in range(len(exons)):
                x = exons[j]
                if i == j :
                    continue
                else:
                    if e[1] >= x[1]:
                        if (e[0] - x[1]) <= spacer:
                            overlapping = 1
                    elif x[1] >= e[1]:
                        if (x[0] - e[1]) <= spacer:
                            ovrelapping = 1
                    if overlapping:
                        # merge exons and add to the exon list
                        exons.append([min(e[0], x[0]), max(e[1], x[1])])
                        # remove the exons e and x
                        exons.remove(e)
                        exons.remove(x)
                        # once an overlapping exon is found, break out of the for loop
                        break
            if overlapping:
                # if an overlapping exon is found, stop this for loop and continue with the
                # while loop with the updated exon list
                break
    exons.sort()
    return exons
def overlap(reg1, reg2, spacer = 0):
    """
    Return overlap between two regions.
    e.g. [10, 30], [20, 40] returns [20, 30]
    """
    regions = sorted([sorted(reg1), sorted(reg2)])
    try:
        if regions[0][1] - regions[1][0] >= spacer:
            return[max([regions[0][0], regions[1][0]]),
                   min([regions[0][1], regions[1][1]])]
        else:
            return []
    except IndexError:
        return []
def remove_overlap(reg1, reg2, spacer = 0):
    """
    Remove overlap between two regions.
    e.g. [10, 30], [20, 40] returns [10, 20], [30, 40]
    """
    regions = sorted([sorted(reg1), sorted(reg2)])
    try:
        if regions[0][1] - regions[1][0] >= spacer:
            coords = sorted(reg1 + reg2)
            return[[coords[0], coords[1]],
                   [coords[2], coords[3]]]
        else:
            return regions
    except IndexError:
        return []
def subtract_overlap (uncovered_regions, covered_regions, spacer = 0):
    """
    Given two sets of regions in the form
    [[start, end], [start, end]], return a set
    of regions that is the second set subtracted from
    the first.
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
        uncovered = [[uncovered_remaining[i-1],
                      uncovered_remaining[i]] for i in \
                     range(1, len(uncovered_remaining))\
                    if uncovered_remaining[i] - uncovered_remaining[i-1]\
                     > 1]
        unc = [uncovered_remaining[0]]
        for u in uncovered:
            unc.extend(u)
        unc.append(uncovered_remaining[-1])
        return [[unc[i], unc[i+1]]for i in range(0, len(unc), 2)               if unc[i+1] - unc[i] > spacer]
    else:
        return []
def trim_overlap (region_list, low = 0.1, high = 0.9, spacer = 0):
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
                            overlapping_region = overlap(reg_i, reg_j, spacer)
                            if len(overlapping_region) > 0:
                                reg_sizes = sorted([reg_i[1] - reg_i[0] + 1,
                                                    reg_j[1] - reg_j[0] + 1])
                                overlap_size = float(overlapping_region[1]                                                     - overlapping_region[0])
                                overlap_ratio = overlap_size/reg_sizes[0]
                                if overlap_ratio <= low:
                                    region_list[i] = "remove"
                                    region_list[j] = "remove"
                                    region_list.extend(remove_overlap(reg_i, reg_j, spacer))
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
                                    print((overlap_ratio, "is outside trim range for ",                                    reg_i, reg_j))

    return region_list
def get_exons (gene_list):
    """Take a list of transcript information in refgene format and return a list of exons
    in the region as [[e1_start, e1_end], [e2_start], [e2_end], ..]. The transcripts must
    belong to the same gene (i.e. have the same gene name).
    Merge overlapping exons. """
    # get start and end coordinates of exons in gene list
    starts = []
    ends = []
    gene_names = []
    gene_ids= []
    #print gene_list
    chrom_list = []
    for gene in gene_list:
        chrom_list.append(gene[2])
    chrom_set = list(set(chrom_list))
    if len(chrom_set) == 0:
        return {}
    while_counter = 0
    while (len(chrom_set)>1) and (while_counter < 50):
        for c in chrom_set:
            if len(c)>5:
                chrom_set.remove(c)
        while_counter += 1
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
                if (i != j) and                    ((e[0] <= x[0] <= e[1]) or (e[0] <= x[1] <= e[1]) or (x[0] <= e[0] <= x[1])):
                    # merge exons and add to the exon list
                    exons.append([min(e[0], x[0]), max(e[1], x[1])])
                    # remove the exons e and x
                    exons.remove(e)
                    exons.remove(x)
                    # change overlapping to 1 so we can stop the outer for loop
                    overlapping = 1
                    # once an overlapping exon is found, break out of the for loop
                    break
            if overlapping:
                # if an overlapping exon is found, stop this for loop and continue with the
                # while loop with the updated exon list
                break
    # get the gene start and end coordinates
    if (len(starts) >= 1) and (len(ends)>=1):
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
        genes = get_snps(region, get_file_locations()[species]["refgene_tabix"])
        for g in genes:
            gene_names.append(g[12])
    except KeyError:
        pass
    return gene_names
def get_gene (gene_name, refgene_file, chrom=None, alternative_chr=1):
    """ Return genomic coordinates of a gene extracted from the refseq genes file.
    Refgene fields are as follows:
    0:bin, 1:name, 2:chrom, 3:strand, 4:txStart, 5:txEnd, 6:cdsStart, 7:cdsEnd, 8:exonCount,
    9:exonStarts, 10:exonEnds, 11:score, 12:name2, 13:cdsStartStat, 14:cdsEndStat, 15:exonFrames.
    Field 12 will be used for name search."""
    # all chromosomes must be included if chromosome of the gene is not provided
    # therefore, chrom cannot be None when alternative_chr is set to 0
    if not (chrom or alternative_chr):
        print("Chromosome of the gene %s must be specified or all chromosomes must be searched.")
        print(("Specify a chromosome or set alternative chromosome to 1." %gene_name))
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
def create_gene_fasta (gene_name_list, wdir, species = "hs", flank=150, multi_file=False):
    """ Get a list of genes, extract exonic sequence + flanking sequence.
    Create fasta files in corresponding directory for each gene if multi_file is True,
    create a single fasta file if False.
    """
    region_list = []
    for gene_name in gene_name_list:
        if gene_name.startswith("chr"):
            coord = get_coordinates(gene_name)
            query = make_region(coord[0], coord[1] - flank, coord[2] + flank)
        else:
            e = get_exons(
                get_gene(gene_name, get_file_locations()[species]["refgene"], alternative_chr=1)
                )
            query = e["chrom"] + ":" + str(e["begin"] - flank) + "-" + str(e["end"] + flank)
        region_list.append(query)
    regions = get_fasta_list(region_list, species)
    if multi_file:
        for i in range(len(region_list)):
            r = region_list[i]
            gene_name = gene_name_list[i]
            filename = wdir + gene_name + ".fa"
            with open(filename, "w") as outfile:
                outfile.write("\n".join())
    else:
        with open(wdir + "multi.fa", "w") as outfile:
            outfile.write("\n".join(region_list))
def get_region_exons(region, species):
    try:
        genes = get_snps(region, get_file_locations()[species]["refgene_tabix"])
    except KeyError:
        genes = []
    return get_exons(genes)
def get_cds(gene_name, species):
    gene_list = get_gene(gene_name,
                         get_file_locations()[species]["refgene"],
                         alternative_chr=1)
    if len(gene_list) > 1:
        print(("More than one refgene entry was found for the gene ", gene_name))
        print("Exons from alternative transcripts will be merged and CDS will be generated from that.")
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
def make_targets_dic(targets_file):
    with open (targets_file) as infile:
        # create a targets dictionary from targets file
        targets = {}
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                name = newline[0]
                if name not in targets:
                    targets[name] = {"must":[]}
                if newline[1] != "target":
                    targets[name][newline[1]] = newline[2]
                else:
                    must_info = newline[2].split(",")
                    must_dic = {
                                "name": must_info[0],
                                "snp_id": must_info[1],
                                "chrom": must_info[2],
                                "begin": must_info[3],
                                "end": must_info[4],
                                "aa_begin":must_info[5],
                                "aa_end":must_info[6],
                                "weight": "none",
                                "segment": "S0:C0"
                                }
                    targets[name]["must"].append(must_dic)
    for t in targets:
        species = targets[t]["species"]
        try:
            gname = targets[t]["gene_name"]
        except KeyError:
            gname = targets[t]["gene_name"] = t
        if species == "pf" and gname == "none":
            try:
                alias_handle= open(get_file_locations()[species]["alias"])
                alias_dic = json.load(alias_handle)
                gname = alias_dic[t]
                alias_handle.close()
            except KeyError:
                gname = t
        flank = int(targets[t]["flank"])
        exons = get_exons(get_gene(gname, get_file_locations()[species]["refgene"], alternative_chr=1))
        cds = get_cds(gname, species)
        try:
            cds_coord = cds["coordinates"]
            targets[t]["cds"] = cds["sequence"]
        except KeyError:
            cds_coord = "none"
            targets[t]["cds"] = "none"
        try:
            targets[t]["orientation"] = exons["orientation"]
        except KeyError:
            targets[t]["orientation"] = "+"
        if (targets[t]["chrom"] == "none") or (targets[t]["begin"] == "none") or              (targets[t]["end"] == "none"):
            try:
                if targets[t]["chrom"] == "none":
                    targets[t]["chrom"] = exons["chrom"]
                targets[t]["begin"] = exons["begin"] - flank + 1
                targets[t]["end"] = exons["end"] + flank
            except KeyError:
                print(("cannot find coordinates for ", t))
                targets[t]["chrom"] = "none"
                targets[t]["begin"] = "none"
                targets[t]["end"] = "none"
        else:
            targets[t]["begin"] = int(targets[t]["begin"])
            targets[t]["end"] = int(targets[t]["end"])
        if species != "hs":
            for m in targets[t]["must"]:
                if m["snp_id"] == "none" and m["begin"] == "none" and m["aa_begin"] != "none":
                    aa_begin = int(m["aa_begin"]) - 1
                    aa_end = int(m["aa_end"])
                    nt_begin = 3 * aa_begin
                    nt_end = 3 * aa_end - 1

                    if targets[t]["orientation"] == "+":
                        m["begin"] = cds_coord[nt_begin] + 1
                        m["end"] = cds_coord[nt_end] + 1
                    elif targets[t]["orientation"] == "-":
                        m["end"] = cds_coord[nt_begin] + 1
                        m["begin"] = cds_coord[nt_end] + 1
                    m["chrom"] = targets[t]["chrom"]


    return targets
def make_subregions (target_region,
                     gene_name,
                     capture_type,
                     capture_targets, maf, mta, refseq, flank):
    """ Take a region with its genomic coordinates, a dictionary of possible targets in the region
    which is the output of targets function. Return subregions depending on the capture type that
    is planned. Capture type can be exons, in which case the sequence is split into
    exons unless the intron between two exons are too small, than those exons will be merged.
    Second option for capture type is targets only, which eliminates sequence regions devoid of
    targets. """

    # define possible capture types and check if specified type is correct
    possible_capture_types = ["exons", "targets", "whole"]
    if capture_type not in possible_capture_types:
        print(("Specified capture type %s is not available." %capture_type))
        return 1
    # convert region to coordinates
    coord = get_coordinates(target_region)
    chrom = coord[0]
    begin = int(coord[1])
    end = int(coord[2])
    # create a list of subregions
    sub = []
    # create an output list
    subregions = []
    # get exon coordinates if the capture type is exons
    # this function is not available at the moment but it should return
    # a list of coordinates in [[e1_begin, e1_end], [e2_begin, e2_end]]
    # subregion will be the target region if all of the capture type is "whole"
    if capture_type == "whole":
        subregions = [[begin, end]]
    elif capture_type == "exons":
        exons = get_exons(get_snps(target_region, hs_refgene_tabix), refseq)["exons"]
        for i in range(len(exons) - 1):
            # define the current exon
            current_exon = exons[i]
            # and its coordinates
            e_start = current_exon[0]
            e_end = current_exon[1]
            # the next exon in the list
            next_exon = exons[i+1]
            # and its coordinates
            n_start = next_exon[0]
            n_end = next_exon[1]
            # the intron between should be at least 2 flanks long so that we can have
            # non-overlapping templates. Append a 1 for standalone exons and 0 for exons to be merged
            if (n_start - e_end - 1) > (2 * flank):
                sub.append(1)
            else:
                # if the intron is too small, current and next exons should be merged
                sub.append(0)
        # merge each exon that has a 0 value in sub list, until no more exons left to be merged
        while 0 in sub:
            for i in range(len(sub)):
                if sub[i] == 0:
                    exons[i] = [exons[i][0], exons[i+1][1]]
                    exons.pop(i+1)
                    sub.pop(i)
                    break
        subregions = exons
    # get the target coordinates if the capture type is "targets"
    elif capture_type == "targets":
        # capture targets should be a target dictionary with "chrx:begin-end" keys.
        # create a list of target coordinates
        target_coordinates = []
        for t in capture_targets:
            # add begin coordinate of the target
            target_coordinates.append([int(capture_targets[t]["begin"]),
                                       int(capture_targets[t]["end"])])
        # sort target coordinates
        target_coordinates.sort()
        # create an end coordinates list to have index numbers coordinates where target region should
        # be split because the next target is too far away
        end_coordinates = []
        for i in range(len(target_coordinates)-1):
            if (target_coordinates[i+1] - target_coordinates[i]) > (2 * flank):
                end_coordinates.append(target_coordinates[i])
                end_coordinates.append(target_coordinates[i+1])
        # initialize the subregion list with the first target
        subregion_list = [target_coordinates[0]]
        subregion_list.extend(end_coordinates)
        # finalize list with the last target
        subregion_list.append(target_coordinates[-1])
        for i in range(0, len(subregion_list), 2):
            subregions.append([subregion_list[i], subregion_list[i+1]])
    regs = []
    counter = 0
    for reg in subregions:
        tar_region = chrom + ":" + str(reg[0]) + "-" + str(reg[1])
        regs.append([gene_name + "_S" + str(counter), tar_region])
        counter += 1
    return regs
def mask_snps (region, snp_file, count_col, min_maf=0.03, min_count=10):
    """lower case mask SNPs in plasmodium genome, return fasta record"""
    # get sequence of the region
    region_fasta = get_fasta(region)
    # extract fasta header
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    # join separate fasta lines into single string
    seq_template_temp = "".join(fasta_list[1:])
    # convert the sequence into list of nucleotides
    seq_template_list = list(seq_template_temp)
    # find snps in the region using get_snps
    region_snps = snp_filter_pf(get_snps(region, snp_file),                  count_col, min_maf, min_count)
    # extract region begin coordinate
    coordinates = get_coordinates(region)
    for snp in region_snps:
        # find snp location in the sequence
        # by subtracting sequence start coordinate from
        # snp coordinate
        snp_index = int(snp[1]) - int(coordinates[1])
        # change the SNP nucleotide to lower case
        seq_template_list[snp_index] = seq_template_list[snp_index].lower()
    # create sequence string from nucleotide list
    fasta_seq = "".join(seq_template_list)
    # add fasta header and create a fasta record of
    # masked sequence
    fasta_rec = fasta_head + "\n" + fasta_seq
    return fasta_rec
def exclude_snps (region, snp_file, maf_col=6, count_col=9, min_maf=0.03, min_count=10):
    """Prepare an exclude list to be used in boulder record.
    Gets SNPs from plasmodium snp file and assumes SNP length
    of 1."""
    # get the region sequence
    region_fasta = get_fasta(region)
    # extract fasta header
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    # join fasta lines in a single string
    seq_template_temp = "".join(fasta_list[1:])
    # convert sequence to list of nucleotides
    seq_template_list = list(seq_template_temp)
    # find region snps filtered for maf
    region_snps = snp_filter_pf(get_snps(region, snp_file),             maf_col, count_col, min_maf, min_count)
    # find sequence start coordinate
    coordinates = get_coordinates (region)
    exclude = []
    for snp in region_snps:
        # find snp location in sequence by
        # subtracting sequence start from snp coordinate
        snp_index = int(snp[1])-int(coordinates[1])
        # create an exclude list for snp index with length 1
        exclude.append([snp_index,1])
    return exclude
def exclude_snps_hla (region, col=23, min_frequency=0.001, min_total_alleles=50):
    """Prepare an exclude list to be used in boulder record.
    Gets SNPs from human snp file and assumes SNP length
    of 1."""
    # get fasta sequence of region using the function get_fasta
    region_fasta = get_fasta(region, species="hs")
    # convert fasta record to one line string without the fasta identifier
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    seq_template_temp = "".join(fasta_list[1:])
    # convert fasta string to list of characters
    seq_template_list = list(seq_template_temp)
    # convert region string to coordinates using get_coordinates function
    coordinates = get_coordinates (region)
    # create a SNP list that is filtered for frequency and allele number
    snp_list = snp_filter_hla(get_snps(region, snp_data_hs), col, min_frequency, min_total_alleles)
    # create a list of excluded snp coordinates
    exclude = []
    # find locations of snps in the sequence
    for snp in snp_list:
        # get snp coordinates relative to the sequence at hand
        snp_index_start = int(snp[2]) + 1 - int(coordinates[1])
        snp_index_end = int(snp[3]) - int(coordinates[1])
        exclude.append([snp_index_start, snp_index_end - snp_index_start + 1])
    # sort excluded snps by their location
    exclude.sort(key=itemgetter(0))
    # merge snps that are close together to reduce excluded region number

    for i in range(len(exclude)-1):
        l_current = exclude[i]
        l_next = exclude[i+1]
        if (l_next[0] - sum(l_current)) < 18:
            l_new = [l_current[0], l_next[0] - l_current[0] + l_next[1]]
            exclude[i+1] = l_new
            exclude[i] = "delete"
    excluded = [x for x in exclude if x != "delete"]
    return excluded
def exclude_snps_hs (region, col=23, min_frequency=0.001, min_total_alleles=50):
    """Prepare an exclude list to be used in boulder record.
    Gets SNPs from human snp file and assumes SNP length
    of 1."""
    snp_data_hs = get_file_locations()["hs"]["snps"]
    # get fasta sequence of region using the function get_fasta
    region_fasta = get_fasta(region, species="hs")
    # convert fasta record to one line string without the fasta identifier
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    seq_template_temp = "".join(fasta_list[1:])
    # convert fasta string to list of characters
    seq_template_list = list(seq_template_temp)
    # convert region string to coordinates using get_coordinates function
    coordinates = get_coordinates (region)
    # create a SNP list that is filtered for frequency and allele number
    snp_list = snp_filter_hs(get_snps(region, snp_data_hs), col, min_frequency, min_total_alleles)
    # create a list of excluded snp coordinates
    exclude = []
    # find locations of snps in the sequence
    for snp in snp_list:
        # get snp coordinates relative to the sequence at hand
        snp_index_start = int(snp[2]) + 1 - int(coordinates[1])
        snp_index_end = int(snp[3]) - int(coordinates[1])
        exclude.append([snp_index_start, snp_index_end - snp_index_start + 1])
    # sort excluded snps by their location
    exclude.sort(key=itemgetter(0))
    # merge snps that are close together to reduce excluded region number

    for i in range(len(exclude)-1):
        l_current = exclude[i]
        l_next = exclude[i+1]
        if (l_next[0] - sum(l_current)) < 18:
            l_new = [l_current[0], l_next[0] - l_current[0] + l_next[1]]
            exclude[i+1] = l_new
            exclude[i] = "delete"
    excluded = [x for x in exclude if x != "delete"]
    return excluded
def mask_snps_hla (region, col=23, min_frequency=0.001, min_total_alleles=50):
    """lowercase mask SNPs in a human genomic region
    return fasta record."""
    snp_data_hs = get_file_locations()["hs"]["snps"]
    # get fasta sequence of region using the function get_fasta
    region_fasta = get_fasta(region, species="hs")
    # convert fasta record to one line string without the fasta identifier
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    seq_template_temp = "".join(fasta_list[1:])
    # convert fasta string to list of characters
    seq_template_list = list(seq_template_temp)
    # convert region string to coordinates using get_coordinates function
    coordinates = get_coordinates (region)
    # create a SNP list that is filtered for frequency and allele number
    snp_list = snp_filter_hla(get_snps(region, snp_data_hs), col, min_frequency, min_total_alleles)
    # convert sequence at SNP locations to lower case
    for snp in snp_list:
        snp_index_start = int(snp[2]) + 1 - int(coordinates[1])
        snp_index_end = int(snp[3]) - int(coordinates[1])
        """
        for i in range(snp_index_start, snp_index_end + 1):
            try:
                seq_template_list[i] = seq_template_list[i].lower()
            except IndexError:
                print i, snp_index_start, snp_index_end, seq_template_temp, snp
         """

        if snp_index_start < 0:
            snp_index_start = 0
        if snp_index_end > len(seq_template_list) - 1:
            snp_index_end = len(seq_template_list) - 1

        """
        if (snp_index_start >= 0) and (snp_index_end <= len(seq_template_list) -1):

            for i in range(snp_index_start, snp_index_end + 1):
                seq_template_list[i] = seq_template_list[i].lower()
        """
        for i in range(snp_index_start, snp_index_end + 1):
            seq_template_list[i] = seq_template_list[i].lower()
    # rebuild the fasta record from modified list
    fasta_seq = "".join(seq_template_list)
    fasta_rec = fasta_head + "\n" + fasta_seq
    return fasta_rec
def mask_snps_hs (region, col=23, min_frequency=0.001, min_total_alleles=50):
    """lowercase mask SNPs in a human genomic region
    return fasta record."""
    snp_data_hs = get_file_locations()["hs"]["snps"]
    # get fasta sequence of region using the function get_fasta
    region_fasta = get_fasta(region, species="hs")
    # convert fasta record to one line string without the fasta identifier
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    seq_template_temp = "".join(fasta_list[1:])
    # convert fasta string to list of characters
    seq_template_list = list(seq_template_temp)
    # convert region string to coordinates using get_coordinates function
    coordinates = get_coordinates (region)
    # create a SNP list that is filtered for frequency and allele number
    snp_list = snp_filter_hs(get_snps(region, snp_data_hs), col, min_frequency, min_total_alleles)
    # convert sequence at SNP locations to lower case
    for snp in snp_list:
        snp_index_start = int(snp[2]) + 1 - int(coordinates[1])
        snp_index_end = int(snp[3]) - int(coordinates[1])
        for i in range(snp_index_start, snp_index_end + 1):
            seq_template_list[i] = seq_template_list[i].lower()
    # rebuild the fasta record from modified list
    fasta_seq = "".join(seq_template_list)
    fasta_rec = fasta_head + "\n" + fasta_seq
    return fasta_rec
def make_boulder (fasta,
                  primer3_input_DIR,
                  exclude_list=[],
                  output_file_name="",
                  sequence_targets=[]):
    """ create a boulder record file in primer3_input_DIR from a given fasta STRING.
    SEQUENCE_ID is the fasta header, usually the genomic region (chrX:m-n)
    exclude_list is [coordinate,length] of any regions primers cannot overlap.
    """
    # parse fasta string, get header and remove remaining nextlines.
    fasta_list=fasta.split("\n")
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
        sequence_target_string = " ".join([",".join(map(str, s))                                         for s in sequence_targets])
    boulder = ("SEQUENCE_ID=" + fasta_head + "\n" +
                "SEQUENCE_TEMPLATE="+seq_template+"\n"+
                "SEQUENCE_TARGET=" + sequence_target_string + "\n" +
                "SEQUENCE_EXCLUDED_REGION="+exclude_region+"\n"+ "=")
    if output_file_name == "":
        outname = fasta_head
    else:
        outname = output_file_name
    with open(primer3_input_DIR + outname, 'w') as outfile:
        outfile.write(boulder)
    return boulder
def snp_masker(wdir,
               output_name,
               region_key,
               species,
               masking= 0,
               maf=0.0001,
               mac=10,
               sequence_targets = []):
    region_snps = get_snps(region_key,
                           get_file_locations()[species]["snps"])
    filtered_snps = snp_filter_hs(region_snps,
                                  min_allele_freq = maf,
                                  min_total_allele = mac)
    begin = get_coordinates(region_key)[1]
    region_fasta = get_fasta(region_key, species).upper()
    fasta_list=region_fasta.split("\n")
    fasta_head = fasta_list[0]
    seq_template_temp = "".join(fasta_list[1:])
    exclude = []
    for d in filtered_snps:
        snp_index_start = int(d[2]) - begin
        snp_index_end = int(d[3]) - begin + 1
        if not masking:
            exclude.append([snp_index_start, snp_index_end - snp_index_start])
        else:
            for i in range(snp_index_start, snp_index_end):
                seq_template_temp[i] = seq_template_temp[i].lower()
            else:
                exclude_ext.append([snp_index_start, snp_index_end - snp_index_start])
    # sort excluded snps by their location
    exclude.sort(key=itemgetter(0))
    # merge snps that are close together to reduce excluded region number
    #print exclude
    for i in range(len(exclude)-1):
        l_current = exclude[i]
        l_next = exclude[i+1]
        if 0 <= l_next[0] - sum(l_current) < 18:
            l_new = [l_current[0], l_next[0] - l_current[0] + l_next[1]]
            exclude[i+1] = l_new
            exclude[i] = "delete"
        elif (sum(l_next) - sum(l_current)) <= 0:
            exclude[i+1] = exclude[i]
            exclude[i] = "delete"
        elif sum(l_next) > sum(l_current) > l_next[0]:
            l_new = [l_current[0], l_next[0] - l_current[0] + l_next[1]]
            exclude[i+1] = l_new
            exclude[i] = "delete"
    excluded = [x for x in exclude if x != "delete"]
    # rebuild the fasta record from modified list
    fasta_seq = "".join(seq_template_temp)
    fasta_rec = fasta_head[1:] + "\n" + fasta_seq
    make_boulder (fasta_rec, wdir,exclude_list=excluded,                      output_file_name=output_name,
                      sequence_targets = sequence_targets)
    return
def make_primers (input_file,  settings, primer3_input_DIR, primer3_output_DIR, output_file="input_file"):
    """ make primers using boulder record file in primer3_input_DIR
    using settings file in primer3_settings_DIR and output as boulder
    record to primer3_output_DIR"""
    file_locations = get_file_locations()
    primer3_settings_DIR = file_locations["all"]["primer3_settings_DIR"]
    # if an output file is specified:
    if output_file != "input_file":
        primer3_out = output_file
    # if no output file is specified, name of the file is the same as input file.
    else:
        primer3_out = input_file
    # call primer3 program using the input and settings file
    primer3_output = subprocess.check_output(["primer3_core", "-p3_settings_file="+primer3_settings_DIR+settings, primer3_input_DIR + input_file])
    # write boulder record to file. Append the settings used to output file name.
    outfile = open (primer3_output_DIR + primer3_out + "_" + settings, 'w')
    outfile.write(primer3_output)
    outfile.close()
    return
def make_primers_worker (l):
    """ A worker function to make primers for multiple regions
    using separate processors. Reads boulder record in given input
    directory and creates primer output files in output directory"""
    file_locations = get_file_locations()
    primer3_settings_DIR = file_locations["all"]["primer3_settings_DIR"]
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
    # call primer3 program using the input and settings file
    primer3_output = subprocess.check_output(["primer3_core",
                                              "-p3_settings_file="+primer3_settings_DIR+settings,
                                              primer3_input_DIR + input_file])
    # write boulder record to file. Append the settings used to output file name.
    outfile = open (primer3_output_DIR + output_file, 'w')
    outfile.write(primer3_output)
    outfile.close()
    return
def make_primers_multi(ext_list, lig_list, pro):

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
def primer_parser3 (input_file,
                    primer3_output_DIR,
                    bowtie2_input_DIR,
                    parse_out, fasta=1, outp=1):
    """ parse a primer3 output file and generate a fasta file in bowtie
    input directory that only contains primer names and sequences to be
    used as bowtie2 input.
    Return a dictionary {sequence_information:{}, primer_information{}}
    first dict has tag:value pairs for input sequence while second dict
    has as many dicts as the primer number returned with primer name keys
    and dicts as values {"SEQUENCE": "AGC..", "TM":"58"...}. Also write
    this dictionary to a json file in primer3_output_DIR. """
    primer_dic = {}
    # all target sequence related information will be placed in
    # sequence_information dictionary.
    primer_dic["sequence_information"] = {}
    # primer information will be kept in primer_information dicts.
    primer_dic["primer_information"] = {}
    # load the whole input file into a list.
    infile = open(primer3_output_DIR+input_file, 'r')
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
        value = pair [1]
        if tag.startswith("SEQUENCE"):
            if tag == "SEQUENCE_ID":
                new_value = value.split(",")[-1].replace("CHR", "chr")
                """
                # if chromosome name is chrx or chry, they should be changed to
                # chrX and chrY.
                temp_coord = get_coordinates(new_value)
                temp_chrom = temp_coord[0]
                if temp_chrom == "chrx":
                    temp_chrom = "chrX"
                elif temp_chrom == "chry":
                    temp_chrom = "chrY"
                # if the chromosome name is coming from an alternative assembly
                # it looks like CHR7_KI270803V1_ALT in primer file while
                # in the original files it looks like chr7_KI270803v1_alt
                # so "ki" must be replaced with "KI"
                elif "ki" in temp_chrom:
                    temp_chrom = temp_chrom.replace("ki", "KI")
                new_value = temp_chrom + ":" + str(temp_coord[1]) + "-" +\
                            str(temp_coord[2])
                """
                primer_dic["sequence_information"][tag]= new_value
            else:
                primer_dic["sequence_information"][tag] = value
    # find how many left primers returned and create empty dictionary
    # for each primer in primer_information dict.
    for pair in lines:
        tag = pair[0]
        value = pair [1]
        if tag == "PRIMER_LEFT_NUM_RETURNED":
            # Add this to sequence information dic because it is sequence specific information
            primer_dic["sequence_information"]["SEQUENCE_LEFT_NUM_RETURNED"] = value
            num_left = value
            # create empty dictionaries with primer name keys
            for i in range(int(value)):
                primer_key = "PRIMER_LEFT_" + str(i)
                primer_dic["primer_information"][primer_key] = {}
    # do the same for right primers found
    for pair in lines:
        tag = pair[0]
        value = pair [1]
        if tag == "PRIMER_RIGHT_NUM_RETURNED":
            primer_dic["sequence_information"]["SEQUENCE_RIGHT_NUM_RETURNED"] = value
            num_right = value
            for i in range(int(value)):
                primer_key = "PRIMER_RIGHT_" + str(i)
                primer_dic["primer_information"][primer_key] = {}
    # get sequence coordinate information to determine genomic coordinates of primers
    # because primer information is relative to template sequence
    sequence_coordinates = get_coordinates(primer_dic["sequence_information"]["SEQUENCE_ID"])
    seq_chr = sequence_coordinates[0]
    """
    if seq_chr == "chrx":
        seq_chr = "chrX"
    elif seq_chr == "chry":
        seq_chr = "chrY"
    """
    seq_start = int(sequence_coordinates[1])
    seq_end = int(sequence_coordinates[2])
    # get primer information from input file and add to primer dictionary
    for pair in lines:
        tag = pair[0]
        value = pair [1]
        if ((tag.startswith("PRIMER_LEFT_") or tag.startswith("PRIMER_RIGHT_")) and             (tag != "PRIMER_LEFT_NUM_RETURNED") and (tag != "PRIMER_RIGHT_NUM_RETURNED")):
            attributes = tag.split('_')
            # primer coordinates tag does not include an attribute value
            # it is only primer name = coordinates, so:
            if len(attributes)>3:
                # then this attribute is not coordinates and should have an attribute value
                # such as TM or HAIRPIN etc.
                primer_name = '_'.join(attributes[0:3])
                attribute_value = '_'.join(attributes[3:])
                primer_dic["primer_information"][primer_name][attribute_value] = value
            else:
                # then this attribute is coordinates and has no attribute value
                # give it an attribute valute "COORDINATES"
                primer_name = '_'.join(attributes[0:3])
                primer_dic["primer_information"][primer_name]['COORDINATES'] = value
                # the coordinates are relative to sequence template
                # find the genomic coordinates
                coordinate_values = value.split(",")
                if tag.startswith("PRIMER_LEFT"):
                    # sequence start is added to primer start to get genomic primer start
                    genomic_start = seq_start + int(coordinate_values[0])
                    # primer len is added "to genomic start because it is a left primer
                    genomic_end = genomic_start + int(coordinate_values[1]) - 1
                    primer_dic["primer_information"][primer_name]['GENOMIC_START'] = genomic_start
                    primer_dic["primer_information"][primer_name]['GENOMIC_END'] = genomic_end
                    primer_dic["primer_information"][primer_name]['CHR'] = seq_chr
                    primer_dic["primer_information"][primer_name]['ORI'] = "forward"
                else:
                    # sequence start is added to primer start to get genomic primer start
                    genomic_start = seq_start + int(coordinate_values[0])
                    # primer len is subtracted from genomic start because it is a right primer
                    genomic_end = genomic_start - int(coordinate_values[1]) + 1
                    primer_dic["primer_information"][primer_name]['GENOMIC_START'] = genomic_start
                    primer_dic["primer_information"][primer_name]['GENOMIC_END'] = genomic_end
                    primer_dic["primer_information"][primer_name]['CHR'] = seq_chr
                    primer_dic["primer_information"][primer_name]['ORI'] = "reverse"
            # add NAME as a key to primer information dictionary
            primer_dic["primer_information"][primer_name]['NAME'] = primer_name
    # if some primers were eliminated from initial primer3 output, remove from dictionary
    for primer in list(primer_dic["primer_information"].keys()):
        if primer_dic["primer_information"][primer] == {}:
            primer_dic["primer_information"].pop(primer)
    # dump the dictionary to json file in primer3_output_DIR if outp parameter is true
    if outp:
        dict_file = open (primer3_output_DIR + parse_out, 'w')
        json.dump(primer_dic, dict_file,indent=1)
        dict_file.close()
    # generate a simple fasta file with primer names
    if fasta:
        outfile = open(bowtie2_input_DIR+parse_out, 'w')
        for primer in primer_dic["primer_information"]:
            # primer name is fasta header and sequence is fasta sequence
            fasta_head = primer
            fasta_line = primer_dic["primer_information"][primer]["SEQUENCE"]
            outfile.write(">" + fasta_head + "\n" + fasta_line +"\n")
        outfile.close()
    return primer_dic
def adjust_worker(chores):
    try:
        p_name = chores[0]
        p_dic = chores[1]
        chroms = chores[2]
        p_copies = chores[3]
        p_tm_diff = chores[4]
        p_end_identity = chores[5]
        p_temp_dir = chores[6]
        p_settings = chores[7]
        p_spec = chores[8]
        Na = float(p_settings["Na"])
        Mg = float(p_settings["Mg"])
        conc = float(p_settings["oligo_conc"])
        alt_arm = int(p_settings["alternative_arms"])
        #chroms = p_coord["C0"]["chromosomes"]
        ref_coord = p_dic["COORDINATES"]
        primer_ori = p_dic["ORI"]
        paralogs = p_dic["PARALOG_COORDINATES"]
        primer_seq = p_dic["SEQUENCE"]
        p_dic["BOWTIE_BINDS"] = []
        try:
            start = paralogs["C0"]["BOWTIE_START"]
            end = paralogs["C0"]["BOWTIE_END"]
        except KeyError:
            return [p_name, p_dic]

        # add reference copy as paralog
        paralogs["C0"]["BOWTIE_BOUND"] = True
        paralogs["C0"]["BOWTIE_SEQUENCE"] = primer_seq
        p_dic["BOWTIE_BINDS"] = ["C0"]
        #if alt_arm:
        paralogs["C0"]["ALT_BOUND"] = True
        paralogs["C0"]["ALT_SEQUENCE"] = primer_seq
        p_dic["ALT_BINDS"] = []

        for i in range(len(p_copies)-1):
            para = paralogs[p_copies[i+1]]
            try:
                para_primer_ori = para["ORI"]
                para_start = para["BOWTIE_START"]
                para_end = para["BOWTIE_END"]
                para_chr = para["CHR"]
                if para_primer_ori == "forward":
                    para_primer_key = para_chr + ":" + str(para_start) + "-" + str(para_end)
                    para_primer_seq = fasta_to_sequence(get_fasta(para_primer_key, p_spec))
                else:
                    para_primer_key = para_chr + ":" + str(para_end) + "-" + str(para_start)
                    para_primer_seq = reverse_complement(fasta_to_sequence(get_fasta(para_primer_key, p_spec)))
                para["BOWTIE_SEQUENCE"] = para_primer_seq
                if para_primer_seq.upper() == primer_seq.upper():
                    p_dic["BOWTIE_BINDS"].append(p_copies[i+1])
                    para["BOWTIE_BOUND"] = True
                else:
                    para["BOWTIE_BOUND"] = False
                    if alt_arm:
                        try:
                            alt_tms = {}
                            ref_tm = hybrid_TM(p_temp_dir, reverse_complement(primer_seq), primer_seq, Na=Na, Mg=Mg, conc=conc)
                            for j in range(-3,4):
                                alt_start = para_start + j
                                if para_primer_ori == "forward":
                                    alt_primer_key = para_chr + ":" + str(alt_start) + "-" + str(para_end)
                                    alt_primer_seq = fasta_to_sequence(get_fasta(alt_primer_key, p_spec))
                                else:
                                    alt_primer_key = para_chr + ":" + str(para_end) + "-" + str(alt_start)
                                    alt_primer_seq = reverse_complement(fasta_to_sequence(
                                                                        get_fasta(alt_primer_key, p_spec)))
                                alt_tm = hybrid_TM(p_temp_dir, alt_primer_seq, reverse_complement(alt_primer_seq),
                                                   Na=Na, Mg=Mg, conc=conc)
                                alt_tms[abs(alt_tm - ref_tm)] = [alt_tm, alt_start, alt_primer_seq]
                            best_alt_tm = min(alt_tms.keys())
                            para["ALT_START"] = alt_tms[best_alt_tm][1]
                            para["ALT_TM"] = alt_tms[best_alt_tm][0]
                            para["ALT_DIFF"] = ref_tm - para["ALT_TM"]
                            para["ALT_SEQUENCE"] = alt_tms[best_alt_tm][2]
                            if best_alt_tm <= p_tm_diff:
                                para["ALT_BOUND"] = True
                                p_dic["ALT_BINDS"].append(p_copies[i+1])
                            else:
                                para["ALT_BOUND"] = False
                        except Exception as e:
                            p_dic["BOWTIE_BINDS"] = [str(e)]
                    else:
                        para["ALT_TM"] = 0
                        para["ALT_TM_DIFF"] = 100
                        para["ALT_BOUND"] = False
            except KeyError:
                para["BOWTIE_TM"] = 0
                para["BOWTIE_TM_DIFF"] = 100
                para["BOWTIE_BOUND"] = False
                para["ALT_TM"] = 0
                para["ALT_TM_DIFF"] = 100
                para["ALT_BOUND"] = False

        return [p_name, p_dic]
    except Exception as e:
        p_dic["error"] = str(e)
        return [p_name, p_dic]
def adjust_paralogs (primer_dic, copies, p_chromosomes, settings,                           primer3_output_DIR, outname, species):
    """ Take a primer dictionary file and add genomic start and end coordinates
    of all its paralogs."""
    # uncomment for using json object instead of dic
    # load the primers dictionary from file
    # with open(primer_file, "r") as infile:
    #     primer_dic = json.load(infile)
    # primer dict consists of 2 parts, sequence_information dict
    # and primer information dict. We wont'change the sequence_info part
    primers = dict(primer_dic["primer_information"])
    # get chromosomes of the paralog copies
    #chromosomes = coordinate_converter["C0"]["chromosomes"]
    # allowed temperature difference for a primer to be considered binding
    # to a paralog. For example, if its tm for ref copy is 60°C, and for C1 59,
    # it will be considered satisfying tm restriction, if tm_diff is 2, but
    # not satisfying if tm_diff is 0.
    tm_diff = float(settings["tm_diff"])
    # A specified number of nucleotides at the 3' of a primer and its paralog primer
    # need to be identical for paralog primers to be considered valid.
    end_identity = -1 * int(settings["3p_identity"])
    # add paralog coordinate information for each primer in primers dic
    chore_list = []
    dics = []
    # create a temp directory in temp_dir
    tmp = primer3_output_DIR + "tmp/"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    for primer in list(primers.keys()):
        chore_list.append([primer, primers[primer], p_chromosomes,                          copies, tm_diff, end_identity, tmp, settings, species])
    num_process = int(settings["processors"])
    #print chore_list[-1]
    p = Pool(num_process)
    p.map_async(adjust_worker, chore_list, callback=dics.extend)
    p.close()
    p.join()
    #print dics[:2]
    #print dics
    for dic in dics:
        primers[dic[0]] = dic[1]
    out_dic = {}
    out_dic["sequence_information"] = primer_dic["sequence_information"]
    out_dic["primer_information"] = primers
    with open(primer3_output_DIR + outname, "w") as outf:
        json.dump(out_dic, outf, indent=1)
    shutil.rmtree(tmp)
    return out_dic
def paralog_primer_worker(chores):
    p_name = chores[0]
    p_dic = chores[1]
    p_coord = chores[2]
    p_copies = chores[3]
    p_tm_diff = chores[4]
    p_end_identity = chores[5]
    p_temp_dir = chores[6]
    p_settings = chores[7]
    p_spec = chores[8]
    Na = float(p_settings["Na"])
    Mg = float(p_settings["Mg"])
    conc = float(p_settings["oligo_conc"])
    chroms = p_coord["C0"]["chromosomes"]
    start = p_dic["GENOMIC_START"]
    end = p_dic["GENOMIC_END"]
    ref_coord = p_dic["COORDINATES"]
    primer_ori = p_dic["ORI"]
    #para_start = p_coord["C0"][start]
    #para_end = p_coord["C0"][end]
    p_dic["PARALOG_COORDINATES"] = {}
    primer_seq = p_dic["SEQUENCE"]
    # add reference copy as paralog
    #ref_tm = hybrid_TM(p_temp_dir, extended_primer_seq, primer_seq, Na=Na, Mg=Mg, conc=conc)
    p_dic["PARALOG_COORDINATES"]["C0"] = {"SEQUENCE": primer_seq,                      "ORI": primer_ori,"CHR":chroms["C0"], "NAME":p_name,                      "GENOMIC_START": start, "GENOMIC_END": end, "COORDINATES":ref_coord,                      }
    for c in p_copies:
        if c!= "C0":
            # check if both ends of the primer has aligned with the reference
            try:
                para_start = p_coord["C0"][c][start]
                para_end = p_coord["C0"][c][end]
            except KeyError:
                # do not add that copy if it is not aligned
                continue
            para_primer_ori = para_start < para_end
            if para_primer_ori:
                para_primer_key = chroms[c] + ":" + str(para_start) + "-" + str(para_end)
                para_primer_seq = fasta_to_sequence(get_fasta(para_primer_key, p_spec))
                p_dic["PARALOG_COORDINATES"][c] = {"SEQUENCE": para_primer_seq,                          "ORI": "forward","CHR":chroms[c], "NAME":p_name,                          "GENOMIC_START": para_start, "GENOMIC_END": para_end,                          "COORDINATES":ref_coord}
            else:
                para_primer_key = chroms[c] + ":" + str(para_end) + "-" + str(para_start)
                para_primer_seq = reverse_complement(                                  fasta_to_sequence(get_fasta(para_primer_key, p_spec)))
                p_dic["PARALOG_COORDINATES"][c] = {"SEQUENCE": para_primer_seq,                          "ORI": "reverse","CHR":chroms[c], "NAME":p_name,                          "GENOMIC_START": para_start, "GENOMIC_END": para_end,                          "COORDINATES":ref_coord}

    return [p_name, p_dic]
def paralog_primers_multi (primer_dic, copies, coordinate_converter, settings,                           primer3_output_DIR, outname, species, outp = 0):
    """ Take a primer dictionary file and add genomic start and end coordinates
    of all its paralogs."""
    # uncomment for using json object instead of dic
    # load the primers dictionary from file
    # with open(primer_file, "r") as infile:
    #     primer_dic = json.load(infile)
    # primer dict consists of 2 parts, sequence_information dict
    # and primer information dict. We wont'change the sequence_info part
    primers = copy.deepcopy(primer_dic["primer_information"])
    # get chromosomes of the paralog copies
    chromosomes = coordinate_converter["C0"]["chromosomes"]
    # allowed temperature difference for a primer to be considered binding
    # to a paralog. For example, if its tm for ref copy is 60°C, and for C1 59,
    # it will be considered satisfying tm restriction, if tm_diff is 2, but
    # not satisfying if tm_diff is 0.
    tm_diff = float(settings["tm_diff"])
    # A specified number of nucleotides at the 3' of a primer and its paralog primer
    # need to be identical for paralog primers to be considered valid.
    end_identity = -1 * int(settings["3p_identity"])
    # add paralog coordinate information for each primer in primers dic
    chore_list = []
    dics = []
    # create a temp directory in temp_dir
    tmp = primer3_output_DIR + "tmp/"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    for primer in list(primers.keys()):
        chore_list.append([primer, primers[primer], coordinate_converter,                          copies, tm_diff, end_identity, tmp, settings, species])
    num_process = int(settings["processors"])
    p = Pool(num_process)
    p.map_async(paralog_primer_worker, chore_list, callback=dics.extend)
    p.close()
    p.join()
    #print dics
    for dic in dics:
        primers[dic[0]] = dic[1]
    out_dic = {}
    out_dic["sequence_information"] = primer_dic["sequence_information"]
    out_dic["primer_information"] = primers
    if outp:
        with open(primer3_output_DIR + outname, "w") as outf:
            json.dump(out_dic, outf, indent=1)
    shutil.rmtree(tmp)
    return out_dic
def alter_worker(chores):
    p_name = chores[0]
    p_dic = chores[1]
    p_coord = chores[2]
    p_copies = chores[3]
    p_tm_diff = chores[4]
    p_end_identity = chores[5]
    p_temp_dir = chores[6]
    p_settings = chores[7]
    p_spec = chores[8]
    Na = float(p_settings["Na"])
    Mg = float(p_settings["Mg"])
    conc = float(p_settings["oligo_conc"])
    chroms = p_coord["C0"]["chromosomes"]
    ref_coord = p_dic["COORDINATES"]
    primer_ori = p_dic["ORI"]
    paralogs = p_dic["PARALOG_COORDINATES"]
    primer_seq = p_dic["SEQUENCE"]
    try:
        start = paralogs["C0"]["BOWTIE_START"]
        end = paralogs["C0"]["BOWTIE_END"]
    except KeyError:
        return [p_name, p_dic]
    # add reference copy as paralog
    ref_tm = hybrid_TM(p_temp_dir, primer_seq, reverse_complement(primer_seq), Na=Na, Mg=Mg, conc=conc)
    paralogs["C0"]["ALT_TM"] = ref_tm
    paralogs["C0"]["ALT_TM_DIFF"] = ref_tm - ref_tm
    paralogs["C0"]["ALT_BOUND"] = True
    p_dic["ALT_BINDS"] = ["C0"]
    for i in range(1, len(p_copies)):
        p_cop = p_copies[i]
        try:
            para = paralogs[p_cop]
        except KeyError:
            continue
        try:
            para_primer_ori = para["ORI"]
            para_start = para["BOWTIE_START"]
            para_end = para["BOWTIE_END"]
            para_chr = chroms[p_cop]
            alt_tms = {}
            for i in range(-3,4):
                alt_start = para_start + i
                if para_primer_ori == "forward":
                    alt_primer_key = para_chr + ":" + str(alt_start) + "-" + str(para_end)
                    alt_primer_seq = fasta_to_sequence(get_fasta(alt_primer_key, p_spec))
                else:
                    alt_primer_key = para_chr + ":" + str(para_end) + "-" + str(alt_start)
                    alt_primer_seq = reverse_complement(fasta_to_sequence(
                                                        get_fasta(alt_primer_key, p_spec)))
                alt_tm = hybrid_TM(p_temp_dir, alt_primer_seq, reverse_complement(alt_primer_seq),
                                   Na=Na, Mg=Mg, conc=conc)
                alt_tms[abs(alt_tm - ref_tm)] = [alt_tm, alt_start, alt_primer_seq]
            best_alt_tm = min(alt_tms.keys())
            para["ALT_START"] = alt_tms[best_alt_tm][1]
            para["ALT_TM"] = alt_tms[best_alt_tm][0]
            para["ALT_DIFF"] = ref_tm - para["ALT_TM"]
            para["ALT_SEQUENCE"] = alt_tms[best_alt_tm][2]
            if best_alt_tm <= p_tm_diff:
                para["ALT_BOUND"] = True
                p_dic["ALT_BINDS"].append(p_cop)

        except KeyError:
            para["ALT_TM"] = 0
            para["ALT_TM_DIFF"] = ref_tm
            para["ALT_BOUND"] = False

    return [p_name, p_dic]
def alternative_paralogs (primer_dic, copies, coordinate_converter, settings,                           primer3_output_DIR, outname, species):
    """ Take a primer dictionary file and add genomic start and end coordinates
    of all its paralogs."""
    # uncomment for using json object instead of dic
    # load the primers dictionary from file
    # with open(primer_file, "r") as infile:
    #     primer_dic = json.load(infile)
    # primer dict consists of 2 parts, sequence_information dict
    # and primer information dict. We wont'change the sequence_info part
    primers = dict(primer_dic["primer_information"])
    # get chromosomes of the paralog copies
    chromosomes = coordinate_converter["C0"]["chromosomes"]
    # allowed temperature difference for a primer to be considered binding
    # to a paralog. For example, if its tm for ref copy is 60°C, and for C1 59,
    # it will be considered satisfying tm restriction, if tm_diff is 2, but
    # not satisfying if tm_diff is 0.
    tm_diff = float(settings["tm_diff"])
    # A specified number of nucleotides at the 3' of a primer and its paralog primer
    # need to be identical for paralog primers to be considered valid.
    end_identity = -1 * int(settings["3p_identity"])
    # add paralog coordinate information for each primer in primers dic
    chore_list = []
    dics = []
    # create a temp directory in temp_dir
    tmp = primer3_output_DIR + "tmp/"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    for primer in list(primers.keys()):
        chore_list.append([primer, primers[primer], coordinate_converter,                          copies, tm_diff, end_identity, tmp, settings, species])
    num_process = int(settings["processors"])
    #print chore_list[-1]
    p = Pool(num_process)
    p.map_async(alter_worker, chore_list, callback=dics.extend)
    p.close()
    p.join()
    #print dics[:2]
    #print dics
    for dic in dics:
        primers[dic[0]] = dic[1]
    out_dic = {}
    out_dic["sequence_information"] = primer_dic["sequence_information"]
    out_dic["primer_information"] = primers
    with open(primer3_output_DIR + outname, "w") as outf:
        json.dump(out_dic, outf, indent=1)
    shutil.rmtree(tmp)
    return out_dic
def fasta_parser(fasta):
    """ Convert a fasta file with multiple sequences
    to a dictionary with fasta headers as keys and sequences
    as values."""
    fasta_dic = {}
    with open(fasta) as infile:
        for line in infile:
            # find the headers
            if line.startswith(">"):
                header = line[1:-1].split(" ")[0]
                if header in fasta_dic:
                    print(("%s occurs multiple times in fasta file" %header))
                fasta_dic[header] = ""
                continue
            try:
                fasta_dic[header] = fasta_dic[header] + line.strip()
            except KeyError:
                fasta_dic[header] = line.strip()
    return fasta_dic
def fasta_to_sequence(fasta):
    """ Convert a fasta sequence to one line sequence"""
    f = fasta.strip().split("\n")
    if len(f) > 0:
        return "".join(f[1:])
    else:
        return ""
def get_sequence(region, species):
    return fasta_to_sequence(get_fasta(region,species))
def bowtie2_run (fasta_file, output_file, bowtie2_input_DIR, bowtie2_output_DIR, species,process_num=4,                 seed_MM=1, mode="-a", seed_len=18, gbar=1, local=0):
    """ Extract primer sequences from the fasta file,
    check alignments for given species genome(s), create
    sam output file(s). Species must be a list!"""
    file_locations = get_file_locations()
    # check if entered species is supported
    genome = file_locations[species]["bowtie2_genome"]
    # determine what type of alignment is wanted
    # local or end-to-end
    if local:
        check_local = "--local"
    else:
        check_local = "--end-to-end"
    dump = subprocess.check_output(["bowtie2",
                                    "-p",
                                    str(process_num),
                                   "-D",
                                    "20",
                                    "-R",
                                    "3",
                                    "-N",
                                    str(seed_MM),
                                    "-L",
                                    str(seed_len),
                                    "-i",
                                    "S,1,0.5",
                                    "--gbar",
                                    str(gbar),
                                    mode,
                                    check_local,
                                   "-x",
                                    genome,
                                    "-f",
                                    bowtie2_input_DIR + fasta_file,
                                   "-S",
                                    bowtie2_output_DIR + output_file])
    return 0
def bowtie(fasta_file, output_file, bowtie2_input_DIR, bowtie2_output_DIR, options,
           species,process_num=4, mode="-a", local=0, fastq = 0):
    """ Extract primer sequences from the fasta file,
    check alignments for given species genome(s), create
    sam output file(s). Species must be a list!"""
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
        com.append("-q " + bowtie2_input_DIR + fasta_file)
    else:
        com.append("-f " + bowtie2_input_DIR + fasta_file)
    com.append("-S " + bowtie2_output_DIR + output_file)
    dump = subprocess.check_output(com)
    return 0
def bwa(fastq_file, output_file, output_type, input_dir,
        output_dir, options, species):
    """ Run bwa alignment on given fastq file using the species bwa indexed genome.
    options should be a list that starts with the command (e.g. mem, aln etc).
    Additional options should be appended as strings of "option value",
    for example, "-t 30" to use 30 threads. Output type can be sam or bam.
    Recommended options ["-t30", "-L500", "-T100"]. Here L500 penalizes clipping
    severely so the alignment becomes end-to-end and T100 stops reporting secondary
    alignments, assuming their score is below 100."""
    genome_file = get_file_locations()[species]["bwa_genome"]
    if output_type == "sam":
        com = ["bwa"]
        com.extend(options)
        com.append(genome_file)
        com.append(input_dir + fastq_file)
        with open(output_dir + output_file, "w") as outfile:
            dump = subprocess.check_call(com, stdout=outfile)
    else:
        com = ["bwa"]
        com.extend(options)
        com.append(genome_file)
        com.append(input_dir + fastq_file)
        sam = subprocess.Popen(com, stdout=subprocess.PIPE)
        bam_com = ["samtools", "view", "-b"]
        with open(output_dir + output_file, "w") as outfile:
            bam = subprocess.Popen(bam_com, stdin=sam.stdout,
                                   stdout=outfile)
def reverse_complement(sequence):
    """ Return reverse complement of a sequence. """
    complement_bases = {
        'g':'c', 'c':'g', 'a':'t', 't':'a', 'n':'n',
        'G':'C', 'C':'G', 'A':'T', 'T':'A', 'N':'N', "-":"-",
        "R":"Y", "Y":"R", "S":"W", "W":"S", "K":"M", "M":"K",
        "B":"V", "V":"B", "D": "H", "H": "D",
        "r":"y", "y":"r", "s":"w", "w":"s", "k":"m", "m":"k",
        "b":"v", "v":"b", "d": "h", "h": "d"
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
def check_TM (s1,s2):
    """ Return melting temparature of two sequences as string.
    Uses ntthal program and current MIP conditions, salt concentration etc.
    ntthal uses Santa Lucia salt correction which yields TMs a few degrees
    higher than the most current R.O. correction. Since it is used for
    determining mispriming, the values it generates are on the safer side.
    For example, a 46C cut off eliminates primers with 43-44 C mispriming.
    """
    # split ntthal command line input to a list of strings
    # to pass into subprocess module
    ntthal_input = ["ntthal", "-r", "-mv", "25", "-dv", "10", "-n", "0.02",                    "-d", "0.4",  "-path",                    "/home/aydemiro/programs/primer3-2.3.6/src/primer3_config/",                    "-s1", s1, "-s2", s2]
    # run ntthal through subprocess
    TM = subprocess.check_output(ntthal_input, stderr=subprocess.STDOUT)
    return TM
def parse_cigar(cigar):
    """ Parse a cigar string which is made up of numbers followed
    by key letters that represent a sequence alignment; return a dictionary
    with alignment keys and number of bases with that alignment key as values.
    Below is some more information about cigar strings.

    2S20M1I2M5D,for, example would mean that the 2 bases are "S"oft clipped from
    5' end of the sequence(read) aligned and it is not part of the alignment;
    following that 2 bases, 20 bases of the read aligns or "Mathces" to the reference sequence,
    mathc here does not mean the bases are identical, just that there is 1 base of
    reference for each base of the read and there are enough similarity between the two
    sequences that they aligned. 1 base following the 20M is an insertion, that is,
    it exists in the read but not in the reference; 5 bases at the end are "D"eletions,
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
    """ Get the length of the reference sequence that a read is aligned to,
    given their cigar string."""
    try:
        # parse cigar string and find out how many insertions are in the alignment
        insertions = parse_cigar(cigar)["I"]
    except KeyError:
        # the key "I" will not be present in the cigar string if there is no insertion
        insertions = 0
    # all the values in the cigar dictionary represent a base in the reference seq,
    # except the insertions, so they should be subtracted
    return sum(parse_cigar(cigar).values()) - insertions
def cleanest_bowtie (primer_file, primer_out,primer3_output_DIR,                   bowtie2_output_DIR, species, settings, host = False,
                    outp = 1):
    """ Take a primer dict with bowtie information added.
    Look at bowtie hits for each primer, determine if they
    are on intended targets or nonspecific. In cases of paralogus
    regions, check all paralogs and determine if the primer
    will bind to any paralog. Create alternative primers if necessary
    and allowed. Get melting temperatures of all hits and add
    all these information to the primer dictionary.
    """
    Na = float(settings["Na"])
    Mg = float(settings["Mg"])
    conc = float(settings["oligo_conc"])
    alt_arm = int(settings["alternative_arms"])
    """
    # read in primer/mip file
    with open (primer3_output_DIR + primer_file, 'r') as handle:
        primers = json.load(handle)

    """
    primers = primer_file
    seq_list = []
    hit_list = []
    bowtie_key ="bowtie_information_" + species
    # read bowtie hits
    for primer_name in primers['primer_information']:
        try:
            primer_seq = primers['primer_information'][primer_name]["SEQUENCE"]
            if not host:
                para = primers['primer_information'][primer_name]["PARALOG_COORDINATES"]
                if "BOWTIE_BINDS" not in primers['primer_information'][primer_name]:
                    primers['primer_information'][primer_name]["BOWTIE_BINDS"] = []
                if "ALT_BINDS" not in primers['primer_information'][primer_name]:
                    primers['primer_information'][primer_name]["ALT_BINDS"] = []
            for bt_hit_name in list(primers['primer_information'][primer_name][bowtie_key].keys()):
                bt_hit = primers['primer_information'][primer_name][bowtie_key][bt_hit_name]
                bt_chrom = bt_hit["chrom"]
                bt_begin = bt_hit["begin"]
                bt_end = bt_hit["end"]
                bt_ori = bt_hit["strand"]
                bt_seq = bt_hit["sequence"]
                if host:
                    seq_list.extend([">" + primer_name, primer_seq])
                    hit_list.extend([">bt_" + bt_hit_name, reverse_complement(bt_seq)])
                    continue
                intended = 0
                # para is a dict like {C0:{"CHR": "chr4", "GENOMIC_START" ..}, C1:{..
                # for non-CNV regions, bowtie mapping should be exactly the same as
                # genomic coordinates, so even if there is 1 bp difference, we'll count
                # this as off target. For CNV regions, a more generous 100 bp
                # padding will be allowed to account for differences in our mapping
                # and bowtie mapping
                map_padding = 1
                if len(para) > 1:
                    map_padding = 100
                for k in para:
                    para_ori = para[k]["ORI"]
                    para_chr = para[k]["CHR"]
                    para_begin = para[k]["GENOMIC_START"]
                    para_end = para[k]["GENOMIC_END"]
                    if (para_ori == bt_ori) and (para_chr == bt_chrom) and                        (abs(para_begin - bt_begin)
                        < map_padding) and (abs(para_end - bt_end) < map_padding):
                        intended = 1
                        para[k]["BOWTIE_END"] = bt_end
                        para[k]["BOWTIE_START"] = bt_begin
                        para[k]["BOWTIE_SEQUENCE"] = bt_seq
                    if intended:
                        if bt_seq.upper() == primer_seq.upper():
                            para[k]["BOWTIE_BOUND"] = True
                            primers['primer_information'][primer_name]["BOWTIE_BINDS"].append(k)
                        else:
                            para[k]["BOWTIE_BOUND"] = False
                            if alt_arm:
                                al = {}
                                al ["ref"] = {"ALT_SEQUENCE": primer_seq}
                                seq_list.extend([">" + primer_name, primer_seq])
                                hit_list.extend([">alt_" + k + "_ref", reverse_complement(
                                                 primer_seq)])
                                for j in range(-3,4):
                                    if j == 0:
                                        continue
                                    alt_start = bt_begin + j
                                    alt_end = bt_end
                                    if para_ori == "forward":
                                        alt_primer_key = create_region(bt_chrom,
                                                                       alt_start,
                                                                       alt_end)
                                        alt_primer_seq = get_sequence(alt_primer_key,
                                                                     species)
                                    else:
                                        alt_primer_key = create_region(bt_chrom,
                                                                       alt_end,
                                                                       alt_start)
                                        alt_primer_seq = get_sequence(alt_primer_key,
                                                                     species)
                                        alt_primer_seq = reverse_complement(alt_primer_seq)
                                    if alt_primer_seq == "":
                                        continue
                                    al[j] = {}
                                    al[j]["ALT_START"] = alt_start
                                    al[j]["ALT_END"] = alt_end
                                    al[j]["ALT_SEQUENCE"] = alt_primer_seq
                                    seq_list.extend([">" + primer_name, alt_primer_seq])
                                    hit_list.extend([">alt_" + k + "_" + str(j),
                                                    reverse_complement(alt_primer_seq)])
                                para[k]["ALTERNATIVES"] = al
                            else:
                                para[k]["ALTERNATIVES"] = {}
                                para[k]["ALT_TM"] = 0
                                para[k]["ALT_TM_DIFF"] = 100
                                para[k]["ALT_BOUND"] = False
                        primers['primer_information'][primer_name][bowtie_key].pop(bt_hit_name)
                        break
                # if the target is intended it should not be added to the cleaned bowtie hits
                # which is for nonspecific targets only
                if not intended:
                    seq_list.extend([">" + primer_name, primer_seq])
                    hit_list.extend([">bt_" + bt_hit_name, reverse_complement(bt_seq)])
                    #print bt_seq
            # end bowtie check
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
                            al = {}
                            al ["ref"] = {"ALT_SEQUENCE": primer_seq}
                            seq_list.extend([">" + primer_name, primer_seq])
                            hit_list.extend([">alt_" + k + "_ref", reverse_complement(
                                             primer_seq)])
                            for j in range(-3,4):
                                if j == 0:
                                    continue
                                alt_start = para_begin + j
                                alt_end = para_end
                                if para_ori == "forward":
                                    alt_primer_key = create_region(para_chr,
                                                                   alt_start,
                                                                   alt_end)
                                    alt_primer_seq = get_sequence(alt_primer_key,
                                                                 species)
                                else:
                                    alt_primer_key = create_region(para_chr,
                                                                   alt_end,
                                                                   alt_start)
                                    alt_primer_seq = get_sequence(alt_primer_key,
                                                                 species)
                                    alt_primer_seq = reverse_complement(alt_primer_seq)
                                if alt_primer_seq == "":
                                    continue
                                al[j] = {}
                                al[j]["ALT_START"] = alt_start
                                al[j]["ALT_END"] = alt_end
                                al[j]["ALT_SEQUENCE"] = alt_primer_seq
                                seq_list.extend([">" + primer_name, alt_primer_seq])
                                hit_list.extend([">alt_" + k + "_" + str(j),
                                                reverse_complement(alt_primer_seq)])
                            para[k]["ALTERNATIVES"] = al
                        else:
                            para[k]["ALTERNATIVES"] = {}
                            para[k]["ALT_TM"] = 0
                            para[k]["ALT_TM_DIFF"] = 100
                            para[k]["ALT_BOUND"] = False

        except KeyError as e:
            #print str(e)

            continue



    #infile.close()
    #outfile.close()
    if len(hit_list) > 0:
        with open(bowtie2_output_DIR + "seq_for_tm", "w") as seq_for_tm, \
            open(bowtie2_output_DIR + "hits_for_tm", "w") as hits_for_tm:
            seq_for_tm.write("\n".join(seq_list))
            hits_for_tm.write("\n".join(hit_list))
        melt_ext = id_generator(6)
        tms = hybrid_TM_files(bowtie2_output_DIR,
                              "seq_for_tm",
                              "hits_for_tm",
                              "melting_" + melt_ext,
                              Na=Na,
                              Mg=Mg,
                              conc=conc)
        tm_list = tms.split("\n")
        tm_count = (len(tm_list) - 1)//2
        for i in range(tm_count):
            items = tm_list[i].split(" ")
            h_tm = float(tm_list[i+tm_count+1].split("\t")[-1])
            p_name = items[2]
            h_type = items[4].split("_")[0]
            h_name = items[4].split("_")[-1][:-1]
            if host:
                bowtie_dic = primers["primer_information"][p_name][bowtie_key]
                for b in bowtie_dic:
                    if h_name == b:
                        bowtie_dic[b]["TM"] = h_tm
            elif h_type == "alt":
                copyname = items[4].split("_")[1]
                try:
                    if h_name == "ref":
                        primers["primer_information"][p_name]["PARALOG_COORDINATES"][copyname]                                   ["ALTERNATIVES"][h_name]["ALT_TM"] = h_tm
                    else:
                        primers["primer_information"][p_name]["PARALOG_COORDINATES"][copyname]                                       ["ALTERNATIVES"][int(h_name)]["ALT_TM"] = h_tm
                except KeyError:
                    pass

            else:
                bowtie_dic = primers["primer_information"][p_name][bowtie_key]
                for b in bowtie_dic:
                    if h_name == b:
                        bowtie_dic[b]["TM"] = h_tm
    # write the (if) updated dictionary to the file
    if outp:
        with open(primer3_output_DIR + primer_out, 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers
def cleanest_bowtie_old (primer_file, primer_out,primer3_output_DIR,                   bowtie2_output_DIR, species, settings):
    """ Take a primer dict with bowtie information added.
    Look at bowtie hits for each primer, determine if they
    are on intended targets or nonspecific. In cases of paralogus
    regions, check all paralogs and determine if the primer
    will bind to any paralog. Create alternative primers if necessary
    and allowed. Get melting temperatures of all hits and add
    all these information to the primer dictionary.
    """
    Na = float(settings["Na"])
    Mg = float(settings["Mg"])
    conc = float(settings["oligo_conc"])
    alt_arm = int(settings["alternative_arms"])
    # read in primer/mip file
    with open (primer3_output_DIR + primer_file, 'r') as handle:
        primers = json.load(handle)
    seq_list = []
    hit_list = []
    bowtie_key ="bowtie_information_" + species
    # read bowtie hits
    for primer_name in primers['primer_information']:
        try:
            primer_seq = primers['primer_information'][primer_name]["SEQUENCE"]
            for bt_hit_name in list(primers['primer_information'][primer_name][bowtie_key].keys()):
                bt_hit = primers['primer_information'][primer_name][bowtie_key][bt_hit_name]
                bt_chrom = bt_hit["chrom"]
                bt_begin = bt_hit["begin"]
                bt_end = bt_hit["end"]
                bt_ori = bt_hit["strand"]
                bt_seq = bt_hit["sequence"][3:]
                intended = 0
                #try:
                para = primers['primer_information'][primer_name]["PARALOG_COORDINATES"]
                if "BOWTIE_BINDS" not in primers['primer_information'][primer_name]:
                    primers['primer_information'][primer_name]["BOWTIE_BINDS"] = []
                if "ALT_BINDS" not in primers['primer_information'][primer_name]:
                    primers['primer_information'][primer_name]["ALT_BINDS"] = []
                # para is a dict like {C0:{"CHR": "chr4", "GENOMIC_START" ..}, C1:{..
                for k in para:
                    para_ori = para[k]["ORI"]
                    para_chr = para[k]["CHR"]
                    para_begin = para[k]["GENOMIC_START"]
                    para_end = para[k]["GENOMIC_END"]
                    if (para_ori == bt_ori) and (para_chr == bt_chrom) and                        (abs(para_begin - bt_begin) < 10) and (abs(para_end - bt_end) < 10):
                        intended = 1
                        para[k]["BOWTIE_END"] = bt_end
                        para[k]["BOWTIE_START"] = bt_begin
                        para[k]["BOWTIE_SEQUENCE"] = bt_seq
                    if intended:
                        if bt_seq.upper() == primer_seq.upper():
                            para[k]["BOWTIE_BOUND"] = True
                            primers['primer_information'][primer_name]["BOWTIE_BINDS"].append(k)
                        else:
                            para[k]["BOWTIE_BOUND"] = False
                            if alt_arm:
                                al = {}
                                al ["ref"] = {"ALT_SEQUENCE": primer_seq}
                                seq_list.extend([">" + primer_name, primer_seq])
                                hit_list.extend([">alt_" + k + "_ref", reverse_complement(
                                                 primer_seq)])
                                for j in range(7):
                                    alt_start = bt_begin + j
                                    alt_primer_seq = bt_hit["sequence"][j:]
                                    al[j] = {}
                                    al[j]["ALT_START"] = alt_start
                                    al[j]["ALT_SEQUENCE"] = alt_primer_seq
                                    seq_list.extend([">" + primer_name, alt_primer_seq])
                                    hit_list.extend([">alt_" + k + "_" + str(j),
                                                    reverse_complement(alt_primer_seq)])
                                para[k]["ALTERNATIVES"] = al
                            else:
                                para[k]["ALTERNATIVES"] = {}
                                para[k]["ALT_TM"] = 0
                                para[k]["ALT_TM_DIFF"] = 100
                                para[k]["ALT_BOUND"] = False
                        primers['primer_information'][primer_name][bowtie_key].pop(bt_hit_name)
                        break
                """
                except KeyError:
                    #print "KE2"
                    primer_start = primers['primer_information']\
                                    [primer_name]["GENOMIC_START"]
                    primer_end = primers['primer_information']\
                                    [primer_name]["GENOMIC_END"]
                    primer_chr = primers['primer_information']\
                                    [primer_name]["CHR"]
                    primer_ori = primers['primer_information']\
                                    [primer_name]["ORI"]
                    if (primer_ori == bt_ori) and (primer_start == bt_begin) and (primer_chr == bt_chrom):
                        intended = 1
                        primers['primer_information'][primer_name][bowtie_key].pop(bt_hit_name)
                """
                # if the target is intended it should not be added to the cleaned bowtie hits
                # which is for nonspecific targets only
                if not intended:
                    seq_list.extend([">" + primer_name, primer_seq])
                    hit_list.extend([">bt_" + bt_hit_name, reverse_complement(bt_seq)])
                    #print bt_seq

        except KeyError as e:
            print((str(e)))

            continue
    #infile.close()
    #outfile.close()
    if len(hit_list) > 0:
        with open(bowtie2_output_DIR + "seq_for_tm", "w") as seq_for_tm,             open(bowtie2_output_DIR + "hits_for_tm", "w") as hits_for_tm:
            seq_for_tm.write("\n".join(seq_list))
            hits_for_tm.write("\n".join(hit_list))
        melt_ext = id_generator(6)
        tms = hybrid_TM_files(bowtie2_output_DIR, "seq_for_tm", "hits_for_tm", "melting_" + melt_ext,
                                Na=Na, Mg=Mg, conc=conc)
        tm_list = tms.split("\n")
        tm_count = (len(tm_list) - 1)//2
        for i in range(tm_count):
            items = tm_list[i].split(" ")
            h_tm = float(tm_list[i+tm_count+1].split("\t")[-1])
            p_name = items[2]
            h_type = items[4].split("_")[0]
            h_name = items[4].split("_")[-1][:-1]
            if h_type == "alt":
                copyname = items[4].split("_")[1]
                try:
                    if h_name == "ref":
                        primers["primer_information"][p_name]["PARALOG_COORDINATES"][copyname]                                   ["ALTERNATIVES"][h_name]["ALT_TM"] = h_tm
                    else:
                        primers["primer_information"][p_name]["PARALOG_COORDINATES"][copyname]                                       ["ALTERNATIVES"][int(h_name)]["ALT_TM"] = h_tm
                except KeyError:
                    print((p_name, h_name, items))

            else:
                bowtie_dic = primers["primer_information"][p_name][bowtie_key]
                for b in bowtie_dic:
                    if h_name == b:
                        bowtie_dic[b]["TM"] = h_tm
    # write the (if) updated dictionary to the file
    with open(primer3_output_DIR + primer_out, 'w') as outfile:
        json.dump(primers, outfile, indent=1)
    return primers
def parse_bowtie (primer_file, bt_file, primer_out,primer3_output_DIR,                   bowtie2_output_DIR, species, settings, outp = 1):
    """ take a bowtie output file and filter top N hits per primer.
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
    infile = open(bowtie2_output_DIR + bt_file, 'r')
    # read in primer/mip file
    """
    with open (primer3_output_DIR + primer_file, 'r') as handle:
        primers = json.load(handle)
    """
    primers = copy.deepcopy(primer_file)
    # create a temp dic to count hits/primer
    counter_dic = {}
    # create a bowtie key that will be used when adding
    # bowtie information to primers
    bowtie_key ="bowtie_information_" + species
    # all bowtie hits that will be used further for TM analysis
    # will need to have sequence information with them
    # region keys for hits (in chrx:begin-end format) will be
    # kept in a list for mass fasta extraction later.
    keys = []
    #
    # read bowtie hits
    for line in infile:
        try:
            if not line.startswith("@"):
                record = line.strip('\n').split('\t')
                primer_name = record [0]
                # increment hit counter for primer
                try:
                    counter_dic[primer_name] +=1
                except KeyError:
                    counter_dic[primer_name] = 1
                # check how many hits have been analyzed for this primer
                # if upper hit limit has been reached, mark primer for removal
                if counter_dic[primer_name] >= M:
                    primers['primer_information'][primer_name]["remove"] = True
                # move on to the next hit if primer hit limit has been reached.
                # no further hits will be added for those primers
                if counter_dic[primer_name] >= N:
                    continue
                flag = record [1]
                # a flag value of 4 means there was no hit, so pass those lines
                if flag == "4":
                    continue
                # chromosome of the bowtie hit
                chrom = record [2]
                # genomic position of bowtie hit
                pos = int(record [3])
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
                    # Primer's 5' is the hit start when the hit is on forward strand
                    # so the nucleotides are added at start position
                    bt_start = hit_start
                    bt_end = hit_end
                    hit_str = "forward"
                    hit_region_key = chrom + ":" + str(hit_start) + "-" + str(hit_end)
                else:
                    bt_start = hit_end
                    bt_end = hit_start
                    hit_str = "reverse"
                    hit_region_key = chrom + ":" + str(hit_start) + "-" + str(hit_end)
                # add region key to keys list for fasta retrieval later
                keys.append(hit_region_key)
                # add all hit information to primer dictionary
                try:
                    primers["primer_information"][primer_name][bowtie_key][str(counter_dic[primer_name])] =                          {"chrom":chrom, "begin": bt_start, "end":bt_end, "key":hit_region_key,
                          "strand":hit_str}
                except KeyError:
                     primers["primer_information"][primer_name][bowtie_key] = {str(counter_dic[primer_name]):                          {"chrom":chrom, "begin": bt_start, "end":bt_end, "key":hit_region_key,
                          "strand":hit_str}}
        except KeyError as ke:
            # in earlier versions of this function the primers with
            # excessive hits were removed during iteration and that lead
            # to keyerrors. Now there should be no key error
            # but print the key if an exception arises.
            # print ke
            continue
    # get the fasta sequences of all hits
    sequence_dic = get_fasta_list(keys, species)
    # remove primers with too many hits
    for p in list(primers["primer_information"].keys()):
        if "remove" in primers["primer_information"][p]:
            primers["primer_information"].pop(p)
        else:
            # add hit sequences to primer dictionary
            # forward strand hits are added directly
            # reverse strand hits are reversed-complemented
            # so the hit is always in the primer orientation and
            # and similar in sequence
            try:
                for h in primers["primer_information"][p][bowtie_key]:
                    if primers["primer_information"][p][bowtie_key][h]["strand"] == "forward":
                        primers["primer_information"][p][bowtie_key][h]["sequence"] =                              sequence_dic[primers["primer_information"][p][bowtie_key][h]["key"]]
                    else:
                        primers["primer_information"][p][bowtie_key][h]["sequence"] =                              reverse_complement(sequence_dic[primers["primer_information"]                             [p][bowtie_key][h]["key"]])
            except KeyError:
                # if there is no bowtie hit for this primer (happens for host species)
                primers["primer_information"][p][bowtie_key] = {}
    # save the updated primers file
    if outp:
        with open(primer3_output_DIR + primer_out, 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers
def alternative(primer_file, output_file, primer3_output_DIR, tm_diff, outp = 1):
    """ Pick the best alternative arm for primers that do not bind
    all paralogs.
    """

    """
    with open (primer3_output_DIR + primer_file) as infile:
        primer_dic = json.load(infile)
        primers = primer_dic["primer_information"]
    """
    primer_dic = primer_file
    primers = primer_dic["primer_information"]
    try:
        for primer_name in primers:
            primer = primers[primer_name]
            para = primer["PARALOG_COORDINATES"]
            for c in para:
                try:
                    alts = para[c]["ALTERNATIVES"]
                    ref_tm = alts["ref"].pop("ALT_TM")
                    alts.pop("ref")
                    sorted_alts = sorted(alts, key=lambda a:                                abs(alts[a]["ALT_TM"] - ref_tm))
                    if abs(alts[sorted_alts[0]]["ALT_TM"] - ref_tm) <= tm_diff:
                        primer["ALT_BINDS"].append(c)
                        para[c].update(alts[sorted_alts[0]])
                    para[c].pop("ALTERNATIVES")
                except KeyError as e:
                    try:
                        para[c].pop("ALTERNATIVES")
                    except KeyError as e:
                        pass
                except IndexError:
                    try:
                        para[c].pop("ALTERNATIVES")
                    except KeyError as e:
                        pass
    except KeyError as e:
        pass
    if outp:
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(primer_dic, outfile, indent=1)
    return primer_dic
def cleaner_bowtie (primer_file, bt_file, bt_out, primer_out,primer3_output_DIR,                   bowtie2_output_DIR, N=50, M=500):
    """ take a bowtie output file and remove all hits on intended
    target. Then, filter top N hits per primer/mip.
    Create a cleaned sam output file. Return nothing."""
    # read in bowtie file
    infile = open(bowtie2_output_DIR + bt_file, 'r')
    # create cleaned bowtie out file
    outfile = open(bowtie2_output_DIR + bt_out, 'w')
    # read in primer/mip file
    with open (primer3_output_DIR + primer_file, 'r') as handle:
        primers = json.load(handle)
    # create a temp dic to count hits/primer or mip
    temp = {}
    # check if primers or mips are cleaned
    # default is a primer file
    mip = 0
    # change mip to 1 if mips are analyzed
    if "pair_information" in list(primers.keys()):
        mip = 1
    counter = 0
    # read bowtie hits
    for line in infile:
        try:
            if not line.startswith("@"):
                record = line.strip('\n').split('\t')
                # another temp dic to
                primer_name = record [0]
                flag = record [1]
                # a flag value of 4 means there was no hit, so pass those lines
                if flag == "4":
                    continue
                # chromosome of the bowtie hit
                chrom = record [2]
                # genomic position of bowtie hit
                pos = int(record [3])
                # get cigar string of alignment
                cigar = record[5]
                # extract which strand is the bowtie hit on
                # true if forward
                strand = ((int(record[1]) % 256) == 0)
                # create an intended variable that will be true if the bt hit is on
                # an intended target
                intended = 0
                # get sequence of mip (if mip) from primer dic
                if mip:
                    if "pairs" in list(primers["pair_information"][primer_name].keys()):
                        para = primers["pair_information"][primer_name]["pairs"]
                        for p in para:
                            mip_start = para[p]["mip_start"] -10
                            mip_end = para[p]["mip_end"] + 10
                            mip_chrom = para[p]["chrom"]
                            if (mip_chrom == chrom) and (mip_start <= pos <= mip_end):
                                intended = 1
                # get sequence of primer from dic
                else:
                    if "PARALOG_COORDINATES" in list(primers['primer_information'][primer_name].keys()):
                        para = primers['primer_information'][primer_name]["PARALOG_COORDINATES"]
                        # para is a dic like {C0:{"CHR": "chr4", "GENOMIC_START" ..}, C1:{..
                        for k in para:
                            primer_ori = para[k]["ORI"]
                            primer_chr = para[k]["CHR"]
                            # expected bowtie hit position depends on the orientation of the primer
                            # orientation of primer and bowtie hit must be the same
                            # and genomic start is the hit position for forward, and end for reverse
                            if primer_ori == "reverse" and not strand and                                (para[k]["GENOMIC_END"]+10) >= pos >= (para[k]["GENOMIC_END"]-10):
                                intended = 1
                                para[k]["BOWTIE_END"] = pos
                                para[k]["BOWTIE_START"] = pos + get_cigar_length(cigar) - 1
                            elif primer_ori == "forward" and strand and                                  (para[k]["GENOMIC_START"]+10) >= pos >= (para[k]["GENOMIC_START"]-10):
                                intended = 1
                                para[k]["BOWTIE_START"] = pos
                                para[k]["BOWTIE_END"] = pos + get_cigar_length(cigar) - 1
                    else:
                        primer_start = primers['primer_information']                                        [primer_name]["GENOMIC_START"]
                        primer_end = primers['primer_information']                                        [primer_name]["GENOMIC_END"]
                        primer_chr = primers['primer_information']                                        [primer_name]["CHR"]
                        primer_ori = primers['primer_information']                                        [primer_name]["ORI"]
                        if (primer_ori == "forward") and strand and (primer_start == pos) and (primer_chr == chrom):
                            intended = 1
                            #print primer_ori, strand, primer_start, pos, primer_chr, chrom
                        elif (primer_ori == "reverse") and (not strand) and (primer_end == pos) and (primer_chr == chrom):
                            intended = 1
                            #print primer_ori, strand, primer_start, pos, primer_chr, chrom
                        elif counter < 50:
                            #print primer_ori, strand, primer_start, pos, primer_chr, chrom
                            counter += 1

                # if the target is intended it should not be added to the cleaned bowtie hits
                # which is for nonspecific targets only
                if intended:
                    continue
                # if it is a nonspecific hit, check how many hits for this primer
                # have already been recorded. If none, create a counter for that primer.
                elif not primer_name in temp:
                    temp[primer_name] = 1
                else:
                    temp[primer_name] += 1
                # if there are less then N hits recorded, write the line to file
                if temp[primer_name] < int(N):
                    outfile.write(line)
                if temp[primer_name] > int(M):
                    if "pair_information" in list(primers.keys()):
                        if primer_name in list(primers["pair_information"].keys()):
                            primers["pair_information"].pop(primer_name)
                    elif primer_name in list(primers["primer_information"].keys()):
                        primers["primer_information"].pop(primer_name)
        except KeyError:
            continue
    infile.close()
    outfile.close()
    # write the (if) updated dictionary to the file
    if primer_out != "":
        # mip dict is not updated
        with open(primer3_output_DIR + primer_out, 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers
def cleaner_host_bowtie (primer_file, bt_file, bt_out,primer3_output_DIR,                   bowtie2_output_DIR, N=50):
    """ take a bowtie output file and remove all hits on intended
    target. Then, filter top N hits per primer/mip.
    Create a cleaned sam output file. Return nothing."""
    # read in bowtie file
    infile = open(bowtie2_output_DIR + bt_file, 'r')
    # create cleaned bowtie out file
    outfile = open(bowtie2_output_DIR + bt_out, 'w')
    # read in primer/mip file
    with open (primer3_output_DIR + primer_file, 'r') as handle:
        primers = json.load(handle)
    # create a temp dic to count hits/primer or mip
    temp = {}
    # check if primers or mips are cleaned
    # default is a primer file
    mip = 0
    # change mip to 1 if mips are analyzed
    if "pair_information" in list(primers.keys()):
        mip = 1
    # read bowtie hits
    for line in infile:
        if not line.startswith("@"):
            record = line.strip('\n').split('\t')
            # another temp dic to
            primer_name = record [0]
            flag = record [1]
            # a flag value of 4 means there was no hit, so pass those lines
            if flag == "4":
                continue
            # chromosome of the bowtie hit
            chrom = record [2]
            # genomic position of bowtie hit
            pos = int(record [3])
            # get cigar string of alignment
            cigar = record[5]
            # extract which strand is the bowtie hit on
            # true if forward
            strand = ((int(record[1]) % 256) == 0)
            if not primer_name in temp:
                temp[primer_name] = 1
            else:
                temp[primer_name] += 1
            # if there are less then N hits recorded, write the line to file
            if temp[primer_name] < int(N):
                outfile.write(line)
    infile.close()
    outfile.close()
    return primers
def add_bt_worker (primer_name, primer_seq, strand, region, mip, fasta_genome, tmp, settings):
    """ Worker function for add_bowtie_multi function.
    Return TM of each bowtie hit passed from parent function."""
    # get the sequence of the region where the bt hit was mapped
    fasta=subprocess.check_output(["samtools",  "faidx", fasta_genome, region])
    # extract 1 line sequence from fasta output
    fasta_list=fasta.split("\n")
    seq_template_temp = "".join(fasta_list[1:])
    Na = float(settings["Na"])
    Mg = float(settings["Mg"])
    conc = float(settings["oligo_conc"])
    # if primer is mapped to reverse strand
    if not strand:
        # then primer will be complementary to forward strand
        # and seq_template_temp is the sequence of forward strand
        seq_template = seq_template_temp
    # if mapped to forward strand
    else:
        # then primer is complementary to reverse strand
        # so we get reverse complement of mapped sequence
        seq_template = reverse_complement(seq_template_temp)
    # ntthal program used for TM checking does not work for >60nt oligos
    # to use the program on MIPS, get subsequences of the MIP and
    # use the maximum temperature. It is not optimal but works fine.
    TM = hybrid_TM(tmp, seq_template, primer_seq, Na=Na, Mg=Mg, conc=conc)
    # add strand information to bowtie dict
    # add TM and region information to bowtie dict
    return TM
def add_bowtie_multi (primer_file, bowtie_out_file, output_file,                      primer3_output_DIR, bowtie2_output_DIR, num_processors,settings, species="pf"):
    """Adds bowtie information to a given primer (or MIP) dictionary file (json object).
    Once done, primer dic still has 2 keys, sequence_information and primer_information.
    Each primer in primer information dic will have a new key: bowtie_information_"species"
    whose value is a list of dictionaries (for each bt hit).
    Mip dictionary will still have 2 keys, sequence_information and pair_information.
    Each primer pair will have bowtie information associated with the mip they make.
    bowtie dictionaries have keys primer_name, flag, chrom, pos, cigar, strand, TM,
    region and other (other parameters not used for now),  with corresponding values.
    if a mip file is given, mip_information has new key bowtie_information_"species"."""
    file_locations = get_file_locations()
    # assign fasta genome file for species
    fasta_genome = file_locations[species]["fasta_genome"]
    # get primer dictionary from file
    with open(primer3_output_DIR + primer_file, "r") as handle:
        primers = json.load(handle)
    # are we adding bowtie information for mips or primers?
    if "pair_information" in list(primers.keys()):
        # then primers already paired and we are looking at MIP data
        mip = 1
    else:
        # then we are looking at individual primers
        mip = 0
    # open bowtie output file
    infile = open (bowtie2_output_DIR + bowtie_out_file, 'r')
    # bowtie_key will be added to primer_information dictionaries for each hit.
    bowtie_key ="bowtie_information_" + species
    # read each line of bowtie output file into list to be fed to worker func
    lines = infile.readlines()
    infile.close()
    # create a temp directory in temp_dir
    tmp = bowtie2_output_DIR + "tmp/"
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    # start multiprocess
    #if __name__ == "__main__":
    # create a pool of "num_process" processors
    p = Pool(num_processors)
    # results of each process will be appended to a list
    worker_results = []
    # for each line in bowtie file, extract alignment info
    for l in lines:
        if not l.startswith("@"):
            record = l.strip('\n').split('\t')
            # create a dic to hold the hit information
            temp_dic = {}
            temp_dic["primer_name"] = record [0]
            temp_dic["flag"] = record [1]
            # a flag value of 4 means there was no alignment
            # so no need to continue with this line/hit
            if temp_dic["flag"] == "4":
                continue
            # chromosome of the bowtie hit
            temp_dic["chrom"] = record [2]
            # genomic position of bowtie hit
            temp_dic["pos"] = record [3]
            # CIGAR string of the hit
            temp_dic["cigar"] = record [5]
            # other information belonging the hit
            # that is not useful for now
            temp_dic["other"] = record [4:5] + record [6:]
            # extract which strand is the bowtie hit on
            # true if forward
            strand = ((int(record[1]) % 256) == 0)
            if strand:
                temp_dic["strand"] = "forward"
            else:
                temp_dic["strand"] = "reverse"
            # get melting temperature of the hit
            primer_name = temp_dic["primer_name"]
            # get start and end coordinates of both primers in mip
            coord = []
            try:
                if mip:
                    primer_seq = primers['pair_information']                                [primer_name]['mip_information']["SEQUENCE"]
                    ext_start = primers['pair_information']                               [primer_name]['extension_primer_information']                               ["GENOMIC_START"]
                    ext_end = primers['pair_information']                             [primer_name]['extension_primer_information']                             ["GENOMIC_END"]
                    lig_start = primers['pair_information']                                [primer_name]['ligation_primer_information']                                ["GENOMIC_START"]
                    lig_end = primers['pair_information']                             [primer_name]['ligation_primer_information']                             ["GENOMIC_END"]
                    primer_chr = primers['pair_information']                             [primer_name]['ligation_primer_information']["CHR"]
                    coord.append(ext_start)
                    coord.append(ext_end)
                    coord.append(lig_start)
                    coord.append(lig_end)
                # get start and end coordinates of primer
                else:
                    primer_seq = primers['primer_information']                                        [primer_name]["SEQUENCE"]
                    primer_start = primers['primer_information']                                          [primer_name]["GENOMIC_START"]
                    primer_end = primers['primer_information']                                        [primer_name]["GENOMIC_END"]
                    primer_chr = primers['primer_information']                                        [primer_name]["CHR"]
                    coord.append(primer_start)
                    coord.append(primer_end)
            except KeyError as e:
                continue
            # check if the hit is on intended target
            # no specific hit should remain after using cleaner_bowtie function
            # so this is just double checking
            #if temp_dic["chrom"] == primer_chr and (temp_dic["pos"] in coord):
            #    continue
            # get the coordinates of bowtie hit +/- 5 bp for hybridization test
            region = temp_dic["chrom"] + ":" + str(int(temp_dic["pos"])-5)            + "-" + str(int(temp_dic["pos"])+len(primer_seq)+5)
            # create TM list to be populated by the worker function
            temp_dic["TM"] = []
            # pick a processor from the pool and pass necessary arguments
            # callback function will append the calculated TM to the list
            res = p.apply_async(add_bt_worker, (primer_name, primer_seq, strand,                                region, mip, fasta_genome, tmp, settings), callback=temp_dic["TM"].append)
            # add the bowtie information dictionary to primer dictionary
            # all bowtie dictionaries for a single mip/primer are appended
            # to a list
            if mip:
                if bowtie_key in list(primers['pair_information']                   [primer_name]['mip_information'].keys()):
                    primers['pair_information'][primer_name]                    ['mip_information'][bowtie_key].append(temp_dic)
                else:
                    primers['pair_information'][primer_name]                    ['mip_information'][bowtie_key] = [temp_dic]
            else:
                if bowtie_key in list(primers['primer_information'][primer_name].keys()):
                    primers['primer_information'][primer_name][bowtie_key].append(temp_dic)
                else:
                    primers['primer_information'][primer_name][bowtie_key] = [temp_dic]
    # close the pool of processor since no other processes will be added
    p.close()
    # wait for each process to finish and then continue with remaining code
    p.join()
    # write the updated dictionary to json file in primer3_output_DIR.
    outfile = open (primer3_output_DIR + output_file, "w")
    json.dump(primers, outfile, indent=1)
    outfile.close()
    # remove files by the tm calculation function
    shutil.rmtree(tmp)
    return primers
def filter_bowtie (primer_file,
                   output_file,
                   primer3_output_DIR,
                   species,
                   TM=46,
                   hit_threshold=0,
                   lower_tm=46,
                   lower_hit_threshold=3,
                  outp = 1):
    """ Check TMs of bowtie hits of given primers,
    on a given genome.Filter the primers with too many
    nonspecific hits."""

    """
    infile = open(primer3_output_DIR + primer_file, 'r')
    primers = json.load(infile)
    infile.close()
    """
    # load primer dictionary
    primers = primer_file
    #print TM, hit_threshold, lower_tm, lower_hit_threshold, species
    #print str(len(primers["primer_information"])) + " primers analyzed!"
    for primer in list(primers["primer_information"].keys()):
        # create a hit count parameter for hits with significant tm
        hc = 0
        lhc = 0
        # check if bowtie information exists in dic
        try:
            bowtie = primers["primer_information"][primer]                     ["bowtie_information_" + species]
            for h in bowtie:
                hit = bowtie[h]
                try:
                    # check if TM information is included in bowtie
                    # for the same TypeError as the mips above
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
                except KeyError as e:
                    continue

            if primer in list(primers["primer_information"].keys()):
                primers["primer_information"][primer].pop("bowtie_information_" + species)
        except KeyError as ke:
            #print ke
            continue
    if outp:
        # write dictionary to file in primer3_output_DIR
        outfile = open (primer3_output_DIR + output_file, 'w')
        json.dump(primers, outfile, indent=1)
        outfile.close()
    return primers
def filter_bowtie_hits (primer_file, output_file, primer3_output_DIR, species, TM=46, hit_threshold=0,
                        lower_tm=46, lower_hit_threshold=3):
    """ Check TMs of bowtie hits of given mips or primers,
    on a given genome.if high, make sure they are the intended
    target. Otherwise, remove primer. Then remove bowtie
    information from dictionary."""
    supported_species = ["pf", "hs"]
    # check if species is supported
    if not species in supported_species:
        return "Species not supported"
    # load primer dictionary
    infile = open(primer3_output_DIR + primer_file, 'r')
    primers = json.load(infile)
    infile.close()
    # check if it is a mip dictionary or primer dictionary
    if "pair_information" in list(primers.keys()):
        # report how many primer pairs are analyzed
        #print str(len(primers["pair_information"])) + " mips analyzed!"
        # remove bowtie hits with TMs higher than given TM
        for primer in list(primers["pair_information"].keys()):
            # create a hit count parameter for hits with significant tm
            hc = 0
            lhc = 0
            # check if bowtie_information key exists for given species
            if ("bowtie_information_" + species) in list(primers["pair_information"]               [primer]["mip_information"].keys()):
                # extract bowtie information
                bowtie = primers["pair_information"][primer]["mip_information"]                         ["bowtie_information_" + species]
                for hit in bowtie:
                    # get TM of the hit
                    # depending on what function had been used to add bowtie
                    # information to the dictionary, TM can be a list or a float
                    # Therefore, we have to check if it is a string or not first.
                    try:
                        # if TM not a list, TypeError is raised
                        if ("TM" in list(hit.keys())) and (hit["TM"][0] >= TM)                            and (primer in list(primers["pair_information"].keys())):
                            hc += 1
                            if hc > hit_threshold:
                                primers["pair_information"].pop(primer)
                                break
                        elif ("TM" in list(hit.keys())) and (hit["TM"][0] >= lower_tm)                            and (primer in list(primers["pair_information"].keys())):
                            lhc += 1
                            if lhc > lower_hit_threshold:
                                primers["pair_information"].pop(primer)
                                break
                    except TypeError as e:
                        if ("TM" in list(hit.keys())) and (hit["TM"] >= TM)                            and (primer in list(primers["pair_information"].keys())):
                            hc += 1
                            if hc > hit_threshold:
                                primers["pair_information"].pop(primer)
                                break
                        elif ("TM" in list(hit.keys())) and (hit["TM"] >= lower_tm)                            and (primer in list(primers["pair_information"].keys())):
                            lhc += 1
                            if lhc > lower_hit_threshold:
                                primers["pair_information"].pop(primer)
                                break
                # if a primer is still in dictionary after filtering for bowtie hits
                # keep the primer and remove the bowtie hits
                if primer in list(primers["pair_information"].keys()):
                    primers["pair_information"][primer]                           ["mip_information"].pop("bowtie_information_" + species)

        # report how many mips are left after filtering
        #print str(len(primers["pair_information"])) +\
        #      " mips with mispriming TM < " + str(TM)
    # if dictionary is for primers only:
    else:
        #print str(len(primers["primer_information"])) + " primers analyzed!"
        for primer in list(primers["primer_information"].keys()):
            # create a hit count parameter for hits with significant tm
            hc = 0
            lhc = 0
            # check if bowtie information exists in dic
            if ("bowtie_information_" + species) in                list(primers["primer_information"][primer].keys()):
                # get bowtie information for the primer
                bowtie = primers["primer_information"][primer]                         ["bowtie_information_" + species]
                for hit in bowtie:
                    # check if TM information is included in bowtie
                    # for the same TypeError as the mips above
                    try:
                        if ("TM" in list(hit.keys())) and (float(hit["TM"][0]) >= TM)                            and (primer in list(primers["primer_information"].keys())):
                            hc += 1
                            if hc > hit_threshold:
                                primers["primer_information"].pop(primer)
                                break
                        elif ("TM" in list(hit.keys())) and (float(hit["TM"][0]) >= lower_tm)                            and (primer in list(primers["primer_information"].keys())):
                            lhc += 1
                            if lhc > lower_hit_threshold:
                                primers["primer_information"].pop(primer)
                                break

                    except TypeError as e:
                        if ("TM" in list(hit.keys())) and (float(hit["TM"]) >= TM)                            and (primer in list(primers["primer_information"].keys())):
                            hc += 1
                            if hc > hit_threshold:
                                primers["primer_information"].pop(primer)
                                break
                        elif ("TM" in list(hit.keys())) and (float(hit["TM"]) >= lower_tm)                            and (primer in list(primers["primer_information"].keys())):
                            lhc += 1
                            if lhc > lower_hit_threshold:
                                primers["primer_information"].pop(primer)
                                break

                if primer in list(primers["primer_information"].keys()):
                    primers["primer_information"][primer].pop("bowtie_information_" + species)

        #print str(len(primers["primer_information"])) + " primers with mispriming TM < " + str(TM)
    # write dictionary to file in primer3_output_DIR
    outfile = open (primer3_output_DIR + output_file, 'w')
    json.dump(primers, outfile, indent=1)
    outfile.close()
    return primers
def filter_bowtie_hits_dup (primer_file, output_file, species, primer3_output_DIR, TM=46):
    """ Check TMs of bowtie hits of given mips or primers,
    on a given genome.if high, make sure they are the intended
    target. Otherwise, remove primer. Then remove bowtie
    information from dictionary.
    This is for use in mips for duplicated regions. It should not filter
    when the hit is on a duplicated copy."""
    supported_species = ["pf", "hs"]
    # check if species is supported
    if not species in supported_species:
        return "Species not supported"
    # load primer dictionary
    infile = open(primer3_output_DIR + primer_file, 'r')
    primers = json.load(infile)
    infile.close()
    # check if it is a mip dictionary or primer dictionary
    if "pair_information" in list(primers.keys()):
        # report how many primer pairs are analyzed
        #print str(len(primers["pair_information"])) + " mips analyzed!"
        # remove bowtie hits with TMs higher than given TM
        for primer in list(primers["pair_information"].keys()):
            # check if bowtie_information key exists for given species
            if ("bowtie_information_" + species) in list(primers["pair_information"]               [primer]["mip_information"].keys()):
                # extract bowtie information
                bowtie = primers["pair_information"][primer]["mip_information"]                         ["bowtie_information_" + species]
                for hit in bowtie:
                    # get TM of the hit
                    # depending on what function had been used to add bowtie
                    # information to the dictionary, TM can be a list or a string
                    # Therefore, we have to check if it is a string or not first.
                    try:
                        if ("TM" in list(hit.keys())) and (float(hit["TM"]) >= TM)                            and (primer in list(primers["pair_information"].keys())):
                            # remove primers with TM higher than specified
                            primers["pair_information"].pop(primer)
                            break
                    except TypeError as e:
                        # if TM is a list and not a string, TypeError is raised
                        if ("TM" in list(hit.keys())) and (float(hit["TM"][0]) >= TM)                            and (primer in list(primers["pair_information"].keys())):
                            primers["pair_information"].pop(primer)
                            break
                # if a primer is still in dictionary after filtering for bowtie hits
                # keep the primer and remove the bowtie hits
                if primer in list(primers["pair_information"].keys()):
                    primers["pair_information"][primer]                           ["mip_information"].pop("bowtie_information_" + species)
        # report how many mips are left after filtering
        #print str(len(primers["pair_information"])) +\
        #      " mips with mispriming TM < " + str(TM)
    # if dictionary is for primers only:
    else:
        #print str(len(primers["primer_information"])) + " primers analyzed!"
        for primer in list(primers["primer_information"].keys()):
            # check if bowtie information exists in dic
            if ("bowtie_information_" + species) in                list(primers["primer_information"][primer].keys()):
                # get bowtie information for the primer
                bowtie = primers["primer_information"][primer]                         ["bowtie_information_" + species]
                for hit in bowtie:
                    # check if TM information is included in bowtie
                    # for the same TypeError as the mips abova
                    # use try/except statements
                    try:
                        if ("TM" in list(hit.keys())) and (float(hit["TM"]) >= TM)                            and (primer in list(primers["primer_information"].keys())):
                            primers["primer_information"].pop(primer)
                            break
                    except TypeError as e:
                        if ("TM" in list(hit.keys())) and (float(hit["TM"][0]) >= TM)                            and (primer in list(primers["primer_information"].keys())):
                            primers["primer_information"].pop(primer)
                            break
                if primer in list(primers["primer_information"].keys()):
                    primers["primer_information"][primer].pop("bowtie_information_" + species)
        #print str(len(primers["primer_information"])) + " primers with mispriming TM < " + str(TM)
    # write dictionary to file in primer3_output_DIR
    outfile = open (primer3_output_DIR + output_file, 'w')
    json.dump(primers, outfile, indent=1)
    outfile.close()
    return primers
def pick_paralog_primer_pairs (ext_file,
                               lig_file,
                               output_file,
                               primer3_output_DIR,
                               min_size, max_size,
                               alternative_arms,
                              region_insertions,
                              subregion_name,
                              outp = 1):
    """ Pick primer pairs from extension and ligation arm candidate
    dictionary files for a given size range. Output dic has 2 keys:
    sequence_information: template sequence information
    pair_information: dictionary with PAIR_NAME:{extenison_primer_information:...
    ligation_primer... orientation:fw/rv}. """
    # load extension and ligation dictionaries
    """
    infile1 = open(primer3_output_DIR + ext_file, 'r')
    infile2 = open(primer3_output_DIR + lig_file, 'r')
    extension = json.load(infile1)
    ligation = json.load(infile2)
    infile1.close()
    infile2.close()
    """
    extension = ext_file
    ligation = lig_file
    # assign primer information dictionaries to a shorter name
    ext = extension["primer_information"]
    lig = ligation["primer_information"]
    # check if extension and ligation dictionaries have primers
    if len(ext) == 0:
        print("There are no extension primers.")
        return 1
    if len(lig) == 0:
        print("There are no ligation primers.")
        return 1
    # assign sequence information dict to shorter name
    ext_seq = extension["sequence_information"]
    lig_seq = ligation["sequence_information"]
    # create a primer pairs dic. This dictionary is similar to primer dic
    primer_pairs = {}
    # has the same sequence_information key:value pairs
    primer_pairs["sequence_information"] = {}
    # has pair information key instead of primer_information
    primer_pairs["pair_information"] = {}
    # populate sequence information (same as extension or ligation)
    primer_pairs["sequence_information"]['SEQUENCE_TEMPLATE'] =       extension["sequence_information"]['SEQUENCE_TEMPLATE']
    primer_pairs["sequence_information"]['SEQUENCE_EXCLUDED_REGION'] =       extension["sequence_information"]['SEQUENCE_EXCLUDED_REGION']
    primer_pairs["sequence_information"]['SEQUENCE_TARGET'] =     extension["sequence_information"]['SEQUENCE_TARGET']
    primer_pairs["sequence_information"]['SEQUENCE_ID'] =     extension["sequence_information"]['SEQUENCE_ID']
    # pick primer pairs
    for e in list(ext.keys()):
        # extension primer information for this mip will be e_info
        e_info = ext[e]
        # get primer coordinates
        ext_start = e_info["GENOMIC_START"]
        ext_end = e_info["GENOMIC_END"]
        # get primer orientation
        ext_ori = ext_end > ext_start
        # if end is greater than start then it is a left(fw) primer, ext_ori is True
        # get coordinates of this primer in paralog copies
        ep_info = e_info["PARALOG_COORDINATES"]
        # paralogs bound by primer
        e_binds = e_info["BOWTIE_BINDS"]
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
                    # ligation end should be greater than extension end for forward pairs
                    position = lig_end > ext_end
                else:
                    # extension end should be greater than ligation end for reverse pairs
                    position = ext_end > lig_end
                # get pair information if relative positions of primers are correct
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
                                    p_coord = [ep_start, ep_end, lp_start, lp_end]
                                    p_coord.sort()
                                    prod_size = p_coord[-1] - p_coord[0] + 1
                                    pairs[p] = {"capture_size":prod_size,"extension_start": ep_start,
                                                "extension_end": ep_end, "ligation_start": lp_start,
                                                "ligation_end": lp_end, "mip_start":p_coord[0],
                                                "mip_end":p_coord[3], "capture_start":p_coord[1]+1,
                                                "capture_end":p_coord[2]-1, "chrom":lp_chrom,
                                                "orientation":pair_ori
                                                 }


                        except KeyError as e:
                            continue
                    # check if any pairs' product is within size limits
                    pair_found = 0
                    captured_copies = []
                    for p in list(pairs.keys()):
                        if not region_insertions.empty:
                            max_insertion_size = region_insertions.loc[
                                (region_insertions["copy_chrom"]
                                 == pairs[p]["chrom"])
                                &(region_insertions["copy_begin"]
                                  > pairs[p]["capture_start"])
                                &(region_insertions["copy_end"]
                                  < pairs[p]["capture_end"]),
                                "max_size"].sum()
                        else:
                            max_insertion_size = 0
                        adjusted_max_size = max((max_size - max_insertion_size),
                                                min_size)
                        adjusted_min_size = min(adjusted_max_size - 30,
                                               min_size)
                        if (adjusted_max_size
                            >= pairs[p]["capture_size"]
                            >= adjusted_min_size):
                            captured_copies.append(p)
                            # only take mips that capture the reference copy
                            # remove the if statement if any copy captured would work
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
                                    &(region_insertions["copy_begin"]
                                      > pairs[p]["capture_start"])
                                    &(region_insertions["copy_end"]
                                      < pairs[p]["capture_end"]),
                                    "max_size"].sum()
                            else:
                                max_insertion_size = 0
                            adjusted_max_size = max((max_size - max_insertion_size),
                                                    min_size)
                            if adjusted_max_size >= pairs[p]["capture_size"] >= 0:
                                captured_copies.append(p)
                        # create a pair name as PAIR_extension primer number
                        # _ligation primer number
                        ext_name = e.split('_')[2]
                        lig_name = l.split('_')[2]
                        pair_name = "PAIR_" + subregion_name + "_"+ ext_name + "_" + lig_name
                        if ext_ori:
                            orientation = "forward"
                        else:
                            orientation = "reverse"
                        primer_pairs["pair_information"][pair_name] = {"pairs":pairs,                                   "extension_primer_information":ext[e],                                   "ligation_primer_information":lig[l], "orientation": orientation,                                   "captured_copies": captured_copies}
                        alt_paralogs = list(
                            (set(l_alt_binds).union(e_alt_binds)
                            ).difference(paralogs)
                        )
                        alts = {}
                        for a in alt_paralogs:
                            try:
                                alt_arms = []
                                p_coord = []
                                if ep_info[a]["BOWTIE_BOUND"]:
                                    ep_start = ep_info[a]["BOWTIE_START"]
                                    ep_end = ep_info[a]["BOWTIE_END"]
                                else:
                                    try:
                                        ep_start = ep_info[a]["ALT_START"]
                                        ep_end = ep_info[a]["ALT_END"]
                                        alt_arms.append("extension")
                                    except KeyError as e:
                                        continue
                                ep_ori = ep_end > ep_start
                                if lp_info[a]["BOWTIE_BOUND"]:
                                    lp_start = lp_info[a]["BOWTIE_START"]
                                    lp_end = lp_info[a]["BOWTIE_END"]
                                else:
                                    try:
                                        lp_start = lp_info[a]["ALT_START"]
                                        lp_end = lp_info[a]["ALT_END"]
                                        alt_arms.append("ligation")
                                    except KeyError as e:
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
                                        p_coord = [ep_start, ep_end, lp_start, lp_end]
                                        p_coord.sort()
                                        prod_size = p_coord[-1] - p_coord[0] + 1
                                        alts[a] = {"capture_size":prod_size,"extension_start": ep_start,
                                                    "extension_end": ep_end, "ligation_start": lp_start,
                                                    "ligation_end": lp_end, "mip_start":p_coord[0],
                                                    "mip_end":p_coord[3], "capture_start":p_coord[1]+1,
                                                    "capture_end":p_coord[2]-1, "chrom":lp_chrom,
                                                    "orientation":pair_ori, "alternative_arms":alt_arms
                                                     }
                            except KeyError as e:
                                continue
                        # check if any pairs' product is within size limits
                        captured_copies = []
                        for a in list(alts.keys()):
                            # does it satisfy arm setting?
                            good_alt = 0
                            if alternative_arms == "any":
                                good_alt = 1
                            elif (len(alts[a]["alternative_arms"]) == 1) and                                 ((alternative_arms == alts[a]["alternative_arms"][0]) or                                   (alternative_arms == "one")):
                                good_alt = 1
                            if good_alt:
                                if not region_insertions.empty:
                                    max_insertion_size = region_insertions.loc[
                                        (region_insertions["copy_chrom"]
                                         == alts[a]["chrom"])
                                        &(region_insertions["copy_begin"]
                                          > alts[a]["capture_start"])
                                        &(region_insertions["copy_end"]
                                          < alts[a]["capture_end"]),
                                        "max_size"].sum()
                                else:
                                    max_insertion_size = 0
                                adjusted_max_size = max((max_size - max_insertion_size),
                                                        min_size)
                                if adjusted_max_size >= alts[a]["capture_size"] >= 0:
                                    captured_copies.append(a)
                                    primer_pairs["pair_information"][pair_name]["pairs"][a] = alts[a]
                        primer_pairs["pair_information"][pair_name]["alt_copies"] = captured_copies
    # return if no pairs found
    if len(primer_pairs["pair_information"]) == 0:
        print("No primer pairs found.")
        return 1
    # print how many primer pairs are found
    #print str(len(primer_pairs["pair_information"])) + " primer pairs found!"
    # write dic to file in primer_output_DIR
    if outp:
        outfile = open(primer3_output_DIR + output_file, 'w')
        json.dump (primer_pairs, outfile, indent=1)
        outfile.close()
    return primer_pairs
def add_capture_sequence(pairs_file, output_file, primer3_output_DIR,species,
                        outp = 1):
    """
    with open(primer3_output_DIR + pairs_file) as infile:
        primer_pairs = json.load(infile)
    """
    primer_pairs = pairs_file
    capture_keys = []
    for p_pair in primer_pairs["pair_information"]:
        pairs = primer_pairs["pair_information"][p_pair]["pairs"]
        for p in pairs:
            paralog_key = pairs[p]["chrom"] + ":" +                           str(pairs[p]["capture_start"]) + "-" +                          str(pairs[p]["capture_end"])
            pairs[p]["capture_key"] = paralog_key
            capture_keys.append(paralog_key)
    capture_sequence_dic = get_fasta_list(capture_keys, species)
    for p_pair in primer_pairs["pair_information"]:
        pairs = primer_pairs["pair_information"][p_pair]["pairs"]
        for p in pairs:
            if pairs[p]["orientation"] == "forward":
                pairs[p]["capture_sequence"] =                    capture_sequence_dic[pairs[p]["capture_key"]]
            else:
                pairs[p]["capture_sequence"] = reverse_complement(
                    capture_sequence_dic[pairs[p]["capture_key"]])
    if outp:
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(primer_pairs, outfile, indent=1)
    return primer_pairs
def make_mips (primer_pairs,
               output_file,
               primer3_output_DIR,
               mfold_input_DIR,
               backbone=None,
              outp = 1):
    """ Make mips from primer pairs by taking reverse complement
    of ligation primer sequence, adding a backbone sequence and
    the extension primer. Standard backbone is used if none
    specified. Add a new key to each primer pair:
    "mip_information" with a dictionary that has SEQUENCE key
    and mip sequence as value."""
    # load primer pair dictionary
    """
    infile = open(primer3_output_DIR + primer_pairs, 'r')
    pairs = json.load(infile)
    infile.close()
    """
    pairs = primer_pairs
    # check if the primer dictionary is empty
    if len(pairs["pair_information"]) == 0:
        print("There are no primer pairs in dictionary")
        return 1
    # use standard backbone if none is specified
    if backbone == None:
        backbone = "NNNNNNCTTCAGCTTCCCGATCCGACGGTAGTGTNNNNNN"
    # get primer sequences for each primer pair
    for primers in pairs["pair_information"]:
        extension_sequence = pairs["pair_information"][primers]        ["extension_primer_information"]["SEQUENCE"]
        ligation_sequence = pairs["pair_information"][primers]        ["ligation_primer_information"]["SEQUENCE"]
        # reverse complement ligation primer
        ligation_rc = reverse_complement(ligation_sequence)
        # add sequences to make the mip
        mip = ligation_rc + backbone + extension_sequence
        # create a dictionary to hold mip information
        mip_dic = {"ref":{'SEQUENCE':mip,
                           "captures": copy.deepcopy(pairs["pair_information"]\
                                                [primers]["captured_copies"])}}
        # add it to the pair dictionary

        if "alt_copies" in list(pairs["pair_information"][primers].keys()):
            alt_sequences = {}
            alt_counter = 0
            alt = pairs["pair_information"][primers]["alt_copies"]
            p_para = pairs["pair_information"][primers]["pairs"]
            e_para = pairs["pair_information"][primers]                    ["extension_primer_information"]["PARALOG_COORDINATES"]
            l_para = pairs["pair_information"][primers]                    ["ligation_primer_information"]["PARALOG_COORDINATES"]
            for a in alt:
                if "extension" in p_para[a]["alternative_arms"]:
                    extension_sequence = e_para[a]["ALT_SEQUENCE"].upper()
                if "ligation" in p_para[a]["alternative_arms"]:
                    ligation_sequence = l_para[a]["ALT_SEQUENCE"].upper()
                value_found = 0
                for key,value in list(alt_sequences.items()):
                    if [extension_sequence, ligation_sequence] == value["sequences"]:
                        value_found = 1
                        value["copies"].append(a)
                        break
                if not value_found:
                    alt_sequences[alt_counter] = {"sequences":[extension_sequence, ligation_sequence],
                                                  "copies":[a]}
                    alt_counter += 1

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
    outfile = open(mfold_input_DIR + output_file, "w")
    for primers in pairs["pair_information"]:
        outline = ">" + primers + "\n" + pairs["pair_information"]        [primers]["mip_information"]["ref"]['SEQUENCE'] + "\n"
        outfile.write(outline)
    outfile.close()
    # write mip dictionary to file in primer3_output_DIR
    if outp:
        outfile = open(primer3_output_DIR + output_file, 'w')
        json.dump(pairs, outfile, indent=1)
        outfile.close()
    return pairs
def check_hairpin (input_fasta, output_json, primer3_output_DIR, bowtie2_input_DIR, mfold_input_DIR, dg=0.8):
    """ Check free energy values of secondary structures formed by mips @60°C.
    Run mfold_quiker script with mip condition parameters, get dG, dH, dS and TM
    values of each mip. TM values are a bit rough according to the author of the program.
    Therefore dG values @ 60°C is used to determine hairpin stability. All mips have
    a secondary structure with 0.82 kcal dG so 0.8 seem like a good cut off value to
    eliminate mips. This roughly means that any mip that is selected will require
    at least 0.8 kcal energy input to form a stable hairpin.
    Add HAIRPIN key is to the mip_information dict, which contains hairpin information:
    minimum dG, dH, dS and TM."""
    # run mfold_quiker script on mip sequences. current working directory is changed when
    # using subprocess.call because the input sequence must be in the directory where the
    # script is run.
    subprocess.check_call(["mfold_quiker", "SEQ="+input_fasta, "T=60"], cwd=mfold_input_DIR, stderr=subprocess.STDOUT)
    # a number of files are generated by the script using the input base name.
    # .det file has the dG TM information we need.
    infile = open (mfold_input_DIR + input_fasta + ".det", 'r')
    lines = []
    for line in infile.readlines():
        newline = line.strip(" ").strip("\n")
        lines.append(newline)
    infile.close()
    # load the mip dictionary from file
    infile = open(primer3_output_DIR + input_fasta, 'r')
    hairpin_dic = json.load(infile)
    infile.close()
    # Lines following the line that starts with the mip name (PAIR..) has the energy information
    for line in lines:
        if line.startswith("PAIR"):
            next_index = (lines.index(line))+1
            # lines are formatted in as follows
            #  dG =      2.86  dH =     -0.10  dS =     -8.88  Tm = -261.9 â
            # remove ° (â) character from line.
            # space is used rather than a tab char, so I split at d, as in dG, dS etc.
            # Also change Tm to dTM to be able to use "d" delimiter.
            newline = lines[next_index].strip("\xe2\x84\x83").replace("Tm", "dTm").split("d")
            # once split, lines look like this: G =      2.86, so I slice the last characters
            # of the string that have the values
            dG = float(newline [1][-10:])
            dH = float(newline [2][-10:])
            dS = float(newline [3][-10:])
            TM = float(newline [4][-8:])
            # if this is the first hairpin information for this mip, add the hairpin information
            if not "HAIRPIN" in list(hairpin_dic["pair_information"][line]["mip_information"]["ref"].keys()):
                hairpin_dic["pair_information"][line]["mip_information"]["ref"]["HAIRPIN"]                = {'dG':dG, 'dH':dH, 'dS':dS, 'TM':TM}
            # if there was already another hairpin for this mip, replace it with this new hairpin
            # if it has a dG vaule that is smaller than the previous mip.
            elif float(hairpin_dic["pair_information"][line]["mip_information"]["ref"]                                  ["HAIRPIN"]['dG']) >= dG:
                hairpin_dic["pair_information"][line]["mip_information"]["ref"]["HAIRPIN"]                            = {'dG':dG, 'dH':dH, 'dS':dS, 'TM':TM}
    # print starting number of mips
    #print str(len(hairpin_dic["pair_information"])) + " mips analyzed!"
    # remove mips with dG<=dg at 60°C
    for mip in list(hairpin_dic["pair_information"].keys()):
        if "HAIRPIN" in list(hairpin_dic["pair_information"][mip]["mip_information"]["ref"].keys()):
            if float(hairpin_dic["pair_information"][mip]                                ["mip_information"]["ref"]["HAIRPIN"]['dG']) <= dg:
                hairpin_dic["pair_information"].pop(mip)
    # print number of mips with reasonable hairpin values
    #print str(len(hairpin_dic["pair_information"])) + " mips left with dG > " + dg + " @ 60°C"
    # create a fasta file for filtered mips to run bowtie2 alignments.
    outfile = open(bowtie2_input_DIR + output_json, 'w')
    for mip in list(hairpin_dic['pair_information'].keys()):
        mip_seq = hairpin_dic["pair_information"][mip]["mip_information"]["ref"]["SEQUENCE"]
        header = mip
        outline = ">" + mip + '\n' + mip_seq + '\n'
        outfile.write(outline)
    outfile.close()
    # save the new dictionary to file
    outfile = open(primer3_output_DIR + output_json, "w")
    json.dump(hairpin_dic, outfile, indent=1)
    outfile.close()
    return hairpin_dic
def hairpin_worker (l):
    file_locations = get_file_locations()
    script_DIR = file_locations["all"]["script_DIR"]
    input_fasta = l[0]
    mip_file = l[1]
    #dg = l[2]
    hairpin_tm = l[2]
    primer3_output_DIR = l[3]
    mfold_input_DIR = l[4]
    mfold_output_DIR = l[5]
    """ Check free energy values of secondary structures formed by mips @60°C.
    Run mfold_quiker script with mip condition parameters, get dG, dH, dS and TM
    values of each mip. TM values are a bit rough according to the author of the program.
    Therefore dG values @ 60°C is used to determine hairpin stability. All mips have
    a secondary structure with 0.82 kcal dG so 0.8 seem like a good cut off value to
    eliminate mips. This roughly means that any mip that is selected will require
    at least 0.8 kcal energy input to form a stable hairpin.
    Add HAIRPIN key is to the mip_information dict, which contains hairpin information:
    minimum dG, dH, dS and TM."""
    # run mfold_quiker script on mip sequences. current working directory is changed when
    # using subprocess.call because the input sequence must be in the directory where the
    # script is run.
    melt_output = subprocess.check_output(["melt.pl", "-n",  "DNA", "--temperature=25",                      "-N", "0.025", "-M", "0.01", script_DIR + mfold_input_DIR + input_fasta ],                     cwd=mfold_input_DIR)
    # a number of files are generated by the script using the input base name.
    # .det file has the dG TM information we need.
    out_list = melt_output.split("\n")[:-1]
    primer_names = []
    tm_values = []
    others = []
    for o in out_list:
        if "PAIR_" in o:
            primer_names.append(o.split(" ")[2][:-1])
        elif not o.startswith("dG"):
            tm_values.append(float(o.split("\t")[-1]))
        else:
            others.append(o)
    """
    with open(primer3_output_DIR + mip_file, 'r') as infile:
        hairpin_dic = json.load(infile)
    for i in xrange(len(primer_names)):
        hairpin_dic["pair_information"][primer_names[i]]["mip_information"]["ref"]["HAIRPIN"] = \
                tm_values[i]


    for mip in hairpin_dic["pair_information"].keys():
        if "HAIRPIN" in hairpin_dic["pair_information"][mip]["mip_information"]["ref"].keys():
            if hairpin_dic["pair_information"][mip]["mip_information"]["ref"]["HAIRPIN"] >= hairpin_tm:
                hairpin_dic["pair_information"].pop(mip)
        else:
            hairpin_dic["pair_information"].pop(mip)

    with open(mfold_output_DIR + input_fasta, 'w') as outfile:
        json.dump(hairpin_dic, outfile, indent=1)
    """
    output_tm = list(zip(primer_names, tm_values))
    with open(mfold_output_DIR + input_fasta, 'w') as outfile:
        outfile.write("\n".join(["\t".join(map(str, o))                                 for o in output_tm]))
    return 0
def make_chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]
def split_file(input_file,
               primer3_output_DIR,
               mfold_input_DIR,
               mfold_output_DIR,
               hairpin_tm,
               line_limit = 600):
    # open the mip fasta file
    # divide the file into smaller, 300 mip sized files
    map_input = []
    with open(mfold_input_DIR + input_file, 'r') as infile:
        counter = 0
        linenum = 0
        outlist = []
        for line in infile:
            linenum += 1
            outlist.append(line.strip())
            if linenum == line_limit:
                with open(mfold_input_DIR + input_file +                       str(counter),  'w') as outfile:
                    outfile.write("\n".join(outlist))
                temp = [input_file + str(counter), input_file, hairpin_tm,                 primer3_output_DIR, mfold_input_DIR, mfold_output_DIR]
                map_input.append(temp)
                counter += 1
                linenum = 0
                outlist = []
        if linenum != 0:
            with open(mfold_input_DIR + input_file +                       str(counter),  'w') as outfile:
                    outfile.write("\n".join(outlist))
            temp = [input_file + str(counter), input_file, hairpin_tm,             primer3_output_DIR, mfold_input_DIR, mfold_output_DIR]
            map_input.append(temp)
    return map_input
def hairpin_multi (input_dict,
                   input_file,
                   output_file, primer3_output_DIR, bowtie2_input_DIR,\
                   mfold_input_DIR, mfold_output_DIR, num_process=20,
                   hairpin_tm=60,
                  outp = 0):
    feeder = split_file(input_file,
                        primer3_output_DIR,
                        mfold_input_DIR,
                        mfold_output_DIR,
                        hairpin_tm)
    #print feeder
    dics = []
    dic = {"pair_information":{}, "sequence_information":{}}
    ################################################################################
    # the following line is commented out because when this script is imported into
    # another, __name__ is no longer __main__ and the multiprocessing does not work.
    # it should work like this in unix systems but may cause problems in windows.
    ################################################################################
    #if __name__ == "__main__":
    # create a pool of "num_process" processors
    p = Pool(num_process)
    p.map_async(hairpin_worker, feeder, callback=dics.append)
    p.close()
    p.join()
    p.join()
    """
        for f in feeder:
        hdic_file = mfold_output_DIR + f[0]
        with open(hdic_file) as hdf:
            d = json.load(hdf)
        dic["sequence_information"] = d["sequence_information"]
        dic["pair_information"].update(d["pair_information"])
        mip_counter += len(dic["pair_information"])
    """
    #print dics
    # check if the worker process was succesful
    if len(dics) == 0:
        # if no results have been generated, then function should return 1
        print("Hairpin analysis was not successful.")
        return 1
    """
    with open(primer3_output_DIR + input_file) as infile:
        hairpin_dict = json.load(infile)
    """
    hairpin_dict = input_dict
    hairpin = hairpin_dict["pair_information"]
    mip_counter = 0
    for f in feeder:
        tm_file = mfold_output_DIR + f[0]
        with open(tm_file) as infile:
            for line in infile:
                newline = line.strip().split("\t")
                mip_name = newline[0]
                mip_tm = float(newline[1])
                if mip_tm < hairpin_tm:
                    hairpin[mip_name]["mip_information"]                    ["ref"]["HAIRPIN"] = mip_tm
                    mip_counter += 1
                else:
                    hairpin.pop(mip_name)

    """
    dic["sequence_information"] = dics[0]["sequence_information"]
    # create a counter for mips remaining in the dictionaries
    mip_counter = 0
    for d in dics:
        dic["pair_information"].update(d["pair_information"])
        mip_counter += len(dic["pair_information"])
    """

    # check if there is at least 1 mip remaining
    if mip_counter == 0:
        # return if no mips left after processing for hairpin
        print("No mips left after hairpin analysis")
        return 1
    if outp:
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(hairpin_dict, outfile, indent=1)
    # create a fasta file for filtered mips to run bowtie2 alignments.
    with open(bowtie2_input_DIR + output_file, 'w') as outfile:
        for mip in list(hairpin.keys()):
            mip_seq = hairpin[mip]["mip_information"]["ref"]["SEQUENCE"]
            header = mip
            outline = ">" + mip + '\n' + mip_seq + '\n'
            outfile.write(outline)
    return hairpin_dict
def score_paralog_primers (primer_file,
                           output_file,
                           primer3_output_DIR,
                           ext,
                           mask_penalty,
                          species,
                          outp = 1):
    """ Score primers in a dictionary according to scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    # load extension and ligation primers from file
    """
    with open (primer3_output_DIR + primer_file, 'r') as infile:
        dic = json.load(infile)

    """
    dic = primer_file
    primers = dic["primer_information"]
    sequence = dic["sequence_information"]
    extension = (ext == "extension")
    # extract template sequence
    seq_template = sequence["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for p in primers:
        # get coordinates of primers in the form of "start_base, len"
        coord = primers[p]["COORDINATES"]
        strand = primers[p]["ORI"]
        if strand == "forward": # then the extension arm is the LEFT primer and ligation arm is RIGHT
            primer_end = int(coord.split(",")[0]) + int(coord.split(",")[1]) - 1
            # 1 is added or subtracted due to sequence index being zero based.
            next_bases = seq_template [(primer_end+1):(primer_end+3)]
        elif strand == "reverse": # then the extension arm is the RIGHT primer and ligation is LEFT
            primer_end = int(coord.split(",")[0]) - int(coord.split(",")[1]) + 1
            next_bases = reverse_complement(seq_template [(primer_end - 2):primer_end])
        # add "NEXT_BASES" key and its value to mip dictionary
        primers[p]["NEXT_BASES"] = next_bases
        #print next_bases
    # define scoring matrices
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    if species.startswith("pf"):
        for i in range(100):
            if i < 15:
                arm_gc_con[i] = worst
            elif i < 20:
                arm_gc_con[i] = low
            elif i < 25:
                arm_gc_con[i] = mid
            elif i < 60:
                arm_gc_con[i] = best
            elif i < 65:
                arm_gc_con[i] = mid
            elif i < 70:
                arm_gc_con[i] = low
            else :
                arm_gc_con[i] = worst
    else:
        for i in range(100):
            if i < 35:
                arm_gc_con[i] = worst
            elif i < 40:
                arm_gc_con[i] = low
            elif i < 45:
                arm_gc_con[i] = mid
            elif i < 60:
                arm_gc_con[i] = best
            elif i < 65:
                arm_gc_con[i] = mid
            elif i < 70:
                arm_gc_con[i] = low
            else :
                arm_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"":worst, "A":worst, "T":worst, "G":worst, "C": worst,
                  "GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # gc clamp matrix: a mid score for g or c clamp and best for extension gg or cc
    best = 200
    mid = 100
    worst = 0
    # using an empirical scoring dictionary here based on initial results
    ext_gc_clamp = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":worst, "CA":worst, "GT":worst, "CT":worst,
                  "AG":mid, "TG":mid, "AC":mid, "TC":mid,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    lig_gc_clamp = {"G": mid, "C": mid, "A": worst, "T": worst}

    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed to get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix [i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix [i] = best
        elif (30 > i > 28):
            extension_len_matrix [i] = low
        elif (i > 29):
            extension_len_matrix [i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (i == 19):
            ligation_len_matrix [i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix [i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix [i] = low
        elif (i > 30):
            ligation_len_matrix [i] = worst
    # score all mips using above matrices
    for p in list(primers.keys()):
        # get arm sequences
        seq = primers[p]["SEQUENCE"]
        # count lower case masked nucleotides
        mask_count = sum(-1 for n in seq if n.islower())
        """
        if mask_count > 2:
            mask_score = -3 * mask_penalty
        else:
            mask_score = mask_count * mask_penalty
        """
        mask_score = mask_count * mask_penalty
        # arm lengths
        if extension:
            len_score = extension_len_matrix[len(seq)]
        else:
            len_score = ligation_len_matrix[len(seq)]

        # gc clamp
        if extension:
            gc_clamp = ext_gc_clamp[seq[-2:].upper()]
        else:
            gc_clamp = lig_gc_clamp[seq[-1].upper()]
        # get gc percentages and convert to int.
        gc = int(float(primers[p]["GC_PERCENT"]))
        # get next base values
        next_b = primers[p]["NEXT_BASES"]
        all_scores = {"arm_len":[len(seq), len_score],                      "arm_gc":[gc, arm_gc_con[gc]],                      "mask_penalty":mask_penalty,                      "gc_clamp": gc_clamp,                      "next_bases":[next_b, next_bases[next_b.upper()]]
                      }
        # calculate total score
        score = len_score + arm_gc_con[gc] + mask_score + next_bases[next_b.upper()]
        # add score to dictionary
        primers[p]["SCORE"] = score
        primers[p]["all_scores"] = all_scores
    if outp:
        # write dictionary to json file
        outfile = open (primer3_output_DIR + output_file, "w")
        json.dump(dic, outfile, indent=1)
        outfile.close()
        #print "hello"
    return dic
def score_primers (ext_file, lig_file, ext_out, lig_out, primer3_output_DIR):
    """ Score primers in a dictionary according to scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    # load extension and ligation primers from file
    with open (primer3_output_DIR + ext_file, 'r') as infile:
        ext_dic = json.load(infile)
    with open (primer3_output_DIR + lig_file, 'r') as infile:
        lig_dic = json.load(infile)
    # extract template sequence
    ext_template = ext_dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    lig_template = lig_dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    ext_temp_lenght = len(ext_template)
    lig_temp_length = len(lig_template)
    # define scoring matrices
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    for i in range(20,81):
        if i < 35:
            arm_gc_con[i] = worst
        elif i < 40:
            arm_gc_con[i] = low
        elif i < 45:
            arm_gc_con[i] = mid
        elif i < 60:
                arm_gc_con[i] = best
        elif i < 65:
                arm_gc_con[i] = mid
        elif i < 70:
            arm_gc_con[i] = low
        else :
            arm_gc_con[i] = worst
    # gc clamp matrix: a mid score for g or c clamp and best for extension gg or cc
    best = 200
    mid = 100
    worst = 0
    ext_gc_clamp = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":worst, "CA":worst, "GT":worst, "CT":worst,
                  "AG":mid, "TG":mid, "AC":mid, "TC":mid,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    lig_gc_clamp = {"G": mid, "C": mid, "A": worst, "T": worst}

    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed to
    # get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix [i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix [i] = best
        elif (30 > i > 28):
            extension_len_matrix [i] = low
        elif (i > 29):
            extension_len_matrix [i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (i == 19):
            ligation_len_matrix [i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix [i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix [i] = low
        elif (i > 30):
            ligation_len_matrix [i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # score extension arm primers using above matrices
    for ext_arm in list(ext_dic["primer_information"].keys()):
        # extract primer sequence
        extension_seq = ext_dic["primer_information"][ext_arm]["SEQUENCE"]
        # extract primer coordinates
        ext_coord = (ext_dic["primer_information"][ext_arm]["COORDINATES"]).split(",")
        ext_start = int(ext_coord[0])
        ext_len = int(ext_coord[1])
        # get gc percentages and convert to int.
        extension_gc = int(float(ext_dic["primer_information"][ext_arm]["GC_PERCENT"]))
        # score gc clamp
        gc_clamp = ext_gc_clamp[extension_seq[-2:].upper()]
        # get primer orientation
        if ext_arm.startswith("PRIMER_LEFT"):
            ext_ori = "forward"
        elif ext_arm.startswith("PRIMER_RIGHT"):
            ext_ori = "reverse"
        # get downstream 2 bases for each primer
        if ext_ori == "forward":
            extension_end = ext_start + ext_len - 1
            # 1 is subtracted due to sequence index being zero based.
            ext_next = ext_template [(extension_end+1):(extension_end+3)]
        else:
            extension_end = ext_start - ext_len + 1
            ext_next = reverse_complement(ext_template [(extension_end - 2):extension_end])
        # add next base tag:value pair to dictionary
        ext_dic["primer_information"][ext_arm]["NEXT_BASES"] = ext_next
        if len(ext_next) < 2:
            score = (extension_len_matrix[ext_len] + arm_gc_con[extension_gc] + gc_clamp)
        else:
            score = (extension_len_matrix[ext_len] + arm_gc_con[extension_gc] +                      next_bases[ext_next.upper()] + gc_clamp)
        # add score to dictionary
        ext_dic["primer_information"][ext_arm]["SCORE"] = score
        # add primer orientation to dic
        ext_dic["primer_information"][ext_arm]["ORIENTATION"] = ext_ori
    # score ligation arm primers using above matrices
    for lig_arm in list(lig_dic["primer_information"].keys()):
        # extract primer sequence
        ligation_seq = lig_dic["primer_information"][lig_arm]["SEQUENCE"]
        # extract primer coordinates
        lig_coord = (lig_dic["primer_information"][lig_arm]["COORDINATES"]).split(",")
        lig_start = int(lig_coord[0])
        lig_len = int(lig_coord[1])
        # get gc percentages and convert to int.
        ligation_gc = int(float(lig_dic["primer_information"][lig_arm]["GC_PERCENT"]))
        # score gc clamp
        gc_clamp = lig_gc_clamp[ligation_seq[-1].upper()]
        # get primer orientation
        if lig_arm.startswith("PRIMER_LEFT"):
            lig_ori = "forward"
        elif lig_arm.startswith("PRIMER_RIGHT"):
            lig_ori = "reverse"
        # get downstream 2 bases for each primer
        if lig_ori == "forward":
            ligation_end = lig_start + lig_len - 1
            # 1 is subtracted due to sequence index being zero based.
            lig_next = lig_template [(ligation_end+1):(ligation_end+3)]
        else:
            ligation_end = lig_start - lig_len + 1
            lig_next = reverse_complement(lig_template [(ligation_end - 2):ligation_end])
        # add next base tag:value pair to dictionary
        lig_dic["primer_information"][lig_arm]["NEXT_BASES"] = lig_next
        if len(lig_next) < 2:
            score = (ligation_len_matrix[lig_len] + arm_gc_con[ligation_gc] + gc_clamp)
        else:
            score = (ligation_len_matrix[lig_len] + arm_gc_con[ligation_gc] +                      next_bases[lig_next.upper()] + gc_clamp)
        # add score to dictionary
        lig_dic["primer_information"][lig_arm]["SCORE"] = score
        # add primer orientation to dic
        lig_dic["primer_information"][lig_arm]["ORIENTATION"] = lig_ori
    # write dictionaries to json file
    with open (primer3_output_DIR + ext_out, "w") as outfile:
     json.dump(ext_dic, outfile, indent=1)
    with open (primer3_output_DIR + lig_out, "w") as outfile:
     json.dump(lig_dic, outfile, indent=1)
    return
def filter_primers (primer_file, output_file,
                    primer3_output_DIR, n, bin_size,
                   outp = 1):
    """ filter primers so that only top n scoring primers
    ending within the same subregion (determined by bin_size)
    remains."""
    # load extension and ligation primers from file
    """
    with open (primer3_output_DIR + primer_file, 'r') as infile:
        dic = json.load(infile)
    """
    dic = primer_file
    template_seq = dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    template_len = len(template_seq)
    forward_bins = {}
    reverse_bins = {}
    best_score = 0
    for i in range(template_len//bin_size + 1):
        forward_bins[i] = []
        reverse_bins[i] = []
    for primer in list(dic["primer_information"].keys()):
        # get primer orientation
        ori = dic["primer_information"][primer]["ORI"]
        # get primer start coordinate
        start = int(dic["primer_information"][primer]["COORDINATES"].split(",")[0])
        primer_len = int(dic["primer_information"][primer]["COORDINATES"].split(",")[1])
        if ori == "forward":
            end = start + primer_len - 1
        elif ori == "reverse":
            end = start - primer_len + 1
        # which bin the start coordinate falls into
        end_bin = end//bin_size
        # get primer score
        score = dic["primer_information"][primer]["SCORE"]
        # append the primer name/score to appropriate bin dic
        if ori == "forward":
            forward_bins[end_bin].append([primer, score])
        elif ori == "reverse":
            reverse_bins[end_bin].append([primer, score])
    best_dic = {}
    best_dic["sequence_information"] = dic["sequence_information"]
    best_dic["primer_information"] = {}
    # find best scoring mips in each forward bin
    for key in forward_bins:
        # sort primers for score
        primer_set = sorted(forward_bins[key], key= itemgetter(1))
        # get best scoring primers (all primers if there are less than n)
        if len(primer_set) < n:
            best_primers = primer_set
        else:
            best_primers = primer_set[-n:]
        # add best primers do dictionary
        for primers in best_primers:
            primer_name = primers[0]
            best_dic["primer_information"][primer_name] = dic["primer_information"][primer_name]
    # find best scoring mips in each reverse bin
    for key in reverse_bins:
        # sort primers for score
        primer_set = sorted(reverse_bins[key], key= itemgetter(1))
        # get best scoring primers (all primers if there are less than n)
        if len(primer_set) < n:
            best_primers = primer_set
        else:
            best_primers = primer_set[-n:]
        # add best primers do dictionary
        for primers in best_primers:
            primer_name = primers[0]
            best_dic["primer_information"][primer_name] = dic["primer_information"][primer_name]
    # write new dic to file
    if outp:
        with open (primer3_output_DIR + output_file, "w") as outfile:
            json.dump(best_dic, outfile, indent=1)
    return best_dic
def filter_mips (mip_dic, bin_size, mip_limit):
    """
    Filter mips in "mip_dic" so that only top scoring mip
    ending within the "bin_size" nucleotides on the same
    strand remain.
    """
    # load extension and ligation primers from file
    shuffled = list(mip_dic.keys())
    random.shuffle(shuffled)
    for m in shuffled:
        if len(mip_dic) <= mip_limit:
            return
        try:
            found = 0
            m_start = mip_dic[m].mip["C0"]["capture_start"]
            m_end = mip_dic[m].mip["C0"]["capture_end"]
            m_func = mip_dic[m].func_score
            m_tech = mip_dic[m].tech_score
            m_ori = mip_dic[m].mip["C0"]["orientation"]
            #shuffled_keys = list(mip_dic.keys())
            #random.shuffle(shuffled_keys)
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
                        if ((abs(n_start - m_start) <= bin_size) or                             (abs(n_end - m_end) <= bin_size)) and                             (m_ori == n_ori):
                            if (m_tech + m_func) >= (n_tech + n_func):
                                mip_dic.pop(n)
                            else:
                                mip_dic.pop(m)
                                break
                except KeyError as e:
                    continue
        except KeyError as e:
            continue
    return
def remove_mips (mip_dic):
    """ filter primers so that only top n scoring primers
    ending within the same subregion (determined by bin_size)
    remains."""
    # load extension and ligation primers from file
    shuffled = list(mip_dic.keys())
    random.shuffle(shuffled)
    for m in shuffled:
        try:
            found = 0
            m_cap_start = mip_dic[m].mip["C0"]["capture_start"]
            m_cap_end = mip_dic[m].mip["C0"]["capture_end"]
            m_start = mip_dic[m].mip["C0"]["mip_start"]
            m_end = mip_dic[m].mip["C0"]["mip_end"]
            m_func = mip_dic[m].func_score
            m_tech = mip_dic[m].tech_score
            m_ori = mip_dic[m].mip["C0"]["orientation"]
            m_score = m_func + m_tech
            for n in shuffled:
                try:
                    if mip_dic[m].name != mip_dic[n].name:
                        n_cap_start = mip_dic[n].mip["C0"]["capture_start"]
                        n_cap_end = mip_dic[n].mip["C0"]["capture_end"]
                        n_start = mip_dic[n].mip["C0"]["mip_start"]
                        n_end = mip_dic[n].mip["C0"]["mip_end"]
                        n_func = mip_dic[n].func_score
                        n_tech = mip_dic[n].tech_score
                        n_ori = mip_dic[n].mip["C0"]["orientation"]
                        n_score = n_func + n_tech
                        remove = False
                        if (m_start <= n_start <= m_end) or                             (n_start <= m_start <= n_end):
                            if m_ori == n_ori:
                                remove = True
                            elif ((n_start < m_start) and                                 (n_cap_start < m_start) and                                 (m_cap_start < n_cap_end) and                                 (n_end < m_cap_end)) or                                 ((m_start < n_start) and                                 (m_cap_start < n_start) and                                 (n_cap_start < m_cap_end) and                                 (m_end < n_cap_end)):
                                remove = False
                            else:
                                remove = True
                        if remove:
                            if m_score >= n_score:
                                mip_dic.pop(n)
                            else:
                                mip_dic.pop(m)
                                break
                except KeyError as e:
                    continue
        except KeyError as e:
            continue
    return
def score_mips (mip_file, primer3_output_DIR, output_file):
    """ Score mips in a dictionary according to a scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    # open mip dictionary from file
    with open (primer3_output_DIR + mip_file, 'r') as infile:
        dic = json.load(infile)
    # add "NEXT_BASES" tag:value pair to dictionary.
    # next_bases are 2 bases immediately downstream of extension primer and ligation primer
    # extract template sequence
    seq_template = dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for mip in dic["pair_information"]:
        # get coordinates of primers in the form of "start_base, len"
        extension_coord = dic["pair_information"][mip]["extension_primer_information"]["COORDINATES"]
        ligation_coord = dic["pair_information"][mip]["ligation_primer_information"]["COORDINATES"]
        # orientation of the mip is used to determine if extension arm or ligation arm is originating
        # from PRIMER_LEFT or PRIMER_RIGTH. When an arm originates from LEFT_PRIMER, it is on the plus
        # strand of DNA and its length is added to the start coordinate to get the end coordinate,
        # while it is subtracted for RIGHT_PRIMERs
        strand = dic["pair_information"][mip]["orientation"]
        if strand == "forward": # then the extension arm is the LEFT primer and ligation arm is RIGHT
            extension_end = int(extension_coord.split(",")[0]) + int(extension_coord.split(",")[1]) - 1
            # 1 is added or subtracted due to sequence index being zero based.
            ligation_end = int(ligation_coord.split(",")[0]) - int(ligation_coord.split(",")[1]) + 1
            ext_next = seq_template [(extension_end+1):(extension_end+3)]
            lig_next = reverse_complement(seq_template [(ligation_end - 2):ligation_end])
        elif strand == "reverse": # then the extension arm is the RIGHT primer and ligation is LEFT
            extension_end = int(extension_coord.split(",")[0]) - int(extension_coord.split(",")[1]) + 1
            ligation_end = int(ligation_coord.split(",")[0]) + int(ligation_coord.split(",")[1]) - 1
            ext_next = reverse_complement(seq_template [(extension_end - 2):extension_end])
            lig_next = seq_template [(ligation_end+1):(ligation_end+3)]
        # add "NEXT_BASES" key and its value to mip dictionary
        dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"] = ext_next
        dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"] = lig_next
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    for i in range(100):
        if i < 35:
            arm_gc_con[i] = worst
        elif i < 40:
            arm_gc_con[i] = low
        elif i < 45:
            arm_gc_con[i] = mid
        elif i < 60:
                arm_gc_con[i] = best
        elif i < 65:
                arm_gc_con[i] = mid
        elif i < 70:
            arm_gc_con[i] = low
        else :
            arm_gc_con[i] = worst
    # capture gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 10000
    mid = 3000
    low = 100
    worse = 0
    worst = -5000
    # define matrix
    cap_gc_con = {}
    for i in range(100):
        if i<20:
            cap_gc_con[i] = worst
        elif i < 35:
            cap_gc_con[i] = worse
        elif i < 40:
            cap_gc_con[i] = low
        elif i < 45:
            cap_gc_con[i] = mid
        elif i < 55:
                cap_gc_con[i] = best
        elif i < 60:
                cap_gc_con[i] = mid
        elif i < 65:
            cap_gc_con[i] = low
        elif i < 80:
            cap_gc_con[i] = worse
        else :
            cap_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # gc clamp matrix: a mid score for g or c clamp and best for extension gg or cc
    best = 200
    mid = 100
    worst = 0
    ext_gc_clamp = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":worst, "CA":worst, "GT":worst, "CT":worst,
                  "AG":mid, "TG":mid, "AC":mid, "TC":mid,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    lig_gc_clamp = {"G": mid, "C": mid, "A": worst, "T": worst}

    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed to get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix [i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix [i] = best
        elif (30 > i > 28):
            extension_len_matrix [i] = low
        elif (i > 29):
            extension_len_matrix [i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (i == 19):
            ligation_len_matrix [i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix [i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix [i] = low
        elif (i > 30):
            ligation_len_matrix [i] = worst
    # score all mips using above matrices
    for mip in list(dic["pair_information"].keys()):
        # get arm sequences
        ligation_seq = dic["pair_information"][mip]["ligation_primer_information"]["SEQUENCE"]
        extension_seq = dic["pair_information"][mip]["extension_primer_information"]["SEQUENCE"]
        # count lower case masked nucleotides
        ligation_mask_penalty = sum(-1000 for n in ligation_seq if n.islower())
        ligation_mask_penalty += sum(-5000 for n in ligation_seq[-5:] if n.islower())
        extension_mask_penalty = sum(-1000 for n in extension_seq if n.islower())
        extension_mask_penalty += sum(-5000 for n in extension_seq[-5:] if n.islower())
        # arm lengths
        ligation_len = len(ligation_seq)
        extension_len = len(extension_seq)
        # find capture gc content
        ligation_start = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_START"])
        ligation_end = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_END"])
        extension_start = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_START"])
        extension_end = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_END"])
        chrom = dic["pair_information"][mip]["extension_primer_information"]["CHR"]
        mip_coord = sorted([ligation_start, ligation_end, extension_end, extension_start])
        capture_key = chrom + ":" + str(mip_coord[1]) + "-" + str(mip_coord[2])
        capture_seq = fasta_to_sequence(get_fasta(capture_key, species="hs"))
        capture_gc = calculate_gc(capture_seq)
        # gc clamp
        gc_clamp = ext_gc_clamp[extension_seq[-2:].upper()] + lig_gc_clamp[ligation_seq[-1]]
        # get gc percentages and convert to int.
        extension_gc = int(float(dic["pair_information"][mip]["extension_primer_information"]["GC_PERCENT"]))
        ligation_gc = int(float(dic["pair_information"][mip]["ligation_primer_information"]["GC_PERCENT"]))
        # get next base values
        extension_next = dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"]
        ligation_next = dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"]
        # desired/undesired copies captured
        copy_bonus = 0
        if "AMP_PARA" in list(dic["pair_information"][mip]["extension_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            ext_copies = dic["pair_information"][mip]["extension_primer_information"]["AMP_PARA"]
            pairs = dic["pair_information"][mip]["pairs"]
        else:
            ext_copies = []
        if "AMP_PARA" in list(dic["pair_information"][mip]["ligation_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            lig_copies = dic["pair_information"][mip]["ligation_primer_information"]["AMP_PARA"]
        else:
            lig_copies = []
        wanted = []
        unwanted = []
        for ec in ext_copies:
            if ec in lig_copies:
                if (ec in desired_copies) and (pairs[ec][-1]):
                    copy_bonus += 1000
                    wanted.append(ec)
                else:
                    copy_bonus -= 500
                    unwanted.append(ec)

        all_scores = {"extension_arm_len":[extension_len, extension_len_matrix[extension_len]],                      "ligation_arm_len":[ligation_len, ligation_len_matrix[ligation_len]],                      "extension_arm_gc":[extension_gc, arm_gc_con[extension_gc]],                      "ligation_arm_gc":[ligation_gc, arm_gc_con[ligation_gc]],                      "ligation_mask_penalty":ligation_mask_penalty,                      "extension_mask_penalty":extension_mask_penalty,                      "capture_gc_content":[capture_gc, cap_gc_con[capture_gc]],                      "copy_bonus": {"intended_copies": wanted, "unintended_copies": unwanted,                                     "desired_copies": desired_copies,"bonus": copy_bonus}, "gc_clamp": gc_clamp,                      "extension_next_bases":[extension_next, next_bases[extension_next.upper()]],                      "ligation_next_bases":[ligation_next, next_bases[ligation_next.upper()]],                         }
        # calculate total score
        score = (extension_len_matrix[extension_len] + ligation_len_matrix[ligation_len] +
                 arm_gc_con[extension_gc] + arm_gc_con[ligation_gc] +
                 next_bases[extension_next.upper()] + next_bases[ligation_next.upper()] +
                 gc_clamp + ligation_mask_penalty + extension_mask_penalty + cap_gc_con[capture_gc] + \
                 copy_bonus)
        # add mip_score to dictionary
        dic["pair_information"][mip]["mip_information"]["mip_score"] = score
        dic["pair_information"][mip]["mip_information"]["all_scores"] = all_scores
        dic["pair_information"][mip]["mip_information"]["capture_seq"] = capture_seq
    # write dictionary to json file
    outfile = open (primer3_output_DIR + output_file, "w")
    json.dump(dic, outfile, indent=1)
    outfile.close()
    return dic
def add_captures (mip_dic, target_dic):
    for mip in list(mip_dic["pair_information"].keys()):
        d = mip_dic["pair_information"][mip]
        # find which diffs are captured by mip
        # extract mip's coortinates
        ext_start = d["extension_primer_information"]["GENOMIC_START"]
        ext_end = d["extension_primer_information"]["GENOMIC_END"]
        lig_start = d["ligation_primer_information"]["GENOMIC_START"]
        lig_end = d["ligation_primer_information"]["GENOMIC_END"]
        coord = [ext_start, ext_end, lig_start, lig_end].sort()
        mip_start = coord[1]
        mip_end = coord[2]
        # create a dictionary for the targets the mip captures
        captured_snps = {}
        for snp in target_snps:
            if mip_end >= target_snps[snp]["begin"] >= mip_start:
                captured_snps[snp] = target_snps[snp]
        # add captured diffs information to mip dictionary
        d["mip_information"]["captured_diffs"] = captured_snps
    return mip_dic
def add_paralog_info(mip_dic, num_para):
    for mip in list(mip_dic["pair_information"].keys()):
        # extract the captured diffs from the mip_dic
        caps = mip_dic["pair_information"][mip]["mip_information"]["captured_diffs"]
        # create a dic for how many paralogs a mip captures
        mip_caps = {}
        # populate the dic with 0 for each paralog key
        for i in range(num_para):
            mip_caps[i] = 0
        # find the diffs captured that identify a paralog
        for diff in caps:
            # paralog specific mip sources are exonic and filtered diffs only
            if (caps[diff]["source"] == "exonic_diffs") or                (caps[diff]["source"] == "filtered_diffs"):
                # the diff is in the form a:1067:CTA:0,1,2 and we want to
                # analyze the last part of it for paralogs identified
                mip_diff = (caps[diff]["diff"].split(":")[-1]).split(",")
                # add the paralogs identified to mip_set paralog caps
                for j in mip_diff:
                    mip_caps[int(j)] += 1
        # calculate how many paralogs identified by the mip
        mip_para = 0
        # for each paralog
        for k in mip_caps:
            # if at least one of the diffs identifies the paralog
            if mip_caps[k] > 0:
                # increment captured paralog number by 1
                mip_para += 1
        # add this information to mip dic
        mip_dic["pair_information"][mip]["mip_information"]["captured_paralog_number"] = mip_para
    return mip_dic
def score_mip_objects (mip_object):
    """ Score mips in a dictionary according to a scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    # open mip dictionary from file
    infile = open (primer3_output_DIR + mip_file, 'r')
    dic = json.load(infile)
    infile.close()
    # add "NEXT_BASES" tag:value pair to dictionary.
    # next_bases are 2 bases immediately downstream of extension primer and ligation primer
    # extract template sequence
    seq_template = dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for mip in dic["pair_information"]:
        # get coordinates of primers in the form of "start_base, len"
        extension_coord = dic["pair_information"][mip]["extension_primer_information"]["COORDINATES"]
        ligation_coord = dic["pair_information"][mip]["ligation_primer_information"]["COORDINATES"]
        # orientation of the mip is used to determine if extension arm or ligation arm is originating
        # from PRIMER_LEFT or PRIMER_RIGTH. When an arm originates from LEFT_PRIMER, it is on the plus
        # strand of DNA and its length is added to the start coordinate to get the end coordinate,
        # while it is subtracted for RIGHT_PRIMERs
        strand = dic["pair_information"][mip]["orientation"]
        if strand == "forward": # then the extension arm is the LEFT primer and ligation arm is RIGHT
            extension_end = int(extension_coord.split(",")[0]) + int(extension_coord.split(",")[1]) - 1
            # 1 is added or subtracted due to sequence index being zero based.
            ligation_end = int(ligation_coord.split(",")[0]) - int(ligation_coord.split(",")[1]) + 1
            ext_next = seq_template [(extension_end+1):(extension_end+3)]
            lig_next = reverse_complement(seq_template [(ligation_end - 2):ligation_end])
        elif strand == "reverse": # then the extension arm is the RIGHT primer and ligation is LEFT
            extension_end = int(extension_coord.split(",")[0]) - int(extension_coord.split(",")[1]) + 1
            ligation_end = int(ligation_coord.split(",")[0]) + int(ligation_coord.split(",")[1]) - 1
            ext_next = reverse_complement(seq_template [(extension_end - 2):extension_end])
            lig_next = seq_template [(ligation_end+1):(ligation_end+3)]
        # add "NEXT_BASES" key and its value to mip dictionary
        dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"] = ext_next
        dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"] = lig_next
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    for i in range(100):
        if i < 35:
            arm_gc_con[i] = worst
        elif i < 40:
            arm_gc_con[i] = low
        elif i < 45:
            arm_gc_con[i] = mid
        elif i < 60:
                arm_gc_con[i] = best
        elif i < 65:
                arm_gc_con[i] = mid
        elif i < 70:
            arm_gc_con[i] = low
        else :
            arm_gc_con[i] = worst
    # capture gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 10000
    mid = 3000
    low = 100
    worse = 0
    worst = -5000
    # define matrix
    cap_gc_con = {}
    for i in range(100):
        if i<20:
            cap_gc_con[i] = worst
        elif i < 35:
            cap_gc_con[i] = worse
        elif i < 40:
            cap_gc_con[i] = low
        elif i < 45:
            cap_gc_con[i] = mid
        elif i < 55:
                cap_gc_con[i] = best
        elif i < 60:
                cap_gc_con[i] = mid
        elif i < 65:
            cap_gc_con[i] = low
        elif i < 80:
            cap_gc_con[i] = worse
        else :
            cap_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # gc clamp matrix: a mid score for g or c clamp and best for extension gg or cc
    best = 200
    mid = 100
    worst = 0
    ext_gc_clamp = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":worst, "CA":worst, "GT":worst, "CT":worst,
                  "AG":mid, "TG":mid, "AC":mid, "TC":mid,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    lig_gc_clamp = {"G": mid, "C": mid, "A": worst, "T": worst}

    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed to get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix [i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix [i] = best
        elif (30 > i > 28):
            extension_len_matrix [i] = low
        elif (i > 29):
            extension_len_matrix [i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (i == 19):
            ligation_len_matrix [i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix [i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix [i] = low
        elif (i > 30):
            ligation_len_matrix [i] = worst
    # score all mips using above matrices
    for mip in list(dic["pair_information"].keys()):
        # get arm sequences
        ligation_seq = dic["pair_information"][mip]["ligation_primer_information"]["SEQUENCE"]
        extension_seq = dic["pair_information"][mip]["extension_primer_information"]["SEQUENCE"]
        # count lower case masked nucleotides
        ligation_mask_penalty = sum(-1000 for n in ligation_seq if n.islower())
        ligation_mask_penalty += sum(-5000 for n in ligation_seq[-5:] if n.islower())
        extension_mask_penalty = sum(-1000 for n in extension_seq if n.islower())
        extension_mask_penalty += sum(-5000 for n in extension_seq[-5:] if n.islower())
        # arm lengths
        ligation_len = len(ligation_seq)
        extension_len = len(extension_seq)
        # find capture gc content
        ligation_start = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_START"])
        ligation_end = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_END"])
        extension_start = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_START"])
        extension_end = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_END"])
        chrom = dic["pair_information"][mip]["extension_primer_information"]["CHR"]
        mip_coord = sorted([ligation_start, ligation_end, extension_end, extension_start])
        capture_key = chrom + ":" + str(mip_coord[1]) + "-" + str(mip_coord[2])
        capture_seq = fasta_to_sequence(get_fasta(capture_key, species="hs"))
        capture_gc = calculate_gc(capture_seq)
        # gc clamp
        gc_clamp = ext_gc_clamp[extension_seq[-2:].upper()] + lig_gc_clamp[ligation_seq[-1]]
        # get gc percentages and convert to int.
        extension_gc = int(float(dic["pair_information"][mip]["extension_primer_information"]["GC_PERCENT"]))
        ligation_gc = int(float(dic["pair_information"][mip]["ligation_primer_information"]["GC_PERCENT"]))
        # get next base values
        extension_next = dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"]
        ligation_next = dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"]
        # desired/undesired copies captured
        copy_bonus = 0
        if "AMP_PARA" in list(dic["pair_information"][mip]["extension_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            ext_copies = dic["pair_information"][mip]["extension_primer_information"]["AMP_PARA"]
            pairs = dic["pair_information"][mip]["pairs"]
        else:
            ext_copies = []
        if "AMP_PARA" in list(dic["pair_information"][mip]["ligation_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            lig_copies = dic["pair_information"][mip]["ligation_primer_information"]["AMP_PARA"]
        else:
            lig_copies = []
        wanted = []
        unwanted = []
        for ec in ext_copies:
            if ec in lig_copies:
                if (ec in desired_copies) and (pairs[ec][-1]):
                    copy_bonus += 1000
                    wanted.append(ec)
                else:
                    copy_bonus -= 500
                    unwanted.append(ec)

        all_scores = {"extension_arm_len":[extension_len, extension_len_matrix[extension_len]],                      "ligation_arm_len":[ligation_len, ligation_len_matrix[ligation_len]],                      "extension_arm_gc":[extension_gc, arm_gc_con[extension_gc]],                      "ligation_arm_gc":[ligation_gc, arm_gc_con[ligation_gc]],                      "ligation_mask_penalty":ligation_mask_penalty,                      "extension_mask_penalty":extension_mask_penalty,                      "capture_gc_content":[capture_gc, cap_gc_con[capture_gc]],                      "copy_bonus": {"intended_copies": wanted, "unintended_copies": unwanted,                                     "desired_copies": desired_copies,"bonus": copy_bonus}, "gc_clamp": gc_clamp,                      "extension_next_bases":[extension_next, next_bases[extension_next.upper()]],                      "ligation_next_bases":[ligation_next, next_bases[ligation_next.upper()]],                         }
        # calculate total score
        score = (extension_len_matrix[extension_len] + ligation_len_matrix[ligation_len] +
                 arm_gc_con[extension_gc] + arm_gc_con[ligation_gc] +
                 next_bases[extension_next.upper()] + next_bases[ligation_next.upper()] +
                 gc_clamp + ligation_mask_penalty + extension_mask_penalty + cap_gc_con[capture_gc] + \
                 copy_bonus)
        # add mip_score to dictionary
        dic["pair_information"][mip]["mip_information"]["mip_score"] = score
        dic["pair_information"][mip]["mip_information"]["all_scores"] = all_scores
        dic["pair_information"][mip]["mip_information"]["capture_seq"] = capture_seq
    # write dictionary to json file
    outfile = open (primer3_output_DIR + output_file, "w")
    json.dump(dic, outfile, indent=1)
    outfile.close()
    return dic
def score_mips_hla (mip_file, primer3_output_DIR, output_file, desired_copies=[]):
    """ Score mips in a dictionary according to a scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    # open mip dictionary from file
    infile = open (primer3_output_DIR + mip_file, 'r')
    dic = json.load(infile)
    infile.close()
    # add "NEXT_BASES" tag:value pair to dictionary.
    # next_bases are 2 bases immediately downstream of extension primer and ligation primer
    # extract template sequence
    seq_template = dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for mip in dic["pair_information"]:
        # get coordinates of primers in the form of "start_base, len"
        extension_coord = dic["pair_information"][mip]["extension_primer_information"]["COORDINATES"]
        ligation_coord = dic["pair_information"][mip]["ligation_primer_information"]["COORDINATES"]
        # orientation of the mip is used to determine if extension arm or ligation arm is originating
        # from PRIMER_LEFT or PRIMER_RIGTH. When an arm originates from LEFT_PRIMER, it is on the plus
        # strand of DNA and its length is added to the start coordinate to get the end coordinate,
        # while it is subtracted for RIGHT_PRIMERs
        strand = dic["pair_information"][mip]["orientation"]
        if strand == "forward": # then the extension arm is the LEFT primer and ligation arm is RIGHT
            extension_end = int(extension_coord.split(",")[0]) + int(extension_coord.split(",")[1]) - 1
            # 1 is added or subtracted due to sequence index being zero based.
            ligation_end = int(ligation_coord.split(",")[0]) - int(ligation_coord.split(",")[1]) + 1
            ext_next = seq_template [(extension_end+1):(extension_end+3)]
            lig_next = reverse_complement(seq_template [(ligation_end - 2):ligation_end])
        elif strand == "reverse": # then the extension arm is the RIGHT primer and ligation is LEFT
            extension_end = int(extension_coord.split(",")[0]) - int(extension_coord.split(",")[1]) + 1
            ligation_end = int(ligation_coord.split(",")[0]) + int(ligation_coord.split(",")[1]) - 1
            ext_next = reverse_complement(seq_template [(extension_end - 2):extension_end])
            lig_next = seq_template [(ligation_end+1):(ligation_end+3)]
        # add "NEXT_BASES" key and its value to mip dictionary
        dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"] = ext_next
        dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"] = lig_next
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    for i in range(100):
        if i < 35:
            arm_gc_con[i] = worst
        elif i < 40:
            arm_gc_con[i] = low
        elif i < 45:
            arm_gc_con[i] = mid
        elif i < 60:
                arm_gc_con[i] = best
        elif i < 65:
                arm_gc_con[i] = mid
        elif i < 70:
            arm_gc_con[i] = low
        else :
            arm_gc_con[i] = worst
    # capture gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 10000
    mid = 3000
    low = 100
    worse = 0
    worst = -5000
    # define matrix
    cap_gc_con = {}
    for i in range(100):
        if i<20:
            cap_gc_con[i] = worst
        elif i < 35:
            cap_gc_con[i] = worse
        elif i < 40:
            cap_gc_con[i] = low
        elif i < 45:
            cap_gc_con[i] = mid
        elif i < 55:
                cap_gc_con[i] = best
        elif i < 60:
                cap_gc_con[i] = mid
        elif i < 65:
            cap_gc_con[i] = low
        elif i < 80:
            cap_gc_con[i] = worse
        else :
            cap_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # gc clamp matrix: a mid score for g or c clamp and best for extension gg or cc
    best = 200
    mid = 100
    worst = 0
    ext_gc_clamp = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":worst, "CA":worst, "GT":worst, "CT":worst,
                  "AG":mid, "TG":mid, "AC":mid, "TC":mid,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    lig_gc_clamp = {"G": mid, "C": mid, "A": worst, "T": worst}

    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed to get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix [i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix [i] = best
        elif (30 > i > 28):
            extension_len_matrix [i] = low
        elif (i > 29):
            extension_len_matrix [i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18,36):
        if (i == 18) or (i == 19):
            ligation_len_matrix [i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix [i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix [i] = low
        elif (i > 30):
            ligation_len_matrix [i] = worst
    # score all mips using above matrices
    for mip in list(dic["pair_information"].keys()):
        # get arm sequences
        ligation_seq = dic["pair_information"][mip]["ligation_primer_information"]["SEQUENCE"]
        extension_seq = dic["pair_information"][mip]["extension_primer_information"]["SEQUENCE"]
        # count lower case masked nucleotides
        ligation_mask_penalty = sum(-1000 for n in ligation_seq if n.islower())
        ligation_mask_penalty += sum(-5000 for n in ligation_seq[-5:] if n.islower())
        extension_mask_penalty = sum(-1000 for n in extension_seq if n.islower())
        extension_mask_penalty += sum(-5000 for n in extension_seq[-5:] if n.islower())
        # arm lengths
        ligation_len = len(ligation_seq)
        extension_len = len(extension_seq)
        # find capture gc content
        ligation_start = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_START"])
        ligation_end = int(dic["pair_information"][mip]["ligation_primer_information"]["GENOMIC_END"])
        extension_start = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_START"])
        extension_end = int(dic["pair_information"][mip]["extension_primer_information"]["GENOMIC_END"])
        chrom = dic["pair_information"][mip]["extension_primer_information"]["CHR"]
        mip_coord = sorted([ligation_start, ligation_end, extension_end, extension_start])
        capture_key = chrom + ":" + str(mip_coord[1]) + "-" + str(mip_coord[2])
        capture_seq = fasta_to_sequence(get_fasta(capture_key, species="hs"))
        capture_gc = calculate_gc(capture_seq)
        # gc clamp
        gc_clamp = ext_gc_clamp[extension_seq[-2:].upper()] + lig_gc_clamp[ligation_seq[-1]]
        # get gc percentages and convert to int.
        extension_gc = int(float(dic["pair_information"][mip]["extension_primer_information"]["GC_PERCENT"]))
        ligation_gc = int(float(dic["pair_information"][mip]["ligation_primer_information"]["GC_PERCENT"]))
        # get next base values
        extension_next = dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"]
        ligation_next = dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"]
        # desired/undesired copies captured
        copy_bonus = 0
        if "AMP_PARA" in list(dic["pair_information"][mip]["extension_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            ext_copies = dic["pair_information"][mip]["extension_primer_information"]["AMP_PARA"]
            pairs = dic["pair_information"][mip]["pairs"]
        else:
            ext_copies = []
        if "AMP_PARA" in list(dic["pair_information"][mip]["ligation_primer_information"].keys())            and "pairs" in list(dic["pair_information"][mip].keys()):
            lig_copies = dic["pair_information"][mip]["ligation_primer_information"]["AMP_PARA"]
        else:
            lig_copies = []
        wanted = []
        unwanted = []
        for ec in ext_copies:
            if ec in lig_copies:
                if (ec in desired_copies) and (pairs[ec][-1]):
                    copy_bonus += 1000
                    wanted.append(ec)
                else:
                    copy_bonus -= 500
                    unwanted.append(ec)

        all_scores = {"extension_arm_len":[extension_len, extension_len_matrix[extension_len]],                      "ligation_arm_len":[ligation_len, ligation_len_matrix[ligation_len]],                      "extension_arm_gc":[extension_gc, arm_gc_con[extension_gc]],                      "ligation_arm_gc":[ligation_gc, arm_gc_con[ligation_gc]],                      "ligation_mask_penalty":ligation_mask_penalty,                      "extension_mask_penalty":extension_mask_penalty,                      "capture_gc_content":[capture_gc, cap_gc_con[capture_gc]],                      "copy_bonus": {"intended_copies": wanted, "unintended_copies": unwanted,                                     "desired_copies": desired_copies,"bonus": copy_bonus}, "gc_clamp": gc_clamp,                      "extension_next_bases":[extension_next, next_bases[extension_next.upper()]],                      "ligation_next_bases":[ligation_next, next_bases[ligation_next.upper()]],                         }
        # calculate total score
        score = (extension_len_matrix[extension_len] + ligation_len_matrix[ligation_len] +
                 arm_gc_con[extension_gc] + arm_gc_con[ligation_gc] +
                 next_bases[extension_next.upper()] + next_bases[ligation_next.upper()] +
                 gc_clamp + ligation_mask_penalty + extension_mask_penalty + cap_gc_con[capture_gc] + \
                 copy_bonus)
        # add mip_score to dictionary
        dic["pair_information"][mip]["mip_information"]["mip_score"] = score
        dic["pair_information"][mip]["mip_information"]["all_scores"] = all_scores
        dic["pair_information"][mip]["mip_information"]["capture_seq"] = capture_seq
    # write dictionary to json file
    outfile = open (primer3_output_DIR + output_file, "w")
    json.dump(dic, outfile, indent=1)
    outfile.close()
    return dic
def score_hs_mips(mip_file, primer3_output_DIR, output_file):
    """ Score mips in a dictionary according to scoring matrix
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then extension arm having
    GC clamp of 2. next_base values are last."""
    # open mip dictionary from file
    infile = open (primer3_output_DIR + mip_file, 'r')
    dic = json.load(infile)
    infile.close()
    # add "NEXT_BASES" tag:value pair to dictionary.
    # next_bases are 2 bases immediately downstream of extension primer and ligation primer
    # extract template sequence
    seq_template = dic["sequence_information"]["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for mip in dic["pair_information"]:
        # get coordinates of primers in the form of "start_base, len"
        extension_coord = dic["pair_information"][mip]["extension_primer_information"]["COORDINATES"]
        ligation_coord = dic["pair_information"][mip]["ligation_primer_information"]["COORDINATES"]
        # orientation of the mip is used to determine if extension arm or ligation arm is originating
        # from PRIMER_LEFT or PRIMER_RIGTH. When an arm originates from LEFT_PRIMER, it is on the plus
        # strand of DNA and its length is added to the start coordinate to get the end coordinate,
        # while it is subtracted for RIGHT_PRIMERs
        strand = dic["pair_information"][mip]["orientation"]
        if strand == "forward": # then the extension arm is the LEFT primer and ligation arm is RIGHT
            extension_end = int(extension_coord.split(",")[0]) + int(extension_coord.split(",")[1]) - 1
            # 1 is added or subtracted due to sequence index being zero based.
            ligation_end = int(ligation_coord.split(",")[0]) - int(ligation_coord.split(",")[1]) + 1
            ext_next = seq_template [(extension_end+1):(extension_end+3)]
            lig_next = reverse_complement(seq_template [(ligation_end - 2):ligation_end])
        elif strand == "reverse": # then the extension arm is the RIGHT primer and ligation is LEFT
            extension_end = int(extension_coord.split(",")[0]) - int(extension_coord.split(",")[1]) + 1
            ligation_end = int(ligation_coord.split(",")[0]) + int(ligation_coord.split(",")[1]) - 1
            ext_next = reverse_complement(seq_template [(extension_end - 2):extension_end])
            lig_next = seq_template [(ligation_end+1):(ligation_end+3)]
        # add "NEXT_BASES" key and its value to mip dictionary
        dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"] = ext_next
        dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"] = lig_next
    # arm gc content score matrix
    # define best to worst values for gc content. Best is very important so it has a huge effect on score.
    best = 1000
    mid = 100
    low = 10
    worst = 0
    # define matrix
    arm_gc_con = {}
    for i in range(20,81):
        if i < 35:
            arm_gc_con[i] = worst
        elif i < 40:
            arm_gc_con[i] = low
        elif i < 45:
            arm_gc_con[i] = mid
        elif i < 60:
                arm_gc_con[i] = best
        elif i < 65:
                arm_gc_con[i] = mid
        elif i < 70:
            arm_gc_con[i] = low
        else :
            arm_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases. This parameter should be important only when comparing
    # mips with equally good gc contents. Therefore the values are low and does not give a mip a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"GG":best, "CC":best, "GC":best, "CG":best,
                  "GA":mid, "CA":mid, "GT":mid, "CT":mid,
                  "AG":low, "TG":low, "AC":low, "TC":low,
                  "AA":worst, "AT":worst, "TA":worst, "TT":worst})
    # score all mips using above matrices
    for mip in list(dic["pair_information"].keys()):
        # get arm sequences
        ligation_seq = dic["pair_information"][mip]["ligation_primer_information"]["SEQUENCE"]
        extension_seq = dic["pair_information"][mip]["extension_primer_information"]["SEQUENCE"]
        # all arms have gc clamp 1. Check if extension arm has gc clamp of 2
        extension_clamp = 0
        if  (extension_seq[-2] == "G") or  (extension_seq[-2] == "C"):
            extension_clamp = 100 # score changed to 100 if arm ends in GG, GC, CG or CC
        # get gc percentages and convert to int.
        extension_gc = int(float(dic["pair_information"][mip]["extension_primer_information"]["GC_PERCENT"]))
        ligation_gc = int(float(dic["pair_information"][mip]["ligation_primer_information"]["GC_PERCENT"]))
        # get next base values
        extension_next = dic["pair_information"][mip]["extension_primer_information"]["NEXT_BASES"]
        ligation_next = dic["pair_information"][mip]["ligation_primer_information"]["NEXT_BASES"]
        # calculate total score
        score = (arm_gc_con[extension_gc] + arm_gc_con[ligation_gc] +
                 next_bases[extension_next.upper()] + next_bases[ligation_next.upper()] +
                 extension_clamp)
        # add mip_score to dictionary
        dic["pair_information"][mip]["mip_information"]["mip_score"] = score
    # write dictionary to json file
    outfile = open (primer3_output_DIR + output_file, "w")
    json.dump(dic, outfile, indent=1)
    outfile.close()
    return dic
# backbone compatible with truseq read 2 primer
hybrid_bb = "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNNNN"
# backbone compatible with truseq read 2 primer, molecular barcodes 10/4 split
hybrid_split = "NNNN" + "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNN"
# backbone compatible with truseq read 2 primer, molecular barcodes 10/4 split
hybrid_split_hp = "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNN"
# backbone compatible with truseq; plus GC buffer 3' downstream of primers
gc_bb = "GCAGATCGGAAGAGCACACCTCGCCAAGCTTTCGGC" + "NNNNNNNNNNNN"
# old backbone with molecular bar codes on one side
slx_bb = "CTTCAGCTTCCCGATCCGACGGTAGTGT" + "NNNNNNNNNNNN"
# backbone dictionary
backbones = {
            "hybrid_bb": "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNNNN",
            "hybrid_split": "NNNN" + "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNN",
            "hybrid_split_hp": "AGATCGGAAGAGCACACGTGACTCGCCAAGCTGAAG" + "NNNNNNNNNN",
            "gc_bb": "GCAGATCGGAAGAGCACACCTCGCCAAGCTTTCGGC" + "NNNNNNNNNNNN",
            "slx_bb": "CTTCAGCTTCCCGATCCGACGGTAGTGT" + "NNNNNNNNNNNN"
            }
def strip_fasta (sequence):
    seq_list = sequence.split('\n')[1:]
    seq_join = "".join(seq_list)
    return seq_join
def calculate_gc (sequence, fasta=0):
    if fasta:
        seq_list = sequence.split('\n')[1:]
        seq_join = "".join(seq_list)
        seq = seq_join.lower()

    else:
        seq = sequence.lower()
    gc_count = seq.count('g') + seq.count('c')
    at_count = seq.count('a') + seq.count('t')
    percent = int (gc_count * 100 / (gc_count + at_count))
    return percent
def translate(sequence, three_letter = False):
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
    gencode3 = {'A': 'Ala',
             'C': 'Cys',
             'D': 'Asp',
             'E': 'Glu',
             'F': 'Phe',
             'G': 'Gly',
             'H': 'His',
             'I': 'Ile',
             'K': 'Lys',
             'L': 'Leu',
             'M': 'Met',
             'N': 'Asn',
             'P': 'Pro',
             'Q': 'Gln',
             'R': 'Arg',
             'S': 'Ser',
             'T': 'Thr',
             'V': 'Val',
             'W': 'Trp',
             'Y': 'Tyr'}
    seq = sequence.upper()
    """Return the translated protein from 'sequence' assuming +1 reading frame"""
    if not three_letter:
        return ''.join([gencode.get(seq[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
    else:
        return ''.join([gencode3.get(gencode.get(seq[3*i:3*i+3],'X'), "X")                                     for i in range(len(sequence)//3)])
def aa_converter(aa_name):
    """
    Output 3 letter and 1 letter amino acid codes for a given
    list of 3 letter or 1 letter amino acid code list.
    """
    gencode3 = {'A': 'Ala',
             'C': 'Cys',
             'D': 'Asp',
             'E': 'Glu',
             'F': 'Phe',
             'G': 'Gly',
             'H': 'His',
             'I': 'Ile',
             'K': 'Lys',
             'L': 'Leu',
             'M': 'Met',
             'N': 'Asn',
             'P': 'Pro',
             'Q': 'Gln',
             'R': 'Arg',
             'S': 'Ser',
             'T': 'Thr',
             'V': 'Val',
             'W': 'Trp',
             'Y': 'Tyr'}
    for a in list(gencode3.keys()):
        gencode3[gencode3[a]] = a
    try:
        return gencode3[aa_name.capitalize()]
    except KeyError as e:
        return "%s is not a valid amino acid name" %a
def compatible_mip_check(m1, m2, overlap_same, overlap_opposite):
    d = m1.mip_dic
    es = ext_start = d["extension_primer_information"]["GENOMIC_START"]
    ee = ext_end = d["extension_primer_information"]["GENOMIC_END"]
    ls = lig_start = d["ligation_primer_information"]["GENOMIC_START"]
    le = lig_end = d["ligation_primer_information"]["GENOMIC_END"]
    # get mip orientation
    ori = d["orientation"]
    m1_set = set(list(range(min([es, ee]), max([es, ee]) + 1))
              + list(range(min([ls, le]), max([ls, le]) + 1)))
    m = m2.mip_dic
    nes = next_ext_start = m["extension_primer_information"]["GENOMIC_START"]
    nee = next_ext_end = m["extension_primer_information"]["GENOMIC_END"]
    nls = next_lig_start = m["ligation_primer_information"]["GENOMIC_START"]
    nle = next_lig_end = m["ligation_primer_information"]["GENOMIC_END"]
    # get mip orientation
    next_ori = m["orientation"]
    m2_set = set(list(range(min([nes, nee]), max([nes, nee]) + 1))
              + list(range(min([nls, nle]), max([nls, nle]) + 1)))
    overlap = len(m1_set.intersection(m2_set))
    if ori == next_ori:
        return overlap <= overlap_same
    else:
        return overlap <= overlap_opposite
def compatible_chains (primer_file, primer3_output_DIR, primer_out, output_file,
                      overlap_same = 0, overlap_opposite = 0, outp = 1):
    try:
        with open (primer3_output_DIR + primer_file, "r") as infile:
            scored_mips = json.load(infile)
    except IOError as e:
        print("Primer file does not exist.")
        return 1
    else:
        # create in/compatibility lists for each mip
        for k in list(scored_mips["pair_information"].keys()):
            # get coordinates of mip arms
            d = scored_mips["pair_information"][k]
            es = ext_start = d["extension_primer_information"]["GENOMIC_START"]
            ee = ext_end = d["extension_primer_information"]["GENOMIC_END"]
            ls = lig_start = d["ligation_primer_information"]["GENOMIC_START"]
            le = lig_end = d["ligation_primer_information"]["GENOMIC_END"]
            # get mip orientation
            ori = d["orientation"]
            #print k, es, ee, ls, le, ori
            # create an incompatibility list
            incompatible = []
            compatible = []
            for mip in list(scored_mips["pair_information"].keys()):
                m = scored_mips["pair_information"][mip]
                nes = next_ext_start = m["extension_primer_information"]                ["GENOMIC_START"]
                nee = next_ext_end = m["extension_primer_information"]                ["GENOMIC_END"]
                nls = next_lig_start = m["ligation_primer_information"]                ["GENOMIC_START"]
                nle = next_lig_end = m["ligation_primer_information"]                ["GENOMIC_END"]
                # get mip orientation
                next_ori = m["orientation"]
                compat = 0
                next_compat = 0
                # check if the two mips are compatible in terms of
                # orientation and coordinates
                if ori == next_ori == "forward":
                    if ((ls < nls) and (ls < nes + overlap_same)) or                      ((ls > nls) and (es + overlap_same> nls)):
                        compat = 1
                elif ori == next_ori == "reverse":
                    if ((ls < nls) and (es < nls + overlap_same)) or                       ((ls > nls) and (ls + overlap_same> nes)):
                        compat = 1
                elif (ori == "forward") and (next_ori == "reverse"):
                    if (ls < nls + overlap_opposite) or                     (es  + overlap_opposite> nes):
                        compat = 1
                    elif (es < nls) and (ee < nls + overlap_opposite) and                         (le  + overlap_opposite> nle) and                         (ls < nee + overlap_opposite):
                        compat = 1
                        next_compat = 1
                    elif (es > nls) and (es  + overlap_opposite> nle) and                          (ee < nee + overlap_opposite) and                         (le + overlap_opposite > nes):
                        compat = 1
                elif (ori == "reverse") and (next_ori == "forward"):
                    if (ls + overlap_opposite > nls) or                     (es < nes + overlap_opposite):
                        compat = 1
                    elif (ls > nes) and (ls + overlap_opposite > nee) and                         (le < nle + overlap_opposite) and                         (ee + overlap_opposite>nls):
                        compat = 1
                    elif (ls < nes) and (le < nes + overlap_opposite) and                          (ee + overlap_opposite > nee) and                         (es < nle + overlap_opposite):
                        compat = 1
                        next_compat = 1
                if not compat:
                    incompatible.append(mip)
                if next_compat:
                    compatible.append(mip)
            d["incompatible"] = incompatible
            d["compatible"] = compatible
        #f = open(primer3_output_DIR + output_file, "w")
        f = []
        def compatible_recurse (l):
            """
            Take a list, l,  of numbers that represent a mip set with
            their corresponding "place" in the mip dictionary, and index
            number, i. Find the subset of mips in the rest of the list
            that are compatible with the mip at index i, using compatibility
            dictionary d. For each mip in the subset, find compatible mips
            in the rest of the list. Recurse until the subset does not have
            any mips. Append each compatible subset to a final result list, f.
            """
            incomp = list(l)
            for il in l:
                incomp.extend(scored_mips["pair_information"][il]["incompatible"])
            comp = set(
                scored_mips["pair_information"][l[-1]]["compatible"]
            ).difference(incomp)
            if len(comp) > 0:
                for n in comp:
                    compatible_recurse(l + [n])
            else:
                #f.write(",".join(l) + "\n")
                """
                l_co = []
                for member in l:
                    d = scored_mips["pair_information"][member]
                    es = d["extension_primer_information"]["GENOMIC_START"]
                    ee = d["extension_primer_information"]["GENOMIC_END"]
                    ls = d["ligation_primer_information"]["GENOMIC_START"]
                    le = d["ligation_primer_information"]["GENOMIC_END"]
                    l_co.extend([es, ee, ls, le])
                l.append([[min(l_co),
                          max(l_co)]])

                """
                f.append(l)
        keys = list(scored_mips["pair_information"].keys())
        for k in keys:
            comp_list = scored_mips["pair_information"][k]["compatible"]
            if len(comp_list) > 0:
                # for each of the mips in the compatibility list,
                for m in comp_list:
                    # create an initial result list to be used by the
                    # compatible_recurse function
                    compatible_recurse([k, m])
            else:
                #comp_list = [k]
                #f.write(k + "\n")
                """
                d = scored_mips["pair_information"][k]
                es = d["extension_primer_information"]["GENOMIC_START"]
                ee = d["extension_primer_information"]["GENOMIC_END"]
                ls = d["ligation_primer_information"]["GENOMIC_START"]
                le = d["ligation_primer_information"]["GENOMIC_END"]
                f.append([k, [[min([es, ee, ls, le]),
                          max([es, ee, ls, le])]]])
                """
                f.append([k])
        fs = [frozenset(fset) for fset in f]
        set_count = len(fs)
        counter = 0
        while((set_count < 10000) and (counter <= 20)):
            counter += 1
            new_fs = set()
            for s1 in fs:
                inc = set()
                for m in s1:
                    inc.update(scored_mips["pair_information"][m]["incompatible"])
                for s2 in fs:
                    if s1 != s2:
                        s3 = s2.difference(inc)
                        if len(s3) > 0:
                            new_fs.add(s1.union(s3))
                        else:
                            new_fs.add(s1)
            fs = new_fs
            if len(fs) < 30000:
                new_fs = set()
                for s1 in fs:
                    for s2 in fs:
                        if s1 != s2:
                            if s1.issubset(s2):
                                break
                    else:
                        new_fs.add(s1)
                fs = new_fs
            set_count = len(fs)
            if 30000 > set_count > 10000:
                fs = [set(s) for s in fs]
                new_fs = set()
                for s1 in fs:
                    inc = set()
                    for m in s1:
                        inc.update(scored_mips["pair_information"][m]["incompatible"])
                    for s2 in fs:
                        if s1 != s2:
                            s1.update(s2.difference(inc))
                    new_fs.add(frozenset(s1))
                fs = new_fs
        if outp:
            with open(primer3_output_DIR + output_file, "w") as outfile:
                outfile.write("\n".join([",".join(s) for s in fs]) + "\n")
        #f.close()
        #print len(output)
        with open(primer3_output_DIR + primer_out, "w") as outfile:
            json.dump(scored_mips, outfile, indent=1)
    return fs
def compatible_mips(primer_file, primer3_output_DIR, primer_out, output_file,
               overlap_same = 0, overlap_opposite = 0):
    try:
        with open (primer3_output_DIR + primer_file, "r") as infile:
            scored_mips = json.load(infile)
    except IOError as e:
        print("Primer file does not exist.")
        return 1
    else:
        # create in/compatibility lists for each mip
        for k in list(scored_mips["pair_information"].keys()):
            # get coordinates of mip arms
            d = scored_mips["pair_information"][k]
            es = ext_start = d["extension_primer_information"]["GENOMIC_START"]
            ee = ext_end = d["extension_primer_information"]["GENOMIC_END"]
            ls = lig_start = d["ligation_primer_information"]["GENOMIC_START"]
            le = lig_end = d["ligation_primer_information"]["GENOMIC_END"]
            # get mip orientation
            ori = d["orientation"]
            #print k, es, ee, ls, le, ori
            # create an incompatibility list
            incompatible = []
            compatible = []
            for mip in list(scored_mips["pair_information"].keys()):
                m = scored_mips["pair_information"][mip]
                nes = next_ext_start = m["extension_primer_information"]                ["GENOMIC_START"]
                nee = next_ext_end = m["extension_primer_information"]                ["GENOMIC_END"]
                nls = next_lig_start = m["ligation_primer_information"]                ["GENOMIC_START"]
                nle = next_lig_end = m["ligation_primer_information"]                ["GENOMIC_END"]
                # get mip orientation
                next_ori = m["orientation"]
                compat = 0
                # check if the two mips are compatible in terms of
                # orientation and coordinates
                if ori == next_ori == "forward":
                    if ((ls < nls) and (ls < nes + overlap_same)) or                      ((ls > nls) and (es + overlap_same> nls)):
                        compat = 1
                elif ori == next_ori == "reverse":
                    if ((ls < nls) and (es < nls + overlap_same)) or                       ((ls > nls) and (ls + overlap_same> nes)):
                        compat = 1
                elif (ori == "forward") and (next_ori == "reverse"):
                    if (ls < nls + overlap_opposite) or                     (es  + overlap_opposite> nes):
                        compat = 1
                    elif (es < nls) and (ee < nls + overlap_opposite) and                         (le  + overlap_opposite> nle) and                         (ls < nee + overlap_opposite):
                        compat = 1
                    elif (es > nls) and (es  + overlap_opposite> nle) and                          (ee < nee + overlap_opposite) and                         (le + overlap_opposite > nes):
                        compat = 1
                elif (ori == "reverse") and (next_ori == "forward"):
                    if (ls + overlap_opposite > nls) or                     (es < nes + overlap_opposite):
                        compat = 1
                    elif (ls > nes) and (ls + overlap_opposite > nee) and                         (le < nle + overlap_opposite) and                         (ee + overlap_opposite>nls):
                        compat = 1
                    elif (ls < nes) and (le < nes + overlap_opposite) and                          (ee + overlap_opposite > nee) and                         (es < nle + overlap_opposite):
                        compat = 1
                if not compat:
                    incompatible.append(mip)
                else:
                    compatible.append(mip)
            d["incompatible"] = incompatible
            d["compatible"] = compatible
        f = open(primer3_output_DIR + output_file, "w")
        #f = []
        def compatible_recurse (l):
            """
            Take a list, l,  of numbers that represent a mip set with
            their corresponding "place" in the mip dictionary, and index
            number, i. Find the subset of mips in the rest of the list
            that are compatible with the mip at index i, using compatibility
            dictionary d. For each mip in the subset, find compatible mips
            in the rest of the list. Recurse until the subset does not have
            any mips. Append each compatible subset to a final result list, f.
            """
            incomp = list(l)
            for il in l:
                incomp.extend(scored_mips["pair_information"][il]["incompatible"])
            comp = set(scored_mips["pair_information"][l[-1]]["compatible"]).difference(incomp)
            if len(comp) > 0:
                for n in comp:
                    compatible_recurse(l + [n])
            else:
                f.write(",".join(l) + "\n")
                #f.append(l)
        keys = list(scored_mips["pair_information"].keys())
        for k in keys:
            comp_list = scored_mips["pair_information"][k]["compatible"]
            if len(comp_list) > 0:
                # for each of the mips in the compatibility list,
                for m in comp_list:
                    # create an initial result list to be used by the
                    # compatible_recurse function
                    compatible_recurse([k, m])
            else:
                comp_list = [k]
                f.write(k + "\n")

        f.close()
        #print len(output)
        with open(primer3_output_DIR + primer_out, "w") as outfile:
            json.dump(scored_mips, outfile, indent=1)
    return
def compatibility (scored_mips, primer3_output_DIR = "", primer_out = "",
                      overlap_same = 0, overlap_opposite = 0):
    # create in/compatibility lists for each mip
    for k in list(scored_mips["pair_information"].keys()):
        # get coordinates of mip arms
        d = scored_mips["pair_information"][k]
        es = ext_start = d["extension_primer_information"]["GENOMIC_START"]
        ee = ext_end = d["extension_primer_information"]["GENOMIC_END"]
        ls = lig_start = d["ligation_primer_information"]["GENOMIC_START"]
        le = lig_end = d["ligation_primer_information"]["GENOMIC_END"]
        # get mip orientation
        ori = d["orientation"]
        #print k, es, ee, ls, le, ori
        # create an incompatibility list
        incompatible = []
        compatible = []
        for mip in list(scored_mips["pair_information"].keys()):
            m = scored_mips["pair_information"][mip]
            nes = next_ext_start = m["extension_primer_information"]            ["GENOMIC_START"]
            nee = next_ext_end = m["extension_primer_information"]            ["GENOMIC_END"]
            nls = next_lig_start = m["ligation_primer_information"]            ["GENOMIC_START"]
            nle = next_lig_end = m["ligation_primer_information"]            ["GENOMIC_END"]
            # get mip orientation
            next_ori = m["orientation"]
            compat = 0
            next_compat = 0
            # check if the two mips are compatible in terms of
            # orientation and coordinates
            if ori == next_ori == "forward":
                if ((ls < nls) and (ls < nes + overlap_same)) or                  ((ls > nls) and (es + overlap_same> nls)):
                    compat = 1
            elif ori == next_ori == "reverse":
                if ((ls < nls) and (es < nls + overlap_same)) or                   ((ls > nls) and (ls + overlap_same> nes)):
                    compat = 1
            elif (ori == "forward") and (next_ori == "reverse"):
                if (ls < nls + overlap_opposite) or                 (es  + overlap_opposite> nes):
                    compat = 1
                elif (es < nls) and (ee < nls + overlap_opposite) and                     (le  + overlap_opposite> nle) and                     (ls < nee + overlap_opposite):
                    compat = 1
                    next_compat = 1
                elif (es > nls) and (es  + overlap_opposite> nle) and                      (ee < nee + overlap_opposite) and                     (le + overlap_opposite > nes):
                    compat = 1
            elif (ori == "reverse") and (next_ori == "forward"):
                if (ls + overlap_opposite > nls) or                 (es < nes + overlap_opposite):
                    compat = 1
                elif (ls > nes) and (ls + overlap_opposite > nee) and                     (le < nle + overlap_opposite) and                     (ee + overlap_opposite>nls):
                    compat = 1
                elif (ls < nes) and (le < nes + overlap_opposite) and                      (ee + overlap_opposite > nee) and                     (es < nle + overlap_opposite):
                    compat = 1
                    next_compat = 1
            if not compat:
                incompatible.append(mip)
            if next_compat:
                compatible.append(mip)
        d["incompatible"] = incompatible
        d["compatible"] = compatible
    """
    with open(primer3_output_DIR + primer_out, "w") as outfile:
        json.dump(scored_mips, outfile, indent=1)
    """

    return scored_mips
def best_mip_set (compatible_mip_sets, compatible_mip_dic, num_para, diff_score_dic, outfile):
    # open file containing compatible mip lists
    with open(primer3_output_DIR + compatible_mip_sets, "r") as infile:
        mip_sets = json.load(infile)
    # load dict file that has the compatibility information and mip information
    with open(primer3_output_DIR + compatible_mip_dic, "r") as infile:
        mip_dic = json.load(infile)
        # if there is any sets in the list
        if len(mip_sets) > 0:
            best_set = []
            best_score = 0
            # create a number to mip name look up dictionary
            # because mip list has only numbers that correspond to a mips place in dict
            num_lookup = {}
            for mip in list(mip_dic["pair_information"].keys()):
                d = mip_dic["pair_information"][mip]
                place = d["place"]
                num_lookup[place] = mip
            # find which diffs each mip set captures
            for mip_set in mip_sets:
                if mip_set != None:
                    # create a dic for diffs captured cumulatively by all mips in the set
                    merged_caps = {}
                    # create a list for mip scores based on mip sequence and not the captured diffs
                    mip_scores = []
                    # create a variable for mip scores based on the diffs captured
                    diff_scores = 0
                    for m in mip_set:
                        # extract the mip name
                        mip_key = num_lookup[m]
                        # extract the captured diffs from the mip_dic and append to capture list
                        caps = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]
                        # get the diff name (e.g. chr1:1000-1001), use it as key for merged caps
                        # using a dict ensures nonredundancy
                        for diff in caps:
                            merged_caps[diff] = caps[diff]
                        # find out how many paralogs the mip identifies
                        mip_para = ["pair_information"][mip_key]["mip_information"]                                   ["captured_paralog_number"]
                        # extract mip score and append to mip_scores
                        ms = mip_dic["pair_information"][mip_key]                                    ["mip_information"]["mip_score"]
                        # add a bonus for mips identifying multiple paralogs
                        cap_bonus = mip_para * cap_coefficient
                        ms += cap_bonus
                        # add total mip score to mip scores
                        mip_scores.append(ms)
                    # to get more evenly (and well) scoring mips rather than having
                    # the top scoring set, create another bonus score
                    bonus_score = 0
                    # get total number of mips in the set
                    mip_count = len(mip_scores)
                    # define a counter for mips scoring above a certain score
                    ms_count = 0
                    for ms in mip_scores:
                        if ms >= lowest_mip_score:
                            ms_count += 1
                    # ideally we want mip_count == ms_count for all mips to be good
                    # but this does not work well for large sets because there is a good
                    # chance one or two will not score well and bonus will always be zero.
                    if (mip_count - ms_count) <= max_poor_mips:
                        bonus_score = uniformity_bonus
                    # create a dict that has how many paralogs are captured by the mip set
                    paralog_caps = {}
                    # populate dict with paralog numbers and an initial value of 0
                    for i in range(num_para):
                        paralog_caps[i] = 0
                    for d in merged_caps:
                        # find out which paralogs the diff uniquely identifies
                        # this information is only relevant for diffs from pdiffs file
                        # and exonic diffs and filtered diffs are the only ones from pdiffs in targets
                        if (merged_caps[d]["source"] == "exonic_diffs") or                            (merged_caps[d]["source"] == "filtered_diffs"):
                            # the diff is in the form a:1067:CTA:0,1,2 and we want to
                            # analyze the last part of it for paralogs identified
                            difference = (merged_caps[d]["diff"].split(":")[-1]).split(",")
                            # add the paralogs identified to mip_set paralog caps
                            for j in difference:
                                paralog_caps[int(j)] += 1
                        # extract the score of the type of difference and increment diff_scores
                        diff_scores += diff_score_dic[merged_caps[d]["source"]]
                    # calculate how many paralogs identified by the set
                    cap_para = 0
                    for k in paralog_caps:
                        if paralog_caps[k] > 0:
                            cap_para += 1
                    # calculate total score of mip set
                    total_score = bonus_score + sum(mip_scores) + (cap_para**2) * paralog_score +                                   diff_scores
                    if total_score > best_score:
                        best_score = total_score
                        best_set = mip_set
            # print the gene name and best set of mips together with total score
            #print gene[0]
            #print best_set, best_score
            temp_dic = {}
            # print scores and diffs captured by each mip in the best set
            for mip in best_set:
                mip_key = num_lookup[mip]
                temp_dic[mip_key] = mip_dic["pair_information"][mip_key]
                #print mip_dic["pair_information"][mip_key]["mip_information"]["mip_score"]
                for diff in mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]:
                    di = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"][diff]
                    #print diff, di["source"], di["diff"]

            with open(primer3_output_DIR + best_mip_sets, "w")  as outfile:
                json.dump(temp_dic, outfile, indent=1)
        # if there are no compatible mipsets, then find the best mip
        elif (len(mip_sets) == 0) and (len(list(mip_dic["pair_information"].keys())) > 0):
            best_score = 0
            best_mip = ""
            for mip_key in list(mip_dic["pair_information"].keys()):
                # extract the captured diffs from the mip_dic and append to capture list
                caps = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]
                # score diffs captured
                diff_scores = 0
                # find diff scores from diff source using diff score dict
                for diff in caps:
                    # extract the score of the type of difference and increment diff_scores
                    diff_scores += diff_score_dic[caps[diff]["source"]]
                # extract mip score
                mip_scores = mip_dic["pair_information"][mip_key]["mip_information"]["mip_score"]
                # find how many paralogs uniquely identified by the mip
                cap_para = mip_dic["pair_information"][mip_key]["mip_information"]                                   ["captured_paralog_number"]
                cap_bonus = cap_para * cap_coefficient
                # calculate total score of mip set
                total_score = mip_scores + (cap_para**2) * paralog_score + diff_scores + cap_bonus
                if total_score > best_score:
                    best_score = total_score
                    best_mip = mip_key
            temp_dic = {}
            temp_dic[best_mip] = mip_dic["pair_information"][best_mip]
            #print gene[0]
            #print [mip_dic["pair_information"][best_mip]["place"]]
            #print best_score
            #print mip_dic["pair_information"][best_mip]["mip_information"]["mip_score"]
            dics = mip_dic["pair_information"][best_mip]["mip_information"]["captured_diffs"]
            #for diff in dics:
                #print dics[diff]["source"], dics[diff]["diff"]
            with open(primer3_output_DIR + scored + "_best_set", "w")  as outfile:
                json.dump(temp_dic, outfile, indent=1)
        else:
            print(("No mips available for target region ", gene[0]))
        return
def get_analysis_settings(settings_file):
    """ Convert analysis settings file to dictionary"""
    settings = {}
    with open(settings_file) as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                value = newline[1].split(",")
                if len(value) == 1:
                    settings[newline[0]] = value[0]
                else:
                    settings[newline[0]] = [v for v in value if v != ""]
    return settings
def filter_mipster (settings):
    """
    Import data from Mipster pipeline, filter and save to filtered_file
    We use a tab separated file coming from Mipster pipeline for data analysis
    This file has a lot of information that will not be used initially
    We will select the columns that will be used.
    Selected colnames are should be specified in settings file
    """
    wdir = settings["workingDir"]
    mipster_file = wdir + settings["mipsterFile"]
    filtered_file = wdir + settings["filteredMipsterFile"]
    # extract data column names from the mipster file
    with open(mipster_file) as infile:
        all_colnames = infile.readline().strip().split("\t")
    # select the relevant columns
    colnames = settings["colNames"]
    # colnames are names from the mipster pipeline and
    # given names are names to be used in this pipeline
    given_names = settings["givenNames"]
    # check if column names are in equal number and match
    col_map = {}
    if len(colnames) == len(given_names):
        for i in range(len(colnames)):
            col_map [colnames[i]] = {"given_name": given_names[i]}
    else:
        print("number of column names and number of given names do not match! Colnames will be used")
        for i in range(len(colnames)):
            col_map [colnames[i]] = {"given_name": colnames[i]}
    # get index of colnames in data file
    for i in range(len(all_colnames)):
        if all_colnames[i] in colnames:
            col_map[all_colnames[i]]["index"] = i
    # extract data of selected columns
    with open(filtered_file, "w") as outfile, open(mipster_file) as all_data:
        outfile.write("#" + "\t".join(given_names) + "\n")
        counter = 0
        for line in all_data:
            if counter > 0 and not line.startswith("s_Sample"):
                newline = line.strip().split("\t")
                outlist = []
                for c in colnames:
                    outlist.append(newline[col_map[c]["index"]])
                outfile.write("\t".join(outlist) + "\n")
            counter += 1
    return
def get_haplotypes(settings):
    """ 1) Extract all haplotypes from new data.
        2) Remove known haplotypes using previous data (if any).
        3) Map haplotypes to species genome to get the best hit(s)
        4) Crosscheck best bowtie hit with the targeted region
        5) Output haplotypes dictionary and off_targets dictionary

        Once this function is called, we will get the new haplotypes present
        in this data set that are on target and where they map on the genome.
        Mapping haplotypes to specific targets/copies is not accomplished here
        """
    wdir = settings["workingDir"]
    mipster_file = wdir + settings["mipsterFile"]
    haplotypes_fq_file = wdir + settings["haplotypesFastqFile"]
    haplotypes_sam_file = wdir + settings["haplotypesSamFile"]
    bwa_options = settings["bwaOptions"]
    sequence_to_haplotype_file = wdir + settings["sequenceToHaplotypeDictionary"]
    call_info_file = settings["callInfoDictionary"]
    species = settings["species"]
    try:
        tol = int(settings["alignmentTolerance"])
    except KeyError as e:
        tol = 50

    #wdir = settings["workingDir"]
    ### DATA EXTRACTION ###
    # if there is no previous haplotype information, an empty dict will be used
    # for instead of the known haplotypes dict
    try:
        with open(sequence_to_haplotype_file) as infile:
            sequence_to_haplotype = json.load(infile)
    except IOError as e:
        sequence_to_haplotype = {}
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    # extract haplotype name and sequence from new data file
    haps = {}
    # extract data column names from the mipster file
    with open(mipster_file) as infile:
        line_number = 0
        for line in infile:
            newline = line.strip().split("\t")
            line_number += 1
            if line_number == 1:
                for i in range(len(newline)):
                    if newline[i] in ["haplotype_ID",
                                      "h_popUID"]:
                        hap_index = i
                    elif newline[i] in ["haplotype_sequence",
                                        'h_seq']:
                        seq_index = i
            else:
                hapname = newline[hap_index]
                hapseq = newline[seq_index]
                hapqual = "H" * len(hapseq)
                # add the haplotype to the dict if it has not been mapped before
                if hapseq not in sequence_to_haplotype:
                    haps[hapname] = {"sequence": hapseq,
                                     "quality": hapqual}

    ### BWA alignment ####
    # create a fastq file for bwa input
    with open(haplotypes_fq_file, "w") as outfile:
        for h in haps:
            outfile.write("@" + h + "\n")
            outfile.write(haps[h]["sequence"] + "\n" + "+" + "\n")
            outfile.write(haps[h]["quality"] + "\n")
    # re-structure haplotypes dictionary and initialize a hits dictionary for all
    # haplotypes that will hold the bowtie hits for each of the haplotypes
    # keys for this dict will be mipnames
    haplotypes = {}
    for h in haps:
        mip_name = h.split(".")[0]
        try:
            haplotypes[mip_name][h] = {"sequence": haps[h]["sequence"]}
        except KeyError as e:
            haplotypes[mip_name] = {h: {"sequence": haps[h]["sequence"]}}
    # run bwa
    bwa(haplotypes_fq_file, haplotypes_sam_file, "sam", "", "",
            bwa_options, species)
    # get best hits from alignment results
    hap_hits = {}
    for h in haps:
        # initialize each haplotype with an empty list and a -5000 score
        hap_hits[h] = [[], [-5000]]
    # find the best bwa hit(s) for each genotype
    with open(haplotypes_sam_file) as infile:
        for line in infile:
            if not line.startswith("@"):
                newline = line.strip().split("\t")
                try:
                    if newline[13].startswith("AS"):
                        score = int(newline[13].split(":")[-1])
                    else:
                        score = -5000
                except IndexError as e:
                    if newline[11].startswith("AS"):
                        score = int(newline[11].split(":")[-1])
                    else:
                        score = -5000
                hapname = newline[0]
                if max(hap_hits[hapname][1]) < score:
                    hap_hits[hapname][0] = [newline]
                    hap_hits[hapname][1] = [score]
                elif max(hap_hits[hapname][1]) == score:
                    hap_hits[hapname][0].append(newline)
                    hap_hits[hapname][1].append(score)
    # update haplotypes dict with bowtie best hit information
    for m in haplotypes:
        for h in haplotypes[m]:
            haplotypes[m][h]["best_hits"] = hap_hits[h]
    # crosscheck the best bwa hit(s) for each haplotype with mip targets
    # mark off target haplotypes
    for m in haplotypes:
        gene_name = m.split("_")[0]
        try:
            call_dict = call_info[gene_name][m]["copies"]
            for h in list(haplotypes[m].keys()):
                haplotypes[m][h]["mapped"] = False
                best_hits = haplotypes[m][h]["best_hits"][0]
                for hit in best_hits:
                    if haplotypes[m][h]["mapped"]:
                        break
                    hit_chrom = hit[2]
                    hit_pos = int(hit[3])
                    for copy_name in call_dict:
                        copy_chrom = call_dict[copy_name]["chrom"]
                        copy_begin = call_dict[copy_name]["capture_start"]
                        copy_end = call_dict[copy_name]["capture_end"]
                        if (
                            (copy_chrom == hit_chrom)
                            and (copy_begin -tol < hit_pos < copy_end + tol)
                        ):
                            haplotypes[m][h]["mapped"] = True
                            break
        except KeyError as e:
            for h in list(haplotypes[m].keys()):
                haplotypes[m][h]["mapped"] = False
    # remove haplotypes that mapped best to an untargeted location on genome
    off_target_haplotypes = {}
    for m in list(haplotypes.keys()):
        for h in list(haplotypes[m].keys()):
            if not haplotypes[m][h]["mapped"]:
                """
                anoter solution to this should be found
                if m.startswith("AMELX"):
                    # AMELX also maps to Y chromosome, which is not off target
                    # this is a quick fix for now but ideally Y chromosome
                    # should be included in the design as a paralogous copy

                    haplotypes[m][h]["mapped"] = True
                else:
                """
                off_target_haplotypes[h] = haplotypes[m].pop(h)
        if len(haplotypes[m]) == 0:
            haplotypes.pop(m)
    hap_file = wdir + settings["tempHaplotypesFile"]
    off_file = wdir + settings["tempOffTargetsFile"]
    with open(hap_file, "w") as out1, open(off_file, "w") as out2:
        json.dump(haplotypes, out1, indent=1)
        json.dump(off_target_haplotypes, out2, indent=1)
    return
def rename_mipster_haplotypes(settings):
    """ 1) Extract all haplotypes from new data.
        2) Remove known haplotypes using previous data (if any).
        3) Map haplotypes to species genome to get the best hit(s)
        4) Crosscheck best bowtie hit with the targeted region
        5) Output haplotypes dictionary and off_targets dictionary

        Once this function is called, we will get the new haplotypes present
        in this data set that are on target and where they map on the genome.
        Mapping haplotypes to specific targets/copies is not accomplished here
        """
    wdir = settings["workingDir"]
    filtered_data = wdir + settings["filteredMipsterFile"]
    renamed_file = filtered_data + "_renamed"
    sequence_to_haplotype_file = wdir +     settings["sequenceToHaplotypeDictionary"]
    # if there is no previous haplotype information, an empty dict will be used
    # for instead of the known haplotypes dict
    with open(sequence_to_haplotype_file) as infile:
            sequence_to_haplotype = json.load(infile)
    # extract haplotype name and sequence from new data file
    problem_sequences = []
    with open(filtered_data) as infile, open(renamed_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                filtered_colnames = line.strip().split("\t")
                filtered_colnames[0] = filtered_colnames[0][1:]
                # find column indexes of relevant fields
                for i in range(len(filtered_colnames)):
                    if filtered_colnames[i] == "haplotype_ID":
                        hap_index = i
                    elif filtered_colnames[i] == "haplotype_sequence":
                        seq_index = i
                    elif filtered_colnames[i] == "haplotype_quality_scores":
                        qual_index = i
                outfile.write(line)
            else:
                newline = line.split("\t")
                hapseq = newline[seq_index]
                try:
                    hapname = sequence_to_haplotype[hapseq]
                    newline[hap_index] = hapname
                    outfile.write("\t".join(newline))
                except KeyError as e:
                    problem_sequences.append(hapseq)
    return problem_sequences
def get_summary_table(settings):
    wdir = settings["workingDir"]
    call_info_file = settings["callInfoDictionary"]
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = []
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.extend(json.load(infile)[psk][1:])
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    with open(unique_haplotype_file) as infile:
        haplotypes = json.load(infile)
    # prepare a list of all mips used in the experiment
    all_mips = []
    for g in call_info:
        for m in call_info[g]:
            if m in used_probes:
                c = "C0"
                temp = [g, m, call_info[g][m]["copies"][c]["chrom"],
                        call_info[g][m]["copies"][c]["capture_start"],
                        call_info[g][m]["copies"][c]["capture_end"]]
                all_mips.append(temp)
    all_mips_srt = sorted(all_mips, key=itemgetter(2,3,4))
    all_mips_keys = {}
    for i in range(len(all_mips_srt)):
        m = all_mips_srt[i]
        k = m[0] + ":" + m[1]
        all_mips_keys[k] = i
    with open(wdir + settings["perSampleResults"]) as infile:
        sample_res = json.load(infile)
    sample_names = sorted(sample_res)
    sample_mip_counts = {}
    for s in sample_names:
        barcode_counts = list(map(int, np.zeros(len(all_mips_srt))))
        for g in sample_res[s]:
            for m in sample_res[s][g]:
                res = sample_res[s][g][m]
                k = ":".join([g, m])
                location = all_mips_keys[k]
                for c in res:
                    filtered_data = res[c]["filtered_data"]
                    for f in filtered_data:
                        hid = f["haplotype_ID"]
                        bc = f["barcode_count"]
                        hap = haplotypes[m][hid]
                        copies = hap["mapped_copies"]
                        for cop in copies:
                            barcode_counts[location] += bc
        sample_mip_counts[s] = barcode_counts
    all_probes = []
    for g in call_info:
        for m in call_info[g]:
            if m in used_probes:
                for c in call_info[g][m]["copies"]:
                    temp = [g, m, c, call_info[g][m]["copies"][c]["chrom"],
                            call_info[g][m]["copies"][c]["capture_start"],
                            call_info[g][m]["copies"][c]["capture_end"]]
                    all_probes.append(temp)
    all_probes_srt = sorted(all_probes, key = itemgetter(3,4,5))
    all_probes_keys = {}
    for i in range(len(all_probes_srt)):
        p = all_probes_srt[i]
        g = p[0]
        m = p[1]
        c = p[2]
        k = ":".join([g, m, c])
        all_probes_keys[k] = i

    bc_header = ["MIP"] + [a[1] for a in all_mips_srt]
    gene_names = ["target_gene_group"] + [a[0] for a in all_mips_srt]
    bc_info = [[s] + sample_mip_counts[s] for s in sorted(sample_names)]
    bc_out = np.transpose([bc_header] + [gene_names] + bc_info)

    """
    with open(wdir + settings["mipCountFile"], "w") as outfile:
        outfile.write("\n".join(["\t".join(b) for b in bc_out]))
    """
    write_list(bc_out, wdir + settings["mipCountFile"])




    min_count = int(settings["minSnpBarcodeCount"])
    min_snp_qual = int(settings["minSnpQuality"])
    min_snp_frac = float(settings["minSnpBarcodeFraction"])
    sample_barcode_counts_uniq = {}
    sample_barcode_counts_multi = {}
    for s in sample_names:
        uniq_barcode_counts = list(map(int, np.zeros(len(all_probes_srt))))
        multi_barcode_counts = []
        for g in sample_res[s]:
            for m in sample_res[s][g]:
                res = sample_res[s][g][m]
                for c in res:
                    filtered_data = res[c]["filtered_data"]
                    copy_barcode_total = 0
                    for f in filtered_data:
                        copy_barcode_total += f["barcode_count"]
                    copy_barcode_total = float(copy_barcode_total)
                    for f in filtered_data:
                        hid = f["haplotype_ID"]
                        bc = f["barcode_count"]
                        if bc >= min_count and bc/copy_barcode_total >= min_snp_frac:
                            hap = haplotypes[m][hid]
                            copies = hap["mapped_copies"]
                            if len(copies) == 1:
                                cop = list(copies.keys())[0]
                                k = ":".join([g, m, cop])
                                location = all_probes_keys[k]
                                uniq_barcode_counts[location] += bc
                            else:
                                location = []
                                for cop in copies:
                                    k = ":".join([g, m, cop])
                                    location.append(all_probes_keys[k])
                                multi_barcode_counts.append([location, bc, hid])
        sample_barcode_counts_uniq[s] = uniq_barcode_counts
        sample_barcode_counts_multi[s] = multi_barcode_counts
    temporary_barcode_counts = np.array([sample_barcode_counts_uniq[s]                                         for s in sorted(sample_names)])
    sample_medians = np.median(temporary_barcode_counts, axis = 1)[:, np.newaxis]
    temporary_sample_norm = temporary_barcode_counts/sample_medians
    probe_medians = np.nanmedian(np.where(temporary_sample_norm != 0,
                                          temporary_sample_norm,
                                          np.nan), axis = 0)[np.newaxis, :]
    temporary_copy_counts = temporary_sample_norm/probe_medians
    nan_mask = np.isnan(temporary_copy_counts)
    inf_mask = np.isinf(temporary_copy_counts)
    temporary_copy_counts[nan_mask] = 0
    temporary_copy_counts[inf_mask] = 0
    sample_barcode_counts = copy.deepcopy(sample_barcode_counts_uniq)
    multi_mapper_ratios = {}
    all_probes_srt_array = np.array(all_probes_srt)
    for sample_index in range(len(sample_names)):
        s = sorted(sample_names)[sample_index]
        multi_barcodes = sample_barcode_counts_multi[s]
        uniq_barcodes = sample_barcode_counts[s]
        multi_mapper_ratios[s] = {}
        for mb in multi_barcodes:
            locs = mb[0]
            bc = mb[1]
            multi_mapper_ratios[s][mb[2]] = {}
            copy_medians_for_loc = []
            for l in locs:
                gene_at_l = all_probes_srt[l][0]
                copy_at_l = all_probes_srt[l][2]
                gene_mask = all_probes_srt_array[:, 0] == gene_at_l
                copy_mask = all_probes_srt_array[:, 2] == copy_at_l
                gene_copy_mask = np.logical_and(gene_mask, copy_mask)
                unique_counts = temporary_copy_counts[sample_index, gene_copy_mask]
                uniq_meds = np.median(unique_counts)
                copy_medians_for_loc.append(uniq_meds)
            copy_ratios = copy_medians_for_loc
            #copy_ratios = np.array(copy_ratios)
            #nan_mask = np.isnan(copy_ratios)
            #copy_ratios[nan_mask] = 0
            copy_totals = float(np.sum(copy_ratios))
            if copy_totals == 0:
                norm_copy_ratios = [1./len(copy_ratios) for cpr in copy_ratios]
            else:
                norm_copy_ratios = [cpr/copy_totals for cpr in copy_ratios]
            for location_index in range(len(locs)):
                l = locs[location_index]
                l_ratio = norm_copy_ratios[location_index]
                uniq_barcodes[l] += bc * l_ratio
                multi_mapper_ratios[s][mb[2]][l] = l_ratio




    mip_header = ["MIP"] + [a[1] for a in all_probes_srt]
    copy_header = ["copy"] + [a[2] + "_" + a[2] for a in all_probes_srt]
    uniq_mip_header = ["MIP_copy"] + [a[1] + "_" + a[2] for a in all_probes_srt]
    gene_names = ["target_gene_group"] + [a[0] for a in all_probes_srt]
    chromosomes = ["chr"] + [a[3] for a in all_probes_srt]
    begins = ["begin"] + [str(a[4]) for a in all_probes_srt]
    ends = ["end"] + [str(a[5]) for a in all_probes_srt]
    bc_info = [[s] + sample_barcode_counts[s] for s in sorted(sample_names)]
    bc_out = np.transpose([uniq_mip_header] + [mip_header] + [copy_header] +                           [gene_names] + [chromosomes] + [begins] +                           [ends] + bc_info)
    write_list(bc_out, wdir + settings["barcodeCountFile"])

    multi_bc_info = [[s] + sample_barcode_counts_multi[s] for s in sorted(sample_names)]
    bc_out = np.transpose([uniq_mip_header] + [mip_header] + [copy_header] +                           [gene_names] + [chromosomes] + [begins] +                           [ends] + multi_bc_info)
    write_list(bc_out, wdir + settings["barcodeCountFile"] + "_multi")

    uniq_bc_info = [[s] + sample_barcode_counts_uniq[s] for s in sorted(sample_names)]
    bc_out = np.transpose([uniq_mip_header] + [mip_header] + [copy_header] +                           [gene_names] + [chromosomes] + [begins] +                           [ends] + uniq_bc_info)
    write_list(bc_out, wdir + settings["barcodeCountFile"] + "_uniq")

    with open(wdir + "temporary_counts", "wb") as outfile:
        pickle.dump(temporary_copy_counts, outfile)
    """
    with open(wdir + settings["perSampleResults"]) as infile:
        sample_res = json.load(infile)
    sample_names = sorted(sample_res.keys())
    """
    ann_keys = settings["annotationKeys"].split(";")
    annotation_id_key = settings["annotationIdKey"]
    sample_diffs = {}
    all_diff_locations = {}
    all_diffs = {}
    for s in sample_names:
        sample_diffs[s] = {}
        for g in sample_res[s]:
            for m in sample_res[s][g]:
                res = sample_res[s][g][m]
                for c in res:
                    filtered_data = res[c]["filtered_data"]
                    copy_barcode_total = 0
                    for f in filtered_data:
                        copy_barcode_total += f["barcode_count"]
                    copy_barcode_total = float(copy_barcode_total)
                    for f in filtered_data:
                        hid = f["haplotype_ID"]
                        bc = f["barcode_count"]
                        if bc < min_count or bc/copy_barcode_total < min_snp_frac:
                            continue
                        try:
                            qual = f["sequence_quality"]
                        except KeyError as e:
                            qual = f["seqence_quality"]
                        hap = haplotypes[m][hid]
                        if not haplotypes[m][hid]["mapped"]:
                            continue
                        copies = hap["mapped_copies"]
                        if len(copies) > 1:
                            multi_mapping = True
                        else:
                            multi_mapping = False
                        for c in copies:
                            k = ":".join([g, m, c])
                            location = all_probes_keys[k]
                            if multi_mapping:
                                bc_ratio = multi_mapper_ratios[s][hid][location]
                                bc = bc * bc_ratio
                            total_depth = sample_barcode_counts[s][location]
                            copy_call = call_info[g][m]["copies"][c]
                            copy_chrom = copy_call["chrom"]
                            copy_start = copy_call["capture_start"]
                            copy_end = copy_call["capture_end"]
                            copy_ori = copy_call["orientation"]
                            ref_seq = copy_call["capture_sequence"]
                            copy_differences = hap["mapped_copies"][c]["differences"]
                            for d in copy_differences:
                                normalized_key = d["vcf_normalized"]
                                var = normalized_key.split(":")
                                try:
                                    annotation_id = d["annotation"][annotation_id_key]
                                except KeyError as e:
                                    annotation_id = "."
                                hap_index = d["hap_index"]
                                start_index = min(hap_index)
                                end_index = max(hap_index) + 1
                                hap_qual_list = []
                                for hi in range(start_index, end_index):
                                    try:
                                        hap_qual_list.append(ord(qual[hi]) - 33)
                                    except IndexError as e:
                                        hap_qual_list = [-1]
                                        break
                                hap_qual = np.mean(hap_qual_list)
                                if hap_qual < min_snp_qual:
                                    continue
                                sam_diff = [normalized_key,
                                            var[0],
                                            int(var[1]),
                                            annotation_id,
                                            var[3],
                                            var[4],
                                            d["clinical"],
                                            d["clinical_id"],
                                            d["psv"],
                                            g,
                                            m,
                                            c,
                                            bc,
                                            total_depth,
                                            [hap_qual],
                                            [bc],
                                            [total_depth],
                                            [multi_mapping],
                                            [hid],
                                            [m],
                                            d["annotation"]]
                                try:
                                    all_diff_locations[normalized_key].append(location)
                                except KeyError as e:
                                    all_diff_locations[normalized_key]  = [location]
                                all_diffs[normalized_key] = sam_diff[:12] +                                                                          [sam_diff[-1]]
                                if not normalized_key in sample_diffs[s]:
                                    sample_diffs[s][normalized_key] = sam_diff
                                else:
                                    sample_diffs[s][normalized_key][12] += sam_diff[12]
                                    if sam_diff[19][0] not in sample_diffs[s][normalized_key][19]:
                                        sample_diffs[s][normalized_key][13] += sam_diff[13]
                                    for i in range(14,20):
                                        sample_diffs[s][normalized_key][i].extend(sam_diff[i])

    try:
        with open(wdir + settings["caseFile"]) as infile:
            case = {line.split("\t")[0] : line.strip().split("\t")[1] for line in infile}
    except IOError as e:
        case = {}

    for a in all_diff_locations:
        all_diff_locations[a] = sorted(list(set(all_diff_locations[a])))
    sorted_diffs = sorted(all_diff_locations,
                          key = lambda a: all_diff_locations[a][0])
    var_header = ["VKEY", "CHROM", "POS", "ID", "REF", "ALT",
                 "TARGET", "TID", "PSV", "PAR", "MIP", "CP"]
    var_header.extend(ann_keys)
    if settings["caseControlStats"] == "True":
        do_stats = True
    else:
        do_stats = False
    if do_stats:
        var_header.extend(["NS", "DP", "AD", "AF", "HO", "HE",
                           "HOCa", "HOCo", "HECa", "HECo",
                           "WTCa", "WTCo", "HOFishOR", "HOFishPval",
                           "HOChiPval", "HEFishOR", "HEFishPval",
                           "HEChiPval", "MTFishOR", "MTFishPval",
                           "MTChiPval", "FORMAT"])
    else:
        var_header.extend(["NS", "DP", "AD", "AF", "HO", "HE", "FORMAT"])
    var_header.extend(sample_names)
    outfile_list = [var_header]
    for d in sorted_diffs:
        diff_locs = all_diff_locations[d]
        diff_list = all_diffs[d][:-1]
        ann = all_diffs[d][-1]
        for k in ann_keys:
            diff_list.append(ann[k])
        wt_count = 0
        wt_case = 0
        wt_control = 0
        het_count = 0
        het_control = 0
        het_case = 0
        hom_count = 0
        hom_case = 0
        hom_control = 0
        diff_depth = 0
        allele_depth = 0
        sdiffs = []
        for s in sample_names:
            try:
                diff = sample_diffs[s][d][12:19]
                try:
                    diff[2] = np.mean(diff[2])
                except ValueError as e:
                    diff[2] = "."
                for i in range(3,7):
                    diff[i] = ",".join(map(str, diff[i]))
                if diff[0] == diff[1]:
                    GT = "1/1"
                    hom_count += 1
                    try:
                        if case[s] == "case":
                            hom_case += 1
                        elif case[s] == "control":
                            hom_control += 1
                    except KeyError as e:
                        pass
                else:
                    GT = "0/1"
                    het_count += 1
                    try:
                        if case[s] == "case":
                            het_case += 1
                        elif case[s] == "control":
                            het_control += 1
                    except KeyError as e:
                        pass
                diff = [GT] + diff
                try:
                    diff.append(case[s])
                except KeyError as e:
                    diff.append(".")
                sam_diff = ":".join(map(str, diff))
                diff_depth += diff[2]
                allele_depth += diff[1]
            except KeyError as e:
                total_depth = 0
                for l in diff_locs:
                    total_depth += sample_barcode_counts[s][l]
                diff = [0, total_depth, ".", 0, total_depth, ".", "."]
                diff_depth += total_depth
                if total_depth > 0:
                    GT = "0/0"
                    wt_count += 1
                    try:
                        if case[s] == "case":
                            wt_case += 1
                        elif case[s] == "control":
                            wt_control += 1
                    except KeyError as e:
                        pass
                else:
                    GT = "./."
                diff = [GT] + diff
                try:
                    diff.append(case[s])
                except KeyError as e:
                    diff.append(".")
                sam_diff = ":".join(map(str, diff))
            sdiffs.append(sam_diff)
        sam_count = sum([hom_count,
                         het_count,
                         wt_count])
        allele_count = hom_count + het_count
        if sam_count > 0:
            allele_freq = round(float(allele_count)/sam_count, 4)
        if do_stats:
            diff_info = [sam_count, diff_depth, allele_depth,
                         allele_freq, hom_count, het_count,
                        hom_case, hom_control,
                       het_case, het_control,
                      wt_case, wt_control]
            diff_info.extend(snp_stats(hom_case, hom_control,
                       het_case, het_control,
                      wt_case, wt_control))
        else:
            diff_info = [sam_count, diff_depth, allele_depth,
                         allele_freq, hom_count, het_count]
        diff_list.extend(diff_info)
        gt_format = ["GT", "AD", "DP", "BQ", "ADS", "DPS", "MM", "HID", "CS"]
        diff_list.append(":".join(gt_format))
        diff_list.extend(sdiffs)
        #outfile_list.append(u"\t".join(map(unicode, diff_list)))
        outfile_list.append(diff_list)
    '''
    sample_diffs_list = {}
    for s in sample_diffs:
        sample_diffs_list[s] = []
        for d in sample_diffs[s]:
            diff = copy.deepcopy(sample_diffs[s][d])
            for i in xrange(14,18):
                diff[i] = ",".join(map(str, diff[i]))
            for k in ann_keys:
                ann = diff.pop(-1)
                diff.append(ann[k])
                diff.append(ann)
            sample_diffs_list[s].append(diff)
        sample_diffs_list[s] = sorted(sample_diffs_list[s], key = itemgetter(1,2,0))
    outfile_list = ["\t".join(["sample",
                               "vcf_key",
                               "chromosome",
                                "position",
                                "rsid",
                                "reference",
                                "alt",
                                "target ID",
                                "PSV",
                                "paralog",
                                "mip",
                                "copy",
                               "allele depth",
                                "total depth",
                                "allele depths",
                                "total depths",
                                "multi map",
                                "haplotype IDs",
                                "targeted SNP"] + ann_keys)]
    for s in sorted(sample_diffs_list):
        var_list = sample_diffs_list[s]
        for diff in var_list:
            outfile_list.append("\t".join([s] + map(str, diff[:-1])))
    '''
    """
    with open(wdir + settings[u"variationTableFile"], u"w") as outfile:
        outfile.write("\n".join(map(map_str, outfile_list)))
    """
    write_list(outfile_list, wdir + settings["variationTableFile"])
    with open(wdir + settings["variationTableFile"] + ".json", "w") as outfile:
        json.dump(outfile_list, outfile)
    with open(wdir + settings["sampleVariationFile"], "w") as outfile:
        json.dump(sample_diffs, outfile)
    return
def filter_variation(variation_json,
                     min_barcode_count,
                     min_barcode_fraction,
                    samples_to_exclude = None,
                    samples_to_include = None):
    """
    Filter SNP variation for total minimum total depth and minimum
    barcode fraction for case control sample sets.
    """
    with open(variation_json) as infile:
        var_table = json.load(infile)
    # each item in the table list is a unique variation
    # except the first item which is the header
    header = copy.deepcopy(var_table[0])
    for h_ind in [0, 6, 7, 8, 9, 10, 11, 12]:
        header[h_ind] = "remove"
    exclude_indexes = []
    # sample IDs are index 40 and above,
    # remove those in the exclude list
    if samples_to_include is not None:
        for i in range(40, len(header)):
            if header[i] not in samples_to_include:
                exclude_indexes.append(i)
                header[i] = "remove"
    if samples_to_exclude is not None:
        for i in range(40, len(header)):
            if header[i] in samples_to_exclude:
                header[i] = "remove"
                exclude_indexes.append(i)
    header = [h for h in header if h != "remove"]
    header = header[:22] + ["HOR", "HPval"] + header[22:]

    filtered_var_table = [header]
    for var in var_table[1:]:
        wt = 0
        wt_case = 0
        wt_control = 0
        het = 0
        het_case = 0
        het_control = 0
        hom = 0
        hom_case = 0
        hom_control = 0
        AD = 0
        DP = 0
        NS = 0
        sample_infos = []
        # first 40 items in each variation list describes
        # general information for the variation such as chromosome
        # and base change, as well as population information such as
        # allele frequency in the sample set. Item at index 40 has
        # the first sample genotype
        for gen_index in range(40, len(var)):
            if gen_index in exclude_indexes:
                continue
            sam = var[gen_index]
            info = sam.split(":")
            # allele depth (mutant)
            ad = int(float(info[1]))
            # total depth at position
            dp= int(float(info[2]))
            # filter for minimum depth
            if dp< min_barcode_count:
                # change genotype to unknown (na)
                info[0] = "./."
            else:
                #AD += ad
                DP += dp
                NS += 1
                af = float(ad)/dp
                if af < min_barcode_fraction:
                    wt += 1
                    info[0] = "0/0"
                    if info[-1] == "case":
                        wt_case += 1
                    elif info[-1] == "control":
                        wt_control += 1
                elif af > 1 - min_barcode_fraction:
                    hom += 1
                    AD += ad
                    info[0] = "1/1"
                    if info[-1] == "case":
                        hom_case += 1
                    elif info[-1] == "control":
                        hom_control += 1
                else:
                    het += 1
                    AD += ad
                    info[0] = "0/1"
                    if info[-1] == "case":
                        het_case += 1
                    elif info[-1] == "control":
                        het_control += 1
            sample_infos.append(":".join(info))
        stats = snp_stats(hom_case, hom_control,
                           het_case, het_control,
                          wt_case, wt_control)
        try:
            fish = list(fisher_exact([[hom_case, hom_control],
                                 [het_case + wt_case,
                                  het_control + wt_control]]))
        except:
            fish = ["na", "na"]
        stats = fish + stats
        updated_var = copy.deepcopy(var[:18])
        for p_ind in [0, 6, 7, 8, 9, 10, 11, 12]:
            updated_var[p_ind] = "remove"
        updated_var = [uv for uv in updated_var if uv != "remove"]
        try:
            AF = (het + hom)/float(NS)
        except ZeroDivisionError as e:
            AF = 0
        updated_var.extend([NS, DP, AD, AF, hom, het,
                           hom_case, hom_control,
                           het_case, het_control,
                           wt_case, wt_control])
        updated_var.extend(stats)
        updated_var.append(var[39])
        updated_var.extend(sample_infos)
        if AF > 0:
            filtered_var_table.append(updated_var)
    write_list(filtered_var_table, variation_json + ".filtered")
    return
def variation_filter(variation_json, min_barcode_count, min_barcode_fraction):
    """
    Filter SNP variation for total minimum total depth and minimum
    barcode fraction for sample sets with no case control comparison.
    """
    with open(variation_json) as infile:
        var_table = json.load(infile)
    # each item in the table list is a unique variation
    # except the first item which is the header
    header = copy.deepcopy(var_table[0])
    for h_ind in [0, 6, 7, 8, 9, 10, 11, 12]:
        header[h_ind] = "remove"
    header = [h for h in header if h != "remove"]
    #header = header[:24] + ["HOR", "HPval"] + header[24:]
    filtered_var_table = [header]
    for var in var_table[1:]:
        wt = 0
        het = 0
        hom = 0
        hom_control = 0
        AD = 0
        DP = 0
        NS = 0
        sample_infos = []
        # first 40 items in each variation list describes
        # general information for the variation such as chromosome
        # and base change, as well as population information such as
        # allele frequency in the sample set. Item at index 40 has
        # the first sample genotype
        for sam in var[25:]:
            info = sam.split(":")
            # allele depth (mutant)
            ad = int(float(info[1]))
            # total depth at position
            dp= int(float(info[2]))
            # filter for minimum depth
            if dp< min_barcode_count:
                # change genotype to unknown (na)
                info[0] = "./."
            else:
                AD += ad
                DP += dp
                NS += 1
                af = float(ad)/dp
                if af < min_barcode_fraction:
                    wt += 1
                    info[0] = "0/0"
                elif af > 1 - min_barcode_fraction:
                    hom += 1
                    info[0] = "1/1"
                else:
                    het += 1
                    info[0] = "0/1"
            sample_infos.append(":".join(info))
        updated_var = copy.deepcopy(var[:18])
        for p_ind in [0, 6, 7, 8, 9, 10, 11, 12]:
            updated_var[p_ind] = "remove"
        updated_var = [uv for uv in updated_var if uv != "remove"]
        try:
            AF = (het + hom)/float(NS)
        except ZeroDivisionError as e:
            AF = 0
        updated_var.extend([NS, DP, AD, AF, hom, het])
        updated_var.append(var[24])
        updated_var.extend(sample_infos)
        if AF > 0:
            filtered_var_table.append(updated_var)
    write_list(filtered_var_table, variation_json + ".filtered.tsv")
    return
def map_str(s):
    try:
        return str(s).encode("utf-8")
    except UnicodeEncodeError as e:
        return s.encode("utf-8")
def write_list(alist, outfile_name):
    """ Convert values of a list to strings and save to file."""
    with open(outfile_name, "w") as outfile:
        outfile.write("\n".encode("utf-8").join(["\t".encode("utf-8").                                    join(map(map_str, l)) for l in alist])                      + "\n".encode("utf-8"))
    return
def snp_stats(hom_case, hom_control,
               het_case, het_control,
              wt_case, wt_control):
    """
    Given case/control genotype numbers in the order:
    1) number of homozygous cases,
    2) homozygous controls,
    3) heterozygous cases,
    4) heterozygous controls
    5) wildtype cases
    6) wildtype controls
    Returns a list of length 9:
    1-3) Homozygous mutant vs wildtype
    1) Odds ratio from Fisher's exact test
    2) P value from Fisher's
    3) P value from chi squared test
    4-6) Heterozygous mutant vs wildtype
    7-9) Mutants combined vs witdtype
    Errors return "na" in place of values
    """
    mut_case = het_case + hom_case
    mut_control = het_control + hom_control
    ho_v_wt = [[hom_case, hom_control],
               [wt_case, wt_control]]
    het_v_wt = [[het_case, het_control],
                [wt_case, wt_control]]
    mut_v_wt = [[mut_case, mut_control],
                [wt_case, wt_control]]
    output = []
    for tbl in [ho_v_wt, het_v_wt, mut_v_wt]:
        try:
            fish = fisher_exact(tbl)
        except:
            fish = ["na", "na"]
        try:
            chi = chi2_contingency(tbl)
        except:
            chi = ["na", "na", "na", "na"]
        output.extend(fish)
        output.append(chi[1])
    return output
def cnv_stats(hom_case, hom_control,
              wt_case, wt_control):
    """
    Given case/control genotype numbers in the order:
    1) number of homozygous cases,
    2) homozygous controls,
    3) heterozygous cases,
    4) heterozygous controls
    5) wildtype cases
    6) wildtype controls
    Returns a list of length 9:
    1-3) Homozygous mutant vs wildtype
    1) Odds ratio from Fisher's exact test
    2) P value from Fisher's
    3) P value from chi squared test
    4-6) Heterozygous mutant vs wildtype
    7-9) Mutants combined vs witdtype
    Errors return "na" in place of values
    """
    ho_v_wt = [[hom_case, hom_control],
               [wt_case, wt_control]]
    output = []
    tbl = ho_v_wt
    try:
        fish = fisher_exact(tbl)
    except:
        fish = ["na", "na"]
    try:
        chi = chi2_contingency(tbl)
    except:
        chi = ["na", "na", "na", "na"]
    output.extend(fish)
    output.append(chi[1])
    return output
def design_mips(design_dir, g):
    print(("Designing MIPs for ", g))
    try:
        Par = mod.Paralog(design_dir + g + "/resources/" + g + ".rinfo")
        Par.run()
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ", Par.paralog_name))
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
        Par = mod.Paralog(design_dir + g + "/resources/" + g + ".rinfo")
        Par.run()
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ", Par.paralog_name))
        if Par.chained_mips:
            print(("All MIPs are chained for paralog ", Par.paralog_name))
        else:
            print(("MIPs are NOT chained for paralog ", Par.paralog_name))
    except Exception as e:
        print((g, str(e), " FAILED!!!"))
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
        fastq_list.append("@" + f )
        fastq_list.append(fasta[f])
        fastq_list.append("+")
        fastq_list.append("H" * len(fasta[f]))
    with open(fastq_file, "w") as outfile:
        outfile.write("\n".join(fastq_list))
    return
def parasight (resource_dir,
               design_info_file,
               designed_gene_list = None,
               extra_extension = ".extra"):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    script_dir = get_file_locations()["all"]["script_DIR"]
    pdf_dir = script_dir + resource_dir + "pdfs/"
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash",
                     "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All             -sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = script_dir + design_info[t]["design_dir"] + t + "/" + t
        backup_name = basename + ".extra"
        filtered_name = basename + "_filtered.pse"
        backup_list.append("scp " + backup_name + " " + backup_name + ".bak")
        backup_list.append("mv " + filtered_name + " " + backup_name)
        psname = basename + ".01.01.ps"
        pdfname = basename + ".pdf"
        gs_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All -sOutputFile="            + pdfname + " " + psname
        if designed_gene_list != None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".pdf")
        else:
            pdf_convert_list.append(t + ".pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".pdf " +                            pdf_dir + t + ".pdf")
        outlist = ["parasight76.pl",
                   "-showseq",
                   basename + ".show",
                   "-extra",
                   basename + extra_extension,
                   "-template",
                   script_dir + "resources/nolabel.pst",
                   "-precode file:" + basename + ".precode" ,
                   "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode" , "w") as outfile:
            outfile.write("$opt{'filename'}='" + t +                          "';&fitlongestline; &print_all (0,'" + basename + "')")
    with open(resource_dir + "backup_commands", "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(resource_dir + "parasight_commands", "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(resource_dir + "gs_commands", "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(resource_dir + "copy_commands", "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(resource_dir + "convert_commands", "w") as outfile:
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
    with open(resource_dir + "visualize", "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return
def parasight_mod (resource_dir,
               design_info_file,
                  species,
               designed_gene_list = None,
               extra_extension = ".extra",
                  maf = 0.1,
                  height = 200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    script_dir = get_file_locations()["all"]["script_DIR"]
    pdf_dir = script_dir + resource_dir + "mod_pdfs/"
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash",
                     "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All             -sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = script_dir + design_info[t]["design_dir"] + t + "/" + t
        showname = basename + ".show"
        try:
            with open(showname) as infile:
                sln = infile.readlines()[-1].strip().split("\t")
                show_region = sln[0] + ":" + sln[2] + "-" + sln[3]
        except IOError as e:
            continue
        reg_snps = get_snps(show_region,
                           get_file_locations()[species]["snps"])
        indels_low = []
        indels_high = []
        snvs_low = []
        snvs_high = []
        for rsnp in reg_snps:
            acs = [a for a in rsnp[23].split(",") if a != ""]
            allele_count = list(map(int, list(map(float, acs))))
            # allele with max count
            max_all = max(allele_count)
            # total alleles
            tot_all = sum(allele_count)
            # minor allele freq
            min_af = (tot_all - max_all)/float(tot_all)
            if "-" in rsnp[22].split(","):
                if min_af >= maf:
                    indels_high.append(rsnp)
                else:
                    indels_low.append(rsnp)
            elif min_af >= maf:
                snvs_high.append(rsnp)
            else:
                snvs_low.append(rsnp)
        backup_name = basename + ".extra"
        try:
            with open(backup_name) as infile, open(backup_name + ".mod", "w") as outfile:
                thirds = 0
                rd_range = list(range(height))
                for line in infile:
                    newline = line.split("\t")
                    if newline[3] in ["all_targets", "target"]:
                        newline[5] = "-10"
                        newline[6] = "4"
                        outfile.write("\t".join(newline))
                    elif newline[3] in ["capture", "extension", "ligation"]:
                        if (thirds % 3)== 0:
                            rd = random.choice(rd_range)
                        thirds += 1
                        newline[5] = str(-30 - rd)
                        outfile.write("\t".join(newline))
                    elif newline[3] == "snp":
                        pass
                    else:
                        outfile.write(line)
                outfile.write("\n")
                for col, snp_list, snp_type, ofs in zip(["pink", "red",
                                     "light green", "dark green"],
                                   [indels_low, indels_high,
                                   snvs_low, snvs_high],
                                                  ["low frequency indel",
                                                   "high frequency indel",
                                                  "low frequency SNP",
                                                  "high frequency SNP"],
                                                       ["-27", "-27",
                                                       "-26", "-26"]):
                    for snp in snp_list:
                        ol = snp[1:4]
                        ol.extend([snp_type, col, ofs, "4",
                                  ""])
                        outfile.write("\t".join(ol) + "\n")
        except IOError as e:
            continue
        psname = basename + ".01.01.ps"
        pdfname = basename + ".mod.pdf"
        gs_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All -sOutputFile="            + pdfname + " " + psname
        if designed_gene_list != None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".mod.pdf")
        else:
            pdf_convert_list.append(t + ".mod.pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".mod.pdf " +                            pdf_dir + t + ".mod.pdf")
        outlist = ["parasight76.pl",
                   "-showseq",
                   basename + ".show",
                   "-extra",
                   basename + extra_extension + ".mod",
                   "-template",
                   script_dir + "resources/nolabel.pst",
                   "-precode file:" + basename + ".precode" ,
                   "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode" , "w") as outfile:
            outfile.write("$opt{'filename'}='" + t +                          "';&fitlongestline; &print_all (0,'" + basename + "')")
    with open(resource_dir + "backup_commands", "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(resource_dir + "parasight_commands", "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(resource_dir + "gs_commands", "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(resource_dir + "copy_commands", "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(resource_dir + "convert_commands", "w") as outfile:
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
    with open(resource_dir + "visualize_mod", "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return
def parasight_shift (resource_dir,
               design_info_file,
                  species,
               designed_gene_list = None,
               extra_extension = ".extra",
                  maf = 0.1,
                  height = 200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    script_dir = get_file_locations()["all"]["script_DIR"]
    pdf_dir = script_dir + resource_dir + "mod_pdfs/"
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash",
                     "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All             -sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = script_dir + design_info[t]["design_dir"] + t + "/" + t
        showname = basename + ".show"
        try:
            with open(showname) as infile:
                sln = infile.readlines()[-1].strip().split("\t")
                show_region = sln[0] + ":" + sln[2] + "-" + sln[3]
        except IOError as e:
            continue
        backup_name = basename + ".extra"
        try:
            with open(backup_name) as infile, open(backup_name + ".mod", "w") as outfile:
                thirds = 0
                rd_range = list(range(height))
                for line in infile:
                    newline = line.split("\t")
                    if newline[3] in ["all_targets", "target"]:
                        newline[5] = "-10"
                        newline[6] = "4"
                        outfile.write("\t".join(newline))
                    elif newline[3] in ["capture", "extension", "ligation"]:
                        if (thirds % 3)== 0:
                            rd = random.choice(rd_range)
                        thirds += 1
                        newline[5] = str(-30 - rd)
                        outfile.write("\t".join(newline))
                    else:
                        outfile.write(line)
                outfile.write("\n")
        except IOError as e:
            continue
        psname = basename + ".01.01.ps"
        pdfname = basename + ".mod.pdf"
        gs_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All -sOutputFile="            + pdfname + " " + psname
        if designed_gene_list != None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".mod.pdf")
        else:
            pdf_convert_list.append(t + ".mod.pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".mod.pdf " +                            pdf_dir + t + ".mod.pdf")
        outlist = ["parasight76.pl",
                   "-showseq",
                   basename + ".show",
                   "-extra",
                   basename + extra_extension + ".mod",
                   "-template",
                   script_dir + "resources/nolabel.pst",
                   "-precode file:" + basename + ".precode" ,
                   "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode" , "w") as outfile:
            outfile.write("$opt{'filename'}='" + t +                          "';&fitlongestline; &print_all (0,'" + basename + "')")
    with open(resource_dir + "backup_commands", "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(resource_dir + "parasight_commands", "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(resource_dir + "gs_commands", "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(resource_dir + "copy_commands", "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(resource_dir + "convert_commands", "w") as outfile:
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
    with open(resource_dir + "visualize_mod", "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return
def parasight_print(gene_list, extra_suffix = ".extra"):
    for g in gene_list:
        print(("cd ../" + g))
        print(("parasight76.pl -showseq "+ g + ".show "
               + "-extra " + g + extra_suffix))
def rescue_mips(design_dir,
                paralog_name,
                redesign_list,
               same_strand_overlap,
               opposite_strand_overlap,
               low,
               high,
                mip_limit,
               chain,
               incomp_score):
    print(("Redesigning MIPs for ", paralog_name))
    paralog_file = design_dir + paralog_name + "/" + paralog_name
    with open(paralog_file, "rb") as infile:
        par = pickle.load(infile)
    subprocess.call(["scp",
                     paralog_file,
                     paralog_file + ".last"])
    par.extra_mips = {}
    redesign_pairs = []
    for segment_name in par.segments:
        seg = par.segments[segment_name]
        """
        seg.rinfo["SELECTION"]["compatibility"]\
        ["low"] = low
        seg.rinfo["SELECTION"]["compatibility"]\
        ["high"] = high
        seg.rinfo["SELECTION"]["compatibility"]\
        ["chain"] = chain

        """

        seg.rinfo["SELECTION"]["compatibility"]        ["low"] = low
        seg.rinfo["SELECTION"]["compatibility"]        ["high"] = high
        seg.rinfo["SELECTION"]["compatibility"]        ["mip_limit"] = mip_limit
        for subregion_name in seg.subregions:
            sub = seg.subregions[subregion_name]
            sub.score_mips()
            temp_scored = copy.deepcopy(sub.primers["hairpin"]["dictionary"])
            best_mipset = sub.mips["best_mipset"]["dictionary"]["mips"]
            scored = list(sub.mips["scored_filtered"]["dictionary"]                            ["pair_information"].keys())
            keep_keys = scored + list(best_mipset.keys())
            for m in list(temp_scored['pair_information'].keys()):
                if m not in keep_keys:
                    temp_scored['pair_information'].pop(m)
            compatible = compatibility(temp_scored,
                                       overlap_same = same_strand_overlap,
                                    overlap_opposite = opposite_strand_overlap)
            alt = sub.mips["best_mipset"]["dictionary"]["alternatives"] = {}
            hairpin = sub.mips["hairpin"]
            compat_info = compatible["pair_information"]
            for m in best_mipset:
                if best_mipset[m].fullname in redesign_list:
                    if m not in redesign_pairs:
                        redesign_pairs.append(m)
                    capture_start = sub.mips["best_mipset"]["dictionary"]                    ["mips"][m].mip["C0"]["capture_start"]
                    capture_end = sub.mips["best_mipset"]["dictionary"]                    ["mips"][m].mip["C0"]["capture_end"]
                    capture_target = [capture_start, capture_end]
                    alt_count = 0
                    alt[m] = {}
                    for p1 in scored:
                        if p1 != m and p1 not in redesign_list:
                            #print p1
                            start1 = hairpin[p1].mip["C0"]["capture_start"]
                            end1 = hairpin[p1].mip["C0"]["capture_end"]
                            capture1 = [start1, end1]
                            if len(overlap(capture_target, capture1)) > 0:
                                remaining_target = subtract_overlap(
                                    [capture_target],
                                    [capture1])
                                if len(remaining_target) == 1:
                                    for p2 in scored:
                                        if p2 != p1 and p2 != m and                                        p1 not in redesign_list and                                        (p2 not in compat_info[p1]["incompatible"]):
                                            start2 = hairpin[p2].mip["C0"]                                            ["capture_start"]
                                            end2 = hairpin[p2].mip["C0"]                                            ["capture_end"]
                                            capture2 = [start2, end2]
                                            if len(subtract_overlap(
                                                [capture_target],
                                                [capture2])) == 1 and\
                                               len(overlap(remaining_target[0],
                                                               capture2)) > 0:
                                                uncovered = subtract_overlap(
                                                    remaining_target,
                                                    [capture2])
                                                if len(uncovered) == 0:
                                                    alt[m][str(alt_count)] = [p1,p2]
                                                    alt_count += 1

            new_set = list(best_mipset.keys())
            for m in best_mipset:
                if best_mipset[m].fullname in redesign_list:
                    print(("Rescue MIPs for ", best_mipset[m].fullname))
                    coverage_needed = [hairpin[m].mip["C0"]                                            ["capture_start"],
                                             hairpin[m].mip["C0"]\
                                            ["capture_end"]]
                    print(("Region to cover is ", coverage_needed))
                    best_rescue = []
                    best_score = -50000
                    for p in alt[m]:
                        p1, p2 = alt[m][p]
                        incomp = len(set(compat_info[p1]["incompatible"] +                                      compat_info[p2]["incompatible"])                        .intersection(new_set))
                        score = 0
                        score += hairpin[p1].tech_score
                        score += hairpin[p2].tech_score
                        score += incomp_score * incomp
                        if score > best_score:
                            best_score = score
                            best_rescue = [p1, p2]
                            if p1 in best_mipset:
                                best_rescue = [p2]
                            if p2 in best_mipset:
                                best_rescue = [p1]
                    resc_capture = []
                    for resc in best_rescue:
                        resc_capture.append([hairpin[resc].mip["C0"]                                            ["capture_start"],
                                             hairpin[resc].mip["C0"]\
                                            ["capture_end"]])
                    print(("Rescue mip coverage is ", merge_overlap(resc_capture)))
                    print(("Rescue score is ", best_score))
                    new_set.extend(best_rescue)
            sub.extra_mips = extra_mips = {}
            for t in new_set:
                if t not in best_mipset:
                    extra_mips[t] = hairpin[t]
            par.extra_mips.update(sub.extra_mips)
    locus_info = par.locus_info
    selected_mips = par.selected_mips
    mip_names_ordered = sorted(selected_mips,
                               key=lambda mip_key: selected_mips[mip_key].\
                                                    mip["C0"]["mip_start"])
    name_counter = int(selected_mips[mip_names_ordered[-1]]                       .fullname.split("_")[-1][3:]) + 1
    extra_names_ordered = sorted(extra_mips,
                               key=lambda mip_key: extra_mips[mip_key].\
                                                    mip["C0"]["mip_start"])
    outfile_list = []
    for mip_name in extra_names_ordered:
        m = extra_mips[mip_name]
        fullname_subregion = m.subregion.fullname
        fullname_mip = fullname_subregion + "_mip" + str(name_counter)
        m.fullname = fullname_mip
        # get coordinate information of the mip
        # reference copy information will be used
        m.mip_start = m.mip["C0"]["mip_start"]
        m.mip_end = m.mip["C0"]["mip_end"]
        m.capture_start = m.mip["C0"]["capture_start"]
        m.capture_end = m.mip["C0"]["capture_end"]
        m.orientation = m.mip["C0"]["orientation"]
        m.chromosome = m.mip["C0"]["chrom"]
        for key in m.mip_dic["mip_information"]:
            if key == "ref":
                fullname_extension = "_ref"
            elif key.startswith("alt"):
                fullname_extension = "_" + key
            else:
                continue
            # ["#pair_name", "mip_name", "chrom", "mip_start", "capture_start",
            # "capture_end", "mip_end", "orientation", "tech_score", "func_score","mip_sequence",
            # "unique_captures", "must_captured", "captured_copies"]
            outlist = [m.name, m.fullname + fullname_extension,
                       m.chromosome, m.mip_start,
                       m.capture_start, m.capture_end,
                       m.mip_end, m.orientation,
                        m.tech_score, m.func_score,
                        m.mip_dic["mip_information"][key]["SEQUENCE"]]
            locus_info["mips"].append(outlist)
            outfile_list.append(outlist)
        name_counter += 1
    for mipname in redesign_pairs:
        m = selected_mips[mipname]
        for key in m.mip_dic["mip_information"]:
            if key == "ref":
                fullname_extension = "_ref"
            elif key.startswith("alt"):
                fullname_extension = "_" + key
            else:
                continue
            outlist = ["#" + m.name, m.fullname + fullname_extension,
                       m.chromosome, m.mip_start,
                           m.capture_start, m.capture_end,
                       m.mip_end, m.orientation,
                            m.tech_score, m.func_score,
                            m.mip_dic["mip_information"][key]["SEQUENCE"]]
            outfile_list.append(outlist)
    par.print_info()
    with open(paralog_file + ".rescued", "wb") as outfile:
        pickle.dump(par, outfile)
    with open(design_dir + paralog_name + "/" + paralog_name +              "_rescue.mips", "w") as outfile:
        outfile.write("\n".join(["\t".join(map(str, outlist)) for                                 outlist in outfile_list]))
    return
def analyze_data_old(settings_file):
    settings = get_analysis_settings(settings_file)
    group_samples(settings)
    update_raw_data(settings)
    get_counts(settings)
    get_unique_probes(settings)
    create_data_table(settings)
    filter_tables(settings)
    get_summary_table(settings)
    return
def process_data_old(wdir, run_ids):
    for rid in run_ids:
        settings_file = wdir + "settings_" + rid
        get_data(settings_file)
        analyze_data_old(settings_file)
    return
def get_data(settings_file):
    settings = get_analysis_settings(settings_file)
    # data extraction
    #filter_mipster (settings)
    get_haplotypes(settings)
    align_haplotypes(settings)
    parse_aligned_haplotypes(settings)
    update_aligned_haplotypes(settings)
    update_unique_haplotypes(settings)
    update_variation(settings)
    get_raw_data(settings)
    return
def analyze_data(settings_file):
    settings = get_analysis_settings(settings_file)
    group_samples(settings)
    update_raw_data(settings)
    get_counts(settings)
    return
def process_data(wdir, run_ids):
    for rid in run_ids:
        settings_file = wdir + "settings_" + rid
        get_data(settings_file)
        analyze_data(settings_file)
    return
def process_haplotypes(settings_file):
    settings = get_analysis_settings(settings_file)
    get_haplotypes(settings)
    align_haplotypes(settings)
    parse_aligned_haplotypes(settings)
    update_aligned_haplotypes(settings)
    update_unique_haplotypes(settings)
    update_variation(settings)
    return
def process_results(wdir,
                   settings_file,
                    sample_sheets = None,
                    meta_files = [],
                    targets_file = None,
                    target_join = "union"
                   ):
    settings = get_analysis_settings(wdir + settings_file)
    if sample_sheets is None:
        sample_sheets = [wdir + "samples.tsv"]
    ##########################################################
    ##########################################################
    # Process 1: use sample sheets, sample sets and meta files
    # to determine which data points from the mipster file
    # should be used, print some statistics
    ##########################################################
    ##########################################################
    # process sample sheets
    run_meta = pd.concat(
        [pd.read_table(s)
         for s in sample_sheets],
         ignore_index = True
    )
    # if only a subset of the run is to be used for this analysis
    # create a sample/probe sets dataframe for used
    # samples and probes , otherwise
    # use all samples in the sample sheets
    run_meta["sample_name"] = (
            run_meta["sample_name"].astype(str)
        )
    run_meta["Sample Name"] = run_meta["sample_name"]
    run_meta["Sample ID"] = run_meta[["sample_name",
                                      "sample_set",
                                      "replicate"]].apply(
        lambda a: "-".join(map(str, a)), axis = 1
    )
    # Sample Set key is reserved for meta data
    # but sometimes erroneously included in the
    # sample sheet. It should be removed.
    try:
        run_meta.drop("Sample Set",
                 inplace = True,
                 axis = 1)
    except (ValueError, KeyError):
        pass
    # drop duplicate values originating from
    # multiple sequencing runs of the same libraries
    run_meta = run_meta.drop_duplicates()
    run_meta = run_meta.groupby(
        ["Sample ID", "Library Prep"]
    ).first().reset_index()
    run_meta.to_csv(wdir + "run_meta.csv")
    # load meta data for samples
    sample_meta_list = []
    try:
        sample_meta = pd.concat(
            [pd.read_table(f) for f in meta_files],
            join = "outer",
            ignore_index = True
        )
    except ValueError as e:
        # if no meta files were provided, create a dummy
        # meta dataframe
        sample_meta = copy.deepcopy(run_meta[["Sample Name"]])
        sample_meta["Meta"] = "Meta"
    sample_meta["Sample Name"] = sample_meta["Sample Name"].astype(str)
    sample_meta = sample_meta.groupby(["Sample Name"]).first().reset_index()
    # Merge Sample meta data and run data
    merged_meta = pd.merge(run_meta, sample_meta,
                           on = "Sample Name",
                           how = "inner")
    merged_meta.to_csv(wdir + "merged_meta.csv")
    print(("{} out of {} samples has meta information and"
           " will be used for analysis.").format(
        merged_meta.shape[0], run_meta.shape[0]
    ))
    # get used sample's ids
    sample_ids = merged_meta["Sample ID"].unique().tolist()
    ##########################################################
    ##########################################################
    # Process 2: extract all observed variants from observed
    # haplotypes and create a variation data frame that will
    # be able to map haplotype IDs to variation.
    ##########################################################
    ##########################################################
    # get all haplotypes and variants observed in the data
    # from the haplotype dictionary
    hap_file = settings["haplotypeDictionary"]
    with open(wdir + hap_file) as infile:
        haplotypes = json.load(infile)
    # keep all variant in all haplotypes in
    # variation list
    variation_list = []
    # keep haplotypes that are the same as reference
    # genome in the reference list
    reference_list = []
    # annotation ID Key specifies if there is and ID field in the vcf
    # which has a database ID of the variation at hand. For example,
    # rsid for variation already defined in dbSNP.
    annotation_id_key = settings["annotationIdKey"]
    unmapped = 0
    for m in haplotypes:
        g = m.split("_")[0]
        for hid in haplotypes[m]:
            hap = haplotypes[m][hid]
            # skip off target haplotypes
            if not hap["mapped"]:
                unmapped += 1
                continue
            copies = hap["mapped_copies"]
            # check if the haplotype is mapping to
            # multiple locations in genome
            if len(copies) > 1:
                multi_mapping = True
            else:
                multi_mapping = False
            for c in copies:
                copy_differences = hap["mapped_copies"][c]["differences"]
                # go through all differences from reference genome
                # get a subset of information included in the
                # haplotype dictionary
                if len(copy_differences) == 0:
                    reference_list.append([hid, c, multi_mapping])
                for d in copy_differences:
                    # all variation is left normalized to reference genome
                    # this is done to align all indels to the same start
                    # to avoid having different locations for the same
                    # indel in a tandem repeat region.
                    # each variation is given a unique key, which is
                    # formed by the first 4 fields of vcf (chr:pos:id:ref:alt)
                    normalized_key = d["vcf_normalized"]
                    var = normalized_key.split(":")
                    raw_key = d["vcf_raw"]
                    raw_var = raw_key.split(":")
                    # get the position of variation prior to
                    # left normalization
                    original_pos = int(raw_var[1])
                    # indels are represented with the preceeding base
                    # like A:AG for an insertion and AG:A for deletion
                    # in both cases, the position is the preceding base
                    # in some cases where the indel is right after probe
                    # arm, we may not actually have coverage in the position
                    # indicated here, so change the position to the next base
                    # where the real change is
                    if len(raw_var[4]) != len(raw_var[3]):
                        original_pos += 1
                    # get the annotation id if any, such as rsID
                    try:
                        annotation_id = d["annotation"][annotation_id_key]
                    except KeyError as e:
                        annotation_id = "."
                    # get the location of variation relative
                    # to haplotype sequence
                    hap_index = d["hap_index"]
                    start_index = min(hap_index)
                    end_index = max(hap_index) + 1
                    temp_list = [normalized_key,
                                            var[0],
                                            int(var[1]),
                                            annotation_id,
                                            var[3],
                                            var[4],
                                            d["psv"],
                                            g, m, c, hid,
                                           raw_key,
                                           original_pos,
                                           start_index,
                                           end_index,
                                           multi_mapping,
                                            ]
                    try:
                        for ak in annotation_keys:
                            temp_list.append(d["annotation"][ak])
                    except NameError as e:
                        annotation_keys = list(d["annotation"].keys())
                        for ak in annotation_keys:
                            temp_list.append(d["annotation"][ak])
                    variation_list.append(temp_list)
    # create pandas dataframes for variants
    colnames = ["VKEY", "CHROM", "POS", "ID", "REF", "ALT",
                "PSV", "Gene", "MIP", "Copy", "Haplotype ID",
                "RAW_VKEY", "Original Position", "Start Index",
                "End Index", "Multi Mapping"]
    clean_annotation_keys = []
    for ak in annotation_keys:
        if ak.startswith("AAChange."):
            clean_annotation_keys.append("AAChangeClean")
        elif ak.startswith("ExonicFunc."):
            clean_annotation_keys.append("ExonicFunc")
        elif ak.startswith("Gene."):
            clean_annotation_keys.append("RefGene")
        else:
            clean_annotation_keys.append(ak)
    colnames = colnames + clean_annotation_keys
    variation_df = pd.DataFrame(variation_list,
                               columns = colnames)
    # create pandas dataframe for reference haplotypes
    reference_df = pd.DataFrame(reference_list,
                               columns = ["Haplotype ID",
                                          "Copy",
                                          "Multi Mapping"])

    # create a dataframe for all mapped haplotypes
    mapped_haplotype_df = pd.concat(
        [variation_df.groupby(["Haplotype ID",
                               "Copy",
                              "Multi Mapping"]).first(
            ).reset_index()[["Haplotype ID",
                               "Copy",
                              "Multi Mapping"]],
         reference_df],
        ignore_index = True)
    print("There are {mh.shape[0]} mapped and {um} unmapped (off target) haplotypes.".format(
        mh = mapped_haplotype_df,
        um = unmapped
    ))
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
    raw_results = pd.read_table(wdir + settings["mipsterFile"])
    raw_results["Haplotype ID"] = raw_results["haplotype_ID"] + "-0"
    # limit the results to the samples intended for this analysis
    raw_results = raw_results.loc[
        raw_results["sample_name"].isin(sample_ids)
    ]
    mapped_results = raw_results.merge(mapped_haplotype_df,
                                         how = "inner")
    print(("There are {rr.shape[0]} data points in raw data,"
            " {mr.shape[0]} are mapped to genome.").format(
        rr = raw_results,
        mr = mapped_results
    ))
    # rename some columns for better visualization
    mapped_results.rename(
        columns = {"sample_name": "Sample ID",
                  "mip_name": "MIP",
                  "gene_name": "Gene",
                  "barcode_count": "Barcode Count",
                  "read_count": "Read Count"},
        inplace = True
    )
    # Try to estimate the distribution of data that is mapping
    # to multiple places in the genome.
    # This is done in 4 steps.
    # 1) Get uniquely mapping haplotypes and barcode counts
    unique_df = mapped_results.loc[~mapped_results["Multi Mapping"]]
    unique_table = pd.pivot_table(unique_df,
               index = "Sample ID",
              columns = ["Gene", "MIP", "Copy"],
              values = ["Barcode Count"],
              aggfunc = np.sum)
    # 2) Estimate the copy number of each paralog gene
    # for each sample from the uniquely mapping data
    try:
        average_copy_count = float(settings["averageCopyCount"])
        norm_percentiles = list(map(float,
                              settings["normalizationPercentiles"]))
    except KeyError as e:
        average_copy_count = 2
        norm_percentiles = [0.4, 0.6]
    unique_df.loc[:, "CA"] = average_copy_count
    unique_df.loc[:, "Copy Average"] = average_copy_count
    unique_df.loc[:, "Adjusted Barcode Count"] = unique_df["Barcode Count"]
    unique_df.loc[:, "Adjusted Read Count"] = unique_df["Read Count"]
    unique_table.fillna(0, inplace = True)
    copy_counts = get_copy_counts(unique_table,
                            average_copy_count,
                            norm_percentiles)
    # 3) Estimate the copy number of each "Gene"
    # from the average copy count of uniquely mapping
    # data for all MIPs within the gene
    cc = copy_counts.groupby(level = ["Gene",
                                 "Copy"], axis = 1).sum()
    gc = copy_counts.groupby(level = ["Gene"], axis = 1).sum()
    ac = cc/gc
    multi_df = mapped_results.loc[mapped_results["Multi Mapping"]]
    # 4) Distribute multi mapping data proportional to
    # Paralog's copy number determined from the
    # uniquely mapping data
    if not multi_df.empty:
        mca = multi_df.apply(
            lambda r: get_copy_average(r, ac),
            axis = 1)
        multi_df.loc[mca.index, "Copy Average"] = mca
        mca = multi_df.groupby(["Sample ID",
                         "Gene"])["Copy Average"].transform(
            normalize_copies
        )
        multi_df.loc[mca.index, "CA"] = mca
        multi_df.loc[: , "Adjusted Barcode Count"] = (
            multi_df.loc[:, "Barcode Count"]
            * multi_df.loc[:, "CA"]
        )
        multi_df.loc[: , "Adjusted Read Count"] = (
            multi_df.loc[:, "Read Count"]
            * multi_df.loc[:, "CA"]
        )
    # Combine unique and multimapping data
    combined_df = pd.concat([unique_df,
                           multi_df],
                           ignore_index = True)
    combined_df.rename(
        columns = {
            "Barcode Count": "Raw Barcode Count",
            "Adjusted Barcode Count": "Barcode Count",
            "Read Count": "Raw Read Count",
            "Adjusted Read Count": "Read Count"
        },
        inplace = True
    )
    # print total read and barcode counts
    print(("Total number of reads and barcodes were {0[0]} and {0[1]}."
          " On target number of reads and barcodes were {1[0]} and {1[1]}.").format(
        raw_results[["read_count",
                     "barcode_count"]].sum(),
         combined_df[["Read Count",
                     "Barcode Count"]].sum().astype(int)
    ))

    ##########################################################
    ##########################################################
    # Process 4: Combine per sample information from process 3
    # with variant and haplotype information from process 2.
    # filter results by given criteria.
    ##########################################################
    ##########################################################
    # Add the statistics for each haplotype to the data
    # such as how many samples had a given haplotype
    # and how many barcodes supported a given haplotype
    # Filter the haplotypes for those criteria to
    # remove possible noise and infrequent haplotypes
    ####### Haplotype Filters #############
    haplotype_min_barcode_filter = int(settings["minHaplotypeBarcodes"])
    haplotype_min_sample_filter = int(settings["minHaplotypeSamples"])
    haplotype_min_sample_fraction_filter = float(settings["minHaplotypeSampleFraction"])
    ####### Haplotype Filters #############
    hap_counts = combined_df.groupby(
        "Haplotype ID"
    )["Barcode Count"].sum().reset_index().rename(
        columns = {"Barcode Count": "Haplotype Barcodes"})
    hap_sample_counts = combined_df.groupby("Haplotype ID")["Sample ID"].apply(
        lambda a: len(set(a))
    ).reset_index(
    ).rename(columns = {"Sample ID": "Haplotype Samples"})
    num_samples = float(combined_df["Sample ID"].unique().size)
    hap_sample_counts["Haplotype Sample Fraction"] = (
        hap_sample_counts["Haplotype Samples"] /num_samples
    )
    hap_counts = hap_counts.merge(hap_sample_counts)
    hap_counts = hap_counts.loc[(hap_counts["Haplotype Samples"]
                                >= haplotype_min_sample_filter)
                               &(hap_counts["Haplotype Sample Fraction"]
                                >= haplotype_min_sample_fraction_filter)
                               &(hap_counts["Haplotype Barcodes"]
                                >= haplotype_min_barcode_filter)]
    variation_df = variation_df.merge(hap_counts,
                                     how = "inner")
    # Rename or remove some columns for downstream analysis
    variation_df["AA Change"] = variation_df["AAChangeClean"].apply(
        split_aa
    )
    variation_df["AA Change Position"] = variation_df["AAChangeClean"].apply(
        split_aa_pos
    )
    try:
        variation_df.drop(["Chr", "Ref", "Alt"],
                         axis = 1,
                         inplace = True)
    except ValueError as e:
        pass

    # if there is a targets file, observed variation can be filtered
    # using the targets. There are 4 ways that the targets data
    # can be added to the variation data.
    # "intersection": filter the observed data to targeted data only,
    # remove targets not observed and data not targeted.
    # "targets": filter the observed data to targeted data only,
    # remove data not targeted but keep targets not observed
    # "data": add the target information to observed data,
    # excluding targets not observed
    # "union"; add the target information to observed data,
    # including the unobserved targets
    # keys ["Vkey", "Chrom", "Pos", "Id", "Ref", "Alt"]
    # must be provided if targets or union method is to be used for
    # merging. These will be used for corresponding columns
    # for targets that are not observed
    if targets_file is not None:
        targets = pd.read_table(targets_file)
        join_dict = {"intersection": "inner",
                "union": "outer",
                "targets": "right",
                "data": "left"}
        targets["Targeted"] = "Yes"
        variation_df = variation_df.merge(
            targets,
            how = join_dict[target_join]
        )
        variation_df["Targeted"].fillna("No",
                                   inplace = True)
        # If a reference genome locus is a mutation of interest
        # such as dhps-437, this information can be supplied
        # in targets file. The rest of the variants will be
        # assinged a False value for this.
        try:
            variation_df["Reference Resistant"].fillna(
                "No", inplace = True
            )
            ref_resistant = True
        except KeyError as e:
            ref_resistant = False
            #variant_counts["Reference Resistant"] = False

        if target_join in ["union", "targets"]:
            data_keys = ["VKEY", "CHROM", "POS", "ID", "REF", "ALT"]
            target_keys = ["Vkey", "Chrom", "Pos", "Id", "Ref", "Alt"]
            for dk, tk in zip(data_keys, target_keys):
                variation_df[dk].fillna(variation_df[tk],
                                       inplace = True)
            variation_df.drop(target_keys,
                             axis = 1,
                             inplace = True)
    else:
        variation_df["Targeted"] = "No"
        ref_resistant = False
    # each variant needs to have a name. This should be provided in
    # the target file. For those variant that are not in the target
    # file, we'll create their names by adding aa-change to gene name
    try:
        variation_df["Mutation Name"].fillna(variation_df["Gene"] + "-"
                                             + variation_df["AA Change"],
                                     inplace = True)
    except KeyError as e:
        variation_df["Mutation Name"] = (variation_df["Gene"] + "-"
                                             + variation_df["AA Change"])
    variation_df.loc[variation_df["AA Change"] == ".",
            "Mutation Name"] = variation_df.loc[
                        variation_df["AA Change"] == "."].apply(
                                                rename_noncoding,
                                                axis = 1)

    # remove columns that will not be used after this point
    variation_df.drop(
        ["RAW_VKEY",
         "Original Position",
         "Haplotype Barcodes",
         "Haplotype Samples",
         "Haplotype Sample Fraction"],
        axis = 1,
        inplace = True
    )

    # Create a chrom, pos tuple to use for the location
    # of a given variant, to be used in coverage calculations
    variation_df["Coverage Key"] = variation_df.apply(
        lambda a: (a["CHROM"], a["POS"]),
        axis = 1
    )

    # load the "call info" dictionary that has
    # all mip information such as the capture coordinates
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    # load the probe set dictionary to extract the
    # probes that were used in this run
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    probe_cop = []
    for m in used_probes:
        g = m.split("_")[0]
        try:
            for c in call_info[g][m]["copies"]:
                probe_cop.append([m,c])
        except KeyError as e:
            continue
    probe_cop = pd.DataFrame(probe_cop,
                            columns = ["MIP", "Copy"])
    # add a place holder column for merging probe information
    probe_cop["Temp"] = "Temp"
    # perform outer merge on results and used probes
    # to include probes that had no coverage in the results
    combined_df = combined_df.merge(probe_cop,
                                    how = "outer").drop("Temp",
                                                       axis = 1)
    # Fill NA values for probes with no coverage
    combined_df["Sample ID"].fillna("Temp",
                                   inplace = True)
    combined_df["Haplotype ID"].fillna(
        combined_df["MIP"] + ".0-0",
        inplace = True
    )
    combined_df["Barcode Count"].fillna(
        0, inplace = True
    )

    # Add sample and barcode depth information for each
    # variant
    variant_counts = combined_df[["Haplotype ID",
                "sequence_quality",
                "Sample ID",
                "Barcode Count",
                "Copy"]].merge(variation_df, how = "right")
    # For unobserved variants, we need a place holder for
    # Sample ID
    variant_counts["Sample ID"].fillna("Temp",
                                       inplace = True)
    variant_counts["Barcode Count"].fillna(0,
                                          inplace = True)

    # Get the sample and barcode depth stats for each variant
    # and filter for given thresholds.
    var_counts = variant_counts.groupby("VKEY").agg(
        {"Sample ID": lambda a: len(set(a)),
        "Barcode Count": "sum",
        "Targeted": lambda a: "Yes" if "Yes" in set(a) else "No"}
    ).rename(
        columns = {"Sample ID": "Variant Samples",
                   "Barcode Count": "Variant Barcodes"}
    ).fillna(0).reset_index()
    var_counts["Variant Sample Fraction"] = var_counts[
        "Variant Samples"
    ].transform(lambda a: a/num_samples)
    # filter variants for specified criteria
    variant_min_barcode_filter = int(settings["minVariantBarcodes"])
    variant_min_sample_filter = int(settings["minVariantSamples"])
    variant_min_sample_fraction_filter = float(settings["minVariantSampleFraction"])
    var_counts = var_counts.loc[((var_counts["Variant Samples"]
                                 >= variant_min_sample_filter)
                               &(var_counts["Variant Barcodes"]
                                 >= variant_min_barcode_filter)
                               &(var_counts["Variant Sample Fraction"]
                                 >= variant_min_sample_fraction_filter))
                               | (var_counts["Targeted"] == "Yes")]
    print("There were {} total and {} unique variants, ".format(
        variant_counts.shape[0],
        len(variant_counts["VKEY"].unique())
    ))
    variant_counts = variant_counts.merge(var_counts,
                                         how = "inner").drop(
        ["Variant Samples",
        "Variant Barcodes",
        "Variant Sample Fraction"],
        axis = 1
    )
    print(("{} total and {} unique variants remain after "
           "filtering variants for "
           "minimum total barcodes of {}, "
           "minimum observed sample number of {}, "
           "and minimum observed sample fraction of {}.").format(
        variant_counts.shape[0],
        len(variant_counts["VKEY"].unique()),
        variant_min_barcode_filter,
        variant_min_sample_filter,
        variant_min_sample_fraction_filter
    ))
    # Add variant sequence quality
    def get_qual(row):
        try:
            # get start of the variation relative to haplotype sequence
            start_index = int(row["Start Index"])
            end_index = int(row["End Index"])
            qual = row["sequence_quality"]
            hap_qual_list = []
            for hi in range(start_index, end_index):
                try:
                    # get phred quality of each base in variation
                    # and convert the phred score to number
                    hap_qual_list.append(ord(qual[hi]) - 33)
                except IndexError as e:
                    continue
                    break
            # calculate quality as the mean for multi base variation
            if len(hap_qual_list) == 0:
                return np.nan
            else:
                return np.mean(hap_qual_list)
        except:
            return np.nan
    variant_counts["Variation Quality"] = variant_counts.apply(
        get_qual, axis = 1
    )

    # filter variants for sequence quality
    variant_min_quality = int(settings["minVariantQuality"])
    variant_counts = variant_counts.loc[(variant_counts["Variation Quality"].isnull())
                      |(variant_counts["Variation Quality"]
                       >= variant_min_quality)]
    print(("{} total and {} unique variants remained after "
           "quality filtering for phred scores >= {}.").format(
        variant_counts.shape[0],
        len(variant_counts["VKEY"].unique()),
        variant_min_quality
    ))
    ##########################################################
    ##########################################################
    # Process 5: Calculate coverage per mip and aslo per
    # variant per sample.
    ##########################################################
    ##########################################################
    # create a position to MIP dictionary for each variant
    # that holds which MIPs cover a given position for
    # coverage calculations
    cpos = variant_counts.groupby(["CHROM", "POS"]).first().reset_index()[["CHROM", "POS"]]
    cpos = cpos.apply(lambda a: (a["CHROM"], a["POS"]), axis = 1).values.tolist()
    position_to_mip = {}
    # go through found variants and add any MIP associated
    # with a given variant
    for m in haplotypes:
        for hid in haplotypes[m]:
            hap = haplotypes[m][hid]
            if not hap["mapped"]:
                continue
            copies = hap["mapped_copies"]
            for c in copies:
                copy_differences = hap["mapped_copies"][c]["differences"]
                # go through all differences from reference genome
                # get a subset of information included in the
                # haplotype dictionary
                for d in copy_differences:
                    # all variation is left normalized to reference genome
                    # this is done to align all indels to the same start
                    # to avoid having different locations for the same
                    # indel in a tandem repeat region.
                    # each variation is given a unique key, which is
                    # formed by the first 4 fields of vcf (chr:pos:id:ref:alt)
                    normalized_key = d["vcf_normalized"]
                    var = normalized_key.split(":")
                    var_pos = (var[0], int(var[1]))
                    if var_pos in cpos:
                        try:
                            position_to_mip[var_pos].add(
                                (m, c)
                            )
                        except KeyError as e:
                            position_to_mip[var_pos] = set()
                            position_to_mip[var_pos].add(
                                (m, c)
                            )
    # add any additional MIP that covers the variant positions
    for var_pos in cpos:
        for g in call_info:
            for m in call_info[g]:
                if m in used_probes:
                    for c in call_info[g][m]["copies"]:
                        ch = call_info[g][m]["copies"][c]["chrom"]
                        cs = call_info[g][m]["copies"][c]["capture_start"]
                        ce = call_info[g][m]["copies"][c]["capture_end"]
                        if ((var_pos[0] == ch)
                            and (cs <= var_pos[1] <= ce)):
                            try:
                                position_to_mip[var_pos].add(
                                    (m, c)
                                )
                            except KeyError as e:
                                position_to_mip[var_pos] = set()
                                position_to_mip[var_pos].add(
                                    (m, c)
                                )
    # Create pivot table of combined barcode counts
    # This is a per MIP per sample barcode count table
    # of the samples with sequencing data
    barcode_counts = pd.pivot_table(combined_df,
                                    index = "Sample ID",
                                    columns = ["MIP",
                                               "Copy"],
                                    values = ["Barcode Count"],
                                    aggfunc = np.sum)
    print("There are {} samples with sequence data".format(
        barcode_counts.shape[0]
    ))
    # barcode count data is only available for samples with data
    # so if a sample has not produced any data, it will be missing
    # these samples should be added with 0 values for each probe
    bc_cols = barcode_counts.columns
    bc_cols = [bc[1:] for bc in bc_cols]
    temp_meta = merged_meta[["Sample ID",
                 "replicate"]].append({"Sample ID": "Temp",
                                     "replicate": 1},
                                     ignore_index = True)
    barcode_counts = pd.merge(temp_meta.set_index("Sample ID"),
                         barcode_counts,
                         left_index = True,
                         right_index = True,
                         how = "left").drop("Temp")
    barcode_counts.drop("replicate", axis = 1, inplace = True)
    barcode_counts.columns = pd.MultiIndex.from_tuples(bc_cols,
                                                      names = ["MIP",
                                                              "Copy"])
    barcode_counts.fillna(0, inplace = True)
    print("There are {} total samples.".format(barcode_counts.shape[0] - 1))

    # create a coverage dictionary for each variant position
    # for each sample
    bc_dict = barcode_counts.to_dict(orient = "index")
    cov_dict = {}
    for ch, po in position_to_mip:
        for m, cp in position_to_mip[(ch,po)]:
            for s in bc_dict:
                try:
                    cov_dict[(s, ch, po)] +=  bc_dict[s][(m,cp)]
                except KeyError as e:
                    cov_dict[(s, ch, po)] =  bc_dict[s][(m,cp)]
    # Include samples with no variants (probably no sequence data)
    # in the variant counts
    variant_counts = variant_counts.merge(
        temp_meta,
        how = "outer"
    ).drop("replicate",
           axis = 1)
    # Columns that will be used for pivoting cannot
    # have NA values, so for samples with no data
    # we need to fill those values with a place holder
    for c in ["CHROM", "POS", "ID", "REF", "ALT"]:
        variant_counts[c].fillna("Temp", inplace = True)
    variant_counts["Barcode Count"].fillna(
        0, inplace = True
    )
    # define a function that returns coverage at given
    # coverage key (sample, chromosome, position)
    def return_coverage(k):
        try:
            return cov_dict[k]
        except KeyError as e:
            return 0
    # create a vcf file for all variants
    # create pivot table with each variant having its own column
    vcf_table = variant_counts.pivot_table(
        index = "Sample ID",
        columns = ["CHROM", "POS", "ID", "REF", "ALT"],
        values = "Barcode Count",
        aggfunc = "sum",
    ).drop("Temp")
    # remove place holder values for samples with no data
    try:
        vcf_table.drop(("Temp", "Temp", "Temp", "Temp", "Temp"),
                      axis = 1,
                      inplace = True)
    except KeyError as e:
        pass
    vcf_table.fillna(0, inplace = True)
    # create coverage table for the vcf_table
    v_cols = vcf_table.columns
    v_index = vcf_table.index
    vcf_co = pd.DataFrame([
        [return_coverage((s, c[0], c[1]))
         for c in v_cols]
        for s in v_index],
        index = v_index,
        columns = v_cols)
    # merge variants on position to get non-reference
    # allele count
    vcf_non_ref = vcf_table.groupby(level = ["CHROM", "POS"],
                     axis = 1).transform("sum")
    # calculate reference allele counts
    vcf_ref = vcf_co - vcf_non_ref
    # get variant qualities
    try:
        vcf_quals = variant_counts.pivot_table(
            index = "Sample ID",
            columns = ["CHROM", "POS", "ID", "REF", "ALT"],
            values = "Variation Quality",
            aggfunc = "mean",
        ).drop("Temp").fillna(".")
    except KeyError:
        vcf_quals = variant_counts.pivot_table(
            index = "Sample ID",
            columns = ["CHROM", "POS", "ID", "REF", "ALT"],
            values = "Variation Quality",
            aggfunc = "mean",
        ).fillna(".")
    # calculate allele frequencies and
    # create genotype calls from frequencies
    # no filtering will be applied here so
    # even low frequency non-refs mixed with ref
    # will be a HET call. This is only for vcf file
    # to be filtered by proper vcf tools later.
    vcf_freq = vcf_table/vcf_co
    vcf_gen = vcf_freq.applymap(lambda a: "0/0" if a == 0
                     else "." if (np.isnan(a) or np.isinf(a))
                      else "0/1" if a <1
                      else "1/1")
    # merge all vcf tables to create the merged vcf
    vcf = (vcf_gen + ":" + vcf_ref.applymap(int).applymap(str)
           + "," + vcf_table.applymap(int).applymap(str)
     + ":" + vcf_co.applymap(int).applymap(str)
     + ":" + vcf_quals.applymap(str)
    ).fillna(".:0,0:0:.").T.sort_index(
        level = ["CHROM", "POS"]
    )
    vcf_samples = vcf.columns.tolist()
    vcf.columns = vcf_samples
    vcf["QUAL"] = "."
    vcf["FILTER"] = "."
    vcf["INFO"] = "."
    vcf["FORMAT"] = "GT:AD:DP:SQ"
    vcf = vcf[["QUAL", "FILTER", "INFO", "FORMAT"] + vcf_samples]
    vcf_header = ["##fileformat=VCFv4.2",
    '##ALT=<ID=NON_REF,Description="Represents a possible alternative allele at this location">"',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in that order">',
     '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
     '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">']
    vcf_file = "variants.vcf"
    with open(wdir + vcf_file, "w") as outfile:
        outfile.write("\n".join(vcf_header) + "\n")
        vcf.reset_index().rename(columns = {"CHROM": "#CHROM"}).to_csv(
            outfile, index = False, sep = "\t"
        )

    # create a variant table from the variant counts
    # this will create a samples by variants table
    # all keys that will be used in the pivot table
    # must be non-NA so fill them with a temporary
    # value
    variant_counts["Mutation Name"].fillna("Temp",
                                      inplace = True)
    variant_counts["CHROM"].fillna("Temp",
                                      inplace = True)
    variant_counts["POS"].fillna("Temp",
                                      inplace = True)
    variant_counts["ID"].fillna("Temp",
                                      inplace = True)
    variant_counts["REF"].fillna("Temp",
                                      inplace = True)
    variant_counts["ALT"].fillna("Temp",
                                      inplace = True)
    variant_counts["Gene"].fillna("Temp",
                                      inplace = True)
    variant_counts["Mutation Name"].fillna("Temp",
                                      inplace = True)
    variant_counts["AA Change Position"].fillna("Temp",
                                      inplace = True)
    variant_counts["ExonicFunc"].fillna("Temp",
                                      inplace = True)
    variant_counts["Targeted"].fillna("No",
                                      inplace = True)

    if ref_resistant:
        variant_counts["Reference Resistant"].fillna("No",
                                                inplace = True)
        # create pivot table for each unique variant
        variant_table = variant_counts.pivot_table(
            index = "Sample ID",
            columns = ["CHROM", "POS", "ID", "REF", "ALT",
                       "Gene", "Mutation Name", "AA Change Position",
                       "ExonicFunc",
                       "Reference Resistant", "Targeted"],
            values = "Barcode Count",
            aggfunc = "sum"
        )
        # drop the temporary sample place holder, if any
        try:
            variant_table.drop("Temp", inplace = True)
        except KeyError:
            pass
        # drop the column with temporary place holders, if any
        try:
            variant_table = variant_table.drop(
                ("Temp", "Temp", "Temp", "Temp", "Temp", "Temp",
                 "Temp", "Temp", "Temp", "No", "No"), axis = 1)
        except KeyError as e:
            pass
        # if a sample did not have a variant, the table value
        # will be NA. Change those to 0.
        variant_table.fillna(0, inplace = True)
        # add amino acid positions and sort table
        col_list = []
        for c in variant_table.columns:
            pos = c[7].split("-")[-1].split("_")[0]
            if pos != ".":
                positions = []
                num_found = False
                for dig in pos:
                    try:
                        int(dig)
                        positions.append(dig)
                        num_found = True
                    except ValueError as e:
                        if num_found:
                            break
                pos = int("".join(positions))
            col_list.append(c[:7] + (pos,) + c[7:])
        column_names = variant_table.columns.names
        new_cols = pd.MultiIndex.from_tuples(
            col_list,
            names = column_names[:7] + ["AA Position"] + column_names[7:]
        )
        variant_table.columns = new_cols
        variant_table.sort_index(level = ["Gene",
                                         "AA Position"],
                                axis = 1)
        variant_table.columns = variant_table.columns.droplevel(level = "AA Position")
        # get coverage table with same indexes as
        # the variant table
        v_cols = variant_table.columns
        v_index = variant_table.index
        variant_cov_df = pd.DataFrame([
            [return_coverage(
                (s, c[0], c[1])
            )
            for c in v_cols]
            for s in v_index],
            index = v_index,
            columns = v_cols)
        # define nonsynonamous changes to aggregate all non-reference
        # amino acids and get to all reference amino acid calls from
        # there
        nonsyn = list(
            set(variant_counts["ExonicFunc"]).difference([".", "synonymous SNV", "Temp"])
        )
        idx = pd.IndexSlice
        # aggregate all nonsyn calls for each amino acid position
        # this is not ideal, as it ignores indels but this is only
        # used for loci where reference genome is actually mutant
        # that is so far only dhps-437 and there are no common indels
        # in the vicinity.
        non_ref_aa_table = variant_table.loc[: , idx[:,:,:,:,:,:,:,:,
                                  nonsyn, :,:]].groupby(level = ["Gene",
                                                         "AA Change Position"],
                                                       axis = 1).transform("sum")
        # non_ref_aa_table loses the synonymous variants in the
        # previous step. We create an all-zero table with
        # variant table indexes by subtracting it from itself
        # than add non_ref table to get the non_ref values
        non_ref_aa_table = (variant_table
                            - variant_table
                            + non_ref_aa_table).fillna(0)
        non_ref_aa_table = non_ref_aa_table.groupby(
            level = ["Gene",
                     "AA Change Position"],
            axis = 1).transform(max)
        non_ref_aa_table = non_ref_aa_table.groupby(
            level = [
                "Gene",
                "Mutation Name",
                "Reference Resistant",
                "Targeted",
                "ExonicFunc"],
            axis = 1).max()
        # create a like-indexed coverage table
        coverage_aa_table = variant_cov_df.groupby(
            level = [
                "Gene",
                "Mutation Name",
                "Reference Resistant",
                "Targeted",
                "ExonicFunc"],
            axis = 1).max()
        # calculate reference amino acid counts
        ref_aa_table = coverage_aa_table - non_ref_aa_table
        # aggregate all variants that lead to the
        # same amino acid change
        mutant_aa_table = variant_table.groupby(
            level = [
                "Gene",
                "Mutation Name",
                 "Reference Resistant",
                 "Targeted",
                 "ExonicFunc"],
            axis = 1).sum()
        # do a sanity check for all the grouping and coverage calculations
        # none of the table values for mutant or reference tables can be
        # larger than the coverage for a given locus.
        if (((mutant_aa_table - coverage_aa_table) > 0).sum().sum()
            + ((ref_aa_table - coverage_aa_table) > 0).sum().sum()) > 0:
            print("Some loci have lower coverage than mutation calls!")
        # where reference is the variant of interest("Reference Resistant")
        # change mutant count to reference count
        mutant_aa_table.loc[:, idx[:, :, "Yes", :]] = ref_aa_table.loc[:, idx[:, :, "Yes", :]]
        mutant_aa_table.columns = mutant_aa_table.columns.droplevel("Reference Resistant")
        coverage_aa_table.columns = coverage_aa_table.columns.droplevel("Reference Resistant")
    else:
        # create pivot table for each unique variant
        variant_table = variant_counts.pivot_table(
            index = "Sample ID",
            columns = ["CHROM", "POS", "ID", "REF", "ALT",
                       "Gene", "Mutation Name", "AA Change Position",
                       "ExonicFunc", "Targeted"],
            values = "Barcode Count",
            aggfunc = "sum"
        ).drop("Temp") # drop the temporary sample place holder
        # drop the column with temporary place holders, if any
        try:
            variant_table = variant_table.drop(
                ("Temp", "Temp", "Temp", "Temp", "Temp", "Temp",
                 "Temp", "Temp", "Temp", "No"), axis = 1)
        except KeyError as e:
            pass
        # if a sample did not have a variant, the table value
        # will be NA. Change those to 0.
        variant_table.fillna(0, inplace = True)
        # get coverage table with same indexes as
        # the variant table
        v_cols = variant_table.columns
        v_index = variant_table.index
        variant_cov_df = pd.DataFrame([
            [return_coverage(
                (s, c[0], c[1])
            )
            for c in v_cols]
            for s in v_index],
            index = v_index,
            columns = v_cols)
        # aggregate all variants that lead to the
        # same amino acid change
        mutant_aa_table = variant_table.groupby(
            level = [
                "Gene",
                "Mutation Name",
                 "Targeted",
                 "ExonicFunc"],
            axis = 1).sum()
        # create a like-indexed coverage table
        coverage_aa_table = variant_cov_df.groupby(
            level = [
                "Gene",
                "Mutation Name",
                "Targeted",
                "ExonicFunc"],
            axis = 1).max()
        # do a sanity check for all the grouping and coverage calculations
        # none of the table values for mutant or reference tables can be
        # larger than the coverage for a given locus.
        if (((mutant_aa_table - coverage_aa_table) > 0).sum().sum()) > 0:
            print("Some loci have lower coverage than mutation calls!")

    combined_df.to_csv(wdir + "haplotype_counts.csv",
                      index = False)
    variant_counts.to_csv(wdir + "variants.csv",
                         index = False)
    barcode_counts.to_csv(wdir + "barcode_counts.csv")
    plot_performance(barcode_counts,
                    wdir = wdir,
                    save = True)
    variant_table.to_csv(wdir + "variant_table.csv")
    variant_cov_df.to_csv(wdir + "variant_coverage_table.csv")
    mutant_aa_table.to_csv(wdir + "mutant_table.csv")
    coverage_aa_table.to_csv(wdir + "mutant_coverage.csv")

    min_mutation_count = int(settings["minMutationCount"])
    min_mutation_fraction = float(settings["minMutationFraction"])
    min_coverage = int(settings["minCoverage"])

    mutant_aa_table = mutant_aa_table.applymap(
        lambda a: 0 if a < min_mutation_count else a
    )
    coverage_aa_table = coverage_aa_table.applymap(
        lambda a: 0 if a < min_coverage else a
    )

    mutant_freq_table = (
        mutant_aa_table/coverage_aa_table
    ).replace(np.inf, np.nan)
    mutant_freq_table.to_csv(wdir + "mutation_frequencies.csv")

    genotypes = mutant_freq_table.applymap(
        lambda x: 2 if x >= (1 - min_mutation_fraction)
        else 1 if x >= min_mutation_fraction else np.nan if np.isnan(x) else 0
    )
    genotypes.to_csv(wdir + "genotypes.csv")
    print(("Per sample mutation frequencies have been "
          "calculated for mutants with at least {} supporting "
          "barcodes and loci with at least {} coverage. "
          "Loci with less coverage will have NA frequencies "
          "and mutants with less supporting barcodes have been "
          "reset to zero frequency.\n Genotypes have been called "
          " using the frequency values. Genotypes with <{} "
           "frequency have been reset to 0 (WT).").format(
        min_mutation_count,
        min_coverage,
        min_mutation_fraction
    ))
    sample_counts = combined_df.groupby("Sample ID")[["Read Count",
                                        "Barcode Count"]].sum()
    target_cov = pd.concat(
        [(barcode_counts>= 1).sum(axis = 1),
         (barcode_counts>= 5).sum(axis = 1),
         (barcode_counts>= 10).sum(axis = 1)],
        axis = 1,
    ).rename(
        columns = {
            0: "targets_with_1_barcodes",
            1: "targets_with_5_barcodes",
            2: "targets_with_10_barcodes"
    }
    )
    sample_counts = sample_counts.merge(target_cov,
                                       how = "outer",
                                       left_index = True,
                                       right_index = True).fillna(0)
    sample_counts.to_csv(wdir + "sample_summary.csv")

    return
def combine_sample_data(gr):
    result = {}
    result["barcode_count"] = gr["barcode_count"].sum()
    result["read_count"] = gr["read_count"].sum()
    result["sequence_quality"] = gr.sort_values("barcode_count",
                       ascending = False)["sequence_quality"].iloc[0]
    result["mip_name"] = gr["mip_name"].iloc[0]
    result["gene_name"] = gr["gene_name"].iloc[0]
    return pd.Series(result)
def combine_info_files(wdir,
                       settings_file,
                       info_files,
                       sample_sheets,
                       combined_file,
                       sample_sets = None):
    settings = get_analysis_settings(wdir + settings_file)
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
        current_run_meta["capital_set"] = current_run_meta[
            "sample_set"
        ].apply(str.upper)
        current_run_meta["Original SID"] = current_run_meta[["sample_name",
                                      "sample_set",
                                      "replicate"]].apply(
            lambda a: "-".join(a), axis = 1
        )
        run_meta.append(current_run_meta)
    run_meta = pd.concat(run_meta,
                        ignore_index = True)
    if sample_sets is not None:
        for s in sample_sets:
            s.append("Temp")
        sps = pd.DataFrame(sample_sets,
                          columns = ["sample_set",
                                    "probe_set",
                                    "Temp"])
    else:
        sps = run_meta.groupby(
            ["sample_set", "probe_set"]
        ).first().reset_index()[["sample_set", "probe_set"]]
        sps["Temp"] = "Temp"
    run_meta = run_meta.merge(sps, how = "inner").drop(
        "Temp", axis = 1
    )
    run_meta_collapsed = run_meta.groupby(
        ["sample_name",
          "capital_set",
          "replicate",
          "Library Prep"]).first().reset_index()[["sample_name",
                                  "capital_set",
                                  "replicate",
                                  "Library Prep"]]
    run_meta_collapsed["new_replicate"] = run_meta_collapsed.groupby(
        "sample_name"
    )["replicate"].transform(lambda g: list(map(str, list(range(1, len(g) + 1)))))
    run_meta = run_meta.merge(run_meta_collapsed)
    run_meta["Sample ID"] = run_meta[["sample_name",
                                  "capital_set",
                                  "new_replicate"]].apply(
        lambda a: "-".join(a), axis = 1
    )
    for i in range(len(info_files)):
        i_file = info_files[i]
        current_run_meta = run_meta.loc[run_meta["sheet_order"] == i]
        current_run_dict = current_run_meta.set_index(
            "Original SID"
        ).to_dict(orient = "index")
        line_number = 0
        try:
            dump = gzip.open(i_file, "rb").readline()
            inf_file = gzip.open(i_file, "rb")
        except IOError as e:
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
                else:
                    ori_sample_id = newline[si_index]
                    try:
                        library = current_run_dict[ori_sample_id]["Library Prep"]
                        sample_id = current_run_dict[ori_sample_id]["Sample ID"]
                        d = ([newline[ci] if ci != si_index else sample_id
                              for ci in col_indexes] + [library])
                        data.append(d)
                    except KeyError as e:
                        continue
    info = pd.DataFrame(data,
                        columns = c_vals + ["Library Prep"])
    info["barcode_count"] = info["barcode_count"].astype(int)
    info["read_count"] = info["read_count"].astype(int)
    info = info.groupby(
        ["sample_name",
        "haplotype_sequence",
        "Library Prep"]
    ).apply(
        combine_sample_data
    ).reset_index()
    m_groups = info.groupby("mip_name")
    h_list = []
    for m, g in m_groups:
        md = pd.DataFrame(g.groupby(["mip_name",
                                    "haplotype_sequence"]).size().sort_values(
        ascending = False).reset_index()).reset_index()
        md["index"] = md["index"].astype(str)
        md["haplotype_ID"] = md["mip_name"] + "." + md["index"]
        h_list.append(md[["haplotype_sequence", "haplotype_ID"]])
    hap_ids = pd.concat(h_list, ignore_index = True)
    info = info.merge(hap_ids)
    info.to_csv(wdir + combined_file,
               index = False,
               sep = "\t")
    run_meta = run_meta.groupby("Sample ID").first().reset_index()
    run_meta = run_meta.drop(["Sample ID",
                              "sample_set",
                             "sheet_order",
                             "replicate"],
                            axis = 1).rename(
            columns = {"capital_set": "sample_set",
                      "new_replicate": "replicate"})
    run_meta.to_csv(wdir + "samples.tsv",
                           sep = "\t",
                          index = False)
def update_probe_sets(mipset_table = "resources/mip_ids/mipsets.csv",
                     mipset_json = "resources/mip_ids/probe_sets.json"):
    mipsets = pd.read_csv(mipset_table)
    mipset_list = mipsets.to_dict(orient="list")
    mipset_dict = {}
    for mipset in mipset_list:
        mlist = mipset_list[mipset]
        mipset_dict[mipset] = [m for m in mlist if not pd.isnull(m)]
    with open(mipset_json, "w") as outfile:
        json.dump(mipset_dict, outfile, indent = 1)
    return
def generate_mock_fastqs(settings_file):
    """
    Generate fastq files for each sample. These files will have stitched and
    barcode corrected reads.
    """
    settings = get_analysis_settings(settings_file)
    wdir = settings["workingDir"]
    sample_results_file = settings["perSampleResults"]
    haplotype_file = settings["haplotypeDictionary"]
    fastq_dir = wdir + "fastq/"
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    with open(wdir + sample_results_file) as infile:
        sample_results = json.load(infile)
    with open(wdir + haplotype_file) as infile:
        haplotypes = json.load(infile)
    for sample in sample_results:
        with gzip.open(fastq_dir + sample + ".fq.gz", "w") as outfile:
            for g in sample_results[sample]:
                for m in sample_results[sample][g]:
                    for c in sample_results[sample][g][m]:
                        filtered_data = sample_results[sample][g][m][c]["filtered_data"]
                        for fd in filtered_data:
                            bc = fd["barcode_count"]
                            hid = fd["haplotype_ID"]
                            qual = fd["sequence_quality"]
                            seq = haplotypes[m][hid]["sequence"]
                            counter = 0
                            for i in range(bc):
                                read_name = "_".join(["@", sample, m, c, str(i)])
                                fastq_lines = "\n".join([read_name, seq, "+", qual]) + "\n"
                                outfile.write(fastq_lines)
    return
def generate_fastqs(wdir, mipster_files, min_bc_count, min_bc_frac):
    """
    Generate fastq files for each sample. These files will have stitched and
    barcode corrected reads.
    """
    fastq_dir = wdir + "fastq/"
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    mipster_dfs = pd.concat([pd.read_table(wdir + mfile,
                                          usecols = [
                                              "s_Sample",
                                              'h_popUID',
                                              "h_seq",
                                              'c_qual',
                                              'c_barcodeCnt',
                                              "c_barcodeFrac"
                                          ])
                             for mfile in mipster_files],
                           axis = 0,
                           ignore_index = True)
    mipster = mipster_dfs.loc[(mipster_dfs["c_barcodeCnt"] >= min_bc_count)
                              &(mipster_dfs["c_barcodeFrac"] >= min_bc_frac)].groupby(
        "s_Sample").apply(lambda x: pd.DataFrame.to_dict(
        x, orient = "index"
    )).to_dict()
    for sample in mipster:
        with gzip.open(fastq_dir + sample + ".fq.gz", "w") as outfile:
            outfile_list = []
            for ind in mipster[sample]:
                row = mipster[sample][ind]
                bc = int(row["c_barcodeCnt"])
                hid = row["h_popUID"]
                qual = row["c_qual"]
                seq = row["h_seq"]
                sample = row["s_Sample"]
                counter = 0
                for i in range(bc):
                    read_name = "_".join(["@", sample, hid, str(ind), str(i)])
                    outfile_list.extend([read_name, seq, "+", qual])
            outfile.write("\n".join(outfile_list) + "\n")
    return
def make_vcf(settings):
    """
    Create a vcf file from variation table to be used
    in vcf friendly programs.
    """
    wdir = settings["workingDir"]
    with open(wdir + settings["variationTableFile"]
         + ".json") as infile:
        var_tab = json.load(infile)
    var_head = var_tab[0]
    var_data = var_tab[1:]
    head_cols = ["CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT"]
    ann_index = list(range(var_head.index("CP") + 1,
                      var_head.index("NS")))
    format_index = var_head.index("FORMAT")
    gen_index = list(range(var_head.index("FORMAT") + 1,
                     len(var_head)))
    info_index = list(range(var_head.index("NS"),
                      var_head.index("HE") + 1))
    col_index = list(range(var_head.index("CHROM"),
                     var_head.index("ALT") + 1))
    mip_index = list(range(var_head.index("TARGET"),
                     var_head.index("PAR") + 1))
    mip_index.append(var_head.index("VKEY"))
    # Allele depth field in the variation table file only shows
    # the alt allele depth. VCF specs want REF, ALT counts.
    # REF count is inferred from total - ALT here, so this
    # will not work for multiallelic sites such as microsatellites.
    updated_var = []
    for v in var_data:
        vd = []
        for i in col_index:
            vd.append(v[i])
        vd.extend([".", "."])
        info = []
        for ind in [info_index, ann_index, mip_index]:
            for i in ind:
                field_name = var_head[i]
                field_value = v[i]
                info.append("=".join(map(str,
                    [field_name, field_value])))
        vd.append(";".join(info))
        format_text = v[format_index]
        vd.append(format_text)
        for i in gen_index:
            gen = v[i].split(":")
            try:
                ad = int(float(gen[1]))
                dp = int(float(gen[2]))
                wd = dp - ad
                ad = str(wd) + "," + str(ad)
                gen[1] = ad
            except ValueError as e:
                pass
            gen = ":".join(gen)
            vd.append(gen)
        updated_var.append(vd)
    # get the meta information lines from ANNOVAR output
    with open(wdir + settings["annotationOutput"]
             + "." + settings["annotationBuildVersion"]
             + "_multianno.vcf") as infile:
        meta_info_lines = []
        for line in infile:
            if line.startswith("##"):
                meta_info_lines.append(line.strip())
            else:
                break
    # add meta information lines for the MIP pipeline fields
    met1 = ['##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##INFO=<ID=AD,Number=A,Type=Integer,Description="Alternate Allele Depth">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, as determined by # heterozygous or homozygous samples divided by NS">',
    '##INFO=<ID=HO,Number=A,Type=Integer,Description="Number of Samples Homozygous for ALT allele">',
    '##INFO=<ID=HE,Number=A,Type=Integer,Description="Number of Samples Heterozygous for ALT Allele">',
           '##INFO=<ID=TARGET,Number=1,Type=String,Description="Whether the position was targeted.">',
    '##INFO=<ID=TID,Number=1,Type=String,Description="Target ID if the position was targeted.">',
    '##INFO=<ID=PSV,Number=1,Type=String,Description="Whether the position was determined as a paralog specific locus.">',
    '##INFO=<ID=PAR,Number=1,Type=String,Description="Gene or region name from MIP pipeline.">',
    '##INFO=<ID=VKEY,Number=1,Type=String,Description="Unique ID given to each variation.">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype.">',
           '##FORMAT=<ID=AD,Number=R,Type=Float,Description="Read depth for each allele, starting with REF.">',
           '##FORMAT=<ID=DP,Number=1,Type=Float,Description=" Read depth at this position for this sample.">',
           '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Base quality for alt allele.">',
           '##FORMAT=<ID=ADS,Number=.,Type=Float,Description="Alt allele read depths when supported by multiple MIPs.">'
           '##FORMAT=<ID=DPS,Number=.,Type=Float,Description="Total read depths when supported by multiple MIPs.">'
           '##FORMAT=<ID=MM,Number=.,Type=String,Description="Whether the allele mapping to multiple locations on reference genome.">',
           '##FORMAT=<ID=HID,Number=.,Type=String,Description="Haplotype ID supporting the alt allele. This can be considered similar to Phase Set. Alleles on the same HID are on the same chromosome.">',
           '##FORMAT=<ID=CS,Number=1,Type=String,Description="Case/control status of sample.">']
    meta_info_lines.extend(met1)
    meta_info_lines = [[m] for m in meta_info_lines]
    up_var_head = ["#CHROM", "POS", "ID", "REF", "ALT",
                   "QUAL", "FILTER", "INFO", "FORMAT"]
    for i in gen_index:
        up_var_head.append(var_head[i])
    meta_info_lines.append(up_var_head)
    meta_info_lines.extend(updated_var)
    #save vcf file
    write_list(meta_info_lines, wdir + "variation.vcf")
def convert_to_int(n):
    """
    Convert values to integers. This is to be used when
    a pandas dataframe converts integers to floats due
    to the presence of NA values and integer values are
    preferred over floats, i.e. string conversion/comparison.
    """
    try:
        return int(n)
    except ValueError as e:
        return np.nan
def get_ternary_genotype(gen):
    """
    Convert a 0/0, 0/1, 1/1 type genotype string to
    0, 1, 2.
    """
    try:
        g = sum(map(int, gen.split(":")[0].split("/")))
    except ValueError as e:
        g = np.nan
    return g
def variation_to_geno(settings, var_file, output_prefix):
    """
    Create PLINK files from variation table file.
    """
    wdir = settings["workingDir"]
    case = {}
    with open(wdir + "case_file") as infile:
        for line in infile:
            newline = line.strip().split("\t")
            case[newline[0]] = newline[1]
    with open (wdir + var_file) as infile:
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
                        ped_sample_info.append(["0",
                                                 sam_name,
                                                 "0",
                                                 "0",
                                                 "0",
                                                 affected])
                        geno_sample_info.append(["0",
                                             sam_name,
                                             sam_name,
                                             affected])
                        samples_used.append(s)
                    except KeyError as e:
                        continue
                used_sample_mask = np.array([s in samples_used for s in sample_ids])
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
                except KeyError as e:
                    genes[gene_name] = [rsid]
                    ordered_genes.append(gene_name)
                genotypes = np.array(newline[sample_start_index:])[used_sample_mask]
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
                   wdir + output_prefix + ".ped")
    write_list(map_snp_ids,
                   wdir + output_prefix + ".map")
    header = ["FAMILY_ID", "INDIVIDUAL_ID", "SAMPLE_ID", "AFFECTION"]
    header.extend([s[1] for s in geno_snp_ids])
    geno_sample_info = [header] + geno_sample_info
    write_list(geno_sample_info, wdir + output_prefix + ".geno")
    write_list([["**", o] + genes[o] for o in ordered_genes],
               wdir + output_prefix + ".hlist")
    return
def absence_presence(col, min_val = 1):
    """
    Given a numerical dataframe column, convert to binary values
    for a minimum threshold.
    This should be used by pandas transform or apply.
    """
    return pd.Series([0 if (c < min_val or np.isnan(c))
                      else 1 for c in col.tolist()])
def plot_performance(barcode_counts,
                    tick_label_size = 8,
                    cbar_label_size = 5,
                    dpi = 300,
                     barcode_threshold = 1,
                    absent_color = "black",
                    present_color = "green",
                    save = False,
                    wdir = None,
                    ytick_freq = None,
                    xtick_freq = None,
                    xtick_rotation = 90,
                    tick_genes = False,
                    gene_name_index = None):
    """
    Plot presence/absence plot for a mip run.
    """
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
        plt.xticks(gene_locs, xlabs,
                   rotation = xtick_rotation,
                  ha = "right")
    for ticklabel in ax.get_xticklabels():
        ticklabel.set_fontsize(tick_label_size)
    for ticklabel in ax.get_yticklabels():
        ticklabel.set_fontsize(tick_label_size)
    ax.set_ylabel("Samples")
    ax.set_xlabel("Probes")
    fig.suptitle("Performance",
                verticalalignment="bottom")
    fig.tight_layout()
    cbar = fig.colorbar(heat, ticks = [0, 1],
                            shrink = 0.2
                       )
    cbar.ax.tick_params(labelsize=cbar_label_size)
    cbar.ax.set_yticklabels(["Absent", "Present"])
    fig.set_dpi(dpi)
    fig.tight_layout()
    if save:
        fig.savefig(wdir + "performance.png",
                    dpi = dpi,
                   bbox_inches='tight')
        plt.close("all")
    else:
        return fig,ax
    return
def plot_coverage(barcode_counts,
                    tick_label_size = 8,
                    cbar_label_size = 5,
                    dpi = 300,
                    log = None,
                    save = False,
                    wdir = None,
                    ytick_freq = None,
                    xtick_freq = None,
                    xtick_rotation = 90,
                    tick_genes = False,
                    gene_name_index = None):
    """
    Plot presence/absence plot for a mip run.
    """
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
        cbar_title = ""
    elif log == 2:
        cbar_title = "log2"
        heat = ax.pcolormesh(barcode_counts.transform(
            lambda a: np.log2(a + 1)
        ))
    elif log == 10:
        cbar_title = "log10"
        heat = ax.pcolormesh(barcode_counts.transform(
            lambda a: np.log10(a + 1)
        ))
    elif log == "ln":
        cbar_title = "log"
        heat = ax.pcolormesh(barcode_counts.transform(
            lambda a: np.log(a + 1)
        ))
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
                   rotation = xtick_rotation,
                  ha = "right")
    for ticklabel in ax.get_xticklabels():
        ticklabel.set_fontsize(tick_label_size)
    for ticklabel in ax.get_yticklabels():
        ticklabel.set_fontsize(tick_label_size)
    ax.set_ylabel("Samples")
    ax.set_xlabel("Probes")
    fig.suptitle("Coverage",
                verticalalignment="bottom")
    fig.tight_layout()
    cbar = fig.colorbar(heat,
                            shrink = 0.5
                       )
    cbar.ax.tick_params(labelsize=cbar_label_size)
    cbar.ax.set_ylabel(cbar_title,
                       fontsize = cbar_label_size,
                      rotation = 90)
    fig.set_dpi(dpi)
    fig.tight_layout()
    if save:
        fig.savefig(wdir + "coverage.png",
                    dpi = dpi,
                   bbox_inches='tight')
        plt.close("all")
    else:
        return fig,ax
    return
def split_aa(aa):
    try:
        return aa.split(";")[0].split(":")[4][2:]
    except IndexError as e:
        return "."
    except AttributeError as e:
        return np.nan
def split_aa_pos(aa):
    try:
        return aa.split(";")[0].split(":")[4][2:-1]
    except IndexError as e:
        return "."
    except AttributeError as e:
        return np.nan
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
    return sorted_counts.loc[sid, idx[:,:,:, :, chrom,
                                      :start,
                                      end:]].sum()
def add_known(group, used_targets):
    group = group.merge(used_targets,
                       how = "outer")
    group["Sample ID"].fillna(method = "ffill", inplace = True)
    group["Sample ID"].fillna(method = "bfill", inplace = True)
    group["Chrom"].fillna(group["CHROM"], inplace = True)
    group["Start"].fillna(group["POS"], inplace = True)
    group["End"].fillna(group["POS"], inplace = True)
    group["CHROM"].fillna(group["Chrom"], inplace = True)
    group["POS"].fillna(group["Start"], inplace = True)
    group["Barcode Count"].fillna(0, inplace = True)
    return group
def find_ref_total(group):
    nr = group.loc[~group["ExonicFunc"].isin(["synonymous SNV",
                                             "."]),
                               "Barcode Count"].sum()
    cov = group["POS Coverage"].max()
    return  cov - nr
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
def genotype(settings,
             vt = None,
             known_targets = None,
            output_file = None):
    """
    Call genotypes, calculate genotype stats.
    """
    wdir = settings["workingDir"]
    # load dataframe containing all counts and mutations
    if vt is None:
        vt = pd.read_csv(wdir + "combined_info.csv")
    # remove reference haplotypes
    vt = vt.loc[~vt["VKEY"].isnull()]
    # split mutation names and positions
    vt["AA Change"] = vt["AAChangeClean"].apply(split_aa)
    vt["AA Change Position"] = vt["AAChangeClean"].apply(split_aa_pos)
    vt.drop(["AAChangeClean"],
       inplace = True,
       axis = 1)
    # combine same variation from different haplotypes within a sample
    vt = pd.DataFrame(vt.groupby(([ 'Sample ID', 'Gene',  'VKEY', 'CHROM',
                      'POS', 'ID',  'REF', 'ALT',
                  'Original Position',  'Multi Mapping',
                     'RefGene', 'ExonicFunc', 'Func.refGene',
                  'GeneDetail.refGene', 'Alt',
                  'AA Change', 'AA Change Position'])).agg({
        "Barcode Count": np.sum,
        "Variation Quality": np.max
    })
                 ).reset_index()

    # get the complete set of variants in the data
    var = vt.groupby([ 'Gene',  'VKEY', 'CHROM',
                      'POS', 'ID',  'REF', 'ALT',
                      'Original Position',  'Multi Mapping',
                      'RefGene', 'ExonicFunc', 'Func.refGene',
                      'GeneDetail.refGene', 'Alt',
                      'AA Change', 'AA Change Position'],
                     as_index = False).first().drop(
        ["Variation Quality",
         "Sample ID",
         "Barcode Count"],
         axis = 1
    )
    # load barcode counts and meta data
    processed_data = load_processed_data(settings)
    barcode_counts = processed_data["Barcode Counts"]
    merged_meta = processed_data["Meta Data"]
    # add known variant information
    if known_targets is not None:
        targets = pd.read_csv(known_targets,
                         sep = "\t")
        targets["known"] = "known"
        used_targets = targets.loc[targets["Gene"].isin(
            barcode_counts.columns.levels[
                barcode_counts.columns.names.index("Gene")
                ].tolist())]
        var = var.merge(used_targets, how = "outer")
        var["Chrom"].fillna(var["CHROM"], inplace = True)
        var["CHROM"].fillna(var["Chrom"], inplace = True)
        var["Original Position"].fillna(var["Start"], inplace = True)
        var["POS"].fillna(var["Start"], inplace = True)
        var["ID"].fillna(".", inplace = True)
        var["VKEY"].fillna(var["CHROM"] + ":.:" + var["POS"].astype(str)
                  + ":" + var["Codon"] + ":" + var["Codon"], inplace = True)
    else:
        var["Chrom"] = var["CHROM"]
        var["Start"] = var["POS"]
    try:
        var["Mutation Name"].fillna(var["AA Change"],
                                     inplace = True)
    except KeyError as e:
        var["Mutation Name"] = var["AA Change"]
    var.loc[var["AA Change"] == ".",
            "Mutation Name"] = var.loc[
                        var["AA Change"] == "."].apply(
                                                rename_noncoding,
                                                axis = 1)
    var["Mutation"] = var[["Gene", "Mutation Name"]].apply(
        lambda a: "-".join(a), axis = 1)
    # merge variant information with count information
    # add sample ids to var so that every sample has a row for every variation
    # even if the sample does not have that variation
    var["Temp Key"] = "Temp"
    merged_meta["Temp Key"] = "Temp"
    var_sam = merged_meta[["Sample ID", "Temp Key"]].drop_duplicates().merge(
        var, how = "outer").drop("Temp Key", axis = 1)
    # deletions have no variation quality
    # combined_df["Variation Quality"].fillna(99, inplace = True)
    # Calculate per base coverage from barcode counts
    mutation_counts = var_sam.merge(vt, how = "outer")
    mutation_counts["Coverage Key"] = mutation_counts.apply(
        lambda a: (a["Sample ID"], a["Chrom"], a["Original Position"]),
        axis = 1
    )
    with open(wdir + "coverage.pkld") as infile:
        coverage = pickle.load(infile)
    mutation_counts["Coverage"] = mutation_counts["Coverage Key"].map(coverage)
    all_mutations = mutation_counts
    all_mutations["Barcode Count"].fillna(0, inplace = True)
    all_mutations["Barcode Fraction"] = (all_mutations["Barcode Count"]
                                     /all_mutations["Coverage"])
    #all_mutations = all_mutations.merge(merged_meta)
    min_snp_qual = int(settings["minSnpQuality"])
    min_snp_frac = float(settings["minSnpBarcodeFraction"])
    qual_filter_mask = ((all_mutations["Variation Quality"] >= min_snp_qual)
                     &(all_mutations["Barcode Fraction"] >= min_snp_frac))
    all_mutations.loc[qual_filter_mask, "Filtered Barcode Count"] = (
        all_mutations.loc[qual_filter_mask, "Barcode Count"]
    )
    all_mutations["Filtered Barcode Count"].fillna(0, inplace = True)
    all_mutations["Mutation Count"] = all_mutations.groupby(
        ["Sample ID", "Mutation"])["Filtered Barcode Count"].transform(sum)
    all_mutations["Filtered Barcode Fraction"] = (
        all_mutations["Filtered Barcode Count"]
        /all_mutations["Coverage"])
    # sum mutation counts per position, then determine reference
    # nucleotide counts for that position
    all_mutations = all_mutations.merge(
        pd.DataFrame(all_mutations.groupby(["Sample ID",
                                            "CHROM",
                                           "POS"])
                    .agg({"Barcode Count": np.sum,
                         "Coverage": np.max})).rename(columns = {
            "Barcode Count": "Non-reference Base Count",
            "Coverage": "POS Coverage"
        }).reset_index()
    )
    all_mutations["Reference Base Count"] = (
        all_mutations["POS Coverage"]
        - all_mutations["Non-reference Base Count"]
    )
    # group all mutations per sample per amino acid position
    # find reference amino acid counts by checking the exonicFunc
    # field. If this is "." (used for noncoding changes)
    # or nonsynonymous SNV, it is considered reference AA.
    all_mutations = all_mutations.merge(
        pd.DataFrame(all_mutations.groupby(["Sample ID",
                                            "Gene",
                                           "AA Change Position"])
                    .apply(find_ref_total)).rename(columns = {
            0: "Reference AA Count"
        }).reset_index()
    )
    all_mutations["Non-reference AA Count"] = (all_mutations["POS Coverage"]
                                            - all_mutations["Reference AA Count"])
    all_mutations["Wildtype AA Count"] = all_mutations["Reference AA Count"]
    try:
        all_mutations["Reference Resistant"]
    except KeyError as e:
        all_mutations["Reference Resistant"] = "No"
    resistant_ref_mask = all_mutations["Reference Resistant"] == "Yes"
    resistant_ref_mask.fillna(False, inplace = True)
    all_mutations.loc[resistant_ref_mask, "Wildtype AA Count"] = (
        all_mutations.loc[resistant_ref_mask, "Non-reference AA Count"])
    all_mutations.loc[
        resistant_ref_mask, "Mutation Count"
    ] = all_mutations.loc[resistant_ref_mask, "Reference AA Count"]
    all_mutations.loc[
        resistant_ref_mask, "Wildtype AA Count"
    ] = all_mutations.loc[resistant_ref_mask, "Non-reference AA Count"]
    all_mutations["Mutation Fraction"] = (all_mutations["Mutation Count"]
                                      /all_mutations["Coverage"])
    min_frac_cutoff = float(settings["minMutationFraction"])
    all_mutations["Genotype"] = all_mutations.apply(get_genotype,
                                                    axis = 1,
                                                   cutoff = min_frac_cutoff)
    all_mutations["Aminotype"] = all_mutations.apply(get_aminotype,
                                                     axis = 1,
                                                   cutoff = min_frac_cutoff)
    all_mutations["Binary Genotype Call"] = all_mutations["Genotype"]
    all_mutations.loc[all_mutations["Genotype"] == "MIX",
                 "Binary Genotype Call"] = "MUT"
    all_mutations["Binary Aminotype Call"] = all_mutations["Aminotype"]
    all_mutations.loc[all_mutations["Aminotype"] == "MIX",
                 "Binary Aminotype Call"] = "MUT"
    if output_file is None:
        output_file = "mutation_summary.tsv"
    all_mutations.to_csv(wdir + output_file,
                        sep = "\t",
                        index = False)
    return
def call_microsats(settings, sim = None, freq_cutoff = 0.005,
                   min_bc_cutoff = 0,
                  use_filtered_mips = True,
                  ref_genome = "Pf3d7"):
    wdir = settings["workingDir"]
    with open(wdir + settings["perSampleResults"]) as infile:
        sample_results = json.load(infile)
    with open(wdir + settings["haplotypeDictionary"]) as infile:
        hap_dict = json.load(infile)
    if sim is None:
        sim = pd.read_csv("resources/pf_MS/simulation.tsv", sep = "\t")
    ref_sim = sim.loc[sim["genome"] == ref_genome]
    strain_freqs = {}
    for sample_name in sample_results:
        sam_res = sample_results[sample_name]
        sam_freq = {}
        for g in sam_res:
            for m in sam_res[g]:
                if use_filtered_mips and (ref_sim.loc[
                    ref_sim["MIP"] == m,
                    "Clean MS MIP"].values[0] == False):
                    continue
                for c in sam_res[g][m]:
                    total_bcs = float(sam_res[g][m][c]["cumulative_data"]["barcode_count"])
                    if total_bcs >= min_bc_cutoff:
                        filtered_data = sam_res[g][m][c]["filtered_data"]
                        ms_types = {}
                        for hd in filtered_data:
                            bcc = hd["barcode_count"]
                            h = hd["haplotype_ID"]
                            h_len = len(hap_dict[m][h]["sequence"])
                            ms_len = int(h_len - ref_sim.loc[ref_sim["MIP"] == m,
                                                         "MS size adjustment"])
                            try:
                                ms_types[ms_len] += bcc
                            except KeyError as e:
                                ms_types[ms_len] = bcc
                        for ml in list(ms_types.keys()):
                            if (ms_types[ml]/total_bcs) < freq_cutoff:
                                ms_types.pop(ml)
                        try:
                            sam_freq[g][m][c] = ms_types
                        except KeyError as e:
                            try:
                                sam_freq[g][m] = {c : ms_types}
                            except KeyError as e:
                                sam_freq[g] = {m : {c : ms_types}}
        strain_freqs[sample_name] = sam_freq
    ms_calls = []
    for s in strain_freqs:
        for g in strain_freqs[s]:
            for m in strain_freqs[s][g]:
                for c in strain_freqs[s][g][m]:
                    for l in strain_freqs[s][g][m][c]:
                        ms_calls.append([s, g, m, c, l, g + "-" + str(int(l)),
                                        strain_freqs[s][g][m][c][l]])
    ms_call_df = pd.DataFrame(ms_calls,
                columns = ["Sample ID",
                          "region",
                          "MIP",
                          "Copy",
                          "length",
                           "haplotype name",
                          "count"]).drop("Copy",
                                        axis = 1)
    merged_calls = pd.DataFrame(ms_call_df.groupby(["Sample ID",
                    "region",
                    "haplotype name",
                    "length",
                    ])["count"].sum())
    merged_calls.reset_index(inplace = True)
    merged_calls["frequency"] = merged_calls.groupby(["Sample ID",
                      "region"])["count"].transform(lambda a: a/a.sum())
    merged_calls.rename(columns = {"length": "MS Length"},
                   inplace = True)
    merged_calls = merged_calls.merge(sim.groupby(
        ["region", "MS Length", "Unique Strain"],
        as_index = False).first()[["region", "MS Length", "Unique Strain"]],
                                     how = "left")
    return {"ms_calls": merged_calls,
            "strain_freqs": strain_freqs}
def get_copy_counts(count_table,
                   average_copy_count = 2,
                   norm_percentiles = [0.4, 0.6]):
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
        lambda a: a/a.sum(), axis = 1)
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
        return ac.loc[r["Sample ID"],
                       (r["Gene"],
                       r["Copy"])]
    except KeyError as e:
        return np.nan
def normalize_copies(a):
    if a.isnull().all():
        a = a.fillna(1)
        return a/a.sum()
    else:
        return a.fillna(0)
def repool(
    wdir,
    data_summary,
    high_barcode_threshold,
    target_coverage_count = None,
    target_coverage_fraction = 0.95,
    target_coverage_key = "targets_with_10_barcodes",
    barcode_coverage_threshold = 10,
    barcode_count_threshold = 100,
    low_coverage_action = "Repool",
    assesment_key = "targets_with_1_barcodes",
    good_coverage_quantile = 0.25,
    output_file = "repool.csv"
):
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
    except KeyError as e:
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
        &(data_summary["total_barcode_count"] < barcode_count_threshold),
        "Status"
    ] = low_coverage_action
    # mark samples with too high barcode coverage
    # these samples will have been sequenced to a high depth but
    # low barcode numbers, so sequencing these more would not make sense.
    # They will be re-captured if more data is needed.
    try:
        data_summary["Barcode Coverage"]
    except KeyError as e:
        data_summary["Barcode Coverage"] = (
            data_summary["total_read_count"]
            /data_summary["total_barcode_count"]
        ).fillna(0)
    data_summary.loc[
        (data_summary["Status"].isnull())
        &(data_summary["Barcode Coverage"] >= barcode_coverage_threshold),
        "Status"
    ] = "Recapture"
    # Zero barcode coverage is presumably due to poor sequencing
    # So low coverage action should be taken.
    data_summary.loc[
        (data_summary["Status"].isnull())
        &(data_summary["Barcode Coverage"] == 0),
        "Status"
    ] = low_coverage_action
    # All remaining samples will be repooled
    data_summary.loc[
        (data_summary["Status"].isnull()),
        "Status"
    ] = "Repool"
    data_summary["Library to Completion"] = ((high_barcode_threshold -
                                              data_summary["total_barcode_count"])
                                             /data_summary["total_barcode_count"])
    # replace inf values with max
    lc_max = data_summary.loc[
        data_summary["Library to Completion"] < np.inf,
        "Library to Completion"
    ].max()
    data_summary.loc[
        data_summary["Library to Completion"] == np.inf,
        "Library to Completion"
    ] = lc_max
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
        /data_summary.loc[data_summary[target_coverage_key] > 0,
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
        (data_summary["Barcodes Per Target Covered"]
          > good_coverage_threshold)
        &(data_summary[assesment_key] < target_coverage_count),
        "Uneven Coverage"
    ] = True
    data_summary.loc[data_summary["Uneven Coverage"].isnull(),
                     "Uneven Coverage"] = False
    try:
        data_summary.to_csv(wdir + output_file,
                           index = False)
    except TypeError as e:
        # in an older version of this function, settings dict
        # was passed instead of wdir, for backwards compatibility
        # we'll catch that error and use wdir from the settings dict
        data_summary.to_csv(wdir["workingDir"] + output_file,
                           index = False)
    print(("Out of %d samples %d are completed, %d will be recaptured and %d repooled" %(
        data_summary.shape[0],
        data_summary.loc[data_summary["Status"] == "Complete"].shape[0],
        data_summary.loc[data_summary["Status"] == "Recapture"].shape[0],
        data_summary.loc[data_summary["Status"] == "Repool"].shape[0])))
    print(("%d samples showed uneven coverage, %d complete, %d to be recaptured, %d repooled"%(
        data_summary.loc[data_summary["Uneven Coverage"]].shape[0],
        data_summary.loc[data_summary["Uneven Coverage"]
                          & (data_summary["Status"] == "Complete")].shape[0],
        data_summary.loc[data_summary["Uneven Coverage"]
                          & (data_summary["Status"] == "Recapture")].shape[0],
        data_summary.loc[data_summary["Uneven Coverage"]
                          & (data_summary["Status"] == "Repool")].shape[0])))
    return
def aa_to_coordinate(gene, species, aa_number, alias = False):
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
        except KeyError as e:
            pass
    cds = get_cds(gene, species)
    if len(cds) == 0:
        return [np.nan, np.nan, np.nan,
                np.nan, np.nan, np.nan]
    ori = cds["orientation"]
    coord = cds["coordinates"]
    chrom = cds["chrom"]
    if ori == "+":
        aa_end = aa_number*3 - 1
        aa_start = aa_number*3 - 3
    else:
        aa_start = aa_number*3 - 1
        aa_end = aa_number*3 - 3
    cds_start = coord[aa_start]
    cds_end = coord[aa_end]
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
    unique_haplotype_file = wdir + settings["haplotypeDictionary"]
    with open(
        unique_haplotype_file) as infile, open(
        unique_haplotype_file + id_generator(6), "w") as outfile:
        haplotypes = json.load(infile)
        json.dump(haplotypes, outfile)
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
                    diffs = haplotypes[m][h]["mapped_copies"][cp]["differences"]
                    aa_changes = {}
                    multi_indels = []
                    for i in range(len(diffs)):
                        d = diffs[i]
                        # get protein change information for the SNP
                        ano = d["annotation"]["AAChangeClean"]
                        # this would look like so
                        # 'mal_mito_3:mal_mito_3:exon1:c.G673A:p.V225I'
                        try:
                            aa = ano.split(":")
                            # get aa change position, e.g. 225 for V225I
                            aa_pos = int(aa[-1].split(".")[-1][1:-1])
                            # this will generate an IndexError or
                            # ValueError when the change is not a SNP
                            # that is a single aa change (such as indels,
                            # or noncoding changes). Those should not be
                            # merged.
                            try:
                                # add the aa change position to the changes dict.
                                aa_changes[aa_pos].append(i)
                            except KeyError as e:
                                aa_changes[aa_pos] = [i]
                        except (IndexError, ValueError) as e:
                            continue
                    # after going through all diffs, look for mutliple diffs
                    # affecting single aminoacid
                    all_merges = []
                    for c in aa_changes:
                        if (len(aa_changes[c]) > 1) and (c not in multi_indels):
                            # break out of loop if indels found
                            mindel = False
                            # merge multiple snps affecting the same aa
                            merge_dict = {}
                            indexes = aa_changes[c]
                            merge = merge_dict[c] = {}
                            # keep positions relative to cDNA in a list
                            c_positions = []
                            # keep positions relative to haplotype in a list
                            h_indexes = []
                            # keep genomic positions of changes in a list
                            g_positions = []
                            # keep the difference between the cDNA and genomic
                            # positions in a list. This will be used to determine
                            # the gene's orientation on the genome.
                            c_offsets = []
                            changes_to_cdna = []
                            changes_to_aa = []
                            for i in indexes:
                                d = diffs[i]
                                # for each diff get the annotation
                                # e.g. 'mal_mito_3:mal_mito_3:exon1:c.G673A:p.V225I'
                                ano = d["annotation"]["AAChangeClean"]
                                aa = ano.split(":")[-1].split(".")[-1]
                                changes_to_aa.append(aa)
                                # get the aa of reference genome (V)
                                aa_ref = aa[0]
                                # get cdna change, e.g. G673A
                                cdna = ano.split(":")[-2].split(".")[-1]
                                changes_to_cdna.append(cdna)
                                # get the mutant base (A)
                                cdna_change = cdna[-1]
                                # compare the sequence on the cDNA with
                                # to the sequence of the haplotype,
                                # to determine the MIP/haplotype's orientation
                                # relative to the cDNA
                                ori = cdna_change == d["hap_base"]
                                try:
                                    cdna_pos = int(cdna[1:-1])
                                    # raises value error if not a SNP
                                except ValueError as e:
                                    multi_indels.append(c)
                                    mindel = True
                                    break
                                # get genomic position of the change
                                diff_start = int(d["annotation"]["Start"])
                                # get the difference between the genomic and
                                # cDNA position of the change, to be used
                                # in determining the gene's orientation
                                pos_diff = diff_start - cdna_pos
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
                            # the gene is on the plus strand, else it is reverse.
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
                                # if the haplotype is in the opposite orientation
                                # as the cDNA
                                h_end_index = h_indexes[-1] + codon_offset + 1
                                h_start_index = h_end_index - 3
                                hap_codon = hap_seq[h_start_index:h_end_index]
                                codon = reverse_complement(hap_codon)
                            # get the genomic position and sequence of the codon
                            if gene_ori:
                                g_start = g_positions[0] - codon_offset
                                g_end = g_start + 2
                            else:
                                g_end = g_positions[-1] + codon_offset
                                g_start = g_end -2
                            # extract the reference codon sequence
                            Ref = get_sequence(create_region(d["chrom"],
                                                                        g_start,
                                                                        g_end), species)
                            if gene_ori:
                                g_codon = Ref
                                Alt = codon
                            else:
                                g_codon = reverse_complement(Ref)
                                Alt = reverse_complement(codon)
                            # calculate merged codon's amino acid
                            merged_aa = translate(codon)
                            # recreate the annotation string for the merge
                            aa_change_base = ano.split(":")[:-2]
                            protein_change = "p." + aa_ref + str(c) + merged_aa
                            coding_change = "c." + g_codon + str(codon_pos) + codon
                            aa_change_base.extend([coding_change, protein_change])
                            AAChange = ":".join(aa_change_base)
                            # determine if the merged change is synonymous
                            if aa_ref == merged_aa:
                                ExonicFunc = "synonymous SNV"
                            else:
                                ExonicFunc = "nonsynonymous SNV"
                            merged_dict = {'annotation': {
                                'AAChangeClean': AAChange,
                                'Alt': Alt,
                                'Chr': d["chrom"],
                                'End': g_end,
                                'ExonicFunc.refGene': ExonicFunc,
                                'Func.refGene': 'exonic',
                                'Gene.refGene': d["annotation"]['Gene.refGene'],
                                'GeneDetail.refGene': d["annotation"]['GeneDetail.refGene'],
                                'Otherinfo': d["annotation"]["Otherinfo"],
                                'Ref': Ref,
                                'Start': g_start
                            },
                                           'base_match': False,
                                            'begin': g_start,
                                            'chrom': d["chrom"],
                                            'clinical': False,
                                            'clinical_id': 'none',
                                            'end': g_end,
                                            'hap_base': hap_codon,
                                            'hap_index': [h_start_index, h_end_index - 1],
                                            'psv': False,
                                            'ref_base': Ref,
                                            'type': 'snp',
                                            'vcf_normalized': ":".join(
                                              [d["chrom"], str(g_start), ".", Ref, Alt]
                                            ),
                                            'vcf_raw': ":".join(
                                              [d["chrom"], str(g_start), ".", Ref, Alt]
                                            ),
                                           "gene_ori": gene_ori,
                                           "ori": ori
                                          }
                            all_merges.append(merged_dict)
                            outlist.append([h, c, merged_aa, aa_ref,
                                           ",".join(changes_to_cdna),
                                            codon, g_codon,
                                           ",".join(changes_to_aa)])
                    # Remove SNPs that were merged, add the merged SNP
                    for c in aa_changes:
                        if (len(aa_changes[c]) > 1) and (c not in multi_indels):
                            indexes = aa_changes[c]
                            for i in indexes:
                                diffs[i] = "remove"
                    diffs.extend(all_merges)
                    diffs = [d for d in diffs if d != "remove"]
                    haplotypes[m][h]["mapped_copies"][cp]["differences"] = diffs
    with open(unique_haplotype_file, "w") as outfile:
        json.dump(haplotypes, outfile, indent = 1)
    # save the report
    write_list(outlist, wdir + "merge_snps_output.txt")
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
    bc = pd.read_csv(wdir + "barcode_counts.csv",
                header = [0,1,2,3,4,5,6],
                index_col = 0)
    bc.T.to_csv(wdir + "barcode_counts.T.csv")
    bc = pd.read_csv(wdir + "barcode_counts.T.csv",
                    index_col = [0,1,2,3,4,5,6]).T
    data_summary = pd.read_csv(wdir + "data_summary.csv",
                          index_col = None)
    merged_meta = pd.read_csv(wdir + "meta_data.csv",
                         index_col = None)
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
    except IOError as e:
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
                    except IndexError as e:
                        info_dict[split_info[0]] = [True]
                # sum of all allele counts will be used as allele count
                ac = info_dict["AC"]
                an = info_dict["AN"][0]
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
                            except KeyError as e:
                                var_info.append(np.nan)
                            except IndexError as e:
                                var_info.append(info_dict[col][0])
                        outlist = outlist + var_info
                        outfile_list.append(outlist)
                # if not a SNP, must be indel
                # indels must have their own line, hence only 1 indel in alt bases
                elif len(alt_bases) > 1:
                    problem_alts.append(newline)
                    break
                # if conforming indel:
                else:
                    alt_base = alt_bases[0]
                    # vcf files have the  indels together with the preceding base
                    # such as REF: TA, ALT: T
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
                            for i in range (2):
                                outlist = [chrom, pos + i, var_id,
                                           "-", alt_bases,
                                          qual, filt]
                                var_info = []
                                for col in info_cols:
                                    try:
                                        var_info.append(info_dict[col][0])
                                    except KeyError as e:
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
                                           ref_base[i], "-",
                                          qual, filt]
                                var_info = []
                                for col in info_cols:
                                    try:
                                        var_info.append(info_dict[col][0])
                                    except KeyError as e:
                                        var_info.append(np.nan)
                                outlist = outlist + var_info
                                outfile_list.append(outlist)
        var_df = pd.DataFrame(outfile_list, columns = ["CHROM",
                                              "POS",
                                              "ID",
                                              "REF",
                                              "ALT",
                                              "QUAL",
                                              "FILTER"]
                     + info_cols)
    var_df = var_df.astype({"AN": int, "AC": int})
    if len(problem_alts) > 0:
        print(("There are %d problematic alleles, see the output list for details"
               %len(problem_alts)))
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
def vcf_to_table(collapsed, output_file):
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
        ref_count = row["AN"] - row ["AC"]
        return str(ref_count) + "," + str(row["AC"]) + ","
    collapsed["AS"] = collapsed.apply(get_allele_strings, axis = 1)
    collapsed["CS"] = collapsed.apply(get_count_strings, axis = 1)
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
    table.to_csv(output_file,
                 sep = "\t",
                 index = False,
                 header = False)
    subprocess.call(["bgzip", "-c", output_file],
                         stdout = open(output_file + ".gz", "w"))
    subprocess.call(["tabix", "-0", "-s 2",
                         "-b 3", "-e 3", output_file + ".gz"])
    return table
def header_to_primer (seq_to_bc_dict,
                     header_string,
                     platform):
    """
    Convert a demultiplexed fastq header to forward and
    reverse primer numbers.
    """
    split_string = header_string.split("+")
    if platform == "miseq":
        try:
            fw = seq_to_bc_dict[split_string[0]]
        except KeyError as e:
            fw = 999
        try:
            rev = seq_to_bc_dict[
                reverse_complement(split_string[1])
        ]
        except KeyError as e:
            rev = 999
    elif platform == "nextseq":
        try:
            fw = seq_to_bc_dict[
                reverse_complement(split_string[1])
            ]
        except KeyError as e:
            fw = 999
        try:
            rev = seq_to_bc_dict[
                reverse_complement(split_string[0])
            ]
        except KeyError as e:
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
        return  fw_seq + "+" + reverse_complement(rev_seq)
def check_stitching(stitch_file):
    """
    Take a stitch log file from MIPWrangler output, return summary datframe.
    """
    with open(stitch_file) as infile:
        stitch = []
        for line in infile:
            newline = line.strip()
            stitch.append(line)
    sti_sum = []
    for l in stitch:
        if '\t"stdOut_" : "[FLASH] Starting FLASH v' in l:
            nl = l.split("\\n")
            #print nl
            #break
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
    sti = pd.DataFrame(sti,
                columns = ["Sample ID",
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
    except IOError as e:
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
    """ Return a list of all possible bases corresponding
    to a given iupac nucleotide code. """
    iupac_dict = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT"
    }
    try:
        return list(iupac_dict[iupac_code.upper()])
    except KeyError as e:
        print((("Non-IUPAC nucleotide code {}."
               " Code must be one of {}").format(
            iupac_code,
            "".join(list(iupac_dict.keys()))
        )))
        return []
def iupac_fasta_converter(header, sequence):
    """
    Given a sequence (header and sequence itself)
    containing iupac characters, return a dictionary with
    all possible sequences converted to ATCG.
    """
    iupac_dict = {
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT"
    }
    iupac_dict = {k: list(iupac_dict[k])
                  for k in list(iupac_dict.keys())}
    if sequence.upper().count("N") >= 10:
        return {header : sequence}
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
    if len (result_list) == 1:
        return {header: result_list[0]}
    else:
        return {header + "-" + str(i) : result_list[i]
                for i in range(len(result_list))}
def save_fasta_dict(fasta_dict, fasta_file, linewidth = 60):
    """ Save a fasta dictionary to file. """
    with open(fasta_file, "w") as outfile:
        for header in fasta_dict:
            outfile.write(">" + header + "\n")
            fasta_seq = fasta_dict[header]
            for i in range(0, len(fasta_seq), linewidth):
                outfile.write(fasta_seq[i: i + linewidth] + "\n")

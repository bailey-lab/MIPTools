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
from sklearn.cluster import MeanShift, DBSCAN
import matplotlib.pyplot as plt
from matplotlib import colors
from sklearn.manifold import TSNE
from scipy.stats import chi2_contingency, fisher_exact
import pysam
import mip_classes as mod
import pandas as pd
import gzip
from primer3 import calcHeterodimerTm
import traceback
from msa_to_vcf import msa_to_vcf as msa_to_vcf

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


def get_file_locations():
    """ All static files such as fasta genomes, snp files, etc. must be listed
    in a file in the working directory. File name is file_locations.
    It is a tab separated text file. First tab has 2 letter species name, or
    "all" for general files used for all species. Second tab is the file name
    and third is the location of the file, either relative to script working
    directory, or the absolute path."""
    file_locations = {}
    with open("/opt/resources/file_locations", "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                if newline[0] not in list(file_locations.keys()):
                    file_locations[newline[0]] = {newline[1]: newline[2]}
                else:
                    file_locations[newline[0]][newline[1]] = newline[2]
    return file_locations


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
            except KeyError:
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


def get_target_coordinates(res_dir, species, capture_size,
                           coordinates_file=None, snps_file=None,
                           genes_file=None, capture_types={}):
    """ Extract MIP target coordinates from provided files. """
    # Get target coordinates specified as genomic coordinates
    if coordinates_file is not None:
        coordinates_file = os.path.join(res_dir, coordinates_file)
        try:
            coord_df = pd.read_table(coordinates_file, index_col=False)
            coord_df.rename(columns={"Name": "name", "Chrom": "chrom",
                            "Start": "begin", "End": "end"}, inplace=True)
            region_coordinates = coord_df.set_index("name").to_dict(
                orient="index")
            # update capture types of targets
            for g in region_coordinates:
                if g not in capture_types:
                    capture_types[g] = region_coordinates[g]["Capture Type"]
        except IOError:
            print(("Target coordinates file {} could not be found.").format(
                (coordinates_file)))
            region_coordinates = {}

    # Get Gene target coordinates
    if genes_file is not None:
        # get the alias file (gene name to gene id mapping) if available
        try:
            with open(get_file_locations()[species]["alias"]) as infile:
                alias = json.load(infile)
        except (KeyError, IOError):
            pass
        try:
            genes_file = os.path.join(res_dir, genes_file)
            genes_df = pd.read_table(genes_file, index_col=False)
            genes = genes_df.set_index("Gene").to_dict(orient="index")
            gene_names = list(genes.keys())
            gene_id_to_gene = {}
            gene_ids = []
            gene_coordinates = {}
            for g in genes:
                try:
                    if np.isnan(genes[g]["Gene ID"]):
                        try:
                            gene_id = alias[g]
                            genes[g]["Gene ID"] = gene_id
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
                gene_ids.append(genes[g]["Gene ID"])
                gene_id_to_gene[genes[g]["Gene ID"]] = g
                capture_types[g] = genes[g]["Capture Type"]
            gene_id_coordinates = gene_to_target(gene_ids, species)
            for gid in gene_id_coordinates:
                gene_coordinates[gene_id_to_gene[gid]] = gene_id_coordinates[
                    gid]
        except IOError:
            print(("Target genes file {} could not be found.").format(
                (genes_file)))
            gene_coordinates = {}
            gene_names = []

    # Get SNP target coordinates
    try:
        snps_file = os.path.join(res_dir, snps_file)
        snp_df = pd.read_table(snps_file, index_col=False)
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
    for c in all_coordinates.keys():
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
    # prioritize gene names over snp or other names
    for t in list(target_names.keys()):
        for n in target_names[t]:
            if n in gene_names:
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
    """ Merge overlapping coordinates for MIP targets.

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
            with open(res_dir + t + ".fa", "w") as outfile:
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
    for f in fasta_sequences.keys():
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
        with open(res_dir + newf + ".fa", "w") as outfile:
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
    """ Parallelize a list of lastz alignments."""
    p = Pool(pro)
    p.map_async(align_region_worker, alignment_list)
    p.close()
    p.join()
    return


def align_region_worker(l):
    """ Worker function for align_region_multi.
    Aligns a single fasta file to a target fasta file.
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
    query_fasta = resource_dir + region_key
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
            "--output=" + resource_dir + output_file,
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
    """ Align sequences given in an alignment dict which contains alignment
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
    with open(resource_dir + output_prefix + ".al", "w") as alignment_file:
        for f in fasta_list:
            fnum = 0
            with open(resource_dir + f + ".al") as alignment:
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
    with open(resource_dir + output_prefix + ".differences", "w") as diff_file:
        for f in fasta_list:
            fnum = 0
            with open(resource_dir + f + ".differences") as diffs:
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
    with open(wdir + name + ".al") as infile:
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
    while overlap_found:
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
    # clean up overlapping region lists by removing duplicates.
    for o in overlaps:
        overlaps[o] = sorted(list(set(overlaps[o])))
    #########################################
    # create a new dictionary for target regions.
    # for each target group in overlaps, we'll have genomic coordinates
    # that will be used as final targets.
    #########################################
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
    # organize target regions, assign region names based on the original
    # target names. Assign a reference target.
    ###########################################
    # sort target regions based on the length of
    # chromosome name and the length of region. Chromosome name is used
    # to distinguish alternate contigs and not use them as reference, but
    # it is not absolutely necessary and it would not behave as expected
    # when chromosome names do not follow that convention, i.e, chr6 and
    # chr6_altXYZ
    for ar in aligned_regions:
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
    for r in region_names:
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
                                max_allowed_indel_size,
                                match_score=1, mismatch_score=5,
                                gap_open_penalty=20, gap_extend_penalty=5
                                ):
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
        ydrop = max_allowed_indel_size * gap_extend_penalty + gap_open_penalty
        alignment_opts = ["--match=" + str(match_score) + "," + str(
            mismatch_score), "--gap=" + str(gap_open_penalty) + "," + str(
            gap_extend_penalty), "--ydrop=" + str(ydrop), "--notransition",
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
    """ Align all regions within a target group to the region selected
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
        with open(resource_dir + t + ".query.fa", "w") as outfile:
            outfile.write(">" + t + "_ref\n")
            outfile.write(get_sequence(query_key, species))
        # create a fasta file that includes all target regions within a group.
        with open(resource_dir + t + ".targets.fa", "w") as outfile:
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
                resource_dir + t + ".targets.fa",
                ["multiple", "unmask", "nameparse=darkspace"],
                ["unmask", "nameparse=darkspace"],
                identity, coverage, gen_out,
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
    return align_region_multi(alignment_commands, num_process)


def intra_alignment_checker(family_name, res_dir, target_regions,
                            region_names):
    """ Following a within group alignment, check if any individual region
    within the group has multiple aligned parts. If found, split that region
    into multiple regions to be re-aligned by intraparalog_aligner.
    """
    alignment_file = family_name + ".aligned"
    new_regions = {}
    with open(res_dir + alignment_file, "r") as alignment:
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
    ret_regions = []
    rnames = []
    for ci in sorted(new_regions):
        ret_regions.extend(sorted(new_regions[ci]))
        if len(new_regions[ci]) > 1:
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
                if overlap(main_region[1:3], ur):
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
                  min_target_size):
    # create fasta files for each target coordinate
    create_target_fastas(res_dir, target_regions, species, flank)

    # add target sequences provided by fasta files
    fasta_targets = add_fasta_targets(res_dir, fasta_files,
                                      fasta_capture_type=fasta_capture_type)
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
    genome_alignment = alignment_parser(res_dir, "merged", spacer=0,
                                        gene_names=gene_names)
    target_regions = copy.deepcopy(genome_alignment[0])
    region_names = copy.deepcopy(genome_alignment[1])
    imperfect_aligners = genome_alignment[2]
    aligned_regions = genome_alignment[3]
    overlaps = genome_alignment[4]

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

    out_dict = {"original_target_regions": original_target_regions,
                "target_regions": target_regions,
                "region_names": region_names,
                "aligned_regions": aligned_regions,
                "capture_types": capture_types,
                "imperfect_aligners": imperfect_aligners,
                "overlaps": overlaps,
                "missed_target_regions": missed_target_regions,
                "missed_target_names": missed_target_names,
                "missed_capture_types": missed_capture_types}
    return out_dict


def alignment_mapper(family_name, res_dir):
    """ Create a coordinate map of within group alignments.
    """
    alignment_file = family_name + ".aligned"
    difference_file = family_name + ".differences"
    with open(
        res_dir + alignment_file, "r"
    ) as alignment, open(res_dir + difference_file, "r") as difference:
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
        with open(res_dir + o, "w") as outfile:
            outfile_list = ["\t".join(["WellPosition", "Name", "Sequence"])]
            plate_mips = order_dict[o]
            for m in plate_mips:
                wp = m[-4] + str(m[-3])
                outfile_list.append("\t".join([wp, m[0], m[-1]]))
            outfile.write("\n".join(outfile_list))
    return
###############################################################
# Data analysis related functions
###############################################################


def get_analysis_settings(settings_file):
    """ Convert analysis settings file to dictionary"""
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
    """ Create a settings file from a settings dictionary."""
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
    sequence_to_haplotype_file = (wdir
                                  + settings["sequenceToHaplotypeDictionary"])
    call_info_file = settings["callInfoDictionary"]
    species = settings["species"]
    try:
        tol = int(settings["alignmentTolerance"])
    except KeyError:
        tol = 50
    # DATA EXTRACTION ###
    # if there is no previous haplotype information, an empty dict will be used
    # for instead of the known haplotypes dict
    try:
        with open(sequence_to_haplotype_file) as infile:
            sequence_to_haplotype = json.load(infile)
    except IOError:
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
                # add the haplotype to the dict
                # if it has not been mapped before
                if hapseq not in sequence_to_haplotype:
                    haps[hapname] = {"sequence": hapseq,
                                     "quality": hapqual}

    # BWA alignment ####
    # create a fastq file for bwa input
    with open(haplotypes_fq_file, "w") as outfile:
        for h in haps:
            outfile.write("@" + h + "\n")
            outfile.write(haps[h]["sequence"] + "\n" + "+" + "\n")
            outfile.write(haps[h]["quality"] + "\n")
    # re-structure haplotypes dictionary and initialize a hits dictionary for
    # all haplotypes that will hold the bowtie hits for each of the haplotypes
    # keys for this dict will be mipnames
    haplotypes = {}
    for h in haps:
        mip_name = h.split(".")[0]
        try:
            haplotypes[mip_name][h] = {"sequence": haps[h]["sequence"]}
        except KeyError:
            haplotypes[mip_name] = {h: {"sequence": haps[h]["sequence"]}}
    # run bwa
    bwa(haplotypes_fq_file, haplotypes_sam_file, "sam", "", "", bwa_options,
        species)
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
                except IndexError:
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
    # create a dataframe from the call info dictionary
    # to determine minimum and maximum coordinates for each gene
    # this will be used for haplotypes that do not map to their intended
    # targets but still are not offtarget haplotypes because they map
    # to one of the regions of interest. A situation like this could be caused
    # by two regions of sufficient similarity to be captured by a single mip
    # even though regions were not determined to be paralogus. Or in any
    # situation where the reads were assigned to the wrong MIP for any reason.
    call_copy = copy.deepcopy(call_info)
    call_df_list = []
    for g in call_copy:
        for m in call_copy[g]:
            mip_number = int(m.split("_")[-1][3:])
            sub_number = int(m.split("_")[-2][3:])
            for c in call_copy[g][m]["copies"]:
                call_dict = call_copy[g][m]["copies"][c]
                try:
                    call_dict.pop("genes")
                except KeyError:
                    pass
                call_dict["gene"] = g
                call_dict["copy"] = c
                call_dict["mip_number"] = mip_number
                call_dict["sub_number"] = sub_number
                call_df_list.append(pd.DataFrame(call_dict, index=[0]))
    call_df = pd.concat(call_df_list)
    gene_df = call_df.groupby(["gene", "copy"]).agg(
        {"chrom": "first",
         "capture_start": np.min,
         "capture_end": np.max,
         "copyname": "first",
         "mip_number": np.max,
         "sub_number": np.max}
    )
    gene_dict = gene_df.to_dict(orient="index")
    gene_df = gene_df.reset_index()
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
                        if ((copy_chrom == hit_chrom) and
                                (copy_begin - tol < hit_pos < copy_end + tol)):
                            haplotypes[m][h]["mapped"] = True
                            break
        except KeyError:
            for h in list(haplotypes[m].keys()):
                haplotypes[m][h]["mapped"] = False
    # remove haplotypes that mapped best to an untargeted location on genome
    off_target_haplotypes = {}
    secondary_haplotypes = []
    secondary_haplotype_dict = {}
    for m in list(haplotypes.keys()):
        for h in list(haplotypes[m].keys()):
            if not haplotypes[m][h]["mapped"]:
                # check other genes/MIPs for off targets that are
                # on other possible targets due to sequence similarity
                best_hits = haplotypes[m][h]["best_hits"][0]
                secondary_haplotype_found = False
                for record in best_hits:
                    if secondary_haplotype_found:
                        break
                    flag = record[1]
                    # a flag value of 4 means there was no hit,
                    # so pass those records
                    if flag == "4":
                        continue
                    hit_chrom = record[2]
                    hit_pos = int(record[3])
                    # get cigar string of alignment
                    cigar = record[5]
                    # extract which strand is the bowtie hit on
                    # true if forward
                    strand = ((int(record[1]) % 256) == 0)
                    # bowtie gives us the start position of the hit
                    # end position is calculated using the cigar string
                    # of the hit
                    hit_end = hit_pos + get_cigar_length(cigar) - 1
                    # create region keys required for sequence retrieval
                    hit_region_key = (hit_chrom + ":" + str(hit_pos)
                                      + "-" + str(hit_end))
                    if strand:
                        orient = "forward"
                    else:
                        orient = "reverse"
                    for k in gene_dict:
                        copy_chrom = gene_dict[k]["chrom"]
                        copy_begin = gene_dict[k]["capture_start"]
                        copy_end = gene_dict[k]["capture_end"]
                        if (((copy_chrom == hit_chrom) and
                                (copy_begin - tol
                                 < hit_pos < copy_end + tol))
                                or ((copy_chrom == hit_chrom) and
                                    (copy_begin - tol
                                     < hit_end < copy_end + tol))):
                            secondary_haplotypes.append(
                                [k[0], k[1], h, hit_chrom, hit_pos, hit_end,
                                 orient, hit_region_key]
                            )
                            haplotypes[m][h]["mapped"] = True
                            secondary_haplotype_dict[h] = (
                                haplotypes[m].pop(h)
                            )
                            secondary_haplotype_found = True
                            break
    if len(secondary_haplotypes) > 0:
        secondary_haplotypes = pd.DataFrame(
            secondary_haplotypes,
            columns=["gene", "copy", "original_hap_ID", "chrom",
                     "capture_start", "capture_end", "orientation",
                     "region_key"]
        )
        secondary_haplotypes = secondary_haplotypes.merge(
            gene_df[["gene", "copy", "copyname", "mip_number", "sub_number"]]
        )
        secondary_haplotypes.to_csv(wdir + "secondary_haplotypes.csv")

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


def align_haplotypes(
        settings, target_actions=["unmask", "multiple"],
        query_actions=["unmask"],
        output_format="general:name1,text1,name2,text2,diff,score",
        alignment_options=["--noytrim"], identity=75, coverage=75
        ):
    """ Get a haplotypes dict and a call_info dict, align each haplotype to
    reference sequences from the call_info dict."""
    wdir = settings["workingDir"]
    haplotypes_file = os.path.join(wdir, settings["tempHaplotypesFile"])
    with open(haplotypes_file) as infile:
        haplotypes = json.load(infile)
    species = settings["species"]
    alignment_dir = wdir + settings["alignmentDir"]
    num_processor = int(settings["processorNumber"])
    command_list = []
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    # create alignment dir if it does not exist
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    for m in haplotypes:
        # create a fasta file for each mip that contains all haplotype
        # sequences for that mip
        haplotype_fasta = alignment_dir + m + ".haps"
        with open(haplotype_fasta, "w") as outfile:
            outfile_list = []
            for h in haplotypes[m]:
                outlist = [">", h, "\n", haplotypes[m][h]["sequence"]]
                outfile_list.append("".join(outlist))
            outfile.write("\n".join(outfile_list))
        haplotype_fasta = m + ".haps"
        # create a reference file for each mip that contains reference
        # sequences for each paralog copy for that mip
        reference_fasta = alignment_dir + m + ".refs"
        with open(reference_fasta, "w") as outfile:
            outfile_list = []
            gene_name = m.split("_")[0]
            for c in call_info[gene_name][m]["copies"]:
                c_ori = call_info[gene_name][m]["copies"][c]["orientation"]
                c_seq = call_info[gene_name][m]["copies"][c][
                    "capture_sequence"]
                if c_ori == "reverse":
                    c_seq = reverse_complement(c_seq)
                outlist = [">", m + "_" + c, "\n", c_seq]
                outfile_list.append("".join(outlist))
            outfile.write("\n".join(outfile_list))
        # name of the alignment output file for the mip
        output_file = m + ".aligned"
        # create the list to be passed to the alignment worker function
        command = [haplotype_fasta, alignment_dir, output_file,
                   reference_fasta, target_actions, query_actions, identity,
                   coverage, output_format, alignment_options, species]
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
    alignment_dir = os.path.join(wdir, settings["alignmentDir"])
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    temp_haplotypes_file = os.path.join(wdir, settings["tempHaplotypesFile"])
    with open(temp_haplotypes_file) as infile:
        haplotypes = json.load(infile)
    alignments = {}
    inverted_alignments = []
    problem_alignments = []
    problem_snps = []
    for m in haplotypes:
        # each mip has all its haplotypes and reference sequences aligned
        # in mipname.aligned file.
        with open(os.path.join(alignment_dir, m + ".aligned")) as al_file:
            for line in al_file:
                problem_al = False
                if not line.startswith("#"):
                    # each line of the alignment file includes an alignment
                    # between the reference copy sequences of a mip
                    # and a haplotype sequence
                    newline = line.strip().split("\t")
                    gene_name = newline[0].split("_")[0]
                    m_name = "_".join(newline[0].split("_")[:-1])
                    ref_copy = newline[0].split("_")[-1]
                    rf_ori = call_info[gene_name][m_name]["copies"][ref_copy][
                        "orientation"]
                    # aligned part of the reference sequence with gaps
                    ref_al = newline[1].upper()
                    if rf_ori == "reverse":
                        ref_al = reverse_complement(ref_al)
                    # aligned part of the reference without gaps
                    ref_used = ref_al.translate(str.maketrans({"-": None}))
                    ref_used = ref_used.upper()
                    hap_name = newline[2]
                    # aligned part of the haplotype with gaps
                    hap_al = newline[3].upper()
                    if rf_ori == "reverse":
                        hap_al = reverse_complement(hap_al)
                    # aligned part of the haplotype without gaps
                    hap_used = hap_al.translate(str.maketrans({"-": None}))
                    hap_used = hap_used.upper()
                    # alignment diff (.for match, : and X mismatch, - gap)
                    diff = newline[4]
                    if rf_ori == "reverse":
                        diff = diff[::-1]
                    score = int(newline[5])
                    # full haplotype sequence
                    hap_seq = haplotypes[m][hap_name]["sequence"].upper()
                    # full reference sequence
                    ref_seq = call_info[gene_name][m]["copies"][ref_copy][
                        "capture_sequence"].upper()
                    # index of where in full reference the alignment begins
                    ref_align_begin = ref_seq.find(ref_used)
                    # index of where in full reference the alignment ends
                    ref_align_end = ref_align_begin + len(ref_used)
                    # index of where in full haplotype sequence the alignment
                    # begins
                    hap_align_begin = hap_seq.find(hap_used)
                    # if the alignment is inverted, i.e. there is a reverse
                    # complement alignment with a significant score, the find
                    # method will not find the haplotype sequence in query or
                    # the target sequence in reference # and return -1. These
                    # alignments have been happening when one copy differs so
                    # much from another, an inverted alignment scores better.
                    # These should be ignored because the real copy the
                    # haplotype comes from will have a better score. However,
                    # there can theoretically be an inversion within a capture
                    # region that produces a legitimate inverted alignment.
                    # Therefore these alignments may be inspected if desired.
                    # We will keep such alignments in a dictionary and save.
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
                    # index of where in full haplotype sequence the alignment
                    # ends
                    hap_align_end = hap_align_begin + len(hap_used)
                    # deal with any existing flanking deletions/insertions
                    # is there any unaligned sequence on the left of alignment
                    left_pad_len = max([hap_align_begin, ref_align_begin])
                    left_pad_diff = abs(hap_align_begin - ref_align_begin)
                    left_pad_ref = ""
                    left_pad_hap = ""
                    left_pad_ref_count = 0
                    left_pad_hap_count = 0
                    # where there are insertions on left, fill the other pad
                    # with gaps
                    for i in range(hap_align_begin - ref_align_begin):
                        # only when ref_align_begin is smaller, we need to pad
                        # left_pad_ref
                        left_pad_ref = "-" + left_pad_ref
                        left_pad_hap = hap_seq[i] + left_pad_hap
                        # counting how many bases from hap_seq is used for
                        # padding
                        left_pad_hap_count += 1
                    # do the same for haplotype sequence
                    for i in range(ref_align_begin - hap_align_begin):
                        # only when ref_align_begin is smaller, we need to pad
                        # left_pad_ref
                        left_pad_hap += "-"
                        left_pad_ref += ref_seq[i]
                        # counting how many bases from ref_seq is used for
                        # padding
                        left_pad_ref_count += 1
                    # add to left_pads the sequences which are there but did
                    # not align
                    for i in range(left_pad_len - left_pad_diff):
                        left_pad_ref += ref_seq[i + left_pad_ref_count]
                        left_pad_hap += hap_seq[i + left_pad_hap_count]
                    # add the left padding info to the alignment
                    for i in range(0, len(left_pad_hap))[::-1]:
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
                        # counting how many bases from hap_seq is used for
                        # padding
                        right_pad_hap_count += 1
                    # do the same for haplotype sequence
                    for i in range(right_pad_ref_len - right_pad_hap_len):
                        right_pad_hap = "-" + right_pad_hap
                        right_pad_ref = ref_seq[-i - 1] + right_pad_ref
                        right_pad_ref_count += 1
                    # add to right the sequences which are there but did not
                    # align
                    for i in range(right_pad_len - right_pad_diff):
                        right_pad_ref = (ref_seq[-i - right_pad_ref_count - 1]
                                         + right_pad_ref)
                        right_pad_hap = (hap_seq[-i - right_pad_hap_count - 1]
                                         + right_pad_hap)
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
                    # hap sequence is accounted for and not just the aligned
                    # part ref_name, ref_copy, ref_seq, ref_al
                    # hap_name, hap_seq, hap_al, diff, score have information
                    # we'll use
                    c_name = ref_copy
                    h_name = hap_name
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
                        # each difference between the hap and ref can be an
                        # indel ("-") or a snp (":" or "x") or the same as
                        # the reference ("."). When dealing with indels, it is
                        # best to call consecutive indels as a cumulative indel
                        # rather than individual indels, i.e. AAA/--- instead
                        # of A/-, A/-, A/- because if we are looking for a
                        # frameshift insertion A/-, having AAA/--- means we
                        # don't observe the frameshift. But if it is kept as
                        # three A/-'s then it looks like the frameshift
                        # mutation is there.
                        if d == "-":
                            # if an indel is encountered, we'll keep track of
                            # it until the end of the indel. That is, when
                            # d != "-"
                            indel_count += 1
                            if hap_al[i] == "-":
                                # if a deletion, hap sequence should have "-"
                                indel_types.append("del")
                                indels.append(ref_al[i])
                                # in cases of deletions, we increment the
                                # genomic pos because the reference has a
                                # nucleotide in this position.
                                genomic_pos += 1
                            elif ref_al[i] == "-":
                                indel_types.append("ins")
                                indels.append(hap_al[i])
                                hap_index += 1
                                # in cases of insertions, we don't increment
                                # the genomic pos because the reference has no
                                # nucleotide in this position insAAA would have
                                # the same start and end positions
                            else:
                                # if neither hap nor ref has "-" at this
                                # position there is a disagreement between the
                                # alignment and the sequences.
                                print(("For the haplotype {} the alignment "
                                       " shows an indel but sequences do not."
                                       " This haplotype  will not have "
                                       "variant calls.").format(h_name))
                                problem_al = True
                                break
                        else:
                            # if the current diff is not an indel,
                            # check if there is preceeding indel
                            if len(indels) > 0:
                                # there should only be a del or ins preceding
                                # this base
                                if len(set(indel_types)) != 1:
                                    # Consecutive insertions and deletions
                                    print(
                                        ("For the haplotype {} there are "
                                         "consecutive insertions and "
                                         "deletions. This haplotype will not"
                                         "have variant calls.").format(h_name)
                                    )
                                    problem_al = True
                                else:
                                    indel_type = list(set(indel_types))[0]
                                    indel_length = len(indels)
                                    # genomic_pos is the current position
                                    # since this position is not an indel,
                                    # indel has ended 1 nucleotide prior to
                                    # this position.
                                    indel_end = genomic_pos - 1
                                    indel_seq = "".join(indels)
                                    buffer_seq = "".join(["-" for j in
                                                         range(indel_length)])
                                    if indel_type == "del":
                                        indel_begin = (genomic_pos
                                                       - indel_length)
                                        ref_base = indel_seq
                                        hap_base = buffer_seq
                                        h_index = [hap_index, hap_index - 1]
                                    else:
                                        # if the preceding indel was an
                                        # insertion the start and end positions
                                        # are the same
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
                    # all differences in between the ref and hap has been found
                    # loop through differences and assign values for
                    # ref_base, hap_base and hap_index.
                    for d in differences:
                        if copy_ori == "reverse":
                            # revert bases to their original strand (-)
                            d["ref_base"] = reverse_complement(d["ref_base"])
                            d["hap_base"] = reverse_complement(d["hap_base"])
                            d["hap_index"] = [-1 * d["hap_index"][0] - 1,
                                              -1 * d["hap_index"][1] - 1]

                    # create a dictionary that holds all the alignment
                    # information for the mip and haplotype
                    al_dict = {"gene_name": gene_name,
                               "mip_name": m,
                               "haplotype_ID": hap_name,
                               "score": score,
                               "differences": differences,
                               "aligned_hap": hap_al,
                               "aligned_ref": ref_al,
                               "diff": diff}
                    # also report alignments that had any problems
                    if problem_al:
                        problem_alignments.append(al_dict)
                    try:
                        alignments[hap_name][ref_copy].append(al_dict)
                    except KeyError:
                        try:
                            alignments[hap_name][ref_copy] = [al_dict]
                        except KeyError:
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
            except KeyError:
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
            except KeyError:
                continue
        # replace the problem dictionary with the updated version
        temp_dict = {}
        for a in probs:
            if a != "remove":
                hap_name = a["haplotype_ID"]
                try:
                    temp_dict[hap_name].append(a)
                except KeyError:
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


def update_aligned_haplotypes(settings):
    """
    Update haplotypes with information from the alignment results.
    Find which paralog copy the haplotype best maps to using the alignment
    scores from lastZ.
    """
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
    # update each haplotype with alignment information
    for m in haplotypes:
        gene_name = m.split("_")[0]
        for h in haplotypes[m]:
            # create a copy dict for each haplotype for each possible
            # paralog gene copy that haplotype may belong to
            copies = haplotypes[m][h]["copies"] = {}
            # get the alignment for this haplotype from alignment dict
            try:
                align = alignments[h]
            except KeyError:
                haplotypes[m][h]["mapped"] = False
                continue
            # update copies dict with alignment information
            for c in align:
                # update haplotype with alignment information
                copies[c] = {"score": align[c]["score"]}
            # sort copies considering alignment scores
            copy_keys_sorted = sorted(copies,
                                      key=lambda a: copies[a]["score"])
            # below code prevents alt contigs to be the best mapping copy
            # unless it is the only mapping copy. Uncomment if needed.
            """
            copy_keys_sorted_temp = copy.deepcopy(copy_keys_sorted)
            # remove alt contigs
            for cop_ind in range(len(copy_keys_sorted)):
                cop = copy_keys_sorted[cop_ind]
                if "alt" in call_info[gene_name][m]["copies"][cop[0]]["chrom"]:
                    copy_keys_sorted[cop_ind] = "remove"
            copy_keys_sorted = [cop_key for cop_key in copy_keys_sorted
                                if cop_key != "remove"]
            if len(copy_keys_sorted) == 0:
                copy_keys_sorted = copy_keys_sorted_temp
            """
            # pick best scoring copies
            # last item in copy_keys_sorted is the best
            best_copy = copy_keys_sorted[-1]
            # create a list of copies that has the best score
            best_copies = [cop for cop in copy_keys_sorted
                           if (copies[cop]["score"]
                               == copies[best_copy]["score"])]
            # create a map dict to be added to haplotype information
            # extract copy keys of best copies
            temp_dic = {}
            for c in best_copies:
                temp_dic[c] = {"copy_name": call_info[gene_name][m]["copies"][
                               c]["copyname"],
                               "differences": alignments[h][c]["differences"],
                               "chrom": call_info[gene_name][m]["copies"][c][
                               "chrom"]}
            haplotypes[m][h]["mapped_copies"] = temp_dic
            # create a single copy name for the haplotype such as HBA1_C0
            # if there are multiple copies that the haplotype matched equally
            # well name will be the compound name of all best mapping copies,
            # e.g. HBA1_C0_HBA2_C1
            mapped_copy_names = []
            for k in sorted(temp_dic.keys()):
                mapped_copy_names.append("_".join([temp_dic[k]["copy_name"],
                                                   k]))
            haplotypes[m][h]["copy_name"] = "_".join(mapped_copy_names)
    temp_mapped_haps_file = wdir + settings["tempMappedHaplotypesFile"]
    with open(temp_mapped_haps_file, "w") as outfile:
        json.dump(haplotypes, outfile, indent=1)
    return


def update_unique_haplotypes(settings):
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
    except IOError:
        unique_haplotypes = {}
    try:
        with open(sequence_to_haplotype_file) as infile:
            sequence_to_haplotype = json.load(infile)
    except IOError:
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
            except KeyError:
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
        except KeyError:
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
    except IOError:
        variation = {}
    var_key_to_uniq_file = wdir + settings["variationKeyToUniqueKey"]
    try:
        with open(var_key_to_uniq_file) as infile:
            var_key_to_uniq = json.load(infile)
    except IOError:
        var_key_to_uniq = {}
    outfile_list = ["##fileformat=VCFv4.1"]
    outfile_list.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT"]))
    temp_variations = []
    for m in haplotypes:
        for h in haplotypes[m]:
            if haplotypes[m][h]["mapped"]:
                try:
                    haplotypes[m][h]["left_normalized"]
                except KeyError:
                    left_normalized = True
                    for c in haplotypes[m][h]["mapped_copies"]:
                        differences = haplotypes[m][h]["mapped_copies"][c][
                            "differences"]
                        for d in differences:
                            var_key = d["vcf_raw"]
                            try:
                                uniq_var_key = var_key_to_uniq[var_key]
                                d["annotation"] = variation[uniq_var_key]
                                d["vcf_normalized"] = uniq_var_key
                            except KeyError:
                                left_normalized = False
                                temp_variations.append(var_key)
                    if left_normalized:
                        haplotypes[m][h]["left_normalized"] = True
    temp_variations = [temp_var.split(":")
                       for temp_var in set(temp_variations)]
    temp_variations = [[v[0], int(v[1])] + v[2:] for v in temp_variations]
    temp_variations = sorted(temp_variations, key=itemgetter(0, 1))
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
                               cwd=wdir, stdout=outfile)
    dump = subprocess.call(["bcftools", "index", "-f", raw_vcf_file + ".gz"],
                           cwd=wdir)
    unmasked_genome = get_file_locations()[species]["unmasked_fasta_genome"]
    dump = subprocess.call(["bcftools", "norm", "-f", unmasked_genome,
                            "-cw", "-w", "0",
                           "-o", norm_vcf_file, raw_vcf_file + ".gz"],
                           cwd=wdir)
    ann_db_dir = get_file_locations()[species]["annotation_db_dir"]
    ann_build = settings["annotationBuildVersion"]
    ann_protocol = settings["annotationProtocol"].replace(";", ",")
    ann_operation = settings["annotationOperation"].replace(";", ",")
    ann_nastring = settings["annotationNaString"]
    ann_out = settings["annotationOutput"]
    try:
        ann_script = settings["annotationScript"]
    except KeyError:
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
    dump = subprocess.check_call(ann_command, cwd=wdir)
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
                if normalized_key not in variation:
                    variation[normalized_key] = {
                        colnames[i]: newline[i] for i in range(len(colnames))
                    }
                line_num += 1
    if line_num != len(temp_variation_keys):
        print("There are more variation keys then annotated variants.")
    for m in haplotypes:
        for h in haplotypes[m]:
            if haplotypes[m][h]["mapped"]:
                try:
                    haplotypes[m][h]["left_normalized"]
                except KeyError:
                    for c in haplotypes[m][h]["mapped_copies"]:
                        differences = haplotypes[m][h]["mapped_copies"][c][
                            "differences"]
                        for d in differences:
                            var_key = d["vcf_raw"]
                            uniq_var_key = var_key_to_uniq[var_key]
                            d["annotation"] = variation[uniq_var_key]
                            d["vcf_normalized"] = uniq_var_key
                            annotation_dict = d["annotation"]
                            for ak in annotation_dict.keys():
                                if ak.startswith("AAChange."):
                                    annotation_dict["AAChangeClean"] = (
                                        annotation_dict.pop(ak)
                                    )
                                elif ak.startswith("ExonicFunc."):
                                    annotation_dict["ExonicFunc"] = (
                                        annotation_dict.pop(ak)
                                    )
                                elif ak.startswith("Gene."):
                                    annotation_dict["GeneID"] = (
                                        annotation_dict.pop(ak)
                                    )
                    haplotypes[m][h]["left_normalized"] = True
    with open(unique_haplotype_file, "w") as outfile:
        json.dump(haplotypes, outfile)
    with open(variation_file, "w") as outfile:
        json.dump(variation, outfile)
    with open(var_key_to_uniq_file, "w") as outfile:
        json.dump(var_key_to_uniq, outfile)
    try:
        m_snps = int(settings["mergeSNPs"])
    except KeyError:
        m_snps = False
    if m_snps:
        dump = merge_snps(settings)
    return


def make_snp_vcf(variant_file, haplotype_file, call_info_file,
                 haplotype_counts_file, vcf_chrom, barcode_count_file,
                 min_cov, min_count, min_freq, vcf_file,
                 settings_file, header_count=11):
    """
    Create a vcf file for SNV only. This will be integrated to process_results
    in the future.
    """
    # Load variant count table
    variant_counts = pd.read_csv(variant_file,
                                 header=list(range(header_count)),
                                 index_col=0)
    # Add variant type to tables, convert position to integer
    cols = variant_counts.columns
    new_index = pd.MultiIndex.from_tuples(
        [(c[0], ) + (int(float(c[1])), ) + c[2:] + ("SNV", )
         if len(c[3]) == len(c[4])
         else (c[0], ) + (int(float(c[1])), ) + c[2:] + ("indel", )
         for c in cols],
        names=cols.names + ["Variant Type"])
    variant_counts.columns = new_index
    # filter indels
    variant_counts = variant_counts.xs("SNV", level="Variant Type", axis=1)
    # load haplotype dict
    with open(haplotype_file) as infile:
        haplotypes = json.load(infile)
    # load call info
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    # Load haplotype counts per sample
    haplotype_counts = pd.read_csv(haplotype_counts_file)
    # Add "copy" information to haplotype counts
    hap_copies = []
    for m in haplotypes:
        g = m.split("_")[0]
        for h in haplotypes[m]:
            try:
                mc = haplotypes[m][h]["mapped_copies"]
                for c in mc:
                    hap_copies.append([h, c,
                                       call_info[g][m]["copies"][c]["chrom"]])
            except KeyError:
                continue
    hap_copies = pd.DataFrame(hap_copies, columns=["Haplotype ID",
                                                   "Copy", "CHROM"])
    # Get all variant positions across the data set
    variant_positions = set()
    cols = variant_counts.columns
    for c in cols:
        pos = c[1]
        ref = c[3]
        alt = c[4]
        len_diff = len(ref) - len(alt)
        if len_diff > 0:
            vp = set(range(pos, pos + len_diff + 1))
        else:
            vp = set([pos])
        variant_positions.update(vp)
    variant_position_set = variant_positions
    # load the probe set dictionary to extract the
    # probes that were used in this run
    settings = get_analysis_settings(settings_file)
    probe_sets_file = settings["mipSetsDictionary"]
    probe_set_keys = settings["mipSetKey"]
    used_probes = set()
    for psk in probe_set_keys:
        with open(probe_sets_file) as infile:
            used_probes.update(json.load(infile)[psk])
    # Create a position_to_mip dictionary that maps each genomic position
    # to all MIPs covering that position
    mip_positions = {}
    position_to_mip = {}
    for g in call_info:
        for m in call_info[g]:
            if m in used_probes:
                for c in call_info[g][m]["copies"]:
                    chrom = call_info[g][m]["copies"][c]["chrom"]
                    if chrom == vcf_chrom:
                        start = call_info[g][m]["copies"][c]["capture_start"]
                        end = call_info[g][m]["copies"][c]["capture_end"]
                        cov_pos = variant_position_set.intersection(
                            range(start, end + 1)
                        )
                        mip_positions[(m, c)] = cov_pos
                        for p in cov_pos:
                            try:
                                position_to_mip[p].add((m, c))
                            except KeyError:
                                position_to_mip[p] = set([(m, c)])
    # Create a dataframe that maps whether a genomic position is the same
    # as the reference or not for each haplotype
    references = []
    for m in haplotypes:
        g = m.split("_")[0]
        cops = call_info[g][m]["copies"]
        copy_keys = list(call_info[g][m]["copies"].keys())
        for c in copy_keys:
            mp = mip_positions[(m, c)]
            capture_chrom = cops[c]["chrom"]
            if capture_chrom == vcf_chrom:
                for h in haplotypes[m]:
                    hap = haplotypes[m][h]
                    refs = {p: 1 for p in mp}
                    try:
                        hap_copy = hap["mapped_copies"][c]
                    except KeyError:
                        continue
                    else:
                        for d in hap_copy["differences"]:
                            vcf_norm = d["vcf_normalized"].split(":")
                            pos = int(vcf_norm[1])
                            r = vcf_norm[3]
                            for j in range(len(r)):
                                refs[pos + j] = 0
                    references.extend([(h, c) + item for item in refs.items()])
    references = pd.DataFrame(
        references, columns=["Haplotype ID", "Copy", "POS", "Reference"]
    )
    # Update the haplotype count dataframe to include each (variant) position's
    # reference and non-reference status.
    references = references.merge(haplotype_counts[
        ["Haplotype ID", "Copy", "Sample ID", "Barcode Count"]
        ])
    # Update the dataframe with Reference base counts by multiplying
    # the haplotype's count by reference count for the haplotype
    references["Ref Count"] = (references["Reference"]
                               * references["Barcode Count"])
    # Create a dictionary that maps each sample's reference base counts at
    # each position of interest
    counts = references.groupby(
        ["Sample ID", "POS"]
    )["Ref Count"].sum().to_dict()
    # Load the barcode counts table, which has the total barcode count for a
    # given MIP for each sample
    barcode_counts = pd.read_csv(barcode_count_file,
                                 header=[0, 1], index_col=0)
    # Create a coverage dictionary from the barcode count table
    cov_dict = barcode_counts.fillna(0).loc[
        :, list(mip_positions.keys())
    ].to_dict()
    # Create reference count and coverage tables corresponding to the variant
    # table we had earlier
    rows = variant_counts.index
    columns = variant_counts.columns
    reference_counts = []
    coverage = []
    for sample in rows:
        sample_refs = []
        sample_cov = []
        for variant in columns:
            pos = variant[1]
            try:
                ref_count = counts[(sample, pos)]
            except KeyError:
                ref_count = 0
            cov = 0
            for k in position_to_mip[pos]:
                cov += cov_dict[k][sample]
            sample_refs.append(ref_count)
            sample_cov.append(cov)
        reference_counts.append(sample_refs)
        coverage.append(sample_cov)
    reference_counts = pd.DataFrame(reference_counts,
                                    columns=columns, index=rows)
    coverage = pd.DataFrame(coverage, columns=columns, index=rows)

    def collapse_snps(g):
        """Take a group of variants on the same position, return a merged
        dataframe which has the allele counts as a comma separated string
        for each  allele in a given position."""
        gv = g.columns.get_level_values
        ref = gv("REF")[0]
        alts = ",".join(gv("ALT"))
        idx = pd.MultiIndex.from_tuples([(".", ref, alts)],
                                        names=["ID", "REF", "ALT"])
        vals = [",".join(map(str, map(int, v))) for v in g.values]
        return pd.DataFrame(vals, columns=idx)
    # group variants on the position to merge multiallelic loci
    collapsed_vars = variant_counts.groupby(level=["CHROM", "POS"],
                                            axis=1).apply(collapse_snps)
    # group coverage and reference counts for multiallelic loci
    collapsed_refs = reference_counts.groupby(
        level=["CHROM", "POS", "ID", "REF"], axis=1
    ).first()
    collapsed_cov = coverage.groupby(level=["CHROM", "POS", "ID", "REF"],
                                     axis=1).first()
    # index of collapsed_vars is lost after groupby operation
    # although the operation does not impact the rows. We'll recover
    # the index from the other data frames which all have the same
    # row indices
    collapsed_vars.index = collapsed_cov.index
    collapsed_vars.sort_index(axis=1, level=["CHROM", "POS"], inplace=True)
    collapsed_refs.sort_index(axis=1, level=["CHROM", "POS"], inplace=True)
    collapsed_cov.sort_index(axis=1, level=["CHROM", "POS"], inplace=True)
    collapsed_refs.columns = collapsed_vars.columns
    collapsed_cov.columns = collapsed_vars.columns
    # merge the variant, reference and coverage count tables to get a "variant
    # string" for each locus and for each sample.
    collapsed_merge = (collapsed_refs.astype(int).astype(str)
                       + "," + collapsed_vars
                       + ":" + collapsed_cov.astype(int).astype(str))

    def call_genotype(s, min_cov, min_count, min_freq):
        """Call genotypes from the variant strings that are in the form:
        ref_count,allele1_count,allele2_count,...:coverage."""
        sp = s.split(":")
        cov = int(sp[-1])
        allele_counts = list(map(int, sp[0].split(",")))
        if cov < min_cov:
            genotypes = []
        else:
            genotypes = []
            for i in range(len(allele_counts)):
                ac = allele_counts[i]
                if ((ac >= min_count) & ((ac/cov) >= min_freq)):
                    genotypes.append(str(i))
        if len(genotypes) == 0:
            genotypes = "."
        else:
            genotypes = "/".join(genotypes) + ":" + s
        return genotypes
    # call genotypes
    vcf = collapsed_merge.applymap(lambda a: call_genotype(
        a, min_cov, min_count, min_freq)
    ).T
    # Update columns of vcf table to remove the column's level name
    vcf_samples = vcf.columns.tolist()
    vcf.columns = vcf_samples
    # Add vcf filler text
    vcf["QUAL"] = "."
    vcf["FILTER"] = "."
    vcf["INFO"] = "."
    vcf["FORMAT"] = "GT:AD:DP"
    vcf = vcf[["QUAL", "FILTER", "INFO", "FORMAT"] + vcf_samples]
    vcf_header = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
        '"Allelic depths for the ref and alt alleles in that order">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
        '"Total read depth (coverage) at this position.>']
    # save vcf file
    with open(vcf_file, "w") as outfile:
        outfile.write("\n".join(vcf_header) + "\n")
        vcf.reset_index().rename(columns={"CHROM": "#CHROM"}).to_csv(
            outfile, index=False, sep="\t"
        )


def make_chrom_vcf(wdir, header_count, min_cov=1, min_count=1, min_freq=0):
    variant_counts = pd.read_csv(
        wdir + "variant_table.csv",
        header=list(range(header_count)), index_col=0
    )
    gb = variant_counts.groupby(level="CHROM", axis=1)
    for chrom in gb.groups.keys():
        gb.get_group(chrom).to_csv(os.path.join(wdir, chrom + ".var.csv"))
    variant_coverage = pd.read_csv(
        wdir + "variant_coverage_table.csv",
        header=list(range(header_count)), index_col=0
    )
    gb = variant_coverage.groupby(level="CHROM", axis=1)
    for chrom in gb.groups.keys():
        gb.get_group(chrom).to_csv(os.path.join(wdir, chrom + ".cov.csv"))

    with open(wdir + "unique_haplotype.dic") as infile:
        haplotypes = json.load(infile)

    with open("/opt/project_resources/mip_ids/call_info.json") as infile:
        call_info = json.load(infile)
    hap_counts = pd.read_csv(wdir + "haplotype_counts.csv")

    chrom_haplotypes = {}
    for m in haplotypes:
        g = m.split("_")[0]
        for h in haplotypes[m]:
            try:
                mc = haplotypes[m][h]["mapped_copies"]
            except KeyError:
                continue
            else:
                for c in mc:
                    chrom = call_info[g][m]["copies"][c]["chrom"]
                    mc[c]["copy_chrom"] = chrom
                    try:
                        chrom_haplotypes[chrom][m][h] = haplotypes[m][h]
                    except KeyError:
                        try:
                            chrom_haplotypes[chrom][m] = {h: haplotypes[m][h]}
                        except KeyError:
                            chrom_haplotypes[chrom] = {
                                m: {h: haplotypes[m][h]}}
    for chrom in chrom_haplotypes:
        with open(os.path.join(wdir, chrom + ".haps.json"), "w") as outfile:
            json.dump(chrom_haplotypes[chrom], outfile)
    hap_copies = []
    for m in haplotypes:
        g = m.split("_")[0]
        for h in haplotypes[m]:
            try:
                mc = haplotypes[m][h]["mapped_copies"]
                for c in mc:
                    hap_copies.append(
                        [h, c, call_info[g][m]["copies"][c]["chrom"]])
            except KeyError:
                continue
    hap_copies = pd.DataFrame(hap_copies,
                              columns=["Haplotype ID", "Copy", "CHROM"])
    haplotype_counts = hap_counts.merge(hap_copies)
    gb = haplotype_counts.groupby("CHROM")
    for chrom in gb.groups.keys():
        gb.get_group(chrom).to_csv(os.path.join(wdir,
                                                chrom + ".hap_counts.txt"))

    call_info_file = "/opt/project_resources/mip_ids/call_info.json"
    chromosomes = set(haplotype_counts["CHROM"])
    for chrom in sorted(chromosomes):
        haplotype_file = wdir + "/" + chrom + ".haps.json"
        variant_file = wdir + "/" + chrom + ".var.csv"
        haplotype_counts_file = wdir + "/" + chrom + ".hap_counts.txt"
        haplotype_file = wdir + "/" + chrom + ".haps.json"
        vcf_chrom = chrom
        barcode_count_file = wdir + "barcode_counts.csv"
        vcf_file = wdir + chrom + ".vcf"
        settings_file = wdir + "settings.txt"
        try:
            make_snp_vcf(vcf_file=vcf_file,
                         haplotype_file=haplotype_file,
                         variant_file=variant_file,
                         haplotype_counts_file=haplotype_counts_file,
                         vcf_chrom=vcf_chrom,
                         barcode_count_file=barcode_count_file,
                         call_info_file=call_info_file,
                         min_cov=min_cov,
                         min_count=min_count,
                         min_freq=min_freq,
                         settings_file=settings_file,
                         header_count=header_count)
        except FileNotFoundError:
            continue


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
                    sample_sheets=None,
                    meta_files=[],
                    targets_file=None,
                    target_join="union"):
    settings = get_analysis_settings(wdir + settings_file)
    if sample_sheets is None:
        sample_sheets = [wdir + "samples.tsv"]
    ##########################################################
    ##########################################################
    # Process 1: use sample sheets, sample sets and meta files
    # to determine which data points from the mipster file
    # should be used, print relevant statistics.
    ##########################################################
    ##########################################################
    # process sample sheets
    run_meta = pd.concat(
        [pd.read_table(s) for s in sample_sheets],
        ignore_index=True
    )
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
    # drop duplicate values originating from
    # multiple sequencing runs of the same libraries
    run_meta = run_meta.drop_duplicates()
    run_meta = run_meta.groupby(
        ["Sample ID", "Library Prep"]
    ).first().reset_index()
    run_meta.to_csv(wdir + "run_meta.csv")
    # load meta data for samples, if given. Use a mock field if not given.
    try:
        sample_meta = pd.concat(
            [pd.read_table(f) for f in meta_files],
            join="outer",
            ignore_index=True
        )
    except ValueError:
        # if no meta files were provided, create a dummy
        # meta dataframe
        sample_meta = copy.deepcopy(run_meta[["Sample Name"]])
        sample_meta["Meta"] = "Meta"
    # Pandas reads sample names that are numbers as numbers
    # these should be converted to string for consistency across samples.
    sample_meta["Sample Name"] = sample_meta["Sample Name"].astype(str)
    sample_meta = sample_meta.groupby(["Sample Name"]).first().reset_index()
    # Merge Sample meta data and run data
    merged_meta = pd.merge(run_meta, sample_meta,
                           on="Sample Name",
                           how="inner")
    merged_meta.to_csv(wdir + "merged_meta.csv")
    print(("{} out of {} samples has meta information and"
           " will be used for analysis.").format(
        merged_meta.shape[0], run_meta.shape[0]
    ))
    # get used sample ids
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
    # keep all variant in all haplotypes in a list
    variation_list = []
    # keep haplotypes that are the same as reference genome
    # in the reference list
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
                copy_chrom = hap["mapped_copies"][c]["chrom"]
                # go through all differences from reference genome
                # get a subset of information included in the
                # haplotype dictionary
                if len(copy_differences) == 0:
                    reference_list.append([hid, c, multi_mapping, copy_chrom])
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
                    except KeyError:
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
                                 g, m, c, hid,
                                 raw_key,
                                 original_pos,
                                 start_index,
                                 end_index,
                                 multi_mapping]
                    try:
                        for ak in annotation_keys:
                            temp_list.append(d["annotation"][ak])
                    except NameError:
                        annotation_keys = list(d["annotation"].keys())
                        for ak in annotation_keys:
                            temp_list.append(d["annotation"][ak])
                    variation_list.append(temp_list)
    # create pandas dataframes for variants
    colnames = ["VKEY", "CHROM", "POS", "ID", "REF", "ALT",
                "Gene", "MIP", "Copy", "Haplotype ID",
                "RAW_VKEY", "Original Position", "Start Index",
                "End Index", "Multi Mapping"]
    colnames = colnames + annotation_keys
    variation_df = pd.DataFrame(variation_list,
                                columns=colnames)
    # create pandas dataframe for reference haplotypes
    reference_df = pd.DataFrame(reference_list,
                                columns=["Haplotype ID",
                                         "Copy",
                                         "Multi Mapping",
                                         "Chrom"])

    # create a dataframe for all mapped haplotypes
    mapped_haplotype_df = pd.concat(
        [variation_df.groupby(
            ["Haplotype ID", "Copy", "Multi Mapping", "CHROM"]
        ).first().reset_index().rename(columns={"CHROM": "Chrom"})[
            ["Haplotype ID", "Copy", "Multi Mapping", "Chrom"]
        ], reference_df], ignore_index=True
    )
    print(
        ("There are {mh.shape[0]} mapped and {um} unmapped (off target)"
         " haplotypes.").format(mh=mapped_haplotype_df, um=unmapped)
    )
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
    # use only the data corresponding to mapped haplotypes
    # filtering the off target haplotypes.
    mapped_results = raw_results.merge(mapped_haplotype_df, how="inner")
    print(("There are {rr.shape[0]} data points in raw data,"
           " {mr.shape[0]} are mapped to genome and their targets.").format(
        rr=raw_results,
        mr=mapped_results
    ))
    # rename some columns for better visualization in tables
    mapped_results.rename(
        columns={"sample_name": "Sample ID",
                 "mip_name": "MIP",
                 "gene_name": "Gene",
                 "barcode_count": "Barcode Count",
                 "read_count": "Read Count"},
        inplace=True
    )
    # Try to estimate the distribution of data that is mapping
    # to multiple places in the genome.
    # This is done in 4 steps.
    # 1) Get uniquely mapping haplotypes and barcode counts
    unique_df = mapped_results.loc[~mapped_results["Multi Mapping"]]
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
    unique_df.loc[:, "CA"] = average_copy_count
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
    ac = cc/gc
    # 4) Distribute multi mapping data proportional to
    # Paralog's copy number determined from the
    # uniquely mapping data
    multi_df = mapped_results.loc[mapped_results["Multi Mapping"]]
    if not multi_df.empty:
        # get the average copy count for the gene the haplotype belongs to
        mca = multi_df.apply(lambda r: get_copy_average(r, ac), axis=1)
        multi_df.loc[mca.index, "Copy Average"] = mca
        mca = multi_df.groupby(
            ["Sample ID", "Gene"]
        )["Copy Average"].transform(normalize_copies)
        multi_df.loc[mca.index, "CA"] = mca
        multi_df.loc[:, "Adjusted Barcode Count"] = (
            multi_df.loc[:, "Barcode Count"]
            * multi_df.loc[:, "CA"]
        )
        multi_df.loc[:, "Adjusted Read Count"] = (
            multi_df.loc[:, "Read Count"]
            * multi_df.loc[:, "CA"]
        )
    # Combine unique and multimapping data
    combined_df = pd.concat([unique_df, multi_df], ignore_index=True)
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
            raw_results[["read_count", "barcode_count"]].sum(),
            combined_df[["Read Count", "Barcode Count"]].sum().astype(int)
        )
    )
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
    ##########################################################
    # Haplotype Filters from the settings file
    haplotype_min_barcode_filter = int(settings["minHaplotypeBarcodes"])
    haplotype_min_sample_filter = int(settings["minHaplotypeSamples"])
    haplotype_min_sample_fraction_filter = float(
        settings["minHaplotypeSampleFraction"]
    )
    # Gather per haplotype data across samples
    hap_counts = combined_df.groupby(
        "Haplotype ID"
    )["Barcode Count"].sum().reset_index().rename(
        columns={"Barcode Count": "Haplotype Barcodes"})
    hap_sample_counts = combined_df.groupby("Haplotype ID")["Sample ID"].apply(
        lambda a: len(set(a))
    ).reset_index(
    ).rename(columns={"Sample ID": "Haplotype Samples"})
    num_samples = float(combined_df["Sample ID"].unique().size)
    hap_sample_counts["Haplotype Sample Fraction"] = (
        hap_sample_counts["Haplotype Samples"] / num_samples
    )
    hap_counts = hap_counts.merge(hap_sample_counts)
    hap_counts = hap_counts.loc[(hap_counts["Haplotype Samples"]
                                 >= haplotype_min_sample_filter)
                                & (hap_counts["Haplotype Sample Fraction"]
                                   >= haplotype_min_sample_fraction_filter)
                                & (hap_counts["Haplotype Barcodes"]
                                   >= haplotype_min_barcode_filter)]
    variation_df = variation_df.merge(hap_counts, how="inner")
    # Rename or remove some columns for downstream analysis
    variation_df["AA Change"] = variation_df["AAChangeClean"].apply(
        split_aa
    )
    variation_df["AA Change Position"] = variation_df["AAChangeClean"].apply(
        split_aa_pos
    )
    try:
        variation_df.drop(["Chr", "Ref", "Alt"], axis=1, inplace=True)
    except KeyError:
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
        targets = pd.read_table(targets_file).drop_duplicates()
        join_dict = {"intersection": "inner",
                     "union": "outer",
                     "targets": "right",
                     "data": "left"}
        targets["Targeted"] = "Yes"
        variation_df = variation_df.merge(
            targets,
            how=join_dict[target_join]
        )
        variation_df["Targeted"].fillna("No", inplace=True)
        # If a reference genome locus is a mutation of interest
        # such as dhps-437, this information can be supplied
        # in targets file. The rest of the variants will be
        # assinged a False value for this.
        try:
            variation_df["Reference Resistant"].fillna("No", inplace=True)
            variation_df["Reference Resistant"] = variation_df[
                "Reference Resistant"
            ].apply(str.capitalize)
            ref_resistant = True
        except KeyError:
            ref_resistant = False
        # if target join method will be "union" or "targets"
        # the target variants will be kept even if they are not observed
        # in the entire data set. Since they are not observed, we do not have
        # information about their properties required in vcf and other formats.
        # These minimum columns must be present in the targets file.
        # They must be prosent as "Title Case" of their corresponding
        # "UPPERCASE" values in vcf: "CHROM" must be supplied as "Chrom"
        if target_join in ["union", "targets"]:
            data_keys = ["VKEY", "CHROM", "POS", "ID", "REF", "ALT"]
            target_keys = ["Vkey", "Chrom", "Pos", "Id", "Ref", "Alt"]
            for dk, tk in zip(data_keys, target_keys):
                variation_df[dk].fillna(variation_df[tk], inplace=True)
            variation_df.drop(target_keys, axis=1, inplace=True)
    else:
        variation_df["Targeted"] = "No"
        ref_resistant = False
    # each variant needs to have a name. This should be provided in
    # the target file. For those variant that are not in the target
    # file, we'll create their names by adding aa-change to gene name
    try:
        variation_df["Mutation Name"].fillna(variation_df["Gene"] + "-"
                                             + variation_df["AA Change"],
                                             inplace=True)
    except KeyError:
        variation_df["Mutation Name"] = (variation_df["Gene"] + "-"
                                         + variation_df["AA Change"])
    # "AA Change" field for noncoding variants are ".", so they will
    # not serve well as unique mutation names. These will be changed
    # by chromosome position and the base change using the "rename_noncoding"
    # function.
    variation_df.loc[
        variation_df["AA Change"] == ".",
        "Mutation Name"
    ] = variation_df.loc[
        variation_df["AA Change"] == "."
    ].apply(rename_noncoding, axis=1)

    # remove columns that will not be used after this point
    variation_df.drop(
        ["RAW_VKEY",
         "Original Position",
         "Haplotype Barcodes",
         "Haplotype Samples",
         "Haplotype Sample Fraction"],
        axis=1,
        inplace=True
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
                probe_cop.append([m, c])
        except KeyError:
            continue
    probe_cop = pd.DataFrame(probe_cop, columns=["MIP", "Copy"])
    # add a place holder column for merging probe information
    probe_cop["Temp"] = "Temp"
    # perform outer merge on results and used probes
    # to include probes that had no coverage in the results
    combined_df = combined_df.merge(probe_cop, how="outer").drop(
        "Temp", axis=1
    )
    # Fill NA values for probes with no coverage in any sample
    combined_df["Sample ID"].fillna("Temp", inplace=True)
    combined_df["Haplotype ID"].fillna(
        combined_df["MIP"] + ".0-0",
        inplace=True
    )
    combined_df["Barcode Count"].fillna(0, inplace=True)

    # Add sample and barcode depth information for each
    # variant
    variant_counts = combined_df[["Haplotype ID",
                                  "sequence_quality",
                                  "Sample ID",
                                  "Barcode Count",
                                  "Copy"]].merge(variation_df, how="right")
    # For unobserved variants, we need a place holder for Sample ID
    variant_counts["Sample ID"].fillna("Temp", inplace=True)
    variant_counts["Barcode Count"].fillna(0, inplace=True)
    variant_counts["Multi Mapping"] = variant_counts["Multi Mapping"].apply(
        lambda a: "Yes" if a is True else "No"
    )
    # Get the sample and barcode depth stats for each variant
    # and filter for given thresholds.
    # First, get the "per variant" statistics
    var_counts = variant_counts.groupby("VKEY").agg(
        {"Sample ID": lambda a: len(set(a)),
         "Barcode Count": "sum",
         "Targeted": lambda a: "Yes" if "Yes" in set(a) else "No"}
    ).rename(
        columns={"Sample ID": "Variant Samples",
                 "Barcode Count": "Variant Barcodes"}
    ).fillna(0).reset_index()
    var_counts["Variant Sample Fraction"] = var_counts[
        "Variant Samples"
    ].transform(lambda a: a/num_samples)
    # filter variants for specified criteria
    variant_min_barcode_filter = int(settings["minVariantBarcodes"])
    variant_min_sample_filter = int(settings["minVariantSamples"])
    variant_min_sample_fraction_filter = float(
        settings["minVariantSampleFraction"]
    )
    var_counts = var_counts.loc[((var_counts["Variant Samples"]
                                  >= variant_min_sample_filter)
                                & (var_counts["Variant Barcodes"]
                                   >= variant_min_barcode_filter)
                                & (var_counts["Variant Sample Fraction"]
                                   >= variant_min_sample_fraction_filter))
                                | (var_counts["Targeted"] == "Yes")]
    print("There were {} total and {} unique variants, ".format(
        variant_counts.shape[0],
        len(variant_counts["VKEY"].unique())
    ))
    # remove "Targeted" column from var counts prior to merge
    var_counts.drop("Targeted", axis=1, inplace=True)
    variant_counts = variant_counts.merge(var_counts,
                                          how="inner").drop(
        ["Variant Samples",
         "Variant Barcodes",
         "Variant Sample Fraction"],
        axis=1
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

    def get_qual(row):
        """ Calculate the sequence quality of a variant from the sequence
        quality of its parent haplotype and the variants position in the
        haplotype.
        """
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
                except IndexError:
                    continue
                    break
            # calculate quality as the mean for multi base variation
            if len(hap_qual_list) == 0:
                return np.nan
            else:
                return np.mean(hap_qual_list)
        except Exception:
            return np.nan

    # calculate variant qualities using the above function
    variant_counts["Variation Quality"] = variant_counts.apply(
        get_qual, axis=1
    )

    # filter variants for sequence quality
    variant_min_quality = int(settings["minVariantQuality"])
    variant_counts = variant_counts.loc[
        (variant_counts["Variation Quality"].isnull())
        | (variant_counts["Variation Quality"] >= variant_min_quality)
    ]
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
    # coverage calculations.
    # First, get all variant positions in the data.
    cpos = variant_counts.groupby(
        ["CHROM", "POS"]
    ).first().reset_index()[["CHROM", "POS"]]
    cpos = cpos.apply(lambda a: (a["CHROM"], a["POS"]), axis=1).values.tolist()
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
                        except KeyError:
                            position_to_mip[var_pos] = set()
                            position_to_mip[var_pos].add(
                                (m, c)
                            )
    # add any additional MIP that covers the variant positions.
    for var_pos in cpos:
        for g in call_info:
            for m in call_info[g]:
                if m in used_probes:
                    for c in call_info[g][m]["copies"]:
                        ch = call_info[g][m]["copies"][c]["chrom"]
                        cs = call_info[g][m]["copies"][c]["capture_start"]
                        ce = call_info[g][m]["copies"][c]["capture_end"]
                        if ((var_pos[0] == ch) and (cs <= var_pos[1] <= ce)):
                            try:
                                position_to_mip[var_pos].add(
                                    (m, c)
                                )
                            except KeyError:
                                position_to_mip[var_pos] = set()
                                position_to_mip[var_pos].add(
                                    (m, c)
                                )
    # Create pivot table of combined barcode counts
    # This is a per MIP per sample barcode count table
    # of the samples with sequencing data
    barcode_counts = pd.pivot_table(combined_df,
                                    index="Sample ID",
                                    columns=["MIP",
                                             "Copy"],
                                    values=["Barcode Count"],
                                    aggfunc=np.sum)
    try:
        barcode_counts.drop("Temp", inplace=True)
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
    all_barcode_counts = pd.merge(merged_meta[
                                    ["Sample ID",
                                     "replicate"]
                                  ].set_index("Sample ID"),
                                  barcode_counts,
                                  left_index=True,
                                  right_index=True,
                                  how="left")
    all_barcode_counts.drop("replicate", axis=1, inplace=True)
    # fix column names
    all_barcode_counts.columns = pd.MultiIndex.from_tuples(
        bc_cols, names=["MIP", "Copy"]
    )
    all_barcode_counts.fillna(0, inplace=True)
    print("There are {} total samples.".format(all_barcode_counts.shape[0]))
    # save barcode and haplotype count files
    combined_df.loc[combined_df["Sample ID"] != "Temp"].to_csv(
        os.path.join(wdir, "haplotype_counts.csv"), index=False
    )
    all_barcode_counts.to_csv(os.path.join(wdir, "all_barcode_counts.csv"))
    # Continue working with the barcode counts that does not include the
    # samples which did not have any data.
    barcode_counts.columns = pd.MultiIndex.from_tuples(bc_cols,
                                                       names=["MIP", "Copy"])
    barcode_counts.fillna(0, inplace=True)
    barcode_counts.to_csv(os.path.join(wdir, "barcode_counts.csv"))
    # Calculate coverage for each variant position for each sample
    bc_dict = barcode_counts.to_dict(orient="index")
    cov_dict = {}
    for ch, po in position_to_mip:
        for m, cp in position_to_mip[(ch, po)]:
            for s in bc_dict:
                try:
                    cov_dict[(s, ch, po)] += bc_dict[s][(m, cp)]
                except KeyError:
                        cov_dict[(s, ch, po)] = bc_dict[s][(m, cp)]

    def return_coverage(k):
        """ Return coverage of a variant position for a sample if the sample
        has any coverage. Return zero if no coverage for that sample.
        """
        try:
            return cov_dict[k]
        except KeyError:
            return 0

    # create a vcf file for all variants
    # create pivot table with each variant having its own column
    vcf_table = variant_counts.pivot_table(
        columns="Sample ID",
        index=["CHROM", "POS", "ID", "REF", "ALT"],
        values="Barcode Count",
        aggfunc="sum",
    )
    # remove place holder sample for unobserved variants
    try:
        vcf_table.drop("Temp", axis=1, inplace=True)
    except KeyError:
        pass
    vcf_table.fillna(0, inplace=True)
    v_cols = vcf_table.columns
    v_index = vcf_table.index
    # Calculate coverage for each variant and position in the vcf table
    vcf_co = pd.DataFrame([
        [return_coverage((s, v[0], v[1]))
         for s in v_cols]
        for v in v_index],
        index=v_index,
        columns=v_cols).fillna(0)
    # merge variants on position to get non-reference allele count.
    # Transforming the groups by the column is extremely slow in python3 for
    # some reason. So we'll transpose the table temporarily.
    vcf_non_ref = vcf_table.groupby(
        level=["CHROM", "POS"],
        axis=0
    ).transform("sum")
    # calculate reference allele counts
    vcf_ref = vcf_co - vcf_non_ref
    # get variant qualities for each variant in each sample.
    variant_counts["Variation Quality"].fillna(-1, inplace=True)
    vcf_quals = variant_counts.pivot_table(
            columns="Sample ID",
            index=["CHROM", "POS", "ID", "REF", "ALT"],
            values="Variation Quality",
            aggfunc="mean"
    ).fillna(-1)
    # convert quality values to string for vcf file.
    vcf_quals = vcf_quals.astype(int).astype(str)
    vcf_quals = vcf_quals.replace("-1", ".")
    try:
        vcf_quals.drop("Temp", axis=1, inplace=True)
    except KeyError:
        pass
    # calculate allele frequencies and create genotype calls from frequencies
    # no filtering will be applied here so even low frequency non-refs mixed
    # with ref will be a HET call. This is only for vcf file to be filtered
    # by proper vcf tools later.
    vcf_freq = vcf_table/vcf_co
    vcf_gen = vcf_freq.applymap(lambda a:
                                "0/0" if a == 0
                                else "." if (np.isnan(a) or np.isinf(a))
                                else "0/1" if a < 1
                                else "1/1")
    # merge all vcf tables to create the merged vcf
    vcf = (vcf_gen + ":" + vcf_ref.astype(int).astype(str)
           + "," + vcf_table.astype(int).astype(str)
           + ":" + vcf_co.astype(int).astype(str)
           + ":" + vcf_quals)
    # Add vcf header
    vcf_samples = vcf.columns.tolist()
    vcf.columns = vcf_samples
    vcf["QUAL"] = "."
    vcf["FILTER"] = "."
    vcf["INFO"] = "."
    vcf["FORMAT"] = "GT:AD:DP:SQ"
    vcf = vcf[["QUAL", "FILTER", "INFO", "FORMAT"] + vcf_samples]
    vcf_header = [
        "##fileformat=VCFv4.2",
        '##ALT=<ID=NON_REF,Description="Represents a possible alternative '
        'allele at this location">"',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
        '"Allelic depths for the ref and alt alleles in that order">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description='
        '"Total read depth (coverage) at this position.>',
        '##FORMAT=<ID=SQ,Number=1,Type=Integer,Description='
        '"Phred scale sequence quality of the variant.">']
    vcf_file = os.path.join(wdir, "variants.vcf")
    with open(vcf_file, "w") as outfile:
        outfile.write("\n".join(vcf_header) + "\n")
        vcf.reset_index().rename(columns={"CHROM": "#CHROM"}).to_csv(
            outfile, index=False, sep="\t"
        )
    # Replace NA values for some fields in the variant counts
    # because they will be used in generating pivot tables
    # they cannot have NA values.
    variant_counts["Gene"].fillna("NA", inplace=True)
    variant_counts["AA Change Position"].fillna("NA", inplace=True)
    variant_counts["ExonicFunc"].fillna("NA", inplace=True)
    # it is possible to filter variants per sample based on their barcode count
    # this can be provided in the settings as minVariantCount
    try:
        min_variant_count = int(settings["minVariantCount"])
    except KeyError:
        min_variant_count = 0
    if ref_resistant:
        variant_counts["Reference Resistant"].fillna("No", inplace=True)
        # create pivot table for each unique variant
        variant_table = variant_counts.pivot_table(
            index="Sample ID",
            columns=["CHROM", "POS", "ID", "REF", "ALT", "Gene",
                     "Mutation Name", "AA Change Position", "ExonicFunc",
                     "Reference Resistant", "Targeted", "Multi Mapping"],
            values="Barcode Count",
            aggfunc="sum"
        )
        # drop the temporary sample place holder, if any
        try:
            variant_table.drop("Temp", inplace=True)
        except KeyError:
            pass
        # if a sample did not have a variant, the table value
        # will be NA. Change those to 0.
        variant_table.fillna(0, inplace=True)
        # Filter based on min count
        variant_table = variant_table.applymap(
            lambda a: a if a >= min_variant_count else 0
        )
        # add amino acid positions and sort table
        # this is to convert an AA Change Position such as Arg59Glu to 59
        # other possible values, for example, 144delGlu for a deletion.
        # we'll just try to get the first number from all possible strings.
        col_list = []
        for c in variant_table.columns:
            pos = c[7].split("-")[-1].split("_")[0]
            # pos is something like Arg59Glu for mutations causing AA changes
            # or "." for intronic or intergenic changes.
            if pos != ".":
                positions = []
                num_found = False
                for dig in pos:
                    try:
                        int(dig)
                        positions.append(dig)
                        num_found = True
                    except ValueError:
                        if num_found:
                            break
                pos = int("".join(positions))
            col_list.append(c[:7] + (pos,) + c[7:])
        column_names = variant_table.columns.names
        new_cols = pd.MultiIndex.from_tuples(
            col_list,
            names=column_names[:7] + ["AA Position"] + column_names[7:]
        )
        variant_table.columns = new_cols
        variant_table.sort_index(level=["Gene", "AA Position"],
                                 axis=1)
        variant_table.columns = variant_table.columns.droplevel(
            level="AA Position"
        )
        # get coverage table with same indexes as
        # the variant table
        v_cols = variant_table.columns
        v_index = variant_table.index
        variant_cov_df = pd.DataFrame(
            [
                [return_coverage((s, c[0], c[1])) for c in v_cols]
                for s in v_index
            ],
            index=v_index,
            columns=v_cols)
        # define nonsynonamous changes to aggregate all non-reference
        # amino acids and get to all reference amino acid calls from there
        nonsyn = list(
            set(variant_counts["ExonicFunc"]).difference(
                [".", "synonymous SNV", "Temp"]
            )
        )
        idx = pd.IndexSlice
        # aggregate all nonsyn calls for each amino acid position
        # this is not ideal, as it ignores indels but this is only
        # used for loci where reference genome is actually mutant
        # that is so far only dhps-437 and there are no common indels
        # in the vicinity. We're also temporarily transpose the
        # dataframe because transform operation is taking much longer
        # on the columns compared to rows.
        variant_table = variant_table.T
        variant_cov_df = variant_cov_df.T
        non_ref_aa_table = variant_table.loc[
            idx[:, :, :, :, :, :, :, :, nonsyn, :, :, :], :
        ].groupby(
            level=["Gene", "AA Change Position"],
            axis=0
        ).transform("sum")
        # non_ref_aa_table loses the synonymous variants in the
        # previous step. We create an all-zero table with
        # variant table indexes by subtracting it from itself
        # than add non_ref table to get the non_ref values
        non_ref_aa_table = (variant_table
                            - variant_table
                            + non_ref_aa_table).fillna(0)
        non_ref_aa_table = non_ref_aa_table.groupby(
            level=["Gene", "AA Change Position"], axis=0
        ).transform(max)
        non_ref_aa_table = non_ref_aa_table.groupby(
            level=[
                "Gene",
                "Mutation Name",
                "Reference Resistant",
                "Targeted",
                "ExonicFunc"],
            axis=0).max()
        # create a like-indexed coverage table
        coverage_aa_table = variant_cov_df.groupby(
            level=[
                "Gene",
                "Mutation Name",
                "Reference Resistant",
                "Targeted",
                "ExonicFunc"],
            axis=0).max()
        # calculate reference amino acid counts
        ref_aa_table = coverage_aa_table - non_ref_aa_table
        # aggregate all variants that lead to the
        # same amino acid change
        mutant_aa_table = variant_table.groupby(
            level=["Gene",
                   "Mutation Name",
                   "Reference Resistant",
                   "Targeted",
                   "ExonicFunc"],
            axis=0).sum()
        # do a sanity check for all the grouping and coverage calculations
        # none of the table values for mutant or reference tables can be
        # larger than the coverage for a given locus.
        if (((mutant_aa_table - coverage_aa_table) > 0).sum().sum()
           + ((ref_aa_table - coverage_aa_table) > 0).sum().sum()) > 0:
            print("Some loci have lower coverage than mutation calls!")
        # Revert transposed dataframes
        variant_table = variant_table.T
        variant_cov_df = variant_cov_df.T
        mutant_aa_table = mutant_aa_table.T
        coverage_aa_table = coverage_aa_table.T
        ref_aa_table = ref_aa_table.T
        # where reference is the variant of interest("Reference Resistant")
        # change mutant count to reference count
        try:
            mutant_aa_table.loc[
                :, idx[:, :, "Yes", :, :],
            ] = ref_aa_table.loc[:, idx[:, :, "Yes", :, :]]
        except KeyError:
            pass
        mutant_aa_table.columns = mutant_aa_table.columns.droplevel(
            "Reference Resistant"
        )
        coverage_aa_table.columns = coverage_aa_table.columns.droplevel(
            "Reference Resistant"
        )
    else:
        # create pivot table for each unique variant
        variant_table = variant_counts.pivot_table(
            index="Sample ID",
            columns=["CHROM", "POS", "ID", "REF", "ALT", "Gene",
                     "Mutation Name", "AA Change Position",
                     "ExonicFunc", "Targeted", "Multi Mapping"],
            values="Barcode Count",
            aggfunc="sum"
        )
        try:
            # drop the temporary sample place holder
            variant_table.drop("Temp", inplace=True)
        except KeyError:
            pass
        # if a sample did not have a variant, the table value will be NA.
        # Change those to 0.
        variant_table.fillna(0, inplace=True)
        # Filter based on min count
        variant_table = variant_table.applymap(
            lambda a: a if a >= min_variant_count else 0
        )
        # get coverage table with same indexes as
        # the variant table
        v_cols = variant_table.columns
        v_index = variant_table.index
        variant_cov_df = pd.DataFrame(
            [
                [return_coverage((s, c[0], c[1])) for c in v_cols]
                for s in v_index
            ],
            index=v_index,
            columns=v_cols)
        # aggregate all variants that lead to the
        # same amino acid change
        variant_table = variant_table.T
        variant_cov_df = variant_cov_df.T
        mutant_aa_table = variant_table.groupby(
            level=[
                "Gene",
                "Mutation Name",
                "Targeted",
                "ExonicFunc"
            ],
            axis=0).sum()
        # create a like-indexed coverage table
        coverage_aa_table = variant_cov_df.groupby(
            level=[
                "Gene",
                "Mutation Name",
                "Targeted",
                "ExonicFunc"],
            axis=0).max()
        variant_table = variant_table.T
        variant_cov_df = variant_cov_df.T
        mutant_aa_table = mutant_aa_table.T
        coverage_aa_table = coverage_aa_table.T
        # do a sanity check for all the grouping and coverage calculations
        # none of the table values for mutant or reference tables can be
        # larger than the coverage for a given locus. However, this is only
        # relevant when the analysis is limited to coding changes.
        if (((mutant_aa_table - coverage_aa_table) > 0).sum().sum()) > 0:
            print(("Some loci have lower coverage than mutation calls!"
                   "This warning is only relevant if the analysis is limited "
                   "to single amino acid changes, excluding indels and "
                   "noncoding changes."))
    variant_counts.to_csv(os.path.join(wdir, "variants.csv"), index=False)
    plot_performance(barcode_counts, wdir=wdir, save=True)
    variant_table.to_csv(os.path.join(wdir, "variant_table.csv"))
    variant_cov_df.to_csv(os.path.join(wdir, "variant_coverage_table.csv"))
    mutant_aa_table.to_csv(os.path.join(wdir, "mutant_table.csv"))
    coverage_aa_table.to_csv(os.path.join(wdir, "mutant_coverage.csv"))

    # Create genotype calls for amino acid changes
    min_mutation_count = int(settings["minMutationCount"])
    min_mutation_fraction = float(settings["minMutationFraction"])
    min_coverage = int(settings["minCoverage"])
    # filter mutants based on the filters from the settings file.
    mutant_aa_table = mutant_aa_table.applymap(
        lambda a: 0 if a < min_mutation_count else a
    )
    coverage_aa_table = coverage_aa_table.applymap(
        lambda a: 0 if a < min_coverage else a
    )
    # Call genotypes from within sample frequency
    mutant_freq_table = (
        mutant_aa_table/coverage_aa_table
    ).replace(np.inf, np.nan)
    mutant_freq_table.to_csv(wdir + "mutation_frequencies.csv")
    genotypes = mutant_freq_table.applymap(
        lambda x: 2 if x >= (1 - min_mutation_fraction)
        else 1 if x >= min_mutation_fraction else np.nan if np.isnan(x) else 0
    )
    genotypes_file = os.path.join(wdir, "genotypes.csv")
    genotypes.to_csv(genotypes_file)
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

    # Create an overview statistics file for samples including
    # total read count, barcode count, and how well they cover each MIP.
    sample_counts = combined_df.groupby("Sample ID")[["Read Count",
                                                      "Barcode Count"]].sum()
    try:
        sample_counts.drop("Temp", inplace=True)
    except KeyError:
        pass
    target_cov = pd.concat(
        [(barcode_counts >= 1).sum(axis=1),
         (barcode_counts >= 5).sum(axis=1),
         (barcode_counts >= 10).sum(axis=1)],
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

    # Create a file with meta data for samples without any data
    no_data = merged_meta.loc[
        ~merged_meta["Sample ID"].isin(sample_counts.index)
    ]
    print(("{} out of {} samples had no data and they were excluded from the"
           " variant calls.").format(no_data.shape[0], merged_meta.shape[0]))
    no_data_file = os.path.join(wdir, "samples_without_data.csv")
    no_data.to_csv(no_data_file)
    make_chrom_vcf(wdir, len(variant_table.columns[0]),
                   min_cov=min_coverage, min_count=min_mutation_count,
                   min_freq=min_mutation_fraction)
    return


###############################################################################
# New contig based analysis for vcf generation
###############################################################################

def get_vcf_haplotypes(settings):
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
    haplotypes_fq_file = wdir + settings["haplotypesFastqFile"]
    haplotypes_sam_file = wdir + settings["haplotypesSamFile"]
    bwa_options = settings["bwaOptions"]
    call_info_file = settings["callInfoDictionary"]
    species = settings["species"]
    try:
        tol = int(settings["alignmentTolerance"])
    except KeyError:
        tol = 200
    # DATA EXTRACTION ###
    # try loading unique haplotypes values. This file is generated
    # in the recent versions of the pipeline but will be missing in older
    # data. If missing, we'll generate it.
    try:
        hap_df = pd.read_csv(wdir + "unique_haplotypes.csv")
    except IOError:
        raw_results = pd.read_table(wdir + settings["mipsterFile"])
        hap_df = raw_results.groupby(
            ["gene_name", "mip_name", "haplotype_ID"])[
            "haplotype_sequence"].first().reset_index()
    # fill in fake sequence quality scores for each haplotype. These scores
    # will be used for mapping only and the real scores for each haplotype
    # for each sample will be added later.
    hap_df["quality"] = hap_df["haplotype_sequence"].apply(
        lambda a: "H" * len(a))
    haps = hap_df.set_index("haplotype_ID").to_dict(orient="index")
    # BWA alignment ####
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
                    call_dict["gene"] = g
                    call_dict["MIP"] = m
                    call_dict["copy"] = c
                    call_dict["mip_number"] = mip_number
                    call_dict["sub_number"] = sub_number
                    call_df_list.append(pd.DataFrame(call_dict, index=[0]))
    call_df = pd.concat(call_df_list)
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

    mapped_haplotypes.to_csv(wdir + "mapped_haplotypes.csv", index=False)
    off_target_haplotypes.to_csv(wdir + "offtarget_haplotypes.csv",
                                 index=False)
    haplotypes.to_csv(wdir + "aligned_haplotypes.csv", index=False)
    haplotype_maps.to_csv(wdir + "all_haplotypes.csv", index=False)
    num_hap = len(set(haplotype_maps["haplotype_ID"]))
    num_off = len(set(off_target_haplotypes["haplotype_ID"]))
    print(("{} of {} haplotypes were off-target, either not mapping to "
           "the reference genome, or best mapping to a region which was "
           "not targeted.").format(num_off, num_hap))
    return


def get_haplotype_counts(settings):
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
    # drop duplicate values originating from
    # multiple sequencing runs of the same libraries
    run_meta = run_meta.drop_duplicates()
    run_meta = run_meta.groupby(
        ["Sample ID", "Library Prep"]
    ).first().reset_index()
    run_meta.to_csv(wdir + "run_meta.csv")
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
    raw_results = pd.read_table(wdir + settings["mipsterFile"])
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
    ac = cc/gc
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
    combined_df = pd.concat([unique_df, multi_df], ignore_index=True)
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
    # Create pivot table of combined barcode counts
    # This is a per MIP per sample barcode count table
    # of the samples with sequencing data
    barcode_counts = pd.pivot_table(combined_df,
                                    index="Sample ID",
                                    columns=["MIP",
                                             "Copy"],
                                    values=["Barcode Count"],
                                    aggfunc=np.sum)
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
    all_barcode_counts.to_csv(os.path.join(wdir, "all_barcode_counts.csv"))
    # Continue working with the barcode counts that does not include the
    # samples which did not have any data.
    barcode_counts.columns = pd.MultiIndex.from_tuples(bc_cols,
                                                       names=["MIP", "Copy"])
    barcode_counts.fillna(0, inplace=True)
    barcode_counts.to_csv(os.path.join(wdir, "barcode_counts.csv"))
    # Create an overview statistics file for samples including
    # total read count, barcode count, and how well they cover each MIP.
    sample_counts = combined_df.groupby("Sample ID")[["Read Count",
                                                      "Barcode Count"]].sum()
    target_cov = pd.concat(
        [(barcode_counts >= 1).sum(axis=1),
         (barcode_counts >= 5).sum(axis=1),
         (barcode_counts >= 10).sum(axis=1)],
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

    # Create a file with meta data for samples without any data
    no_data = run_meta.loc[
        ~run_meta["Sample ID"].isin(sample_counts.index)
    ]
    print(("{} out of {} samples had no data and they were excluded from the"
           " variant calls.").format(no_data.shape[0], run_meta.shape[0]))
    no_data_file = os.path.join(wdir, "samples_without_data.csv")
    no_data.to_csv(no_data_file)
    return


def split_contigs(settings):
    """ Get a haplotypes dict and a call_info dict, align each haplotype to
    reference sequences from the call_info dict."""
    wdir = settings["workingDir"]
    # load mapped haplotypes dataframe
    mapped_haplotypes = pd.read_csv(os.path.join(wdir,
                                                 "mapped_haplotypes.csv"))

    # get all sample IDs to use in variant files.
    sample_ids = pd.read_csv(os.path.join(wdir, "run_meta.csv"))[
        "Sample ID"].tolist()
    # split haplotypes to contigs where any overlapping haplotype will be
    # added to the contig. There is also going to be a buffer of 35 bp, that
    # is, any haplotypes that are closer than 35 bp to a contig will also
    # end up in the contig even if there is no overlap.

    def get_contig(g):
        intervals = zip(g["capture_start"], g["capture_end"])
        return pd.DataFrame(merge_overlap(
            [list(i) for i in intervals], spacer=35))
    contigs = mapped_haplotypes.groupby("Chrom").apply(get_contig)
    contigs = contigs.reset_index()
    contigs.rename(columns={"level_1": "contig", 0: "contig_capture_start",
                            1: "contig_capture_end"}, inplace=True)
    contigs["contig_start"] = contigs["contig_capture_start"] - 20
    contigs["contig_end"] = contigs["contig_capture_end"] + 20

    # get target SNPs if targets are to be included in the variants results
    # even if they are not observed. This is specified by targetJoin parameter
    # in settings
    if int(settings["targetJoin"]):
        targets = pd.read_table("/opt/project_resources/targets.tsv")
        targets["target_length"] = targets["Ref"].apply(len)
        targets["End"] = targets["Pos"] + targets["target_length"] - 1
        targets = targets.merge(contigs)
        targets = targets.loc[(targets["contig_start"] <= targets["Pos"])
                              & (targets["contig_end"] >= targets["End"])]
        targets["capture_start"] = targets["Pos"]
        targets["capture_end"] = targets["End"]
        targets["orientation"] = "forward"
        targets["forward_sequence"] = targets["Alt"]
    # create a contig dictionary in this format:
    # {chromX: {contig#: {contig_start: .., contig_end: ..}}}
    c_dict = contigs.set_index(["Chrom", "contig"]).to_dict(orient="index")
    contig_dict = {}
    for key, value in c_dict.items():
        try:
            contig_dict[key[0]][key[1]] = value
        except KeyError:
            contig_dict[key[0]] = {key[1]: value}
    # assign each haplotype to its contig
    merge_contigs = mapped_haplotypes.merge(contigs)
    mapped_haplotypes = merge_contigs.loc[
        (merge_contigs["contig_capture_start"]
         <= merge_contigs["capture_start"])
        & (merge_contigs["contig_capture_end"]
           >= merge_contigs["capture_end"])]
    species = settings["species"]

    haplotype_counts = pd.read_csv(os.path.join(wdir, "haplotype_counts.csv"))
    contig_list = []
    contigs_dir = os.path.join(wdir, "contigs")
    if not os.path.exists(contigs_dir):
        os.makedirs(contigs_dir)
    gb = mapped_haplotypes.groupby(["Chrom", "contig"])
    contig_info_dict = {}
    for k in gb.groups.keys():
        contig_info = {}
        contig_haplotypes = gb.get_group(k)
        contig_info["chrom"] = k[0]
        contig_info["contig"] = k[1]
        contig_name = contig_info["chrom"] + "_" + str(contig_info["contig"])
        contig_info["contig_name"] = contig_name
        contig_info["contig_start"] = int(
            contig_haplotypes.iloc[0]["contig_start"])
        contig_info["contig_end"] = int(
            contig_haplotypes.iloc[0]["contig_end"])
        contig_info["contigs_dir"] = contigs_dir
        contig_counts = haplotype_counts.merge(
            contig_haplotypes[["haplotype_ID", "Copy"]])
        contig_counts_file = os.path.join(contigs_dir,
                                          contig_name + "_counts.csv")
        contig_info["contig_counts_file"] = contig_counts_file
        contig_counts.to_csv(contig_counts_file, index=False)
        contig_info["species"] = species
        contig_haplotypes_file = os.path.join(contigs_dir,
                                              contig_name + "_haps.csv")
        contig_info["contig_haplotypes_file"] = contig_haplotypes_file
        contig_haplotypes.to_csv(contig_haplotypes_file, index=False)
        contig_info["min_coverage"] = int(settings["minVariantCoverage"])
        contig_info["min_count"] = int(settings["minVariantCount"])
        contig_info["min_wsaf"] = int(settings["minVariantWsaf"])
        contig_info["sample_ids"] = sample_ids
        contig_info["aligner"] = settings["multipleSequenceAligner"]
        contig_info["msa_to_vcf"] = settings["msaToVcf"]
        contig_info["snp_only"] = int(settings["snpOnlyVcf"])
        try:
            contig_info["contig_targets"] = targets.loc[
                (targets["Chrom"] == contig_info["chrom"])
                & (targets["contig"] == contig_info["contig"])]
            if contig_info["contig_targets"].empty:
                contig_info["contig_targets"] = None
        except NameError:
            contig_info["contig_targets"] = None
        contig_list.append(contig_info)
        try:
            contig_info_dict[k[0]][k[1]] = contig_info
        except KeyError:
            contig_info_dict[k[0]] = {k[1]: contig_info}

    with open(os.path.join(wdir, "contig_info.pkl"), "wb") as outfile:
        pickle.dump(contig_info_dict, outfile)
    results = []
    pro = int(settings["processorNumber"])
    p = NoDaemonProcessPool(pro)
    p.map_async(process_contig, contig_list, callback=results.extend)
    p.close()
    p.join()

    with open(os.path.join(wdir, "contig_process_results.json"),
              "w") as outfile:
        json.dump(results, outfile)

    return (contig_info_dict, results)


def merge_contigs(settings, contig_info_dict, results):
    # merge contig vcfs for each chromosome
    wdir = settings["workingDir"]
    species = settings["species"]
    genome_fasta = get_file_locations()[species]["fasta_genome"]
    for chrom in contig_info_dict:
        chrom_vcf_list = os.path.join(wdir, chrom + "_vcf_files.txt")
        chrom_vcf_file = os.path.join(wdir, chrom + ".vcf")
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
        subprocess.call(["bcftools", "concat", "-f", chrom_vcf_list,
                         "-o", chrom_vcf_file])

        split_vcf_file = os.path.join(wdir, chrom + ".split.vcf")
        subprocess.call(["bcftools", "norm", "-m-both", "-N",
                         chrom_vcf_file, "-o", split_vcf_file])

        filt_vcf_file = os.path.join(wdir, chrom + ".split.filt.vcf")

        minVariantBarcodes = settings["minVariantBarcodes"]
        minVariantSamples = settings["minVariantSamples"]
        minVariantSampleFraction = settings["minVariantSampleFraction"]
        minVariantSampleTotal = settings["minVariantSampleTotal"]
        minVariantMeanQuality = settings["minVariantMeanQuality"]
        minVariantMeanWsaf = settings["minVariantMeanWsaf"]
        minMipCountFraction = settings["minMipCountFraction"]

        filter_expressions = [
            "((INFO/AD[1] >= " + minVariantBarcodes + ")",
            "(INFO/AC[1] >= " + minVariantSamples + ")",
            "(INFO/AF[1] >= " + minVariantSampleFraction + ")",
            "(INFO/AN >= " + minVariantSampleTotal + ")",
            "(INFO/QS[1] >= " + minVariantMeanQuality + ")",
            "(INFO/WSAF[1] >= " + minVariantMeanWsaf + ")",
            "(INFO/MCF[1] >= " + minMipCountFraction + ")"]

        filter_expressions = " & ".join(filter_expressions)
        filter_expressions = filter_expressions + ') | (OT !=".")'

        subprocess.call(["bcftools", "view", "-i", filter_expressions,
                         split_vcf_file, "-o", filt_vcf_file])

        merged_vcf_file = os.path.join(wdir, chrom + ".merged.filt.vcf")
        subprocess.call(["bcftools", "norm", "+m-any", "-N",
                         filt_vcf_file, "-o", merged_vcf_file])

        norm_vcf_file = os.path.join(wdir, chrom + ".norm.vcf")
        subprocess.call(["bcftools", "norm", "-m-both", "-f", genome_fasta,
                         filt_vcf_file, "-o", norm_vcf_file])

        # annotate with snpEff
        ann_db = settings["snpEffDb"]
        ann = subprocess.Popen(["java", "-Xmx10g", "-jar",
                                "/opt/species_resources/snpEff/snpEff.jar",
                                ann_db, norm_vcf_file], stdout=subprocess.PIPE)
        annotated_vcf_file = os.path.join(wdir, chrom + ".norm.ann.vcf")
        with open(annotated_vcf_file, "wb") as outfile:
            outfile.write(ann.communicate()[0])


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
            subprocess.call(["muscle", "-in", fasta_file, "-out",
                             alignment_file])
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
            return
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
            (vcf_merge["capture_start"] <= vcf_merge["END"])
            & (vcf_merge["capture_end"] >= vcf_merge["POS"]))
        vcf_merge.loc[~vcf_merge["covered"], "genotype"] = np.nan
        vcf_clean = vcf_merge.loc[~vcf_merge["genotype"].isnull()]
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
                return "."
            genotypes = []
            for i in range(allele_count):
                if (allele_depths[i] >= min_count) and (wsaf[i] >= min_wsaf):
                    genotypes.append(i)
            if len(genotypes) == 0:
                return "."
            else:
                gt = "/".join(map(str, sorted(genotypes)))
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
                             ",".join(map(str, wsaf.round(3)))])

        collapsed_vcf = pd.DataFrame(combined_vcf.groupby(
            ["CHROM", "POS", "REF", "ALT", "Sample ID"]).apply(collapse_vcf)
        ).reset_index()
        vcf_table = collapsed_vcf.pivot_table(
            index=["CHROM", "POS", "REF", "ALT"],
            columns="Sample ID", aggfunc="first")
        vcf_table.fillna(".", inplace=True)

        def get_var_summary(row):
            val = row.values
            ad = []
            quals = []
            wsafs = []
            mip_counts = []
            hap_counts = []
            for v in val:
                if v != ".":
                    ad.append(list(map(int, v.split(":")[1].split(","))))
                    quals.append(v.split(":")[3].split(","))
                    mip_counts.append(list(map(
                        int, v.split(":")[4].split(","))))
                    hap_counts.append(list(map(
                        int, v.split(":")[5].split(","))))
                    wsafs.append(list(map(float, v.split(":")[6].split(","))))
            if len(ad) == 0:
                return "."
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
                "AC=" + ",".join(map(str, (np.array(ad) >= min_count).sum(
                    axis=0))),
                "AF=" + ",".join(map(str, ((np.array(ad) >= min_count).sum(
                    axis=0)/len(ad)).round(5))),
                "AN=" + str(len(ad)),
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
        vcf_table = vcf_table.loc[:, samples].fillna(".")
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
            "##INFO=<ID=QS,Number=R,Type=Float,Description="
            '"Average sequence quality per allele.">',
            "##INFO=<ID=AN,Number=1,Type=Integer,Description="
            '"Number of samples with genotype calls.">',
            "##INFO=<ID=AC,Number=R,Type=Integer,Description="
            '"Number of samples carrying the allele.">',
            "##INFO=<ID=AF,Number=R,Type=Float,Description="
            '"Frequency of samples carrying the allele.">',
            "##INFO=<ID=WSAF,Number=R,Type=Float,Description="
            '"Average nonzero WithinSampleAlleleFrequency.">',
            "##INFO=<ID=MC,Number=R,Type=Float,Description="
            '"Average number of MIPs supporting the allele (when called).">',
            "##INFO=<ID=HC,Number=R,Type=Float,Description="
            '"Average number of haplotypes supporting the allele'
            ' (when called).">',
            "##INFO=<ID=MCF,Number=R,Type=Float,Description="
            '"MC expressed as the fraction of MAX MC.">',
            "##INFO=<ID=OT,Number=.,Type=String,Description="
            '"Variant position overlaps with a target.">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description='
            '"Allelic depths for the ref and alt alleles in that order.">',
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
        print(("Caught exception in worker thread for contig {}.").format(
            contig_name))
        traceback.print_exc()
        print()
        raise e


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
        except Exception:
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
        print("Region list is empty.")
        return
    file_locations = get_file_locations()
    genome_fasta = file_locations[species]["fasta_genome"]
    region_file = "/tmp/regions_" + id_generator(10) + ".txt"
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
    with open(output_file, "w") as outfile:
        outfile.write(get_fasta(region, species))


def get_snps(region, snp_file):
    """ Take a region string and a  tabix'ed snp file,
    return a list of snps which are lists of
    tab delimited information from the snp file. """
    # extract snps using tabix, in tab separated lines
    snp_temp = subprocess.check_output(["tabix", snp_file, region])
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
                                        region, snp_file])
    # split the lines (each SNP)
    snps_split = snp_temp.split("\n")[:-1]
    # add each snp in the region to a list
    # as lists of
    snps = []
    for line in snps_split:
        snp = line.split('\t')[:8]
        snps.append(snp)
    return snps


def targets(must_file, diff_list):
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


def merge_overlap(intervals, spacer=0):
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
                if i == j:
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


def overlap(reg1, reg2, spacer=0):
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
    using separate processors. Read boulder record in given input
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
    primer3_output = subprocess.check_output(
        ["primer3_core",
         "-p3_settings_file="+primer3_settings_DIR+settings,
         primer3_input_DIR + input_file]
    )
    # write boulder record to file.
    with open(primer3_output_DIR + output_file, 'w') as outfile:
        outfile.write(primer3_output.decode("UTF-8"))
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


def primer_parser3(input_file, primer3_output_DIR, bowtie2_input_DIR,
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
        dict_file = open(primer3_output_DIR + parse_out, 'w')
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


def paralog_primer_worker(chores):
    p_name = chores[0]
    p_dic = chores[1]
    p_coord = chores[2]
    p_copies = chores[3]
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
            # check if both ends of the primer has aligned with the reference
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
                p_dic["PARALOG_COORDINATES"][c] = {"ORI": "forward",
                                                   "CHR": chroms[c],
                                                   "NAME": p_name,
                                                   "GENOMIC_START": para_start,
                                                   "GENOMIC_END": para_end,
                                                   "COORDINATES": ref_coord,
                                                   "KEY": para_primer_key}
            else:
                para_primer_key = chroms[c] + ":" + str(para_end) + "-" + str(
                   para_start)
                p_dic["PARALOG_COORDINATES"][c] = {"ORI": "reverse",
                                                   "CHR": chroms[c],
                                                   "NAME": p_name,
                                                   "GENOMIC_START": para_start,
                                                   "GENOMIC_END": para_end,
                                                   "COORDINATES": ref_coord,
                                                   "KEY": para_primer_key}

    return [p_name, p_dic]


def paralog_primers_multi(primer_dict, copies, coordinate_converter, settings,
                          primer3_output_DIR, outname, species, outp=0):
    """ Take a primer dictionary file and add genomic start and end coordinates
    of all its paralogs."""
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
        with open(primer3_output_DIR + outname, "w") as outf:
            json.dump(primer_dict, outf, indent=1)
    return primer_dict


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
                    print(("%s occurs multiple times in fasta file" % header))
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
    return fasta_to_sequence(get_fasta(region, species))


def bowtie2_run(fasta_file, output_file, bowtie2_input_DIR,
                bowtie2_output_DIR, species, process_num=4,
                seed_MM=1, mode="-a", seed_len=18, gbar=1, local=0):
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


def parse_cigar(cigar):
    """ Parse a cigar string which is made up of numbers followed
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


def parse_bowtie(primer_dict, bt_file, primer_out, primer3_output_DIR,
                 bowtie2_output_DIR, species, settings, outp=1):
    """ Take a bowtie output file and filter top N hits per primer.
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
        with open(primer3_output_DIR + primer_out, 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers


def process_bowtie(primers, primer_out, primer3_output_DIR,
                   bowtie2_output_DIR, species, settings, host=False, outp=1):
    """ Take a primer dict with bowtie information added.
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
                                    if (alt_start < 0) or (alt_end < 0):
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
                                if (alt_start < 0) or (alt_end < 0):
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
        with open(primer3_output_DIR + primer_out, 'w') as outfile:
            json.dump(primers, outfile, indent=1)
    return primers


def filter_bowtie(primers, output_file, primer3_output_DIR, species, TM=46,
                  hit_threshold=0, lower_tm=46, lower_hit_threshold=3, outp=1):
    """ Check TMs of bowtie hits of given primers, on a given genome.
    Filter the primers with too many nonspecific hits."""
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
        outfile = open(primer3_output_DIR + output_file, 'w')
        json.dump(primers, outfile, indent=1)
        outfile.close()
    return primers


def alternative(primer_dic, output_file,
                primer3_output_DIR, tm_diff, outp=1):
    """ Pick the best alternative arm for primers that do not bind all
    paralogs. This is done picking the alternative primer with melting
    temperature that is closest to the original primer.
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
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(primer_dic, outfile, indent=1)
    return primer_dic


def score_paralog_primers(primer_dict, output_file, primer3_output_DIR,
                          ext, mask_penalty, species, outp=1):
    """ Score primers in a dictionary according to a scoring matrix.
    Scoring matrices are somewhat crude at this time.
    Arm GC content weighs the most, then arms GC clamp and arm length
    Next_base values are last."""
    primers = primer_dict["primer_information"]
    sequence = primer_dict["sequence_information"]
    extension = (ext == "extension")
    # extract template sequence
    seq_template = sequence["SEQUENCE_TEMPLATE"]
    # find the coordinates of next bases
    for p in primers:
        # get coordinates of primers in the form of "start_base, len"
        coord = primers[p]["COORDINATES"]
        strand = primers[p]["ORI"]
        if strand == "forward":
            primer_end = (int(coord.split(",")[0])
                          + int(coord.split(",")[1]) - 1)
            # 1 is added or subtracted due to sequence index being zero based.
            next_bases = seq_template[(primer_end+1):(primer_end+3)]
        elif strand == "reverse":
            primer_end = (int(coord.split(",")[0])
                          - int(coord.split(",")[1]) + 1)
            next_bases = reverse_complement(seq_template[
                                               (primer_end - 2):primer_end
                                               ])
        # add "NEXT_BASES" key and its value to mip dictionary
        primers[p]["NEXT_BASES"] = next_bases
    # define scoring matrices
    # arm gc content score matrix
    # define best to worst values for gc content.
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
            else:
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
            else:
                arm_gc_con[i] = worst
    # next base score matrix
    # define best to worst values for next bases.
    # This parameter should be important only when comparing mips with equally
    # good gc contents. Therefore, the values are low and does not give a mip
    # a big +.
    best = 10
    mid = 5
    low = 2
    worst = 0
    # define matrix
    next_bases = ({"": worst, "A": worst, "T": worst, "G": worst, "C": worst,
                   "GG": best, "CC": best, "GC": best, "CG": best,
                   "GA": mid, "CA": mid, "GT": mid, "CT": mid,
                   "AG": low, "TG": low, "AC": low, "TC": low,
                   "AA": worst, "AT": worst, "TA": worst, "TT": worst})
    # gc clamp matrix
    # gc clamp will be updated taking into account that although G/C base
    # pairs are the most stable, G/X mismatches are also the most stable.
    # mismatch stability order is G>T>A>C with C being the most discriminating
    # base.
    ext_gc_clamp = {"G": 0, "C": 200, "A": 50, "T": 100}
    lig_gc_clamp = {"G": 0, "C": 200, "A": 50, "T": 100}
    # extension arm lengths score matrix
    # this is added for plasmodium since length is more relaxed
    # to get best possible mips with higher TMs
    # which sometimes leads to very long arms.
    best = 50
    mid = 25
    low = 5
    worst = 0
    extension_len_matrix = {}
    for i in range(18, 36):
        if (i == 18) or (25 <= i <= 28):
            extension_len_matrix[i] = mid
        elif (19 <= i <= 24):
            extension_len_matrix[i] = best
        elif (30 > i > 28):
            extension_len_matrix[i] = low
        elif (i > 29):
            extension_len_matrix[i] = worst
    # ligation arm lengths score matrix
    best = 50
    mid = 25
    low = 10
    worst = 0
    ligation_len_matrix = {}
    for i in range(18, 36):
        if (i == 18) or (i == 19):
            ligation_len_matrix[i] = mid
        elif (20 <= i <= 26):
            ligation_len_matrix[i] = best
        elif (27 <= i <= 30):
            ligation_len_matrix[i] = low
        elif (i > 30):
            ligation_len_matrix[i] = worst
    # score all mips using above matrices
    for p in list(primers.keys()):
        # get arm sequences
        seq = primers[p]["SEQUENCE"]
        # count lower case masked nucleotides
        mask_count = sum(-1 for n in seq if n.islower())
        mask_score = mask_count * mask_penalty
        # arm lengths
        if extension:
            len_score = extension_len_matrix[len(seq)]
        else:
            len_score = ligation_len_matrix[len(seq)]
        # gc clamp
        if extension:
            gc_clamp = ext_gc_clamp[seq[-1].upper()]
        else:
            gc_clamp = lig_gc_clamp[seq[-1].upper()]
        # get gc percentages and convert to int.
        gc = int(float(primers[p]["GC_PERCENT"]))
        # get next base values
        next_b = primers[p]["NEXT_BASES"]
        all_scores = {"arm_len": [len(seq), len_score],
                      "arm_gc": [gc, arm_gc_con[gc]],
                      "mask_penalty": mask_penalty,
                      "gc_clamp": gc_clamp,
                      "next_bases": [next_b, next_bases[next_b.upper()]]}
        # calculate total score
        score = (len_score + arm_gc_con[gc] + mask_score
                 + next_bases[next_b.upper()])
        # add score to dictionary
        primers[p]["SCORE"] = score
        primers[p]["all_scores"] = all_scores
    if outp:
        # write dictionary to json file
        outfile = open(primer3_output_DIR + output_file, "w")
        json.dump(primer_dict, outfile, indent=1)
        outfile.close()
    return primer_dict


def filter_primers(primer_dict, output_file,
                   primer3_output_DIR, n, bin_size, outp=1):
    """ Filter primers so that only top n scoring primers
    ending within the same subregion (determined by bin_size) remain.
    For example, bin_size=3 and n=1 would chose the best scoring primer
    among primers that end within 3 bps of each other."""
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
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(best_primer_dict, outfile, indent=1)
    return best_primer_dict


def pick_paralog_primer_pairs(extension, ligation, output_file,
                              primer3_output_DIR, min_size, max_size,
                              alternative_arms, region_insertions,
                              subregion_name, outp=1):
    """ Pick primer pairs from extension and ligation arm candidate
    dictionary files for a given size range."""
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
                        adjusted_max_size = max((max_size
                                                 - max_insertion_size),
                                                min_size)
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
                            adjusted_max_size = max((max_size
                                                     - max_insertion_size),
                                                    min_size)
                            if (adjusted_max_size
                                    >= pairs[p]["capture_size"] >= 0):
                                captured_copies.append(p)
                        # create a pair name as
                        # PAIR_extension primer number_ligation primer number
                        ext_name = e.split('_')[2]
                        lig_name = l.split('_')[2]
                        pair_name = ("PAIR_" + subregion_name + "_" + ext_name
                                     + "_" + lig_name)
                        if ext_ori:
                            orientation = "forward"
                        else:
                            orientation = "reverse"
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
                                adjusted_max_size = max((max_size
                                                         - max_insertion_size),
                                                        min_size)
                                if (adjusted_max_size
                                        >= alts[a]["capture_size"] >= 0):
                                    captured_copies.append(a)
                                    primer_pairs["pair_information"][
                                        pair_name]["pairs"][a] = alts[a]
                        primer_pairs["pair_information"][pair_name][
                            "alt_copies"] = captured_copies
    # return if no pairs found
    if len(primer_pairs["pair_information"]) == 0:
        print("No primer pairs found.")
        return 1
    # write dict to file in primer_output_DIR
    if outp:
        with open(primer3_output_DIR + output_file, 'w') as outfile:
            json.dump(primer_pairs, outfile, indent=1)
    return primer_pairs


def add_capture_sequence(primer_pairs, output_file, primer3_output_DIR,
                         species, outp=1):
    """
    Extract the sequence between primers using the genome sequence and
    primer coordinates.
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
        with open(primer3_output_DIR + output_file, "w") as outfile:
            json.dump(primer_pairs, outfile, indent=1)
    return primer_pairs


def make_mips(pairs, output_file, primer3_output_DIR, mfold_input_DIR,
              backbone, outp=1):
    """ Make mips from primer pairs by taking reverse complement
    of ligation primer sequence, adding a backbone sequence and
    the extension primer. Standard backbone is used if none
    specified. Add a new key to each primer pair:
    "mip_information" with a dictionary that has SEQUENCE key
    and mip sequence as value."""
    # check if the primer dictionary is empty
    if len(pairs["pair_information"]) == 0:
        print("There are no primer pairs in dictionary")
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
    with open(mfold_input_DIR + output_file, "w") as outfile:
        for primers in pairs["pair_information"]:
            outline = (">" + primers + "\n" + pairs["pair_information"]
                       [primers]["mip_information"]["ref"]['SEQUENCE'] + "\n")
            outfile.write(outline)
    # write mip dictionary to file in primer3_output_DIR
    if outp:
        outfile = open(primer3_output_DIR + output_file, 'w')
        json.dump(pairs, outfile, indent=1)
        outfile.close()
    return pairs


def check_hairpin(pairs, output_file, settings, output_dir, outp=1):
    """ Check possible hiybridization between the MIP arms themselves or
    between the MIP arms and the probe backbone. Remove MIPs with likely
    hairpins.
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
    backbone_tm = float(settings["ligation"]["hairpin_tm"])
    backbone_name = settings["mip"]["backbone"]
    backbone = mip_backbones[backbone_name]
    # go through mips and calculate hairpins
    # we will calculate hairpins by looking at TMs between arm sequences
    # and backbone sequences since the whole MIP sequence is too long
    # for nearest neighbor calculations (at least or primer3 implementation).
    for p in pairs["pair_information"].keys():
        pair_dict = pairs["pair_information"][p]
        mip_dict = pair_dict["mip_information"]
        # for each primer pair we can have a number of mips due to paralog
        # copies having alternative mips. We'll go through each mip.
        for m in mip_dict.keys():
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


def make_chunks(l, n):
    """ Yield successive n-sized chunks from list l.
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
                with open(mfold_input_DIR + input_file +
                                       str(counter),  'w') as outfile:
                    outfile.write("\n".join(outlist))
                temp = [input_file + str(counter), input_file, hairpin_tm,
                                 primer3_output_DIR, mfold_input_DIR, mfold_output_DIR]
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
                except KeyError:
                    continue
        except KeyError:
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
                except KeyError:
                    continue
        except KeyError:
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
        if i< 20:
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
    percent = int(gc_count * 100 / (gc_count + at_count))
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
    except KeyError:
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


def compatible_chains(primer_file, primer3_output_DIR, primer_out, output_file,
                      overlap_same=0, overlap_opposite=0, outp=1):
    try:
        with open(primer3_output_DIR + primer_file, "r") as infile:
            scored_mips = json.load(infile)
    except IOError:
        print("Primer file does not exist.")
        return 1
    else:
        # create in/compatibility lists for each mip
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
            for mip in list(scored_mips["pair_information"].keys()):
                m = scored_mips["pair_information"][mip]
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
                    incompatible.add(mip)
                if next_compat:
                    compatible.add(mip)
            d["incompatible"] = incompatible
            d["compatible"] = compatible

        mip_sets = set()

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
                "compatible"].difference(incomp)
            # if there are mips that can be added, call compatible_recurse
            # function for each of those mips
            if len(comp) > 0:
                for n in comp:
                    compatible_recurse(l + [n])
            # stop recursing when the mip chain cannot be elongated
            else:
                mip_sets.add(frozenset(l))

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
                mip_sets.add(frozenset([k]))
        # the initial mip sets only contain mip chains. We can expand each
        # such set by merging with other sets after removing incompatible
        # mips from the second set.
        set_count = len(mip_sets)
        counter = 0
        expanded_mipset = True
        while((set_count < 10000) and (counter <= 20) and expanded_mipset):
            counter += 1
            new_mip_sets = set()
            expanded_mipset = False
            for s1 in mip_sets:
                inc = set()
                for m in s1:
                    inc.update(scored_mips["pair_information"][m][
                        "incompatible"])
                for s2 in mip_sets.difference(s1):
                    s3 = s2.difference(inc).difference(s1)
                    if len(s3) > 0:
                        new_mip_sets.add(frozenset(s1.union(s3)))
                        expanded_mipset = True
            mip_sets = new_mip_sets
            new_mip_sets = set()
            for s1 in mip_sets:
                for s2 in mip_sets.difference(s1):
                    if s1.issubset(s2):
                        break
                else:
                    new_mip_sets.add(s1)
            mip_sets = new_mip_sets
            set_count = len(mip_sets)
        if outp:
            with open(primer3_output_DIR + output_file, "w") as outfile:
                outfile.write("\n".join([",".join(s) for s in mip_sets])
                              + "\n")
        with open(primer3_output_DIR + primer_out, "wb") as outfile:
            pickle.dump(scored_mips, outfile)
    return mip_sets


def compatible_mips(primer_file, primer3_output_DIR, primer_out, output_file,
               overlap_same = 0, overlap_opposite = 0):
    try:
        with open (primer3_output_DIR + primer_file, "r") as infile:
            scored_mips = json.load(infile)
    except IOError:
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


def write_list(alist, outfile_name):
    """ Convert values of a list to strings and save to file."""
    with open(outfile_name, "w") as outfile:
        outfile.write("\n".join(["\t".join(map(str, l))
                                for l in alist]) + "\n")
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
        Par.run_paralog()
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ",
                   Par.paralog_name))
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
        Par.run_paralog()
        if Par.copies_captured:
            print(("All copies were captured for paralog ", Par.paralog_name))
        else:
            print(("Some copies were NOT captured for paralog ",
                   Par.paralog_name))
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
        fastq_list.append("@" + f)
        fastq_list.append(fasta[f])
        fastq_list.append("+")
        fastq_list.append("H" * len(fasta[f]))
    with open(fastq_file, "w") as outfile:
        outfile.write("\n".join(fastq_list))
    return


def parasight(resource_dir,
              design_info_file,
              designed_gene_list=None,
              extra_extension=".extra"):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
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


def parasight_mod(resource_dir,  design_info_file, species,
                  designed_gene_list=None, extra_extension=".extra",
                  maf=0.1, height=200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    pdf_dir = os.path.join(resource_dir, "mod_pdfs")
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
        basename = os.path.join(design_info[t]["design_dir"], t, t)
        showname = basename + ".show"
        try:
            with open(showname) as infile:
                sln = infile.readlines()[-1].strip().split("\t")
                show_region = sln[0] + ":" + sln[2] + "-" + sln[3]
        except IOError:
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
        except IOError:
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


def parasight_shift(resource_dir, design_info_file, species,
                    designed_gene_list=None, extra_extension=".extra",
                    height=200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    pdf_dir = os.path.join(resource_dir, "mod_pdfs")
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
        basename = os.path.join(design_info[t]["design_dir"], t, t)
        backup_name = basename + ".extra"
        try:
            with open(backup_name) as infile, open(
                    backup_name + ".mod", "w") as outfile:
                thirds = 0
                rd_range = list(range(height))
                for line in infile:
                    newline = line.split("\t")
                    if newline[3] in ["all_targets", "target"]:
                        newline[5] = "-10"
                        newline[6] = "4"
                        outfile.write("\t".join(newline))
                    elif newline[3] in ["capture", "extension", "ligation"]:
                        if (thirds % 3) == 0:
                            rd = random.choice(rd_range)
                        thirds += 1
                        newline[5] = str(-30 - rd)
                        outfile.write("\t".join(newline))
                    else:
                        outfile.write(line)
                outfile.write("\n")
        except IOError:
            continue
        psname = basename + ".01.01.ps"
        pdfname = basename + ".mod.pdf"
        gs_command = ("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                      + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                      "-sOutputFile=" + pdfname + " " + psname)
        if designed_gene_list is not None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".mod.pdf")
        else:
            pdf_convert_list.append(t + ".mod.pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".mod.pdf "
                        + pdf_dir + t + ".mod.pdf")
        outlist = ["parasight76.pl",
                   "-showseq", basename + ".show",
                   "-extra", basename + extra_extension + ".mod",
                   "-template", "/opt/resources/nolabel.pst",
                   "-precode file:" + basename + ".precode", "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode", "w") as outfile:
            outfile.write("$opt{'filename'}='" + t
                          + "';&fitlongestline; &print_all (0,'"
                          + basename + "')")
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


def parasight_print(gene_list, extra_suffix=".extra"):
    for g in gene_list:
        print(("cd ../" + g))
        print(("parasight76.pl -showseq " + g + ".show "
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


def get_data(settings_file):
    settings = get_analysis_settings(settings_file)
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


def combine_sample_data(gr):
    """ Combine data from multiple sequencing runs for the same sample.

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
    run_meta_collapsed = run_meta.groupby(
        ["sample_name", "capital_set", "replicate", "Library Prep"]
    ).first().reset_index()[["sample_name", "capital_set",
                             "replicate", "Library Prep"]]
    run_meta_collapsed["new_replicate"] = run_meta_collapsed.groupby(
        "sample_name"
    )["replicate"].transform(
        lambda g: list(map(str, list(range(1, len(g) + 1))))
    )
    run_meta = run_meta.merge(run_meta_collapsed)
    run_meta["Sample ID"] = run_meta[["sample_name",
                                      "capital_set",
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
    info.to_csv(wdir + combined_file, index=False, sep="\t")
    info.groupby(["gene_name", "mip_name", "haplotype_ID"])[
        "haplotype_sequence"].first().reset_index().to_csv(
            wdir + "unique_haplotypes.csv", index=False)
    run_meta = run_meta.groupby("Sample ID").first().reset_index()
    run_meta = run_meta.drop(["Sample ID",
                              "sample_set",
                              "sheet_order",
                              "replicate"],
                             axis=1).rename(
        columns={"capital_set": "sample_set",
                 "new_replicate": "replicate"}
    )
    run_meta.to_csv(wdir + "samples.tsv", sep="\t", index=False)


def update_probe_sets(mipset_table = "/opt/resources/mip_ids/mipsets.csv",
                     mipset_json = "/opt/resources/mip_ids/probe_sets.json"):
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
    Generate fastq files for each sample. These files will have stitched and
    barcode corrected reads.
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


def convert_to_int(n):
    """
    Convert values to integers. This is to be used when
    a pandas dataframe converts integers to floats due
    to the presence of NA values and integer values are
    preferred over floats, i.e. string conversion/comparison.
    """
    try:
        return int(n)
    except ValueError:
        return np.nan


def get_ternary_genotype(gen):
    """
    Convert a 0/0, 0/1, 1/1 type genotype string to
    0, 1, 2.
    """
    try:
        g = sum(map(int, gen.split(":")[0].split("/")))
    except ValueError:
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
                    except KeyError:
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
                except KeyError:
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
    except IndexError:
        return "."
    except AttributeError:
        return np.nan
def split_aa_pos(aa):
    try:
        return aa.split(";")[0].split(":")[4][2:-1]
    except IndexError:
        return "."
    except AttributeError:
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
                            except KeyError:
                                ms_types[ms_len] = bcc
                        for ml in list(ms_types.keys()):
                            if (ms_types[ml]/total_bcs) < freq_cutoff:
                                ms_types.pop(ml)
                        try:
                            sam_freq[g][m][c] = ms_types
                        except KeyError:
                            try:
                                sam_freq[g][m] = {c : ms_types}
                            except KeyError:
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
    except KeyError:
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
        &(data_summary["total_barcode_count"] < barcode_count_threshold),
        "Status"
    ] = low_coverage_action
    # mark samples with too high barcode coverage
    # these samples will have been sequenced to a high depth but
    # low barcode numbers, so sequencing these more would not make sense.
    # They will be re-captured if more data is needed.
    try:
        data_summary["Barcode Coverage"]
    except KeyError:
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
    except TypeError:
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
        except KeyError:
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
                                # add the aa change position to the changes
                                # dict.
                                aa_changes[aa_pos].append(i)
                            except KeyError:
                                aa_changes[aa_pos] = [i]
                        except (IndexError, ValueError):
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
                                except ValueError:
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
                                'ExonicFunc': ExonicFunc,
                                'Func.refGene': 'exonic',
                                'GeneID': d["annotation"]['GeneID'],
                                'GeneDetail.refGene': d["annotation"]['GeneDetail.refGene'],
                                'Otherinfo': d["annotation"]["Otherinfo"],
                                'Ref': Ref,
                                'Start': g_start
                            },
                                           'begin': g_start,
                                            'chrom': d["chrom"],
                                            'end': g_end,
                                            'hap_base': hap_codon,
                                            'hap_index': [h_start_index, h_end_index - 1],
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
                            except KeyError:
                                var_info.append(np.nan)
                            except IndexError:
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
                                           ref_base[i], "-",
                                          qual, filt]
                                var_info = []
                                for col in info_cols:
                                    try:
                                        var_info.append(info_dict[col][0])
                                    except KeyError:
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


def iupac_fasta_converter(header, sequence):
    """
    Given a sequence (header and sequence itself) containing iupac characters,
    return a dictionary with all possible sequences converted to ATCG.
    """
    iupac_dict = {"R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
                  "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
                  "N": "ACGT"}
    iupac_dict = {k: list(iupac_dict[k])
                  for k in list(iupac_dict.keys())}
    if sequence.upper().count("N") >= 10:
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
    """ Save a fasta dictionary to file. """
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
    """
    Create a sample sheet to be used by bcl2fasq file from sample list.
    """
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

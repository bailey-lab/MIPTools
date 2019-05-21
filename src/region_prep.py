import sys
import mip_functions as mip
import mip_classes as mod
import pickle
import json
import copy
import os
import subprocess
import pandas as pd
import argparse

# Read input arguments
parser = argparse.ArgumentParser(description=""" Run MIP design pipeline.""")
parser.add_argument("-n", "--design-name",
                    help="A short name for the current design.",
                    required=True)
parser.add_argument("-s", "--species",
                    help="Target species name.",
                    required=True)
parser.add_argument("-o", "--host-species",
                    help="Host species name.",
                    Default=None)
parser.add_argument("-p", "--processor-number",
                    help="Number of available processors.",
                    type=int,
                    Default=7)
parser.add_argument("--parallel-processes",
                    help="Number of designs  carried out in parallel.",
                    type=int,
                    Default=1)
parser.add_argument("--flank",
                    help="Number of bases to flank target sites on each side.",
                    type=int,
                    Default=250)
parser.add_argument("--single-mip-threshold",
                    help=("Target regions smaller than this will have a "
                          "single MIP designed."),
                    type=int,
                    Default=0)
parser.add_argument("--min-target-size",
                    help="Length threshold to eliminate un-aligned targets.",
                    type=int,
                    Default=150)
parser.add_argument("--max-indel-size",
                    help=("Indel size limit between paralogs, allowing this"
                          "size gapped alignment within paralogs."),
                    type=int,
                    Default=50)
parser.add_argument("--coordinates-file",
                    help="Path to file containing target coordinates.",
                    Default=None)
parser.add_argument("--genes-file",
                    help="Path to file containing target genes.",
                    Default=None)
parser.add_argument("--snps-file",
                    help="Path to file containing SNP coordinates.",
                    Default=None)
parser.add_argument("--fasta-files",
                    help=("Path(s) to fasta file(s) containing target "
                          "sequences."),
                    nargs="*")
parser.add_argument("--fasta-capture-type",
                    help="Capture type for targets specified in fasta files.",
                    Default="whole")
parser.add_argument("--genome-coverage",
                    help="Minimum alignment length to reference genome.",
                    type=int,
                    Default=None)
parser.add_argument("--genome-identity",
                    type=int,
                    help=("Minimum sequence identity (0-100) for genomic "
                          "alignments."),
                    Default=100)
parser.add_argument("--local-coverage",
                    type=int,
                    help="Minimum alignment length to reference paralog.",
                    Default=None)
parser.add_argument("--local-identity",
                    type=int,
                    help=("Minimum sequence identity (0-100) for paralog "
                          "alignments."),
                    Default=100)
parser.add_argument("--design-dir",
                    help="Path to location to output design files.",
                    Default="/opt/analysis")
parser.add_argument("--resource-dir",
                    help="Path to location to output prep files.",
                    Default="/opt/project_resources")
parser.add_argument("--targets-rinfo-template",
                    help="Path to rinfo template for 'targets' capture type.",
                    Default=("/opt/project_resources/rinfo_templates/"
                             "targets_template.rinfo"))
parser.add_argument("--exons-rinfo-template",
                    help="Path to rinfo template for 'exons' capture type.",
                    Default=("/opt/project_resources/rinfo_templates/"
                             "exonic_template.rinfo"))
parser.add_argument("--whole-rinfo-template",
                    help="Path to rinfo template for 'whole' capture type.",
                    Default=("/opt/project_resources/rinfo_templates/"
                             "whole_template.rinfo"))
parser.add_argument("--output-file",
                    help="Base name to save region prep results.",
                    Default=("/opt/project_resources/design_info"))
# parse arguments from command line
args = vars(parser.parse_args())
design_dir = args["design_dir"]
res_dir = args["resource_dir"]
design_name = args["design_name"]
species = args["species"]
host_species = args["host_species"]
num_process = args["processor_number"]
parallel_processes = args["parallel_processes"]
processor_per_region = int(num_process/parallel_processes)
flank = args["flank"]  # set flank to 150 for maximum capture size of 170
single_mip_threshold = args["single_mip_threshold"]
genes_file = args["genes_file"]
snps_file = args["snps_file"]
coordinates_file = args["coordinates_file"]
fasta_files = args["fasta_files"]
fasta_capture_type = args["fasta_capture_type"]
max_allowed_indel_size = args["max_indel_size"]
min_target_size = args["min_target_size"]
genome_coverage = args["genome_coverage"]
genome_identity = args["genome_identity"]
intra_coverage = args["local_coverage"]
intra_identity = args["local_identity"]
targets_rinfo_template = args["targets_rinfo_template"]
exonic_rinfo_template = args["exonic_rinfo_template"]
whole_rinfo_template = args["whole_rinfo_template"]
output_file = args["output_file"]
# extract target coordinates from target files
target_coordinates = mip.get_target_coordinates(
    res_dir, species, capture_size=flank, coordinates_file=coordinates_file,
    snps_file=snps_file, genes_file=genes_file)
target_regions = target_coordinates["target_regions"]
capture_types = target_coordinates["capture_types"]
gene_names = target_coordinates["gene_names"]
target_names = target_coordinates["target_names"]

# align targets to reference genome
target_alignments = mip.align_targets(
    res_dir, target_regions, species, flank, fasta_files, fasta_capture_type,
    genome_identity, genome_coverage, num_process, gene_names,
    max_allowed_indel_size, intra_identity, intra_coverage, capture_types,
    min_target_size)

# get alignment results. Run alignments again in case there are remaining
# unaligned targets.
final_target_regions = target_alignments["target_regions"]
final_region_names = target_alignments["region_names"]
imperfect_aligners = target_alignments["imperfect_aligners"]
overlaps = target_alignments["overlaps"]
final_capture_types = target_alignments["capture_types"]
target_regions = target_alignments["missed_target_regions"]
capture_types = target_alignments["missed_capture_types"]
extra_target_alignments = mip.align_targets(
    res_dir, target_regions, species, flank, [], fasta_capture_type,
    genome_identity, genome_coverage, num_process, gene_names,
    max_allowed_indel_size, intra_identity, intra_coverage, capture_types,
    min_target_size)

# update initial alignment results with the new alignment results
final_target_regions.update(target_alignments["target_regions"])
final_region_names.update(target_alignments["region_names"])
imperfect_aligners.extend(target_alignments["imperfect_aligners"])
overlaps.update(target_alignments["overlaps"])
unaligned_target_regions = target_alignments["missed_target_regions"]
if len(unaligned_target_regions) > 0:
    unaligned_targets_file = os.path.join(res_dir, "unaligned_targets.json")
    print(("{} regions remain un-aligned after genomic alignments."
           " See {} and decide if an additional alignment is needed.").format(
           len(unaligned_target_regions), unaligned_targets_file))
    with open(unaligned_targets_file, "w") as outfile:
        json.dump(unaligned_target_regions, outfile, indent=1)

# perform cross mapping of coordinates between paralog copies
alignments = {}
for t in final_target_regions:
    alignments[t] = mip.alignment_mapper(t, res_dir)

# check if any SNP target provided has failed the region prep prpocess
snp_coordinates = target_coordinates["snp_coordinates"]
snps_in_targets = {}
snps_missed = list(snp_coordinates.keys())
for s in snp_coordinates:
    snp = snp_coordinates[s]
    snp_begin = snp["begin"]
    snp_end = snp["end"]
    snp_chrom = snp["chrom"]
    snp["name"] = s
    snp["size_difference"] = snp["InsertionLength"]
    for t in final_target_regions:
        for i in range(len(final_target_regions[t])):
            r = final_target_regions[t][i]
            region_chrom = r[0]
            region_start = r[1]
            region_end = r[2]
            if region_chrom == snp_chrom and region_start <= \
                    snp_begin <= snp_end <= region_end:
                snp_copy = snp["copy"] = "C" + str(i)
                try:
                    snps_missed.remove(s)
                except ValueError:
                    pass
                try:
                    snps_in_targets[t][s] = snp
                except KeyError:
                    snps_in_targets[t] = {s: snp}
if len(snps_missed) > 0:
    print(("Some SNP targets are missing from target regions. {}").format(
        snps_missed))

# find which of the original target is now part of each final target group
included_targets = {}
for t in final_target_regions:
    overlapping_region_list = overlaps[t]
    for o in overlapping_region_list:
        try:
            regions_included = target_names[o]
        except KeyError:
            regions_included = [o]
        try:
            included_targets[t].extend(regions_included)
        except KeyError:
            included_targets[t] = regions_included

# set capture template files for each capture type for all regions
capture_priority = ["whole", "exons", "targets"]
capture_rinfos = {"whole": whole_rinfo_template,
                  "exons": exonic_rinfo_template,
                  "targets": targets_rinfo_template}
region_capture_types = {}
for t in final_target_regions:
    inc = included_targets[t]
    included_capture_types = []
    for i in inc:
        included_capture_types.append(final_capture_types[i])
    for c in capture_priority:
        if c in included_capture_types:
            region_capture_types[t] = c
            break
    else:
        print(("Capture type not found for {}. "
               "Setting capture type to 'whole'").format(t))

# create a rinfo file for each target using the rinfo template
# save the target_dict to be used in the design module
for g in final_target_regions:
    target_dict = {"regions": final_target_regions[g],
                   "alignments": alignments[g],
                   "genes": overlaps[g],
                   "copy_names": final_region_names[g],
                   "included_targets": included_targets[g]}
    single_mip = True
    for fr in final_target_regions[g]:
        if fr[-1] > single_mip_threshold:
            single_mip = False
    cap_type = region_capture_types[g]
    rinfo_file = capture_rinfos[cap_type]
    if not os.path.exists(design_dir + g + "/resources/"):
        os.makedirs(design_dir + g + "/resources/")
    with open(design_dir + g + "/resources/" + g + ".pkl", "wb") as outfile:
        pickle.dump(target_dict, outfile)
    with open(rinfo_file) as infile:
        outfile = design_dir + g + "/resources/" + g + ".rinfo"
        outfile_list = []
        for line in infile:
            newline = line.strip().split("\t")
            if line.startswith("NAME"):
                newline[1] = g
                outfile_list.append(newline)
            elif line.startswith("DESIGN_DIR"):
                newline[1] = design_dir
                outfile_list.append(newline)
            elif line.startswith("SPECIES"):
                newline[1] = species
                outfile_list.append(newline)
            elif line.startswith("HOST_SPECIES"):
                newline[1] = host_species
                outfile_list.append(newline)
            elif (line.startswith("CAPTURE")
                  and (newline[1] == "max_mips")
                  and single_mip):
                newline[2] = 1
                outfile_list.append(newline)
            elif (line.startswith("SETTINGS")
                  and (newline[1] == "processors")):
                newline[2] = newline[3] = newline[4] = processor_per_region
                outfile_list.append(newline)
            else:
                outfile_list.append(newline)
        try:
            target_snps = snps_in_targets[g].keys()
            rinfo_keys = ["name", "copy", "chrom", "begin", "end",
                          "size_difference"]
            for rk in rinfo_keys:
                outlist = ["MUST", rk]
                for s in target_snps:
                    outlist.append(snps_in_targets[g][s][rk])
                outfile_list.append(outlist)
        except KeyError:
            pass
        mip.write_list(outfile_list, outfile)

# save a dictionary containing all region prep data
design_info = {}
design_info = {"resource_dir": res_dir,
               "design_dir": design_dir,
               "final_target_regions": final_target_regions,
               "final_region_names": final_region_names,
               "included_targets": included_targets,
               "alignments": alignments,
               "overlaps": overlaps,
               "region_capture_types": region_capture_types}
with open(output_file + ".pkld", "wb") as outfile:
    pickle.dump(design_info, outfile)

# create a summary of target regions for visual evaluation
target_info = {}
for g in final_target_regions:
    target_info[g] = {"resource_dir": res_dir,
                      "design_dir": design_dir,
                      "final_target_regions": final_target_regions[g],
                      "included_targets": included_targets[g],
                      "region_capture_types": region_capture_types[g]}

with open(output_file + ".json", "w") as outfile:
    json.dump(target_info, outfile, indent=1)

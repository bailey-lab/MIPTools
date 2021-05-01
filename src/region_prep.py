"""Prepare target genomic regions for MIP design."""
import mip_functions as mip
import pickle
import json
import os
import argparse
import subprocess
import pandas as pd


def region_prep(design_name, species, host_species=None, processor_number=7,
                parallel_processes=1, flank=150, merge_distance=None,
                design_dir="/opt/analysis",
                resource_dir="/opt/project_resources",
                single_mip_threshold=0, genes_file=None, snps_file=None,
                coordinates_file=None, fasta_files=[],
                fasta_capture_type="whole",
                max_indel_size=50, min_target_size=150, genome_coverage=250,
                genome_identity=100, local_coverage=100, local_identity=100,
                targets_rinfo_template=("/opt/resources/templates/"
                                        "rinfo_templates/template_rinfo.txt"),
                exons_rinfo_template=("/opt/resources/templates/"
                                      "rinfo_templates/template_rinfo.txt"),
                whole_rinfo_template=("/opt/resources/templates/"
                                      "rinfo_templates/template_rinfo.txt"),
                output_file="/opt/project_resources/design_info"):
    """Prepare target genomic regions for MIP design."""
    if host_species is None:
        host_species = "none"

    processor_per_region = int(processor_number/parallel_processes)

    if merge_distance is None:
        merge_distance = flank
    else:
        merge_distance = int(merge_distance/2)

    # extract target coordinates from target files
    target_coordinates = mip.get_target_coordinates(
        resource_dir, species, capture_size=merge_distance,
        coordinates_file=coordinates_file, snps_file=snps_file,
        genes_file=genes_file)
    target_regions = target_coordinates["target_regions"]
    capture_types = target_coordinates["capture_types"]
    gene_names = target_coordinates["gene_names"]
    target_names = target_coordinates["target_names"]

    # align targets to reference genome
    target_alignments = mip.align_targets(
        resource_dir, target_regions, species, flank, fasta_files,
        fasta_capture_type, genome_identity, genome_coverage, processor_number,
        gene_names, max_indel_size, local_identity, local_coverage,
        capture_types, min_target_size, merge_distance,
        "alignment_results.json")

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
        resource_dir, target_regions, species, flank, [], fasta_capture_type,
        genome_identity, genome_coverage, processor_number, gene_names,
        max_indel_size, local_identity, local_coverage, capture_types,
        min_target_size, merge_distance, "extra_alignment_results.json")

    # update initial alignment results with the new alignment results
    final_target_regions.update(extra_target_alignments["target_regions"])
    final_region_names.update(extra_target_alignments["region_names"])
    imperfect_aligners.extend(extra_target_alignments["imperfect_aligners"])
    overlaps.update(extra_target_alignments["overlaps"])
    final_capture_types.update(extra_target_alignments["capture_types"])
    unaligned_target_regions = extra_target_alignments["missed_target_regions"]
    if len(unaligned_target_regions) > 0:
        unaligned_targets_file = os.path.join(resource_dir,
                                              "unaligned_targets.json")
        print(("{} regions remain un-aligned after genomic alignments."
               " See {} and decide if an additional alignment is "
               "needed.").format(
               len(unaligned_target_regions), unaligned_targets_file))
        with open(unaligned_targets_file, "w") as outfile:
            json.dump(unaligned_target_regions, outfile, indent=1)

    # perform cross mapping of coordinates between paralog copies
    alignments = {}
    for t in final_target_regions:
        alignments[t] = mip.alignment_mapper(t, resource_dir)

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
        snp["size_difference"] = snp["MaxInsertionLength"]
        for t in final_target_regions:
            for i in range(len(final_target_regions[t])):
                r = final_target_regions[t][i]
                region_chrom = r[0]
                region_start = r[1]
                region_end = r[2]
                if region_chrom == snp_chrom and region_start <= \
                        snp_begin <= snp_end <= region_end:
                    snp["copy"] = "C" + str(i)
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
                      "exons": exons_rinfo_template,
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
            region_capture_types[t] = "whole"

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
        r_dir = os.path.join(design_dir, g, "resources")
        if not os.path.exists(r_dir):
            os.makedirs(r_dir)
        with open(os.path.join(r_dir, g + ".pkl"), "wb") as outfile:
            pickle.dump(target_dict, outfile)
        with open(rinfo_file) as infile:
            outfile = os.path.join(r_dir, g + ".rinfo")
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
                elif (line.startswith("SELECTION")
                      and (newline[1] == "chain")):
                    if cap_type == "targets":
                        newline[2] = 0
                    else:
                        newline[2] = 1
                    outfile_list.append(newline)
                elif (line.startswith("CAPTURE")
                      and (newline[1] == "single_mip")
                      and single_mip):
                    newline[2] = 1
                    outfile_list.append(newline)
                elif (line.startswith("CAPTURE")
                      and (newline[1] == "flank")):
                    newline[2] = flank
                    outfile_list.append(newline)
                elif (line.startswith("CAPTURE")
                      and (newline[1] == "capture_type")):
                    newline[2] = cap_type
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
        # copy primer settings files
        primer_templates_dir = "/opt/resources/primer3_settings"
        ext_primer_template = os.path.join(primer_templates_dir,
                                           "extension_primer_settings.txt")
        lig_primer_template = os.path.join(primer_templates_dir,
                                           "ligation_primer_settings.txt")
        subprocess.call(["rsync", ext_primer_template,
                         lig_primer_template, r_dir])
    # save a dictionary containing all region prep data
    design_info = {}
    design_info = {"resource_dir": resource_dir,
                   "design_dir": design_dir,
                   "final_target_regions": final_target_regions,
                   "final_region_names": final_region_names,
                   "included_targets": included_targets,
                   "alignments": alignments,
                   "overlaps": overlaps,
                   "region_capture_types": region_capture_types,
                   "parallel_processes": parallel_processes}
    with open(output_file + ".pkld", "wb") as outfile:
        pickle.dump(design_info, outfile)

    # create a summary of target regions for visual evaluation
    target_info = {}
    for g in final_target_regions:
        target_info[g] = {"resource_dir": resource_dir,
                          "design_dir": design_dir,
                          "final_target_regions": final_target_regions[g],
                          "included_targets": included_targets[g],
                          "region_capture_types": region_capture_types[g],
                          "parallel_processes": parallel_processes}

    with open(output_file + ".json", "w") as outfile:
        json.dump(target_info, outfile, indent=1)

    # create a table for visual evaluation
    tr_list = []
    for g in target_info:
        tinfo = target_info[g]
        inc = tinfo["included_targets"]
        tr_regs = tinfo["final_target_regions"]
        cap_type = tinfo["region_capture_types"]
        for i in range(len(tr_regs)):
            temp_list = [g]
            temp_list.append("C" + str(i))
            temp_list.extend(tr_regs[i][:3])
            temp_list.append(",".join(inc))
            temp_list.append(cap_type)
            tr_list.append(temp_list)
    tr_df = pd.DataFrame(tr_list, columns=["Name", "Copy", "Chrom", "Start",
                                           "End", "IncludedTargets",
                                           "CaptureType"])
    tr_df.to_csv(output_file + ".tsv", sep="\t", index=False)


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser(
        description=""" Run MIP design pipeline.""")
    parser.add_argument("-n", "--design-name",
                        help="A short name for the current design.",
                        required=True)
    parser.add_argument("-s", "--species",
                        help="Target species name.",
                        required=True)
    parser.add_argument("-o", "--host-species",
                        help="Host species name.",
                        default=None)
    parser.add_argument("-p", "--processor-number",
                        help="Number of available processors.",
                        type=int,
                        default=7)
    parser.add_argument("--parallel-processes",
                        help="Number of designs  carried out in parallel.",
                        type=int,
                        default=1)
    parser.add_argument("--flank",
                        help=("Number of bases to flank target sites on "
                              "each side."),
                        type=int,
                        default=150)
    parser.add_argument("--merge-distance",
                        help=("Targets closer than this will be placed in the "
                              "same target region. Same distance will be used "
                              "for merging aligned targets. If no argument is "
                              "given, it will default to 'flank' value."),
                        type=int,
                        default=None)
    parser.add_argument("--single-mip-threshold",
                        help=("Target regions smaller than this will have a "
                              "single MIP designed."),
                        type=int,
                        default=0)
    parser.add_argument("--min-target-size",
                        help=("Length threshold to eliminate un-aligned "
                              "targets."),
                        type=int,
                        default=150)
    parser.add_argument("--max-indel-size",
                        help=("Indel size limit between paralogs, allowing "
                              "this size gapped alignment within paralogs."),
                        type=int,
                        default=50)
    parser.add_argument("--coordinates-file",
                        help="Path to file containing target coordinates.",
                        default=None)
    parser.add_argument("--genes-file",
                        help="Path to file containing target genes.",
                        default=None)
    parser.add_argument("--snps-file",
                        help="Path to file containing SNP coordinates.",
                        default=None)
    parser.add_argument("--fasta-files",
                        help=("Path(s) to fasta file(s) containing target "
                              "sequences."),
                        nargs="*")
    parser.add_argument("--fasta-capture-type",
                        help=("Capture type for targets specified "
                              "in fasta files."),
                        default="whole")
    parser.add_argument("--genome-coverage",
                        help="Minimum alignment length to reference genome.",
                        type=int,
                        default=1000)
    parser.add_argument("--genome-identity",
                        type=int,
                        help=("Minimum sequence identity (0-100) for genomic "
                              "alignments."),
                        default=100)
    parser.add_argument("--local-coverage",
                        type=int,
                        help="Minimum alignment length to reference paralog.",
                        default=100)
    parser.add_argument("--local-identity",
                        type=int,
                        help=("Minimum sequence identity (0-100) for paralog "
                              "alignments."),
                        default=100)
    parser.add_argument("--design-dir",
                        help="Path to location to output design files.",
                        default="/opt/analysis")
    parser.add_argument("--resource-dir",
                        help="Path to location to output prep files.",
                        default="/opt/project_resources")
    parser.add_argument("--targets-rinfo-template",
                        help=("Path to rinfo template for 'targets' "
                              "capture type."),
                        default=("/opt/resources/templates/rinfo_templates/"
                                 "template_rinfo.txt"))
    parser.add_argument("--exons-rinfo-template",
                        help=("Path to rinfo template for 'exons' "
                              "capture type."),
                        default=("/opt/resources/templates/rinfo_templates/"
                                 "template_rinfo.txt"))
    parser.add_argument("--whole-rinfo-template",
                        help=("Path to rinfo template for 'whole' "
                              "capture type."),
                        default=("/opt/resources/templates/rinfo_templates/"
                                 "template_rinfo.txt"))
    parser.add_argument("--output-file",
                        help="Base name to save region prep results.",
                        default=("/opt/project_resources/design_info"))

    # parse arguments from command line
    args = vars(parser.parse_args())
    # run region_prep function with command line arguments
    region_prep(**args)

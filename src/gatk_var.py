import argparse
import json
import subprocess


# Read input arguments
parser = argparse.ArgumentParser(
    description="""Call variants using GATK pipeline.""")

parser.add_argument("-d", "--design-info",
                    help=("Path to design info file created by the region "
                          "prep module."),
                    default="/opt/project_resources/design_info.json")

parser.add_argument("-r", "--resource-dir",
                    help=("Path to directory where output files are saved."),
                    default="/opt/project_resources/")

parser.add_argument("-f", "--filter-type",
                    help=(("If a selected MIP file is provided, how should "
                           "the listed MIPs be handled.")),
                    default="remove",
                    choices=["remove", "keep"])

parser.add_argument("-n", "--design-name",
                    help=("A name given this MIP design project."),
                    required=True)

args = vars(parser.parse_args())
genome_dict_file =
call_info_file

def create_intervals(call_info_file, interval_bed_file, interval_list_file,
                     padding=200):
    with open(call_info_file) as infile:
        call_info = json.load(infile)
    with open(interval_bed_file, "w") as outfile:
        for g in call_info:
            for m in call_info[g]:
                for c in call_info[g][m]["copies"]:
                    copy_dict = call_info[g][m]["copies"][c]
                    capture_start = copy_dict["capture_start"] - padding
                    if capture_start < 0:
                        capture_start = 0
                    capture_end = str(int(copy_dict["capture_end"] + padding))
                    capture_start = str(int(capture_start))
                    outfile.write("\t".join([copy_dict["chrom"],
                                             capture_start, capture_end]
                                            ) + "\n")

    subprocess.call(["gatk", "BedToIntervalList", "-I", interval_bed_file,
                     "-O", interval_list_file, "-SD", genome_dict_file])


def haplotype_caller_worker(l):
    bam = l["bam"]
    ref = l["reference"]
    jav = l["java_options"]

import mip_functions as mip
import pandas as pd
import argparse
import os
import subprocess

# Read input arguments
parser = argparse.ArgumentParser(
    description="""Call variants using freebayes.""")
parser.add_argument("-f", "--fastq",
                    help="Set this flag to create fastq files from MIP data.",
                    action="store_true")
parser.add_argument("-d", "--fastq")
parser.add_argument("-a", "--align",
                    help="Set this flag to align fasta files to genome.",
                    action="store_true")
parser.add_argument("-e", "--extra-options",
                    help=("Additional freebayes options to pass directly "
                          "to freebayes. Options must have + in place of -. "
                          "For example, ++pooled+continuous if you want to "
                          "pass --pooled-continuous"),
                    action="append",
                    default=[])

# parse arguments from command line
args = vars(parser.parse_args())

gatk = args["gatk"]

if freebayes:
freebayes_options = args["freebayes_options"]
map_haplotypes = args["map_haplotypes"]

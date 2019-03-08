import sys
sys.path.append("/opt/src/")
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
parser.add_argument("-o", "--host",
                    help="Host species name.",
                    Default=None)
parser.add_argument("-p", "--processor-number",
                    help="Number of available processors.",
                    type=int,
                    Default=7)
parser.add_argument("-z", "--parallel-designs",
                    help="Number of designs  carried out in parallel.",
                    type=int,
                    Default=None)
parser.add_argument("-z", "--parallel-designs",
                    help="Number of designs  carried out in parallel.",
                    type=int,
                    Default=None)
parser.add_argument("-f", "--flank-size",
                    help="Number of bases to flank target sites.",
                    type=int,
                    Default=250)
parser.add_argument("-g", "--single-mip-treshold",
                    help=("Target regions smaller than this will have a "
                          "single MIP designed."),
                    type=int,
                    Default=0)
parser.add_argument("-y", "--min-target-size",
                    help="Size threshold to eliminate un-aligned targets.",
                    type=int,
                    Default=None)
parser.add_argument("-", "--coordinates-file",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
parser.add_argument("-y", "--genes-file",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
parser.add_argument("-y", "--snps-file",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
parser.add_argument("-y", "--fasta-files",
                    help="Size threshold to eliminate un-aligned targets.",
                    nargs="*")
parser.add_argument("-y", "--fasta-capture-type",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default="whole")
parser.add_argument("-y", "--genomic-coverage-threshold",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
parser.add_argument("-y", "--genomic-identity-threshold",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
parser.add_argument("-y", "--local-coverage-threshold",
                    help="Size threshold to eliminate un-aligned targets.",
                    Default=None)
# parse arguments from command line
args = vars(parser.parse_args())

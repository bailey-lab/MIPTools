import argparse
from Bio import SeqIO
from itertools import chain
import os
import numpy as np
from random import sample
import re
import subprocess

# Parse input arguments
parser = argparse.ArgumentParser(
    description="""Downsample the number of UMIs sequenced per MIP."""
)
parser.add_argument(
    "-t",
    "--downsample-threshold",
    help="The threshold at which UMIs will be downsampled.",
    default=2000,
    type=int,
)
parser.add_argument(
    "-w",
    "--weighted",
    action="store_true",
    help="Whether to apply a weight when randomly sampling UMIs.",
)
parser.add_argument(
    "file",
    nargs="+",
    help="The files on which to downsample the UMIs.",
)
args = vars(parser.parse_args())
threshold = int(args["threshold"])

# Remove empty first element from list
if args["file"][0] == "":
	args["file"] = args["file"][1:]

# Iterate over all the files input
for file in args["file"]:
    # Unzip the file
    subprocess.run(["gzip", "-df", file])
    unzipped = os.path.splitext(file)[0]

    # Read the file
    records = list(SeqIO.parse(unzipped, "fastq"))

    # Randomly select a certain number of records. Either weigh by the read
    # count, or just randomly select a sample.
    if len(records) > threshold:
        if args["weighted"]:
            # Find the read counts for each UMI
            read_cnts = []
            for r in records:
                read_cnts.append(re.findall("readCnt=(\\d+)", r.id))

            # Flatten the list, convert to int, and make list sum to one
            read_cnts = list(chain.from_iterable(read_cnts))
            read_cnts = [int(x) for x in read_cnts]
            weights = [x / sum(read_cnts) for x in read_cnts]

            # Subset the list using weights
            subset = [
                records[i]
                for i in np.random.choice(
                    len(records), threshold, False, weights
                )
            ]
        else:
            subset = sample(records, threshold)

        # Write the subsetted file
        SeqIO.write(subset, unzipped, "fastq")

    # Zip the file
    subprocess.run(["gzip", "-f", unzipped])

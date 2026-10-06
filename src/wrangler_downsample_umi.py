import argparse
from Bio import SeqIO
from itertools import chain
import os
import numpy as np
import random
import re
import subprocess
import zlib
from multiprocessing import Pool

# Parse input arguments
parser = argparse.ArgumentParser(
    description="""Downsample the number of UMIs sequenced per MIP."""
)
parser.add_argument(
    "-c",
    "--cpu-count",
    help="The number of available processors to use.",
    default=1,
    type=int,
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
    "-s",
    "--seed",
    help="""Random seed for reproducible downsampling. Each file is seeded
    independently using this value combined with the file path, so the UMIs
    kept for a file do not depend on how files are distributed among worker
    processes. If not provided, selection is random and not reproducible.""",
    default=None,
    type=int,
)
parser.add_argument(
    "file",
    nargs="+",
    help="The files on which to downsample the UMIs.",
)
args = vars(parser.parse_args())
cpu_count = args["cpu_count"]
downsample_threshold = int(args["downsample_threshold"])
weighted = args["weighted"]
seed = args["seed"]

# Remove empty first element from list
if args["file"][0] == "":
    args["file"] = args["file"][1:]


def downsammple_fastq(file, downsample_threshold, weighted, seed=None):
    """Downsamples a FASTQ file by removing UMIs.

    Args:
        file (str): The path of the FASTQ file.
        downsample_threshold (int): The threshold at which UMIs will be
            downsampled.
        weighted (bool): Whether to downsample, weighing by the read count of
            each UMI.
        seed (int): Base random seed. If given, a per-file seed is derived from
            it and the file path, so results are reproducible regardless of
            worker scheduling. If None, the global (unseeded) generators are
            used.
    """
    # Set up random number generators
    if seed is None:
        py_rng = random
        np_rng = np.random
    else:
        file_seed = (seed + zlib.crc32(os.path.normpath(file).encode())) % 2**32
        py_rng = random.Random(file_seed)
        np_rng = np.random.default_rng(file_seed)

    # Unzip the file
    subprocess.run(["gzip", "-df", file])
    unzipped = os.path.splitext(file)[0]

    # Read the file
    records = list(SeqIO.parse(unzipped, "fastq"))

    # Randomly select a certain number of records. Either weigh by the read
    # count, or just randomly select a sample.
    if len(records) > downsample_threshold:
        if weighted:
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
                for i in np_rng.choice(
                    len(records), downsample_threshold, False, weights
                )
            ]
        else:
            subset = py_rng.sample(records, downsample_threshold)

        # Write the subsetted file
        SeqIO.write(subset, unzipped, "fastq")

    # Zip the file
    subprocess.run(["gzip", "-f", unzipped])


# Create a pool of workers to parallelize process
p = Pool(cpu_count)

# Iterate over all the files input
for file in args["file"]:
    p.apply_async(
        downsammple_fastq, [file, downsample_threshold, weighted, seed]
    )

# Close pool and wait for all child processes to terminate
p.close()
p.join()

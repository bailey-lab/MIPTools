import mip_functions as mip
import pickle
import os
import pandas as pd
import argparse
# Read input arguments
parser = argparse.ArgumentParser(
    description=""" Create a QC report following a bcl demultiplex operateon.
    """)
parser.add_argument(
    "-p", "--platform",
    help=("Sequencing platform."),
    required=True,
    choices=["nextseq", "miseq"]
)
# parse arguments from command line
args = vars(parser.parse_args())
platform = args["platform"]
stats_dir = "/opt/analysis/Stats"
bc_dict = "/opt/resources/barcode_dict.json"
# load barcode dict to be passed to the header-primer conversion function
with open(bc_dict, "rb") as infile:
    bc_dict = pickle.load(infile)
# fastq summary fies created for each lane contains raw read numbers
# and reads passing filter for tile and sample
fsums = []
# demultiplexing summary files have the information with the most popular
# unindexed sample barcodes.
dfiles = []
# scan the stats dir and extract information from the relevant files
for f in os.scandir(stats_dir):
    if f.name.startswith("FastqSummary"):
        fd = pd.read_table(f.path)
        lane = "Lane" + f.name.split(".")[0][-1]
        fd["Lane"] = lane
        fsums.append(fd)
    elif f.name.startswith("DemuxSummary"):
        lane = "Lane" + f.name.split(".")[0][-1]
        dfiles.append([f.path, lane])
dsums = []
for d, l in dfiles:
    with open(d) as infile:
        start = False
        for line in infile:
            if start:
                counts = line.strip().split("\t")
                counts.append(l)
                dsums.append(counts)
            elif line.startswith("### Columns: Index_Sequence Hit_Count"):
                start = True
# create summary dataframe for read counts
fsums = pd.concat(fsums)[["NumberOfReadsRaw", "NumberOfReadsPF", "Lane"]]
fsums = fsums.groupby("Lane").sum().reset_index()
print(("Total number of raw reads and reads passing filter were "
      "{0[NumberOfReadsRaw]} and {0[NumberOfReadsPF]}, respectively.").format(
    fsums.sum()
))
fsums.to_csv("/opt/analysis/Stats/ReadSummary.csv")
# create a dataframe with unindexed read information
dsums = pd.DataFrame(dsums, columns=["Header", "Read Count", "Lane"])
dsums["Read Count"] = dsums["Read Count"].astype(int)
# get primer indexes for corresponding headers
dsums["Fw,Rev"] = dsums["Header"].apply(
    lambda a: mip.header_to_primer(bc_dict, a, platform)
)
dsums["Fw"] = dsums["Fw,Rev"].apply(lambda a: a[0])
dsums["Rev"] = dsums["Fw,Rev"].apply(lambda a: a[1])
dsums["Index Difference"] = dsums["Rev"] - dsums["Fw"]
# separate 999 values which do not correspond to our indexes
caught = dsums.loc[(dsums["Fw"] != 999)
                   & (dsums["Rev"] != 999)]
print(("There were {} undetermined reads. {} of these belong to "
      "possible primer pairs.").format(dsums["Read Count"].sum(),
                                       caught["Read Count"].sum()))
dsums.to_csv("/opt/analysis/Stats/UndeterminedIndexSummary.csv")
caught.to_csv("/opt/analysis/Stats/UndeterminedPrimerSummary.csv")

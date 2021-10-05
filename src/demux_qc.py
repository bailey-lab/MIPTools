"""Generate demultiplexing statistics after a sequencine run."""
import mip_functions as mip
import pickle
import os
import pandas as pd
import argparse


def main(platform, stats_dir):
    """Generate demultiplexing statistics after a sequencine run."""
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
    fsums = pd.concat(fsums)

    # create per sample read count summary if sample sheet is provided
    try:
        sample_sheet = os.path.join(stats_dir, "SampleSheet.csv")
        sample_sheet_list = []
        with open(sample_sheet) as infile:
            parse_line = False
            for line in infile:
                if line.startswith("Sample_ID"):
                    parse_line = True
                if parse_line:
                    sample_sheet_list.append(line.split(","))
        sample_sheet = pd.DataFrame(sample_sheet_list[1:],
                                    columns=sample_sheet_list[0])
        sample_sheet["Sample_ID"] = sample_sheet["Sample_ID"].astype(int)
        sample_sheet = sample_sheet.rename(
            columns = {"Sample_ID": "SampleNumber", "Sample_Name": "Sample ID"})
        sample_sums = fsums.groupby(["SampleNumber", "Lane"], as_index=False)[
            ["NumberOfReadsRaw", "NumberOfReadsPF"]].sum()
        sample_sums = sample_sums.merge(sample_sheet)
        sample_sums.to_csv(
            os.path.join(stats_dir, "PerSampleReadCounts.csv"),
            index = False)
    except IOError:
        pass
    # Print out read number summary.
    fsums = fsums[["NumberOfReadsRaw", "NumberOfReadsPF", "Lane"]].groupby(
        "Lane").sum().reset_index()
    print(("Total number of raw reads and reads passing filter were "
           "{0[NumberOfReadsRaw]:,} and {0[NumberOfReadsPF]:,}, "
           "respectively.").format(
        fsums.sum()
    ))
    fsums.to_csv(os.path.join(stats_dir, "ReadSummary.csv"))
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
    print(("There were {:,} undetermined reads. {:,} of these belong to "
          "possible primer pairs.").format(dsums["Read Count"].sum(),
                                           caught["Read Count"].sum()))
    dsums.to_csv(os.path.join(stats_dir, "UndeterminedIndexSummary.csv"))
    caught.to_csv(os.path.join(stats_dir, "UndeterminedPrimerSummary.csv"))


if __name__ == "__main__":
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
    parser.add_argument(
        "-d", "--stats-dir",
        help=("Path to directory where demultiplexing stats are saved."),
        default="/opt/analysis/Stats"
    )

    # parse arguments from command line
    args = vars(parser.parse_args())
    main(platform=args["platform"], stats_dir=args["stats_dir"])

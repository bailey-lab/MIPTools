import pandas as pd
import argparse
import os
import numpy as np


def main(capture_384, sample_plates, sheets_96, output_file,
         wdir="/opt/analysis",
         quadrants="/opt/resources/sample_prep/quadrants.csv",
         forward_plates="/opt/resources/sample_prep/forward_plates.csv",
         reverse_plates="/opt/resources/sample_prep/reverse_plates.csv"):
    quad = pd.read_csv(quadrants)
    forward_plates = pd.read_csv(forward_plates)
    reverse_plates = pd.read_csv(reverse_plates)
    if capture_384 is not None:
        capture_384 = [os.path.join(wdir, p) for p in capture_384]
        capture_384 = [pd.read_table(p) for p in capture_384]
        capture_384 = pd.concat(capture_384, ignore_index=True, axis=0)
        sample_plates = [os.path.join(wdir, p) for p in sample_plates]
        sample_plates = [pd.read_table(p) for p in sample_plates]
        sample_plates = pd.concat(sample_plates, ignore_index=True, axis=0)
        sample_plates["row"] = sample_plates["sample_well"].apply(
            lambda a: a[0])
        sample_plates["column"] = sample_plates["sample_well"].apply(
            lambda a: int(a[1:]))
        plating_cols = ["sample_name", "sample_plate", "row", "column"]
        sample_plates = sample_plates.loc[:, plating_cols]
        plates_without_samples = set(capture_384["sample_plate"]).difference(
            sample_plates["sample_plate"])
        if len(plates_without_samples) > 0:
            print(("{} does not have corresponding sample plates.").format(
                plates_without_samples))
        samples_without_plates = set(sample_plates["sample_plate"]).difference(
            capture_384["sample_plate"])
        if len(samples_without_plates) > 0:
            print(("{} does not have corresponding capture plates.").format(
                samples_without_plates))
        captures = capture_384.merge(sample_plates)
        captures = captures.merge(forward_plates)
        captures = captures.merge(reverse_plates)
        captures = captures.merge(quad)
        captures = captures.drop_duplicates()

    if sheets_96 is not None:
        sheets_96 = [os.path.join(wdir, p) for p in sheets_96]
        sheets_96 = [pd.read_table(p) for p in sheets_96]
        sheets_96 = pd.concat(sheets_96, ignore_index=True, axis=0)
        sheets_96 = sheets_96.drop_duplicates()

    capture_columns = ["sample_name", "sample_set", "probe_set",
                       "replicate", "fw", "rev", "owner", "Library Prep",
                       "sample_plate", "capture_plate",
                       "capture_plate_row", "capture_plate_column",
                       "quadrant", "FW_plate", "REV_plate", "384 Column"]

    if (capture_384 is not None) and (sheets_96 is not None):
        com = pd.concat([captures, sheets_96], axis=0, ignore_index=True,
                        sort=False)
        com = com.loc[:, capture_columns]
    elif capture_384 is not None:
        # "replicate" column is only available in 96 well format
        # this will need to be set to NaN value if no 96 well sample sheet
        # was used.
        captures["replicate"] = np.nan
        com = captures.loc[:, capture_columns]
    elif sheets_96 is not None:
        com = sheets_96
    else:
        print("At least one sample sheet must be provided.")
        return

    com = com.drop_duplicates()

    def assign_replicate(replicates):
        replicates = list(replicates)
        group_size = len(replicates)
        reps_available = set(range(1, group_size + 1))
        reps_used = set(replicates)
        reps_available = reps_available.difference(reps_used)
        reps_available = sorted(reps_available, reverse=True)
        reps_used = set()
        for i in range(group_size):
            rep = replicates[i]
            if np.isnan(rep) or (rep in reps_used):
                rep = int(reps_available.pop())
                replicates[i] = rep
                reps_used.add(rep)
            else:
                replicates[i] = int(rep)
                reps_used.add(rep)
        return pd.Series(replicates)
    try:
        com["Replicate"] = com.groupby(["sample_name", "sample_set"])[
            "replicate"].transform(assign_replicate).astype(int)
        com.drop("replicate", inplace=True, axis=1)
        com.rename(columns={"Replicate": "replicate"}, inplace=True)
    except ValueError:
        print("Error in assigning replicates. Please make sure "
              "the 'sample_name' and 'sample_set' fields have "
              "valid, non-empty values in all provided files.")

    if com.shape[0] != (com.groupby(["fw", "rev"]).first().shape[0]):
        size_file = os.path.join(wdir, "repeating_primers.csv")
        print(("There are repeating forward/reverse primer pairs.\n"
               "Sample sheet will be split based on the probe sets used.\n"
               "Inspect {} for repeating primer information.").format(
                   size_file))
        c_size = com.groupby(["fw", "rev"]).size().sort_values(
            ascending=False)
        c_size = c_size.loc[c_size > 1].reset_index()
        com.merge(c_size).to_csv(size_file, index=False)
        gb = com.groupby("probe_set")
        for group_key in gb.groups:
            g = gb.get_group(group_key)
            if g.shape[0] != (g.groupby(["fw", "rev"]).first().shape[0]):
                size_file = os.path.join(wdir, group_key + "_repeating.csv")
                print(("There are repeating forward/reverse primer pairs "
                       "within probe set {}. Inspect {} and correct the "
                       "sample sheet before proceeding with demultiplexing."
                       ).format(group_key, size_file))
                g_size = g.groupby(["fw", "rev"]).size().sort_values(
                    ascending=False)
                g_size = g_size.loc[g_size > 1].reset_index()
                g.merge(g_size).to_csv(size_file, index=False)
            g.to_csv(os.path.join(wdir, group_key + "_" + output_file),
                     index=False, sep="\t")
    else:
        com.to_csv(os.path.join(wdir, output_file), index=False, sep="\t")

    for sample_id in com["sample_name"]:
        if not sample_id.replace("-", "").isalnum():
            print(("Sample names can only contain "
                   "alphanumeric characters and '-'. "
                   "{} has invalid characters. "
                   "This sample will not be processed.").format(sample_id))


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser(
        description=""" Create a sample sheet given sample prep files using
        384 well format as well as regular sample sheets.
        """)
    parser.add_argument("-c", "--capture-plates",
                        help=("Capture prep plate file(s)."),
                        nargs="*")
    parser.add_argument("-s", "--sample-plates",
                        help=("Sample plate file(s)."),
                        nargs="*")
    parser.add_argument("-t", "--sample-sheets",
                        help=("Finished sample sheet file(s)."),
                        nargs="*")
    parser.add_argument("-o", "--output-file",
                        help=("Output file name."),
                        required=True)
    parser.add_argument("-w", "--working-directory",
                        help=("Path to working directory."),
                        default="/opt/analysis")
    parser.add_argument("-q", "--quadrants",
                        help=("A file defining quadrants of a 384 plate."),
                        default=("/opt/resources/sample_prep/"
                                 "quadrants.csv"))
    parser.add_argument("-f", "--forward-plates",
                        help=("Forward  primer plate file."),
                        default=("/opt/resources/sample_prep/"
                                 "forward_plates.csv"))
    parser.add_argument("-r", "--reverse-plates",
                        help=("Reverse primer plate file."),
                        default=("/opt/resources/sample_prep/"
                                 "reverse_plates.csv"))
    args = vars(parser.parse_args())

    main(args["capture_plates"], args["sample_plates"], args["sample_sheets"],
         args["output_file"], args["working_directory"], args["quadrants"],
         args["forward_plates"], args["reverse_plates"])

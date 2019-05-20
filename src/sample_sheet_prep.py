import pandas as pd
import argparse
import os
import numpy as np


def main(capture_384, sample_plates, sheets_96, output_file,
         wdir="/opt/analysis",
         quadrants="/opt/base_resources/sample_prep/quadrants.csv",
         forward_plates="/opt/base_resources/sample_prep/forward_plates.csv",
         reverse_plates="/opt/base_resources/sample_prep/reverse_plates.csv"):
    quad = pd.read_csv(quadrants)
    forward_plates = pd.read_csv(forward_plates)
    reverse_plates = pd.read_csv(reverse_plates)
    if len(capture_384) > 0:
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

    if len(sheets_96) > 0:
        sheets_96 = [os.path.join(wdir, p) for p in sheets_96]
        sheets_96 = [pd.read_table(p) for p in sheets_96]
        sheets_96 = pd.concat(sheets_96, ignore_index=True, axis=0)
        sheets_96 = sheets_96.drop_duplicates()

    capture_columns = ["sample_name", "sample_set", "probe_set",
                       "replicate", "fw", "rev", "owner", "Library Prep",
                       "sample_plate", "capture_plate",
                       "capture_plate_row", "capture_plate_column",
                       "quadrant", "FW_plate", "REV_plate", "384 Column"]

    if (len(capture_384) > 0) and (len(sheets_96) > 0):
        com = pd.concat([captures, sheets_96], axis=0, ignore_index=True)
        com = com.loc[:, capture_columns]
    elif len(capture_384) > 0:
        com = captures.loc[:, capture_columns]
    elif len(sheets_96) > 0:
        com = sheets_96
    else:
        print("At least one sample sheet needs to be provided.")
        return

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

    com["Replicate"] = com.groupby(["sample_name", "sample_set"])[
        "replicate"].transform(assign_replicate).astype(int)
    com.drop("replicate", inplace=True, axis=1)
    com.rename(columns={"Replicate": "replicate"}, inplace=True)

    if com.shape[0] != (com.groupby(["fw", "rev"]).first().shape[0]):
        print("There are repeating forward/reverse primer pairs.\n"
              "Sample sheet will be split based on the probe sets used.")
        gb = com.groupby("probe_set")
        for group_key in gb.groups:
            g = gb.get_group(group_key)
            if g.shape[0] != (g.groupby(["fw", "rev"]).first().shape[0]):
                print(("There are repeating forward/reverse primer pairs "
                       "within probe set {}. Correct sample sheet before "
                       "proceeding with demultiplexing.").format(group_key))
            g.to_csv(os.path.join(wdir, group_key, "_", output_file))
    else:
        com.to_csv(os.path.join(wdir, output_file))


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
    parser.add_argument("-q", "--quadrants",
                        help=("A file defining quadrants of a 384 plate."),
                        default=("/opt/base_resources/sample_prep/"
                                 "quadrants.csv"))
    parser.add_argument("-f", "--forward-plates",
                        help=("Forward  primer plate file."),
                        default=("/opt/base_resources/sample_prep/"
                                 "forward_plates.csv"))
    parser.add_argument("-r", "--reverse-plates",
                        help=("Reverse primer plate file."),
                        default=("/opt/base_resources/sample_prep/"
                                 "reverse_plates.csv"))
    args = vars(parser.parse_args())

    main(args["capture_plates"], args["sample_plates"], args["sheets_96"],
         args["output_file"], args["quadrants"], args["forward_plates"],
         args["reverse_plates"])

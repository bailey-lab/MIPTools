import pandas as pd
import argparse
import os
import numpy as np
import subprocess
import pickle


def sample_sheet_prep(
        capture_plates, sample_plates, legacy_sheets, output_file,
        working_directory="/opt/analysis",
        quadrants="/opt/resources/sample_prep/quadrants.csv",
        forward_plates="/opt/resources/sample_prep/forward_plates.csv",
        reverse_plates="/opt/resources/sample_prep/reverse_plates.csv",
        barcode_dictionary="/opt/resources/sample_prep/barcode_dict.pickle",
        platform="nextseq",
        template_dir="/opt/resources/templates/sample_sheet_templates/"):
    quad = pd.read_csv(quadrants)
    forward_plates = pd.read_csv(forward_plates)
    reverse_plates = pd.read_csv(reverse_plates)
    wdir = working_directory
    if capture_plates is not None:
        capture_paths = [os.path.join(wdir, p) for p in capture_plates]
        capture_plates = []
        for p in capture_paths:
            try:
                capture_plates.append(pd.read_table(p).rename(
                    columns={"Library Prep": "library_prep"}))
            except IOError:
                raise Exception(("Error: Capture plate file {} does not exist "
                       "in the run directory {}").format(p, wdir))
        if len(capture_plates) > 0:
            capture_plates = pd.concat(capture_plates, ignore_index=True,
                                       axis=0)
            sample_paths = [os.path.join(wdir, p) for p in sample_plates]
            sample_plates = []
            for p in sample_paths:
                try:
                    sample_plates.append(pd.read_table(p))
                except IOError:
                    raise Exception(("Error: Sample plate file {} does not "
                           "exist in the run directory {}").format(p, wdir))
            if len(sample_plates) > 0:
                sample_plates = pd.concat(sample_plates, ignore_index=True,
                                          axis=0)
                sample_plates["row"] = sample_plates["sample_well"].apply(
                    lambda a: a[0])
                sample_plates["column"] = sample_plates["sample_well"].apply(
                    lambda a: int(a[1:]))
                plating_cols = ["sample_name", "sample_plate", "row", "column"]
                sample_plates = sample_plates.loc[:, plating_cols]
                # if a sample plate is referenced in capture plates but
                # sample plate information is not provided, raise an error.
                plates_without_samples = set(
                    capture_plates["sample_plate"]).difference(
                        sample_plates["sample_plate"])
                if len(plates_without_samples) > 0:
                    raise Exception(("Error: Sample plate(s) {} do not "
                           "exist.").format(plates_without_samples))
                # if a sample plate is provided but not referenced in the
                # captures, print a warning.
                samples_without_plates = set(
                    sample_plates["sample_plate"]).difference(
                        capture_plates["sample_plate"])
                if len(samples_without_plates) > 0:
                    print(("Warning: {} does not have corresponding capture "
                           "plates.").format(samples_without_plates))
                captures = capture_plates.merge(sample_plates)
                captures = captures.merge(forward_plates)
                captures = captures.merge(reverse_plates)
                captures = captures.merge(quad)
                captures = captures.drop_duplicates()
                captures["replicate"] = np.nan
            else:
                raise Exception("Error: sample plates {} were not found."
                    ).format(sample_plates)
        else:
            raise Exception("Error: capture plates {} were not found."
                ).format(capture_plates)
    else:
        capture_plates = []

    if legacy_sheets is not None:
        legacy_paths = [os.path.join(wdir, p) for p in legacy_sheets]
        legacy_sheets = []
        for p in legacy_paths:
            try:
                legacy_sheets.append(pd.read_table(p).rename(
                    columns={"Library Prep": "library_prep"}))
            except IOError:
                raise Exception(("Error: Sample sheet file {} does not exist "
                       "in the run directory {}").format(p, wdir))
        legacy_sheets = pd.concat(legacy_sheets, ignore_index=True, axis=0)
        legacy_sheets = legacy_sheets.drop_duplicates()
    else:
        legacy_sheets = []

    if (len(capture_plates) > 0) and (len(legacy_sheets) > 0):
        sample_sheet = pd.concat([captures, legacy_sheets], axis=0,
            ignore_index=True, sort=False)
    elif len(capture_plates) > 0:
        # "replicate" column is only available in 96 well format
        # this will need to be set to NaN value if no 96 well sample sheet
        # was used.
        captures["replicate"] = np.nan
        sample_sheet = captures
    elif len(legacy_sheets) > 0:
        sample_sheet = legacy_sheets
    else:
        raise Exception("At least one sample sheet must be provided.")

    sample_sheet = sample_sheet.drop_duplicates()

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
        sample_sheet["Replicate"] = sample_sheet.groupby(
            ["sample_name", "sample_set"])["replicate"].transform(
                assign_replicate).astype(int)
        sample_sheet.drop("replicate", inplace=True, axis=1)
        sample_sheet.rename(columns={"Replicate": "replicate"}, inplace=True)
    except ValueError:
        raise Exception("Error in assigning replicates. Please make sure "
              "the 'sample_name' and 'sample_set' fields have "
              "valid, non-empty values in all provided files.")

    # load barcode dictionary file to add sample barcode sequences
    with open(barcode_dictionary, "rb") as bc_file:
        barcodes = pickle.load(bc_file)

    # create the sample id as the combination of sample name, sample set
    # and replicate number
    sample_sheet["sample_id"] = sample_sheet[
        ["sample_name", "sample_set", "replicate"]].apply(
            lambda a: "-".join(map(str, a)), axis=1)

    # check sample ids for formatting.
    sample_sheet["valid_sample_id"] = sample_sheet["sample_id"].map(
        lambda a: a.replace("-", "").isalnum())
    invalid_samples = sample_sheet.loc[~sample_sheet["valid_sample_id"]]
    if invalid_samples.shape[0] > 0:
        invalid_samples_file = os.path.join(wdir, "invalid_samples.csv")
        invalid_samples.to_csv(invalid_samples_file)
        raise Exception(("Sample IDs can only contain alphanumeric characters"
               " and '-'. There are samples with invalid characters. "
               "Please correct the sample ids saved in {}")).format(
               invalid_samples_file)

    # get i5 and i7 index (sample barcode) sequences for each sample.
    # orientation of i5 sequences are different for miseq and nextseq
    sample_sheet["i7"] = sample_sheet["rev"].map(
        lambda a: barcodes[a]["index_sequence"])
    if platform == "nextseq":
        sample_sheet["i5"] = sample_sheet["fw"].map(
            lambda a: barcodes[a]["index_sequence"])
    elif platform == "miseq":
        sample_sheet["i5"] = sample_sheet["fw"].map(
            lambda a: barcodes[a]["sequence"])

    # reset the index to make sure index starts from zero
    sample_sheet.reset_index(inplace=True)
    # Create barcode names. This naming scheme is following earlier notation
    # used but probably any unique name would work.
    sample_sheet["reverse_index"] = "S" + sample_sheet["rev"].astype(str)
    sample_sheet["forward_index"] = "N" + sample_sheet["fw"].astype(str)
    # there are 4 fields that can be left empty in the sample sheet but
    # we still need to have those columns as empty strings
    empty_cols = ["empty" + str(i) for i in range(4)]
    for c in empty_cols:
        sample_sheet[c] = ""

    # check whether all required columns have valid values
    required_columns = ["sample_name", "sample_set", "probe_set",
                        "replicate", "fw", "rev", "library_prep",
                        "sample_id", "i5", "i7", "reverse_index",
                        "forward_index"]
    missing_values = sample_sheet[required_columns].isnull().any()
    missing_values = missing_values.loc[missing_values].index.to_list()
    if len(missing_values) > 0:
        bad_sample_sheet = os.path.join(wdir, "bad_sample_sheet.csv")
        sample_sheet.to_csv(bad_sample_sheet)
        raise Exception(("Error: Required column(s) {} cannot have missing "
            "values. Please inspect the file {}").format(
                ", ".join(missing_values), bad_sample_sheet))

    # create a small function to save files
    def make_sample_sheet(s_sheet, pfix):
        """Save sample sheet to temporary file and concatenate to template."""
        # select the columns needed in the final sample sheet
        cols = ["sample_id", "empty0", "empty1", "reverse_index", "i7",
                "forward_index", "i5", "empty2", "empty3"]
        # save sample sheet to a temporary file
        sample_sheet_tail = os.path.join(wdir, "temp_samples.csv")
        s_sheet.loc[:, cols].to_csv(sample_sheet_tail, index=True,
            header=False)
        # the sample sheet we generated needs some lines from the template file
        # we'll cat that file before the sample sheet tail just saved.
        sample_sheet_head = os.path.join(
            template_dir, platform + "_sample_sheet_template.csv")

        sample_sheet_file = os.path.join(wdir, pfix + "SampleSheet.csv")
        with open(sample_sheet_file, "w") as final_sample_sheet:
            res = subprocess.run(["cat", sample_sheet_head, sample_sheet_tail],
                                 stdout=final_sample_sheet,
                                 stderr=subprocess.PIPE)
            if res.stderr == b"":
                subprocess.run(["rm", sample_sheet_tail])
            else:
                raise Exception((
                    "Error creating final sample sheet file: {}").format(
                        res.stderr))
        # save the entire dataframe as _samples.tsv file
        s_sheet.to_csv(os.path.join(wdir, pfix + output_file),
            index=False, sep="\t")
        return


    # check if there are non-unique primer pairs
    if sample_sheet.shape[0] != (
            sample_sheet.groupby(["fw", "rev"]).first().shape[0]):
        size_file = os.path.join(wdir, "repeating_primers.csv")
        # nonunique primer pairs may be allowed if they belong to
        # different probe sets becouse the data can be separated based
        # on the probe sequences.
        print(("There are repeating forward/reverse primer pairs.\n"
               "Sample sheet will be split based on the probe sets used.\n"
               "Inspect {} for repeating primer information.").format(
                   size_file))
        c_size = sample_sheet.groupby(["fw", "rev"]).size().sort_values(
            ascending=False)
        c_size = c_size.loc[c_size > 1].reset_index()
        sample_sheet.merge(c_size).to_csv(size_file, index=False)

        # create a separate sample sheet file for each probe set
        gb = sample_sheet.groupby("probe_set")
        for group_key in gb.groups:
            g = gb.get_group(group_key)
            # raise error if non-unique primers exist within probesets
            if g.shape[0] != (g.groupby(["fw", "rev"]).first().shape[0]):
                size_file = os.path.join(wdir, group_key + "_repeating.csv")
                raise Exception(("There are repeating forward/reverse primer "
                       "pairs within probe set {}. Inspect {} and correct the "
                       "sample sheet.").format(group_key, size_file))
            make_sample_sheet(g, group_key + "_")
    else:
        # save a single sample sheet when all primer pairs are unique
        make_sample_sheet(sample_sheet, "")


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
    parser.add_argument("-t", "--legacy-sheets",
                        help=("Legacy sample sheet file(s)."),
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
    parser.add_argument("-b", "--barcode-dictionary",
                        help="Path to sample barcode dictionary.",
                        default=("/opt/resources/sample_prep/"
                                 "barcode_dict.pickle"))
    parser.add_argument("-p", "--platform",
                        help="Sequencing platform",
                        default="nextseq",
                        choices=["nextseq", "miseq"])
    parser.add_argument("-d", "--template-dir",
                        help="Directory containing sample sheet headers.",
                        default=("/opt/resources/templates/"
                                 "sample_sheet_templates/"))

    args = vars(parser.parse_args())

    sample_sheet_prep(**args)

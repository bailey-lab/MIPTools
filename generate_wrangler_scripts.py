import json
import pickle
import os
from itertools import zip_longest
import pandas as pd
import numpy as np
import argparse

# Read input arguments
parser = argparse.ArgumentParser(
    description=""" Generate bash scripts to be used for processing
    after a MIP sequencing run.
    """)
parser.add_argument(
    "-e", "--experiment-id",
    help=("A Unique id given to each sequencing run by the user."),
    required=True
)
parser.add_argument("-p", "--platform",
                    help="Sequencing platform.",
                    choices=["nextseq", "miseq"],
                    required=True)
parser.add_argument("-n", "--nextseq-id",
                    help="A unique id given by nextseq machine.",
                    default="")
parser.add_argument("-c", "--cpu-count",
                    type=int,
                    help="Number of available processors to use.",
                    default=1)
parser.add_argument("-s", "--server-num",
                    type=int,
                    help="Starting number for MIP server.",
                    default=1)
parser.add_argument("-x", "--access-token-file",
                    help="Basespace access token for user.",
                    default="/opt/resources/access_token.txt")
parser.add_argument("-d", "--raw-data-dir",
                    help=("Absolute path to base directory where sequencing "
                          "files should be saved to."),
                    default="/opt/work")
parser.add_argument("--processed-data-dir",
                    help=("Absolute path to base directory where "
                          "MIPWrangler output file should be copied to."),
                    default="/opt/work")
parser.add_argument("-a", "--analysis-data-dir",
                    help=("Absolute path to base directory for "
                          "MIPWrangler working directory."),
                    default="/opt/work")
parser.add_argument("-w", "--cluster-script",
                    help=("MIPWrangler script name. Absolute path"
                          "if not in $PATH."),
                    default="runMIPWranglerNoCutoffCurrent.sh")
parser.add_argument("-r", "--project-resource-dir",
                    help=("Path to directory where project specific resources "
                          "such as probe sets used, mip arm info etc. are"),
                    default="/opt/project_resources")
parser.add_argument("-r", "--base-resource-dir",
                    help=("Path to directory where general resources such as "
                          "barcode dictionary, sample sheet "
                          "templates etc. are."),
                    default="/opt/base_resources")
parser.add_argument("-l", "--sample-list",
                    help=("File providing a list of samples with associated "
                          "information."),
                    required=True)
# parse arguments from command line
args = vars(parser.parse_args())
experiment_id = args["experiment_id"]
platform = args["platform"]
nextseq_id = args["nextseq_id"]
cluster_script = args["cluster_script"]
cpu_count = args["cpu_count"]
server_num = args["server_num"]
access_token_file = args["access_token_file"]
raw_data_dir = os.path.abspath(args["raw_data_dir"])
analysis_data_dir = os.path.abspath(args["analysis_data_dir"])
processed_data_dir = os.path.abspath(args["processed_data_dir"])
project_resource_dir = os.path.abspath(args["project_resource_dir"])
base_resource_dir = os.path.abspath(args["base_resource_dir"])
sample_list_file = args["sample_list"]
new_mip_set = args["new_mip_set"]
new_mip_arms = args["new_mip_arms"]
experiment_name = experiment_id + "_" + platform
raw_dir = os.path.join(raw_data_dir, experiment_name)
sample_sheet_template = os.path.join(
    base_resource_dir,
    "templates",
    "sample_sheet_templates",
    platform + "_sample_sheet_template.csv"
)
raw_mip_ids_dir = os.path.join(raw_dir, "mip_ids")
sample_sheet = os.path.join(raw_mip_ids_dir, "SampleSheet.csv")
fastq_dir = os.path.join(raw_dir, "fastq")
analysis_dir = os.path.join(analysis_data_dir, experiment_name)
barcode_dict_file = os.path.join(
    base_resource_dir, "barcode_dict.json")
# create dirs if they do not exist
for d in [raw_mip_ids_dir, fastq_dir]:
    if not os.path.exists(d):
        os.makedirs(d)
# Nextseq data needs to be downloaded from BaseSpace and converted to
# fastq prior to demultiplexing.
if platform == "nextseq":
    download_commands = [
        ["cd", raw_dir],
        ["python /usr/bin/BaseSpaceRunDownloader_v2.py -r",
            nextseq_id, "-a", '"$(< ' + access_token_file + ')"']
        ]
    demux_commands = [
        ["cd", os.path.join(raw_dir, nextseq_id)],
        ["ulimit -n 9999"],
        [("nohup bcl2fastq -o ../fastq --sample-sheet "
          "../mip_ids/SampleSheet.csv --no-lane-splitting")]
    ]
# MiSeq run data is received as fastq or bcl (preferred), no download necessary
# When bcl is received, fastq conversion and demultiplexing is required.
else:
    download_commands = []
    demux_commands = [
        ["cd", os.path.join(raw_dir, "run_folder")],
        ["ulimit -n 9999"],
        [("nohup bcl2fastq -o ../fastq --sample-sheet"
            " ../mip_ids/SampleSheet.csv --no-lane-splitting")],
    ]
# First part of the MIPWrangler process is to extract the sequences and
# stitch forward and reverse reads. This is done with runGzExtractStitch
stitch_base = "nohup MIPWrangler runGzExtractStitch"
# stitch_base = "nohup MIPWrangler mipSetupAndExtractByArm"
stitch_commands = {}
wrangler_commands = {}
# Load the barcode dictionary which contains sequences of sample barcodes
with open(barcode_dict_file, "rb") as in1:
    barcode_dic = pickle.load(in1)
# read in sample information
sample_names = []
sample_info = {}
with open(sample_list_file) as infile:
    linenum = 0
    for line in infile:
        newline = line.strip().split("\t")
        if linenum == 0:
            colnames = newline
            linenum += 1
        else:
            sample_dict = {colname: colvalue for colname, colvalue
                           in zip(colnames, newline)}
            sample_set = sample_dict["sample_set"]
            sample_name = sample_dict["sample_name"]

            replicate_number = sample_dict["replicate"]
            forward_index = sample_dict["fw"]
            reverse_index = sample_dict["rev"]
            sample_id = "-".join([sample_name, sample_set, replicate_number])
            if sample_id in sample_info:
                print("Repeating sample name ", sample_id)
            if not sample_id.replace("-", "").isalnum():
                print(("Sample IDs can only contain "
                       "alphanumeric characters and '-'. "
                       "{} has invalid characters.").format(sample_id))
                continue
            if platform == "nextseq":
                sample_dict.update(
                    {"i7": barcode_dic[reverse_index]["index_sequence"],
                     "i5": barcode_dic[forward_index]["index_sequence"]})
            elif platform == "miseq":
                sample_dict.update(
                    {"i7": barcode_dic[reverse_index]["index_sequence"],
                     "i5": barcode_dic[forward_index]["sequence"]})
            sample_dict["sample_index"] = linenum
            linenum += 1
            sample_info[sample_id] = sample_dict
            sample_names.append(sample_id)
# Check for samples sharing one or both barcodes. One barcode sharing is
# allowed but a warning can be printed if desired by setting the warning
#  to True. If both barcodes are shared among two samples, those samples
# will be ignored and a message will be broadcast.
warnings = False
samples_sharing = []
for s1 in sample_info:
    for s2 in sample_info:
        if s1 != s2:
            if ((sample_info[s1]["fw"] == sample_info[s2]["fw"])
                    and (sample_info[s1]["rev"] == sample_info[s2]["rev"])):
                samples_sharing.append([s1, s2])
            elif warnings and (
                (sample_info[s1]["fw"] == sample_info[s2]["fw"])
                    or (sample_info[s1]["rev"] == sample_info[s2]["rev"])
            ):
                print("Samples %s and %s share a barcode" % (s1, s2))
samples_sharing_set = []
if len(samples_sharing) > 0:
    for s in samples_sharing:
        samples_sharing_set.extend(s)
    samples_sharing_set = set(samples_sharing_set)
    print("There are %d samples sharing the same barcode pair"
          % len(samples_sharing_set))
    pd.DataFrame(samples_sharing).to_csv(
        os.path.join(analysis_data_dir, "samples_sharing_barcodes.tsv"),
        sep="\t"
    )
# create sample sheet
with open(sample_sheet_template) as infile, open(sample_sheet, "w") as outfile:
    outfile_list = infile.readlines()
    outfile_list = [o.strip() for o in outfile_list]
    for sample_id in sample_names:
        if sample_id in samples_sharing_set:
            continue
        reverse_index = sample_info[sample_id]["rev"]
        forward_index = sample_info[sample_id]["fw"]
        sample_index = str(sample_info[sample_id]["sample_index"])
        outlist = [sample_index, sample_id, "", "",
                   "S" + reverse_index,
                   sample_info[sample_id]["i7"],
                   "N" + forward_index,
                   sample_info[sample_id]["i5"], "", ""]
        outfile_list.append(",".join(outlist))
    outfile.write("\n".join(outfile_list))
with open(os.path.join(raw_mip_ids_dir, "samples.dic"), "w") as outfile:
    json.dump(sample_info, outfile, indent=1)
# Sequencing runs are organized into sample sets which should be analyzed
# separately. The sample sets are listed in the sample list file.
# Probe sets column in the sample list indicates the probes used for a given
# sample. We'll separate samples into probe and sample sets  for analysis.
sample_sets = {}
for s in sample_info:
    sam = sample_info[s]
    s_set = sam["sample_set"]
    try:
        sample_sets[s_set]["sample_names"].append(s)
    except KeyError:
        sample_sets[s_set] = {"sample_names": [s],
                              "probe_sets": [],
                              "probe_set_strings": []}
    for k in sam:
        if k.startswith("probe_set"):
            sample_sets[s_set]["probe_set_strings"].extend(sam[k].split(";"))
for s_set in sample_sets:
    sample_sets[s_set]["probe_set_strings"] = list(
        set(sample_sets[s_set]["probe_set_strings"]))
    for pss in sample_sets[s_set]["probe_set_strings"]:
        sample_sets[s_set]["probe_sets"].append(pss.split(","))
# If a new mip set is used, update the mipsets.csv file and run the following
mipset_table = os.path.join(project_resource_dir, "mip_ids", "mipsets.csv")
mipsets = pd.read_csv(mipset_table)
mipset_list = mipsets.to_dict(orient="list")
# convert the mip sets dataframe to dict for easy access
all_probes = {}
# keep mip arm files for each mip set in a dictionary
mip_arms_dict = {}
for mipset in mipset_list:
    list_m = mipset_list[mipset]
    # the file name should be the second line in the mipsets.csv
    mip_arms_dict[mipset] = list_m[0]
    # rest of the lines have probe names in the set
    set_m = set(list_m[1:])
    set_m.remove(np.nan)
    all_probes[mipset] = set_m
# save probe dictionary to file
probe_set_file = os.path.join(project_resource_dir,
                              "mip_ids", "probe_sets.json")
with open(probe_set_file, "w") as outfile:
    json.dump(all_probes, outfile)
subset_names = []
# For each sample and probe set create
# 1) MIPWrangler input files (samples etc.)
# 2) Scripts for MIPWrangler Part I (extract + stitch)
# 3) Scripts for MIPWrangler Part II (clustering)
for s_set in sample_sets:
    sample_subset = sample_sets[s_set]["sample_names"]
    probe_sets = sample_sets[s_set]["probe_sets"]
    for pset_names in probe_sets:
        probes = set()
        mip_arms_list = []
        for p_name in pset_names:
            probes.update(all_probes[p_name])
            arm_file = os.path.join(project_resource_dir,
                                    "mip_ids",
                                    mip_arms_dict[p_name])
            try:
                with open(arm_file) as infile:
                    mip_arms_list.append(pd.read_table(infile))
            except IOError:
                print(("MIP arm file {} is required but missing for "
                      "the probe set {}").format(arm_file, p_name))
        mip_arms_table = pd.concat(mip_arms_list,
                                   ignore_index=True).drop_duplicates()
        mip_arms_table = mip_arms_table.loc[
            mip_arms_table["mip_id"].isin(probes)
        ]
        mip_family_names = mip_arms_table["mip_family"].tolist()
        subset_name = s_set + "_" + "_".join(pset_names)
        # Create MIPWrangler Input files
        with open(
            os.path.join(
                raw_mip_ids_dir,
                subset_name + "_allMipsSamplesNames.tab.txt"
            ), "w"
        ) as outfile:
            outfile_list = ["\t".join(["mips", "samples"])]
            mips_samples = zip_longest(
                mip_family_names, sample_subset, fillvalue=""
            )
            for ms in mips_samples:
                outfile_list.append("\t".join(ms))
            outfile.write("\n".join(outfile_list))
            pd.DataFrame(mip_arms_table).groupby(
                "mip_id").first().reset_index().dropna(
                    how="all", axis=1
                    ).to_csv(
                        os.path.join(
                            raw_mip_ids_dir,
                            subset_name + "_mipArms.txt"
                        ), sep="\t",
                        index=False
                    )
        # Create MIPWrangler part I script commands
        stitch_commands[subset_name] = [
            ["mkdir -p", os.path.join(analysis_dir, subset_name)],
            ["cd", os.path.join(analysis_dir,  subset_name)],
            [
                stitch_base, "--mipArmsFilename",
                os.path.join(raw_mip_ids_dir, subset_name + "_mipArms.txt"),
                "--mipSampleFile",
                os.path.join(
                    raw_mip_ids_dir,
                    subset_name + "_allMipsSamplesNames.tab.txt"
                ),
                "--numThreads",
                str(cpu_count),
                "--masterDir analysis",
                "--dir",
                fastq_dir,
            ]
        ]
        # Create MIPWrangler part II script commands
        wrangler_commands[subset_name] = [
            ["cd",
             "analysis"],
            ["nohup",
             cluster_script,
             str(server_num),
             str(cpu_count)],
            ["rsync -a " + os.path.join(
                analysis_dir,  subset_name +
                "/analysis/populationClustering/allInfo.tab.txt"
            ), os.path.join(processed_data_dir,
                            experiment_id + "_" + subset_name + ".txt")
            ]
        ]
        server_num += 1
        if subset_name in subset_names:
            print("%s is already in subset_names!" % subset_name)
        subset_names.append(subset_name)
# Save all scripts to files.
all_commands = download_commands + demux_commands
with open(os.path.join(raw_dir, "download_and_demux.sh"), "w") as outfile:
    outfile.write("\n".join([" ".join(a) for a in all_commands]) + "\n")
for s in subset_names:
    with open(os.path.join(raw_dir, s + ".sh"), "w") as outfile:
        outfile.write("\n".join(
            [" ".join(c) for c in stitch_commands[s]]) + "\n")
        outfile.write("\n".join(
            [" ".join(c) for c in wrangler_commands[s]]) + "\n")
    all_commands.extend(stitch_commands[s])
    all_commands.extend(wrangler_commands[s])
with open(os.path.join(raw_dir, "wrangle_combine.sh"), "w") as outfile:
    outfile.write("\n".join([" ".join(a) for a in all_commands]) + "\n")

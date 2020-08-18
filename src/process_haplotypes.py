import mip_functions as mip
import pandas as pd
import argparse
import os
import subprocess

# Read input arguments
parser = argparse.ArgumentParser(description="""Process MIP haplotypes.""")
parser.add_argument("-d", "--data-dir",
                    help="Path to directory containing the data.",
                    default="/opt/data")
parser.add_argument("-a", "--analysis-dir",
                    help="Path to directory to use for analysis.",
                    default="/opt/analysis")
parser.add_argument("-f", "--files",
                    help="MIPWrangler output file to be processed.",
                    action="append")
parser.add_argument("-s", "--sample-sheets",
                    help="Sample sheet belonging to MIPWrangler output.",
                    action="append")
parser.add_argument("-g", "--sample-groups",
                    help="Sample group to analyze.",
                    action="append",
                    default=None)
parser.add_argument("-p", "--probe-sets",
                    help="Probe sets to process.",
                    action="append",
                    default=None)
parser.add_argument("-u", "--probe-sets-used",
                    help="Comma separated probe sets to analyze.",
                    default=None)
parser.add_argument("-t", "--threads",
                    help="Number of CPU threads to use.",
                    default="1")
parser.add_argument("-e", "--species",
                    help="Species/genome to use (pf, hg38, etc.)",
                    required=True)
parser.add_argument("-w", "--wrangler-process",
                    help="Set this flag to process MIPWrangler output.",
                    action="store_true")
parser.add_argument("-m", "--map-haplotypes",
                    help="Set this flag to map haplotypes to genome.",
                    action="store_true")
parser.add_argument("-x", "--min-haplotype-umi",
                    help=("Minimum UMI requirement for a haplotype "
                          "across samples."),
                    type=int,
                    default=1)
parser.add_argument("-y", "--min-haplotype-sample",
                    help=("Minimum observed sample requirement for a "
                          "haplotype."),
                    type=int,
                    default=1)
parser.add_argument("-z", "--min-haplotype-fraction",
                    help=("Minimum observed sample fraction for a "
                          "haplotype."),
                    type=float,
                    default=0.0001)

# parse arguments from command line
args = vars(parser.parse_args())
data_dir = args["data_dir"]
wdir = args["analysis_dir"]
info_files = args["files"]
sample_sheets = args["sample_sheets"]
groups = args["sample_groups"]
probes = args["probe_sets"]
species = args["species"]
probe_sets_used = args["probe_sets_used"]
wrangler_process = args["wrangler_process"]
processorNumber = args["threads"]

bwaExtra = ["-t", processorNumber]

# copy the template settings file
temp_settings_file = (
    "/opt/resources/templates/analysis_settings_templates/settings.txt")
subprocess.run(["scp", temp_settings_file,
                "/opt/analysis/template_settings.txt"])

# extract the settings template
temp_settings = mip.get_analysis_settings(
    "/opt/analysis/template_settings.txt")

# update bwa settings with the options set above
bwaOptions = temp_settings["bwaOptions"]
try:
    bwaOptions.extend(bwaExtra)
except AttributeError:
    bwaOptions = [bwaOptions]
    bwaOptions.extend(bwaExtra)


# create a dictionary for which settings should be updated
# using the user specified parameters.
update_keys = {"processorNumber": processorNumber,
               "bwaOptions": bwaOptions,
               "species": species}
# update the settings
for k, v in update_keys.items():
    temp_settings[k] = v
# create a settings file in the analysis directory.
settings_file = "settings.txt"
settings_path = os.path.join(wdir, settings_file)
mip.write_analysis_settings(temp_settings, settings_path)
settings = mip.get_analysis_settings(settings_path)

if wrangler_process:
    if (info_files is None) or (sample_sheets is None):
        print(("Please provide MIPWrangler output and sample sheets for "
               "processing or do not set -w/--wrangler-process flag to "
               "skip MIPWrangler output processing."))
    else:
        try:
            info_files = [os.path.join(data_dir, i) for i in info_files]
            sample_sheets = [os.path.join(data_dir, s) for s in sample_sheets]
        except IOError:
            print(("Please provide the path for MIPWrangler output and sample "
                   "sheets relative to the data directory for processing or "
                   "do not set -w/--wrangler-process flag to "
                   "skip MIPWrangler output processing."))
        else:
            if len(info_files) != (len(sample_sheets)):
                print("Please provide one sample sheet per MIPWrangler "
                      "output.")
            else:
                if (groups is not None) and (probes is not None):
                    sample_groups = list(zip(groups, probes))
                elif (groups is None) and (probes is not None):
                    print(("Probes used supplied but sample sets were not. "
                           "Analysis will be carried out on all data."))
                elif (groups is None) and (probes is not None):
                    print(("sample sets supplied but probes used were not. "
                           "Analysis will be carried out on all data."))
                else:
                    sample_groups = None

                if probe_sets_used is None:
                    if probes is None:
                        probes = []
                        ss = pd.concat(pd.read_table(s for s in sample_sheets),
                                       ignore_index=True)
                        for p in ss["probe_set"].unique():
                            probes.extend(p.split(";"))
                    probe_sets_used = []
                    for p in probes:
                        probe_sets_used.extend(p.split(","))
                else:
                    probe_sets_used = probe_sets_used.split(",")
                probe_sets_used.append("")
                # Pass the probe_sets to be analyzed to  settings
                settings["mipSetKey"] = list(set(probe_sets_used))
                mip.write_analysis_settings(settings, settings_path)
                settings = mip.get_analysis_settings(settings_path)

                if len(info_files) > 1:
                    mip.combine_info_files(wdir,
                                           settings_file,
                                           info_files,
                                           sample_sheets,
                                           settings["mipsterFile"],
                                           sample_sets=sample_groups)
                else:
                    mip.process_info_file(wdir,
                                          settings_file,
                                          info_files,
                                          sample_sheets,
                                          settings["mipsterFile"],
                                          sample_sets=sample_groups)

if args["map_haplotypes"]:
    # filter haplotype sequences based on the number of total supporting UMIs
    settings["minHaplotypeBarcodes"] = args["min_haplotype_umi"]
    # filter haplotype sequences based on the number of samples they were
    # observed in
    settings["minHaplotypeSamples"] = args["min_haplotype_sample"]
    # filter haplotype sequences based on the fraction of samples they were
    # observed in
    settings["minHaplotypeSampleFraction"] = args["min_haplotype_fraction"]
    mip.get_vcf_haplotypes(settings)
    mip.get_haplotype_counts(settings)

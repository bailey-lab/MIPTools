import argparse
import pickle
import mip_functions as mip

# Read input arguments
parser = argparse.ArgumentParser(
    description=""" Design MIP probes.""")
parser.add_argument("-t", "--design-info",
                    help=("Path to design info file created by the region "
                          "prep module."),
                    Default="/opt/project_resources/design_info.json")
parser.add_argument("-e", "--exclude",
                    help=("Target names to exclude"),
                    Default=None,
                    nargs="*")
parser.add_argument("-i", "--include",
                    help=("Target names to include."),
                    Default=None,
                    nargs="*")
# parse arguments from command line
args = vars(parser.parse_args())
design_info_file = args["design_info"]
include = args["include"]
exclude = args["exclude"]
with open(design_info_file, "rb") as infile:
    design_info = pickle.load(infile)
design_dir = design_info[design_info.keys()[0]]["design_dir"]
parallel_processes = design_info[design_info.keys()[0]][
    "parallel_processes"]
if include is not None:
    targets = set(design_info.keys()).intersection(include)
    if len(targets) == 0:
        print(("None of {} is in targets.").format(include))
elif exclude is not None:
    targets = set(design_info.keys()).difference(exclude)
    if len(targets) == 0:
        print(("No target was left after removing excluded targets: {}."
               ).format(exclude))
else:
    targets = set(design_info.keys())

if len(targets) == 0:
    print("No targets remain for MIP design, exiting.")
else:
    mip.design_mips_multi(design_dir, targets, parallel_processes)
    resource_dir = design_info[design_info.keys()[0]]["resource_dir"]
    mip.parasight(resource_dir, design_info_file,
                  designed_gene_list=None, extra_extension=".extra",
                  use_json=True)

import argparse
import pickle
import json
import os
import probe_summary_generator
# Read input arguments
parser = argparse.ArgumentParser(
    description=""" Order MIP probes.""")
parser.add_argument("-d", "--design-info",
                    help=("Path to design info file created by the region "
                          "prep module."),
                    Default="/opt/project_resources/design_info.json")
parser.add_argument("-r", "--resource-dir",
                    help=("Path to directory where output files are saved."),
                    Default="/opt/project_resources/")
parser.add_argument("-f", "--filter-type",
                    help=(("If a selected MIP file is provided, how should "
                           "the listed MIPs be handled.")),
                    Default=None,
                    choices=["remove", "keep", None])
parser.add_argument("-i", "--include",
                    help=("Target names to include."),
                    Default=None,
                    nargs="*")

args = vars(parser.parse_args())
design_info_file = args["design_info"]
with open(design_info_file) as infile:
    design_info = json.load(infile)
filter_type = args["filter_type"]
mip_info = {}
call_info = {}
finished_genes = []
unfinished_genes = []


for gene_name in design_info:
    try:
        des = design_info[gene_name]["design_dir"]
        par_file = os.path.join(des, gene_name, gene_name)
        with open(par_file, "rb") as infile:
            Par = pickle.load(infile)
        res = Par.order_mips(filter_type=filter_type)
        mip_info[gene_name] = res["mip_info"]
        call_info[gene_name] = res["call_info"]
        finished_genes.append(gene_name)
    except Exception as e:
        unfinished_genes.append([gene_name, e])

resource_dir = args["resource_dir"]
mip_ids_dir = os.path.join(resource_dir, "mip_ids")
if not os.path.exists(mip_ids_dir):
    os.makedirs(mip_ids_dir)

mip_info_file = os.path.join(mip_ids_dir, "mip_info.json")
with open(mip_info_file, "w") as outfile:
    json.dump(mip_info, outfile, indent=1)

call_info_file = os.path.join(mip_ids_dir, "call_info.json")
with open(call_info_file, "w") as outfile:
    json.dump(call_info, outfile, indent=1)

unf_file = os.path.join(resource_dir, "failed_targets.json")
with open(unf_file, "w") as outfile:
    json.dump(unfinished_genes, outfile, indent=1)

fin_file = os.path.join(resource_dir, "completed_targets.json")
with open(fin_file, "w") as outfile:
    json.dump(finished_genes, outfile, indent=1)

probe_file = os.path.join(mip_ids_dir, "probe_info.csv")
probe_summary_generator.get_probe_info(mip_info_file=mip_info_file,
                                       output_file=probe_file)

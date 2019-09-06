import argparse
import json
import mip_functions as mip

# Read input arguments
parser = argparse.ArgumentParser(
    description=""" Design MIP probes.""")
parser.add_argument("-d", "--design-info",
                    help=("Path to design info file created by the region "
                          "prep module."),
                    default="/opt/project_resources/design_info.json")

# parse arguments from command line
args = vars(parser.parse_args())
design_info_file = args["design_info"]
with open(design_info_file) as infile:
    design_info = json.load(infile)
random_key = list(design_info.keys())[0]
design_dir = design_info[random_key]["design_dir"]
resource_dir = design_info[random_key]["resource_dir"]
mip.parasight(resource_dir, design_info_file,
              designed_gene_list=None, extra_extension=".extra",
              use_json=True)
mip.parasight_print(resource_dir, design_dir, design_info_file,
                    designed_gene_list=None, extra_extension=".extra",
                    use_json=True, print_out=False)

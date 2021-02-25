import argparse
import json
import mip_functions as mip


def design(design_info_file="/opt/project_resources/design_info.json",
           exclude=None, include=None):
    """Design MIPs following region prep."""
    with open(design_info_file) as infile:
        design_info = json.load(infile)
    random_key = list(design_info.keys())[0]
    design_dir = design_info[random_key]["design_dir"]
    parallel_processes = design_info[random_key]["parallel_processes"]
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
        resource_dir = design_info[random_key]["resource_dir"]
        mip.parasight(resource_dir, design_info_file,
                      designed_gene_list=None, extra_extension=".extra",
                      use_json=True)
        mip.parasight_print(resource_dir, design_dir, design_info_file,
                            designed_gene_list=None, extra_extension=".extra",
                            use_json=True, print_out=False)


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser(
        description=""" Design MIP probes.""")
    parser.add_argument("-d", "--design-info-file",
                        help=("Path to design info file created by the region "
                              "prep module."),
                        default="/opt/project_resources/design_info.json")
    parser.add_argument("-e", "--exclude",
                        help=("Target names to exclude"),
                        default=None,
                        nargs="*")
    parser.add_argument("-i", "--include",
                        help=("Target names to include."),
                        default=None,
                        nargs="*")
    # parse arguments from command line
    args = vars(parser.parse_args())

    design(**args)

"""Prepare finalized files for MIP order and MIP use."""
import argparse
import pickle
import json
import os
import probe_summary_generator
import pandas as pd
from itertools import product
import mip_functions as mip


def order_mips(design_name,
               design_info_file="/opt/project_resources/design_info.json",
               resource_dir="/opt/project_resources/",
               filter_type="remove"):
    """Prepare finalized files for MIP order and MIP use."""
    with open(design_info_file) as infile:
        design_info = json.load(infile)

    mip_info = {}
    call_info = {}
    finished_genes = []
    no_mip_genes = []
    unfinished_genes = []

    for gene_name in design_info:
        try:
            des = design_info[gene_name]["design_dir"]
            par_file = os.path.join(des, gene_name, gene_name)
            with open(par_file, "rb") as infile:
                Par = pickle.load(infile)
            res = Par.order_mips(filter_type=filter_type)
            if len(res["call_info"]) > 0:
                call_info[gene_name] = res["call_info"]
                mip_info[gene_name] = res["mip_info"]
                finished_genes.append(gene_name)
            else:
                no_mip_genes.append(gene_name)
        except Exception as e:
            unfinished_genes.append([gene_name, str(e)])

    mip_ids_dir = os.path.join(resource_dir, "mip_ids")
    if not os.path.exists(mip_ids_dir):
        os.makedirs(mip_ids_dir)

    mip_info_file = os.path.join(mip_ids_dir, "mip_info.json")
    with open(mip_info_file, "w") as outfile:
        json.dump(mip_info, outfile, indent=1)

    call_info_file = os.path.join(mip_ids_dir, "call_info.json")
    with open(call_info_file, "w") as outfile:
        json.dump(call_info, outfile, indent=1)

    unf_file = os.path.join(resource_dir, "failed_pipeline.json")
    with open(unf_file, "w") as outfile:
        json.dump(unfinished_genes, outfile, indent=1)

    fin_file = os.path.join(resource_dir, "completed_targets.json")
    with open(fin_file, "w") as outfile:
        json.dump(finished_genes, outfile, indent=1)

    no_mip_file = os.path.join(resource_dir, "failed_designs.json")
    with open(no_mip_file, "w") as outfile:
        json.dump(no_mip_genes, outfile, indent=1)

    print(("{} target designs were succesfully finished, with at least one "
           "MIP. See {} for a list of completed targets.").format(
               len(finished_genes), fin_file))

    print(("{} target designs went through the pipeline without errors, but "
           "did not yield any MIPs See {} for a list of these "
           "targets.").format(
               len(no_mip_genes), no_mip_file))

    print(("{} target designs failed to go through the pipeline."
           " See {} for a list of failed targets.").format(
               len(unfinished_genes), unf_file))

    probe_file = os.path.join(mip_ids_dir, "probe_info.csv")
    probe_summary_generator.get_probe_info(mip_info_file=mip_info_file,
                                           output_file=probe_file)

    probe_df = pd.read_csv(probe_file)[
        ["probe_id", "Probe Sequence", "MIP", "Gene"]]
    probe_df = probe_df.groupby(
        ["Probe Sequence", "MIP"]).first().reset_index()
    probe_count = len(probe_df)
    plate_count = probe_count // 96
    plate_count += 1
    plates = [design_name + "_" + str(c) for c in range(plate_count)]
    n_rep = "(N)"
    first_n = "(N:25252525)"
    probe_df["Sequence"] = probe_df["Probe Sequence"].apply(
        lambda a: a.replace("N", n_rep).replace(n_rep, first_n, 1))
    probe_df["probe_number"] = probe_df["MIP"].apply(
        lambda a: int(a.split("_")[-1][3:]))

    probe_df = probe_df.sort_values(["Gene", "probe_number"])
    wells = [well[0] + str(well[1]) for well in
             product(["A", "B", "C", "D", "E", "F", "G", "H"], range(1, 13))]
    wells = product(plates, wells)
    ind = pd.MultiIndex.from_tuples(
        list(wells)[: len(probe_df)], names=["Plate", "WellPosition"])
    probe_df.index = ind
    probe_df = probe_df.reset_index()
    probe_df["Name"] = probe_df["MIP"] + "_" + probe_df["probe_id"]
    gb = probe_df[["Plate", "WellPosition", "Name", "Sequence"]].groupby(
        "Plate")

    with pd.ExcelWriter("/opt/project_resources/probe_order.xlsx") as outfile:
        for gname in gb.groups:
            gr = gb.get_group(gname)
            gr.to_excel(outfile, sheet_name=gname, index=False)

    mip.parasight(resource_dir, design_info_file,
                  designed_gene_list=None, extra_extension=".extra",
                  use_json=True)
    random_key = list(design_info.keys())[0]
    design_dir = design_info[random_key]["design_dir"]
    mip.parasight_print(resource_dir, design_dir, design_info_file,
                        designed_gene_list=None, extra_extension=".extra",
                        use_json=True, print_out=False)


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser(
        description="""Order MIP probes.""")

    parser.add_argument("-d", "--design-info-file",
                        help=("Path to design info file created by the region "
                              "prep module."),
                        default="/opt/project_resources/design_info.json")

    parser.add_argument("-r", "--resource-dir",
                        help="Path to directory where output files are saved.",
                        default="/opt/project_resources/")

    parser.add_argument("-f", "--filter-type",
                        help=(("If a selected MIP file is provided, how "
                               "should the listed MIPs be handled.")),
                        default="remove",
                        choices=["remove", "keep"])

    parser.add_argument("-n", "--design-name",
                        help=("A name given this MIP design project."),
                        required=True)

    args = vars(parser.parse_args())

    order_mips(**args)

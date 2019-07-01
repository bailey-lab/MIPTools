import pandas as pd
from itertools import product
import json
import os


def generate_probe_order(probe_file, design_name, filter_probes,
                         mip_info_file, call_info_file):
    probe_df = pd.read_csv(probe_file)
    if filter_probes:
        try:
            probe_df = probe_df.loc[probe_df["Filter"] == "PASS"]
            dirname = os.path.dirname(probe_file)
            basename = os.path.basename(probe_file)
            filtered_name = os.path.join(dirname, "filtered_" + basename)
            probe_df.to_csv(filtered_name)
        except KeyError:
            print("Probe file does not have 'Filter' field."
                  " All probes will be used for the order.")
            filter_probes = False
    probe_df = probe_df[["probe_id", "Probe Sequence", "MIP", "Gene"]]
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
            gr[["WellPosition", "Name", "Sequence"]].to_excel(
                outfile, sheet_name=gname, index=False)

    # uptade mip_info and call info files to remove filtered MIPs
    if filter_probes:
        filtered_mips = set(probe_df["MIP"])
        try:
            with open(call_info_file) as infile:
                call_info = json.load(infile)
            for g in list(call_info.keys()):
                for m in list(call_info[g].keys()):
                    if m not in filtered_mips:
                        call_info[g].pop(m)
                if len(call_info[g]) == 0:
                    call_info.pop(g)
            dirname = os.path.dirname(call_info_file)
            basename = os.path.basename(call_info_file)
            filtered_name = os.path.join(dirname, "filtered_" + basename)
            with open(filtered_name, "w") as outfile:
                json.dump(call_info, outfile, indent=1)
        except IOError:
            print(("call_info file {} does not exist, so it will not be "
                   "updated for the filtered MIPs.").format(call_info_file))

        try:
            with open(mip_info_file) as infile:
                mip_info = json.load(infile)
            for g in list(mip_info.keys()):
                for m in list(mip_info[g]["mips"].keys()):
                    if m not in filtered_mips:
                        mip_info[g]["mips"].pop(m)
                if len(mip_info[g]["mips"]) == 0:
                    mip_info.pop(g)
            dirname = os.path.dirname(mip_info_file)
            basename = os.path.basename(mip_info_file)
            filtered_name = os.path.join(dirname, "filtered_" + basename)
            with open(filtered_name, "w") as outfile:
                json.dump(mip_info, outfile, indent=1)
        except IOError:
            print(("mip_info file {} does not exist, so it will not be updated"
                   " for the filtered MIPs.").format(mip_info_file))


if __name__ == "__main__":
    import argparse
    # Read input arguments
    parser = argparse.ArgumentParser(
        description=("""Create a .xlsx file for probe order from a
                     probe summary file."""))

    parser.add_argument("-m", "--mip-info-file",
                        help=("Path to mip info file."),
                        default=None)

    parser.add_argument("-c", "--call-info-file",
                        help=("Path to call info file."),
                        default=None)

    parser.add_argument("-f", "--filter-probes",
                        help=("Should the probes be filtered based on ."
                              "'Filter' field in probe info file."),
                        action="store_true")

    parser.add_argument("-p", "--probe-info-file",
                        help=("Path to probe info file."),
                        required=True)

    parser.add_argument("-n", "--name",
                        help=("A name for the probe order."),
                        required=True)

    args = vars(parser.parse_args())

    filter_probes = args["filter_probes"]
    mip_info_file = args["mip_info_file"]
    call_info_file = args["call_info_file"]
    design_name = args["name"]
    probe_file = args["probe_info_file"]
    if filter_probes and ((mip_info_file is None)
                          or (call_info_file is None)):
        print("mip_info_file and call_info_file must be specified if "
              "'-f' flag is used to filter probes.")
    else:
        generate_probe_order(probe_file, design_name, filter_probes,
                             mip_info_file, call_info_file)

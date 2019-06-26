import pandas as pd
from itertools import product


def generate_probe_order(probe_file, design_name):
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
            gr[["WellPosition", "Name", "Sequence"]].to_excel(
                outfile, sheet_name=gname, index=False)


if __name__ == "__main__":
    import argparse
    # Read input arguments
    parser = argparse.ArgumentParser(
        description=("""Create a .xlsx file for probe order from a
                     probe summary file."""))

    parser.add_argument("-f", "--probe-info-file",
                        help=("Path to probe info file."),
                        default=("/opt/project_resources/mip_ids/"
                                 "probe_info.csv"))

    parser.add_argument("-n", "--name",
                        help=("A name for the probe order."),
                        required=True)

    args = vars(parser.parse_args())
    generate_probe_order(args["probe_info_file"], args["name"])

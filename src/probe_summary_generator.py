import json
import pandas as pd
import mip_functions as mip
import os

def get_probe_info(probe_set_keys=None, output_file=None,
                   probe_sets_file=None, mip_info_file=None):
    used_probes = set()
    if probe_sets_file is None:
        probe_sets_file = "/opt/project_resources/mip_ids/probe_sets.json"
    if mip_info_file is None:
        mip_info_file = "/opt/project_resources/mip_ids/mip_info.json"
    if probe_set_keys is not None:
        for psk in probe_set_keys:
            with open(probe_sets_file) as infile:
                used_probes.update(json.load(infile)[psk])

    with open(mip_info_file) as infile:
        mip_info = json.load(infile)
    mip_list = []
    info_list = []
    for g in mip_info:
        for m in mip_info[g]["mips"]:
            if (probe_set_keys is None) or (m in used_probes):
                minfo = mip_info[g]["mips"][m]["mip_dic"]["mip_information"]
                try:
                    cap_info = mip_info[g]["mips"][m]["capture_info"][
                        "captured_targets"]["must"]
                except KeyError:
                    cap_info = {}
                for probe in minfo:
                    for c in minfo[probe]["captures"]:
                        captured = []
                        if c in cap_info:
                            for cpt in cap_info[c]:
                                captured.append(cpt)
                        mip_list.append([g, m, probe, c, minfo[probe][
                            "SEQUENCE"], ",".join(captured)])
                minfo = mip_info[g]["mips"][m]["mip_info"]
                for c in minfo:
                    copy_dict = minfo[c]
                    try:
                        # lists and dicts cause problems in dataframe
                        # generation, so remove potential lists etc.
                        copy_dict.pop("alternative_arms")
                    except KeyError:
                        pass
                    try:
                        # lists and dicts cause problems in dataframe
                        # generation, so remove potential lists etc.
                        copy_dict.pop("haplotypes")
                    except KeyError:
                        pass
                    copy_dict["Gene"] = g
                    copy_dict["MIP"] = m
                    copy_dict["Copy"] = c
                    info_list.append(pd.DataFrame(copy_dict, index=[0]))
    mip_list = pd.DataFrame(mip_list)
    mip_list.columns = ["Gene", "MIP", "probe_id", "Copy", "Probe Sequence",
                        "Targets"]
    info_df = pd.concat(info_list, ignore_index=True)
    info_df = info_df.merge(mip_list)
    if output_file is not None:
        info_df.to_csv(output_file, index=False)
    return info_df


def get_probe_call_info(probe_set_keys=None, output_file=None,
                        probe_sets_file=None, call_info_file=None):
    used_probes = set()
    if probe_sets_file is None:
        probe_sets_file = "/opt/project_resources/mip_ids/probe_sets.json"
    if call_info_file is None:
        call_info_file = "/opt/project_resources/mip_ids/call_info.json"
    if probe_set_keys is not None:
        for psk in probe_set_keys:
            with open(probe_sets_file) as infile:
                used_probes.update(json.load(infile)[psk])

    with open(call_info_file) as infile:
        mip_info = json.load(infile)
    info_list = []
    for g in mip_info:
        for m in mip_info[g]:
            if (probe_set_keys is None) or (m in used_probes):
                minfo = mip_info[g][m]["copies"]
                for c in minfo:
                    copy_dict = minfo[c]
                    copy_dict["Gene"] = g
                    copy_dict["MIP"] = m
                    copy_dict["Copy"] = c
                    try:
                        copy_dict.pop("genes")
                    except KeyError:
                        pass
                    info_list.append(pd.DataFrame(copy_dict, index=[0]))
    info_df = pd.concat(info_list, ignore_index=True)
    if output_file is not None:
        info_df.to_csv(output_file, index=False)
    return info_df


def merge_call_infos(file_list, output_file=None):
    """Merge multiple call_info files in json format.

    Call info files are hierarchically organized into
    {gene: "mips": {mip1: ..}} structure.
    If two files have overlapping gene names, preceeding info file's gene
    will be updated with the following file. If mip names are shared within
    genes, the mips from the latter file will replace the preceeding file's
    mip.

    Save the merged call info file to specified file, if given.
    """
    info = {}
    for f in file_list:
        with open(f) as infile:
            temp_info = json.load(infile)
        for g in temp_info:
            try:
                info[g].update(temp_info[g])
            except KeyError:
                info[g] = temp_info[g]
    if output_file is not None:
        with open(output_file, "w") as outfile:
            json.dump(info, outfile, indent=1)
    return info


def merge_mip_infos(file_list, output_file=None):
    """Merge multiple mip_info files in json format.

    Mip info files are hierarchically organized into
    {gene: "mips": {mip1: ..}} structure.
    If two files have overlapping gene names, preceeding info file's gene
    will be updated with the following file. If mip names are shared within
    genes, the mips from the latter file will replace the preceeding file's
    mip.

    Save the merged mip info file to specified file, if given.
    """
    info = {}
    for f in file_list:
        with open(f) as infile:
            temp_info = json.load(infile)
        for g in temp_info:
            try:
                info[g]["mips"].update(temp_info[g]["mips"])
            except KeyError:
                info[g] = {"mips": temp_info[g]["mips"]}
    if output_file is not None:
        with open(output_file, "w") as outfile:
            json.dump(info, outfile, indent=1)
    return info


def update_probe_sets(
        mipset_table="/opt/project_resources/mip_ids/mipsets.csv",
        mipset_json="/opt/project_resources/mip_ids/probe_sets.json"):
    mipsets = pd.read_csv(mipset_table)
    mipset_list = mipsets.to_dict(orient="list")
    mipset_dict = {}
    for mipset in mipset_list:
        mlist = mipset_list[mipset]
        mipset_dict[mipset] = [m for m in mlist if not pd.isnull(m)]
    with open(mipset_json, "w") as outfile:
        json.dump(mipset_dict, outfile, indent=1)
    return


def generate_mip_arms_file(probe_set_key, probe_sets_file=None,
                           mip_info_file=None):
    if probe_sets_file is None:
        probe_sets_file = "/opt/project_resources/mip_ids/mipsets.csv"
    mipsets = pd.read_csv(probe_sets_file)
    mipset_list = mipsets.to_dict(orient="list")[probe_set_key]
    mipset_list = [m for m in mipset_list if not pd.isnull(m)]
    mip_arms_file = mipset_list[0]
    used_probes = set(mipset_list[1:])

    if mip_info_file is None:
        mip_info_file = "/opt/project_resources/mip_ids/mip_info.json"
    with open(mip_info_file) as infile:
        mip_info = json.load(infile)
    mip_list = []
    for g in mip_info:
        for m in mip_info[g]["mips"]:
            if m in used_probes:
                minfo = mip_info[g]["mips"][m]["mip_dic"]["mip_information"]
                for probe in minfo:
                    probe_seq = minfo[probe]["SEQUENCE"]
                    nlocs = []
                    nfound = False
                    for i in range(len(probe_seq)):
                        if probe_seq[i] == "N":
                            if not nfound:
                                nlocs.append(i)
                                nfound = True
                                continue
                        elif nfound:
                            nlocs.append(i)
                            nfound = False
                    if len(nlocs) != 4:
                        print(("Two groups of Ns were not found for MIP "
                               "{}. MIP arms will not be created for it."
                               ).format(m))
                    lig = mip.reverse_complement(probe_seq[:nlocs[0]])
                    ext = probe_seq[nlocs[-1]:]
                    lig_umi_len = nlocs[1] - nlocs[0]
                    ext_umi_len = nlocs[3] - nlocs[2]
                    temp = [g, m, m + "_" + probe, ext, lig, ext_umi_len,
                            lig_umi_len, probe_set_key]
                    mip_list.append(temp)

    mip_arms = pd.DataFrame(mip_list)
    mip_arms.columns = ["gene_name", "mip_family", "mip_id", "extension_arm",
                        "ligation_arm", "extension_barcode_length",
                        "ligation_barcode_length", "mipset"]

    mip_arms.to_csv(os.path.join("/opt/project_resources/mip_ids",
                                 mip_arms_file), index=False, sep="\t")
    update_probe_sets()
    return

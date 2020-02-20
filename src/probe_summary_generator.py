import json
import pandas as pd


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

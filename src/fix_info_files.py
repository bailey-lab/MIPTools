import json
import subprocess
import mip_functions as mip


def fix(call_file, mip_file, species, verbose=True):
    with open(call_file) as infile:
        call_info = json.load(infile)
    with open(mip_file) as infile:
        mip_info = json.load(infile)

    pairs_keys = ['capture_size', 'mip_start', 'mip_end', 'ligation_start',
                  'ligation_end', 'extension_start', 'extension_end',
                  'capture_start', 'capture_end']
    cap_keys = ["begin", "end", "copy_begin", "copy_end"]
    cap_copy_keys = ["begin", "end"]

    for g in mip_info:
        for m in mip_info[g]["mips"]:
            md = mip_info[g]["mips"][m]["mip_dic"]
            pairs = md["pairs"]
            for c in pairs:
                for k in pairs_keys:
                    pairs[c][k] = int(pairs[c][k])
                for k in pairs_keys:
                    if k != "capture_end":
                        if pairs[c][k] == pairs[c]["capture_end"]:
                            if verbose:
                                print(m, c, k)
                            pairs[c]["capture_end"] -= 1
                chrom = pairs[c]["chrom"]
                try:
                    ck = pairs[c]["capture_key"]
                except KeyError:
                    ck = "none"
                capture_key = mip.create_region(
                    chrom, pairs[c]["capture_start"], pairs[c]["capture_end"])
                if ck != capture_key:
                    if verbose:
                        print(ck, capture_key, m, c)
                    pairs[c]["capture_sequence"] = mip.get_sequence(
                        capture_key, species)
                    pairs[c]["capture_key"] = capture_key
                mip_info[g]["mips"][m]["mip_info"][c] = pairs[c]
            for pt in ["extension_primer_information",
                       "ligation_primer_information"]:
                prim = md[pt]["PARALOG_COORDINATES"]
                for c in prim:
                    prim[c]["GENOMIC_START"] = int(prim[c]["GENOMIC_START"])
                    prim[c]["GENOMIC_END"] = int(prim[c]["GENOMIC_END"])
                    ori = prim[c]["ORI"]
                    chrom = prim[c]["CHR"]
                    if ori == "forward":
                        prim_key = mip.create_region(
                            chrom, prim[c]["GENOMIC_START"],
                            prim[c]["GENOMIC_END"])
                        prim[c]["SEQUENCE"] = mip.get_sequence(prim_key,
                                                               species)
                    elif ori == "reverse":
                        prim_key = mip.create_region(
                            chrom, prim[c]["GENOMIC_END"],
                            prim[c]["GENOMIC_START"])
                        prim[c]["SEQUENCE"] = mip.reverse_complement(
                            mip.get_sequence(prim_key, species))
                    prim[c]["KEY"] = prim_key
                    if prim[c]["SEQUENCE"] == "":
                        if verbose:
                            print("Empty sequence", m, c)
            capture_info = mip_info[g]["mips"][m]["capture_info"][
                "captured_targets"]
            for t in capture_info:
                for ck in cap_keys:
                    try:
                        capture_info[t][ck] = int(capture_info[t][ck])
                    except KeyError:
                        continue
                for copyname in capture_info[t]["copies"]:
                    for ck in cap_copy_keys:
                        capture_info[t]["copies"][copyname][ck] = int(
                            capture_info[t]["copies"][copyname][ck])

    for g in call_info:
        for m in call_info[g]:
            snps = call_info[g][m]["snps"]
            for s in snps:
                snps[s]["begin"] = int(snps[s]["begin"])
                snps[s]["end"] = int(snps[s]["end"])
                cop = snps[s]["copies"]
                for c in cop:
                    cop[c]["begin"] = int(cop[c]["begin"])
                    cop[c]["end"] = int(cop[c]["end"])
            copies = call_info[g][m]["copies"]
            for c in copies:
                md = mip_info[g]["mips"][m]["mip_info"][c]
                for k in copies[c]:
                    if k in md:
                        copies[c][k] = md[k]

    subprocess.call(["scp", mip_file, mip_file + ".bak2"])
    with open(mip_file, "w") as outfile:
        json.dump(mip_info, outfile, indent=1)

    subprocess.call(["scp", call_file, call_file + ".bak2"])
    with open(call_file, "w") as outfile:
        json.dump(call_info, outfile, indent=1)

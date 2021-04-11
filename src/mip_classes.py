
"""
Created on Tue Apr 29 16:50:42 2014

@author: oz
"""
import json
import pickle
import os
import random
from operator import itemgetter
import mip_functions as mip
import copy
import numpy as np
import pandas as pd
import traceback
import primer3
print("Classes reloading.")


class Locus:
    """
    Locus class provides a construct for creating objects representing
    any DNA locus that is to be used for MIP designs.
    Multiple loci can and should be represented by a single
    object when there is high sequence similarity, i.e. paralogus sequences.
    A locus must have at least one sequence that is designated as the
    reference when multiple loci represented in a single Locus object.
    Locus objects will be instantiated objects of the Paralog class.
    """

    def __init__(self, paralog, rinfo_file):
        # connect the paralog object creating the locus so that the locus
        # "knows" what paralog object it belongs to
        self.paralog = paralog
        # rifo file has capture specifics for the mips
        self.rinfo_file = rinfo_file
        # convert rinfo file to dictionary
        self.rinfo = self.rinfo_parser(self.rinfo_file)
        # get species and host species from rinfo dict
        self.species = self.rinfo["SPECIES"]
        # host species can be none, or any contaminating species
        # it is created for human, when designing parasite mips
        self.host = self.rinfo["HOST_SPECIES"]
        # get the design directory relative to root
        self.design_dir = self.rinfo["DESIGN_DIR"]
        # get file location dictionary for locations of static files
        # such as genome fasta files, snp files etc.
        self.file_locations = mip.get_file_locations()
        # give the object a name
        self.paralog_name = self.rinfo["NAME"]
        self.segment_name = "S0"
        self.cwd = self.paralog.cwd + self.segment_name + "/"
        if not os.path.exists(self.cwd):
            os.makedirs(self.cwd)
        self.resource_dir = self.paralog.resource_dir
        self.regions = self.paralog.regions
        self.alignments = self.paralog.alignments
        self.target_genes = self.paralog.target_genes
        self.segment_dic = self.paralog.segment_dic
        self.begin = self.paralog.begin
        self.end = self.paralog.end
        self.chrom = self.paralog.chrom
        self.extras_dic = self.paralog.extras_dic
        self.copies = self.paralog.copies
        self.zcoordinates = self.paralog.zcoordinates
        self.pcoordinates = self.paralog.pcoordinates
        self.pdiffs = self.paralog.pdiffs
        self.extra_snps = self.paralog.extra_snps
        self.snps = self.paralog.snps
        self.must = self.paralog.must
        self.exons = self.paralog.exons
        self.get_variants()
        self.subregions = {}
        self.designed = False
        self.failed = False

    def get_segment_dic(self):
        """ Take a rinfo dictionary and extract the segments in it.
        Return a dictionary in the form:
        {'S0': ['C0', 'C1', 'C2', 'S0']}
        """
        # get region information
        regions = self.regions
        # get alignment information
        alignments = self.alignments
        copies = []
        extras = []
        # get copies from alignment dict
        for k in alignments:
            key_id = k.split("_")[-1]
            copy_num = int(key_id[1:])
            # sequences that originate from reference genome have ids like C0,
            # C1. Extra sequences that has no known reference genome positions
            # have ids such as X0, X1.
            if key_id.startswith("C"):
                copies.append(copy_num)
            elif key_id.startswith("X"):
                extras.append(copy_num)
        segment_copies = []
        copies.sort()
        extras.sort()
        segment_dic = {"S0": {}}
        seg = segment_dic["S0"]
        for c in copies:
            copy_id = "C" + str(c)
            segment_copies.append(copy_id)
            ori = alignments[self.paralog_name + "_" + copy_id]["strand2"]
            if ori == "+":
                orientation = "forward"
            elif ori == "-":
                orientation = "reverse"
            seg[copy_id] = {"chrom": regions[c][0],
                            "begin": regions[c][1],
                            "end": regions[c][2],
                            "length": regions[c][2]-regions[c][1] + 1,
                            "orientation": orientation,
                            "copy_key": self.paralog_name + "_" + copy_id}
        extras_dic = {}
        for e in extras:
            copy_id = "X" + str(e)
            ori = alignments[self.paralog_name + "_" + copy_id]["strand2"]
            if ori == "+":
                orientation = "forward"
            elif ori == "-":
                orientation = "reverse"
            extras_dic[copy_id] = {"orientation": orientation,
                                   "copy_key": self.paralog_name
                                   + "_" + copy_id}
        return {"segment_dic": segment_dic,
                "extras_dic": extras_dic,
                "segment_copies": segment_copies}

    def segment_converter(self):
        """
        Convert 0-offset coordinates between paralogs for a given segment.
        """
        segment_dic = self.segment_dic["S0"]
        coordinates = {"C0": {"chromosomes": {}}}
        for c in segment_dic:
            coordinates["C0"]["chromosomes"][c] = segment_dic[c]["chrom"]
            co = self.alignments[segment_dic[c]["copy_key"]]["coordinates"]
            rev_co = self.alignments[segment_dic[c]["copy_key"]][
                "reverse_coordinates"
            ]
            ref_start = segment_dic["C0"]["begin"] - 1
            copy_start = segment_dic[c]["begin"] - 1
            if c != "C0":
                coordinates["C0"][c] = {}
                coordinates[c] = {}
                for i in co:
                    # convert .5 coordinates to lower whole number
                    j = int(np.floor(i))
                    coordinates["C0"][c][j + ref_start] = int(
                        np.floor(co[i])
                    ) + copy_start
                for i in rev_co:
                    j = int(np.floor(i))
                    coordinates[c][j + copy_start] = int(
                        np.floor(rev_co[i])
                    ) + ref_start
        return coordinates

    def base_converter(self):
        """ Convert 0-offset pcoordinates to 1-offset."""
        zcoordinates = self.zcoordinates
        segment_dic = self.segment_dic["S0"]
        pcoordinates = {}
        for c in segment_dic:
            pcoordinates[c] = {}
            if c == "C0":
                for co in zcoordinates[c]:
                    if co != "chromosomes":
                        pcoordinates[c][co] = {}
                        for k in zcoordinates[c][co]:
                            pcoordinates[c][co][k+1] = (zcoordinates
                                                        [c][co][k] + 1)
                    else:
                        pcoordinates[c][co] = zcoordinates[c][co]
            else:
                pcoordinates[c] = {}
                for k in zcoordinates[c]:
                    pcoordinates[c][k+1] = zcoordinates[c][k] + 1
        return pcoordinates

    def get_pdiffs(self):
        """ Return 1-offset pdiffs and extra snps. """
        alignments = self.alignments
        segment_dic = self.segment_dic["S0"]
        diffs = {}
        for c in segment_dic:
            copy_diffs = diffs[c] = {}
            copy_key = segment_dic[c]["copy_key"]
            snps = alignments[copy_key]["snps"]
            for s in snps:
                query_begin = int(s.split("-")[0])
                query_end = int(s.split("-")[1])
                # get zero-offset begin coordinates for regions
                genomic_ref_begin = segment_dic["C0"]["begin"] - 1
                ref_chrom = segment_dic["C0"]["chrom"]
                begin = genomic_ref_begin + query_begin
                end = genomic_ref_begin + query_end
                ref_ori = snps[s]["query_orientation"]
                ref_base = snps[s]["query_base"]
                target_begin = snps[s]["target_begin"]
                target_end = snps[s]["target_end"]
                genomic_copy_begin = segment_dic[c]["begin"] - 1
                copy_chrom = segment_dic[c]["chrom"]
                copy_begin = genomic_copy_begin + target_begin
                copy_end = genomic_copy_begin + target_end
                copy_base = snps[s]["target_base"]
                diff_ori = "forward"
                if ref_ori == "-":
                    diff_ori = "reverse"
                    ref_base = mip.reverse_complement(ref_base)
                    copy_base = mip.reverse_complement(copy_base)
                # Insertions and single changes will be represented by a single
                # diff. Insertions will have the smaller coordinate as the
                # 1-offset coordinate so no changes to end and begin
                # coordinates for the reference but the insertion begin
                # coordinate will be incremented by 1.
                if end == begin:
                    d = {"chrom": ref_chrom,
                         "begin": begin,
                         "end": end + 1,
                         "type": "insertion",
                         "base": ref_base,
                         "orientation": diff_ori,
                         "copy_chrom": copy_chrom,
                         "copy_begin": copy_begin + 1,
                         "copy_end": copy_end,
                         "copy_base": copy_base,
                         "size_difference": copy_end - copy_begin}
                elif end - begin == 1:
                    d = {"chrom": ref_chrom,
                         "begin": begin + 1,
                         "end": end,
                         "type": "SNV",
                         "base": ref_base,
                         "orientation": diff_ori,
                         "copy_chrom": copy_chrom,
                         "copy_begin": copy_begin + 1,
                         "copy_end": copy_end,
                         "copy_base": copy_base,
                         "size_difference": 0}
                # multiple changes will be split
                # if it is not a deletion
                elif "-" not in copy_base:
                    for i in range(end - begin):
                        if diff_ori == "forward":
                            d = {"chrom": ref_chrom,
                                 "begin": begin + i + 1,
                                 "end": begin + i + 1,
                                 "type": "SNV",
                                 "base": ref_base[i],
                                 "orientation": diff_ori,
                                 "copy_chrom": copy_chrom,
                                 "copy_begin": copy_begin + i + 1,
                                 "copy_end": copy_begin + i + 1,
                                 "copy_base": copy_base[i],
                                 "size_difference": 0}
                        else:
                            d = {"chrom": ref_chrom,
                                 "begin": begin + i + 1,
                                 "end": begin + i + 1,
                                 "type": "SNV",
                                 "base": ref_base[i],
                                 "orientation": diff_ori,
                                 "copy_chrom": copy_chrom,
                                 "copy_begin": copy_end - i,
                                 "copy_end": copy_end - i,
                                 "copy_base": copy_base[-i-1],
                                 "size_difference": 0}
                else:  # if deletion
                    if diff_ori == "forward":
                        d = {"chrom": ref_chrom,
                             "begin": begin + 1,
                             "end": end,
                             "type": "deletion",
                             "base": ref_base,
                             "orientation": diff_ori,
                             "copy_chrom": copy_chrom,
                             "copy_begin": copy_begin,
                             "copy_end": copy_begin,
                             "copy_base": copy_base,
                             "size_difference": begin - end}
                    else:
                        d = {"chrom": ref_chrom,
                             "begin": begin + 1,
                             "end": end,
                             "type": "deletion",
                             "base": ref_base,
                             "orientation": diff_ori,
                             "copy_chrom": copy_chrom,
                             "copy_begin": copy_end,
                             "copy_end": copy_end,
                             "copy_base": copy_base,
                             "size_difference": begin - end}
                copy_diffs[d["copy_chrom"]
                           + ":" + str(d["copy_begin"])
                           + d["copy_base"] + ":"
                           + ref_base] = d

        extras_dic = self.extras_dic
        extra_snps = {}
        for c in extras_dic:
            copy_key = extras_dic[c]["copy_key"]
            snps = alignments[copy_key]["snps"]
            copy_snps = extra_snps[c] = {}
            for s in snps:
                query_begin = int(s.split("-")[0])
                query_end = int(s.split("-")[1])
                genomic_ref_begin = segment_dic["C0"]["begin"] - 1
                ref_chrom = segment_dic["C0"]["chrom"]
                begin = genomic_ref_begin + query_begin
                end = genomic_ref_begin + query_end
                ref_ori = snps[s]["query_orientation"]
                ref_base = snps[s]["query_base"]
                copy_base = snps[s]["target_base"]
                diff_ori = "forward"
                if ref_ori == "-":
                    diff_ori = "reverse"
                    ref_base = mip.reverse_complement(ref_base)
                    copy_base = mip.reverse_complement(copy_base)
                # insertions and single changes will be represented by a single
                # diff
                if end == begin:
                    d = {"chrom": ref_chrom,
                         "begin": begin,
                         "end": end,
                         "type": "insertion",
                         "base": ref_base,
                         "orientation": diff_ori,
                         "copy_base": copy_base,
                         "size_difference": len(copy_base)}
                elif end - begin == 1:
                    d = {"chrom": ref_chrom,
                         "begin": begin + 1,
                         "end": end,
                         "type": "SNV",
                         "base": ref_base,
                         "orientation": diff_ori,
                         "copy_base": copy_base,
                         "size_difference": 0}
                elif "-" not in copy_base:
                    for i in range(end - begin):
                        if diff_ori == "forward":
                            d = {"chrom": ref_chrom,
                                 "begin": begin + i + 1,
                                 "end": begin + i + 1,
                                 "type": "SNV",
                                 "base": ref_base[i],
                                 "orientation": diff_ori,
                                 "copy_base": copy_base[i],
                                 "size_difference": 0}
                        else:
                            d = {"chrom": ref_chrom,
                                 "begin": begin + i + 1,
                                 "end": begin + i + 1,
                                 "type": "SNV",
                                 "base": ref_base[i],
                                 "orientation": diff_ori,
                                 "copy_base": copy_base[-i-1],
                                 "size_difference": 0}
                else:
                    if diff_ori == "forward":
                        d = {"chrom": ref_chrom,
                             "begin": begin + 1,
                             "end": end,
                             "type": "deletion",
                             "base": ref_base,
                             "orientation": diff_ori,
                             "copy_base": copy_base,
                             "size_difference": begin - end}
                    else:
                        d = {"chrom": ref_chrom,
                             "begin": begin + 1,
                             "end": end,
                             "type": "deletion",
                             "base": ref_base,
                             "orientation": diff_ori,
                             "copy_base": copy_base,
                             "size_difference": begin - end}
                copy_snps[c + ":" + str(d["begin"])
                          + d["copy_base"] + ":"
                          + ref_base] = d
        return {"pdiffs": diffs,
                "extra_snps": extra_snps}

    def get_snps_from_table(self):
        """Return SNPs from a UCSC style SNP table.

        Return SNPs within a Locus object using a UCSC style SNP table.
        This is deprecated in favor of using get_snps function which takes
        a vcf file instead.
        """
        pcoordinates = self.pcoordinates
        chrom = pcoordinates["C0"]["chromosomes"]["C0"]
        segment_dic = self.segment_dic["S0"]
        snps = {}
        snp_file = mip.get_file_locations()[self.species]["snps"]
        for c in segment_dic:
            snps[c] = {}
            copy_chrom = segment_dic[c]["chrom"]
            begin = segment_dic[c]["begin"]
            end = segment_dic[c]["end"]
            copy_key = copy_chrom + ":" + str(begin) + "-" + str(end)
            snp_list = mip.get_snps(copy_key, snp_file)
            for s in snp_list:
                alleles = s[22].split(",")
                allele_counts = s[23].split(",")
                if len(alleles) > 1:
                    alleles.pop(-1)
                    allele_counts.pop(-1)
                    allele_counts = list(map(int, list(map(float,
                                                           allele_counts))))
                    reference_allele_count = allele_counts[0]
                    total_allele_count = sum(allele_counts)
                    maf = float(total_allele_count
                                - reference_allele_count)/total_allele_count
                else:
                    maf = 0
                snp_ori = "forward"
                ref_base = s[8]
                observed = s[9].split("/")
                if s[6] == "-":
                    snp_ori = "reverse"
                    ref_base = mip.reverse_complement(ref_base)
                    observed = list(map(mip.reverse_complement, observed))
                d = {"copy_chrom": s[1],
                     "copy_begin": int(s[2]),
                     "copy_end": int(s[3]),
                     "snp_id": s[4],
                     "orientation": snp_ori,
                     "copy_base": ref_base,
                     "observed": observed,
                     "class": s[11],
                     "snp_function": s[15].split(","),
                     "exceptions": s[18],
                     "alleles": alleles,
                     "allele_counts": allele_counts,
                     "minor_allele_freq": maf}
                # snps from ucsc are 0-offset. change to 1-offset
                if d["copy_begin"] != d["copy_end"]:
                    d["copy_begin"] += 1
                if c != "C0":
                    try:
                        begin = pcoordinates[c][d["copy_begin"]]
                        end = pcoordinates[c][d["copy_end"]]
                    except KeyError:
                        # if the snp is in a region of the copy that is not
                        # aligned to the ref, we will not add the SNP
                        continue
                else:
                    begin = d["copy_begin"]
                    end = d["copy_end"]
                d["chrom"] = chrom
                if begin > end:
                    d["begin"] = end
                    d["end"] = begin
                else:
                    d["begin"] = begin
                    d["end"] = end
                snps[c][d["chrom"] + ":" + str(d["begin"]) + "-" + str(
                    d["end"])] = d
        return snps

    def get_snps(self):
        pcoordinates = self.pcoordinates
        chrom = pcoordinates["C0"]["chromosomes"]["C0"]
        segment_dic = self.segment_dic["S0"]
        snps = {}
        try:
            af_name = self.rinfo["CAPTURE"]["S0"]["allele_frequency_name"]
        except KeyError:
            af_name = "AF"
        try:
            snp_file = mip.get_file_locations()[self.species]["snps"]
        except KeyError:
            return {}
        for c in segment_dic:
            snps[c] = {}
            copy_chrom = segment_dic[c]["chrom"]
            begin = segment_dic[c]["begin"]
            end = segment_dic[c]["end"]
            copy_key = copy_chrom + ":" + str(begin) + "-" + str(end)
            snp_list = mip.get_vcf_snps(copy_key, snp_file)
            for s in snp_list:
                schr = s[0]
                spos = int(s[1])
                sref = s[3]
                salts = s[4].split(",")
                sinfo = s[7].split(";")
                sdict = {}
                for inf in sinfo:
                    try:
                        split_inf = inf.split("=")
                        sdict[split_inf[0]] = split_inf[1].split(",")
                    except IndexError:
                        sdict[inf] = True
                try:
                    sdict["AF"] = ["0" if af == "." else af
                                   for af in sdict[af_name]]
                except KeyError:
                    pass
                d = {"copy_chrom": schr,
                     "copy_begin": spos,
                     "copy_base": sref,
                     "alleles": salts,
                     "info": sdict}
                if c != "C0":
                    try:
                        begin = pcoordinates[c][d["copy_begin"]]
                    except KeyError:
                        # if the snp is in a region of the copy
                        # that is not aligned to the ref, we will not add
                        # the SNP
                        continue
                else:
                    begin = d["copy_begin"]
                d["chrom"] = chrom
                d["begin"] = begin
                snp_key = schr + ":" + str(spos) + ":" + sref + ":" + s[4]
                snps[c][snp_key] = d
        return snps

    def get_must(self):
        must_dic = self.rinfo["MUST"]
        pcoordinates = self.pcoordinates
        chrom = pcoordinates["C0"]["chromosomes"]["C0"]
        must = {}
        for m in must_dic:
            c = must_dic[m]["copy"]
            try:
                copy_must = must[c]
            except KeyError:
                copy_must = must[c] = {}
            try:
                sz_diff = int(must_dic[m]["size_difference"])
            except (KeyError, ValueError):
                sz_diff = 0
            d = {"name": m,
                 "copy_begin": int(must_dic[m]["begin"]),
                 "copy_end": int(must_dic[m]["end"]),
                 "copy_chrom": must_dic[m]["chrom"],
                 "size_difference": sz_diff}
            if d["size_difference"] == 0:
                d["type"] = "SNV"
            elif d["size_difference"] < 0:
                d["type"] = "deletion"
            elif d["size_difference"] > 0:
                d["type"] = "insertion"
            if c != "C0":
                    begin = pcoordinates[c][d["copy_begin"]]
                    end = pcoordinates[c][d["copy_end"]]
            else:
                begin = d["copy_begin"]
                end = d["copy_end"]
            d["chrom"] = chrom
            if begin > end:
                d["begin"] = end
                d["end"] = begin
            else:
                d["begin"] = begin
                d["end"] = end
            copy_must[d["name"]] = d
        return must

    def get_variants(self):
        columns = ["chrom", "begin", "end", "base",
                   "copy_chrom", "copy_begin", "copy_end", "copy_base",
                   "size_difference", "orientation",
                   "AF", "type", "variant_type", "copy_id"]
        coordinate_converter = self.pcoordinates
        chrom = coordinate_converter["C0"]["chromosomes"]["C0"]
        variant_dicts = {"pdiffs": self.pdiffs,
                         "snps": self.snps,
                         "extra_snps": self.extra_snps,
                         "must": self.must}
        all_variants = []
        for variant_type in variant_dicts:
            vd = variant_dicts[variant_type]
            for copy_id in vd:
                for s in vd[copy_id]:
                    snp = vd[copy_id][s]
                    if variant_type != "snps":
                        var = []
                        for col in columns[:-2]:
                            try:
                                var.append(snp[col])
                            except KeyError:
                                var.append(np.nan)
                        var.append(variant_type)
                        var.append(copy_id)
                        all_variants.append(var)
                    else:
                        copy_chrom = snp["copy_chrom"]
                        copy_base = snp["copy_base"]
                        copy_begin = snp["copy_begin"]
                        alleles = snp["alleles"]
                        try:
                            allele_freqs = list(map(float, snp["info"]["AF"]))
                            af_start_index = int(self.rinfo["CAPTURE"]["S0"][
                                "af_start_index"])
                            allele_freqs = allele_freqs[af_start_index:]
                        except KeyError:
                            allele_freqs = [1 for i in range(len(alleles))]
                        for allele_index in range(len(alleles)):
                            split_snp = {"copy_chrom": copy_chrom,
                                         "chrom": chrom}
                            a = alleles[allele_index]
                            a_freq = allele_freqs[allele_index]
                            # check if allele is insertion/deletion or snp
                            a_len = len(a) - len(copy_base)
                            if a_len == 0:  # snv allele
                                split_snp["copy_begin"] = copy_begin
                                split_snp["copy_end"] = copy_begin
                                split_snp["copy_base"] = copy_base
                                if copy_id != "C0":
                                    ref_coordinates = [coordinate_converter[
                                        copy_id][copy_begin],
                                                       coordinate_converter[
                                        copy_id][copy_begin]
                                                      ]
                                    begin = min(ref_coordinates)
                                    end = max(ref_coordinates)
                                else:
                                    begin = copy_begin
                                    end = copy_begin
                                split_snp["begin"] = begin
                                split_snp["end"] = end
                                split_snp["base"] = a
                                split_snp["size_difference"] = 0
                                split_snp["AF"] = a_freq
                                split_snp["type"] = "SNV"
                            elif a_len < 0:  # deletion allele
                                # deletions in vcf are represented so:
                                # 100 ATTT AT meaning the bases at
                                # positions 102 and 103 are deleted
                                # the position we have as begin is 100
                                copy_deletion_begin = (copy_begin
                                                       + len(a) - 1)
                                # this corresponds to the last undeleted base
                                # before deletion, 101
                                copy_deletion_end = (copy_begin
                                                     + len(copy_base))
                                # this corresponds to the first undeleted base
                                # after deletion, 104
                                if copy_id != "C0":
                                    deletion_coord = [
                                        coordinate_converter[copy_id]
                                        [copy_deletion_begin],
                                        coordinate_converter[copy_id]
                                        [copy_deletion_end]
                                    ]
                                    deletion_begin = min(deletion_coord)
                                    deletion_end = max(deletion_coord)
                                else:
                                    deletion_begin = copy_deletion_begin
                                    deletion_end = copy_deletion_end
                                split_snp["copy_begin"] = (copy_deletion_begin
                                                           + 1)
                                split_snp["copy_end"] = copy_deletion_end - 1
                                split_snp["copy_base"] = copy_base[len(a):]
                                split_snp["begin"] = deletion_begin + 1
                                split_snp["end"] = deletion_end - 1
                                split_snp["base"] = -a_len * "-"
                                split_snp["size_difference"] = a_len
                                split_snp["AF"] = a_freq
                                split_snp["type"] = "deletion"
                            elif a_len > 0:  # insertion allele
                                # insertions are represented so: 100 AAA AAATT
                                # showing a 2 bp insertion between positions
                                # 102 and 103, the begin position we have is
                                # 100
                                copy_insertion_begin = copy_begin + len(
                                    copy_base) - 1
                                # this position corresponds to the position
                                # before the insertion, 102
                                copy_insertion_end = copy_insertion_begin + 1
                                if copy_id != "C0":
                                    insertion_coord = [
                                        coordinate_converter[copy_id]
                                        [copy_insertion_begin],
                                        coordinate_converter[copy_id]
                                        [copy_insertion_end]
                                    ]
                                    insertion_begin = min(insertion_coord)
                                else:
                                    insertion_begin = copy_insertion_begin
                                split_snp["copy_begin"] = copy_insertion_begin
                                split_snp["copy_end"] = copy_insertion_begin
                                split_snp["copy_base"] = copy_base[-1]
                                split_snp["begin"] = insertion_begin
                                split_snp["end"] = insertion_begin
                                split_snp["base"] = a[-a_len - 1:]
                                split_snp["size_difference"] = a_len
                                split_snp["AF"] = a_freq
                                split_snp["type"] = "insertion"
                            var = []
                            for col in columns[:-2]:
                                try:
                                    var.append(split_snp[col])
                                except KeyError:
                                    var.append(np.nan)
                            var.append(variant_type)
                            var.append(copy_id)
                            all_variants.append(var)
        self.all_variants = pd.DataFrame(all_variants, columns=columns)
        snp_df = self.all_variants.loc[self.all_variants["variant_type"]
                                       == "snps"]
        # we need per position variants for masking purposes.
        # i.e. deletions represented as AAA->--- should be A->- in 3 positions.
        # a small function to divide SNPs is provided by split_snps

        def split_snps(r):
            r_chrom = r["chrom"]
            r_begin = r["begin"]
            r_end = r["end"]
            af = r["AF"]
            res = []
            for i in range(r_begin, r_end + 1):
                res.append([r_chrom, i, af])
            return pd.DataFrame(res, columns=["chrom", "position", "AF"])
        split_result_list = []
        for i in snp_df.index:
            split_result_list.append(
                split_snps(snp_df.loc[i])
            )
        try:
            split = pd.concat(split_result_list, ignore_index=True)
        except ValueError:
            split = pd.DataFrame(split_result_list,
                                 columns=["chrom", "position", "AF"])
        # if there are no snps, skip groupby
        # to avoid losing chrom and position columns
        if not split.empty:
            split = split.groupby(["chrom", "position"]).aggregate(
                {"AF": sum}).reset_index()
        maf_filter = float(self.rinfo["CAPTURE"]["S0"]["maf_for_arm_design"])
        split = split.loc[split["AF"] >= maf_filter]
        self.split_snps = split

        # A similar split is needed on non-snp variants (pdiffs, etc.)
        # we don't have allele counts for these, so it'll be only
        # chromosomes and positions in the split diffs.
        diff_df = self.all_variants.loc[self.all_variants["variant_type"]
                                        != "snps"]
        split_diff_list = [
            pd.DataFrame([[diff_df.loc[diff_index]["chrom"],
                           i] for i in range(
                diff_df.loc[diff_index]["begin"],
                diff_df.loc[diff_index]["end"] + 1)],
                     columns=["chrom", "position"])
            for diff_index in diff_df.index]
        try:
            split_diff_df = pd.concat(split_diff_list,
                                      ignore_index=True).drop_duplicates()
        except ValueError:
            split_diff_df = pd.DataFrame(split_diff_list,
                                         columns=["chrom", "position"])
        self.split_pdiffs = split_diff_df.reset_index(drop=True)

        ######################################################################
        # prepare a dataframe with potential size variation in the region
        # this will be used in picking primer pairs and not in masking
        # for each position we'll try to find the largest indel which
        # has a minor allele frequency of >= insertion_filter
        ######################################################################

        insertion_df = self.all_variants.loc[
            (self.all_variants["type"] == "insertion")
            & (self.all_variants["variant_type"] == "snps")]
        insertion_filter = float(self.rinfo["CAPTURE"]["S0"]["maf_for_indels"])
        filt_insertions = insertion_df.loc[(insertion_df["AF"]
                                           >= insertion_filter)]
        if not filt_insertions.empty:
            max_sizes = filt_insertions.groupby(
                ["copy_chrom", "copy_begin"])["size_difference"].max()
            max_sizes = max_sizes.reset_index().rename(
                columns={"size_difference": "max_size"})
            insertion_df = insertion_df.merge(max_sizes, how="inner")
            # Collapse on position and select the largest insertion
            collapsed_insertions = insertion_df.groupby(["copy_chrom",
                                                         "copy_begin"]).agg(
                {"copy_end": "first", "AF": sum, "max_size": max}
            ).reset_index()
            # filter for maf
            filt_insertions = collapsed_insertions.loc[collapsed_insertions[
                "AF"] >= maf_filter]
        self.insertions = filt_insertions
        return

    def get_exons(self):
        copies = self.segment_dic["S0"]
        for c in copies:
            copy_key = (copies[c]["chrom"] + ":" + str(copies[c]["begin"])
                        + "-" + str(copies[c]["end"]))
            copies[c]["genes"] = mip.get_region_exons(copy_key, self.species)
        return

    def get_projected_exons(self):
        """
        Create a dictionary of segment exons projected on reference copy.

        keys: copy_names: [list of exons of this copy projected on reference
        copy coordinates], "merged": [list of all exons in segment projected on
        ref copy, overlapping exons merged]}
        """
        copies = self.segment_dic["S0"]
        try:
            reference_exons = copies["C0"]["genes"]["exons"]
        except KeyError:
            reference_exons = []
        reference_ori = copies["C0"]["orientation"]
        projected_exons = {}
        projected_exons["C0"] = [e for e in reference_exons if copies["C0"][
            "begin"] <= e[0] <= e[1] <= copies["C0"]["end"]]
        con = self.pcoordinates
        for c in copies:
            if c != "C0":
                try:
                    copy_exons = copies[c]["genes"]["exons"]
                    copy_exons = [e for e in copy_exons if copies[c]["begin"]
                                  <= e[0] <= e[1] <= copies[c]["end"]]
                except KeyError:
                    continue
                copy_ori = copies[c]["orientation"]
                # create a list of projected exons for the given copy
                projected_copy_exons = []
                if reference_ori == copy_ori:
                    for ce in copy_exons:
                        # add copy exon if it was aligned with ref
                        # otherwise we'll get a key error from
                        # trying to convert the exon coordinate to ref
                        try:
                            projected_exon_start = con[c][ce[0]]
                            projected_exon_end = con[c][ce[1]]
                            projected_copy_exons.append([projected_exon_start,
                                                         projected_exon_end])
                        except KeyError:
                            continue
                else:
                    for ce in reversed(copy_exons):
                        try:
                            projected_exon_start = con[c][ce[1]]
                            projected_exon_end = con[c][ce[0]]
                            projected_copy_exons.append([projected_exon_start,
                                                         projected_exon_end])
                        except KeyError:
                            continue
                projected_exons[c] = projected_copy_exons
        all_exons = []
        for p in projected_exons:
            all_exons.extend(projected_exons[p])
        projected_exons["merged"] = sorted(mip.merge_overlap(all_exons),
                                           key=itemgetter(0))
        return projected_exons

    def make_subregions(self):
        """
        Create subregions of the target region based on the capture purpose.

        Subregions are the main operational unit in the design
        pipeline, i.e. MIPs are designed for each subregion, MIPs of each
        subregion are "aware" of each other in terms of compatibility etc.
        For example, if capture purpose is to sequence exons only
        and this is specified in the capture type as "exons", then generate
        a subregion for each exon provided that exons are separated by enough
        distance; otherwise two exons may be included in a single subregion.
        """
        capture_info = self.rinfo["CAPTURE"]["S0"]
        capture_type = capture_info["capture_type"]
        flank = int(capture_info["flank"])
        # create a list of subregions
        subregions = []
        must = self.must
        pdiffs = self.pdiffs
        extra_snps = self.extra_snps
        target_diffs = capture_info["target_diffs"].split(",")
        # create targets attribute for the segment
        self.targets = capture_targets = {"must": must}
        if "pdiffs" in target_diffs:
            capture_targets["pdiffs"] = pdiffs
        if "extra_snps" in target_diffs:
            capture_targets["extra_snps"] = extra_snps
        if capture_type == "whole":
            subregions = [[self.begin + flank, self.end - flank]]
        elif capture_type == "exons":
            exons = copy.deepcopy(self.exons)
            exonic_subregions = mip.merge_overlap(exons, 2*flank)
            # check if all regions in "must capture" list are covered in exons
            uncaptured = copy.deepcopy(must)
            for e in exonic_subregions:
                for c in must:
                    for m in must[c]:
                        if (e[0] <= must[c][m]["begin"]
                                <= must[c][m]["end"] <= e[1]):
                            uncaptured[c].pop(m)
            for c in list(uncaptured.keys()):
                if len(uncaptured[c]) == 0:
                    uncaptured.pop(c)
            # if there are targets falling outside of exons, these will be
            # identified and the exons will be extended, or new intervals
            # for those will be added to include all targets.
            uncaptured_coordinates = []
            for c in uncaptured:
                for u in uncaptured[c]:
                    uncaptured_coordinates.append([uncaptured[c][u]["begin"],
                                                  uncaptured[c][u]["end"]])
            complete_targets = copy.deepcopy(exons)
            complete_targets.extend(uncaptured_coordinates)
            subregions = mip.merge_overlap(complete_targets, 2 * flank)
        # get the target coordinates if the capture type is "targets"
        elif capture_type == "targets":
            # create a list of target coordinates
            target_coordinates = []
            for ct in capture_targets:
                for c in capture_targets[ct]:
                    for t in capture_targets[ct][c]:
                        target_coordinates.append(
                            [capture_targets[ct][c][t]["begin"],
                             capture_targets[ct][c][t]["end"]])
            subregions = mip.merge_overlap(target_coordinates,
                                           spacer=2*flank)

        return subregions

    def rinfo_parser(self, rinfo_file):
        """Parse .rinfo file for a Paralog object."""
        settings_dict = {}
        must_list = []
        field_names = {}
        with open(rinfo_file) as infile:
            for line in infile:
                if not line.startswith("#"):
                    newline = line.strip().split("\t")
                    if len(newline) == 2:
                        settings_dict[newline[0]] = newline[1]
                    elif newline[0] == "MUST":
                        must_list.append(newline[1:])
                    else:
                        key = newline[0]
                        field = newline[1]
                        values = newline[2:]
                        if field == "settings_for":
                            settings_dict[key] = {}
                            field_names[key] = {}
                            for i in range(len(values)):
                                field_names[key][i] = values[i]
                        else:
                            for i in range(len(values)):
                                fname = field_names[key][i]
                                try:
                                    settings_dict[key][fname][field] = values[
                                        i]
                                except KeyError:
                                    settings_dict[key][fname] = {
                                        field: values[i]}
        must_dict = {}
        must_list = np.transpose(must_list).tolist()
        try:
            f_names = must_list[0]
            for i in range(1, len(must_list)):
                m_dict = dict(list(zip(f_names, must_list[i])))
                must_dict[m_dict["name"]] = m_dict
        except IndexError:
            pass
        settings_dict["MUST"] = must_dict
        return settings_dict

    def get_gene_names(self, rinfo):
        """Currently Unused."""
        gene_names = []
        region_dic = rinfo["REGION"]
        for r in region_dic:
            gene_names.append([region_dic[r]["copyname"],
                              region_dic[r]["chrom"]])
        r = []
        for g in gene_names:
            if g not in r:
                r.append(g)
        return r

    def add_subregion(self, name, subregion):
        """Add a subregion object to subregions dictionary of a segment."""
        self.subregions[name] = Subregion(self, name, subregion)

    def make_primers(self):
        """Design primers for extension and ligation arms."""
        ext_list = []
        lig_list = []
        ext_set = self.rinfo["SETTINGS"]["extension"]["settings_file"]
        lig_set = self.rinfo["SETTINGS"]["ligation"]["settings_file"]
        processors = int(self.rinfo["SETTINGS"]["mip"]["processors"])
        for s in self.subregions:
            dir_in = self.subregions[s].primer3_input_DIR
            extension_boulder = self.subregions[s].fullname + "_ext"
            ligation_boulder = self.subregions[s].fullname + "_lig"
            dir_out = self.subregions[s].primer3_output_DIR
            extension_out = self.subregions[s].fullname + "_ext"
            ligation_out = self.subregions[s].fullname + "_lig"
            primer_settings_dir = self.paralog.resource_dir
            ext_list.append([extension_boulder, ext_set, extension_out, dir_in,
                             dir_out, primer_settings_dir,
                             s, self.paralog_name, "extension"])
            lig_list.append([ligation_boulder, lig_set, ligation_out, dir_in,
                             dir_out, primer_settings_dir,
                             s, self.paralog_name, "ligation"])
            self.subregions[s].primers = {"original": {
                "extension": extension_out, "ligation": ligation_out}}
        mip.make_primers_multi(ext_list, lig_list, processors)
        return

    def set_segment_rinfo(self):
        # this function is only to be called from
        # set_rinfo function of the paralog class
        self.rinfo = self.paralog.rinfo
        for s in self.subregions:
            self.subregions[s].update_subregion()


class Paralog(Locus):
    """
    Represent a gene or a paralogus gene family.

    Requires a region information (rinfo) file which is generated by the
    upstream region prep module. Instantiates a Locus object using the sequence
    information in the rinfo file.
    """

    def __init__(self, rinfo_file):
        # rifo file has capture specifics for the mips
        self.rinfo_file = rinfo_file
        # convert rinfo file to dictionary
        self.rinfo = self.rinfo_parser(self.rinfo_file)
        # get species and host species from rinfo dict
        self.species = self.rinfo["SPECIES"]
        # host species can be none, or any contaminating species
        # it is created for human, when designing parasite mips
        self.host = self.rinfo["HOST_SPECIES"]
        # get the design directory relative to root
        self.design_dir = self.rinfo["DESIGN_DIR"]
        # check if chained overlapping MIPs are desired
        self.chain_mips = int(
            self.rinfo["SELECTION"]["compatibility"]["chain"])
        # get file location dictionary for locations of static files
        # such as genome fasta files, snp files etc.
        self.file_locations = mip.get_file_locations()
        # give the object a name
        self.paralog_name = self.rinfo["NAME"]
        # resource_dir is employed to indicate the location of rinfo and other
        # files. Files should be in a directory named as the paralog family
        # name which is set up by the upstream region prep module.
        self.cwd = os.path.join(self.design_dir, self.paralog_name) + "/"
        self.resource_dir = os.path.join(self.cwd, "resources") + "/"
        # load the target dictionary, generated by the upstream region prep
        # module. It holds the information regarding capture targets such as
        # genomic coordinates.
        target_dict_file = os.path.join(self.resource_dir,
                                        self.paralog_name + ".pkl")
        with open(target_dict_file, "rb") as infile:
            self.target_dict = pickle.load(infile)
        self.regions = self.target_dict["regions"]
        self.alignments = self.target_dict["alignments"]
        self.target_genes = self.target_dict["genes"]
        self.segment_dic = self.get_segment_dic()["segment_dic"]
        self.begin = self.segment_dic["S0"]["C0"]["begin"]
        self.end = self.segment_dic["S0"]["C0"]["end"]
        self.chrom = self.segment_dic["S0"]["C0"]["chrom"]
        self.extras_dic = self.get_segment_dic()["extras_dic"]
        self.copies = self.get_segment_dic()["segment_copies"]
        self.zcoordinates = self.segment_converter()
        self.pcoordinates = self.base_converter()
        self.pdiffs = self.get_pdiffs()["pdiffs"]
        self.extra_snps = self.get_pdiffs()["extra_snps"]
        self.snps = self.get_snps()
        self.must = self.get_must()
        self.get_exons()
        self.exons = self.get_projected_exons()["merged"]
        self.segments = {}
        self.add_segments()
        self.designed = False
        self.failed = False

    def add_segments(self):
        """Create a locus object for each segment in the segment_dic."""
        for s in self.segment_dic:
            self.segments[s] = Locus(self, self.rinfo_file)

    def set_rinfo(self):
        """
        Update rinfo for the paralog, segment and gene objects.

        This is used when settings are changed in the rinfo file.
        """
        self.rinfo = self.rinfo_parser(self.rinfo_file)
        for s in self.segments:
            self.segments[s].set_segment_rinfo()

    def check_chained(self):
        """Check if all mips for this paralog are chained."""
        seg = self.segments["S0"]
        self.chained_mips = True
        for subregion in seg.subregions:
            try:
                if not seg.subregions[subregion].chained_mips:
                    self.chained_mips = False
                    return
            except AttributeError:
                self.chained_mips = False
                return

    def check_copies(self):
        """Check if all MIPs are capturing all copies in the paralog."""
        for mipname in self.mips:
            m = self.mips[mipname]
            cc = m.captured_copies
            if set(cc) != set(self.copies):
                self.copies_captured = False
                break
        else:
            self.copies_captured = True

    def run_segments(self, seg):
        """Make subregions and primers for each segment/gene."""
        if seg.designed:
            return
        else:
            s = seg.make_subregions()
            if len(s) == 0:
                print(("Nothing interesting in segment %s of paralog %s"
                       % (seg.segment_name, seg.paralog.paralog_name)))
                return "failed"
            s.sort(key=itemgetter(0))
            counter = 0
            for r in s:
                seg.add_subregion("Sub" + str(counter), r)
                counter += 1
            seg.make_primers()
            seg.designed = True
            with open(self.cwd + self.paralog_name, "wb") as savefile:
                pickle.dump(self, savefile)
            return

    def run_subregions(self, sub):
        """Run mip design pipeline for the given subregion."""
        if sub.designed:
            return
        elif not sub.failed:
            try:
                sub.parse_primers()
                sub.bowtie2_run()
                sub.score_primers()
                picked_pairs = sub.pick_primer_pairs()
                if picked_pairs == 1:
                    print(sub.fullname, " has no primer pairs.")
                    return
                sub.make_mips()

                with open(self.cwd + self.paralog_name, "wb") as savefile:
                    pickle.dump(self, savefile)

                sub.hairpin()
                sub.score_mips()
                sub.compatible()
                sub.best_mipset()
                sub.designed = True
                with open(self.cwd + self.paralog_name, "wb") as savefile:
                    pickle.dump(self, savefile)
            except Exception as e:
                sub.failed = True
                print(sub.fullname, " failed due to error ", str(e))
                traceback.print_exc()
        else:
            print(sub.fullname, " design has failed before",
                  " and will not be repeated.")
            return

    def run_paralog(self):
        """
        Top function that runs other functions in the design pipeline.

        Runs run_segments for segments and genes, run_subregions for
        subregions.
        """
        if self.designed:
            return
        else:
            for s in self.segments:
                seg = self.segments[s]
                if not seg.designed:
                    run_result = self.run_segments(seg)
                    if run_result == "failed":
                        continue
                for subregion in seg.subregions:
                    sub = seg.subregions[subregion]
                    self.run_subregions(sub)
            self.check_chained()
            self.summarize_mips()
            self.check_copies()
            self.print_info()
            self.sort_mips()
            self.order_mips()
            self.designed = True
            with open(self.cwd + self.paralog_name, "wb") as savefile:
                pickle.dump(self, savefile)

    def summarize_mips(self):
        """Extract summary information from mips in the paralog."""
        # output for listing brief mip information such as name and sequence
        self.mipfile = self.cwd + self.paralog_name + ".mips"
        outfile = open(self.mipfile, "w")
        summary_file = self.cwd + self.paralog_name + "_summary.txt"
        summary_list = []
        # create a list to hold all lines to be output
        outfile_list = []
        # create mips attribute for the paralog, a dictionary to hold
        # mip objects
        self.mips = {}
        # create locus info dictionary that holds coordinate information about
        # everything in the locus. This will be used in drawing by Parasight.
        locus_info = {"segments": {}, "exons": [], "subregions": {},
                      "snps": {}, "primers": {}, "all_targets": {}}
        # populate the locus info dictionary
        paralog_segments = self.segments
        for s in paralog_segments:
            # make sure each segment is designed and there was no failure
            if paralog_segments[s].designed and not paralog_segments[s].failed:
                # add segment coordinates
                locus_info["segments"][paralog_segments[s].segment_name] = {
                    "chrom": paralog_segments[s].chrom,
                    "begin": paralog_segments[s].begin,
                    "end": paralog_segments[s].end,
                    "copies": ",".join(paralog_segments[s].copies)}
                # add exons in the segment
                locus_info["exons"].extend(list(paralog_segments[s].exons))
                # add targets in the segment (snps, must capture regions,
                # pdiffs etc.)
                locus_info["all_targets"].update(paralog_segments[s].targets)
                # add subregion information
                for sub in paralog_segments[s].subregions:
                    # make sure all subregions run through without error
                    if (
                        paralog_segments[s].subregions[sub].designed
                        and not paralog_segments[s].subregions[sub].failed
                    ):
                        # add each mip from this subregion to the paralog
                        self.mips.update(paralog_segments[s].subregions[
                            sub].mips["best_mipset"]["dictionary"]["mips"])
                    # add subregion coordinates
                    locus_info["subregions"][paralog_segments[s].subregions[
                        sub].fullname] = {
                        "chrom": paralog_segments[s].subregions[sub].chrom,
                        "begin": paralog_segments[s].subregions[sub].begin,
                        "end": paralog_segments[s].subregions[sub].end}
                    # all potential primers in the subregion
                    locus_info["primers"][paralog_segments[s].subregions[
                        sub].fullname] = (
                       {"extension": paralog_segments[s].subregions[
                           sub].primers["parsed"]["extension"]["dictionary"],
                        "ligation": paralog_segments[s].subregions[
                            sub].primers["parsed"]["ligation"]["dictionary"]}
                    )
                    # number of primers designed and filtered
                    sub_pr = paralog_segments[s].subregions[sub].primers
                    for h_arm in ["extension", "ligation"]:
                        p_designed = len(sub_pr["parsed"][h_arm]
                                         ["dictionary"]["primer_information"])
                        p_hit = sub_pr["bt_added"][h_arm]["primers_remaining"]
                        p_hit_host = sub_pr["bt_added"][h_arm][
                            "primers_remaining_host"]
                        p_tm = sub_pr["bowtied"][h_arm]["primers_remaining"]
                        p_tm_host = sub_pr["bowtied"][h_arm][
                            "primers_remaining_host"]
                        p_score_filt = len(sub_pr["filtered"][h_arm]
                                           ["dictionary"]
                                           ["primer_information"])
                        sum_line = ("Paralog {} subregion {} {} arm had {} "
                                    "primers designed. {} remained after "
                                    "primers with too many alignments, "
                                    " {} after primers with too many "
                                    " alignments to host genome (if any) "
                                    "{} after primers with high TM alignments"
                                    " {} after primers with high TM alignments"
                                    " to host genome (if any), {} after "
                                    "primers with lower scores compared to "
                                    "competing primers were removed.").format(
                                        self.paralog_name, sub, h_arm,
                                        p_designed, p_hit, p_hit_host, p_tm,
                                        p_tm_host, p_score_filt
                                    )
                        summary_list.append(sum_line)
                    try:
                        pair_num = len(sub_pr["pairs"]["dictionary"]
                                       ["pair_information"])
                        hp_num = len(paralog_segments[s].subregions[sub].mips[
                                     "hairpin"])
                    except KeyError:
                        pair_num = hp_num = 0
                    sum_line = ("{} primer pairs were found for Paralog {} "
                                "subregion {} corresponding to the same "
                                "number of MIPs. {} MIPs remained after "
                                "filtering for hairpins.").format(
                                    pair_num, self.paralog_name, sub, hp_num)
                    summary_list.append(sum_line)
        # add copies in the segment to locus info
        segment_dic = self.segment_dic["S0"]
        copy_dic = {}
        for c in segment_dic:
            try:
                copy_name = "_".join(segment_dic[c]["genes"]["names"])
            except KeyError:
                copy_name = "none"
            if c != "C0":
                copy_begin = min(self.pcoordinates[c].values())
                copy_end = max(self.pcoordinates[c].values())

                copy_dic[c] = {"chrom": self.chrom,
                               "begin": copy_begin,
                               "end": copy_end,
                               "name": copy_name}
            else:
                copy_dic[c] = {"chrom": self.chrom,
                               "begin": segment_dic[c]["begin"],
                               "end": segment_dic[c]["end"],
                               "name": copy_name}
        locus_info["copies"] = copy_dic
        # add must-capture targets
        locus_info["targets"] = {}
        for c in self.must:
            for t in self.must[c]:
                locus_info["targets"][self.must[c][t]["name"]] = {
                    "chrom": self.must[c][t]["chrom"],
                    "begin": self.must[c][t]["begin"],
                    "end": self.must[c][t]["end"]
                }
        # header for mip order output
        outline = ["#pair_name", "mip_name", "chrom", "mip_start",
                   "capture_start", "capture_end", "mip_end", "orientation",
                   "tech_score", "func_score", "mip_sequence",
                   "unique_captures", "must_captured", "captured_copies"]
        outfile_list.append("\t".join(outline))
        must_cap = copy.deepcopy(self.must)
        # order mip names in the paralog object according to their genomic
        # start location
        mip_names_ordered = sorted(self.mips, key=lambda mip_key: self.mips[
            mip_key].mip["C0"]["mip_start"])
        locus_info["mips"] = []
        # rename mips showing their place on the genome
        name_counter = 0
        # are there any targets
        targets_present = False
        # record total number of MIPs
        sum_line = ("Total number of MIPs selected for Paralog {} was {}. "
                    ).format(self.paralog_name, len(mip_names_ordered))
        summary_list.append(sum_line)
        for mip_name in mip_names_ordered:
            m = self.mips[mip_name]
            fullname_subregion = m.subregion.fullname
            fullname_mip = fullname_subregion + "_mip" + str(name_counter)
            m.fullname = fullname_mip
            # get coordinate information of the mip
            # reference copy information will be used
            m.mip_start = m.mip["C0"]["mip_start"]
            m.mip_end = m.mip["C0"]["mip_end"]
            m.capture_start = m.mip["C0"]["capture_start"]
            m.capture_end = m.mip["C0"]["capture_end"]
            m.orientation = m.mip["C0"]["orientation"]
            m.chromosome = m.mip["C0"]["chrom"]
            caps = m.captures
            unique_caps = m.capture_info["unique_captures"]
            ms = []
            for c in self.must:
                for t in self.must[c]:
                    # set targets_present to True if there are any targets
                    targets_present = True
                    try:
                        if t in caps:
                            line = ("Target " + self.must[c][t]["name"]
                                    + " is captured.")
                            outfile_list.append("#" + line)
                            # add the name of captured target to list
                            ms.append(self.must[c][t]["name"])
                            sum_line = ("Target {} for paralog {} was "
                                        "captured.").format(
                                            self.must[c][t]["name"],
                                            self.paralog_name)
                            summary_list.append(sum_line)
                            # remove captured target from target dict
                            must_cap[c].pop(t)
                    except KeyError:
                        continue
            # if none of the must-targets captured by this mip, add none to ms
            if len(ms) == 0:
                ms = ["none"]
            # create extensions for mips to account for alternative mips for
            # different copies. A mip that captures reference copy will have
            # _ref extension, alternatives will have _alt#
            for key in m.mip_dic["mip_information"]:
                if key == "ref":
                    fullname_extension = "_ref"
                elif key.startswith("alt"):
                    fullname_extension = "_" + key
                else:
                    continue
                outlist = [m.name, m.fullname + fullname_extension,
                           m.chromosome, m.mip_start, m.capture_start,
                           m.capture_end, m.mip_end, m.orientation,
                           m.tech_score, m.func_score,
                           m.mip_dic["mip_information"][key]["SEQUENCE"],
                           ",".join(unique_caps), ",".join(ms),
                           ",".join(m.captured_copies)]
                locus_info["mips"].append(outlist)
                outline = "\t".join(map(str, outlist))
                outfile_list.append(outline)
            name_counter += 1
        # look for must-targets that have not been captured, print if any
        target_remaining = False
        for c in must_cap:
            if len(must_cap[c]) > 0:
                target_remaining = True
                for t in must_cap[c]:
                    line = ("Target" + " " + self.must[c][t]["name"]
                            + " is NOT captured.")
                    outfile_list.append("#" + line)
                    sum_line = ("Target {} for paralog {} was NOT "
                                "captured.").format(self.must[c][t]["name"],
                                                    self.paralog_name)
                    summary_list.append(sum_line)
        if targets_present and not target_remaining:
            line = "All targets have been captured."
            self.must_captured = True
            outfile_list.append("#" + line + "\n")
            sum_line = ("All targets for paralog {} were captured.").format(
                self.paralog_name)
            summary_list.append(sum_line)
        else:
            self.must_captured = False
        if self.chain_mips:
            if self.chained_mips:
                line = "All mips are chained"
            else:
                line = "Mips are not chained"
            outfile_list.append("#" + line + "\n")
            sum_line = (line + " for paralog {}.").format(
                self.paralog_name)
            summary_list.append(sum_line)
        outfile.write("\n".join(outfile_list))
        outfile.close()
        with open(summary_file, "w") as outfile:
            outfile.write("\n".join(summary_list) + "\n")
        print("\n".join(summary_list) + "\n")
        self.locus_info = locus_info
        return

    def print_info(self):
        """Create visualization files for parasight."""
        # set drawing parameters for each label: color, vertical placement,
        # thickness
        info = {"exons": {"color": "steel blue", "offset": -15, "width": 4},
                "copies": {"color": ["dark slate blue", "slate blue"],
                           "offset": [-10, 1.5], "width": 2},
                "extension": {"color": ["cyan4"], "offset": [-45, -10, -30],
                              "width": 4},
                "ligation": {"color": ["maroon"], "offset": [-45, -10, -30],
                             "width": 4},
                "capture": {"color": ["black", "gray"],
                            "offset": [-45, -10, -30], "width": 2},
                "subregions": {"color": "rosy brown", "offset": -21,
                               "width": 4},
                "targets": {"color": "gold", "offset": -26, "width": 8},
                "snps": {"color": "red", "offset": -26, "width": 4},
                "pdiffs": {"color": "dark green", "offset": -26, "width": 4},
                "primers": {"color": {"forward_extension": "cyan4",
                                      "reverse_extension": "maroon",
                                      "forward_ligation": "maroon",
                                      "reverse_ligation": "cyan4"},
                            "offset": {"forward_extension": -22.5,
                                       "reverse_extension": -20.5,
                                       "forward_ligation": -22,
                                       "reverse_ligation": -20},
                            "width": {"forward_extension": 0.5,
                                      "reverse_extension": 0.5,
                                      "forward_ligation":  0.5,
                                      "reverse_ligation": 0.5}}}

        outlist = ["seq", "begin", "end", "label", "color", "offset", "width",
                   "name", "uniq_copies", "targets_capped", "copies_bound",
                   "tech_score"]

        begin_coordinates = []
        end_coordinates = []
        outfile_list = ["\t".join(map(str, outlist))]
        copies = self.locus_info["copies"]
        subregions = self.locus_info["subregions"]
        exons = self.locus_info["exons"]
        mips = self.locus_info["mips"]
        targets = self.locus_info["targets"]
        primers = self.locus_info["primers"]
        snps = self.segments["S0"].split_snps.to_dict(orient="index")
        pdiffs = self.segments["S0"].split_pdiffs.to_dict(orient="index")
        for c in sorted(copies):
            chrom = copies[c]["chrom"]
            copy_number = int(c[1:])
            width = info["copies"]["width"]
            begin_coordinates.append(copies[c]["begin"])
            end_coordinates.append(copies[c]["end"])
            outlist = [copies[c]["chrom"], copies[c]["begin"],
                       copies[c]["end"], "copy",
                       info["copies"]["color"][copy_number % 2],
                       (info["copies"]["offset"][0] +
                        width * copy_number * info["copies"]["offset"][1]),
                       width, copies[c]["name"], c]
            outfile_list.append("\t".join(map(str, outlist)))
        for s in subregions:
            outlist = [subregions[s]["chrom"], subregions[s]["begin"],
                       subregions[s]["end"], "subregion",
                       info["subregions"]["color"],
                       info["subregions"]["offset"],
                       info["subregions"]["width"], s]
            outfile_list.append("\t".join(map(str, outlist)))
        for s in exons:
            outlist = [chrom, s[0], s[1], "exon", info["exons"]["color"],
                       info["exons"]["offset"], info["exons"]["width"], s]
            outfile_list.append("\t".join(map(str, outlist)))
        for s in targets:
            outlist = [targets[s]["chrom"], targets[s]["begin"],
                       targets[s]["end"] + 1,  "target",
                       info["targets"]["color"], info["targets"]["offset"],
                       info["targets"]["width"], s]
            outfile_list.append("\t".join(map(str, outlist)))
        # information in mips is organized in a list
        # ["#pair_name", "mip_name", "chrom", "mip_start", "capture_start",
        # "capture_end", "mip_end", "orientation", "tech_score", "func_score",
        # "mip_sequence", "unique_captures", "must_captured",
        # "captured_copies"]
        for i in mips:
            coordinates = [[i[2], i[3], i[4] - 1],
                           [i[2], i[4], i[5]],
                           [i[2], i[5] + 1, i[6]]]
            # check if the mip covers the reference (ending with _ref)
            # or an alternative copy (ending with _alt#)
            try:
                mip_number = int(i[1].split("_")[-1][3:])
            except ValueError:
                mip_number = -1
            mip_number += 1
            jitter = random.choice(list(range(10)))
            if i[7] == "forward":
                this_arm = "extension"
                next_arm = "capture"
                last_arm = "ligation"
                color_index = 0
                offset_index = 0
            else:
                this_arm = "ligation"
                next_arm = "capture"
                last_arm = "extension"
                color_index = -1
                offset_index = -1
            arms = [this_arm, next_arm, last_arm]
            for j in range(len(arms)):
                prt = []
                armtype = arms[j]
                prt.extend(coordinates[j])
                prt.append(armtype)
                prt.append(info[armtype]["color"][color_index])
                prt.append(info[armtype]["offset"][offset_index]
                           + info[armtype]["offset"][1] * mip_number
                           - jitter)
                prt.append(info[armtype]["width"])
                prt.append(i[1])
                if armtype == "capture":
                    prt.append(i[1])
                else:
                    prt.append("")
                try:
                    prt.append(i[11])
                    prt.append(i[13])
                except IndexError:
                    prt.extend(["", ""])
                if armtype == "capture":
                    prt.append(i[8])
                else:
                    prt.append("")

                outfile_list.append("\t".join(map(str, prt)))
        for s in snps:
            outlist = [snps[s]["chrom"], snps[s]["position"],
                       snps[s]["position"] + 1, "snp", info["snps"]["color"],
                       info["snps"]["offset"], info["snps"]["width"], s]
            outfile_list.append("\t".join(map(str, outlist)))
        for s in pdiffs:
            outlist = [pdiffs[s]["chrom"], pdiffs[s]["position"],
                       pdiffs[s]["position"] + 1, "pdiffs",
                       info["pdiffs"]["color"], info["pdiffs"]["offset"],
                       info["pdiffs"]["width"], s]
            outfile_list.append("\t".join(map(str, outlist)))
        for sbr in primers:
            ext_primers = primers[sbr]["extension"]["primer_information"]
            lig_primers = primers[sbr]["ligation"]["primer_information"]
            for p in ext_primers:
                primer_ori = ext_primers[p]["ORI"]
                if primer_ori == "forward":
                    outlist = [ext_primers[p]["CHR"],
                               ext_primers[p]["GENOMIC_START"],
                               ext_primers[p]["GENOMIC_END"], p,
                               info["primers"]["color"]["forward_extension"],
                               info["primers"]["offset"]["forward_extension"],
                               info["primers"]["width"]["forward_extension"]]
                    outfile_list.append("\t".join(map(str, outlist)))
                else:
                    outlist = [ext_primers[p]["CHR"],
                               ext_primers[p]["GENOMIC_END"],
                               ext_primers[p]["GENOMIC_START"], p,
                               info["primers"]["color"]["reverse_extension"],
                               info["primers"]["offset"]["reverse_extension"],
                               info["primers"]["width"]["reverse_extension"]]
                    outfile_list.append("\t".join(map(str, outlist)))
            for p in lig_primers:
                primer_ori = lig_primers[p]["ORI"]
                if primer_ori == "forward":
                    outlist = [lig_primers[p]["CHR"],
                               lig_primers[p]["GENOMIC_START"],
                               lig_primers[p]["GENOMIC_END"], p,
                               info["primers"]["color"]["reverse_ligation"],
                               info["primers"]["offset"]["reverse_ligation"],
                               info["primers"]["width"]["reverse_ligation"]]
                    outfile_list.append("\t".join(map(str, outlist)))
                else:
                    outlist = [lig_primers[p]["CHR"],
                               lig_primers[p]["GENOMIC_END"],
                               lig_primers[p]["GENOMIC_START"], p,
                               info["primers"]["color"]["forward_ligation"],
                               info["primers"]["offset"]["forward_ligation"],
                               info["primers"]["width"]["forward_ligation"]]
                    outfile_list.append("\t".join(map(str, outlist)))
        self.showfile = self.cwd + self.paralog_name + ".show"
        self.extrafile = self.cwd + self.paralog_name + ".extra"
        with open(self.showfile, "w") as show, open(
                self.extrafile, "w") as extra:
            show.write("header\n")
            show.write("\t".join(map(str, [chrom,  max(end_coordinates),
                                           min(begin_coordinates),
                                           max(end_coordinates)])))
            extra.write("\n".join(outfile_list))
        # Below attributes will be assigned by sort_mips function
        # when some mips are selected to be included or excluded from
        # the design. If no such selection is used, we assign them here
        # to the "unselected" values.
        self.selected_mips = self.mips
        self.selectedfile = self.extrafile
        self.selected_mipfile = self.mipfile
        return

    def sort_mips(self, filter_type="remove"):
        purple = "#a0462009f0e4"  # purple in parasight
        blue = "#00830000ffff"
        colors = [blue, purple, "blue", "purple"]
        self.selectedfile = self.cwd + self.paralog_name + "_selected.pse"
        parasight_file = self.selectedfile
        try:
            with open(parasight_file) as infile:
                para_list = infile.readlines()
        except IOError:
            return
        if filter_type not in ["keep", "remove"]:
            print(("Filter type must be either 'keep' or 'remove'."
                   " No filtering will be applied for type {}")).format(
                       filter_type)
        if filter_type == "keep":
            select = set()
        else:
            select = set([self.mips[m].fullname for m in self.mips])
        for i in para_list:
            fields = i.strip().split("\t")
            label_name = fields[7]
            # check if line describes a mip
            if label_name in ["extension", "ligation", "capture"]:
                mip_name = "_".join(fields[8].split("_")[:-1])
                if filter_type == "keep":
                    if fields[3] in colors:
                        select.add(mip_name)
                elif filter_type == "remove":
                    if fields[3] in colors:
                        select.remove(mip_name)
        outlist = []
        for i in para_list:
            fields = i.strip().split("\t")
            label_name = fields[7]
            # check if line describes a mip
            if label_name in ["extension", "ligation", "capture"]:
                mip_name = "_".join(fields[8].split("_")[:-1])
                if mip_name in select:
                    outlist.append(i)
            else:
                outlist.append(i)
        self.selected_mips = {}
        for m in self.mips:
            if self.mips[m].fullname in select:
                self.selected_mips[m] = self.mips[m]
        self.selectedfile = self.cwd + self.paralog_name + "_filtered.pse"
        outfile = open(self.selectedfile, "w")
        outfile.write("".join(outlist))
        infile.close()
        outfile.close()
        outfile_list = []
        must_cap = copy.deepcopy(self.must)
        mip_names_ordered = sorted(self.selected_mips,
                                   key=lambda mip_key: self.selected_mips[
                                       mip_key].mip["C0"]["mip_start"])
        for mip_name in mip_names_ordered:
            m = self.selected_mips[mip_name]
            m.mip_start = m.mip["C0"]["mip_start"]
            m.mip_end = m.mip["C0"]["mip_end"]
            m.capture_start = m.mip["C0"]["capture_start"]
            m.capture_end = m.mip["C0"]["capture_end"]
            m.orientation = m.mip["C0"]["orientation"]
            m.chromosome = m.mip["C0"]["chrom"]
            caps = m.captures
            unique_caps = m.capture_info["unique_captures"]
            ms = []
            for c in self.must:
                for t in self.must[c]:
                    try:
                        if t in caps:
                            line = ("Target " + self.must[c][t]["name"]
                                    + " is captured.")
                            print(("Target %s of %s is captured."
                                   % (self.must[c][t]["name"],
                                      self.paralog_name)))
                            outfile_list.append("#" + line)
                            # add the name of captured target to list
                            ms.append(self.must[c][t]["name"])
                            # remove captured target from target dict
                            must_cap[c].pop(t)
                    except KeyError:
                        continue
            # if none of the must-targets captured by this mip,
            # add none to ms list
            if len(ms) == 0:
                ms = ["none"]
            # create extensions for mips to account for alternative mips for
            # different copies. A mip that captures reference copy will have
            # _ref extension, alternatives will have _alt#
            for key in m.mip_dic["mip_information"]:
                if key == "ref":
                    fullname_extension = "_ref"
                elif key.startswith("alt"):
                    fullname_extension = "_" + key
                else:
                    continue
                # ["#pair_name", "mip_name", "chrom", "mip_start",
                # "capture_start", "capture_end", "mip_end", "orientation",
                # "tech_score", "func_score","mip_sequence",
                # "unique_captures", "must_captured", "captured_copies"]
                outlist = [m.name, m.fullname + fullname_extension,
                           m.chromosome, m.mip_start, m.capture_start,
                           m.capture_end, m.mip_end, m.orientation,
                           m.tech_score, m.func_score,
                           m.mip_dic["mip_information"][key]["SEQUENCE"],
                           ",".join(unique_caps), ",".join(ms),
                           ",".join(m.captured_copies)]
                outline = "\t".join(map(str, outlist))
                outfile_list.append(outline)
        # look for must-targets that have not been captured, print if any
        target_remaining = False
        for c in must_cap:
            if len(must_cap[c]) > 0:
                target_remaining = True
                for t in must_cap[c]:
                    line = ("Target" + " " + self.must[c][t]["name"]
                            + " is NOT captured.")
                    outfile_list.append("#" + line)
                    print(("Target %s of %s is NOT captured."
                           % (self.must[c][t]["name"], self.paralog_name)))
        if not target_remaining:
            line = "All targets have been captured."
            outfile_list.append("#" + line + "\n")
            print("All targets have been captured for %s" % self.paralog_name)
        else:
            print("Some targets have NOT been captured for %s"
                  % self.paralog_name)
        self.selected_mipfile = self.cwd + self.paralog_name + "_order"
        outfile = open(self.selected_mipfile, "w")
        outfile.write("\n".join(outfile_list))
        outfile.close()
        with open(self.cwd + self.paralog_name, "wb") as savefile:
            pickle.dump(self, savefile)

    def order_mips(self, filter_type="remove"):
        """
        Create call_info and mip_info files for mip order.

        These files are used by the order_mips.py script to prepare the
        actual probe order. The files are also used in the data analysis
        pipeline.
        """
        mip_info = {"mips": {}, "paralog_info": {},
                    "mip_names": {"mip_to_pair": {}, "pair_to_mip": {}}}
        # use sort_mips function to remove any filtered MIPs, sort MIPs
        # according to their genomic position and rename according to that
        # position.
        self.sort_mips(filter_type)
        desired_attr = ['target_genes', 'exons', 'pcoordinates', 'extra_snps',
                        'segment_dic', 'extras_dic', 'copies', 'regions',
                        'failed', 'begin', 'copies_captured', 'snps',
                        'target_dict', 'chained_mips', 'zcoordinates',
                        'designed', 'end', 'chrom', 'must', 'pdiffs',
                        'rinfo', 'paralog_name', 'alignments', "species"]
        for d in desired_attr:
            mip_info["paralog_info"][d] = self.__dict__[d]
        for pair_name in self.selected_mips:
            m = self.selected_mips[pair_name]
            mip_name = m.fullname
            mip_info["mip_names"]["mip_to_pair"][mip_name] = pair_name
            mip_info["mip_names"]["pair_to_mip"][pair_name] = mip_name
            mip_info["mips"][mip_name] = {"mip_info": m.mip,
                                          "capture_info": m.capture_info,
                                          "mip_dic": m.mip_dic}
            ligation = mip_info["mips"][mip_name]["mip_dic"][
                "ligation_primer_information"]['PARALOG_COORDINATES']
            extension = mip_info["mips"][mip_name]["mip_dic"][
                "extension_primer_information"]['PARALOG_COORDINATES']
            mipdic = mip_info["mips"][mip_name]["mip_info"]
            for c in ligation:
                if c not in mipdic:
                    try:
                        ext_start = extension[c]['ALT_START']
                        ext_end = extension[c]['ALT_END']
                    except KeyError:
                        try:
                            ext_start = extension[c]['BOWTIE_START']
                            ext_end = extension[c]['BOWTIE_END']
                        except KeyError:
                            try:
                                ext_start = extension[c]['GENOMIC_START']
                                ext_end = extension[c]['GENOMIC_END']
                            except KeyError:
                                continue
                    try:
                        lig_start = ligation[c]['ALT_START']
                        lig_end = ligation[c]['ALT_END']
                    except KeyError:
                        try:
                            lig_start = ligation[c]['BOWTIE_START']
                            lig_end = ligation[c]['BOWTIE_END']
                        except KeyError:
                            try:
                                lig_start = ligation[c]['GENOMIC_START']
                                lig_end = ligation[c]['GENOMIC_END']
                            except KeyError:
                                continue
                    chrom = ligation[c]["CHR"]
                    orientation = extension[c]["ORI"]
                    mip_start = sorted([ext_start, ext_end,
                                        lig_start, lig_end])[0]
                    mip_end = sorted([ext_start, ext_end,
                                      lig_start, lig_end])[3]
                    capture_start = sorted([ext_start, ext_end,
                                            lig_start, lig_end])[1] + 1
                    capture_end = sorted([ext_start, ext_end,
                                          lig_start, lig_end])[2] - 1
                    capture_size = capture_end - capture_start + 1
                    capture_key = chrom + ":" + str(capture_start) + "-" + str(
                        capture_end)
                    mipdic[c] = {"chrom": chrom,
                                 "capture_start": capture_start,
                                 "capture_end": capture_end,
                                 "mip_start": mip_start,
                                 "mip_end": mip_end,
                                 "capture_size": capture_size,
                                 "capture_key": capture_key,
                                 "captured": False,
                                 "orientation": orientation,
                                 "ligation_start": lig_start,
                                 "ligation_end": lig_end,
                                 "extension_start": ext_start,
                                 "extension_end": ext_end}
                else:
                    mipdic[c]["captured"] = True
        call_info = {}
        for mip_name in mip_info["mips"]:
            minfo = mip_info["mips"][mip_name]["mip_info"]
            call_info[mip_name] = {"copies": {}}
            copies = call_info[mip_name]["copies"]
            for c in minfo:
                chrom = minfo[c]["chrom"]
                capture_start = minfo[c]["capture_start"]
                capture_end = minfo[c]["capture_end"]
                ori = minfo[c]["orientation"]
                captured = minfo[c]["captured"]
                try:
                    genes = mip_info["paralog_info"]["segment_dic"][
                        "S0"][c]["genes"]["names"]
                except KeyError:
                    genes = []
                try:
                    cap_seq = minfo[c]["capture_sequence"]
                except KeyError:
                    cap_seq = mip.get_sequence(minfo[c]["capture_key"],
                                               self.species)
                    if ori == "reverse":
                        cap_seq = mip.reverse_complement(cap_seq)
                minfo[c]["capture_sequence"] = cap_seq
                copies[c] = {"capture_sequence": cap_seq,
                             "chrom": chrom,
                             "capture_start": capture_start,
                             "capture_end": capture_end,
                             "orientation": ori,
                             "captured": captured,
                             "copyname": self.paralog_name + "-" + c,
                             "genes": genes}
        for m in call_info:
            copies = call_info[m]["copies"]
            md = mip_info["mips"][m]
            caps = md["capture_info"]["captured_targets"]
            for target_type in caps:
                for copy_id in caps[target_type]:
                    try:
                        copies[copy_id]["variants"][target_type] = caps[
                            target_type][copy_id]
                    except KeyError:
                        copies[copy_id]["variants"] = {target_type: caps[
                            target_type][copy_id]}
        return {"call_info": call_info, "mip_info": mip_info}


class Subregion(Locus):
    """
    Subregion class is used to create objects representing DNA loci.

    This DNA locus is a subregion of interest of a paralog object. For example,
    if exonic regions of a gene is to be targeted for MIP capture,
    each exon would be a subregion. MIPs are designed at the subregion level.
    """

    def __init__(self, locus, name, intervals):
        self.subregion_name = name
        self.locus = locus
        self.intervals = intervals
        self.designed = False
        self.failed = False
        self.update_subregion()

    def update_subregion(self):
        self.chrom = self.locus.chrom
        self.snps = self.locus.snps
        self.capture = self.locus.rinfo["CAPTURE"]["S0"]
        if int(self.capture["single_mip"]):
            self.single = True
        else:
            self.single = False
        self.settings = self.locus.rinfo["SETTINGS"]
        self.flank = int(self.capture["flank"])
        self.segment_begin = self.locus.begin
        self.segment_end = self.locus.end
        if (self.intervals[0] - self.segment_begin) >= self.flank:
            self.begin = self.intervals[0] - self.flank
        else:
            self.begin = self.segment_begin
        if (self.segment_end - self.intervals[1]) >= self.flank:
            self.end = self.intervals[1] + self.flank
        else:
            self.end = self.segment_end
        self.region_key = (self.chrom + ":" + str(self.begin + 1)
                           + "-" + str(self.end))
        self.fullname = (self.locus.paralog.paralog_name
                         + "_" + self.locus.segment_name
                         + "_" + self.subregion_name)

        self.cwd = self.locus.cwd + self.subregion_name + "/"
        self.species = self.locus.species
        self.host = self.locus.host
        self.create_dirs()
        self.targets = self.get_targets()
        self.masking()
        self.update_scoring()

    def get_targets(self):
        """Get targets of parent Paralog object falling into the subregion."""
        targets = {}
        for target_type in self.locus.targets:
            for copy_id in self.locus.targets[target_type]:
                for t in self.locus.targets[target_type][copy_id]:
                    if (self.begin
                            <= self.locus.targets[target_type]
                            [copy_id][t]["begin"] <= self.end):
                        try:
                            targets[target_type][copy_id][t] = (
                                self.locus.targets[target_type][copy_id][t]
                            )
                        except KeyError:
                            try:
                                targets[target_type][copy_id] = {
                                    t: self.locus.targets
                                    [target_type][copy_id][t]
                                }
                            except KeyError:
                                targets[target_type] = {
                                    copy_id: {
                                        t: self.locus.targets
                                        [target_type][copy_id][t]
                                    }
                                }
        return targets

    def create_dirs(self):
        """Create standard subdirectories for MIP design."""
        if not os.path.exists(self.cwd):
            # create a list of subdirectories, for input and output files
            # to be used throughout the program for this subregion.
            env = mip.create_dirs(self.cwd)
            # create the main dir
            os.makedirs(self.cwd)
            # create subdirs
            for path in env:
                os.makedirs(path)
        # assign the directory variables.
        self.primer3_input_DIR = self.cwd + "primer3_input_files/"
        self.primer3_output_DIR = self.cwd + "primer3_output_files/"
        self.bowtie2_input_DIR = self.cwd + "bowtie2_input/"
        self.bowtie2_output_DIR = self.cwd + "bowtie2_output/"
        self.mfold_input_DIR = self.cwd + "mfold_input/"
        self.mfold_output_DIR = self.cwd + "mfold_output/"
        return

    def masking(self):
        """
        Prepare primer design templates.

        Prepare the design templates in boulder record format to be used
        by primer3.
        """
        # get masking criteria from settings
        mask_diffs_lig = int(self.capture["mask_diffs_lig"])
        mask_diffs_ext = int(self.capture["mask_diffs_ext"])
        mask_snps_lig = int(self.capture["mask_snps_lig"])
        mask_snps_ext = int(self.capture["mask_snps_ext"])
        # get region sequence
        region_fasta = mip.get_fasta(self.region_key, self.species)
        # convert fasta record to one line string without the fasta identifier
        fasta_list = region_fasta.split("\n")
        fasta_head = fasta_list[0]
        seq_template_temp = "".join(fasta_list[1:]).upper()
        # convert fasta string to list of characters
        seq_template_list_lig = list(seq_template_temp)
        seq_template_list_ext = list(seq_template_temp)
        # keep positions to be excluded from primer design in a list
        exclude_lig = []
        exclude_ext = []
        # get locus diffs (paralogus differences, target SNPs etc.)
        region_diffs = copy.deepcopy(self.locus.split_pdiffs)
        # diff positions are genomic but we need local positions
        # get the local (index) of diffs by subtracting the genomic region
        # start position from diff position.
        region_diffs["index_position"] = region_diffs["position"] - self.begin
        # limit diffs to those within this subregion
        region_diffs = region_diffs.loc[
            (region_diffs["position"] >= self.begin)
            & (region_diffs["position"] <= self.end)
        ]
        diff_positions = region_diffs["index_position"].tolist()
        # two options for excluding positions from designs are mask and
        # exclude. Masking is simple lowercase masking which can be used
        # by primer3 to avoid placing 3' end of a primer by setting
        # lowercase mask option. Additionally, lowercase masking is a feature
        # that is propagated down the pipeline which can be used in scoring
        # each primer according to lowercase positions as opposed to letting
        # primer3 make a final decision to design the primer or not.
        # Exclude option prevents any primer overlapping with a given position
        # and thus is more strict.
        if mask_diffs_lig:
            for i in diff_positions:
                seq_template_list_lig[i] = seq_template_list_lig[i].lower()
        else:
            for i in diff_positions:
                exclude_lig.append([i, i])
        if mask_diffs_ext:
            for i in diff_positions:
                seq_template_list_ext[i] = seq_template_list_ext[i].lower()
        else:
            for i in diff_positions:
                exclude_ext.append([i, i])
        # carry out masking procedure for region SNPs
        region_snps = copy.deepcopy(self.locus.split_snps)
        region_snps["index_position"] = region_snps["position"] - self.begin
        region_snps = region_snps.loc[(region_snps["position"] >= self.begin)
                                      & (region_snps["position"] <= self.end)]
        snp_positions = region_snps["index_position"].tolist()
        if mask_snps_lig:
            for i in snp_positions:
                seq_template_list_lig[i] = seq_template_list_lig[i].lower()
        else:
            for i in snp_positions:
                exclude_lig.append([i, i])
        if mask_snps_ext:
            for i in snp_positions:
                seq_template_list_ext[i] = seq_template_list_ext[i].lower()
        else:
            for i in snp_positions:
                exclude_ext.append([i, i])
        # There is a limit to how many positions can be specified in an exclude
        # list for primer 3. We reduce the number we provide by collapsing
        # positions that are closer than 17 nucleotides assuming the minimum
        # primer length is 17  and there cannot be a primer between two
        # excluded positions that are closer than 17  nt.
        exclude_ext = mip.merge_overlap(exclude_ext, spacer=17)
        exclude_lig = mip.merge_overlap(exclude_lig, spacer=17)
        # boulder records specifies the positions as a space separated list
        # of comma separated values as index1,length1 index2,length2 ...
        # so we convert the index values (index_start, index_end) to that
        # format.
        exclude_ext = [[e[0], e[1] - e[0] + 1] for e in exclude_ext]
        exclude_lig = [[l[0], l[1] - l[0] + 1] for l in exclude_lig]
        # rebuild the fasta record from modified list of nucleotides.
        fasta_seq_lig = "".join(seq_template_list_lig)
        fasta_rec_lig = (">" + self.fullname + "_lig" + "," + fasta_head[1:]
                         + "\n" + fasta_seq_lig)
        fasta_seq_ext = "".join(seq_template_list_ext)
        fasta_rec_ext = (">" + self.fullname + "_ext" + "," + fasta_head[1:]
                         + "\n" + fasta_seq_ext)
        lig_output_name = self.fullname + "_lig"
        ext_output_name = self.fullname + "_ext"
        # create boulder records using the masked sequences.
        mip.make_boulder(fasta_rec_lig, self.primer3_input_DIR,
                         exclude_list=exclude_lig,
                         output_file_name=lig_output_name)
        mip.make_boulder(fasta_rec_ext, self.primer3_input_DIR,
                         exclude_list=exclude_ext,
                         output_file_name=ext_output_name)
        return

    def update_scoring(self):
        """
        Get scoring parameters from the settings (rinfo) file.

        It is possible to assign scores to each SNP/diff captured, although
        this is not used extensively. More documentation should be available
        for each scoring parameter where they are used.
        """
        self.scoring = {"snp_scores": {}, "diff_scores": {}}
        snps = self.capture["target_snp_functions"].split(",")
        if snps != ["none"]:
            snp_scores = list(map(int,
                              self.capture["score_snp_functions"].split(",")))
            for i in range(len(snps)):
                self.scoring["snp_scores"][snps[i]] = snp_scores[i]
        diffs = self.capture["target_diffs"].split(",")
        if diffs != ["none"]:
            diff_scores = list(map(int,
                               self.capture["score_target_diffs"].split(",")))
            for i in range(len(diffs)):
                self.scoring["diff_scores"][diffs[i]] = diff_scores[i]
        self.scoring["mask_penalty"] = int(self.capture["mask_penalty"])
        self.scoring["unique_copy_bonus"] = int(
            self.capture["unique_copy_bonus"]
        )
        self.scoring["alternative_copy_penalty"] = int(
            self.capture["alternative_copy_penalty"]
        )
        self.scoring["technical_score_coefficient"] = float(
            self.capture["technical_score_coefficient"]
        )
        self.scoring["chain_bonus"] = int(self.capture["chain_bonus"])
        self.scoring["chain_coverage"] = float(self.capture["chain_coverage"])
        self.scoring["must_bonus"] = int(self.capture["must_bonus"])
        self.scoring["set_copy_bonus"] = int(self.capture["set_copy_bonus"])
        return

    def parse_primers(self):
        """Parse primer file generated by primer3. Return a dict."""
        self.primers = {"original": {}, "parsed": {}}
        for arm in ["extension", "ligation"]:
            name = self.fullname + "_" + arm[:3]
            self.primers["original"][arm] = {"filename": name}
            self.primers["parsed"][arm] = {"filename": name + "_parsed"}
            parsed = mip.primer_parser3(
                self.primers["original"][arm]["filename"],
                self.primer3_output_DIR, self.bowtie2_input_DIR,
                self.primers["parsed"][arm]["filename"], fasta=1, outp=0
            )
            par = mip.paralog_primers(
                parsed, self.locus.copies, self.locus.pcoordinates,
                self.settings[arm], self.primer3_output_DIR,
                self.primers["parsed"][arm]["filename"],
                self.locus.species,
                outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"]))
            self.primers["parsed"][arm]["dictionary"] = par
        return

    def bowtie2_run(self):
        """
        Align primers to reference (and host) genomes.

        Use the alignments to decide which primers are on target, which have
        too many off targets and remove those with too many off targets.
        When there are multiple paralogs, the primers are designed on one
        paralog. This function also determines if a given primer is likely to
        bind to the other paralogs. If not, and if it is allowed, it will
        design alternative primers for those paralogs that are not likely to
        be bound by the original primer. The alternative primer is based on
        the original primer when possible, that is, when the paralog copy
        shows up as in the bowtie alignment, in which case, the
        bowtie alinged region of the paralog is extended or clipped from the
        5' end to match the melting temperature of the original primers. When
        bowtie alignment does not catch the paralog copy, the original
        lastz alignment used in region preps are used to get the alternative
        primer. Similarly, the sequence from the alignment will be adjusted
        from the 5' end to match the melting temperature.
        """
        species = self.species
        host = self.host
        self.bowtie = {}
        self.bowtie["input"] = {}
        self.bowtie["output"] = {}
        self.bowtie["cleaned"] = {}
        self.primers["bt_added"] = {}
        self.primers["bowtied"] = {}
        self.primers["alternative"] = {}
        self.primers["bt_filtered"] = {}
        output_level = int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
        # carry out everything for extension and ligation arms
        for arm in ["extension", "ligation"]:
            name = self.fullname + "_" + arm[:3]
            # parsed primer3 output will be used as bowtie input.
            # primer_parser creates a fasta file of primers in
            # the bowtie2_input_DIR.
            self.bowtie["input"][arm] = {
                "filename": self.primers["parsed"][arm]["filename"]
            }
            # output files will be generated for both the target and the
            # host species in the bowtie2_output_DIR.
            self.bowtie["output"][arm] = {"filename": name + "_" + species,
                                          "hostfile": name + "_" + host}
            # run bowtie2 for the species
            mip.bowtie2_run(self.bowtie["input"][arm]["filename"],
                            self.bowtie["output"][arm]["filename"],
                            self.bowtie2_input_DIR, self.bowtie2_output_DIR,
                            species,
                            process_num=int(self.settings[arm]["processors"]),
                            mode=self.settings[arm]["bowtie_mode"],
                            seed_len=int(self.settings[arm]["seed_len"]),
                            local=int(self.settings[arm]["local"]))
            # run bowtie2 for host species
            if host != "none":
                mip.bowtie2_run(
                    self.bowtie["input"][arm]["filename"],
                    self.bowtie["output"][arm]["hostfile"],
                    self.bowtie2_input_DIR, self.bowtie2_output_DIR, host,
                    process_num=int(self.settings[arm]["processors"]),
                    mode=self.settings[arm]["bowtie_mode"],
                    seed_len=int(self.settings[arm]["seed_len"]),
                    local=int(self.settings[arm]["local"]))
            # parse bowtie hits, remove primers with excessive hits
            # and add bowtie sequence information to primer dict
            # first for the species, create the filename for primers dictionary
            # that has the bowtie information for the species
            self.primers["bt_added"][arm] = {
                "filename": name + "_" + species + "_bt_added",
                "hostfile": name + "_" + host + "_bt_added"
            }
            bt_added = mip.parse_bowtie(
                self.primers["parsed"][arm]["dictionary"],
                self.bowtie["output"][arm]["filename"],
                self.primers["bt_added"][arm]["filename"],
                self.primer3_output_DIR,
                self.bowtie2_output_DIR,
                self.species, self.settings[arm], outp=output_level
            )
            # add number of remaining primers after primers with excessive
            # alignments were removed
            self.primers["bt_added"][arm]["primers_remaining"] = len(
                bt_added["primer_information"])
            self.primers["bt_added"][arm]["primers_remaining_host"] = len(
                bt_added["primer_information"])
            # add bowtie information for host species primer
            # dict to be used here is species bt_added dict, not the original
            # and the bowtie file to be used is that of the host species
            # primers with excessive alignments (> upper_hit_limit) are
            # also removed at this step.
            if host != "none":
                bt_added_host = mip.parse_bowtie(
                    bt_added,
                    self.bowtie["output"][arm]["hostfile"],
                    self.primers["bt_added"][arm]["hostfile"],
                    self.primer3_output_DIR, self.bowtie2_output_DIR,
                    self.host, self.settings[arm], outp=output_level
                )
                # add number of remaining primers after primers with excessive
                # alignments to host genome were removed
                self.primers["bt_added"][arm]["primers_remaining_host"] = len(
                    bt_added["primer_information"])
            # add temperature information of bowtie hits and create alternative
            # arms if necessary and allowed (for paralogus genes only).
            self.primers["bowtied"][arm] = {
                "filename": name + "_" + species + "_bowtied",
                "hostfile": name + "_" + host + "_bowtied"
            }
            if host != "none":
                # if there is a host species, then the most up to date primers
                # are in the hostfile, rather than the species file. That will
                # be the starting input dict.
                # First, add the species bowtie information
                bowtied_host = mip.process_bowtie(
                    bt_added_host, self.primers["bowtied"][arm]["filename"],
                    self.primer3_output_DIR, self.bowtie2_output_DIR,
                    self.species, self.settings[arm],
                    outp=output_level)
                # then add the host information, using species-bowtied
                # dictionary as input
                bowtied = mip.process_bowtie(
                    bowtied_host, self.primers["bowtied"][arm]["hostfile"],
                    self.primer3_output_DIR, self.bowtie2_output_DIR,
                    self.host, self.settings[arm], True, outp=output_level
                )
            else:
                # if there is no host species, add bt information for species
                # here the primer file is that of the species
                bowtied = mip.process_bowtie(
                    bt_added, self.primers["bowtied"][arm]["filename"],
                    self.primer3_output_DIR, self.bowtie2_output_DIR,
                    self.species, self.settings[arm], outp=output_level
                )
            # filter primers with nonspecific binding
            self.primers["bt_filtered"][arm] = {
                "filename": name + "_" + species + "_bt_filtered",
                "hostfile": name + "_" + host + "_bt_filtered"}

            if host != "none":
                # the primer dict with both host and species bt information
                # will be used first
                filtered_host = mip.filter_bowtie(
                    bowtied, self.primers["bt_filtered"][arm]["hostfile"],
                    self.primer3_output_DIR, self.host,
                    TM=float(self.settings[arm]["filter_tm"]),
                    hit_threshold=int(self.settings[arm]["hit_threshold"]),
                    lower_tm=float(self.settings[arm]["lower_filter_tm"]),
                    lower_hit_threshold=int(self.settings[arm]
                                            ["lower_hit_threshold"]),
                    outp=output_level
                )
                # filter the species hits next
                filtered = mip.filter_bowtie(
                    filtered_host,
                    self.primers["bt_filtered"][arm]["filename"],
                    self.primer3_output_DIR,
                    self.species,
                    TM=float(self.settings[arm]["filter_tm"]),
                    hit_threshold=int(self.settings[arm]["hit_threshold"]),
                    lower_tm=float(self.settings[arm]["lower_filter_tm"]),
                    lower_hit_threshold=int(self.settings[arm]
                                            ["lower_hit_threshold"]),
                    outp=output_level
                )
                # add number of remaining primers after primers with high TM
                # off target alignments were removed
                self.primers["bowtied"][arm]["primers_remaining"] = len(
                    filtered["primer_information"])
                self.primers["bowtied"][arm]["primers_remaining_host"] = len(
                    filtered_host["primer_information"])
            else:
                # if there is no host species, filter the species bowtie hits
                filtered = mip.filter_bowtie(
                    bowtied, self.primers["bt_filtered"][arm]["filename"],
                    self.primer3_output_DIR,
                    self.species,
                    TM=float(self.settings[arm]["filter_tm"]),
                    hit_threshold=int(self.settings[arm]["hit_threshold"]),
                    lower_tm=float(self.settings[arm]["lower_filter_tm"]),
                    lower_hit_threshold=int(self.settings[arm]
                                            ["lower_hit_threshold"]))
                # add number of remaining primers after primers with high TM
                # off target alignments were removed
                self.primers["bowtied"][arm]["primers_remaining"] = len(
                    filtered["primer_information"])
                self.primers["bowtied"][arm]["primers_remaining_host"] = len(
                    filtered["primer_information"])

            # Pick the best alternative primer according to the TM information,
            # if necessary. when there is no need for alternatives, the
            # alternative function will return the same dictionary as the
            # input.
            try:
                tm_diff = float(self.settings[arm]["alt_tm_diff"])
            except KeyError:
                tm_diff = float(self.settings[arm]["tm_diff"])
            self.primers["alternative"][arm] = {"filename": name
                                                + "_alternative"}
            alternative = mip.alternative(
                filtered, self.primers["alternative"][arm]["filename"],
                self.primer3_output_DIR, tm_diff, outp=1)
            self.primers["alternative"][arm]["dictionary"] = alternative
        return

    def score_primers(self):
        """
        Score and filter primers.

        Scoring and filtering parameters are specified in the rinfo file.
        """
        self.primers["scored"] = {}
        self.primers["filtered"] = {}
        mask_penalty = self.scoring["mask_penalty"]
        for arm in ["extension", "ligation"]:
            name = self.fullname + "_" + arm[:3]
            self.primers["scored"][arm] = {"filename": name + "_scored"}
            self.primers["filtered"][arm] = {"filename": name + "_filtered"}
            scored = mip.score_paralog_primers(
                self.primers["alternative"][arm]["dictionary"],
                self.primers["scored"][arm]["filename"],
                self.primer3_output_DIR, arm, mask_penalty, self.species,
                mip.mip_backbones[self.settings["mip"]["backbone"]],
                outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
            )
            filtered = mip.filter_primers(
                scored, self.primers["filtered"][arm]["filename"],
                self.primer3_output_DIR, int(self.settings[arm]["pick_size"]),
                int(self.settings[arm]["bin_size"]),
                outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
            )
            self.primers["filtered"][arm]["dictionary"] = filtered
        return

    def pick_primer_pairs(self):
        """Pick primer pairs satisfying min/max product size criteria."""
        self.primers["paired"] = {"filename": self.fullname + "_paired"}
        settings = self.settings["mip"]
        self.primers["pairs"] = {"filename": self.fullname + "_pairs"}
        # determine product size range by using the planned sequence  read
        # lengths, UMI lengths and desired read overlap.
        extension_read_len = int(self.settings["extension"]["read_length"])
        ligation_read_len = int(self.settings["ligation"]["read_length"])
        extension_umi_len = int(self.settings["extension"]["umi_length"])
        ligation_umi_len = int(self.settings["ligation"]["umi_length"])
        minimum_read_overlap = int(self.settings["mip"]
                                   ["minimum_read_overlap"])
        maximum_read_overlap = int(self.settings["mip"]
                                   ["maximum_read_overlap"])
        size_min = (extension_read_len + ligation_read_len - extension_umi_len
                    - ligation_umi_len - maximum_read_overlap)
        size_max = (extension_read_len + ligation_read_len - extension_umi_len
                    - ligation_umi_len - minimum_read_overlap)
        paired = mip.pick_paralog_primer_pairs(
            self.primers["filtered"]["extension"]["dictionary"],
            self.primers["filtered"]["ligation"]["dictionary"],
            self.primers["paired"]["filename"],
            self.primer3_output_DIR,
            size_min,
            size_max,
            settings["alternative_arms"],
            self.locus.insertions,
            self.subregion_name,
            outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
        )
        if paired == 1:
            return 1
        pairs = mip.add_capture_sequence(
            paired,
            self.primers["pairs"]["filename"],
            self.primer3_output_DIR,
            self.species)
        self.primers["pairs"]["dictionary"] = pairs
        return 0

    def make_mips(self):
        """Create MIP sequences by combining arms and the backbone."""
        self.mips = {"original": {}}
        self.primers["original"]["mip"] = {}
        self.primers["original"]["mip"]["filename"] = self.fullname + "_mips"
        mips_made = mip.make_mips(
            self.primers["pairs"]["dictionary"],
            self.primers["original"]["mip"]["filename"],
            self.primer3_output_DIR,
            self.mfold_input_DIR,
            mip.mip_backbones[self.settings["mip"]["backbone"]],
            outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
        )
        self.primers["original"]["mip"]["dictionary"] = mips_made
        mip_bb = mip.mip_backbones[self.settings["mip"]["backbone"]]
        for m in mips_made["pair_information"]:
            temp = mips_made["pair_information"][m]
            self.mips["original"][m] = Mip(self, m, temp, mip_bb)
        return

    def hairpin(self):
        """Calculate MIP hairpin TMs."""
        self.primers["hairpin"] = {}
        self.primers["hairpin"]["filename"] = self.fullname + "_hp"
        self.mips["hairpin"] = {}
        hp = mip.check_hairpin(
            self.primers["original"]["mip"]["dictionary"],
            self.primers["hairpin"]["filename"],
            self.settings,
            self.primer3_output_DIR,
            outp=int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"])
        )
        self.primers["hairpin"]["dictionary"] = hp
        for h in hp["pair_information"]:
            hairpin_tms = hp["pair_information"][h]["hairpin"]
            self.mips["original"][h].hairpin = hairpin_tms
            self.mips["hairpin"][h] = self.mips["original"][h]
        return

    def score_mips(self):
        """Score each designed MIP object."""
        for m in list(self.mips["hairpin"].keys()):
            self.mips["hairpin"][m].add_capture_info()
        for m in list(self.mips["hairpin"].keys()):
            self.mips["hairpin"][m].score_mip_object()
        try:
            selection_skip = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["skip"])
        except KeyError:
            selection_skip = False
        if selection_skip:
            self.skip_compatibility = True
        else:
            self.skip_compatibility = False
        temp_scored = copy.deepcopy(self.primers["hairpin"]["dictionary"])
        self.mips["scored_filtered"] = {"dictionary": temp_scored,
                                        "filename": self.fullname
                                        + "_scored_filtered"}

        with open(self.primer3_output_DIR + self.mips["scored_filtered"][
                "filename"], "w") as infile:
            json.dump(self.mips["scored_filtered"]["dictionary"], infile,
                      indent=1)
        return

    def compatible(self):
        """Find compatible MIP sets."""
        try:
            overlap_same = int(self.locus.rinfo["SELECTION"]["compatibility"]
                               ["same_strand_overlap"])
            overlap_opposite = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["opposite_strand_overlap"])
        except KeyError:
            overlap_same = overlap_opposite = 0
        compatible_chains = int(self.locus.rinfo["SELECTION"][
            "compatibility"]["chain"])
        if compatible_chains:
            self.chain_mips = True
        else:
            self.chain_mips = False
        self.primers["compatible"] = {
            "filename": self.fullname + "_compat",
            "listname": self.fullname + "_compatibles"
        }
        if not self.skip_compatibility:
            bin_size = int(self.locus.rinfo["SELECTION"]["compatibility"][
                "bin_size"])
            trim_increment = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["trim_increment"])
            trim_limit = int(self.locus.rinfo["SELECTION"]["compatibility"][
                "trim_limit"])
            mip_limit = int(self.locus.rinfo["SELECTION"]["compatibility"][
                "mip_limit"])
            compatible_sets = mip.compatible_chains(
                self.mips["scored_filtered"]["filename"],
                self.mips["hairpin"],
                self.primer3_output_DIR,
                self.primers["compatible"]["filename"],
                self.primers["compatible"]["listname"],
                int(self.scoring["must_bonus"]),
                int(self.scoring["set_copy_bonus"]),
                overlap_same,
                overlap_opposite,
                int(self.locus.rinfo["CAPTURE"]["S0"]["output_level"]),
                bin_size, trim_increment, trim_limit, mip_limit,
                self.chain_mips, self.intervals)
            self.primers["compatible"]["list"] = compatible_sets
        return

    def best_mipset(self):
        """Find best scoring MIP set."""
        self.mips["best_mipset"] = {"filename": self.fullname + "_best_set",
                                    "dictionary": {"mips": {}}}
        try:
            overlap_same = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["same_strand_overlap"])
            overlap_opposite = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["opposite_strand_overlap"])
        except KeyError:
            overlap_same = overlap_opposite = 0
        try:
            max_overlap_same = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["max_same_strand_overlap"])
            max_overlap_opposite = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["max_opposite_strand_overlap"])
        except KeyError:
            max_overlap_same = overlap_same
            max_overlap_opposite = overlap_opposite
        intervals = [self.begin + self.flank, self.end - self.flank]
        subregion_size = intervals[1] - intervals[0] + 1
        best_set_chained = False
        chain_bonus = self.scoring["chain_bonus"]
        coverage_threshold = self.scoring["chain_coverage"]
        must_bonus = self.scoring["must_bonus"]
        set_copy_bonus = self.scoring["set_copy_bonus"]
        scored_filtered = (self.mips["scored_filtered"]
                           ["dictionary"]["pair_information"])
        if not self.skip_compatibility:
            # load dict file that has the compatibility information and mip
            # information
            with open(self.primer3_output_DIR + self.primers["compatible"][
                    "filename"], "rb") as infile:
                mip_dic = pickle.load(infile)
            # if there is any sets in the list
            if not self.single:
                all_mip_sets = [list(s) for s in self.primers[
                    "compatible"]["list"]]
                for mip_set in all_mip_sets:
                    # create a dic for diffs captured cumulatively by all
                    # mips in the set
                    merged_caps = []
                    # create a list for mip scores based on mip sequence and
                    # not the captured diffs
                    mip_scores = []
                    # create a list for what is captured by the set (only must
                    # captures)
                    must_captured = []
                    # create a list for other targets captured
                    targets_captured = []
                    # a list for mip coordinates
                    capture_coordinates = []
                    for mip_key in mip_set:
                        # extract the mip name
                        # extract the captured diffs from the mip_dic and
                        # append to capture list
                        mip_obj = self.mips["hairpin"][mip_key]
                        uniq = mip_obj.capture_info["unique_captures"]
                        merged_caps.extend(uniq)
                        must_captured.extend(mip_obj.captures)
                        targets_captured.extend(mip_obj.captured_targets)
                        if ((mip_obj.tech_score > 0)
                                and (mip_obj.func_score > 0)):
                            mip_scores.append(
                                float(mip_obj.tech_score * mip_obj.func_score)
                                / 1000
                            )
                        else:
                            mip_scores.append(
                                float(mip_obj.tech_score + mip_obj.func_score)
                                / 1000
                            )
                        mcoord = sorted(
                            [mip_obj.extension["C0"]["GENOMIC_START"],
                             mip_obj.ligation["C0"]["GENOMIC_START"],
                             mip_obj.extension["C0"]["GENOMIC_END"],
                             mip_obj.ligation["C0"]["GENOMIC_END"]]
                        )
                        capture_coordinates.append([mcoord[1] + 1,
                                                    mcoord[2] - 1])
                    merged_capture_coordinates = mip.merge_overlap(
                        capture_coordinates
                    )
                    coverage_size = 0
                    for covered_region in merged_capture_coordinates:
                        target_overlap = mip.overlap(intervals, covered_region)
                        try:
                            coverage_size += (target_overlap[1]
                                              - target_overlap[0] + 1)
                        except IndexError:
                            coverage_size += 0
                    chain_coverage = float(coverage_size) / subregion_size
                    cb = chain_coverage * chain_bonus
                    scp = len(set(merged_caps)) * set_copy_bonus
                    must_set = list(set(must_captured))
                    mb = len(must_set) * must_bonus
                    total_score = mb + cb + scp + sum(mip_scores)
                    try:
                        if total_score > best_score:
                            best_score = total_score
                            best_set = mip_set
                            best_caps = list(set(merged_caps))
                    except NameError:
                        best_score = total_score
                        best_set = mip_set
                        best_caps = list(set(merged_caps))
                # check if mips need to be chained and
                # run the while loop for chained mips
                best_set_chained = False
                cover_attempt = 0
                new_mip_found = False
                while (self.chain_mips
                       and not best_set_chained
                       and cover_attempt <= 500):
                    # check coverage of best set found
                    capture_coordinates = []
                    for mip_key in best_set:
                        mip_obj = self.mips["hairpin"][mip_key]
                        mcoord = sorted(
                            [mip_obj.extension["C0"]["GENOMIC_START"],
                             mip_obj.ligation["C0"]["GENOMIC_START"],
                             mip_obj.extension["C0"]["GENOMIC_END"],
                             mip_obj.ligation["C0"]["GENOMIC_END"]])
                        capture_coordinates.append([mcoord[1] + 1,
                                                    mcoord[2] - 1])
                    merged_capture_coordinates = mip.merge_overlap(
                        capture_coordinates
                    )
                    coverage_size = 0
                    for covered_region in merged_capture_coordinates:
                        target_overlap = mip.overlap(intervals, covered_region)
                        try:
                            coverage_size += (target_overlap[1]
                                              - target_overlap[0] + 1)
                        except IndexError:
                            coverage_size += 0
                    chain_coverage = float(coverage_size) / subregion_size
                    if chain_coverage >= coverage_threshold:
                        best_set_chained = True
                    else:
                        if not new_mip_found:
                            if overlap_same < max_overlap_same:
                                overlap_same += 2
                            elif overlap_opposite < max_overlap_opposite:
                                overlap_opposite += 2
                            else:
                                break
                        cover_attempt += 1
                        best_copy = copy.deepcopy(best_set)
                        # check coverage of best set found
                        capture_coordinates_copy = []
                        for mip_key in best_copy:
                            mip_obj = self.mips["hairpin"][mip_key]
                            mcoord = sorted(
                                [mip_obj.extension["C0"]["GENOMIC_START"],
                                 mip_obj.ligation["C0"]["GENOMIC_START"],
                                 mip_obj.extension["C0"]["GENOMIC_END"],
                                 mip_obj.ligation["C0"]["GENOMIC_END"]])
                            capture_coordinates_copy.append([mcoord[1] + 1,
                                                             mcoord[2] - 1])
                        merged_capture_coordinates_copy = mip.merge_overlap(
                            capture_coordinates_copy
                        )
                        new_mip_found = False
                        for new_mip_key in scored_filtered:
                            if new_mip_key in best_copy:
                                continue
                            new_mip_obj = self.mips["hairpin"][new_mip_key]
                            # check this mip for compatibility to
                            # other mips in the best mipset
                            for bm in best_copy:
                                bm_obj = self.mips["hairpin"][bm]
                                if not mip.compatible_mip_check(
                                    new_mip_obj, bm_obj, overlap_same,
                                    overlap_opposite
                                ):
                                    break
                            else:
                                mcoord = sorted(
                                    [new_mip_obj.extension["C0"]
                                     ["GENOMIC_START"],
                                     new_mip_obj.ligation["C0"]
                                     ["GENOMIC_START"],
                                     new_mip_obj.extension["C0"]
                                     ["GENOMIC_END"],
                                     new_mip_obj.ligation["C0"]
                                     ["GENOMIC_END"]]
                                )
                                cap_coord = [mcoord[1] + 1,
                                             mcoord[2] - 1]
                                so = mip.subtract_overlap(
                                    [cap_coord],
                                    merged_capture_coordinates_copy
                                )
                                if len(so) < 1:
                                    continue
                                mip_set = copy.deepcopy(best_copy)
                                mip_set.append(new_mip_key)
                                # create a dic for diffs captured cumulatively
                                # by all mips in the set
                                merged_caps = []
                                # create a list for mip scores based on mip
                                # sequence and not the captured diffs
                                mip_scores = []
                                # create a list for what is captured by the
                                # set (only must captures)
                                must_captured = []
                                # create a list for other targets captured
                                targets_captured = []
                                # a list for mip coordinates
                                capture_coordinates = []
                                for mip_key in mip_set:
                                    mip_obj = self.mips["hairpin"][mip_key]
                                    uniq = mip_obj.capture_info[
                                        "unique_captures"]
                                    merged_caps.extend(uniq)
                                    must_captured.extend(mip_obj.captures)
                                    targets_captured.extend(
                                        mip_obj.captured_targets)
                                    if (mip_obj.tech_score > 0
                                            and mip_obj.func_score > 0):
                                        mip_scores.append(float(
                                            mip_obj.tech_score
                                            * mip_obj.func_score)/1000)
                                    else:
                                        mip_scores.append(float(
                                            mip_obj.tech_score
                                            + mip_obj.func_score)/1000)
                                    mcoord = sorted(
                                        [mip_obj.extension["C0"][
                                         "GENOMIC_START"],
                                         mip_obj.ligation["C0"][
                                         "GENOMIC_START"],
                                         mip_obj.extension["C0"][
                                         "GENOMIC_END"],
                                         mip_obj.ligation["C0"][
                                         "GENOMIC_END"]])
                                    capture_coordinates.append(
                                        [mcoord[1] + 1, mcoord[2] - 1])
                                merged_capture_coordinates = mip.merge_overlap(
                                    capture_coordinates)
                                coverage_size = 0
                                for covered_region in merged_capture_coordinates:
                                    target_overlap = mip.overlap(
                                        intervals, covered_region)
                                    try:
                                        coverage_size += (target_overlap[1]
                                                          - target_overlap[0]
                                                          + 1)
                                    except IndexError:
                                        coverage_size += 0
                                chain_coverage = float(
                                    coverage_size) / subregion_size
                                cb = chain_coverage * chain_bonus
                                scp = len(set(merged_caps)) * set_copy_bonus
                                must_set = list(set(must_captured))
                                mb = len(must_set) * must_bonus
                                total_score = mb + cb + scp + sum(mip_scores)
                                if total_score > best_score:
                                    new_mip_found = True
                                    best_score = total_score
                                    best_set = mip_set
                                    best_caps = list(set(merged_caps))
                for mip_key in best_set:
                    mip_obj = self.mips["hairpin"][mip_key]
                    self.mips["best_mipset"]["dictionary"]["mips"][mip_key] = (
                        mip_obj)
                self.mips["best_mipset"]["dictionary"]["caps"] = best_caps
                self.mips["best_mipset"]["dictionary"]["score"] = best_score
            # if a single MIP should be selected instead of a set of MIPs
            # select the best scoring MIP
            elif (self.single
                  and (len(list(mip_dic["pair_information"].keys())) > 0)):
                for mip_key in list(mip_dic["pair_information"].keys()):
                    mip_obj = self.mips["hairpin"][mip_key]
                    uniq = mip_obj.capture_info["unique_captures"]
                    if mip_obj.tech_score != 0 and mip_obj.func_score != 0:
                        total_score = mip_obj.tech_score * mip_obj.func_score
                    else:
                        total_score = mip_obj.tech_score + mip_obj.func_score
                    try:
                        if total_score > best_score:
                            best_score = total_score
                            best_mip = mip_key
                            best_caps = uniq
                    except NameError:
                        best_score = total_score
                        best_mip = mip_key
                        best_caps = uniq
                self.mips["best_mipset"]["dictionary"]["mips"][best_mip] = (
                    self.mips["original"][best_mip])
                self.mips["best_mipset"]["dictionary"]["caps"] = best_caps
                self.mips["best_mipset"]["dictionary"]["score"] = best_score
                best_set_chained = False
            else:
                print("No mips available for target region %s" % self.fullname)
        # if compatibility test should be skipped
        else:
            with open(self.primer3_output_DIR + self.mips["scored_filtered"][
                    "filename"], "r") as infile:
                scored_dic = json.load(infile)
            mip_dic = scored_dic["pair_information"]
            # create a dic for diffs captured cumulatively by all mips in the
            # set
            merged_caps = []
            # create a list for mip scores based on mip sequence and not the
            # captured diffs
            mip_scores = []
            mip_coordinates = []
            for mip_key in mip_dic:
                # extract the captured diffs from the mip_dic and append to
                # capture list
                mip_obj = self.mips["original"][mip_key]
                uniq = mip_obj.capture_info["unique_captures"]
                merged_caps.extend(uniq)
                if mip_obj.tech_score != 0 and mip_obj.func_score != 0:
                    mip_scores.append(float(mip_obj.tech_score
                                      * mip_obj.func_score)/1000)
                else:
                    mip_scores.append(float(mip_obj.tech_score
                                      + mip_obj.func_score)/1000)
                mip_scores.append(mip_obj.tech_score * mip_obj.func_score)
                mcoor = sorted([mip_obj.extension["C0"]["GENOMIC_START"],
                                mip_obj.extension["C0"]["GENOMIC_END"],
                                mip_obj.ligation["C0"]["GENOMIC_START"],
                                mip_obj.ligation["C0"]["GENOMIC_END"]])
                mip_coordinates.append([mcoor[1]+1, mcoor[2]+1])
            # calculate total score of mip set
            total_score = ((len(set(merged_caps))**2) + 1) * sum(mip_scores)
            self.mips["best_mipset"]["dictionary"]["mips"] = {}
            for smip in scored_dic["pair_information"]:
                self.mips["best_mipset"]["dictionary"]["mips"][smip] = (
                    self.mips["hairpin"][smip])
            self.mips["best_mipset"]["dictionary"]["caps"] = list(set(
                merged_caps))
            self.mips["best_mipset"]["dictionary"]["score"] = total_score
            # get the overlap of mip coordinates
            # the overlap should have only a single coordinate
            # if all the mips in the set is overlapping
            chained_mips = mip.merge_overlap(mip_coordinates)
            if (len(chained_mips) == 1) and (len(mip_coordinates) > 1):
                best_set_chained = True
            else:
                best_set_chained = False
        self.chained_mips = best_set_chained
        # remove MIPs that do not capture any targets
        # only if capture type is "targets"
        try:
            selection_skip = int(self.locus.rinfo["SELECTION"][
                "compatibility"]["skip"])
        except KeyError:
            selection_skip = False
        if (self.capture["capture_type"] == "targets") and not selection_skip:
            best_mipset = self.mips["best_mipset"]["dictionary"]["mips"]
            for mip_key in list(best_mipset.keys()):
                if len(best_mipset[mip_key].capture_info[
                        "captured_targets"]) == 0:
                    best_mipset.pop(mip_key)
        return


class Mip():
    """
    Mip class provides a construct to represent a MIP probe.

    It is instantiated by a subregion object.
    """

    def __init__(self, subregion, name, mip_dic, backbone):
        try:
            self.extension = mip_dic["extension_primer_information"][
                "PARALOG_COORDINATES"]
            self.ligation = mip_dic["ligation_primer_information"][
                "PARALOG_COORDINATES"]
        except KeyError:
            self.extension = {"C0": mip_dic["extension_primer_information"]}
            self.ligation = {"C0": mip_dic["ligation_primer_information"]}
        self.mip_dic = mip_dic
        self.mip = mip_dic["pairs"]
        self.name = name
        self.backbone = backbone
        self.sequence = mip_dic["mip_information"]["ref"]["SEQUENCE"]
        self.subregion = subregion
        self.species = self.subregion.species
        self.host = self.subregion.host
        try:
            self.alt_copies = mip_dic["alt_copies"]
        except KeyError:
            self.alt_copies = []
        self.captured_copies = mip_dic["captured_copies"]
        self.captured_copies.extend(self.alt_copies)

    def add_capture_info(self):
        # find what is captured by a mip
        self.capture_info = {"captured_targets": {}, "unique_captures": []}
        cap_targets = self.capture_info["captured_targets"]
        # get target dictionary
        targets = self.subregion.targets
        extension_read_len = int(self.subregion.settings["extension"][
            "read_length"])
        ligation_read_len = int(self.subregion.settings["ligation"][
            "read_length"])
        extension_min_trim = int(self.subregion.settings["extension"][
            "minimum_trim"])
        ligation_min_trim = int(self.subregion.settings["ligation"][
            "minimum_trim"])
        extension_umi_len = int(self.subregion.settings["extension"][
            "umi_length"])
        ligation_umi_len = int(self.subregion.settings["ligation"][
            "umi_length"])
        try:
            arms = self.subregion.capture["arms"]
        except KeyError:
            arms = "any"
        if arms not in ["extension",  "ligation",  "any", "both", "capture"]:
            arms = "capture"
        insertions = self.subregion.locus.insertions
        for target_type in targets:
            for copy_id in targets[target_type]:
                if copy_id in self.mip:
                    copy_mip = self.mip[copy_id]
                    cs = copy_mip["capture_start"]
                    ce = copy_mip["capture_end"]
                    cc = copy_mip["chrom"]
                    co = copy_mip["orientation"]
                    if co == "forward":
                        ext_cs = cs
                        ext_ce = (
                            cs + extension_read_len
                            - extension_umi_len - extension_min_trim)
                        lig_ce = ce
                        lig_cs = (
                            ce - ligation_read_len
                            + ligation_umi_len + ligation_min_trim)
                        if not insertions.empty:
                            ext_insertion_len = insertions.loc[
                                (insertions["copy_chrom"] == cc)
                                & (insertions["copy_begin"] > ext_cs)
                                & (insertions["copy_end"] < ext_ce),
                                "max_size"].sum()
                            lig_insertion_len = insertions.loc[
                                (insertions["copy_chrom"] == cc)
                                & (insertions["copy_begin"] > lig_cs)
                                & (insertions["copy_end"] < lig_ce),
                                "max_size"].sum()
                        else:
                            ext_insertion_len = lig_insertion_len = 0
                    else:
                        lig_cs = cs
                        lig_ce = (
                            cs + ligation_read_len
                            - ligation_umi_len - ligation_min_trim)
                        ext_ce = ce
                        ext_cs = (
                            ce - extension_read_len
                            + extension_umi_len + extension_min_trim)
                        if not insertions.empty:
                            lig_insertion_len = insertions.loc[
                                (insertions["copy_chrom"] == cc)
                                & (insertions["copy_begin"] > lig_cs)
                                & (insertions["copy_end"] < lig_ce),
                                "max_size"].sum()
                            ext_insertion_len = insertions.loc[
                                (insertions["copy_chrom"] == cc)
                                & (insertions["copy_begin"] > ext_cs)
                                & (insertions["copy_end"] < ext_ce),
                                "max_size"].sum()
                        else:
                            ext_insertion_len = lig_insertion_len = 0
                    for t in targets[target_type][copy_id]:
                        tar = targets[target_type][copy_id][t]
                        tbeg = tar["copy_begin"]
                        try:
                            tend = tar["copy_end"]
                        except KeyError:
                            tend = tbeg
                        tchrom = tar["copy_chrom"]
                        try:
                            size_diff = int(tar["size_difference"])
                        except KeyError:
                            size_diff = 0
                        ext_size_diff = max([size_diff, ext_insertion_len])
                        lig_size_diff = max([size_diff, lig_insertion_len])
                        if co == "forward":
                            ext_coverage_start = ext_cs
                            ext_coverage_end = ext_ce - ext_size_diff
                            lig_coverage_start = lig_cs + lig_size_diff
                            lig_coverage_end = lig_ce
                        else:
                            ext_coverage_start = ext_cs + ext_size_diff
                            ext_coverage_end = ext_ce
                            lig_coverage_start = lig_cs
                            lig_coverage_end = lig_ce - lig_size_diff
                        covered_by_capture = False
                        covered_by_ext = False
                        covered_by_lig = False
                        if tchrom == cc:
                            if cs <= tbeg <= tend <= ce:
                                covered_by_capture = True
                            if (ext_coverage_start <= tbeg <= tend
                                    <= ext_coverage_end):
                                covered_by_ext = True
                            if (lig_coverage_start <= tbeg <= tend
                                    <= lig_coverage_end):
                                covered_by_lig = True
                        covered_by_both = covered_by_ext and covered_by_lig
                        covered_by_any = covered_by_ext or covered_by_lig
                        if (((arms == "capture") and covered_by_capture)
                                or ((arms == "extension") and covered_by_ext)
                                or ((arms == "ligation") and covered_by_lig)
                                or (((arms == "both") and covered_by_both))
                                or ((arms == "any") and covered_by_any)):
                            try:
                                cap_targets[
                                    target_type][copy_id][t] = tar
                            except KeyError:
                                try:
                                    cap_targets[
                                        target_type][copy_id] = {t: tar}
                                except KeyError:
                                    cap_targets[target_type] = {
                                        copy_id: {t: tar}
                                    }
        all_captured_sequences = []
        for c in self.captured_copies:
            all_captured_sequences.append(
                self.mip[c]["capture_sequence"]
            )
        for c in self.captured_copies:
            copy_seq = self.mip[c]["capture_sequence"]
            if all_captured_sequences.count(copy_seq) == 1:
                self.capture_info["unique_captures"].append(c)
        return

    def score_mip_object(self):
        """Score a mip object technically and functionally.

        Technical scores estimates how well we think the mip will work.
        Functional scores  reflect how informative the capture will be,
        assuming the mip works.
        """
        # TECHNICAL SCORING
        # technical scoring coefficients were calculated based on
        # linear models of various parameters and provided as a dict
        with open("/opt/resources/mip_scores.dict", "rb") as infile:
            linear_coefs = pickle.load(infile)
        # the model was developed using specific reaction conditions as below.
        # actual conditions may be different from these but we'll use these
        # for the model.
        na = 25  # Sodium concentration
        mg = 10  # magnesium concentration
        conc = 0.04  # oligo concentration
        # arms are scored on the primer level previously
        # get extension arm score
        extension_score = self.mip_dic["extension_primer_information"]["SCORE"]
        # get ligation arm score
        ligation_score = self.mip_dic["ligation_primer_information"]["SCORE"]
        # Calculate probe scores
        # calculate free energy between extension and ligation arms
        extension_arm = self.mip_dic["extension_primer_information"][
            "SEQUENCE"]
        ligation_arm = self.mip_dic["ligation_primer_information"]["SEQUENCE"]
        probe_dg = primer3.calcHeterodimer(
            extension_arm, mip.reverse_complement(ligation_arm),
            mv_conc=na, dv_conc=mg, dna_conc=conc, dntp_conc=0).dg
        # calculate gc content of the captured sequences
        gc_values = []
        for c in self.mip:
            capture_gc = mip.calculate_gc(self.mip[c]["capture_sequence"])
            gc_values.append(capture_gc)
        capture_gc = sum(gc_values)/len(gc_values)
        # create a mip parameter dict
        score_features = {"probe_dg": probe_dg,
                          "capture_gc": capture_gc}
        # calculate technical score using the linear model provided
        tech_score = 0
        for feature in score_features:
            degree = linear_coefs[feature]["degree"]
            mip_feature = score_features[feature]
            poly_feat = [pow(mip_feature, i) for i in range(degree + 1)]
            tech_score += sum(linear_coefs[feature]["coef"] * poly_feat)
            tech_score += linear_coefs[feature]["intercept"]
        tech_score += extension_score
        tech_score += ligation_score
        # The score calculated in the model is generally between 0 and 1
        # older scores were between 0 and 2000. So we'll multiply the new
        # score to be comparable to the old.
        tech_score = round(tech_score * 2000)
        # FUNCTIONAL SCORING
        func_score = 0
        diff_scores = self.subregion.scoring["diff_scores"]
        snp_scores = self.subregion.scoring["snp_scores"]
        diff_scores.update(snp_scores)
        self.captures = []
        self.captured_targets = []
        capped = self.capture_info["captured_targets"]
        for target_type in capped:
            for copy_id in capped[target_type]:
                for t in capped[target_type][copy_id]:
                    func_score += diff_scores[target_type]
                    if target_type == "must":
                        self.captures.append(t)
                    else:
                        self.captured_targets.append(t)
        unique_copy_bonus = self.subregion.scoring["unique_copy_bonus"]
        alternative_copy_penalty = self.subregion.scoring[
            "alternative_copy_penalty"]
        copy_bonus = (len(self.capture_info["unique_captures"])
                      * unique_copy_bonus + (len(self.captured_copies)
                      - len(self.alt_copies)/2.) * alternative_copy_penalty)
        func_score += copy_bonus
        self.tech_score = self.subregion.scoring[
            "technical_score_coefficient"] * tech_score
        self.func_score = func_score
        self.capture_gc = capture_gc
        return

def parasight_mod(resource_dir,  design_info_file, species,
                  designed_gene_list=None, extra_extension=".extra",
                  maf=0.1, height=200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    pdf_dir = os.path.join(resource_dir, "mod_pdfs")
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash", "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                        + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                        "-sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = os.path.join(design_info[t]["design_dir"], t, t)
        showname = basename + ".show"
        try:
            with open(showname) as infile:
                sln = infile.readlines()[-1].strip().split("\t")
                show_region = sln[0] + ":" + sln[2] + "-" + sln[3]
        except IOError:
            continue
        reg_snps = get_snps(show_region,
                            get_file_locations()[species]["snps"])
        indels_low = []
        indels_high = []
        snvs_low = []
        snvs_high = []
        for rsnp in reg_snps:
            acs = [a for a in rsnp[23].split(",") if a != ""]
            allele_count = list(map(int, list(map(float, acs))))
            # allele with max count
            max_all = max(allele_count)
            # total alleles
            tot_all = sum(allele_count)
            # minor allele freq
            min_af = (tot_all - max_all)/float(tot_all)
            if "-" in rsnp[22].split(","):
                if min_af >= maf:
                    indels_high.append(rsnp)
                else:
                    indels_low.append(rsnp)
            elif min_af >= maf:
                snvs_high.append(rsnp)
            else:
                snvs_low.append(rsnp)
        backup_name = basename + ".extra"
        try:
            with open(backup_name) as infile, open(backup_name + ".mod", "w") as outfile:
                thirds = 0
                rd_range = list(range(height))
                for line in infile:
                    newline = line.split("\t")
                    if newline[3] in ["all_targets", "target"]:
                        newline[5] = "-10"
                        newline[6] = "4"
                        outfile.write("\t".join(newline))
                    elif newline[3] in ["capture", "extension", "ligation"]:
                        if (thirds % 3)== 0:
                            rd = random.choice(rd_range)
                        thirds += 1
                        newline[5] = str(-30 - rd)
                        outfile.write("\t".join(newline))
                    elif newline[3] == "snp":
                        pass
                    else:
                        outfile.write(line)
                outfile.write("\n")
                for col, snp_list, snp_type, ofs in zip(["pink", "red",
                                     "light green", "dark green"],
                                   [indels_low, indels_high,
                                   snvs_low, snvs_high],
                                                  ["low frequency indel",
                                                   "high frequency indel",
                                                  "low frequency SNP",
                                                  "high frequency SNP"],
                                                       ["-27", "-27",
                                                       "-26", "-26"]):
                    for snp in snp_list:
                        ol = snp[1:4]
                        ol.extend([snp_type, col, ofs, "4",
                                  ""])
                        outfile.write("\t".join(ol) + "\n")
        except IOError:
            continue
        psname = basename + ".01.01.ps"
        pdfname = basename + ".mod.pdf"
        gs_command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite " +             "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All -sOutputFile="            + pdfname + " " + psname
        if designed_gene_list != None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".mod.pdf")
        else:
            pdf_convert_list.append(t + ".mod.pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".mod.pdf " +                            pdf_dir + t + ".mod.pdf")
        outlist = ["parasight76.pl",
                   "-showseq",
                   basename + ".show",
                   "-extra",
                   basename + extra_extension + ".mod",
                   "-template",
                   script_dir + "resources/nolabel.pst",
                   "-precode file:" + basename + ".precode" ,
                   "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode" , "w") as outfile:
            outfile.write("$opt{'filename'}='" + t +                          "';&fitlongestline; &print_all (0,'" + basename + "')")
    with open(resource_dir + "backup_commands", "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(resource_dir + "parasight_commands", "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(resource_dir + "gs_commands", "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(resource_dir + "copy_commands", "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(resource_dir + "convert_commands", "w") as outfile:
        outfile.write("\n".join(pdf_merge_list))
    visualization_list = ["#!/usr/bin/env bash"]
    visualization_list.append("chmod +x backup_commands")
    visualization_list.append("./backup_commands")
    visualization_list.append("chmod +x parasight_commands")
    visualization_list.append("./parasight_commands")
    visualization_list.append("chmod +x gs_commands")
    visualization_list.append("./gs_commands")
    visualization_list.append("chmod +x copy_commands")
    visualization_list.append("./copy_commands")
    visualization_list.append("chmod +x convert_commands")
    visualization_list.append("./convert_commands")
    with open(resource_dir + "visualize_mod.sh", "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return


def parasight_shift(resource_dir, design_info_file, species,
                    designed_gene_list=None, extra_extension=".extra",
                    height=200):
    with open(design_info_file, "rb") as infile:
        design_info = pickle.load(infile)
    output_list = ["#!/usr/bin/env bash"]
    pdf_dir = os.path.join(resource_dir, "mod_pdfs")
    backup_list = ["#!/usr/bin/env bash"]
    gs_list = ["#!/usr/bin/env bash"]
    pdf_list = ["#!/usr/bin/env bash"]
    pdf_merge_list = ["#!/usr/bin/env bash", "cd " + pdf_dir]
    pdf_convert_list = ["gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                        + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                        "-sOutputFile=merged.pdf"]
    if not os.path.exists(pdf_dir):
        os.makedirs(pdf_dir)
    for t in design_info:
        basename = os.path.join(design_info[t]["design_dir"], t, t)
        backup_name = basename + ".extra"
        try:
            with open(backup_name) as infile, open(
                    backup_name + ".mod", "w") as outfile:
                thirds = 0
                rd_range = list(range(height))
                for line in infile:
                    newline = line.split("\t")
                    if newline[3] in ["all_targets", "target"]:
                        newline[5] = "-10"
                        newline[6] = "4"
                        outfile.write("\t".join(newline))
                    elif newline[3] in ["capture", "extension", "ligation"]:
                        if (thirds % 3) == 0:
                            rd = random.choice(rd_range)
                        thirds += 1
                        newline[5] = str(-30 - rd)
                        outfile.write("\t".join(newline))
                    else:
                        outfile.write(line)
                outfile.write("\n")
        except IOError:
            continue
        psname = basename + ".01.01.ps"
        pdfname = basename + ".mod.pdf"
        gs_command = ("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite "
                      + "-dPDFSETTINGS=/prepress -dAutoRotatePages=/All "
                      "-sOutputFile=" + pdfname + " " + psname)
        if designed_gene_list is not None:
            if t in designed_gene_list:
                pdf_convert_list.append(t + ".mod.pdf")
        else:
            pdf_convert_list.append(t + ".mod.pdf")
        gs_list.append(gs_command)
        pdf_list.append("cp " + basename + ".mod.pdf "
                        + pdf_dir + t + ".mod.pdf")
        outlist = ["parasight76.pl",
                   "-showseq", basename + ".show",
                   "-extra", basename + extra_extension + ".mod",
                   "-template", "/opt/resources/nolabel.pst",
                   "-precode file:" + basename + ".precode", "-die"]
        output_list.append(" ".join(outlist))
        with open(basename + ".precode", "w") as outfile:
            outfile.write("$opt{'filename'}='" + t
                          + "';&fitlongestline; &print_all (0,'"
                          + basename + "')")
    with open(resource_dir + "backup_commands", "w") as outfile:
        outfile.write("\n".join(backup_list))
    with open(resource_dir + "parasight_commands", "w") as outfile:
        outfile.write("\n".join(output_list))
    with open(resource_dir + "gs_commands", "w") as outfile:
        outfile.write("\n".join(gs_list))
    with open(resource_dir + "copy_commands", "w") as outfile:
        outfile.write("\n".join(pdf_list))
    pdf_merge_list.append(" ".join(pdf_convert_list))
    with open(resource_dir + "convert_commands", "w") as outfile:
        outfile.write("\n".join(pdf_merge_list))
    visualization_list = ["#!/usr/bin/env bash"]
    visualization_list.append("chmod +x backup_commands")
    visualization_list.append("./backup_commands")
    visualization_list.append("chmod +x parasight_commands")
    visualization_list.append("./parasight_commands")
    visualization_list.append("chmod +x gs_commands")
    visualization_list.append("./gs_commands")
    visualization_list.append("chmod +x copy_commands")
    visualization_list.append("./copy_commands")
    visualization_list.append("chmod +x convert_commands")
    visualization_list.append("./convert_commands")
    with open(resource_dir + "visualize_mod.sh", "w") as outfile:
        outfile.write("\n".join(visualization_list))
    return


def rescue_mips(design_dir,
                paralog_name,
                redesign_list,
               same_strand_overlap,
               opposite_strand_overlap,
               low,
               high,
                mip_limit,
               chain,
               incomp_score):
    print(("Redesigning MIPs for ", paralog_name))
    paralog_file = design_dir + paralog_name + "/" + paralog_name
    with open(paralog_file, "rb") as infile:
        par = pickle.load(infile)
    subprocess.call(["scp",
                     paralog_file,
                     paralog_file + ".last"])
    par.extra_mips = {}
    redesign_pairs = []
    for segment_name in par.segments:
        seg = par.segments[segment_name]
        """
        seg.rinfo["SELECTION"]["compatibility"]\
        ["low"] = low
        seg.rinfo["SELECTION"]["compatibility"]\
        ["high"] = high
        seg.rinfo["SELECTION"]["compatibility"]\
        ["chain"] = chain

        """

        seg.rinfo["SELECTION"]["compatibility"]        ["low"] = low
        seg.rinfo["SELECTION"]["compatibility"]        ["high"] = high
        seg.rinfo["SELECTION"]["compatibility"]        ["mip_limit"] = mip_limit
        for subregion_name in seg.subregions:
            sub = seg.subregions[subregion_name]
            sub.score_mips()
            temp_scored = copy.deepcopy(sub.primers["hairpin"]["dictionary"])
            best_mipset = sub.mips["best_mipset"]["dictionary"]["mips"]
            scored = list(sub.mips["scored_filtered"]["dictionary"]                            ["pair_information"].keys())
            keep_keys = scored + list(best_mipset.keys())
            for m in list(temp_scored['pair_information'].keys()):
                if m not in keep_keys:
                    temp_scored['pair_information'].pop(m)
            compatible = compatibility(temp_scored,
                                       overlap_same = same_strand_overlap,
                                    overlap_opposite = opposite_strand_overlap)
            alt = sub.mips["best_mipset"]["dictionary"]["alternatives"] = {}
            hairpin = sub.mips["hairpin"]
            compat_info = compatible["pair_information"]
            for m in best_mipset:
                if best_mipset[m].fullname in redesign_list:
                    if m not in redesign_pairs:
                        redesign_pairs.append(m)
                    capture_start = sub.mips["best_mipset"]["dictionary"]                    ["mips"][m].mip["C0"]["capture_start"]
                    capture_end = sub.mips["best_mipset"]["dictionary"]                    ["mips"][m].mip["C0"]["capture_end"]
                    capture_target = [capture_start, capture_end]
                    alt_count = 0
                    alt[m] = {}
                    for p1 in scored:
                        if p1 != m and p1 not in redesign_list:
                            #print p1
                            start1 = hairpin[p1].mip["C0"]["capture_start"]
                            end1 = hairpin[p1].mip["C0"]["capture_end"]
                            capture1 = [start1, end1]
                            if len(overlap(capture_target, capture1)) > 0:
                                remaining_target = subtract_overlap(
                                    [capture_target],
                                    [capture1])
                                if len(remaining_target) == 1:
                                    for p2 in scored:
                                        if p2 != p1 and p2 != m and                                        p1 not in redesign_list and                                        (p2 not in compat_info[p1]["incompatible"]):
                                            start2 = hairpin[p2].mip["C0"]                                            ["capture_start"]
                                            end2 = hairpin[p2].mip["C0"]                                            ["capture_end"]
                                            capture2 = [start2, end2]
                                            if len(subtract_overlap(
                                                [capture_target],
                                                [capture2])) == 1 and\
                                               len(overlap(remaining_target[0],
                                                               capture2)) > 0:
                                                uncovered = subtract_overlap(
                                                    remaining_target,
                                                    [capture2])
                                                if len(uncovered) == 0:
                                                    alt[m][str(alt_count)] = [p1,p2]
                                                    alt_count += 1

            new_set = list(best_mipset.keys())
            for m in best_mipset:
                if best_mipset[m].fullname in redesign_list:
                    print(("Rescue MIPs for ", best_mipset[m].fullname))
                    coverage_needed = [hairpin[m].mip["C0"]                                            ["capture_start"],
                                             hairpin[m].mip["C0"]\
                                            ["capture_end"]]
                    print(("Region to cover is ", coverage_needed))
                    best_rescue = []
                    best_score = -50000
                    for p in alt[m]:
                        p1, p2 = alt[m][p]
                        incomp = len(set(compat_info[p1]["incompatible"] +                                      compat_info[p2]["incompatible"])                        .intersection(new_set))
                        score = 0
                        score += hairpin[p1].tech_score
                        score += hairpin[p2].tech_score
                        score += incomp_score * incomp
                        if score > best_score:
                            best_score = score
                            best_rescue = [p1, p2]
                            if p1 in best_mipset:
                                best_rescue = [p2]
                            if p2 in best_mipset:
                                best_rescue = [p1]
                    resc_capture = []
                    for resc in best_rescue:
                        resc_capture.append([hairpin[resc].mip["C0"]                                            ["capture_start"],
                                             hairpin[resc].mip["C0"]\
                                            ["capture_end"]])
                    print(("Rescue mip coverage is ", merge_overlap(resc_capture)))
                    print(("Rescue score is ", best_score))
                    new_set.extend(best_rescue)
            sub.extra_mips = extra_mips = {}
            for t in new_set:
                if t not in best_mipset:
                    extra_mips[t] = hairpin[t]
            par.extra_mips.update(sub.extra_mips)
    locus_info = par.locus_info
    selected_mips = par.selected_mips
    mip_names_ordered = sorted(selected_mips,
                               key=lambda mip_key: selected_mips[mip_key].\
                                                    mip["C0"]["mip_start"])
    name_counter = int(selected_mips[mip_names_ordered[-1]]                       .fullname.split("_")[-1][3:]) + 1
    extra_names_ordered = sorted(extra_mips,
                               key=lambda mip_key: extra_mips[mip_key].\
                                                    mip["C0"]["mip_start"])
    outfile_list = []
    for mip_name in extra_names_ordered:
        m = extra_mips[mip_name]
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
        for key in m.mip_dic["mip_information"]:
            if key == "ref":
                fullname_extension = "_ref"
            elif key.startswith("alt"):
                fullname_extension = "_" + key
            else:
                continue
            # ["#pair_name", "mip_name", "chrom", "mip_start", "capture_start",
            # "capture_end", "mip_end", "orientation", "tech_score", "func_score","mip_sequence",
            # "unique_captures", "must_captured", "captured_copies"]
            outlist = [m.name, m.fullname + fullname_extension,
                       m.chromosome, m.mip_start,
                       m.capture_start, m.capture_end,
                       m.mip_end, m.orientation,
                        m.tech_score, m.func_score,
                        m.mip_dic["mip_information"][key]["SEQUENCE"]]
            locus_info["mips"].append(outlist)
            outfile_list.append(outlist)
        name_counter += 1
    for mipname in redesign_pairs:
        m = selected_mips[mipname]
        for key in m.mip_dic["mip_information"]:
            if key == "ref":
                fullname_extension = "_ref"
            elif key.startswith("alt"):
                fullname_extension = "_" + key
            else:
                continue
            outlist = ["#" + m.name, m.fullname + fullname_extension,
                       m.chromosome, m.mip_start,
                           m.capture_start, m.capture_end,
                       m.mip_end, m.orientation,
                            m.tech_score, m.func_score,
                            m.mip_dic["mip_information"][key]["SEQUENCE"]]
            outfile_list.append(outlist)
    par.print_info()
    with open(paralog_file + ".rescued", "wb") as outfile:
        pickle.dump(par, outfile)
    with open(design_dir + paralog_name + "/" + paralog_name +              "_rescue.mips", "w") as outfile:
        outfile.write("\n".join(["\t".join(map(str, outlist)) for                                 outlist in outfile_list]))
    return
def plot_clusters(settings, cluster_method="dbscan"):
    wdir = settings["workingDir"]
    if cluster_method == "dbscan":
        cluster_output_file = os.path.join(wdir, settings["dbScanOutputFile"])
    else :
        cluster_output_file = os.path.join(wdir, settings["clusterOutputFile"])
    with open(cluster_output_file, "rb") as infile:
        cnv_calls = pickle.load(infile)
    filtered_tables_file = os.path.join(wdir, settings["filteredTablesFile"])
    with open(filtered_tables_file, "rb") as infile:
        all_tables = pickle.load(infile)
    figsize = tuple(map(int, settings["figsize"]))
    ymax = int(settings["ymax"])
    if cluster_method == "dbscan":
        image_dir = os.path.join(wdir, "dbscan_images")
    else:
        image_dir = os.path.join(wdir, "cluster_images")
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    for gene_name in cnv_calls:
        clusters = cnv_calls[gene_name]["clusters"]
        tables = all_tables[gene_name]
        probe_table = tables["probe_information"]
        copies = [probe_table[0,2]]
        locations = [0]
        for i in range(len(probe_table) -1):
            current_copy = probe_table[i+1, 2]
            if current_copy != copies[-1]:
                copies.append(current_copy)
                locations.extend([i, i+1])
        locations.append(len(probe_table))
        locations = [[locations[i], locations[i+1]] for i in range(0, len(locations), 2)]
        upper_limit = cnv_calls[gene_name]["plotting"]["upper_limit"]
        lower_limit = cnv_calls[gene_name]["plotting"]["lower_limit"]
        fig1 = plt.figure(figsize=figsize)
        fig2 = plt.figure(figsize=figsize)
        fig3 = plt.figure(figsize=figsize)
        for c in range(len(clusters)):
            ax1 = fig1.add_subplot(len(clusters), 1, c)
            ax2 = fig2.add_subplot(len(clusters), 1, c)
            ax3 = fig3.add_subplot(len(clusters), 1, c)
            #ax1 = axarray[c]
            #fig, ax1 = plt.subplots()
            plot_data = clusters[c]["cluster_data"]
            median_data = clusters[c]["cluster_medians"]
            collapsed_median = clusters[c]["collapsed_medians"]
            odds = clusters[c]["cluster_stats"]["fisher"]["OR"]
            fish_p = clusters[c]["cluster_stats"]["fisher"]["pvalue"]
            chi_p = clusters[c]["cluster_stats"]["chi"]["pvalue"]
            cluster_size = len(clusters[c]["cluster_samples"])
            cluster_summary_list = ["cluster size " + str(cluster_size),
                                    "odds ratio " + str(odds),
                                    "fisher's exact pvalue " + str(fish_p),
                                    "chi squared pvalue " + str(chi_p)]
            cluster_summary = "\n".join(cluster_summary_list)
            for d in plot_data:
                ax1.plot(d, lw = 0.2)
            ax1.axhline(1.5, lw = 2, c = "r")
            ax1.axhline(2.5, lw = 2, c = "g")
            ax2.plot(median_data, lw = 1, c = "k")
            ax2.axhline(1.5, lw = 2, c = "r")
            ax2.axhline(2.5, lw = 2, c = "g")
            ax3.plot(collapsed_median, lw = 1, c = "k")
            ax3.plot(upper_limit, lw = 2, c = "g")
            ax3.plot(lower_limit, lw = 2, c = "r")
            colors = ["y", "k"]
            for i in range(len(copies)):
                copyname = copies[i]
                loc = locations[i]
                ax1.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 8, c = colors[i%2])
                ax1.text(x = np.mean(loc), y = - 0.5, s = copyname, fontsize = 10,
                         horizontalalignment = "right", rotation = "vertical")
                ax2.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 2, c = colors[i%2])
                ax2.plot(np.arange(loc[0]-1, loc[1]+1), np.zeros(loc[1]- loc[0]+2),
                         lw = 8, c = colors[i%2])
                ax2.text(x = np.mean(loc), y = - 0.5, s = copyname, fontsize = 10,
                         horizontalalignment = "right", rotation = "vertical")
            ax1.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            ax1.set_title("Gene " + gene_name + " cluster " + str(c))
            ax1.set_ylim([0,ymax])
            ax2.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            #ax1.annotate(cluster_summary, xy=(-12, -12), xycoords='axes points',
            #    size=10, ha='right', va='top',bbox=dict(boxstyle='round', fc='w'))
            #ax2.annotate(cluster_summary, xy=(-12, -12), xycoords='axes points',
            #    size=8, ha='right', va='top',bbox=dict(boxstyle='round', fc='w'))
            ax2.set_title("Gene " + gene_name + " cluster " + str(c))
            ax2.set_ylim([0,ymax])
            ax3.legend(cluster_summary_list, loc= "upper right", fontsize=8)
            ax3.set_title("Gene " + gene_name + " cluster " + str(c))
            ax3.set_ylim([0,2*ymax])
        # plot clusters
        if cluster_method == "dbscan":
            ms = cnv_calls[gene_name]["dbscan"]
            Y_plot = cnv_calls[gene_name]["tsne"]["Y"]
            T_plot = cnv_calls[gene_name]["tsne"]["V"]
        else:
            ms = cnv_calls[gene_name]["mean_shift"]
            Y_plot = cnv_calls[gene_name]["tsne"]["Y"]
            T_plot = cnv_calls[gene_name]["tsne"]["T"]
        cluster_labels = ms.labels_
        labels_unique = np.unique(cluster_labels)
        colors = plt.cm.Spectral(np.linspace(0, 1, len(labels_unique)))
        n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
        fig4 = plt.figure()
        if cluster_method != "dbscan":
            ax4 = fig4.add_subplot(2,1,2)
            ax5 = fig4.add_subplot(2,1,1)
            cluster_centers = ms.cluster_centers_
            for k, col in zip(labels_unique, colors):
                if k == -1:
                    # Black used for noise.
                    col = 'k'
                my_members = cluster_labels == k
                cluster_center = cluster_centers[k]
                ax5.plot(T_plot[my_members, 0], T_plot[my_members, 1],
                         'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=5)
                ax5.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor="w",
                         markeredgecolor='k', markersize=14)
                ax5.annotate(str(k), xy = cluster_center, ha="center",
                             va = "center", fontsize = 12)
        else:
            ax4 = fig4.add_subplot(1,1,1)
        for k, col in zip(labels_unique, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
            my_members = cluster_labels == k
            ax4.plot(Y_plot[my_members, 0], Y_plot[my_members, 1],
                     'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=5)
            ax4.plot(np.mean(Y_plot[my_members, 0]), np.mean(Y_plot[my_members, 1]),
                     'o', markerfacecolor="w",
                     markeredgecolor='k', markersize=14)
            ax4.annotate(str(k), xy = (np.mean(Y_plot[my_members, 0]),
                                       np.mean(Y_plot[my_members, 1])),
                         ha="center",
                         va = "center", fontsize = 12)
        if cluster_method == "dbscan":
            ax4.set_title(gene_name + " DBSCAN clusters")
        else:
            ax4.set_title(gene_name + " meanshift clusters")
            ax5.set_title(gene_name + " tSNE clusters")
        # save figures
        fig1.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig2.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig3.tight_layout(pad = 1, w_pad = 1, h_pad = 10)
        fig1.savefig(os.path.join(image_dir, gene_name), dpi = 96)
        fig2.savefig(os.path.join(image_dir, gene_name) + "_medians", dpi = 96)
        fig3.savefig(os.path.join(image_dir, gene_name) + "_collapsed", dpi = 96)
        fig4.savefig(os.path.join(image_dir, gene_name) + "_meanshift", dpi = 96)
        plt.close("all")
    return

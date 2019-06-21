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


def best_mip_set (compatible_mip_sets, compatible_mip_dic, num_para, diff_score_dic, outfile):
    # open file containing compatible mip lists
    with open(os.path.join(primer3_output_dir, compatible_mip_sets), "r") as infile:
        mip_sets = json.load(infile)
    # load dict file that has the compatibility information and mip information
    with open(os.path.join(primer3_output_dir, compatible_mip_dic), "r") as infile:
        mip_dic = json.load(infile)
        # if there is any sets in the list
        if len(mip_sets) > 0:
            best_set = []
            best_score = 0
            # create a number to mip name look up dictionary
            # because mip list has only numbers that correspond to a mips place in dict
            num_lookup = {}
            for mip in list(mip_dic["pair_information"].keys()):
                d = mip_dic["pair_information"][mip]
                place = d["place"]
                num_lookup[place] = mip
            # find which diffs each mip set captures
            for mip_set in mip_sets:
                if mip_set != None:
                    # create a dic for diffs captured cumulatively by all mips in the set
                    merged_caps = {}
                    # create a list for mip scores based on mip sequence and not the captured diffs
                    mip_scores = []
                    # create a variable for mip scores based on the diffs captured
                    diff_scores = 0
                    for m in mip_set:
                        # extract the mip name
                        mip_key = num_lookup[m]
                        # extract the captured diffs from the mip_dic and append to capture list
                        caps = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]
                        # get the diff name (e.g. chr1:1000-1001), use it as key for merged caps
                        # using a dict ensures nonredundancy
                        for diff in caps:
                            merged_caps[diff] = caps[diff]
                        # find out how many paralogs the mip identifies
                        mip_para = ["pair_information"][mip_key]["mip_information"]                                   ["captured_paralog_number"]
                        # extract mip score and append to mip_scores
                        ms = mip_dic["pair_information"][mip_key]                                    ["mip_information"]["mip_score"]
                        # add a bonus for mips identifying multiple paralogs
                        cap_bonus = mip_para * cap_coefficient
                        ms += cap_bonus
                        # add total mip score to mip scores
                        mip_scores.append(ms)
                    # to get more evenly (and well) scoring mips rather than having
                    # the top scoring set, create another bonus score
                    bonus_score = 0
                    # get total number of mips in the set
                    mip_count = len(mip_scores)
                    # define a counter for mips scoring above a certain score
                    ms_count = 0
                    for ms in mip_scores:
                        if ms >= lowest_mip_score:
                            ms_count += 1
                    # ideally we want mip_count == ms_count for all mips to be good
                    # but this does not work well for large sets because there is a good
                    # chance one or two will not score well and bonus will always be zero.
                    if (mip_count - ms_count) <= max_poor_mips:
                        bonus_score = uniformity_bonus
                    # create a dict that has how many paralogs are captured by the mip set
                    paralog_caps = {}
                    # populate dict with paralog numbers and an initial value of 0
                    for i in range(num_para):
                        paralog_caps[i] = 0
                    for d in merged_caps:
                        # find out which paralogs the diff uniquely identifies
                        # this information is only relevant for diffs from pdiffs file
                        # and exonic diffs and filtered diffs are the only ones from pdiffs in targets
                        if (merged_caps[d]["source"] == "exonic_diffs") or                            (merged_caps[d]["source"] == "filtered_diffs"):
                            # the diff is in the form a:1067:CTA:0,1,2 and we want to
                            # analyze the last part of it for paralogs identified
                            difference = (merged_caps[d]["diff"].split(":")[-1]).split(",")
                            # add the paralogs identified to mip_set paralog caps
                            for j in difference:
                                paralog_caps[int(j)] += 1
                        # extract the score of the type of difference and increment diff_scores
                        diff_scores += diff_score_dic[merged_caps[d]["source"]]
                    # calculate how many paralogs identified by the set
                    cap_para = 0
                    for k in paralog_caps:
                        if paralog_caps[k] > 0:
                            cap_para += 1
                    # calculate total score of mip set
                    total_score = bonus_score + sum(mip_scores) + (cap_para**2) * paralog_score +                                   diff_scores
                    if total_score > best_score:
                        best_score = total_score
                        best_set = mip_set
            # print the gene name and best set of mips together with total score
            #print gene[0]
            #print best_set, best_score
            temp_dic = {}
            # print scores and diffs captured by each mip in the best set
            for mip in best_set:
                mip_key = num_lookup[mip]
                temp_dic[mip_key] = mip_dic["pair_information"][mip_key]
                #print mip_dic["pair_information"][mip_key]["mip_information"]["mip_score"]
                for diff in mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]:
                    di = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"][diff]
                    #print diff, di["source"], di["diff"]

            with open(os.path.join(primer3_output_DIR, best_mip_sets), "w")  as outfile:
                json.dump(temp_dic, outfile, indent=1)
        # if there are no compatible mipsets, then find the best mip
        elif (len(mip_sets) == 0) and (len(list(mip_dic["pair_information"].keys())) > 0):
            best_score = 0
            best_mip = ""
            for mip_key in list(mip_dic["pair_information"].keys()):
                # extract the captured diffs from the mip_dic and append to capture list
                caps = mip_dic["pair_information"][mip_key]["mip_information"]["captured_diffs"]
                # score diffs captured
                diff_scores = 0
                # find diff scores from diff source using diff score dict
                for diff in caps:
                    # extract the score of the type of difference and increment diff_scores
                    diff_scores += diff_score_dic[caps[diff]["source"]]
                # extract mip score
                mip_scores = mip_dic["pair_information"][mip_key]["mip_information"]["mip_score"]
                # find how many paralogs uniquely identified by the mip
                cap_para = mip_dic["pair_information"][mip_key]["mip_information"]                                   ["captured_paralog_number"]
                cap_bonus = cap_para * cap_coefficient
                # calculate total score of mip set
                total_score = mip_scores + (cap_para**2) * paralog_score + diff_scores + cap_bonus
                if total_score > best_score:
                    best_score = total_score
                    best_mip = mip_key
            temp_dic = {}
            temp_dic[best_mip] = mip_dic["pair_information"][best_mip]
            #print gene[0]
            #print [mip_dic["pair_information"][best_mip]["place"]]
            #print best_score
            #print mip_dic["pair_information"][best_mip]["mip_information"]["mip_score"]
            dics = mip_dic["pair_information"][best_mip]["mip_information"]["captured_diffs"]
            #for diff in dics:
                #print dics[diff]["source"], dics[diff]["diff"]
            with open(os.path.join(primer3_output_dir, scored + "_best_set"), "w")  as outfile:
                json.dump(temp_dic, outfile, indent=1)
        else:
            print(("No mips available for target region ", gene[0]))
        return


def compatibility (scored_mips, primer3_output_DIR = "", primer_out = "",
                      overlap_same = 0, overlap_opposite = 0):
    # create in/compatibility lists for each mip
    for k in list(scored_mips["pair_information"].keys()):
        # get coordinates of mip arms
        d = scored_mips["pair_information"][k]
        es = ext_start = d["extension_primer_information"]["GENOMIC_START"]
        ee = ext_end = d["extension_primer_information"]["GENOMIC_END"]
        ls = lig_start = d["ligation_primer_information"]["GENOMIC_START"]
        le = lig_end = d["ligation_primer_information"]["GENOMIC_END"]
        # get mip orientation
        ori = d["orientation"]
        #print k, es, ee, ls, le, ori
        # create an incompatibility list
        incompatible = []
        compatible = []
        for mip in list(scored_mips["pair_information"].keys()):
            m = scored_mips["pair_information"][mip]
            nes = next_ext_start = m["extension_primer_information"]            ["GENOMIC_START"]
            nee = next_ext_end = m["extension_primer_information"]            ["GENOMIC_END"]
            nls = next_lig_start = m["ligation_primer_information"]            ["GENOMIC_START"]
            nle = next_lig_end = m["ligation_primer_information"]            ["GENOMIC_END"]
            # get mip orientation
            next_ori = m["orientation"]
            compat = 0
            next_compat = 0
            # check if the two mips are compatible in terms of
            # orientation and coordinates
            if ori == next_ori == "forward":
                if ((ls < nls) and (ls < nes + overlap_same)) or                  ((ls > nls) and (es + overlap_same> nls)):
                    compat = 1
            elif ori == next_ori == "reverse":
                if ((ls < nls) and (es < nls + overlap_same)) or                   ((ls > nls) and (ls + overlap_same> nes)):
                    compat = 1
            elif (ori == "forward") and (next_ori == "reverse"):
                if (ls < nls + overlap_opposite) or                 (es  + overlap_opposite> nes):
                    compat = 1
                elif (es < nls) and (ee < nls + overlap_opposite) and                     (le  + overlap_opposite> nle) and                     (ls < nee + overlap_opposite):
                    compat = 1
                    next_compat = 1
                elif (es > nls) and (es  + overlap_opposite> nle) and                      (ee < nee + overlap_opposite) and                     (le + overlap_opposite > nes):
                    compat = 1
            elif (ori == "reverse") and (next_ori == "forward"):
                if (ls + overlap_opposite > nls) or                 (es < nes + overlap_opposite):
                    compat = 1
                elif (ls > nes) and (ls + overlap_opposite > nee) and                     (le < nle + overlap_opposite) and                     (ee + overlap_opposite>nls):
                    compat = 1
                elif (ls < nes) and (le < nes + overlap_opposite) and                      (ee + overlap_opposite > nee) and                     (es < nle + overlap_opposite):
                    compat = 1
                    next_compat = 1
            if not compat:
                incompatible.append(mip)
            if next_compat:
                compatible.append(mip)
        d["incompatible"] = incompatible
        d["compatible"] = compatible

    return scored_mips

def targets(must_file, diff_list):
    """ Take a file with snps or regions that must be captured by mips and a
    list of other variations of interest (created in ucsc table format) and
    return a target dictionary"""
    print(("Extracting target snps from file " + must_file))
    # create a dict of positions that should be targeted
    targets = {}
    # top priority will be given to a must have list of rs numbers
    # or genomic coordinates. Read from a tab separated file that has
    # id, chr, beg, end (None if not specified). rs numbers are not used at this time.
    # so, genomic coordinates must be provided
    with open(must_file, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                newline = line.strip().split("\t")
                snp_name = newline[0]
                snp_chr = newline[1]
                snp_begin = newline[2]
                snp_end = newline[3]
                # from the coordinates provided, create the dict key
                key = snp_chr + ":" + snp_begin + "-" + snp_end
                # add snp to dictionary
                targets[key] = {"chrom":snp_chr, "begin":snp_begin, "end":snp_end,                                 "name": snp_name, "diff": "na", "source": "must"}
    # add snps from the diff_list
    for diff in diff_list:
        snp_name = diff[4]
        snp_chr = diff[1]
        snp_begin = diff[2]
        snp_end = diff[3]
        src = diff[15]
        # from the coordinates, create the dict key
        key = snp_chr + ":" + snp_begin + "-" + snp_end
        # check if the key is already in the dict
        if not key in list(targets.keys()):
            # add snp to dictionary
            targets[key] = {"chrom":snp_chr, "begin":snp_begin, "end":snp_end,                             "name": snp_name, "diff": "na", "source": src}
    # return targets dictionary
    return targets


def remove_mips(mip_dic):
    """ filter primers so that only top n scoring primers
    ending within the same subregion (determined by bin_size)
    remains."""
    # load extension and ligation primers from file
    shuffled = list(mip_dic.keys())
    random.shuffle(shuffled)
    for m in shuffled:
        try:
            found = 0
            m_cap_start = mip_dic[m].mip["C0"]["capture_start"]
            m_cap_end = mip_dic[m].mip["C0"]["capture_end"]
            m_start = mip_dic[m].mip["C0"]["mip_start"]
            m_end = mip_dic[m].mip["C0"]["mip_end"]
            m_func = mip_dic[m].func_score
            m_tech = mip_dic[m].tech_score
            m_ori = mip_dic[m].mip["C0"]["orientation"]
            m_score = m_func + m_tech
            for n in shuffled:
                try:
                    if mip_dic[m].name != mip_dic[n].name:
                        n_cap_start = mip_dic[n].mip["C0"]["capture_start"]
                        n_cap_end = mip_dic[n].mip["C0"]["capture_end"]
                        n_start = mip_dic[n].mip["C0"]["mip_start"]
                        n_end = mip_dic[n].mip["C0"]["mip_end"]
                        n_func = mip_dic[n].func_score
                        n_tech = mip_dic[n].tech_score
                        n_ori = mip_dic[n].mip["C0"]["orientation"]
                        n_score = n_func + n_tech
                        remove = False
                        if (m_start <= n_start <= m_end) or                             (n_start <= m_start <= n_end):
                            if m_ori == n_ori:
                                remove = True
                            elif ((n_start < m_start) and                                 (n_cap_start < m_start) and                                 (m_cap_start < n_cap_end) and                                 (n_end < m_cap_end)) or                                 ((m_start < n_start) and                                 (m_cap_start < n_start) and                                 (n_cap_start < m_cap_end) and                                 (m_end < n_cap_end)):
                                remove = False
                            else:
                                remove = True
                        if remove:
                            if m_score >= n_score:
                                mip_dic.pop(n)
                            else:
                                mip_dic.pop(m)
                                break
                except KeyError:
                    continue
        except KeyError:
            continue
    return



def add_captures (mip_dic, target_dic):
    for mip in list(mip_dic["pair_information"].keys()):
        d = mip_dic["pair_information"][mip]
        # find which diffs are captured by mip
        # extract mip's coortinates
        ext_start = d["extension_primer_information"]["GENOMIC_START"]
        ext_end = d["extension_primer_information"]["GENOMIC_END"]
        lig_start = d["ligation_primer_information"]["GENOMIC_START"]
        lig_end = d["ligation_primer_information"]["GENOMIC_END"]
        coord = [ext_start, ext_end, lig_start, lig_end].sort()
        mip_start = coord[1]
        mip_end = coord[2]
        # create a dictionary for the targets the mip captures
        captured_snps = {}
        for snp in target_snps:
            if mip_end >= target_snps[snp]["begin"] >= mip_start:
                captured_snps[snp] = target_snps[snp]
        # add captured diffs information to mip dictionary
        d["mip_information"]["captured_diffs"] = captured_snps
    return mip_dic


def add_paralog_info(mip_dic, num_para):
    for mip in list(mip_dic["pair_information"].keys()):
        # extract the captured diffs from the mip_dic
        caps = mip_dic["pair_information"][mip]["mip_information"]["captured_diffs"]
        # create a dic for how many paralogs a mip captures
        mip_caps = {}
        # populate the dic with 0 for each paralog key
        for i in range(num_para):
            mip_caps[i] = 0
        # find the diffs captured that identify a paralog
        for diff in caps:
            # paralog specific mip sources are exonic and filtered diffs only
            if (caps[diff]["source"] == "exonic_diffs") or                (caps[diff]["source"] == "filtered_diffs"):
                # the diff is in the form a:1067:CTA:0,1,2 and we want to
                # analyze the last part of it for paralogs identified
                mip_diff = (caps[diff]["diff"].split(":")[-1]).split(",")
                # add the paralogs identified to mip_set paralog caps
                for j in mip_diff:
                    mip_caps[int(j)] += 1
        # calculate how many paralogs identified by the mip
        mip_para = 0
        # for each paralog
        for k in mip_caps:
            # if at least one of the diffs identifies the paralog
            if mip_caps[k] > 0:
                # increment captured paralog number by 1
                mip_para += 1
        # add this information to mip dic
        mip_dic["pair_information"][mip]["mip_information"]["captured_paralog_number"] = mip_para
    return mip_dic





def filter_data(data_file, filter_field, comparison, criteria, output_file):
    """
    Data fields to filter: sample_name, gene_name, mip_name, barcode_count etc.
    Comparison: for numeric values: gt (greater than), lt, gte (greater than or equal to), lte
    criteria must be numeric as well.
    for string values: in, nin (not in), criteria must be a list to include or exclude.
    for string equality assessments, a list with single value can be used.
    """
    filtered_data = []
    with open(data_file) as infile:
        data = json.load(infile)
    for d in data:
        field_value = d[filter_field]
        if comparison == "gt":
            if field_value > criteria:
                filtered_data.append(d)
        elif comparison == "gte":
            if field_value >= criteria:
                filtered_data.append(d)
        elif comparison == "lt":
            if field_value < criteria:
                filtered_data.append(d)
        elif comparison == "lte":
            if field_value <= criteria:
                filtered_data.append(d)
        elif comparison == "in":
            if field_value in criteria:
                filtered_data.append(d)
        elif comparison == "nin":
            if field_value not in criteria:
                filtered_data.append(d)
    with open(output_file, "w") as outfile:
        json.dump(filtered_data, outfile, indent = 1)



def create_data_table (settings):
    wdir = settings["workingDir"]
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    uniq_file = os.path.join(wdir, settings["uniqueProbeFile"])
    with open(uniq_file) as infile:
        uniq_dict = json.load(infile)
    sample_results_file = os.path.join(wdir, settings["perSampleResults"])
    with open(sample_results_file) as infile:
        counts = json.load(infile)
    sample_names = list(counts.keys())
    big_table = []
    all_tables = {}
    for gene in call_info:
        gene_info = call_info[gene]
        try:
            all_probes = uniq_dict["unique_probes"][gene]
        except KeyError:
            all_probes = {}
        table = []
        for p in all_probes:
            for c in all_probes[p]:
                split_copies = c.split("_")
                if len(split_copies) == 2:
                    c_id = c.split("_")[-1]
                    probe_info = gene_info[p]["copies"][c_id]
                    try:
                        chrom = int(probe_info["chrom"][3:])
                    except ValueError:
                        chrom = 23
                    start = probe_info["capture_start"]
                    end = probe_info["capture_end"]
                    ori = probe_info["orientation"]
                    gc_frac = calculate_gc(probe_info["capture_sequence"])
                elif len(split_copies) > 2:
                    gc_list = []
                    for i in range(0, len(split_copies), 2):
                        c_id = split_copies[i+1]
                        probe_info = gene_info[p]["copies"][c_id]
                        gc_list.append(calculate_gc(probe_info["capture_sequence"]))
                    gc_frac = np.mean(gc_list)
                    start = "na"
                    end = "na"
                    ori = "na"
                    chrom = 99
                s_counts = []
                for s in sample_names:
                    try:
                        bc = counts[s][gene][p][c]["cumulative_data"]                             ["sample_normalized_barcode_count"]
                    except KeyError:
                        bc = 0

                    s_counts.append(bc)
                p_counts = [gene, p, c, chrom, start, end, gc_frac, ori]
                p_counts.extend(s_counts)
                table.append(p_counts)
                big_table.append(p_counts)
        table = sorted(table, key = itemgetter(3, 4, 5))
        all_tables[gene] = table
    big_table = sorted(big_table, key = itemgetter(3, 4, 5))
    result = {"tables": all_tables,
            "big_table" : big_table,
            "sample_names": sample_names}
    tables_file = os.path.join(wdir, settings["tablesFile"])
    with open(tables_file, "w") as outfile:
        json.dump(result, outfile)
    return


def filter_tables (settings):
    wdir = settings["workingDir"]
    tables_file = os.path.join(wdir, settings["tablesFile"])
    with open(tables_file) as infile:
        tables_dic = json.load(infile)
    tables = copy.deepcopy(tables_dic["tables"])
    sample_info_file = os.path.join(wdir, settings["sampleInfoFile"])
    with open(sample_info_file) as infile:
        sample_info = json.load(infile)["samples"]
    min_probe_median = int(settings["minimumProbeMedian"])
    min_sample_median = int(settings["minimumSampleMedian"])
    filtered_tables = {}
    #normalized_tables = {}
    #barcode_tables = {}
    sample_names = copy.deepcopy(tables_dic["sample_names"])
    sample_normalizers = np.asarray([sample_info[s]["sample_normalizer"] for s in sample_names])
    sample_array = np.asarray(sample_names)
    #filtered_samples = sample_array[sample_filter]
    for g in tables:
        t_array = np.asarray(tables[g])
        split_array = np.hsplit(t_array, [8])
        t_info = split_array[0]
        t_data = np.array(split_array[1], float)
        #s_filtered = t_data[:, sample_filter]
        barcode_table = np.transpose(t_data)/sample_normalizers[:, np.newaxis]
        try:
            probe_filter = np.median(t_data, axis = 1) >= min_probe_median
        except IndexError:
            continue
        sample_filter = np.median(barcode_table, axis = 1) >= min_sample_median
        info_filtered = t_info[probe_filter, :]
        data_filtered = t_data[:, sample_filter]
        data_filtered = data_filtered[probe_filter, :]
        median_normalized = 2*(data_filtered / np.median(data_filtered, axis = 1)[:,np.newaxis])
        barcode_filtered = np.transpose(barcode_table)[:, sample_filter]
        barcode_filtered = barcode_filtered[probe_filter, :]
        #filtered_data = np.hstack((info_filtered, data_filtered))
        #normalized_data = np.hstack((info_filtered, median_normalized))
        #barcode_data = np.hstack((info_filtered, barcode_filtered))
        filtered_tables[g] = {"probe_information":info_filtered,
                              "barcode_counts": barcode_filtered,
                              "sample_normalized_barcode_counts": data_filtered,
                              "median_normalized_copy_counts": median_normalized,
                              "sample_names": sample_array[sample_filter]}
        #normalized_tables[g] = normalized_data
        #barcode_tables [g] = barcode_data
    filtered_tables_file = os.path.join(wdir, settings["filteredTablesFile"])
    with open(filtered_tables_file, "wb") as outfile:
        pickle.dump(filtered_tables, outfile)
    return


def generate_clusters(settings):
    cluster_output = {}
    problem_clusters = []
    wdir = settings["workingDir"]
    filtered_tables_file = os.path.join(wdir, settings["filteredTablesFile"])
    with open(filtered_tables_file, 'rb') as infile:
        tables = pickle.load(infile)
    case_file = os.path.join(wdir, settings["caseFile"])
    case_status = {}
    try:
        with open(case_file) as infile:
            for line in infile:
                newline = line.strip().split("\t")
                case_status[newline[0]] = newline[1]
    except IOError:
        pass
    for g in tables:
        try:
            probes = list(tables[g]["probe_information"][:, 1])
            uniq_probes = []
            for p in probes:
                if p not in uniq_probes:
                    uniq_probes.append(p)
            sample_names = tables[g]["sample_names"]
            labels = []
            for s in sample_names:
                try:
                    labels.append(case_status[s])
                except KeyError:
                    labels.append("na")
            copy_counts = tables[g]["median_normalized_copy_counts"]
            barcode_counts = np.transpose(tables[g]["sample_normalized_barcode_counts"])
            collapsed_counts = np.zeros((len(uniq_probes), len(sample_names)))
            for i in range(len(uniq_probes)):
                u = uniq_probes[i]
                for j in range(len(probes)):
                    p = probes[j]
                    if u == p:
                        collapsed_counts[i] += copy_counts[j]
            repeat_counts = []
            for u in uniq_probes:
                repeat_counts.append(probes.count(u))
            upper_limit = np.asarray(repeat_counts)*2 + 0.5
            lower_limit = np.asarray(repeat_counts)*2 - 0.5
            copy_counts = np.transpose(copy_counts)
            collapsed_counts = np.transpose(collapsed_counts)
            ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=100)
            #ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=30)
            #Y  = ts.fit_transform(copy_counts)
            Y  = ts.fit_transform(barcode_counts)
            ts = TSNE(n_components=2, init = "pca", random_state=0, perplexity=50)
            V  = ts.fit_transform(Y)
            U  = ts.fit_transform(V)
            T  = ts.fit_transform(U)
            tsne_keys = {"Y": Y, "V": V, "U": U, "T": T}
            ms = MeanShift(bandwidth = 5, cluster_all = True)
            ms.fit(tsne_keys[settings["tsneKey"]])
            cluster_labels = ms.labels_
            cluster_centers = ms.cluster_centers_
            labels_unique = np.unique(cluster_labels)
            n_clusters_ = len(labels_unique)
            cluster_dict = {"cluster_labels" : cluster_labels,
                            "unique_probes": uniq_probes,
                           "all_probes": probes,
                           "plotting": {"upper_limit": upper_limit,
                                        "lower_limit": lower_limit},
                           "mean_shift": ms,
                           "tsne": {"median": copy_counts,
                                    "collapsed_median": collapsed_counts,
                                    "T": T,
                                    "U": U,
                                    "V": V,
                                    "Y": Y},
                           "clusters": {}}
            for k in range(n_clusters_):
                my_members = cluster_labels == k
                cluster_center = cluster_centers[k]
                cluster_case_count = list(np.asarray(labels)[my_members]).count("case")
                cluster_control_count = list(np.asarray(labels)[my_members]).count("control")
                other_case_count = labels.count("case") - cluster_case_count
                other_control_count = labels.count("control") - cluster_control_count
                cluster_table = [[cluster_case_count, cluster_control_count],
                                 [other_case_count, other_control_count]]
                try:
                    cluster_fish = fisher_exact(cluster_table)
                except:
                    cluster_fish = ["na", "na"]
                try:
                    cluster_chi = chi2_contingency(cluster_table)
                except:
                    cluster_chi = ["na", "na", "na", "na"]
                cluster_dict["clusters"][k] = {"cluster_table": cluster_table,
                           "cluster_stats": {"fisher": {"OR": cluster_fish[0],
                                                        "pvalue": cluster_fish[1]},
                                             "chi": {"pvalue": cluster_chi[1],
                                                     "expected": cluster_chi[3]}},
                           "cluster_data": copy_counts[my_members],
                           "cluster_medians": np.median(copy_counts[my_members], axis=0),
                           "cluster_samples": sample_names[my_members],
                           "cluster_members": my_members,
                           "collapsed_data": collapsed_counts[my_members],
                           "collapsed_medians": np.median(collapsed_counts[my_members], axis=0)
                           }
            cluster_output[g] = cluster_dict
        except Exception:
            problem_clusters.append([g, e])
    cluster_output_file = os.path.join(wdir, settings["clusterOutputFile"])
    with open(cluster_output_file, "wb") as outfile:
        pickle.dump(cluster_output, outfile)
    return problem_clusters


def dbscan_clusters(settings):
    cluster_output = {}
    problem_clusters = []
    wdir = settings["workingDir"]
    cluster_output_file = os.path.join(wdir, settings["clusterOutputFile"])
    with open(cluster_output_file, "rb") as infile:
        cnv_calls = pickle.load(infile)
    filtered_tables_file = os.path.join(wdir, settings["filteredTablesFile"])
    with open(filtered_tables_file, "rb") as infile:
        tables = pickle.load(infile)
    case_file = os.path.join(wdir, settings["caseFile"])
    case_status = {}
    try:
        with open(case_file) as infile:
            for line in infile:
                newline = line.strip().split("\t")
                case_status[newline[0]] = newline[1]
    except IOError:
        pass
    tsne_key = settings["tsneKey"]
    tsne_plot_key = settings["tsnePlotKey"]
    min_samples = int(settings["minClusterSamples"])
    max_unclustered_frac = float(settings["maxUnclusteredFrac"])
    max_cluster_count = int(settings["maxClusterCount"])
    #tsne_keys = {"Y": Y, "V": V, "U": U, "T": T}
    for g in tables:
        try:
            probes = list(tables[g]["probe_information"][:, 1])
            uniq_probes = []
            for p in probes:
                if p not in uniq_probes:
                    uniq_probes.append(p)
            sample_names = tables[g]["sample_names"]
            labels = []
            for s in sample_names:
                try:
                    labels.append(case_status[s])
                except KeyError:
                    labels.append("na")
            copy_counts = tables[g]["median_normalized_copy_counts"]
            collapsed_counts = np.zeros((len(uniq_probes), len(sample_names)))
            for i in range(len(uniq_probes)):
                u = uniq_probes[i]
                for j in range(len(probes)):
                    p = probes[j]
                    if u == p:
                        collapsed_counts[i] += copy_counts[j]
            repeat_counts = []
            for u in uniq_probes:
                repeat_counts.append(probes.count(u))
            upper_limit = np.asarray(repeat_counts)*2 + 0.5
            lower_limit = np.asarray(repeat_counts)*2 - 0.5
            copy_counts = np.transpose(copy_counts)
            collapsed_counts = np.transpose(collapsed_counts)
            Y = cnv_calls[g]["tsne"][tsne_key]
            P = cnv_calls[g]["tsne"][tsne_plot_key]
            cluster_count = []
            eps_range = np.arange(2., 11.)/2
            for eps in eps_range:
                db = DBSCAN( eps = eps, min_samples = min_samples).fit(Y)
                cluster_labels = db.labels_
                unclustered = float(list(cluster_labels).count(-1))/len(cluster_labels)
                # Number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
                cluster_count.append([n_clusters_, unclustered])
            # find the best clustering
            indexes_PF = []
            for i in range(len(cluster_count)):
                cc = cluster_count[i][0]
                un = cluster_count[i][1]
                if cc <= max_cluster_count and un <= max_unclustered_frac :
                    indexes_PF.append(i)
            if len(indexes_PF) == 0:
                problem_clusters.append([g, "no clusters found!"])
                continue
            else:
                if len(indexes_PF) == 1:
                    best_index = indexes_PF[0]
                else:
                    best_index = indexes_PF[0]
                    for i in indexes_PF[1:]:
                        cc = cluster_count[i][0]
                        best_cc = cluster_count[best_index][0]
                        if best_cc == cc:
                            best_index = i
                db = DBSCAN( eps = eps_range[best_index], min_samples = min_samples).fit(Y)
                cluster_labels = db.labels_
                unclustered = float(list(cluster_labels).count(-1))/len(cluster_labels)
                # Number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
            #cluster_labels = ms.labels_
            #cluster_centers = ms.cluster_centers_
            labels_unique = np.unique(cluster_labels)
            cluster_dict = {"cluster_labels" : cluster_labels,
                            "unique_probes": uniq_probes,
                           "all_probes": probes,
                           "plotting": {"upper_limit": upper_limit,
                                        "lower_limit": lower_limit},
                           "mean_shift": cnv_calls[g]["mean_shift"],
                           "dbscan": db,
                           "tsne": cnv_calls[g]["tsne"],
                           "clusters": {}}
            for k in range(n_clusters_):
                my_members = cluster_labels == k
                #cluster_center = cluster_centers[k]
                cluster_case_count = list(np.asarray(labels)[my_members]).count("case")
                cluster_control_count = list(np.asarray(labels)[my_members]).count("control")
                other_case_count = labels.count("case") - cluster_case_count
                other_control_count = labels.count("control") - cluster_control_count
                cluster_table = [[cluster_case_count, cluster_control_count],
                                 [other_case_count, other_control_count]]
                try:
                    cluster_fish = fisher_exact(cluster_table)
                except:
                    cluster_fish = ["na", "na"]
                try:
                    cluster_chi = chi2_contingency(cluster_table)
                except:
                    cluster_chi = ["na", "na", "na", "na"]
                cluster_dict["clusters"][k] = {"cluster_table": cluster_table,
                           "cluster_stats": {"fisher": {"OR": cluster_fish[0],
                                                        "pvalue": cluster_fish[1]},
                                             "chi": {"pvalue": cluster_chi[1],
                                                     "expected": cluster_chi[3]}},
                           "cluster_data": copy_counts[my_members],
                           "cluster_medians": np.median(copy_counts[my_members], axis=0),
                           "cluster_samples": sample_names[my_members],
                           "cluster_members": my_members,
                           "collapsed_data": collapsed_counts[my_members],
                           "collapsed_medians": np.median(collapsed_counts[my_members], axis=0)
                           }
                cluster_output[g] = cluster_dict
        except Exception as e:
            problem_clusters.append([g, e])
    db_output_file = os.path.join(wdir, settings["dbScanOutputFile"])
    with open(db_output_file, "wb") as outfile:
        pickle.dump(cluster_output, outfile)
    return problem_clusters

def get_unique_probes(settings):
    wdir = settings["workingDir"]
    unique_haplotype_file = os.path.join(wdir, settings["haplotypeDictionary"])
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    with open(settings["callInfoDictionary"]) as infile:
        call_info = json.load(infile)
    # keep which copies a mip can be mapped to in copy_names dictionary
    copy_names = {}
    problem_copy_names = []
    gene_names = []
    for m in unique_haplotypes:
        g = m.split("_")[0]
        gene_names.append(g)
        for h in unique_haplotypes[m]:
            if unique_haplotypes[m][h]["mapped"]:
                try:
                    cn = unique_haplotypes[m][h]["copy_name"]
                except KeyError:
                    problem_copy_names.append(h)
                try:
                    copy_names[g][m].append(cn)
                except KeyError:
                    try:
                        copy_names[g][m] = [cn]
                    except KeyError:
                        copy_names[g] = {m: [cn]}
    for g in copy_names:
        for m in copy_names[g]:
            copy_names[g][m] = sorted(list(set(copy_names[g][m])))
    # keep all copy names that have been observed for a gene in all_copies dict
    all_copy_names = {}
    all_copy_names_list = []
    for g in copy_names:
        cn = []
        for m in copy_names[g]:
            cn.extend(copy_names[g][m])
        cn = sorted(list(set(cn)))
        all_copy_names[g] = cn
        all_copy_names_list.extend(cn)

    # keep probe names that always maps to specific paralog copies in unique_keys
    # uniq_keys[g][m1] = ["HBA1_C0", "HBA2_C1"]
    uniq_keys = {}
    # keep probe names that never maps to specific paralog copies in non_uniq_keys
    # non_uniq_keys[g][m2] = ["HBA1_C0_HBA2_C1"]
    non_uniq_keys = {}
    # keep probe names that maps to specific paralog copies for some copies only
    # in semi_uniq_keys [g][m3] = ["HBA1_C0"] (only maps C0 specifically)
    #semi_uniq_keys = {}
    # if a mip does not map to any target copy (always off target or never works)
    no_copy_keys = {}
    for g in all_copy_names:
        copies = all_copy_names[g]
        keys = sorted(list(call_info[g].keys()),
                      key= lambda a: call_info[g][a]["copies"]["C0"]["capture_start"])
        for k in keys:
            try:
                key_copies = copy_names[g][k]
                nucopies = []
                toss_copies = []
                try:
                    for c in key_copies:
                        split_copies = c.split("_")
                        nucopies.append([split_copies[i+1]                                              for i in range(0, len(split_copies), 2)])
                except IndexError:
                    print(split_copies)
                    break
                for i in range(len(nucopies)):
                    query_list = nucopies[i]
                    for j in range(len(nucopies)):
                        target_list = nucopies[j]
                        if i != j:
                            for c in query_list:
                                if c in target_list:
                                    toss_copies.append(i)
                toss_copies = list(set(toss_copies))
                uniq_copies = [key_copies[i]                                for i in range(len(key_copies)) if i not in toss_copies]
                non_uniq_copies = [key_copies[i]                                for i in range(len(key_copies)) if i in toss_copies]
                if len(uniq_copies) > 0:
                    try:
                        uniq_keys[g][k] = uniq_copies
                    except KeyError:
                        uniq_keys[g] = {k: uniq_copies}
                if len(non_uniq_copies) > 0:
                    try:
                        non_uniq_keys[g][k] = non_uniq_copies
                    except KeyError:
                        non_uniq_keys[g] = {k: non_uniq_copies}
            except KeyError:
                try:
                    no_copy_keys[g].append(k)
                except KeyError:
                    no_copy_keys[g] = [k]
    gene_names = sorted(gene_names)
    result =  {"copy_names": copy_names,
            "all_copy_names": all_copy_names,
            "problem_copy_names": problem_copy_names,
            "gene_names": gene_names,
            "unique_probes": uniq_keys,
            "non_unique_probes": non_uniq_keys,
            "no_copy_probes": no_copy_keys}
    uniq_file = os.path.join(wdir, settings["uniqueProbeFile"])
    with open(uniq_file, "w") as outfile:
        json.dump(result, outfile, indent = 1)
    return



def get_data(settings_file):
    settings = get_analysis_settings(settings_file)
    get_haplotypes(settings)
    align_haplotypes(settings)
    parse_aligned_haplotypes(settings)
    update_aligned_haplotypes(settings)
    update_unique_haplotypes(settings)
    update_variation(settings)
    get_raw_data(settings)
    return


def analyze_data(settings_file):
    settings = get_analysis_settings(settings_file)
    group_samples(settings)
    update_raw_data(settings)
    get_counts(settings)
    return


def process_data(wdir, run_ids):
    for rid in run_ids:
        settings_file = os.path.join(wdir, "settings_" + rid)
        get_data(settings_file)
        analyze_data(settings_file)
    return



###############################################################################
# older analysis functions.
###############################################################################


def get_raw_data(settings):
    """ Extract raw data from filtered_data file. If there is data from a
    previous run, new data can be added to old data. If this sample set is new,
    or being analyzed separately, than existing data should be "na".
    Return a list of data point dictionaries. One dict for each
    haplotype/sample. Write this list to disk.
    """
    wdir = settings["workingDir"]
    unique_haplotype_file = os.path.join(wdir, settings["haplotypeDictionary"])
    sequence_to_haplotype_file = os.path.join(
        wdir, settings["sequenceToHaplotypeDictionary"])
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    with open(sequence_to_haplotype_file) as infile:
        sequence_to_haplotype = json.load(infile)
    problem_data = []
    existing_data_file = os.path.join(wdir, settings["existingData"])
    try:
        with open(existing_data_file) as infile:
            raw_data = json.load(infile)
    except IOError:
        raw_data = []
    mipster_file = os.path.join(wdir, settings["mipsterFile"])
    colnames = dict(list(zip(settings["colNames"],
                         settings["givenNames"])))
    with open(mipster_file) as infile:
        # filteredData only contains relevant fields for this analysis
        # the field names are given in the settings dict and
        # kept in given_names list
        run_ID = settings["runID"]
        line_number = 0
        for line in infile:
            newline = line.strip().split("\t")
            if not line_number == 0:
                line_number += 1
                col_indexes = {}
                try:
                    for ck in list(colnames.keys()):
                        col_indexes[
                            newline.index(ck)
                        ] = {"name": colnames[ck]}
                        if ck == 'c_barcodeCnt':
                            bc_index = newline.index(ck)
                        elif ck == 'c_readCnt':
                            rc_index = newline.index(ck)
                except ValueError:
                    for ck in list(colnames.values()):
                        col_indexes[newline.index(ck)] = {"name": ck}
                        if ck == 'barcode_count':
                            bc_index = newline.index(ck)
                        elif ck == 'read_count':
                            rc_index = newline.index(ck)
                # each data point will be a dict
                data_dic = {}
                for i in list(col_indexes.keys()):
                    # check data type and convert data appropriately
                    if i in [bc_index, rc_index]:
                        data_point = int(newline[i])
                    else:
                        data_point = newline[i]
                    data_dic[col_indexes[i]["name"]] = data_point
                data_dic["run_ID"] = run_ID
                # once data dict is created, check the haplotype sequence
                # in the sequence to haplotype dict and update the data dict
                # with values in the haplotype dict
                seq = data_dic.pop("haplotype_sequence")
                try:
                    uniq_id = sequence_to_haplotype[seq]
                    data_dic.pop("haplotype_quality_scores")
                    data_dic["haplotype_ID"] = uniq_id
                    if unique_haplotypes[data_dic["mip_name"]][uniq_id][
                            "mapped"]:
                        raw_data.append(data_dic)
                except KeyError:
                    problem_data.append(data_dic)
                    continue
    # dump the raw_data list to a json file
    raw_data_file = os.path.join(wdir, settings["rawDataFile"])
    with open(raw_data_file, "w") as outfile:
        json.dump(raw_data, outfile)
    problem_data_file = os.path.join(wdir, settings["rawProblemData"])
    with open(problem_data_file, "w") as outfile:
        json.dump(problem_data, outfile, indent=1)
    return


def group_samples (settings):
    wdir = settings["workingDir"]
    raw_data_file = wdir  + settings["rawDataFile"]
    with open(raw_data_file) as infile:
        raw_data = json.load(infile)
    samples = {}
    problem_samples = {}
    mip_names = {}
    #diploid_mipnames = []
    #all_mipnames = []
    copy_stable_genes = settings["copyStableGenes"]
    for c in raw_data:
        gene_name = c["gene_name"]
        sample_name = c["sample_name"]
        mip_name = c["mip_name"]
        if (copy_stable_genes == "na") or (gene_name in copy_stable_genes):
            try:
                samples[sample_name]["diploid_barcode_count"] += c["barcode_count"]
                samples[sample_name]["diploid_read_count"] += c["read_count"]
                samples[sample_name]["diploid_mip_names"].append(mip_name)
                #diploid_mipnames.append(mip_name)
                """
                except KeyError:
                    try:
                        samples[sample_name]["diploid_barcode_count"] = c["barcode_count"]
                        samples[sample_name]["diploid_read_count"] = c["read_count"]
                        samples[sample_name]["diploid_mip_names"] = [mip_name]
                        diploid_mipnames.append(mip_name)
                """
            except KeyError:
                samples[sample_name] = {"diploid_barcode_count": c["barcode_count"],
                                        "diploid_read_count": c["read_count"],
                                        "diploid_mip_names": [mip_name],
                                        "total_barcode_count": c["barcode_count"],
                                        "total_read_count": c["read_count"],
                                        "all_mip_names": [mip_name]
                                        }
                #diploid_mipnames.append(mip_name)
        try:
            samples[sample_name]["total_barcode_count"] += c["barcode_count"]
            samples[sample_name]["total_read_count"] += c["read_count"]
            samples[sample_name]["all_mip_names"].append(mip_name)
            #all_mipnames.append(mip_name)
            """
            except KeyError:
                try:
                    samples[sample_name]["total_barcode_count"] = c["barcode_count"]
                    samples[sample_name]["total_read_count"] = c["read_count"]
                    samples[sample_name]["all_mip_names"] = [mip_name]
                    all_mipnames.append(mip_name)
            """
        except KeyError:
            samples[sample_name] = {"diploid_barcode_count": 0,
                                        "diploid_read_count": 0,
                                        "diploid_mip_names": [],
                                        "total_barcode_count": c["barcode_count"],
                                        "total_read_count": c["read_count"],
                                        "all_mip_names": [mip_name]
                                        }
            #all_mipnames.append(mip_name)


    for s in list(samples.keys()):
        try:
            samples[s]["diploid_mip_names"] = list(set((samples[s]["diploid_mip_names"])))
            samples[s]["diploid_mip_number"] = len(samples[s]["diploid_mip_names"])
            samples[s]["all_mip_names"] = list(set((samples[s]["all_mip_names"])))
            samples[s]["total_mip_number"] = len(samples[s]["all_mip_names"])
            mip_names[s] = {}
            mip_names[s]["diploid_mip_names"] = samples[s].pop("diploid_mip_names")
            mip_names[s]["all_mip_names"] = samples[s].pop("all_mip_names")
            try:
                samples[s]["average_diploid_barcode_count"] = round(
                samples[s]['diploid_barcode_count']
                 /samples[s]['diploid_mip_number'], 2)
                samples[s]["sample_normalizer"] = round(
                100/(samples[s]["average_diploid_barcode_count"]), 2)
            except ZeroDivisionError:
                problem_samples[s] = samples.pop(s)
        except KeyError:
            problem_samples[s] = samples.pop(s)

    results =  {"samples": samples,
            "problem_samples": problem_samples,
            "mip_names": mip_names}
    sample_info_file = os.path.join(wdir, settings["sampleInfoFile"])
    with open(sample_info_file, "w") as outfile:
        json.dump(results,outfile, indent=1)
    return


def update_raw_data (settings):
    wdir = settings["workingDir"]
    raw_data_file = wdir  + settings["rawDataFile"]
    with open(raw_data_file) as infile:
        raw_data = json.load(infile)
    unique_haplotype_file = os.path.join(wdir, settings["haplotypeDictionary"])
    with open(unique_haplotype_file) as infile:
        unique_haplotypes = json.load(infile)
    sample_info_file = os.path.join(wdir, settings["sampleInfoFile"])
    with open(sample_info_file) as infile:
        samples = json.load(infile)["samples"]
    normalized_data = []
    problem_data = []
    for r in raw_data:
        try:
            sample_name = r["sample_name"]
            r["sample_normalizer"] = samples[sample_name]["sample_normalizer"]
            r["sample_normalized_barcode_count"] = r["sample_normalizer"] * r["barcode_count"]
            hid = r["haplotype_ID"]
            mid = r["mip_name"]
            r["copy_name"] = unique_haplotypes[mid][hid]["copy_name"]
            normalized_data.append(r)
        except KeyError:
            problem_data.append(r)
    normalized_data_file = os.path.join(wdir, settings["normalizedDataFile"])
    with open(normalized_data_file, "w") as outfile:
        json.dump(normalized_data, outfile)
    problem_data_file = os.path.join(wdir, settings["normalizedProblemData"])
    with open(problem_data_file, "w") as outfile:
        json.dump(problem_data, outfile, indent=1)
    return


def get_counts (settings):
    wdir = settings["workingDir"]
    data_file = os.path.join(wdir, settings["normalizedDataFile"])
    with open(data_file) as infile:
        data = json.load(infile)
    min_barcode_count = int(settings["minBarcodeCount"])
    min_barcode_fraction = float(settings["minBarcodeFraction"])
    counts = {}
    samples = {}
    # get sample data across probes
    for t in data:
        g = t["gene_name"]
        m = t["mip_name"]
        c = t["copy_name"]
        s = t["sample_name"]
        try:
            samples[s][g][m][c]["raw_data"].append(t)
        except KeyError:
            try:
                samples[s][g][m][c] = {"raw_data": [t]}
            except KeyError:
                try:
                    samples[s][g][m] = {c: {"raw_data": [t]}}
                except KeyError:
                    try:
                        samples[s][g] = {m: {c: {"raw_data": [t]}}}
                    except KeyError:
                        samples[s] = {g: {m: {c: {"raw_data": [t]}}}}
    # filter/merge individual data points for each sample/paralog copy
    # data for the same haplotype, possibly from separate runs will be merged
    # data that does not pass a low barcode threshold will be filtered
    merge_keys = ["sample_normalized_barcode_count",
                  "barcode_count",
                  "read_count"]
    for s in samples:
        for g in samples[s]:
            for m in samples[s][g]:
                for c in samples[s][g][m]:
                    data_points = samples[s][g][m][c]["raw_data"]
                    grouped_data = samples[s][g][m][c]["grouped_data"] = []
                    merged_data = samples[s][g][m][c]["merged_data"] = []
                    filtered_data = samples[s][g][m][c]["filtered_data"] = []
                    cumulative_data = samples[s][g][m][c]["cumulative_data"] = {                            "sample_normalized_barcode_count": 0,
                            "barcode_count": 0,
                            "read_count": 0,
                            "haplotypes": []}
                    temp_data_points = copy.deepcopy(data_points)
                    # group data points with same haplotype_ID together
                    for i in range(len(temp_data_points)):
                        di = temp_data_points[i]
                        if di != "remove":
                            same_haps = [di]
                            for j in range(len(temp_data_points)):
                                if i != j:
                                    dj = temp_data_points[j]
                                    if dj != "remove" and (di["haplotype_ID"] ==                                                            dj["haplotype_ID"]):
                                        same_haps.append(dj)
                                        temp_data_points[j] = "remove"
                            grouped_data.append(same_haps)
                    for dl in grouped_data:
                        # find data point with most barcodes
                        temp_barcode_count = 0
                        best_data_index = 0
                        for i in range(len(dl)):
                            b = dl[i]["barcode_count"]
                            if b > temp_barcode_count:
                                temp_barcode_count = b
                                best_data_index = i
                        # use data point with more barcodes as base
                        dm = copy.deepcopy(dl[best_data_index])
                        for i in range(len(dl)):
                            if i != best_data_index:
                                dt = dl[i]
                                for k in merge_keys:
                                    dm[k] += dt[k]
                        merged_data.append(dm)
                    # filter haplotypes with insufficient count/frequency
                    temp_barcode_counts = []
                    for d in merged_data:
                        temp_barcode_counts.append(d["barcode_count"])
                    frac_threshold = min_barcode_fraction * sum(temp_barcode_counts)
                    for d in merged_data:
                        if d["barcode_count"] >= frac_threshold and                           d["barcode_count"] >= min_barcode_count:
                            filtered_data.append(d)
                    for d in filtered_data:
                        for k in merge_keys:
                            cumulative_data[k] += d[k]
                        cumulative_data["haplotypes"].append(d["haplotype_ID"])
    # get probe information across samples
    template_dict = {"sample_normalized_barcode_count": [],
                     "barcode_count": [],
                     "read_count" : [],
                     "haplotypes": [],
                     "sample_names": []}
    for s in samples:
        for g in samples[s]:
            for m in samples[s][g]:
                for c in samples[s][g][m]:
                    filtered_data = samples[s][g][m][c]["filtered_data"]
                    cumulative_data = samples[s][g][m][c]["cumulative_data"]
                    try:
                        norm_counts = counts[g][m][c]["sample_normalized_barcode_count"]
                    except KeyError:
                        try:
                            counts[g][m][c] = copy.deepcopy(template_dict)
                        except KeyError:
                            try:
                                counts[g][m] = {c: copy.deepcopy(template_dict)}
                            except KeyError:
                                counts[g] = {m: {c: copy.deepcopy(template_dict)}}
                    for k in merge_keys:
                        counts[g][m][c][k].append(cumulative_data[k])
                    counts[g][m][c]["sample_names"].append(s)
                    counts[g][m][c]["haplotypes"].append(filtered_data)
    sample_results_file = os.path.join(wdir, settings["perSampleResults"])
    with open(sample_results_file, "w") as outfile:
        json.dump(samples, outfile, indent=1)
    probe_results_file = os.path.join(wdir, settings["perProbeResults"])
    with open(probe_results_file, "w") as outfile:
        json.dump(counts, outfile, indent=1)
    return



###############################################################################
# END older analysis functions.
###############################################################################

##########################################################
# Possibly deprecated functions
##########################################################


def cnv_stats(hom_case, hom_control,
              wt_case, wt_control):
    """
    Given case/control genotype numbers in the order:
    1) number of homozygous cases,
    2) homozygous controls,
    3) heterozygous cases,
    4) heterozygous controls
    5) wildtype cases
    6) wildtype controls
    Returns a list of length 9:
    1-3) Homozygous mutant vs wildtype
    1) Odds ratio from Fisher's exact test
    2) P value from Fisher's
    3) P value from chi squared test
    4-6) Heterozygous mutant vs wildtype
    7-9) Mutants combined vs witdtype
    Errors return "na" in place of values
    """
    ho_v_wt = [[hom_case, hom_control],
               [wt_case, wt_control]]
    output = []
    tbl = ho_v_wt
    try:
        fish = fisher_exact(tbl)
    except:
        fish = ["na", "na"]
    try:
        chi = chi2_contingency(tbl)
    except:
        chi = ["na", "na", "na", "na"]
    output.extend(fish)
    output.append(chi[1])
    return output


def snp_stats(hom_case, hom_control,
               het_case, het_control,
              wt_case, wt_control):
    """
    Given case/control genotype numbers in the order:
    1) number of homozygous cases,
    2) homozygous controls,
    3) heterozygous cases,
    4) heterozygous controls
    5) wildtype cases
    6) wildtype controls
    Returns a list of length 9:
    1-3) Homozygous mutant vs wildtype
    1) Odds ratio from Fisher's exact test
    2) P value from Fisher's
    3) P value from chi squared test
    4-6) Heterozygous mutant vs wildtype
    7-9) Mutants combined vs witdtype
    Errors return "na" in place of values
    """
    mut_case = het_case + hom_case
    mut_control = het_control + hom_control
    ho_v_wt = [[hom_case, hom_control],
               [wt_case, wt_control]]
    het_v_wt = [[het_case, het_control],
                [wt_case, wt_control]]
    mut_v_wt = [[mut_case, mut_control],
                [wt_case, wt_control]]
    output = []
    for tbl in [ho_v_wt, het_v_wt, mut_v_wt]:
        try:
            fish = fisher_exact(tbl)
        except:
            fish = ["na", "na"]
        try:
            chi = chi2_contingency(tbl)
        except:
            chi = ["na", "na", "na", "na"]
        output.extend(fish)
        output.append(chi[1])
    return output

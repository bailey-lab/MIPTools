"""Contains a suite of functions for copy number variation"""

import copy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import probe_summary_generator
import mip_functions as mip
import allel
import subprocess
plt.style.use('ggplot')
dnacopy = importr("DNAcopy")


def filter_samples(barcode_counts, settings, sample_threshold,
                   probe_threshold):
    """Filter a UMI count table based on per sample and per probe thresholds.

    First, samples failing average UMI threshold are removed. Then, probes
    failing the average UMI count are removed.
    """
    filter_level = settings["copyStableLevel"]
    filter_values = settings["copyStableValues"]
    col = barcode_counts.columns
    all_values = col.get_level_values(filter_level)
    if filter_values == "none":
        filter_values = all_values
    indexer = [c in filter_values for c in all_values]
    cns_counts = barcode_counts.loc[:, indexer]
    cns_medians = cns_counts.median(axis=1)
    sample_mask = cns_medians.loc[
        cns_medians >= sample_threshold].index.tolist()
    probe_medians = barcode_counts.loc[sample_mask].median()
    probe_mask = probe_medians.loc[
        probe_medians >= probe_threshold].index.tolist()
    masked_counts = barcode_counts.loc[sample_mask, probe_mask]
    return masked_counts


def sample_normalize(masked_counts, settings):
    """Normalize a UMI count table sample-wise."""
    filter_level = settings["copyStableLevel"]
    filter_values = settings["copyStableValues"]
    col = masked_counts.columns
    all_values = col.get_level_values(filter_level)
    if filter_values == "none":
        filter_values = all_values
    indexer = [c in filter_values for c in all_values]
    cns_counts = masked_counts.loc[:, indexer]
    cns_totals = cns_counts.sum(axis=1)
    sample_normalized = masked_counts.div(cns_totals, axis=0)
    return sample_normalized


def probe_normalize(sample_normalized, settings):
    """Probe-wise normalize a (sample normalized) count table.

    Probe-wise normalize a (sample normalized) count table assuming
    the average (median or other specified value) probe represents a certain
    copy number. For human samples, for example, assumes the average sample
    will have 2 copies of each target.
    """
    average_copy_count = int(settings["averageCopyCount"])
    norm_percentiles = list(map(float, settings["normalizationPercentiles"]))
    copy_counts = sample_normalized.transform(
        lambda a: average_copy_count * a/(a.quantile(norm_percentiles).mean()))
    return copy_counts


def call_copy_numbers(copy_counts, settings):
    """Call integer copy numbers based on normalized count tables.

    Define breakpoints using DNACopy algorithm.
    """
    copy_calls = {}
    problem_genes = []
    try:
        diploid_upper = float(settings["upperNormalPloidy"])
        diploid_lower = float(settings["lowerNormalPloidy"])
        ploidy = int(settings["ploidy"])
    except KeyError:
        diploid_upper = 2
        diploid_lower = 2
        ploidy = 2
    # Set analysis_level.
    analysis_level = settings["copyNumberCallLevel"]
    # Set cluster method.
    # Options are kmeans and tsne.
    cluster_method = settings["cnvClusterMethod"]
    n_clusters = int(settings["cnvClusterNumber"])
    # segmentation parameters
    # number of probes threshold for large genes
    large_gene_threshold = int(settings["largeGroupTreshold"])
    min_segment_size_large = int(settings["minSegmentSizeLarge"])
    min_segment_size_small = int(settings["minSegmentSizeSmall"])
    sdundo_large = float(settings["sdUndoLarge"])
    sdundo_small = float(settings["sdUndoSmall"])
    # Ecludian distance to re-assign a sample to a larger copy state
    merge_distance = float(settings["reassignDistance"])
    # get analysis units
    if analysis_level != "none":
        analysis_units = set(copy_counts.columns.get_level_values(
            analysis_level))
    else:
        analysis_units = [None]
    # Copy number analysis can and should be performed on data that is
    # separated into logical units.
    # Most intuitively the analysis unit is the gene group,
    # called gene here.
    # For example if we have data on two groups such as
    # Glycophorins and Alphahemoglobins, copy numbers
    # should be called separately on these groups.
    # Sometimes it may be useful to analyze on different units,
    # when target groups contain multiple gene names, for example.
    # The unit will be referred to as gene whether or not it actually
    # is a gene.
    for g in analysis_units:
        # get copy counts for the unit
        try:
            # drop samples where all values are NA
            # drop probes with any NA values
            if g is not None:
                N = copy_counts.xs(
                    g, level=analysis_level, axis=1, drop_level=False).dropna(
                    axis=0, how="all").dropna(axis=1, how="any")
            else:
                N = copy_counts.dropna(axis=0, how="all").dropna(
                    axis=1, how="any")
                g = "all"
            # Cluster samples based on copy counts  accross the analysis unit.
            if cluster_method == "tsne":
                ts = TSNE(init="pca", perplexity=50, random_state=0)
                T = ts.fit_transform(N)
                # apply kmeans clustering to copy counts for this gene
                kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(T)
            else:
                kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(N)
                T = None
            copy_calls[g] = {"copy_counts": N, "tsne": T, "kmeans": kmeans}
            clusters = copy_calls[g]["kmeans_clusters"] = {}
            labels = kmeans.labels_
            # next step is to find break points for each
            # cluster found in kmeans the algorithm removes
            # break points that create segments which are
            # not significantly different, based on SD.
            # Lower number of SDs seem to work better for
            # larger genes. So there will be two different
            # undo_SD vaules depending on the size of the gene,
            # in terms of # of probes.
            if N.shape[1] > large_gene_threshold:
                undo_SD = sdundo_large
                min_segment_size = min_segment_size_large
            else:
                undo_SD = sdundo_small
                min_segment_size = min_segment_size_small
            # find break points for each kmeans cluster center
            for i in range(n_clusters):
                # create an array of the median of copy numbers
                cluster_mask = labels == i
                a_list = np.median(N.loc[cluster_mask], axis=0)
                # Add a negligible value to the array to get rid of
                # zeros (they will lead to an error when taking log)
                a_vec = np.array(a_list) + 1e-20
                # get log2 values of copy numbers
                b_vec = robjects.FloatVector(np.log2(a_vec))
                maploc = N.columns.get_level_values(level="begin")
                res = dnacopy.segment(
                    dnacopy.smooth_CNA(
                        dnacopy.CNA(
                            genomdat=b_vec, maploc=robjects.IntVector(maploc),
                            chrom=robjects.StrVector(
                                ["x" for j in range(len(a_list))])),
                        smooth_region=2),
                    undo_splits="sdundo", undo_SD=undo_SD, min_width=2)
                clusters[i] = res
            # Clean up segmentation
            # All break points found for gene will be merged into
            # a list to be used in determining whether a break point
            # should be used.
            break_points = []
            for clu in clusters:
                # Get segmentation results for the cluster
                # from DNACopy output
                res = clusters[clu]
                seg_sizes = np.array(res[1][4],
                                     dtype=int)
                # Convert log2 copy values to real values
                segment_cns = 2 ** np.array(res[1][5], dtype=float)
                # create a list of segment break points,
                # starting from first probe
                bp_starts = [0]
                for s in seg_sizes:
                    bp_starts.append(bp_starts[-1] + s)
                # extend the cumulative break point list
                # with the bps from this cluster
                if bp_starts[-1] > min_segment_size:
                    for i in range(len(bp_starts) - 1):
                        if seg_sizes[i] < min_segment_size:
                            # if this is the last segment,
                            # merge with the previous
                            if i == len(bp_starts) - 2:
                                bp_starts[i] = "remove"
                            # if this is the first segment,
                            # merge with the next
                            else:
                                if i == 0:
                                    bp_starts[i+1] = "remove"
                                # if this is a middle segment
                                else:
                                    # if copy numbers of flanking
                                    # segments are the same,
                                    # remove the segment entirely
                                    if (int(segment_cns[i-1])
                                            == int(segment_cns[i+1])):
                                        bp_starts[i] = "remove"
                                        bp_starts[i+1] = "remove"
                                    # if copy numbers of flanking
                                    # segments are the same,
                                    # merge segment with the larger
                                    # flanking segment
                                    if seg_sizes[i-1] > seg_sizes[i+1]:
                                        bp_starts[i] = "remove"
                                    else:
                                        bp_starts[i+1] = "remove"
                # recalculate break points and segment sizes
                bp_starts = [b for b in bp_starts if b != "remove"]
                break_points.extend(bp_starts)
            # remove recurrent break points
            break_points = sorted(set(break_points))
            seg_sizes = [break_points[i+1] - break_points[i]
                         for i in range(len(break_points) - 1)]
            # remove break points leading to small segments
            # this will cause smaller segments to be merged with
            # the larger of its flanking segments
            if break_points[-1] > min_segment_size:
                for i in range(len(break_points) - 1):
                    if seg_sizes[i] < min_segment_size:
                        if i == len(break_points) - 2:
                            break_points[i] = "remove"
                        else:
                            if i == 0:
                                break_points[i+1] = "remove"
                            else:
                                if seg_sizes[i-1] > seg_sizes[i+1]:
                                    break_points[i] = "remove"
                                else:
                                    break_points[i+1] = "remove"
            # recalculate break points and segment sizes
            break_points = [b for b in break_points if b != "remove"]
            seg_sizes = [break_points[i+1] - break_points[i]
                         for i in range(len(break_points) - 1)]
            copy_calls[g]["break_points"] = break_points
            copy_calls[g]["segment_sizes"] = seg_sizes
            segment_calls = copy_calls[g]["segment_calls"] = {}
            for i in range(len(break_points) - 1):
                seg_start = break_points[i]
                seg_end = break_points[i + 1]
                # get copy count data for the segment
                seg_counts = N.iloc[:, seg_start:seg_end]
                # Filter noisy probes.
                # If a probe in the segment has much higher
                # std compared to the median std in all probes
                # in the segment, remove the noisy probe
                seg_stds = seg_counts.std()
                med_std = seg_stds.median()
                std_mask = seg_stds <= 2 * med_std
                seg_filtered = seg_counts.loc[:, std_mask]
                seg_rounded = seg_filtered.median(
                    axis=1).round()
                delta_diploid_upper = seg_filtered - diploid_upper
                delta_diploid_lower = seg_filtered - diploid_lower
                delta_rounded = seg_filtered.transform(
                    lambda a: a - seg_rounded)
                distance_to_diploid_lower = delta_diploid_lower.apply(
                    np.linalg.norm, axis=1)
                distance_to_diploid_upper = delta_diploid_upper.apply(
                    np.linalg.norm, axis=1)
                distance_to_diploid_rounded = delta_rounded.apply(
                    np.linalg.norm, axis=1)
                diploid_mask = ((distance_to_diploid_upper
                                <= distance_to_diploid_rounded)
                                | (distance_to_diploid_lower
                                   <= distance_to_diploid_rounded))
                diploid = seg_filtered.loc[diploid_mask].transform(
                    lambda a: a - a + ploidy)
                other_ploid = seg_filtered.loc[~diploid_mask].transform(
                    lambda a: a - a + seg_rounded.loc[~diploid_mask])
                seg_ploid = pd.concat([diploid, other_ploid])
                segment_calls[i] = seg_ploid
            copy_calls[g]["copy_calls"] = pd.concat(segment_calls, axis=1)
            copy_calls[g]["copy_calls"].columns = (
                copy_calls[g]["copy_calls"].columns.droplevel(level=0))
            copy_calls[g]["copy_calls"] = copy_calls[g]["copy_calls"].reindex(
                N.index)
            grouped = pd.DataFrame(copy_calls[g]["copy_calls"].groupby(
                by=copy_calls[g]["copy_calls"].columns.tolist(),
                as_index=False).apply(len))
            grouped = grouped.reset_index().sort_values(by=0, ascending=False)
            grouped.pop(0)
            grouped.columns = pd.MultiIndex.from_tuples(
                grouped.columns,  names=copy_counts.columns.names)
            grouped.index = range(len(grouped))
            copy_calls[g]["collapsed_copy_calls"] = copy.deepcopy(
                copy_calls[g]["copy_calls"])
            for count_i in range(len(N)):
                n_val = N.iloc[count_i]
                for group_i in range(len(grouped)):
                    g_val = grouped.iloc[group_i]
                    c_val = copy_calls[g]["collapsed_copy_calls"].iloc[count_i]
                    g_distance = np.linalg.norm((n_val - g_val).dropna())
                    c_distance = np.linalg.norm((n_val - c_val).dropna())
                    if g_distance - c_distance < merge_distance:
                        copy_calls[g]["collapsed_copy_calls"].iloc[count_i] = (
                            g_val)
                        break
        except Exception as e:
            problem_genes.append([g, e])
    return [copy_calls, problem_genes]


def update_count_columns(count_table):
    """Add genomic coordinates to UMI count table."""
    info_df = probe_summary_generator.get_probe_call_info()
    info_dict = info_df.groupby(["MIP", "Copy"]).first().to_dict(
        orient="index")
    cols = count_table.columns
    new_cols = [c + (info_dict[c]["chrom"], info_dict[c]["capture_start"],
                info_dict[c]["capture_end"], info_dict[c]["Gene"])
                for c in cols]
    new_cols = pd.MultiIndex.from_tuples(
        new_cols, names=["MIP", "Copy", "Chrom", "begin", "end", "Gene"])
    count_table.columns = new_cols
    count_table = count_table.sort_index(level=["Chrom", "begin"], axis=1)
    return


def control_normalize(control_cnv_vcf, sample_normalized_counts, settings):
    """Probe-wise normalize a (sample normalized) count table.

    Probewise normalization based on known copy number of control samples.
    This is specifically designed for use with 1KG project control samples.
    """
    # Create a cnv table from cnv vcf file. This is based on 1KG project's
    # copy number call vcf file. It will not work for any other file.
    ######
    # Extract variants overlapping MIP captures
    with open("regions.txt", "w") as outfile:
        for c in sample_normalized_counts.columns:
            outfile.write("\t".join(map(str, c[2:5])) + "\n")

    res = subprocess.run(["bcftools", "view", "-R", "regions.txt",
                          "-o", "regional.vcf", control_cnv_vcf],
                         stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr)
    kg_vcf = allel.read_vcf("regional.vcf",
                            fields=["CHROM", "POS", "END", "ALT",
                                    "samples", "calldata/*"],
                            alt_number=20)

    alts = []
    for variant in kg_vcf["variants/ALT"]:
        new_alt = [1]
        for v in variant:
            try:
                new_alt.append(int(v[3:-1]))
            except (ValueError, IndexError):
                new_alt.append(np.nan)
        new_alt.append(np.nan)
        alts.append(new_alt)
    alts = np.array(alts)
    copy_counts = []
    for i in range(len(alts)):
        alt_alleles = alts[i]
        sample_genotypes = kg_vcf["calldata/GT"][i]
        variant_copy_counts = []
        for s in sample_genotypes:
            try:
                sample_cn = alt_alleles[s[0]] + alt_alleles[s[1]]
            except IndexError:
                sample_cn = np.nan
            variant_copy_counts.append(sample_cn)
        copy_counts.append(variant_copy_counts)
    cc_array = np.array(copy_counts)

    variant_index = np.stack([kg_vcf["variants/CHROM"],
                              kg_vcf["variants/POS"],
                              kg_vcf["variants/END"]], axis=1)

    vindex = pd.MultiIndex.from_arrays(variant_index.T)

    vcf_df = pd.DataFrame(cc_array, index=vindex, columns=kg_vcf["samples"])
    vcf_df = vcf_df.sort_index()
    control_cnv_table = vcf_df.loc[(pd.IndexSlice[:, :, 0:]), :].T

    id_to_name = {}
    sample_names = set()
    for s in sample_normalized_counts.index:
        sn = "-".join(s.split("-")[:-2])
        sample_names.add(sn)
        id_to_name[s] = sn

    control_sample_names = sample_names.intersection(control_cnv_table.index)
    control_sample_ids = [s for s in sample_normalized_counts.index
                          if id_to_name[s] in control_sample_names]
    control_sample_index = [id_to_name[s] for s in control_sample_ids]
    cont_sn = sample_normalized_counts.loc[control_sample_ids]
    cont_cc = control_cnv_table.loc[control_sample_index]
    cont_cc.index = control_sample_ids

    def get_non_diploid(r):
        if r.max() > 2:
            if r.min() >= 2:
                return r.max()
            else:
                return np.nan
        else:
            return r.min()
    df_list = []
    for c in cont_sn.columns:
        chrom = c[2]
        interval = [c[3], c[4]]
        variant_list = []
        for i in cont_cc.columns:
            if (chrom == i[0]) and (len(mip.overlap(i[1:3], interval)) > 0):
                variant_list.append(cont_cc.loc[:, i])
        if len(variant_list) == 1:
            df_list.extend(variant_list)
        elif len(variant_list) > 1:
            var_df = pd.concat(variant_list, axis=1)
            var_df = var_df.apply(get_non_diploid, axis=1)
            df_list.append(var_df)
        else:
            df_list.append(pd.Series(
                [2 for j in range(len(control_sample_ids))],
                index=control_sample_ids))

    cont_cc = pd.concat(df_list, axis=1)
    cont_cc.columns = sample_normalized_counts.columns
    single_copy = cont_sn / cont_cc
    norm_percentiles = list(map(float, settings["normalizationPercentiles"]))
    normalizer = single_copy.quantile(norm_percentiles).mean()
    control_normalized = sample_normalized_counts / normalizer
    return control_normalized, cont_cc


def simple_copy_caller(gene_name, count_table, probe_threshold,
                       sd_multiplier=4):
    """Make simple copy calls based on UMI counts for a gene."""
    sample_normalized = count_table.div(count_table.mean(axis=1), axis=0)
    geno_counts = count_table.xs(
        gene_name, level="Gene", drop_level=False, axis=1)
    probe_medians = geno_counts.median(axis=1)
    filtered_samples = probe_medians.loc[
        probe_medians >= probe_threshold].index
    normalized = sample_normalized.xs(
        gene_name, level="Gene", drop_level=False, axis=1)
    normalized = normalized.loc[filtered_samples]
    probe_medians = normalized.median()
    normalized = normalized.div(probe_medians, axis=1)
    std_base = normalized.median(axis=1).std()
    normalized = normalized.loc[
        :, normalized.std() < (sd_multiplier * std_base)]
    return {"median_copy_count": normalized.median(axis=1).sort_values(),
            "mean_copy_count": normalized.mean(axis=1).sort_values(),
            "normalized_coverage": normalized}

"""Contains functions required for variant calling and PCA analysis"""

import numpy as np
from sklearn import decomposition
from itertools import cycle
import matplotlib.pyplot as plt


def call_genotypes(count_table, coverage_table, min_count,
                   min_coverage, min_freq):
    """Calls genotypes filtered by tallies and calculates prevalences""" 
    # filter mutation counts for minimum count parameter
    # by setting counts to zero if it is below threshold
    filtered_mutation_counts = count_table.applymap(
        lambda x: 0 if x < min_count else x)
    # filter loci without enough coverage by setting
    # coverage to zero if it is below threshold
    filtered_mutation_coverage = coverage_table.applymap(
        lambda x: 0 if x < min_coverage else x)
    # calculate within sample frequency
    freq = filtered_mutation_counts / filtered_mutation_coverage
    freq.replace(np.inf, np.nan, inplace=True)
    # call genotypes using the minimum within sample
    # allele frequency parameter
    genotypes = freq.applymap(
        lambda x: np.nan if (np.isnan(x) or np.isinf(x))
        else 0 if x == 0
        else 0 if x < min_freq
        else 1 if x < (1 - min_freq)
        else 2).sort_index(axis=1)
    # calculate prevalences
    prev = freq.applymap(
        lambda x: np.nan if (np.isnan(x) or np.isinf(x))
        else 0 if x == 0
        else 0 if x < min_freq
        else 1).sort_index(axis=1)
    return {"genotypes": genotypes, "prevalences": prev, "wsaf": freq,
            "filtered_mutation_counts": filtered_mutation_counts,
            "filtered_mutation_coverage": filtered_mutation_coverage}


def filter_variants(variant_table, sample_drop=None, variant_drop=None,
                    variant_filters=None, impute_func=None, remove_zero=False):
    variant_table = variant_table.sort_index(axis=1)
    """Filters variants, imputes missing values, and prints summary"""
    # filter variants that have 0 frequency across all samples
    if remove_zero:
        allele_sums = variant_table.sum(axis=0)
        allele_filter = allele_sums > 0
        variant_table = variant_table.loc[:, allele_filter]

    # filter data based on specified variant properties
    if variant_filters is not None:
        vf = []
        vl = []
        for lev in variant_filters:
            vf.append(variant_filters[lev])
            vl.append(lev)
        tab = variant_table.xs(vf, level=vl, axis=1, drop_level=False)
    else:
        tab = variant_table
    # filter samples that have less than % of loci called
    if sample_drop is not None:
        table_filt = tab.dropna(axis=0,
                                thresh=tab.shape[1] * (1 - sample_drop))
    else:
        table_filt = tab
    # filter variants that are missing in > % of the samples
    if variant_drop is not None:
        table_filt = table_filt.dropna(
            axis=1, thresh=table_filt.shape[0] * (1 - variant_drop))
    # impute missing values with given function
    if impute_func is not None:
        if impute_func == "mode":
            filled_table = table_filt.fillna(table_filt.mode().iloc[0])
        elif impute_func == "mean":
            filled_table = table_filt.fillna(table_filt.mean())
        elif impute_func == "median":
            filled_table = table_filt.fillna(table_filt.median())
        else:
            print(("{} is not an available imputation function. "
                   "Imputation is not done. Please provide one of "
                   "'mode', 'mean' or 'median'").format(impute_func))
    else:
        filled_table = table_filt
    # print a summary of changing sample and variant numbers with filters
    print(("Started with {} samples and {} variants.\n"
           "{} variants were left after the variants were selected using {}"
           " criteria. \n"
           "After filtering samples and variants for at most {} and {} "
           "missingness, respectively, {} samples and {} variants were left. "
           "Remaining missing values were imputed with the '{}' of the "
           "missing variant across all samples.").format(
        variant_table.shape[0], variant_table.shape[1], tab.shape[1],
        variant_filters, sample_drop, variant_drop, filled_table.shape[0],
        filled_table.shape[1], impute_func))
    return filled_table


def plot_PCA(variant_table, n_components=3, meta_data=None, hue_column=None,
             all_colors=('tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                         'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
                         'tab:olive', 'tab:cyan'), scatter_size=15,
             fontsize=16, fig_size=(4, 12), fig_dpi=150):
    """Perform principal component analysis (PCA) and plots results""" 
    # PCA and plot first 3 components
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(variant_table)
    X = pca.transform(variant_table)
    sample_names = variant_table.index.tolist()
    if (meta_data is not None) and (hue_column is not None):
        sites = meta_data.set_index("Sample ID").loc[sample_names,
                                                     hue_column].values
        unique_sites = set(sites)
        site_color_dict = dict(zip(unique_sites, cycle(all_colors)))
    else:
        sites = np.ones(len(sample_names))
        site_color_dict = {1: "k"}
    fig, axes = plt.subplots(3, 1)
    ax = axes[0]
    for ctry, color in site_color_dict.items():
        ctry_mask = sites == ctry
        if len(X[ctry_mask, 0]) > 0:
            ax.scatter(X[ctry_mask, 0], X[ctry_mask, 1],
                       c=color, label=ctry, s=scatter_size)
    if meta_data is not None:
        ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC1 (%0.1f%%)" % pca.explained_variance_[0],
                  fontsize=fontsize)
    ax.set_ylabel("PC2 (%0.1f%%)" % pca.explained_variance_[1],
                  fontsize=fontsize)
    ax = axes[1]
    for ctry, color in site_color_dict.items():
        ctry_mask = sites == ctry
        if len(X[ctry_mask, 0]) > 0:
            ax.scatter(X[ctry_mask, 1], X[ctry_mask, 2],
                       c=color, label=ctry, s=scatter_size)
    if meta_data is not None:
        ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC2 (%0.1f%%)" % pca.explained_variance_[1],
                  fontsize=fontsize)
    ax.set_ylabel("PC3 (%0.1f%%)" % pca.explained_variance_[2],
                  fontsize=fontsize)
    ax = axes[2]
    for ctry, color in site_color_dict.items():
        ctry_mask = sites == ctry
        if len(X[ctry_mask, 0]) > 0:
            ax.scatter(X[ctry_mask, 0], X[ctry_mask, 2],
                       c=color, label=ctry, s=scatter_size)
    if meta_data is not None:
        ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC1 (%0.1f%%)" % pca.explained_variance_[0],
                  fontsize=fontsize)
    ax.set_ylabel("PC3 (%0.1f%%)" % pca.explained_variance_[2],
                  fontsize=fontsize)
    fig.set_size_inches(*fig_size)
    fig.set_dpi(fig_dpi)
    return pca, X

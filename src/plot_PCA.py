from itertools import cycle
from sklearn import decomposition
import numpy as np
import matplotlib.pyplot as plt


def plot_PCA(var_table, cov_table, min_variant_depth, min_coverage,
             min_within_sample_freq, sample_filter, variant_filter,
             filter_levels=[], filter_values=[],
             meta_data=None, hue_column=None,
             all_colors=('tab:blue', 'tab:orange', 'tab:green', 'tab:red',
                         'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray',
                         'tab:olive', 'tab:cyan')):
    # filter variants
    var_filt = var_table.applymap(lambda a: a if a >= min_variant_depth else 0)
    cov_filt = cov_table.applymap(lambda a: a if a >= min_coverage else 0)
    # calculate within sample frequencies of the variants
    freqs = var_filt/cov_filt
    print(freqs.shape)
    freqs.replace(np.inf, np.nan, inplace=True)
    # filter data based on specified variant properties
    if len(filter_values) > 0:
        freqs = freqs.xs(filter_values, level=filter_levels, axis=1)
    print(freqs.shape)
    # filter variants that have 0 frequency across all samples
    allele_sums = freqs.sum(axis=0)
    allele_filter = allele_sums > 0
    freqs = freqs.loc[:, allele_filter]
    print(freqs.shape)
    # call biallelic genotypes from frequencies,
    # setting a minimum within sample frequency threshold to call
    # an allele present
    geno = freqs.applymap(lambda a: np.nan if np.isnan(a)
                          else 0 if a <= min_within_sample_freq else 1)
    geno = geno.sort_index(axis=1)
    tab = geno.T
    # filter samples that have less than % of loci called
    table_filt = tab.dropna(axis=1, thresh=tab.shape[0] * sample_filter)
    print(table_filt.shape)
    # filter variants that are missing in > % of the samples
    table_filt = table_filt.dropna(
        axis=0, thresh=table_filt.shape[1] * variant_filter)
    print(table_filt.shape)
    filled_table = table_filt.T.fillna(table_filt.T.mode().iloc[0]).T
    # PCA and plot first 3 components
    pca = decomposition.PCA(n_components=3)
    pca.fit(filled_table.T)
    X = pca.transform(filled_table.T)
    sample_names = filled_table.columns.tolist()
    if (meta_data is not None) and (hue_column is not None):
        sites = meta_data.set_index("Sample ID").loc[
            sample_names, hue_column].values
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
                       c=color, label=ctry, s=10)
    ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC1 (%0.1f%%)" % pca.explained_variance_[0], fontsize=16)
    ax.set_ylabel("PC2 (%0.1f%%)" % pca.explained_variance_[1], fontsize=16)
    ax = axes[1]
    for ctry, color in site_color_dict.items():
        ctry_mask = sites == ctry
        if len(X[ctry_mask, 0]) > 0:
            ax.scatter(X[ctry_mask, 1], X[ctry_mask, 2],
                       c=color, label=ctry, s=15)
    ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC2 (%0.1f%%)" % pca.explained_variance_[1], fontsize=16)
    ax.set_ylabel("PC3 (%0.1f%%)" % pca.explained_variance_[2], fontsize=16)
    ax = axes[2]
    for ctry, color in site_color_dict.items():
        ctry_mask = sites == ctry
        if len(X[ctry_mask, 0]) > 0:
            ax.scatter(X[ctry_mask, 0], X[ctry_mask, 2],
                       c=color, label=ctry, s=15)
    ax.legend()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel("PC1 (%0.1f%%)" % pca.explained_variance_[0], fontsize=16)
    ax.set_ylabel("PC3 (%0.1f%%)" % pca.explained_variance_[2], fontsize=16)
    fig.set_size_inches(4, 12)
    fig.set_dpi(150)

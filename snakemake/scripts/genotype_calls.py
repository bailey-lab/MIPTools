mutation_counts=snakemake.input.mutation_counts
mutation_coverage=snakemake.input.mutation_coverage
min_count=snakemake.params.min_count
min_coverage=snakemake.params.min_coverage
min_freq=snakemake.params.min_freq
num_samples_wsaf=snakemake.params.num_samples_wsaf
min_wsaf=snakemake.params.min_wsaf
num_samples_umi=snakemake.params.num_samples_umi
min_umi=snakemake.params.min_umi

filtered_cov=snakemake.output.filtered_cov
filtered_alt=snakemake.output.filtered_alt
wsaf=snakemake.output.wsaf
filtered_gen=snakemake.output.filtered_gen
prevalences_inp=snakemake.output.prevalences_inp
final_filtered_gen=snakemake.output.final_filtered_gen
filtered_prev_inp=snakemake.output.filtered_prev_inp

import sys
sys.path.append("/opt/src")
import PCA
import pandas as pd
import os
wdir = "/opt/analysis/"

mutation_counts=pd.read_csv(mutation_counts, header=list(range(6)), index_col=0)
mutation_coverage=pd.read_csv(mutation_coverage, header=list(range(6)), index_col=0)

print('counts are', mutation_counts, 'cov is', mutation_coverage, 'min count is', min_count, 'min_cov is', min_coverage, 'min freq is', min_freq)

gt_calls = PCA.call_genotypes(mutation_counts, mutation_coverage,
                              min_count, min_coverage, min_freq)

filtered_mutation_counts = gt_calls["filtered_mutation_counts"]
filtered_mutation_counts.to_csv(os.path.join(wdir, filtered_alt))

filtered_mutation_coverage = gt_calls["filtered_mutation_coverage"]
filtered_mutation_coverage.to_csv(os.path.join(wdir, filtered_cov))

freq = gt_calls["wsaf"]
freq.to_csv(os.path.join(wdir, wsaf))

genotypes = gt_calls["genotypes"]
genotypes.to_csv(os.path.join(wdir, filtered_gen))

prevalences = gt_calls["prevalences"]
prevalences.to_csv(os.path.join(wdir, prevalences_inp))

wsaf_filter = ((freq > min_wsaf).sum()) >= num_samples_wsaf
print((f"{wsaf_filter.sum()} of {freq.shape[1]} variants will remain after the "
        "wsaf filter"))

umi_filter = ((filtered_mutation_counts >= min_umi).sum()) > num_samples_umi
print((f"{umi_filter.sum()} of {freq.shape[1]} variants will remain after the "
        "UMI filter"))

targ = freq.columns.get_level_values("Targeted") == "Yes"

variant_mask = targ | (wsaf_filter & umi_filter)
print((f"{variant_mask.sum()} variants will remain in the final call set.\n"
       f"{targ.sum()} variants were targeted and will be kept, and"
       f"{(wsaf_filter & umi_filter).sum()} will be removed by the combined UMI"
       "and WSAF filters."))

filtered_genotypes = genotypes.loc[:, variant_mask]
filtered_genotypes.to_csv(os.path.join(wdir, final_filtered_gen))

filtered_prevalences = prevalences.loc[:, variant_mask]
filtered_prevalences.to_csv(os.path.join(wdir, filtered_prev_inp))
mutation_counts=snakemake.input.mutation_counts
mutation_coverage=snakemake.input.mutation_coverage
min_count=snakemake.params.min_count
min_coverage=snakemake.params.min_coverage
min_freq=snakemake.params.min_freq

filtered_cov=snakemake.output.filtered_cov
filtered_alt=snakemake.output.filtered_alt
wsaf=snakemake.output.wsaf
filtered_gen=snakemake.output.filtered_gen
prevalences_inp=snakemake.output.prevalences_inp

import sys
sys.path.append("/opt/src")
import PCA

gt_calls = PCA.call_genotypes(mutation_counts, mutation_coverage,
                              min_count, min_coverage, min_freq)
gt_calls.keys()

filtered_mutation_counts = gt_calls["filtered_mutation_counts"]
filtered_mutation_counts.to_csv(os.path.join(
        wdir, filtered_alt))
filtered_mutation_counts.head()


filtered_mutation_coverage = gt_calls["filtered_mutation_coverage"]
filtered_mutation_coverage.to_csv(os.path.join(
        wdir, filtered_cov))
filtered_mutation_coverage.head()



freq = gt_calls["wsaf"]
freq.to_csv(os.path.join(
        wdir, wsaf))
freq.head()



genotypes = gt_calls["genotypes"]
genotypes.to_csv(os.path.join(
        wdir, filtered_gen))
genotypes.head()

prevalences = gt_calls["prevalences"]
prevalences.to_csv(os.path.join(wdir, prevalences_inp))
prevalences.head()
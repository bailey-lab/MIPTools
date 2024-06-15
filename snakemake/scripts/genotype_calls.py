mutation_counts=snakemake.input.mutation_counts
mutation_coverage=snakemake.input.mutation_coverage
min_count=snakemake.params.min_count
min_coverage=snakemake.params.min_coverage
min_freq=snakemake.params.min_freq

import sys
sys.path.append("/opt/src")
import PCA

gt_calls = PCA.call_genotypes(mutation_counts, mutation_coverage,
                              min_count, min_coverage, min_freq)
gt_calls.keys()

filtered_mutation_counts = gt_calls["filtered_mutation_counts"]
filtered_mutation_counts.to_csv(os.path.join(
        wdir, "filtered_alternate_AA_table.csv"))
filtered_mutation_counts.head()


filtered_mutation_coverage = gt_calls["filtered_mutation_coverage"]
filtered_mutation_coverage.to_csv(os.path.join(
        wdir, "filtered_coverage_AA_table.csv"))
filtered_mutation_coverage.head()



freq = gt_calls["wsaf"]
freq.to_csv(os.path.join(
        wdir, "within_sample_allele_frequencies.csv"))
freq.head()



genotypes = gt_calls["genotypes"]
genotypes.to_csv(os.path.join(
        wdir, "filtered_genotypes_table.csv"))
genotypes.head()



prevalences = gt_calls["prevalences"]
prevalences.to_csv(os.path.join(wdir, "prevalences_input_table.csv"))
prevalences.head()
import sys
sys.path.append("/opt/src")
import subprocess
import mip_functions as mip
import pandas as pd

wdir='/opt/analysis/'

high_barcode_threshold=snakemake.params.high_barcode_threshold
low_coverage_action=snakemake.params.low_coverage_action
target_coverage_count=snakemake.params.target_coverage_count
target_coverage_fraction=snakemake.params.target_coverage_fraction
target_coverage_key=snakemake.params.target_coverage_key
barcode_coverage_threshold=snakemake.params.barcode_coverage_threshold
barcode_count_threshold=snakemake.params.barcode_count_threshold
assessment_key=snakemake.params.assessment_key
good_coverage_quantile=snakemake.params.good_coverage_quantile
repool_csv=snakemake.params.repool_csv

sample_summary = pd.read_csv(wdir+'sample_summary.csv')
meta = pd.read_csv(wdir+'run_meta.csv')
data_summary = pd.merge(sample_summary, meta)

mip.repool(wdir, data_summary, high_barcode_threshold, 
target_coverage_count=target_coverage_count,
target_coverage_fraction=target_coverage_fraction,
target_coverage_key=target_coverage_key,
barcode_coverage_threshold=barcode_coverage_threshold,
barcode_count_threshold=barcode_count_threshold, 
low_coverage_action=low_coverage_action,
assesment_key=assessment_key,
good_coverage_quantile=good_coverage_quantile,
output_file=repool_csv)

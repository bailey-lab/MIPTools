import os
version = os.environ['VERSION']
import tomllib
import subprocess
from box import Box
from pathlib import Path
with open(f"/opt/user/config.toml", "rb") as f:
    config = Box(tomllib.load(f))

config_threads = config.universal_settings.general_cpu_count
config_bwa_extra = config.check_run_stats.bwa_extra
config_wrangler_info_file = config.variant_calling_inputs.wrangler_info_file_name
config_species = config.check_run_stats.species
config_probe_set = config.check_run_stats.get('probe_set')
config_sample_set = config.check_run_stats.get('sample_set')
config_freebayes_cpu_count = config.prevalence_calling.freebayes_cpu_count
config_min_haplotype_barcodes = config.check_run_stats.min_haplotype_barcodes
config_min_haplotype_samples = config.check_run_stats.min_haplotype_samples
config_min_haplotype_sample_fraction = config.check_run_stats.min_haplotype_sample_fraction

config_high_UMI_theshold = config['check_run_stats']['high_UMI_threshold']
config_low_coverage_action = config['check_run_stats']['low_coverage_action']
config_target_coverage_count = config['check_run_stats'].get('target_coverage_count')
config_target_coverage_fraction = config['check_run_stats']['target_coverage_fraction']
config_target_coverage_key = config['check_run_stats']['target_coverage_key']
config_UMI_coverage_threshold = config['check_run_stats']['UMI_coverage_threshold']
config_UMI_count_threshold = config['check_run_stats']['UMI_count_threshold']
config_assessment_key = config['check_run_stats']['assessment_key']
config_good_coverage_quantile = config['check_run_stats']['good_coverage_quantile']

output_folder = "/opt/user/stats_and_variant_calling"
log_folder = output_folder + "/run_settings"
base_resources = "/opt/resources"
snakemake_directory = "/opt/snakemake"

subprocess.call(f"mkdir -p {log_folder}", shell=True)


rule all:
	input:
		repool_csv=output_folder + "/repool.csv",
		UMI_counts=output_folder + "/UMI_counts.csv",
		output_graph=output_folder + "/umi_heatmap.html",


rule copy_params:
	"""
	copies snakemake file, config file, profile, and python scripts to output
	folder
	"""
	input:
		snakefile=snakemake_directory + "/02_check_run_stats.smk",
		configfile=f"/opt/user/config.toml",
		scripts=snakemake_directory + "/scripts",
	output:
		snakefile=log_folder + "/02_check_run_stats.smk",
		configfile=log_folder + f"/config.toml",
		scripts=directory(log_folder + "/scripts"),
	resources:
		log_dir=log_folder,
	shell:
		"""
		cp {input.snakefile} {output.snakefile}
		cp {input.configfile} {output.configfile}
		cp -r {input.scripts} {output.scripts}
		"""


rule modify_ozkan_settings:
	"""
	copies Ozkan's default settings, plus any user updated settings, to an
	output folder alongside the data for later reference.
	"""
	input:
		snakefile=log_folder + "/02_check_run_stats.smk",
		configfile=log_folder + f"/config.toml"
	params:
		template_settings=base_resources + "/templates/analysis_settings_templates/settings.txt",
		processor_number=config_threads,
		bwa_extra=config_bwa_extra,
		species=config_species,
		probe_set=config_probe_set.strip(),
		freebayes_threads=config_freebayes_cpu_count,
		min_haplotype_barcodes=config_min_haplotype_barcodes,
		min_haplotype_samples=config_min_haplotype_samples,
		min_haplotype_sample_fraction=config_min_haplotype_sample_fraction,
		wdir=output_folder,
	output:
		user_settings=output_folder + "/settings.txt",
	resources:
		log_dir=log_folder,
	script:
		"scripts/modify_ozkan_settings.py"


rule parse_info_file:
	"""
	parses the original info file into multiple sub-files
	"""
	input:
		user_settings=output_folder + "/settings.txt",
	output:
		data=output_folder + "/data.tsv",
		samples=output_folder + "/samples.tsv",
		unique_haplotypes=output_folder + "/unique_haplotypes.csv",
	params:
		wdir=output_folder,
		settings_file="settings.txt",
		info_files=[f"/opt/wrangled_data/{config_wrangler_info_file}"],
		sample_sheets=f"/opt/wrangled_data/sample_sheet.tsv",
		sample_set=config_sample_set.strip(),
		probe_set=config_probe_set.strip(),
	resources:
		log_dir=log_folder,
	script:
		"scripts/parse_info_file.py"


rule map_haplotypes:
	"""
	maps haplotypes against the reference genome and outputs several tables
	showing these mappings and whether they are on target.
	"""
	input:
		data=output_folder + "/data.tsv",
		samples=output_folder + "/samples.tsv",
		unique_haplotypes=output_folder + "/unique_haplotypes.csv",
	params:
		wdir=output_folder,
		settings_file="settings.txt",
	output:
		fastq_haps=output_folder + "/haplotypes.fq",
		haps_sam=output_folder + "/haplotypes_bwa.sam",
		aligned_haps=output_folder + "/aligned_haplotypes.csv",
		all_haps=output_folder + "/all_haplotypes.csv",
		mapped_haps=output_folder + "/mapped_haplotypes.csv",
		offtarget_haps=output_folder + "/offtarget_haplotypes.csv",
		metadata=output_folder + "/run_meta.csv",
		UMI_counts=output_folder + "/UMI_counts.csv",
		haplotype_counts=output_folder + "/haplotype_counts.csv",
		sample_summary=output_folder + "/sample_summary.csv",
	# resources below are currently not utilized - haven't figured out a way to
	# get singularity profile, slurm profile, and high ulimits all at once.
	resources:
		mem_mb=200000,
		time_min=4320,
		nodes=20,
		log_dir=log_folder,
	script:
		"scripts/map_haplotypes.py"


rule graph_UMIs:
	"""
	graphs the UMIs that worked and the UMIs that failed
	"""
	input:
		UMI_counts=output_folder + "/UMI_counts.csv",
		sample_summary_csv=output_folder + "/sample_summary.csv"
	params:
		wdir=output_folder,
	output:
		output_graph=output_folder + "/umi_heatmap.html",
		umi_vs_probe_graph = output_folder + "/umi_count_vs_probe_coverage.html"
	resources:
		log_dir=log_folder,
	script:
		"scripts/graph_UMIs.py"


rule make_repool_table:
	"""
	creates a table that recommends (for each sample) whether it needs to be
	repooled or recaptured based on some user-defined thresholds
	"""
	input:
		output_folder + "/run_meta.csv",
	params:
		high_UMI_threshold=config_high_UMI_theshold,
		low_coverage_action=config_low_coverage_action,
		target_coverage_count=config_target_coverage_count,
		target_coverage_fraction=config_target_coverage_fraction,
		target_coverage_key=config_target_coverage_key,
		UMI_coverage_threshold=config_UMI_coverage_threshold,
		UMI_count_threshold=config_UMI_count_threshold,
		assessment_key=config_assessment_key,
		good_coverage_quantile=config_good_coverage_quantile,
		repool_csv= output_folder + "/repool.csv",
		wdir = output_folder,
	resources:
		log_dir=log_folder,
	output:
		repool_csv=output_folder + "/repool.csv",
	script:
		"scripts/make_repool_table.py"

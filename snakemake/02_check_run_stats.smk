for line in open('/opt/build.sh','r'): 
	if 'mip_version=\"' in line: exec(line)
configfile: f'/opt/config/config_v{mip_version}.yaml'


output_folder = "/opt/user/stats_and_variant_calling"
log_folder = output_folder + "/run_settings"
base_resources = "/opt/resources"
snakemake_directory = "/opt/snakemake"
wrangler_folder = "/opt/user/wrangled_data"
import subprocess

subprocess.call(f"mkdir -p {log_folder}", shell=True)


rule all:
	input:
		snakefile=log_folder + "/02_check_run_stats.smk",
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
		configfile=f"/opt/config/config_v{mip_version}.yaml",
		scripts=snakemake_directory + "/scripts",
	output:
		snakefile=log_folder + "/02_check_run_stats.smk",
		configfile=log_folder + f"/config_v{mip_version}.yaml",
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
	params:
		template_settings=base_resources + "/templates/analysis_settings_templates/settings.txt",
		processor_number=config["general_cpu_count"],
		bwa_extra=config["bwa_extra"],
		species=config["species"],
		probe_set=config["probe_set"].strip(),
		freebayes_threads=config["freebayes_cpu_count"],
		min_haplotype_barcodes=config["min_haplotype_barcodes"],
		min_haplotype_samples=config["min_haplotype_samples"],
		min_haplotype_sample_fraction=config["min_haplotype_sample_fraction"],
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
		info_files=[wrangler_folder + '/' + config["wrangler_file"]],
		sample_sheets=wrangler_folder + "/sample_sheet.tsv",
		sample_set=config["sample_set"].strip(),
		probe_set=config["probe_set"].strip(),
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
		high_UMI_threshold=config["high_UMI_threshold"],
		low_coverage_action=config["low_coverage_action"],
		target_coverage_count=config["target_coverage_count"],
		target_coverage_fraction=config["target_coverage_fraction"],
		target_coverage_key=config["target_coverage_key"],
		UMI_coverage_threshold=config["UMI_coverage_threshold"],
		UMI_count_threshold=config["UMI_count_threshold"],
		assessment_key=config["assessment_key"],
		good_coverage_quantile=config["good_coverage_quantile"],
		repool_csv= output_folder + "/repool.csv",
		wdir = output_folder,
	resources:
		log_dir=log_folder,
	output:
		repool_csv=output_folder + "/repool.csv",
	script:
		"scripts/make_repool_table.py"

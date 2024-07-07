configfile: "/opt/config/config.yaml"


output_folder = "/opt/user/wrangled_data"
fastq_folder = "/opt/fastq_dir"
snakemake_folder = "/opt/snakemake"
project_resources_dir = "/opt/project_resources"


rule all:
	"""
	creates a mip_ids folder and an allMipsSamplesNames.tab.txt file. extracts mips,
	corrects mips, and generates files that can be used to determine sample names as
	well as sample names that had extractable data
	"""
	input:
		setup_finished=output_folder + "/setup_finished.txt",
		output_configfile=output_folder + "/snakemake_params/config.yaml",


rule copy_files:
	input:
		setup_snakefile=snakemake_folder + "/wrangler_by_sample_setup.smk",
		finish_snakefile=snakemake_folder + "/wrangler_by_sample_finish.smk",
		input_configfile="/opt/config/config.yaml",
		in_scripts=snakemake_folder + "/scripts",
	output:
		setup_snakefile=output_folder + "/snakemake_params/setup_run.smk",
		finish_snakefile=output_folder + "/snakemake_params/finish_run.smk",
		output_configfile=output_folder + "/snakemake_params/config.yaml",
		out_scripts=directory(output_folder + "/snakemake_params/scripts"),
	shell:
		"""
		cp {input.setup_snakefile} {output.setup_snakefile}
		cp {input.finish_snakefile} {output.finish_snakefile}
		cp {input.input_configfile} {output.output_configfile}
		cp -r {input.in_scripts} {output.out_scripts}
		"""


rule generate_mip_files:
	"""
	given that I'm repackaging miptools wrangler (so wrangler.sh is not needed)
	and that the existing generate_wrangler_scripts.py seems unnecessarily
	convoluted and that only two files are needed by subsequent steps
	(mipArms.txt and allMipsSamplesNames.tab.txt) I wrote my own
	script for this. Input is an arms file and a sample sheet. Output is an arms
	file with rearranged columns and a two column file with names of all mips
	and names of all samples (with no pairing between columns of any given row).
	"""
	input:
		arms_file=project_resources_dir + "/mip_ids/mip_arms.txt",
		sample_sheet="/opt/input_sample_sheet_directory/"
		+ config["input_sample_sheet"].split("/")[-1],
		fastq_folder=fastq_folder,
	params:
		sample_set=config["sample_set"],
		probe_sets=config["probe_set"],
	output:
		mip_arms=output_folder + "/mip_ids/mipArms.txt",
		sample_file=output_folder + "/mip_ids/allMipsSamplesNames.tab.txt",
		sample_sheet=output_folder + "/sample_sheet.tsv",
	script:
		"scripts/generate_mip_files.py"


rule setup:
	input:
		mip_arms=output_folder + "/mip_ids/mipArms.txt",
		sample_file=output_folder + "/mip_ids/allMipsSamplesNames.tab.txt",
	params:
		output_dir = output_folder + "/analysis",
		project_resources = project_resources_dir,
		fastq_dir=fastq_folder,
	output:
		setup_finished=output_folder + "/setup_finished.txt",
	threads: config["general_cpu_count"]
	shell:
		"""
		MIPWrangler mipSetup \
		  --mipArmsFilename {output_folder}/mip_ids/mipArms.txt \
		  --mipSampleFile {output_folder}/mip_ids/allMipsSamplesNames.tab.txt \
		  --numThreads {threads} \
		  --masterDir {params.output_dir} \
		  --dir {fastq_folder} --mipServerNumber 1
		touch {output.setup_finished}
		"""

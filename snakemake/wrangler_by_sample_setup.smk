import os
from pathlib import Path
from box import Box
version = os.environ['VERSION']

import tomllib
with open(f"/opt/user/config.toml", "rb") as f:
    config = Box(tomllib.load(f))

sample_set = config.wrangler_settings.sample_set
probe_set = config.wrangler_settings.probe_set
thread_count = config.wrangler_settings.cpu_count
output_directory_title = config.wrangler_settings.output_title
input_sample_sheet = Path("/opt") / Path(config.wrangler_inputs.input_sample_sheet).name
output_folder = f"/opt/user/{output_directory_title}"
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
		output_configfile=output_folder + f"/snakemake_params/config.toml",


rule copy_files:
	input:
		input_configfile=f"/opt/user/config.toml",
	output:
		output_configfile=os.path.join(output_folder, "snakemake_params", "config.toml"),
	shell:
		"""
		cp {input.input_configfile} {output.output_configfile}
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
		sample_sheet=input_sample_sheet,
		fastq_folder=fastq_folder,
	params:
		sample_set=sample_set,
		probe_sets=probe_set,
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
	threads: thread_count
	shell:
		"""
		rm -rf {params.output_dir}
		MIPWrangler mipSetup \
		  --mipArmsFilename {output_folder}/mip_ids/mipArms.txt \
		  --mipSampleFile {output_folder}/mip_ids/allMipsSamplesNames.tab.txt \
		  --numThreads {threads} \
		  --masterDir {params.output_dir} \
		  --dir {fastq_folder} --mipServerNumber 1
		touch {output.setup_finished}
		"""

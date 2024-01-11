'''
creates a mip_ids folder and an allMipsSamplesNames.tab.txt file. extracts mips,
corrects mips, and generates files that can be used to determine sample names as
well as sample names that had extractable data.
'''

configfile: 'wrangler_by_sample.yaml'
output=config['output_folder']

rule all:
	input:
		setup_finished=output+'/setup_finished.txt',
#		good_samples=output+'/successfully_extracted_samples.txt',
		output_configfile=output+'/snakemake_params/wrangler_by_sample.yaml'

rule copy_files:
	input:
		setup_snakefile='setup_run.smk',
		finish_snakefile='finish_run.smk',
		input_configfile='wrangler_by_sample.yaml',
		in_scripts='scripts'
	output:
		setup_snakefile=output+'/snakemake_params/setup_run.smk',
		finish_snakefile=output+'/snakemake_params/finish_run.smk',
		output_configfile=output+'/snakemake_params/wrangler_by_sample.yaml',
		out_scripts=directory(output+'/snakemake_params/scripts')
	shell:
		'''
		cp {input.setup_snakefile} {output.setup_snakefile}
		cp {input.finish_snakefile} {output.finish_snakefile}
		cp {input.input_configfile} {output.output_configfile}
		cp -r {input.in_scripts} {output.out_scripts}
		'''

rule generate_mip_files:
	'''
	given that I'm repackaging miptools wrangler (so wrangler.sh is not needed)
	and that the existing generate_wrangler_scripts.py seems unnecessarily
	convoluted and that only two files are needed by subsequent steps
	(mipArms.txt and allMipsSamplesNames.tab.txt) I wrote my own
	script for this. Input is an arms file and a sample sheet. Output is an arms
	file with rearranged columns and a two column file with names of all mips
	and names of all samples (with no pairing between columns of any given row).
	'''
	input:
		arms_file=config['project_resources']+'/mip_ids/mip_arms.txt',
		sample_sheet=config['input_sample_sheet'],
		fastq_folder=config['fastq_dir']
	params:
		sample_set=config['sample_set_used'],
		probe_sets=config['probe_sets_used']
	output:
		mip_arms=output+'/mip_ids/mipArms.txt',
		sample_file=output+'/mip_ids/allMipsSamplesNames.tab.txt',
		sample_sheet=output+'/sample_sheet.tsv'
	script:
		'scripts/generate_mip_files.py'

rule setup:
	input:
		mip_arms=output+'/mip_ids/mipArms.txt',
		sample_file=output+'/mip_ids/allMipsSamplesNames.tab.txt'
	params:
		output_dir='/opt/analysis/analysis',
		project_resources=config['project_resources'],
		wrangler_dir=output,
		sif_file=config['miptools_sif'],
		fastq_dir=config['fastq_dir']
	output:
		setup_finished=output+'/setup_finished.txt'
	threads: config['cpu_count']
	shell:
		'''
		singularity exec \
		-B {params.project_resources}:/opt/project_resources \
		-B {params.wrangler_dir}:/opt/analysis \
		-B {params.fastq_dir}:/opt/data \
		{params.sif_file} \
		MIPWrangler mipSetup --mipArmsFilename /opt/analysis/mip_ids/mipArms.txt --mipSampleFile /opt/analysis/mip_ids/allMipsSamplesNames.tab.txt --numThreads {threads} --masterDir {params.output_dir} --dir /opt/data --mipServerNumber 1
		touch {output.setup_finished}
		'''

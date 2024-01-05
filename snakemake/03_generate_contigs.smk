configfile: 'miptools_analysis_no_jupyter.yaml'
#singularity: config['sif_file']
output_folder=config['output_directory']
log_folder=config['output_directory']+'/run_settings/generate_contigs_and_run_freebayes'
import subprocess
subprocess.call(f'mkdir {log_folder}', shell=True)

rule all:
	input:
		freebayes_command_dict=output_folder+'/freebayes_command_dict.yaml',
		snakefile=log_folder+'/03_generate_contigs.smk'

rule copy_params:
	'''
	copies snakemake file, config file, profile, and python scripts to output
	folder
	'''
	input:
		generate_contigs_snakefile='03_generate_contigs.smk',
		run_freebayes_snakefile = '04_run_freebayes.smk',
		configfile='miptools_analysis_no_jupyter.yaml',
		profile='singularity_profile',
		scripts='scripts'
	output:
		generate_contigs_snakefile=log_folder+'/03_generate_contigs.smk',
		run_freebayes_snakefile = log_folder+'/04_run_freebayes.smk',
		configfile=log_folder+'/miptools_analysis_no_jupyter.yaml',
		profile=directory(log_folder+'/singularity_profile'),
		scripts=directory(log_folder+'/scripts')
	resources:
		log_dir=log_folder
	shell:
		'''
		cp {input.generate_contigs_snakefile} {output.generate_contigs_snakefile}
		cp {input.run_freebayes_snakefile} {output.run_freebayes_snakefile}
		cp {input.configfile} {output.configfile}
		cp -r {input.profile} {output.profile}
		cp -r {input.scripts} {output.scripts}
		'''

rule generate_contigs:
	'''
	generates padded bams and padded fastqs and a list of commands to run freebayes
	'''
	input:
		output_folder+'/aligned_haplotypes.csv'
	output:
		#contig_vcfs=directory(output_folder+'/contig_vcfs'),
		padded_bams=directory(output_folder+'/padded_bams'),
		padded_fastqs=directory(output_folder+'/padded_fastqs'),
		freebayes_command_dict=output_folder+'/freebayes_command_dict.yaml'
		#variants_index=output_folder+'/variants.vcf.gz.csi',
		#variants=output_folder+'/variants.vcf.gz',
		#unfixed_variants=output_folder+'/unfixed.vcf.gz',
		#new_header=output_folder+'/new_vcf_header.txt',
		#warnings=output_folder+'/freebayes_warnings.txt',
		#errors=output_folder+'/freebayes_errors.txt',
		#targets_index=output_folder+'/targets.vcf.gz.tbi',
		#targets_vcf=output_folder+'/targets.vcf.gz'
	params:
		targets_file=config['target_aa_annotation'],
		freebayes_settings=config['freebayes_settings'],
		wdir='/opt/analysis',
		settings_file='settings.txt'
	#resources below are currently not utilized - haven't figured out a way to
	#get singularity profile, slurm profile, and high ulimits all at once.
	resources:
		mem_mb=200000,
		nodes=16,
		time_min=5760,
		log_dir=log_folder
	singularity: config['sif_file']
	script:
		'scripts/generate_contigs.py'

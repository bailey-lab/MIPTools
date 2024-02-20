configfile: 'variant_calling.yaml'
output_folder='/opt/analysis'
log_folder='/opt/analysis/run_settings'
import subprocess
subprocess.call(f'mkdir {log_folder}', shell=True)

if config['target_aa_annotation']:
	targeting=config['target_aa_annotation']
elif config['target_nt_annotation']:
	targeting=config['target_nt_annotation']
elif config['target_aa_annotation'] and config['target_nt_annotation:
	print("can't set both target_aa_annotation and target_nt_annotation, one of"
	"these needs to be false")
	exit()
else:
	targeting=None
	
	
else:
	print("

targets_file_choice=config['target_aa_annotation']

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
		generate_contigs_snakefile='/opt/snakemake/03_generate_contigs.smk',
		run_freebayes_snakefile = '/opt/snakemake/04_run_freebayes.smk',
		scripts='/opt/snakemake/scripts'
	output:
		generate_contigs_snakefile=log_folder+'/03_generate_contigs.smk',
		run_freebayes_snakefile = log_folder+'/04_run_freebayes.smk',
		scripts=directory(log_folder+'/scripts')
	resources:
		log_dir=log_folder
	shell:
		'''
		cp {input.generate_contigs_snakefile} {output.generate_contigs_snakefile}
		cp {input.run_freebayes_snakefile} {output.run_freebayes_snakefile}
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
		targets_file=targeting,
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
	script:
		'scripts/generate_contigs.py'

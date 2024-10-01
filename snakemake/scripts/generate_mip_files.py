import pandas as pd
import subprocess
import os
import time

arms_file=snakemake.input.arms_file
input_fastq_folder=snakemake.input.fastq_folder
input_sample_sheet=snakemake.input.sample_sheet
desired_sample_set=snakemake.params.sample_set.replace(' ','').split(',')
desired_sample_set=set([sample_set.upper() for sample_set in desired_sample_set])
desired_probe_sets=snakemake.params.probe_sets.replace(' ', '').strip().split(',')
mip_arms=snakemake.output.mip_arms
sample_file=open(snakemake.output.sample_file, 'w')
output_sample_sheet=snakemake.output.sample_sheet
subprocess.call(f'cp {input_sample_sheet} {output_sample_sheet}', shell=True)

#grab only selected columns from original arms file and output them to new arms file
arms_df=pd.read_table(arms_file)
arms_df=arms_df[['mip_id', 'mip_family', 'extension_arm', 'ligation_arm', 'extension_barcode_length', 'ligation_barcode_length', 'gene_name', 'mipset']]
arms_df.to_csv(mip_arms, index=False, sep='\t')
sequenced_samples=[sample.split('_')[0] for sample in os.listdir(input_fastq_folder)]

print('sequenced samples are', sequenced_samples)

samples_used=set([])
for line_number, line in enumerate(open(input_sample_sheet)):
	line=line.strip().split('\t')
	if line_number==0:
		sample_name_c, sample_set_c, replicate_c, probe_set_c=line.index('sample_name'), line.index('sample_set'), line.index('replicate'), line.index('probe_set')
	else:
		probe_sets=line[probe_set_c].replace(' ', '').strip().split(',')
		probe_sets=[entry.upper() for entry in probe_sets]
		sample_set=line[sample_set_c].replace(' ', '').strip()
		for desired_probe_set in desired_probe_sets:
			new_sample_name=f'{line[sample_name_c]}-{line[sample_set_c]}-{line[replicate_c]}'
			if new_sample_name in sequenced_samples:
				if desired_probe_set.upper() in probe_sets and sample_set.upper() in desired_sample_set:
					samples_used.add(new_sample_name)

family_df=arms_df[['mip_family']]
family_dict=family_df.to_dict()
family_list=sorted([family_dict['mip_family'][row] for row in family_dict['mip_family']])
sample_list=sorted(list(samples_used))

print(f'number of MIPs is: {len(family_list)}.')
print(f'number of samples is: {len(sample_list)}.')
print('if this is not the number of MIPs you were expecting or not the number '
	'of samples you were expecting, examine your input files carefully and try '
	'again.')
time.sleep(10)

if len(sample_list)<1:
	raise Exception('no samples were retrieved from the sample '
		'sheet. This likely means that the probe_set column of your sample '
		'sheet does not match the probe_set value you provided in the config '
		'file, or that your sample sheet does not match the sample_set value '
		'you provided in the config file. Carefully review both the sample '
		'sheet and config files in a plain text editor and in a spreadsheet; '
		'spreadsheet is good for viewing columns, while text editor is good '
		'for viewing hidden characters like double quotes. Subtle mistakes in '
		'the sample sheet, for example empty columns, double quotes in '
		'probe_set or sample_set columns, and columns with no header, can also '
		'cause this issue')

if len(family_list)<1:
	raise Exception('no MIPs were in the mip_arms.txt file. This '
		'likely means that the mip_arms.txt file located in the mip_ids '
		'subfolder of the project_resources folder you provided in the '
		'config file is either corrupted or does not exist. Create a mip_ids '
		'subfolder inside your project_resources folder and put a correctly '
		'formatted mip_arms.txt file inside it, and try again. This '
		'mip_arms.txt file should only contain the mips found in your mip '
		'panel')


bigger_size=max(len(sample_list), len(family_list))
sample_list=sample_list+['']*(bigger_size-len(sample_list))
family_list=family_list+['']*(bigger_size-len(family_list))

sample_file.write('mips\tsamples\n')
for entry in range(bigger_size):
	sample_file.write(f'{family_list[entry]}\t{sample_list[entry]}\n')

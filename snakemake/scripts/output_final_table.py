'''
still needs to be written
'''
import os
import gzip
#from natsort import natsorted
all_targets=snakemake.params.all_targets
prefix=snakemake.params.prefix
suffix=snakemake.params.suffix
final_table=gzip.open(snakemake.output.final_table, mode='wt')

full_list=[]
header=''
for target in all_targets:
	file_path=prefix+target+suffix
	if os.path.exists(file_path):
		for line_number, line in enumerate(gzip.open(file_path, mode='rt')):
			if line_number>0:
				line=line.strip().split('\t')
#				line[18]='not_shown'
				full_list.append(line)
			elif len(header)==0:
				header=line
#sorted_list=natsorted(full_list)
full_list.sort()

final_table.write(header)
for line in full_list:
	final_table.write('\t'.join(line)+'\n')

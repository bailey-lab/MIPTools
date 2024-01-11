sample_file=snakemake.input.sample_file
good_samples=snakemake.output.good_samples
nested_output=snakemake.params.nested_output

samples=[]
for line in open(sample_file):
	line=line.strip().split('\t')
	if len(line)>1:
		samples.append(line[1])

samples=samples[1:]
summary_dict={}
for sample in samples:
	summary_file=nested_output+f'/analysis/{sample}/{sample}_mipExtraction/extractInfoSummary.txt'
	for line_number, line in enumerate(open(summary_file)):
		if line_number>0:
			good_read_count=int(line.strip().split('\t')[6].split('(')[0])
			summary_dict[sample]=good_read_count

if len(summary_dict)==len(samples):
	#only make output file if all extractions gave data - double checks that all
	#extractions really did complete for all samples
	good_samples=open(good_samples, 'w')
	for sample in summary_dict:
		if summary_dict[sample]>0:
			good_samples.write(sample+'\n')

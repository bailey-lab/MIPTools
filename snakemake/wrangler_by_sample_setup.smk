configfile: 'wrangler_by_sample.yaml'
output=config['output_folder']

all_samples, all_targets=[],[]

for line_number, line in enumerate(open(output+'/mip_ids/allMipsSamplesNames.tab.txt')):
	if line_number>0:
		line=line.rstrip().split('\t')
		if len(line)>1 and len(line[1])>0:
			all_samples.append(line[1])
		if len(line[0])>0:
			all_targets.append(line[0])

final_dict={1: expand(output+'/analysis/{sample}/{sample}_mipExtraction/log.txt', sample=all_samples),
			2: expand(output+'/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt', sample=all_samples),
			3: output+'/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json',
			4: expand(output+'/clustering_status/{sample}_mip_clustering_finished.txt', sample=all_samples),
			5: expand(output+'/analysis/populationClustering/{target}/analysis/log.txt', target=all_targets),
			6: output+'/allInfo.tsv.gz'}
output_choice=config['output_choice']
final_out=final_dict[output_choice]

rule all:
	input:
		final_out

rule extract_by_arm:
	input:
	params:
		output_dir='/opt/analysis/analysis',
		wrangler_dir=output,
		sif_file=config['miptools_sif'],
		fastq_dir=config['fastq_dir']
	resources:
		time_min=240
	output:
		output+'/analysis/{sample}/{sample}_mipExtraction/log.txt'
	shell:
		'''
		singularity exec \
		-B {params.fastq_dir}:/opt/data \
		-B {params.wrangler_dir}:/opt/analysis \
		{params.sif_file} \
		MIPWrangler mipExtractByArm --masterDir {params.output_dir} --sample {wildcards.sample} --overWriteDirs --minCaptureLength=30
		'''
if config['downsample_umi_count']<2**32:
	rule mip_barcode_correction:
		input:
			good_samples=expand(output+'/analysis/{sample}/{sample}_mipExtraction/log.txt', sample=all_samples)
		params:
			output_dir='/opt/analysis/analysis',
			wrangler_dir=output,
			sif_file=config['miptools_sif'],
			downsample_seed=config['downsample_seed'],
			downsample_amount=config['downsample_umi_count']
		resources:
			mem_mb=config['memory_mb_per_step'],
			time_min=20
		output:
			barcode_corrections_finished=output+'/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt'
		shell:
			'''
			singularity exec \
			-B {params.wrangler_dir}:/opt/analysis \
			{params.sif_file} \
			MIPWrangler mipBarcodeCorrection --masterDir {params.output_dir} \
			--downSampleAmount {params.downsample_amount} --downSampleSeed \
			{params.downsample_seed} --overWriteDirs --sample {wildcards.sample}
			'''
else:
	rule mip_barcode_correction:
		input:
			good_samples=expand(output+'/analysis/{sample}/{sample}_mipExtraction/log.txt', sample=all_samples)
		params:
			output_dir='/opt/analysis/analysis',
			wrangler_dir=output,
			sif_file=config['miptools_sif'],
			downsample_seed=config['downsample_seed'],
		resources:
			mem_mb=config['memory_mb_per_step'],
			time_min=20
		output:
			barcode_corrections_finished=output+'/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt'
		shell:
			'''
			singularity exec \
			-B {params.wrangler_dir}:/opt/analysis \
			{params.sif_file} \
			MIPWrangler mipBarcodeCorrection --masterDir {params.output_dir} \
			--doNotDownSample --downSampleSeed \
			{params.downsample_seed} --overWriteDirs --sample {wildcards.sample}
			'''


rule correct_for_same_barcode_contam:
	input:
		all_corrected=expand(output+'/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt', sample=all_samples)
	params:
		output_dir='/opt/analysis/analysis',
		wrangler_dir=output,
		sif_file=config['miptools_sif'],
	resources:
		mem_mb=40000,
		time_min=1440,
		nodes=20
	threads: 20
	output:
		#name is controlled by --logFile
		corrected_barcode_marker=output+'/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json'
	shell:
		'''
		singularity exec \
		-B {params.wrangler_dir}:/opt/analysis \
		{params.sif_file} \
		MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir {params.output_dir} --numThreads {threads} --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1
		'''

rule mip_clustering:
	input:
		corrected_barcode_marker=output+'/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json',
		#sample_dir=output+'/analysis/{sample}'
	params:
		output_dir='/opt/analysis/analysis',
		wrangler_dir=output,
		sif_file=config['miptools_sif']
	resources:
		mem_mb=config['memory_mb_per_step'],
		time_min=60,
	output:
		mip_clustering=output+'/clustering_status/{sample}_mip_clustering_finished.txt'
	shell:
		'''
		singularity exec \
		-B {params.wrangler_dir}:/opt/analysis \
		{params.sif_file} \
		MIPWrangler mipClustering --masterDir {params.output_dir} --overWriteDirs --par /opt/resources/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps --sample {wildcards.sample}
		touch {output.mip_clustering}
		'''

rule pop_cluster_target:
	input:
		mip_cluster_files=expand(output+'/clustering_status/{sample}_mip_clustering_finished.txt', sample=all_samples)
	params:
		output_dir='/opt/analysis/analysis',
		wrangler_dir=output,
		sif_file=config['miptools_sif']
	resources:
		mem_mb=config['memory_mb_per_step'],
		time_min=60,
	output:
		pop_clustering=output+'/analysis/populationClustering/{target}/analysis/log.txt'
	shell:
		'''
		singularity exec \
		-B {params.wrangler_dir}:/opt/analysis \
		{params.sif_file} \
		MIPWrangler mipPopulationClustering --keepIntermediateFiles --masterDir {params.output_dir} --overWriteDirs --cutoff 0 --countEndGaps --fraccutoff 0.005 --mipName {wildcards.target}
		touch {output.pop_clustering}
		'''

rule output_final_table:
	'''
	cat together output files of previous step into a final file, do a "natural
	sort" to sort things similar to how Nick's are output. gzip it
	'''
	input:
		pop_clustering=expand(output+'/analysis/populationClustering/{target}/analysis/log.txt', target=all_targets)
#		final_sample_outputs=expand('/path/to/sample/outputs/{sample}.something', sample=sample_list)
	params:
		all_targets=all_targets,
		prefix=output+'/analysis/populationClustering/',
		suffix='/analysis/selectedClustersInfo.tab.txt.gz'
	resources:
		mem_mb=20000,
		time_min=480		
	output:
		final_table=output+'/allInfo.tsv.gz'
	script:
		'scripts/output_final_table.py'

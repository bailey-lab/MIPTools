for line in open('/opt/build.sh','r'): 
	if 'mip_version=\"' in line: exec(line)
configfile: f'/opt/config/config_v{mip_version}.yaml'


output_folder = "/opt/user/wrangled_data"
base_resources = "/opt/resources"

all_samples, all_targets = [], []

for line_number, line in enumerate(
	open(output_folder + "/mip_ids/allMipsSamplesNames.tab.txt")
):
	if line_number > 0:
		line = line.rstrip().split("\t")
		if len(line) > 1 and len(line[1]) > 0:
			all_samples.append(line[1])
		if len(line[0]) > 0:
			all_targets.append(line[0])

final_dict = {
	1: expand(
		output_folder + "/analysis/{sample}/{sample}_mipExtraction/log.txt",
		sample=all_samples,
	),
	2: expand(
		output_folder
		+ "/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt",
		sample=all_samples,
	),
	3: output_folder + "/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json",
	4: expand(
		output_folder + "/clustering_status/{sample}_mip_clustering_finished.txt",
		sample=all_samples,
	),
	5: expand(
		output_folder + "/analysis/populationClustering/{target}/analysis/log.txt",
		target=all_targets,
	),
	6: [output_folder+'/saved_panels/'+config['probe_set']+'.json', output_folder + "/extractInfoSummary.tsv.gz"]

}
output_choice = config["output_choice"]
final_out = final_dict[output_choice]


rule all:
	input:
		final_out,


rule extract_by_arm:
	params:
		output_dir = output_folder + "/analysis",
	resources:
		time_min=240,
	output:
		output_folder + "/analysis/{sample}/{sample}_mipExtraction/log.txt",
	shell:
		"""
		MIPWrangler mipExtractByArm --masterDir {params.output_dir} --sample {wildcards.sample} --overWriteDirs --minCaptureLength=30
		"""


if config["downsample_umi_count"] < 2**32:

	rule mip_barcode_correction:
		input:
			good_samples=expand(
				output_folder + "/analysis/{sample}/{sample}_mipExtraction/log.txt",
				sample=all_samples,
			),
		params:
			output_dir=output_folder + "/analysis",
			downsample_seed=config["downsample_seed"],
			downsample_amount=config["downsample_umi_count"],
		resources:
			mem_mb=config["memory_mb_per_step"],
			time_min=20,
		output:
			barcode_corrections_finished=output_folder
			+ "/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt",
		shell:
			"""
			MIPWrangler mipBarcodeCorrection --masterDir {params.output_dir} \
			--downSampleAmount {params.downsample_amount} --downSampleSeed \
			{params.downsample_seed} --overWriteDirs --sample {wildcards.sample}
			"""

else:

	rule mip_barcode_correction:
		input:
			good_samples=expand(
				output_folder + "/analysis/{sample}/{sample}_mipExtraction/log.txt",
				sample=all_samples,
			),
		params:
			output_dir=output_folder + "/analysis",
			downsample_seed=config["downsample_seed"],
		resources:
			mem_mb=config["memory_mb_per_step"],
			time_min=20,
		output:
			barcode_corrections_finished=output_folder
			+ "/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt",
		shell:
			"""
			MIPWrangler mipBarcodeCorrection --masterDir {params.output_dir} \
			--doNotDownSample --downSampleSeed \
			{params.downsample_seed} --overWriteDirs --sample {wildcards.sample}
			"""


rule correct_for_same_barcode_contam:
	input:
		all_corrected=expand(
			output_folder
			+ "/analysis/{sample}/{sample}_mipBarcodeCorrection/barcodeFilterStats.tab.txt",
			sample=all_samples,
		),
	params:
		output_dir=output_folder + "/analysis",
	resources:
		mem_mb=40000,
		time_min=1440,
		nodes=20,
	threads: 20
	output:
		#name is controlled by --logFile
		corrected_barcode_marker=output_folder
		+ "/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json",
	shell:
		"""
		MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir {params.output_dir} --numThreads {threads} --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1
		"""


rule mip_clustering:
	input:
		corrected_barcode_marker=output_folder
		+ "/analysis/logs/mipCorrectForContamWithSameBarcodes_run1.json",
	params:
		output_dir=output_folder + "/analysis",
	resources:
		mem_mb=config["memory_mb_per_step"],
		time_min=60,
	output:
		mip_clustering=output_folder
		+ "/clustering_status/{sample}_mip_clustering_finished.txt",
	shell:
		"""
		MIPWrangler mipClustering --masterDir {params.output_dir} --overWriteDirs --par {base_resources}/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps --sample {wildcards.sample}
		touch {output.mip_clustering}
		"""


rule pop_cluster_target:
	input:
		mip_cluster_files=expand(
			output_folder + "/clustering_status/{sample}_mip_clustering_finished.txt",
			sample=all_samples,
		),
	params:
		output_dir=output_folder + "/analysis",
	resources:
		mem_mb=config["memory_mb_per_step"],
		time_min=60,
	output:
		pop_clustering=output_folder
		+ "/analysis/populationClustering/{target}/analysis/log.txt",
	shell:
		"""
		MIPWrangler mipPopulationClustering --keepIntermediateFiles --masterDir {params.output_dir} --overWriteDirs --cutoff 0 --countEndGaps --fraccutoff 0.005 --mipName {wildcards.target}
		touch {output.pop_clustering}
		"""


rule output_final_table:
	"""
	cat together output files of previous step into a final file, do a "natural
	sort" to sort things similar to how Nick's are output. gzip it
	"""
	input:
		pop_clustering=expand(
			output_folder + "/analysis/populationClustering/{target}/analysis/log.txt",
			target=all_targets,
		),
	params:
		all_targets=all_targets,
		prefix=output_folder + "/analysis/populationClustering/",
		suffix="/analysis/selectedClustersInfo.tab.txt.gz",
	resources:
		mem_mb=20000,
		time_min=480,
	output:
		final_table=output_folder + "/allInfo.tsv.gz",
	script:
		"scripts/output_final_table.py"

rule concatenate_summary_files:
	input:
		rules.output_final_table.output,
	output:
		extract_info_summary = output_folder + "/extractInfoSummary.tsv.gz",
		extract_info_by_target = output_folder + "/extractInfoByTarget.tsv.gz",
		stitch_info_by_target = output_folder + "/stitchInfoByTarget.tsv.gz",

	shell:
		r"""
		sed -r '1d;s/(\s+)?\S+//2' {output_folder}/analysis/resources/sampleInputFiles.tab.txt |
			awk '$2=$1' |
			sed "s/ /\//g;s/$/_mipExtraction\/extractInfoSummary.txt/g;s/^/\/opt\/user\/wrangled_data\/analysis\//g" \
			| xargs cat \
			| sed '1!{{/Sample/d}}' \
			| pigz > {output.extract_info_summary}
		
		sed -r '1d;s/(\s+)?\S+//2' {output_folder}/analysis/resources/sampleInputFiles.tab.txt |
			awk '$2=$1' |
			sed "s/ /\//g;s/$/_mipExtraction\/extractInfoByTarget.txt/g;s/^/\/opt\/user\/wrangled_data\/analysis\//g" \
			| xargs cat \
			| sed '1!{{/Sample/d}}' \
			| pigz > {output.extract_info_by_target}
		
		sed -r '1d;s/(\s+)?\S+//2' {output_folder}/analysis/resources/sampleInputFiles.tab.txt |
			awk '$2=$1' |
			sed "s/ /\//g;s/$/_mipExtraction\/stitchInfoByTarget.txt/g;s/^/\/opt\/user\/wrangled_data\/analysis\//g" \
			| xargs cat \
			| sed '1!{{/Sample/d}}' \
			| pigz > {output.stitch_info_by_target}
		"""

rule convert_pmo:
	input:
		final_table=output_folder + "/allInfo.tsv.gz",
		mip_arms=output_folder + "/mip_ids/mipArms.txt"
	params:
		species_ID=config['species_ID'],
		genome_URL=config['genome_URL'],
		gff_URL=config['gff_URL'],
		panel_ID=config['probe_set'],
		genome_name=config['genome_name'],
		genome_version=config['genome_version'],
		wrangled_name=config['run_ID']
	output:
		panel_pmo=output_folder+'/saved_panels/'+config['probe_set']+'.json',
		microhaplotype_pmo=output_folder+'/saved_microhaps/'+config['run_ID']+'.json',
		combined_pmo=output_folder+'/final_PMO/'+config['probe_set']+'_'+config['run_ID']+'.json'
	script:
		'scripts/convert_pmo.py'
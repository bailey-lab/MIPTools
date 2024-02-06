import tqdm_thing

sample_fastq=snakemake.input.sample_fastq
sample_output=snakemake.output.sample_output
sample=snakemake.wildcards.sample


#see /nfs/jbailey5/baileyweb/asimkin/miptools/miptools_by_sample_prototyping/output/analysis/analysis/D10-JJJ-28/D10-JJJ-28_mipExtraction/k13_S0_Sub0_mip3_ref/k13_S0_Sub0_mip3_ref.fastq.gz for example UMIs
#wrangler_downsample_umi(wrangler_downsample_list)

mip_barcode_correction(barcode_correction_list)
#command: MIPWrangler mipBarcodeCorrection --keepIntermediateFiles --masterDir {params.output_dir} --numThreads {threads} --overWriteDirs --sample {sample}
#marker file: /nfs/jbailey5/baileyweb/asimkin/miptools/miptools_by_sample_prototyping/output/analysis/analysis/D10-JJJ-2/D10-JJJ-2_mipBarcodeCorrection/barcodeFilterStats.tab.txt

mip_correct_for_contam_with_same_barcodes(correct_for_contam_list)
#

mip_clustering(mip_clustering)


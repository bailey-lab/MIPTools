corrected_barcode_marker=snakemake.input.corrected_barcode_marker
mip_cluster_finished=snakemake.output.mip_cluster_finished
sample=snakemake.wildcards.sample

MIPWrangler mipClustering --masterDir {params.output_dir} --numThreads {threads} --overWriteDirs --par /opt/resources/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps --sample sample
subprocess.call(f'touch {mip_cluster_finished}', shell=True)


#!/usr/bin/env bash

if [[ $# -ne 5 && $# -ne 6 ]]; then
    msg="Illegal number of parameters. Needs five or six arguments:\n"
    msg="${msg}1) The name of the MIP server number.\n"
    msg="${msg}2) The number of threads to use.\n"
    msg="${msg}3) The population clustering fraction cutoff.\n"
    msg="${msg}4) The threshold for downsampling the UMI count.\n"
    msg="${msg}5) A flag indicating if downsmapling should be weighted.\n"
    msg="${msg}   Either an empty string or the -w flag as a string.\n"
    msg="${msg}6) (Optional) The random seed option for downsampling as a string,\n"
    msg="${msg}   e.g. '--seed 7', or an empty string to use the downsampler's default seed."
    echo ${msg} >&2
    exit 2
fi

# Correct barcodes
MIPWrangler mipBarcodeCorrectionMultiple                --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipBarcodeCorrecting_run1 --allowableErrors 6
MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1

# Downsample UMI counts
find . -type f -path '*mipBarcodeCorrection/*.fastq.gz' -exec python /opt/src/wrangler_downsample_umi.py --cpu-count $2 --downsample-threshold $4 $5 ${6:-} {} +

# Cluster barcodes and MIPs
MIPWrangler mipClusteringMultiple                       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipClustering_run1 --par /opt/resources/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps
MIPWrangler mipPopulationClusteringMultiple             --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipPopClustering_run1 --cutoff 0 --countEndGaps --fraccutoff $3

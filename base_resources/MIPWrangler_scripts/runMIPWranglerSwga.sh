#!/usr/bin/env bash

if [[ $# -ne 3 ]]; then
    msg="Illegal number of parameters. Needs three arguments:\n"
    msg="${msg}1) The name of the MIP server number\n"
    msg="${msg}2) The number of threads to use\n"
    msg="${msg}3) The population clustering fraction cutoff."
    echo ${msg} >&2
    exit 2
fi

# Correct barcodes
MIPWrangler mipBarcodeCorrectionMultiple                --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipBarcodeCorrecting_run1 --allowableErrors 6
MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1

# Cluster barcodes and MIPs
MIPWrangler mipClusteringMultiple                       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipClustering_run1 --par /opt/resources/clustering_pars/illumina_swga.pars.txt --countEndGaps
MIPWrangler mipPopulationClusteringMultiple             --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipPopClustering_run1 --cutoff 0 --countEndGaps

#!/usr/bin/env bash

if [[ $# -ne 4 ]]; then
    msg="Illegal number of parameters. Needs three arguments:\n"
    msg="${msg}1) The name of the MIP server number\n"
    msg="${msg}2) The number of threads to use\n"
    msg="${msg}3) The population clustering fraction cutoff\n"
    msg="${msg}4) The threshold for downsampling the UMI count."
    echo ${msg} >&2
    exit 2
fi

# Correct barcodes
MIPWrangler mipBarcodeCorrectionMultiple                --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipBarcodeCorrecting_run1 --allowableErrors 6
MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1

# Downsample UMI counts
for file in $(find -type f -path '*mipBarcodeCorrection/*.fastq.gz'); do
    # Unzip file
    gzip -df ${file}

    # Find name of unzipped file
    unzipped=$(echo ${file} | sed -e 's/\.[^./]*$//')

    # Replace values greater than the threshold with the threshold
    awk -v threshold="$4" -F 'readCnt=|]' '/^@/ {
        if ( $2 > threshold )
            downsample=threshold;
        else
            downsample=$2;
        sub(/readCnt=[[:digit:]]+/,"readCnt="downsample);
    } 1' ${unzipped} > ${file}.tmp && mv ${file}.tmp ${unzipped}

    # Rezip unzip files
    gzip -f ${unzipped}
done

# Cluster barcodes and MIPs
MIPWrangler mipClusteringMultiple                       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipClustering_run1 --par /opt/resources/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps
MIPWrangler mipPopulationClusteringMultiple             --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipPopClustering_run1 --cutoff 0 --countEndGaps --fraccutoff $3

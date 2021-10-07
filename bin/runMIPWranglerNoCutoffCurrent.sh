#!/usr/bin/env bash
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters. Needs 2 arguments: 1) name of mip server number, 2) num of threads to use." >&2
    exit 2
fi

MIPWrangler mipBarcodeCorrectionMultiple                --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipBarcodeCorrecting_run1 --allowableErrors 6
MIPWrangler mipCorrectForContamWithSameBarcodesMultiple --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipCorrectForContamWithSameBarcodes_run1
MIPWrangler mipClusteringMultiple                       --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipClustering_run1 --par /opt/resources/clustering_pars/illumina_collapseHomoploymers.pars.txt --countEndGaps
MIPWrangler mipPopulationClusteringMultiple             --masterDir $(realpath ./)  --numThreads $2 --overWriteDirs --overWriteLog --logFile mipPopClustering_run1 --cutoff 0 --countEndGaps --fraccutoff 0
#nohup MIPWrangler mav			                        --masterDir $(realpath ./)  --numThreads $2 --port $((10000+$1)) --name mip$1  &

# read arguments from command line
args <- commandArgs(TRUE)
fas = args[1]
output_file = args[2]

# load the DECIPHER library in R
library(DECIPHER)

# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet(fas)

# perform the alignment
aligned <- AlignSeqs(seqs)

# write the alignment to a new FASTA file
writeXStringSet(aligned, file=output_file)

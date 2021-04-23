# This script merges fastq files in different directories, possibly from
# multiple sequencing runs, belonging to same libraries.
set -e
set -u
# set defaults
fastq_dirs=""
combined_dir=""
parallel=1
while getopts f:c:p: OPT; do
    case "$OPT" in
        f)
          fastq_dirs="$fastq_dirs $OPTARG";;
        c)
          combined_dir="$OPTARG";;

        p)
          parallel="$OPTARG";;

        *)
          echo "Invalid option: \"${OPT}\". \
          Usage: 'bash -f fastq_dir1 [[-f fastq_dir2]...] -c combined_dir \
          -p number_of_cpus"
          exit
    esac
done
mkdir $combined_dir || (echo "Could not create \"$combined_dir\". \
      Make sure an absolute path provided with -c option AND that the \
      Path provided does not already exist." && false)
cd $combined_dir
find $fastq_dirs -name '*_R1*' -printf '%f\n'|cut -d_ -f1 |sort|uniq > \
     sample_ids.txt


if command -v parallel &> /dev/null
then
  cat sample_ids.txt |parallel -I samplename -j $parallel find \
      $fastq_dirs -name samplename'*_R1_*' -exec cat '{}' + '>' \
      samplename_R1_001.fastq.gz
  cat sample_ids.txt |parallel -I samplename -j $parallel find \
      $fastq_dirs -name samplename'*_R2_*' -exec cat '{}' + '>' \
      samplename_R2_001.fastq.gz
else
  if [ $parallel > 1 ]
    then
      echo "'parallel' program is not available. A single cpu will be used!"
  fi
  while read samplename; do
      find $fastq_dirs -name ${samplename}'*_R1_*' |xargs cat > \
          ${samplename}_R1_001.fastq.gz
      find $fastq_dirs -name ${samplename}'*_R2_*' |xargs cat > \
          ${samplename}_R2_001.fastq.gz
  done < sample_ids.txt
fi
find $fastq_dirs -name '*_samples.tsv' -exec scp '{}' . ';'
cat *_samples.tsv> "concatenated_samples.tsv"
cat -n "concatenated_samples.tsv" | sort -uk2 |sort -n | cut -f2- > \
    "combined_samples.tsv"
a=$(grep sample_name combined_samples.tsv |grep sample_set | \
    grep replicate |wc -l)
if [ $a -ne 1 ]; then  echo "Warning: combined_samples.tsv file does not \
    contain a single header line which means the ..._samples.tsv files used \
    for each fastq directory had different headers. combined_samples.tsv file \
    created is not a valid sample sheet as it is now. You will need to create \
    the sample sheet file to use it downstream."
else echo "A combined sample sheet 'combined_samples.tsv' has \
 been created to be used with the merged fastq files in downstream analysis."
fi

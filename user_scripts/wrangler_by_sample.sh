#################################################
# set the ulimit high (necessary for big datasets)
#################################################
ulimit -n $(ulimit -Hn)

#################################################
# set the home directory as the current working directory
#################################################
newhome=$(pwd -P)

###############################################
# function to parse the yaml file edited by the user
# pulls out the location of the sif file, output directory, etc.
############################################
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}


#parse_yaml wrangler_by_sample.yaml
eval $(parse_yaml wrangler_by_sample.yaml)

input_sample_sheet_directory="$(dirname "${input_sample_sheet}")"





############################
# setup the run
##########################

# create output directory if it doesn't exist
mkdir -p $output_folder

#replace leading and trailing whitespace in variables (If I learn more unix I'll wrap this in a function or add to the yaml parser above):
project_resources="$(echo -e "${project_resources}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
output_folder="$(echo -e "${output_folder}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
input_sample_sheet_directory="$(echo -e "${input_sample_sheet_directory}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
fastq_dir="$(echo -e "${fastq_dir}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="-B $project_resources:/opt/project_resources
 -B $output_folder:/opt/analysis
 -B $input_sample_sheet_directory:/opt/input_sample_sheet_directory
 -B $fastq_dir:/opt/data
 -H $newhome"
 
snakemake_args="--cores $cpu_count --keep-going --rerun-incomplete --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec $singularity_bindings $miptools_sif snakemake \
   -s /opt/snakemake/wrangler_by_sample_setup.smk --unlock 

   singularity exec $singularity_bindings $miptools_sif snakemake \
   -s /opt/snakemake/wrangler_by_sample_finish.smk --unlock
}

#parse command line arguments to do the unlocking
while getopts "u" opt; do
        case ${opt} in
            u) unlock
               exit 1 ;;
        esac
    done

##################################
# Step 1: Set Up Wrangler Run
#################################
singularity exec \
 $singularity_bindings \
 $miptools_sif snakemake -s /opt/snakemake/wrangler_by_sample_setup.smk $snakemake_args

 ##################################
# Step 2: Finish Wrangler Run
#################################
singularity exec \
 $singularity_bindings \
 $miptools_sif snakemake -s /opt/snakemake/wrangler_by_sample_finish.smk $snakemake_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n
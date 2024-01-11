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


eval $(parse_yaml wrangler_by_sample.yaml)

function parse_sample_sheet_directory {
    readarray -d "/" -t strarr <<< "$input_sample_sheet"
    for (( n=1; n < ${#strarr[*]}-1; n++))
        do
         input_sample_sheet_directory+="/${strarr[n]}"
        done
    echo $input_sample_sheet_directory
}
input_sample_sheet_directory=$(parse_sample_sheet_directory)


############################
# setup the run
##########################

# create output directory if it doesn't exist
mkdir -p $output_folder

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="-B $project_resources:/opt/project_resources
 -B $output_folder:/opt/analysis
 -B $input_sample_sheet_directory:/opt/input_sample_sheet_directory
 -B $fastq_dir:/opt/data
 -B /home/charlie/projects/MIPTools_wrangler_in_sif/snakemake:/opt/snakemake
 -H $newhome"
 
snakemake_args="--cores $cpu_count --keep-going --rerun-incomplete --latency-wait 60"

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

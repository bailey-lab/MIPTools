########################################################
# README
# this file uses variant_calling.yaml for its parameters
#########################################################


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

eval $(parse_yaml variant_calling.yaml)

############################
# setup the run
##########################

# create output directory if it doesn't exist
mkdir -p $output_directory

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="-B $project_resources:/opt/project_resources
 -B $species_resources:/opt/species_resources
 -B $wrangler_directory:/opt/data
 -B $output_directory:/opt/analysis
 -H $newhome"
 
snakemake_args="--cores $processor_number --keep-going --rerun-incomplete --use-conda --latency-wait 60"

##################################
# Step 1: Check Run Stats
#################################
singularity exec \
 $singularity_bindings \
 $sif_file snakemake -s /opt/snakemake/02_check_run_stats.smk $snakemake_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n

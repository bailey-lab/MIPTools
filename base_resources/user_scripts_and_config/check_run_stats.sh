#!/bin/bash
#################################################
# set the ulimit high (necessary for big datasets)
#################################################
ulimit -n $(ulimit -Hn)

##############################################################
# set the home directory as the absolute (non-softlinked)
#current working directory and change directory to this folder
##############################################################
newhome=$(pwd -P)
cd $newhome

###################################
# import variables from yaml function
################################
yml (){
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

rmwt () {
   no_hash=$(echo -e $1 | sed -e 's/\#.*$//')
   echo -e $no_hash | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

eval $(yml config.yaml)

############################
# setup the run
##########################

# create output directory if it doesn't exist
mkdir -p $(rmwt $variant_calling_folder)

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="
 -B $(rmwt $project_resources):/opt/project_resources
 -B $(rmwt $species_resources):/opt/species_resources
 -B $(rmwt $wrangler_folder):/opt/user/wrangled_data
 -B $(rmwt $variant_calling_folder):/opt/user/stats_and_variant_calling
 -B /home/charlie/data/MIPTools/snakemake:/opt/snakemake
 -B /home/charlie/data/MIPTools/base_resources:/opt/resources
 -B /home/charlie/data/MIPTools/src:/opt/src
 -B $newhome:/opt/config"
 
snakemake_args="--cores $(rmwt $general_cpu_count) --keep-going --rerun-incomplete --use-conda --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec \
   $singularity_bindings \
   $(rmwt $miptools_sif) \
   snakemake -s /opt/snakemake/02_check_run_stats.smk --unlock 
}

#parse command line arguments to do the unlocking
while getopts "u" opt; do
        case ${opt} in
            u) unlock
               exit 1 ;;
        esac
    done

##################################
# Step 1: Check Run Stats
#################################
singularity exec \
  $singularity_bindings \
  $(rmwt $miptools_sif) \
  snakemake -s /opt/snakemake/02_check_run_stats.smk \
  $snakemake_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n

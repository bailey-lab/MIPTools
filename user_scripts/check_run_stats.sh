#################################################
# set the ulimit high (necessary for big datasets)
#################################################
ulimit -n $(ulimit -Hn)

###################################
# import variables from yaml function
################################
yml () {
   echo $(eval ./yq .$1 < config.yaml) 
}

############################
# setup the run
##########################

# create output directory if it doesn't exist
mkdir -p $(yml 'variant_calling_folder')

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="
 -B $(yml 'project_resources'):/opt/project_resources
 -B $(yml 'species_resources'):/opt/species_resources
 -B $(yml 'wrangler_folder'):/opt/data
 -B $(yml 'variant_calling_folder'):/opt/analysis
 -B /d/MIPTools/snakemake:/opt/snakemake
 -B $(pwd -P):/opt/config"
 
snakemake_args="--cores $(yml 'general_cpu_count') --keep-going --rerun-incomplete --use-conda --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
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
  $(yml 'miptools_sif') \
  snakemake -s /opt/snakemake/02_check_run_stats.smk \
  $snakemake_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n
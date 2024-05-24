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
mkdir -p $(yml 'wrangler_folder')

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="
 -B $(yml 'project_resources'):/opt/project_resources
 -B $(yml 'wrangler_folder'):/opt/analysis
 -B $(dirname $(yml 'input_sample_sheet')):/opt/input_sample_sheet_directory
 -B $(yml 'fastq_dir'):/opt/data
 -B /d/MIPTools/snakemake:/opt/snakemake
 -B /d/MIPTools/src:/opt/src
 -B $(pwd -P):/opt/config"
 
snakemake_args="--cores $(yml 'general_cpu_count') --keep-going --rerun-incomplete --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/wrangler_by_sample_setup.smk \
   --unlock 

   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/wrangler_by_sample_finish.smk \
   --unlock
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
 $(yml 'miptools_sif') \
 snakemake -s /opt/snakemake/wrangler_by_sample_setup.smk \
 $snakemake_args

 ##################################
# Step 2: Finish Wrangler Run
#################################
singularity exec \
 $singularity_bindings \
 $(yml 'miptools_sif') \
 snakemake -s /opt/snakemake/wrangler_by_sample_finish.smk \
 $snakemake_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n

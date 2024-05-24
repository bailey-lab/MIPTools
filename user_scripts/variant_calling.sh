#################################################
# set the ulimit high (necessary for big datasets)
#################################################
ulimit -n $(ulimit -Hn)

#################################################
# set the home directory as the current working directory
#################################################
cwd=$(pwd -P)

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
 -B $(yml 'wrangled_folder'):/opt/data
 -B $(yml 'variant_calling_folder'):/opt/analysis
 -B /d/MIPTools/snakemake/:/opt/snakemake
 -B $cwd:/opt/config"
 
snakemake_args="--cores $(yml 'general_cpu_count') --keep-going --rerun-incomplete --use-conda --latency-wait 60"
freebayes_args="--cores $(yml 'freebayes_cpu_count') --keep-going --rerun-incomplete --use-conda --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/02_check_run_stats.smk --unlock 

   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/03_generate_contigs.smk --unlock

   singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/04_run_freebayes.smk --unlock
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
  
###############################
# Step 2: Generate Contigs
##############################
singularity exec \
  $singularity_bindings \
  $(yml 'miptools_sif') \
  snakemake -s /opt/snakemake/03_generate_contigs.smk \
  $snakemake_args

###############################  
# Step 3: Run Freebayes
###############################
singularity exec \
   $singularity_bindings \
   $(yml 'miptools_sif') \
   snakemake -s /opt/snakemake/04_run_freebayes.smk \
   $freebayes_args

#################################
# confirm the ulimit settings #
################################
echo 'ulimit is' 
ulimit -n

#####################################################
# remove two unnecessary files that were generated #
#####################################################
rm -f snpEff_summary.html
rm -f snpEff_genes.txt

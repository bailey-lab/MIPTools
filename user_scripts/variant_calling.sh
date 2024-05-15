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

#replace leading and trailing whitespace in variables (If I learn more unix I'll wrap this in a function or add to the yaml parser above):
project_resources="$(echo -e "${project_resources}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
species_resources="$(echo -e "${species_resources}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
wrangler_directory="$(echo -e "${wrangler_directory}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
output_directory="$(echo -e "${output_directory}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

# define singularity bindings and snakemake arguments to be used each time snakemake is called
singularity_bindings="-B $project_resources:/opt/project_resources
 -B $species_resources:/opt/species_resources
 -B $wrangler_directory:/opt/data
 -B $output_directory:/opt/analysis
 -H $newhome"
 
snakemake_args="--cores $processor_number --keep-going --rerun-incomplete --use-conda --latency-wait 60"
freebayes_args="--cores $freebayes_threads --keep-going --rerun-incomplete --use-conda --latency-wait 60"

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   echo "unlocking"
   singularity exec $singularity_bindings $miptools_sif snakemake \
   -s /opt/snakemake/02_check_run_stats.smk --unlock 

   singularity exec $singularity_bindings snakemake \
   -s /opt/snakemake/03_generate_contigs.smk --unlock

   singularity exec $singularity_bindings snakemake \
   -s /opt/snakemake/04_run_freebayes.smk --unlock
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
 $sif_file snakemake -s /opt/snakemake/02_check_run_stats.smk $snakemake_args
  
###############################
# Step 2: Generate Contigs
##############################
singularity exec \
  $singularity_bindings \
  $sif_file snakemake -s /opt/snakemake/03_generate_contigs.smk $snakemake_args

###############################  
# Step 3: Run Freebayes
###############################
singularity exec \
  $singularity_bindings \
  $sif_file snakemake -s /opt/snakemake/04_run_freebayes.smk $freebayes_args

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
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

# set the mip_version to avoid any conflicts between the shell script and the
# version of miptools in the sif file
mip_version=v0.5.0
check_for_sif(){
    no_sif=false
    if [[ ! -f $miptools_sif ]]; then
        echo ""
        echo "error: the path to the sif in the config file cannot be found, please check on it"
        no_sif=true
    fi
    if [[ ! $miptools_sif == *"$mip_version"*  ]]; then
        echo ""
        echo "it looks like you do not have a version $mip_version sif selected in your config file"
        echo "please edit the config file to choose a sif file version $mip_version"
        no_sif=true
    fi
}

###################################
# import variables from yaml function
###################################
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

establish_binds () {
    miptools_sif=$(rmwt $miptools_sif)
    project_resources=$(rmwt $project_resources)
    species_resources=$(rmwt $species_resources)
    input_sample_sheet_directory=$(rmwt $(dirname $input_sample_sheet))
    fastq_dir=$(rmwt $fastq_dir)
    wrangler_folder=$(rmwt $wrangler_folder)
    variant_calling_folder=$(rmwt $variant_calling_folder)
    prevalence_metadata=$(rmwt $prevalence_metadata)

    singularity_bindings="-B $newhome:/opt/config"
    if [ ! -z $project_resources ]; then singularity_bindings="$singularity_bindings
        -B $project_resources:/opt/project_resources"; fi
    if [ ! -z $species_resources ]; then singularity_bindings="$singularity_bindings
        -B $species_resources:/opt/species_resources"; fi
    if [ ! -z $input_sample_sheet_directory ]; then singularity_bindings="$singularity_bindings
        -B $input_sample_sheet_directory:/opt/input_sample_sheet_directory"; fi
    if [ ! -z $fastq_dir ]; then singularity_bindings="$singularity_bindings
        -B $fastq_dir:/opt/fastq_dir"; fi
    if [ ! -z $wrangler_folder ]; then singularity_bindings="$singularity_bindings
        -B $wrangler_folder:/opt/user/wrangled_data"
        mkdir -p $wrangler_folder; fi
    if [ ! -z $variant_calling_folder ]; then singularity_bindings="$singularity_bindings
        -B $variant_calling_folder:/opt/user/stats_and_variant_calling"
        mkdir -p $variant_calling_folder; fi
    if [ ! -z $prevalence_metadata ]; then singularity_bindings="$singularity_bindings
        -B $prevalence_metadata:/opt/user/prevalence_metadata"; fi
}

eval $(yml config_v*.yaml)
establish_binds

##########################################
# optional: unlock a crashed snakemake run
##########################################

unlock() {
   singularity run \
     --app unlock_snakemake \
     $singularity_bindings \
     $miptools_sif
}

#parse command line arguments to do the unlocking
while getopts "u" opt; do
        case ${opt} in
            u) unlock
               exit 1 ;;
        esac
    done

check_for_sif
if [[ $no_sif = true ]]; then exit; fi

singularity run \
  --app wrangler \
  $singularity_bindings \
  $miptools_sif \
  -c $general_cpu_count

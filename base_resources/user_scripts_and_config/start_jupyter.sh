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
mkdir -p $(rmwt $wrangler_folder)

# define singularity bindings and snakemake arguments to be used each time snakemake is called
if [ -z "$prevalence_metadata" ]; then # don't include prevalence data if user has left it blank
   singularity_bindings="
    -B $(rmwt $project_resources):/opt/project_resources
    -B $(rmwt $species_resources):/opt/species_resources
    -B $(rmwt $wrangler_folder):/opt/user/wrangled_data
    -B $(rmwt $variant_calling_folder):/opt/user/stats_and_variant_calling
    -B $newhome:/opt/config"
else
   singularity_bindings="
    -B $(rmwt $project_resources):/opt/project_resources
    -B $(rmwt $species_resources):/opt/species_resources
    -B $(rmwt $wrangler_folder):/opt/user/wrangled_data
    -B $(rmwt $variant_calling_folder):/opt/user/stats_and_variant_calling
    -B $(rmwt $prevalence_metadata):/opt/user/prevalence_metadata
    -B $newhome:/opt/config"
fi

singularity run \
  $singularity_bindings \
  --app jupyter $(rmwt $miptools_sif) -d /opt/user

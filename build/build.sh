# set the version number. this will be used to update the version number of the output 
# sif file and user scripts and config files 
mip_version="0.5.0"

# - change the config file name in the bash scripts
# - change the sif file name in the config.yaml

# - change the name of the user scripts to match the current version number
scripts_path="../base_resources/user_scripts_and_config"
mv $scripts_path/check_run_stats*\.sh "$scripts_path/check_run_stats_v$mip_version.sh"
mv $scripts_path/start_jupyter*\.sh "$scripts_path/start_jupyter_v$mip_version.sh"
mv $scripts_path/variant_calling*\.sh "$scripts_path/variant_calling_v$mip_version.sh"
mv $scripts_path/wrangler_by_sample*\.sh "$scripts_path/wrangler_by_sample_v$mip_version.sh"
mv $scripts_path/config*\.yaml "$scripts_path/config_v$mip_version.yaml"
mv $scripts_path/run_miptools*\.sh "$scripts_path/run_miptools_v$mip_version.sh"

# - change the snakemake file name in the config file
sed -i "s/miptools_v.*\.sif/miptools_v$mip_version\.sif/g" $scripts_path/config_v$mip_version.yaml
sed -i "s/mip_version=.*/mip_version=v$mip_version/g" $scripts_path/run_miptools_v$mip_version.sh
sed -i "s/mip_version=.*/mip_version=v$mip_version/g" $scripts_path/wrangler_by_sample_v$mip_version.sh
sed -i "s/mip_version=.*/mip_version=v$mip_version/g" $scripts_path/check_run_stats_v$mip_version.sh
sed -i "s/mip_version=.*/mip_version=v$mip_version/g" $scripts_path/variant_calling_v$mip_version.sh
sed -i "s/mip_version=.*/mip_version=v$mip_version/g" $scripts_path/start_jupyter_v$mip_version.sh
sed -i "s#\/opt\/config\/config.*\.yaml#\/opt\/config\/config_v$mip_version\.yaml#g" ../base_resources/jupyter_notebooks/prevalence_plotting.ipynb


build_conda() {
    if [ ! -f conda.sif ]; then
        sudo singularity build conda.sif conda.def
    fi
}

build_wrangler() {
    if [ ! -f wrangler.sif ]; then
        sudo singularity build wrangler.sif wrangler.def
    fi
}

build_conda &
build_wrangler &

wait

sudo singularity build miptools_v$mip_version.sif MIPTools.def

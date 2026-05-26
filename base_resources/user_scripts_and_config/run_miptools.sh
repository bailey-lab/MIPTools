#!/usr/bin/env bash

# do not edit version here, the build script will apply the version defined there at build time
VERSION=dev

# set the ulimit high in case there are a very large number of files
ulimit -n $(ulimit -Hn)

# change directory to the absolute path of the current working directory to
# avoid binding issues with symlinks
newhome=$(pwd -P)
cd $newhome

# import bindmounts and ensure that they are absolute paths
current_section=""
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*($|#) ]] && continue

    if [[ "$line" =~ ^[[:space:]]*\[([^]]+)\] ]]; then
        current_section="${BASH_REMATCH[1]//./_}_"
        continue
    fi

    line=$(sed -E 's/("[^"]*"|[0-9]+)[[:space:]]*#.*/\1/' <<< "$line")

    if [[ "$line" =~ ^[[:space:]]*([a-zA-Z_][a-zA-Z0-9_]*)[[:space:]]*=[[:space:]]*(\"(.*)\"|[0-9]+)[[:space:]]*$ ]]; then
        declare "${current_section}${BASH_REMATCH[1]}=${BASH_REMATCH[3]:-${BASH_REMATCH[2]}}"
    fi
done < "config.toml"

miptools_sif="$(readlink -f "$universal_inputs_miptools_sif")"
project_resources="$(readlink -f "$universal_inputs_project_resources_directory")"
input_sample_sheet="$(readlink -f "$wrangler_inputs_input_sample_sheet")"
fastq_dir="$(readlink -f "$wrangler_inputs_fastq_dir")"
species_resources="$(readlink -f "$variant_calling_inputs_species_resources")"
prevalence_metadata_file="$(readlink -f "$prevalence_summary_inputs_prevalence_metadata_file")"
wrangler_directory="$(readlink -f "$variant_calling_inputs_wrangler_directory")"

singularity_options=(-B "$newhome:/opt/user")
    [[ -n "$project_resources" && -d "$project_resources" ]] &&
        singularity_options+=(-B "$project_resources:/opt/project_resources")
    [[ -n "$species_resources" && -d "$species_resources" ]] &&
        singularity_options+=(-B "$species_resources:/opt/species_resources")
    [[ -n "$input_sample_sheet" && -f "$input_sample_sheet" ]] &&
        singularity_options+=(-B "$input_sample_sheet:/opt/$(basename "$input_sample_sheet")")
    [[ -n "$fastq_dir" && -d "$fastq_dir" ]] &&
        singularity_options+=(-B "$fastq_dir:/opt/fastq_dir")
    [[ -n "$wrangler_directory" && -d "$wrangler_directory" ]] &&
        singularity_options+=(-B "$wrangler_directory:/opt/wrangled_data")
    [[ -n "$prevalence_metadata_file" && -f "$prevalence_metadata_file" ]] &&
        singularity_options+=(-B "$prevalence_metadata_file:/opt/$(basename "$prevalence_metadata_file")")

    singularity_options+=(-B "$HOME/data/MIPTools/snakemake:/opt/snakemake")
    singularity_options+=(-B "$HOME/data/MIPTools/src:/opt/src")


# set the version to avoid any conflicts between the shell script and the
# version of miptools in the sif file
check_for_sif(){
    if [[ ! -e $miptools_sif ]]; then
        echo ""
        echo "error: the path to the sif in the config file cannot be found, please check on it"
        exit 1
    fi
    if [[ ! $(singularity exec $miptools_sif printenv VERSION) == "$VERSION"  ]]; then
        echo ""
        echo "it looks like you do not have a version $VERSION sif selected in your config file"
        echo "please edit the config file to choose a sif file version $VERSION"
        exit 1
    fi
}

# give user options to edit config or run different pipelines
main_menu (){
    echo ""
    echo "Enter a number to select one of the following actions"
    PS3='Choose an option: '
    options=(
        "edit config" "run wrangler" "check run stats" \
        "variant calling" "start jupyter" "unlock snakemake" \
    )
    select opt in "${options[@]}"
    do
        case $opt in
            "edit config")
                ./micro config_$VERSION.yaml
                break
                ;;
            "run wrangler")
                check_for_sif
                singularity run \
                    --app wrangler \
                    "${singularity_options[@]}" \
                    "$miptools_sif" \
                    -c "$wrangler_settings_cpu_count"
                break
                ;;
            "check run stats")
                check_for_sif
                singularity run \
                    --app check_run_stats \
                    "${singularity_options[@]}" \
                    "$miptools_sif" \
                    -c "$check_run_stats_cpu_count"
                break
                ;;
            "variant calling")
                check_for_sif
                singularity run \
                    --app variant_calling \
                    "${singularity_options[@]}" \
                    "$miptools_sif" \
                    -c "$check_run_stats_cpu_count" \
                    -f "$prevalence_calling_freebayes_cpu_count"
                break
                ;;
            "start jupyter")
                check_for_sif
                singularity run \
                    --app jupyter \
                    --env prevalence_metadata_file="$prevalence_metadata_file" \
                    "${singularity_options[@]}" \
                    "$miptools_sif" \
                    -d /opt/user
                break
                ;;
            "unlock snakemake")
                singularity run \
                    --app unlock_snakemake \
                    "${singularity_options[@]}" \
                    "$miptools_sif"
                break
                ;;
            *) echo "invalid option $REPLY";;
        esac
    done
}

main_menu

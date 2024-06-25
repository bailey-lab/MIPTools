sif_file=miptools_dev.sif
project_resources=input_data/DR23K_project_resources
species_resources=input_data/pf_species_resources
prevalence_metadata=input_data/metadata_files
wrangler_folder=wrangled_data
variant_output=stats_and_variant_calling

singularity run \
  -B $project_resources:/opt/project_resources \
  -B $species_resources:/opt/species_resources \
  -B $wrangler_folder:/opt/data \
  -B $variant_output:/opt/analysis \
  -B $prevalence_metadata:/opt/prevalence_metadata \
  --app jupyter $sif_file
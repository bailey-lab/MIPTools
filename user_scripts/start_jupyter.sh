#!/bin/bash

project_resources='/home/charlie/resources/miptools_test-data/DR1_project_resources'
species_resources='/home/charlie/resources/miptools_test-data/pf_species_resources'
wrangler_folder='/home/charlie/resources/miptools_test-data_wrangler'
variant_output='/home/charlie/projects/MIPTools/user_scripts/local_files/variant_jupyter'
sif_file='/home/charlie/resources/mip.sif'
prevalence_metadata='/home/charlie/resources/choropleth_inputs'

singularity run \
  -B $project_resources:/opt/project_resources \
  -B $species_resources:/opt/species_resources \
  -B $wrangler_folder:/opt/data \
  -B $variant_output:/opt/analysis \
  -B $prevalence_metadata:/opt/prevalence_metadata \
  --app jupyter $sif_file

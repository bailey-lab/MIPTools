#!/bin/bash

project_resources='/home/charlie/resources/miptools_test-data/DR1_project_resources'
species_resources='/home/charlie/resources/miptools_test-data/pf_species_resources'
wrangler_folder='/home/charlie/projects/user_scripts/wrangler'
variant_output='/home/charlie/projects/user_scripts/variant_jupyter'
sif_file='/home/charlie/projects/MIPTools/MIP1200.sif'

singularity run \
  -B $project_resources:/opt/project_resources \
  -B $species_resources:/opt/species_resources \
  -B $wrangler_folder:/opt/data \
  -B $variant_output:/opt/analysis \
  --app jupyter $sif_file

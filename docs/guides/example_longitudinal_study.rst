=============================
An example longitudinal study
=============================

.. note:: 
	
	This guide is currently under active development and should not be used by
	new users of this software until it has been finished and validated.

This dataset consists of cherrypicked samples from the Conrad et al. article
"Evolution of Partial Resistance to Artemisinins in Malaria Parasites in
Uganda." This article was published in 2023 in the New England Journal of
Medicine with PMID 37611122.

Although the original study consists of thousands of samples from sixteen
districts of Uganda across seven years, we needed a much smaller dataset for
this tutorial, that would be easy to download, analyze, and interpret. We
therefore chose three years (2016, 2019, and 2022) from four districts (Agago,
Amolatar, Kole, and Kabale), using only the informative reads from 220 samples
chosen to yield prevalences that reproduce the main findings of the original
study. We also chose to include only MIPs covering Kelch13, dhps, dhfr, crt,
and mdr1. The breakdown of these 220 samples is as follows:

- 20 samples from each of the four districts in 2016
- 17 samples from Agago from 2019
- 13 samples from Amolatar from 2019
- 14 samples from Kabale from 2019
- 20 samples from Kole from 2019
- 20 samples from Agago, Kabale, and Kole from 2022
- 16 samples from Amolatar from 2022

The example dataset can be obtained here:
(Download link not built yet)

This data is organized into five main folders:
	- **pf_species_resources:** This folder includes an indexed copy of the
	  Pf3D7 falciparum genome, gene annotations, common SNPs, and a directory of
	  key files.

	- **barebones_DR_project_resources:** This is a special project resources
	  folder that includes MIPs for K13, dhps, dhfr, crt, and mdr1 genes only.

	- **metadata_files:** These include 2016_metadata.tsv, 2019_metadata.tsv,
	  and 2022_metadata.tsv, with collection site information for each sample.
	  These also include a geojson file, which we will use later for mapping
	  mutation prevalences.

	- **fastq_files:** These are the raw paired end sequencing reads associated
	  with each sample. Because there are 220 samples and two read files (one
	  forward and one reverse) associated with each sample, you should see 440
	  fastq files in this folder.

	- **sample_sheets:** This folder contains sample sheets for each year,
	  including 2016_sample_sheet.tsv, 2019_sample_sheet.tsv, and
	  2022_sample_sheet.tsv

Wrangling
=========

This step converts raw reads into error-corrected haplotypes, and collapses
multiple reads that are duplicates of the same original sampled molecule into a
single representative consensus sequence. This step also reports the number of
unique molecular identifiers (UMIs) associated with every haplotype for every
MIP for every sample.

| We've provided a bash script for converting the yaml settings into instructions for the wrangler. You can obtain the bash script for wrangling with this command (put it in the same folder as the settings yaml file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/wrangler_by_sample.sh`

| After editing the config.yaml file, you can execute the wrangler script with:
| :code:`bash wrangler_by_sample.sh`

Interpreting the wrangler output
--------------------------------

When the wrangler is finished, you should open the folder labeled  

Checking Run Statistics
=======================

| You can obtain the script for checking run stats here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/check_run_stats.sh`

| After editing the relevant config.yaml file sections, you can execute the check_run_stats script with:
| :code:`bash check_run_stats.sh`


Variant Calling
===============
| You can obtain the script for variant calling here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/variant_calling.sh`

| After editing the relevant config.yaml file sections, you can execute the variant_calling script with:
| :code:`bash variant_calling.sh`

prevalence Calling
==================
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
https://baileylab.brown.edu/MIPTools/download/test_data_new_tutorial.tar.gz

This data is organized into four main folders:
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

There is also a sample sheet (cherrypicked_sample_sheet.tsv) that has
information about which samples are associated with each 'experiment'. In this
case the data is separated into three experiments (known as sample_sets):
- PRX-00 (2016 dataset)
- PRX-04 (2019 dataset)
- PRX-07 (2022 dataset)

Editing Settings
================
| For convenience, settings can be passed in to all steps via a single shared
| yaml file. Later in the tutorial, we'll show some more advanced usage options
| available for troubleshooting and passing more customizable inputs. For now,
| you can obtain an example simple settings file with this command:
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/config.yaml`
| After downloading, open the file for editing with a text editor and make sure
| to **follow the instructions in this file**, editing it to contain the
| correct paths to the files you downloaded above (including project resources,
| species resources, sample sheet, and sif files), as well as the locations
| where you'd like the output to be sent.

| Make sure to edit the settings for all steps that you intend to run before
| running them. If you open the config file, you should see which settings (in
| the config.yaml file downloaded above) pertain to each of the steps below.


Wrangling
=========

This step converts raw reads into error-corrected haplotypes, and collapses
multiple reads that are duplicates of the same original sampled molecule into a
single representative consensus sequence. This step also reports the number of
unique molecular identifiers (UMIs) associated with every haplotype for every
MIP for every sample.

For this tutorial, you'll need to wrangle three sample_sets: The 2016 dataset
(PRX-00), the 2019 dataset (PRX-04), and the 2022 dataset (PRX-07). This
pipeline can analyze multiple sample sets at once, so you'll be inputting
PRX-00,PRX-04,PRX-07 in the sample_set variable of config.yaml.

| We've provided a bash script for converting the yaml settings into
| instructions for the wrangler. You can obtain the bash script for wrangling
| with this command (put it in the same folder as the settings yaml file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/wrangler_by_sample.sh`

| After editing the config.yaml file, you can execute the wrangler script with:
| :code:`bash wrangler_by_sample.sh`

Interpreting the wrangler output
--------------------------------



Checking Run Statistics
=======================

| You can obtain the script for checking run stats here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/check_run_stats.sh`

| After editing the relevant config.yaml file sections, you can execute the check_run_stats script with:
| :code:`bash check_run_stats.sh`

Interpreting the run statistics
-------------------------------


Variant Calling
===============
| You can obtain the script for variant calling here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/variant_calling.sh`

| After editing the relevant config.yaml file sections, you can execute the variant_calling script with:
| :code:`bash variant_calling.sh`

Interpreting the variant calling
--------------------------------

prevalence Calling
==================

Interpreting the prevalence calling
-----------------------------------
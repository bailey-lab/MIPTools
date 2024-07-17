=============================
An example longitudinal study
=============================

.. note:: 
	
	This guide is currently under active development and should not be used by
	new users of this software until it has been finished and validated. This
	tutorial assumes that you have access to a linux computer that has a copy
	of singularity installed on it, and that you have basic familiarity with
	running commands in a terminal and editing plain text files.


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
study. We also chose to include only MIPs covering essential drug resistance
genes (Kelch13, dhps, dhfr, crt, mdr1, PM2, and PM3). The breakdown of these
220 samples is as follows:

- 20 samples from each of the four districts in 2016
- 17 samples from Agago in 2019
- 13 samples from Amolatar in 2019
- 14 samples from Kabale in 2019
- 20 samples from Kole in 2019
- 20 samples from Agago, Kabale, and Kole in 2022
- 16 samples from Amolatar in 2022

| The tutorial dataset can be downloaded from here:
| https://baileylab.brown.edu/MIPTools/download/tutorial_dataset.tar.gz
| The dataset can be extracted with this command:
| :code:`tar -xvzf tutorial_dataset.tar.gz`

| You can obtain a copy of our latest sif file from here:
| https://baileylab.brown.edu/MIPTools/download/miptools_dev.sif
| This includes all executable programs needed for analysis

For this tutorial to work without modifying any of the settings, you'll need to
move the sif file into the tutorial_dataset folder.

| To obtain a copy of the config files used for this tutorial, cd into the
 tutorial_dataset folder and execute this command:
| :code:`singularity run -B $(pwd -P):/opt/config miptools_dev.sif`
 (remember that a prerequisite for this tutorial is an installed copy of
 singularity on your computer).

| In general, when you analyze any dataset, you should cd into a folder and run
 the
| :code:`singularity run -B $(pwd -P):/opt/config /path/to/your/downloaded/miptools_dev.sif`
 step to download all config files. Unlike in the tutorial, you'll then need to
 modify the config.yaml file to point to your input files and parameters. If the
 computer you're on is a remote computer, you can edit the config.yaml file
 using the text editor "micro" using this command:
| :code:`./micro config.yaml` 
| micro offers mouse support (so you can click a field of text to start editing
 it) and allows you to copy with ctrl-C, paste with ctrl-V, save with ctrl-S,
 and quit with ctrl-Q.

Understanding the input data
----------------------------
If you cd to the tutorial_dataset folder, you'll see a subfolder called
input_data. This folder is organized into four main folders:

- **pf_species_resources:** This folder includes an indexed copy of the
  Pf3D7 falciparum genome, gene annotations, common SNPs, and a directory of
  key files.
- **DR23K_project_resources:** Contains information about the mip panel (in
  this case a panel called DR23K). Bound internally to :code:`/opt/project_resources`.
  Includes:
    
  - A targets.tsv file with the genomic coordinates of any protein-coding
    mutations that are of particular interest.
  - A geneid_to_genename.tsv file that converts gene IDs (from  into common gene names.
  - A mip_ids folder that contains the mip arms of MIPs that target regions of the
    genome that are of interest.
- **metadata_files:** These include 2016_metadata.tsv, 2019_metadata.tsv,
  and 2022_metadata.tsv, with collection site information for each sample.
  These also include a geojson file, which we will use later for mapping
  mutation prevalences.
- **fastq_files:** These are the raw paired end sequencing reads associated
  with each sample. Because there are 220 samples and two read files (one
  forward and one reverse) associated with each sample, you should see 440
  fastq files in this folder.

The input_data folder also contains a sample sheet called
cherrypicked_sample_sheet.tsv that has information about which samples are
associated with a given project (because sometimes multiple projects are
sequenced at the same time on a sequencer). These are known as sample_sets.
For this dataset, there are three sample_set designations:

- PRX-00 (2016 dataset)
- PRX-04 (2019 dataset)
- PRX-07 (2022 dataset)

Editing Settings
================
For convenience, settings can be passed in to all steps via a single shared
yaml file, called config.yaml. We've edited these settings to run with this
tutorial dataset, but we highly recommend opening this file for editing with a
text editor and reading the comments thoroughly - this file specifies inputs
and outputs and controls all aspects of the behavior of the program.

Wrangling
=========
This step converts raw reads into error-corrected haplotypes, and collapses
multiple reads that are duplicates of the same original sampled molecule into a
single representative consensus sequence. Consensus sequences that are
identical to each other are named as haplotypes. This step also reports the
number of unique molecular identifiers (UMIs) associated with every haplotype
for every MIP for every sample. This UMI count is equivalent to the number of
times each type of genetic sequence was seen in each original sample (prior to
PCR amplification).


| After changing directory to tutorial_dataset, you can execute the wrangler
 script with:
| :code:`bash wrangler_by_sample.sh`


Interpreting the wrangler output
--------------------------------
In the pre-configured settings, output of the wrangling step will go to a
folder called "wrangled_data." This is controlled by the wrangler_folder
variable in the config.yaml file.  If you'd like to see the "raw" outputs of
the wrangler, the main output file is called allInfo.tsv.gz and it can be
unzipped for reading in tabular format. Each row gives UMI counts, genetic
sequence, and statistics associated with a single haplotype associated with a
particular MIP of a particular sample.

Later steps will parse this table into graphical formats that will be easier to
interpret.

If you'd like to learn more about how to directly interpret the wrangler
output, you can check out the
:ref:`advanced_wrangler_interpretation` page.

Checking Run Statistics
=======================

This step converts the wrangler output data into graphs and tables that tell a
user which samples and mips succeeded and which may need to be run again.

| While in the folder tutorial_dataset, you can execute the check_run_stats
 command with:
| :code:`bash check_run_stats.sh`

| Alternatively, you can create a folder in the tutorial_dataset folder called
 "stats_and_variant_calling" and then run this jupyter script:
| :code:`bash start_jupyter.sh`
 Click the "analysis" folder, and click the "check_run_stats.ipynb" file. Follow
 the instructions in the notebook.



Interpreting the run statistics
-------------------------------
In the pre-configured settings, output of the check_run_stats step will go to a
folder called "stats_and_variant_calling." This is controlled by the
variant_calling_folder variable in the config.yaml file. There are a few key
output files that are useful to examine:

- **umi_heatmap.html**: This file can be downloaded and opened with a web
  browser. It includes The names of all samples (y-axis) and the names of all
  MIPs (x-axis). In the tutorial dataset, DR23K has 121 mips, and in the
  tutorial dataset, there are 220 samples. Not all of these samples are
  visible, but if you zoom in (by clicking and dragging) you can see all
  labels. By hovering over a box on the heatmap, you can see how many UMIs are
  associated with each sample and each MIP.
  
  - If you look for bright rows in this dataset, you can see that some samples,
    such as KO-07-001-PRX-07-1, performed extremely well across almost all MIPs,
    with UMI counts >2^12 for almost all MIPs, while if you look for dim rows,
    you might notice that other samples, such as AM-07-89-PRX-07-1, performed
    very poorly with UMI counts <2^4 for almost all MIPs.
  - Similarly, if you look for bright columns in this dataset, you might notice
    that most MIPs perform relatively well, while a few have very dim columns
    and perform poorly across all samples (e.g. crt_S0_Sub0_mip9).

- **umi_count_vs_probe_coverage.html**: This file is also meant to be
  downloaded and opened with a web browser. The x-axis represents total UMIs
  for a sample, while the y-axis represents number of MIPs having at least 10
  UMIs within that sample. By hovering over individual points, you can see which
  samples have a large number of MIPs that have more than 10 UMIs (indicating
  that they are well-sampled) and which do not. A "good" dataset will show a few
  points forming a vertical line along the y-axis line near x=10*UMI_count.
  Since we have 121 MIPs, our vertical line should occur at x=1,210). In a
  "good" dataset, almost all samples would have 10 UMIs for almost all MIPs, and
  the vast majority of points should form a horizontal line with a y-value near
  the number of MIPs (121 in our case). For the tutorial dataset, MIPs are not
  performing very well - most samples appear along the vertical line, and the
  vertical line extends well past x=1,210, indicating uneven coverage. Even as
  UMIs increase well past the theoretical minimum, this is not enough to
  saturate most MIPs with 10 UMIs. The line doesn't become horizontal until
  x=50,000, indicating that 50,000 UMIs are needed to start having good UMI
  coverage for nearly all MIPs. Hardly any samples approach the y=121 line. The
  best performing samples retrieve 118 MIPs (out of 121), so there is no sample
  that recovered all 121 MIPs. Many of these samples should be redone (either
  repooled or re-captured).
- **repool.csv**: This file gives recommendations regarding which samples are
  "Complete" (if at least 95% of MIPs have at least 10 UMIs), which should be
  "Repooled" (if the sample is not "Complete" and the number of reads is
  similar to the number of UMIs) and which should be "Recaptured" (if the
  sample is not "Complete" and the number of reads is much higher than the
  number of UMIs). Thresholds for these recommendations are based on the repool
  spreadsheet settings from the config.yaml file. In the tutorial dataset, 21
  of the samples are "Complete", 53 of the samples are "Recapture" and 146 of
  the samples are "Repool". Out of 8,904,984 reads, 6,119,806 reads, or 68.7%,
  came from the 21 "Complete" samples. The "Complete" samples monopolized the
  sequencing reads, and used 68.7% of the reads despite making up only 17.2% of
  the samples. The "Recapture" samples have plenty of sequencing reads for each
  UMI, but they all come from only a few UMIs. By repeating the MIP capture
  reactions for these samples, hopefully more UMIs will be recovered. After
  repeating the MIP capture reactions on the "Recapture" samples, by
  re-sequencing a pool of the 199 samples that are not "Complete", 68.7% of the
  reads should be freed up to give more sequencing depth to the remaining
  samples. This process can be repeated until almost all samples are "Complete".
  Reads from earlier runs can be pooled with reads from later runs so that reads
  from samples that are not "Complete" are not wasted.

Variant Calling
===============
This step takes haplotypes (from the Wrangling step) and maps them to the
reference genome (in this case 3D7). This step uses an annotation file and a
list of mutations of interest to name all of the mutations that were seen in the
dataset, as well as count the number of UMIs that were associated with the
reference genome and the number of UMIs that were associated with the mutant in
each sample.

| After editing the relevant config.yaml file sections you can execute the
 variant_calling script (while in the tutorial_dataset folder) with:
| :code:`bash variant_calling.sh`

| Alternatively, you can create a folder in the tutorial_dataset folder called
 "stats_and_variant_calling" and then run this jupyter script:
| :code:`bash start_jupyter.sh`
 Click the "analysis" folder, and click the "variant_calling.ipynb" file. Follow
 the instructions in the notebook.


Interpreting the variant calling
--------------------------------
In the pre-configured settings, output of the check_run_stats step will go to a
folder called "stats_and_variant_calling." This is controlled by the
variant_calling_folder variable in the config.yaml file. There are a few key
output files that are useful to examine:

- **variants.vcf.gz**: Each row of this file is a genomic position. Each column
  is an individual sample. For the rows that have mutations, the codes
  (described in the header) show various statistics for each mutation, such as
  number of UMIs supporting the mutation, number of UMIs that covered the
  region, and the confidence of the variant caller (in this case Freebayes) that
  the mutation is real. This file can be used by many downstream applications
  (such as Identity by Descent) that expect VCF files as inputs.
- **AA tables files**: This tutorial dataset examines drug resistance mutations.
  The files below describe the number of UMI counts associated with each
  mutation. Every column is a different mutation, and every row is a sample.

  - *coverage_AA_table.csv* - The number of total UMIs associated with the
    region of the genome covered by the mutation in a sample.
  - *reference_AA_table.csv* - The number of UMIs associated with the reference
    allele in a sample.
  - *alternate_AA_table.csv* - The number of UMIs associated with the mutated
    allele in a sample.

The within sample allele frequency of a mutation is obtained by dividing the
alternate UMI count in a sample be the coverage UMI count of the sample, and the
prevalence of a mutation is obtained by counting the number of samples that meet
some minimum coverage UMI count and that have an alternate UMI count greater
than some minimum level. By setting a minimum UMI coverage of three and a
minimum UMI alternate count of one, we can see how many samples meet these
criteria. As two examples:

- The crt-Asn75Glu mutation (column BG) has 183 samples that have values of at
  least 3 in the coverage_AA_table, and 11 of these samples have values of at
  least 1 in the alternate_AA_table. The overall prevalence of the crt-Asn75Glu
  mutation at these coverage and alternate thresholds is 11/183 or 6%.
- The dhfr-ts-Cys59Arg mutation (column D) has 199 samples that have values of at
  least 3 in the coverage_AA_table, and 193 of these samples have values of at
  least 1 in the alternate_AA_table. The overall prevalence of the
  dhfr-ts-Cys59Arg mutation at these coverage and alternate thresholds is
  193/197, or 97%

In the next section, we'll use metadata files to perform a more detailed
prevalence calling for individual regions and individual years.

prevalence Calling
==================
For this step, you'll need to open a Jupyter notebook. If you change directory
to the tutorial_dataset folder, you can launch the jupyter notebook with this
command:

| :code:`bash start_jupyter.sh`

After launching the jupyter notebook, leave the terminal window open. If you're
running the Jupyter notebook on a remote server, you may need to use port
forwarding to view the output Jupyter notebook. The command for this is shown at
the top of the Jupyter notebook output screen, and needs to be executed on your
local computer. After executing this, you can click one of the links on the
running Jupyter notebook screen. The link will open on your web-browser. Click
on the "analysis" folder link, and then click the link labeled
"prevalence_plotting." Follow the instructions in the notebook.
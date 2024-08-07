====================
Abbreviated tutorial
====================
This abbreviated tutorial provides a minimal version of instructions for
people who want to get a "working version" up and running quickly to explore
the data for themselves. A little background is still needed.

Background
==========
MIPTools supports several computational steps, including:
	- MIP design: This is for choosing resions of the genome that you'd like to
	  target.

	- Wrangling: This is for finding what haplotypes are associated with each
	  targeted region and the number of times each haplotype was seen in each
	  sample. This can be thought of as the "core" purpose of MIPTools. The output
	  is a tab delimited file called **allInfo.tsv.gz**

	- stat-checking: This is for figuring out which samples and which targeted
	  regions (aka MIPs) performed well and which did not. There are several
	  graphical outputs and comma separated files produced at this stage. Among
	  the most important are **UMI_counts.csv** (how many times each targeted
	  region of the genome was seen in each sample) and **repool.csv** (which MIPs
	  need to be re-sequenced or repooled in each sample). **umi_heatmap.html** is
	  also useful for visualizing the performance of each MIP in each sample
	  (open this with a web browser).

	- variant calling: This is for estimating the number of times that each MIP
	  was associated with a mutation in each sample. The outputs are a VCF file
	  containing mutated genomic positions for each sample and several tables
	  that annotate mutations with associated amino acid changes for each sample.
	  For the tutorial, the main outputs of interest are three tables that can be
	  used to infer which samples contain each mutation:

	  - *coverage_AA_table.csv*: how many times the mutation was sequenced

	  - *reference_AA_table.csv*: how many times the reference allele was seen in each sample

	  - *alternate_AA_table.csv*: how many times the alternate (mutant) allele was seen in each sample

Input File Structure
--------------------
A few directories are required for most operations.

	- **species_resources:** Contains information about the genome you targeted MIPs against.
	  Bound internally to :code:`/opt/species_resources`. Includes:

		- *fasta file*: Genome reference sequence in fasta format.

	  	- *bowtie2_genome*: The reference genome indexed using bowtie2.

  		- *bwa_genome*: The reference genome indexed using bwa.

  		- *SNPs*: VCF formatted locations of known SNPs in the reference genome.
		  Useful during MIP design to avoid targeting a polymorphic region of the genome.

		- *refgene*: RefGen style gene/gene prediction table in GenePred format.
  		  These are available at http://genome.ucsc.edu under Tools/Table Browser
		  for most species. If you have gff3/gtf formatted files, they can be
		  converted to GenePred format using Jim Kent's programs
		  `gff3ToGenePred <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred>`_
  		  and `gtfToGenePred <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred>`_.

  		- *refgene_tabix*: An index of the refgene file, created using tabix.

		- *file_locations.tsv*: This file is required for all operations. It is a
	  	  tab separated text file showing where each required file will be
	  	  located in the container. Each line corresponds to one file. First
	  	  field states the species for the file, second field states what kind of
	  	  file it is and the last field is the absolute path to the file within the
		  singularity container

	- **project_resources:** Contains information about your mip panel. Bound internally to 
	  :code:`/opt/project_resources`. Includes:
		- A targets.tsv file with the genomic coordinates of any protein-coding mutations
		  that are of particular interest.
		- a mip_ids folder that contains the mip arms of MIPs that target regions of the
		  genome that are of interest.
		- A few redundant json and csv files for easy access to MIP information. Our goal
		  is to remove these in the future.

Setting up your environment
===========================
Now that you know what steps will be performed and how files are organized, you
can set up your computer for analysis. This analysis can be done on any linux
computer that has singularity installed.

| You can obtain a copy of our latest sif file from here:
| https://baileylab.brown.edu/MIPTools/download/miptools_dev.sif
| This includes all executable programs needed for analysis


| In general, when you analyze any dataset, you should cd into a folder where
you'd like your analysis to go and run this command to download all relevant
scripts (plus a little text editor):
| :code:`singularity run -B $(pwd -P):/opt/config /path/to/your/downloaded/miptools_dev.sif`

Editing Settings
----------------
For convenience, settings can be passed in to all steps via a single shared
yaml file, called config.yaml. After downloading, open the file for editing
with a text editor and make sure to **follow the instructions in this file**,
editing it to contain the correct paths to your files (including project
resources, species resources, sample sheet, and sif files), as well as the
locations where you'd like the output to be sent.

|If the computer you're on is a remote computer, you can edit the config.yaml
 file using the text editor "micro" using this command:
| :code:`./micro config.yaml` 
| micro offers mouse support (so you can click a field of text to start editing
 it) and allows you to copy with ctrl-C, paste with ctrl-V, save with ctrl-S,
 and quit with ctrl-Q.

Analyzing an example dataset
----------------------------
To assist you, we've created a hypothetical dataset. This dataset contains 56
samples sequenced with 57 targeted genomic regions corresponding to known drug
resistance mutations of the P. falciparum genome. It assumes that MIPs have
already been designed, and that these MIPs have been used to target our regions
of interest, and that samples have already been pooled together, sequenced with
illumina paired end reads, and demultiplexed.

| You can download the tutorial dataset from here:
| https://baileylab.brown.edu/MIPTools/download/test-data.tar.gz
The dataset includes a project_resources folder, a species_resources folder, a
sample sheet, and a fastq directory with demultiplexed illumina paired end
reads as output.

| The downloaded tutorial dataset can be extracted with this command:
| :code:`tar -xvzf tutorial_dataset.tar.gz`

Wrangling
---------
See the "background" section above for what wrangling is and what output files
it produces.

| After editing the config.yaml file using the instructions in the yaml, you can
 execute the wrangler script with:
| :code:`bash wrangler_by_sample.sh`

Checking run stats
------------------
See the "background" section above for what check_run_stats is and what output
files it produces.

| After editing the config.yaml file using the instructions in the yaml, you can
 execute the check_run_stats script with:
| :code:`bash check_run_stats.sh`


Variant Calling
---------------
See the "background" section above for what variant calling is and what output
files it produces.

| After editing the config.yaml file using the instructions in the file, you
can execute the variant_calling script with:
| :code:`bash variant_calling.sh`

Resource Requirements
=====================
If you use the default processor counts, wrangling and variant calling should complete in approximately
five minutes each for the tutorial dataset, with checking run stats completing substantially faster.

More generally, resources required vary widely depending on the project. Wrangling and variant calling
require the most RAM and processing power, and both of these steps can be parallelized across multiple
processors. Wrangling with ~7,000 samples can take up to three days to complete, and some variant calling
steps on datasets this large can take a little over a week. The more processors (also known as CPUs or
threads) you ask for, the faster the job will run, the more RAM will be required, and the higher the
probability that the job will crash if RAM is insufficient. With a dataset containing 7,000 samples, a
single processor might require up to 150 GB of RAM in the variant calling step. Internally, MIPTools uses
snakemake so that if a job crashes partway through, you can rerun it and MIPTools will pick up where it
left off. Therefore, you might consider running a job once and requesting a large number of processors
(e.g. 15) so that most of the steps finish quickly. Then, if the job crashes, you might edit the settings
file to request fewer processors (e.g. 4 or even 2 or 1) so that any remaining particularly tricky steps
can be run with a lower likelihood of crashing.

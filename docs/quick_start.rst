====================
Abbreviated tutorial
====================
This abbreviated tutorial provides a minimal version of instructions for
people who want to get a "working version" up and running quickly to explore
the data for themselves. A little background is still needed. The tutorial is
designed to be run on a machine that has a local terminal, and requires a
little bit of familiarity with navigating files in unix, as well as access to a
machine that has singularity installed on it.

Background
==========
MIPTools supports several computational steps, including:
	- MIP design: This is for choosing resions of the genome that you'd like to
	  target.

	- Wrangling: This is for finding what haplotypes are associated with each
	  targeted region and the number of times each haplotype was seen in each
	  sample. This can be thought of as the "core" purpose of MIPTools. The output
	  is a tab delimited file called **allInfo.tsv.gz**

	- check_run_stats: After wrangling is finished, this is for figuring out
	  which samples and which targeted regions (aka MIPs) performed well and
	  which did not. There are several graphical outputs and comma separated
	  files produced at this stage. Among the most important are
	  **UMI_counts.csv** (how many times each targeted region of the genome was
	  seen in each sample) and **repool.csv** (which MIPs need to be
	  re-sequenced or repooled in each sample). **umi_heatmap.html** is
	  also useful for visualizing the performance of each MIP in each sample
	  (open this with a web browser).

	- variant calling: After wrangling is finished, this is for estimating the
	  number of times that each MIP was associated with a mutation in each
	  sample. If check_run_stats has not been run already, this step
	  automatically gathers statistics that form inputs for the variant calling
	  step. The outputs are a VCF file containing mutated genomic positions
	  for each sample and several tables that annotate mutations with
	  associated amino acid changes for each sample. For the tutorial, the main
	  outputs of interest are three tables that can be used to infer which
	  samples contain each mutation:

	  - *coverage_AA_table.csv*: how many times the mutation was sequenced
	  - *reference_AA_table.csv*: how many times the reference allele was seen
	    in each sample
	  - *alternate_AA_table.csv*: how many times the alternate (mutant) allele
	    was seen in each sample
    
    - prevalence calling: This is used to call the prevalence of mutations, and
      to graph these prevalences in a few different ways. This is performed
      using Jupyter notebooks (described below).

Input File Structure
--------------------
A few directories are required for most operations.

	- **species_resources:** Contains information about the genome you targeted MIPs against.
	  Includes:

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

	- **project_resources:** Contains information about your mip panel.
	  Includes:

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


In general, when you analyze any dataset, you should cd into a folder where
you'd like your analysis to go and run this command to download all relevant
scripts (plus a little text editor):

.. code-block:: console

    singularity run -B $(pwd -P):/opt/config /path/to/your/downloaded/miptools_dev.sif

Setting up a run
----------------

We've provided a script that allows you to edit your settings, wrangle your
data, generate run statistics, and perform variant calling, all from one
terminal interface. To launch this script, run the script that begins with
'run_miptools' using bash. An example (using version 0.5.0) is below.

.. code-block:: console

    bash run_miptools_v0.5.0.sh

For convenience, settings can be passed in to all steps via a single shared
yaml file, called config.yaml. After launching the run_miptools script, you can
edit the file by selecting option 1. 
**Carefully follow the instructions in this file**,
editing it to contain the correct paths to your files (including project
resources, species resources, sample sheet, and sif files, all described above),
as well as the locations where you'd like the output to be sent. The file
you're editing is called config and ends in .yaml. You can edit it in a
different text editor if you prefer, but the run_miptools script uses 'micro'.

When you're finished editing the file, type ctrl-S to save and then ctrl-Q to
quit.

run_miptools options
--------------------

After you've finished editing the config file, you can wrangle by selecting
option 2, or check run statistics on a finished wrangling job by selecting
option 3, or perform variant calling on a finished wrangling job by selecting
option 4. Wrangling, stat checking, and variant calling are all described in
the "background" section above. Option 5 allows you to launch a jupyter
notebook for an alternative interface for analyzing wrangled data (described
below). Option 6 is for unlocking snakemake. This is useful when a job crashes
partway through, in which case the program snakemake, which we use internally
for scheduling jobs, will lock the working directory. This option will unlock
the working directory again. Option 7 is used to quit out of the run_miptools
script.

.. _jupyter_instructions:

jupyter notebooks
-----------------
If you choose option 5 in the run_miptools script, this will launch a jupyter
notebook that you can access from your web browser. Some advanced commands are
available here, as well as detailed descriptions of each step associated with
statistics generation and variant calling. Jupyter notebooks are also used to
allow you to select mutations of interest and visualize prevalences.

If you're running this on a machine other than your physical machine (e.g. on a
server or a cluster that you're logging into with ssh) you'll need to scroll to
the beginning of the jupyter notebook launch message and copy the code into a
second terminal window **on your local machine**. The code you're looking for
should look something like this:

.. code-block:: console

    ssh -fNL localhost:####:###.###.###.##:#### your_username@server

This command allows you to make a connection between the remote server and your
local web browser. Because this code needs to be run from your local machine,
users accessing a server from a web browser (e.g. from open on demand browser
windows) rather than a local terminal window may struggle with this step.

After running this command (or ignoring it if running miptools on a local
machine) you can choose one of the http links to launch jupyter notebooks in
your web browser (depending on your system, it may be the second or third link
that works).

Once the browser launches, choose 'stats_and_variant_calling' and choose one of
the following notebooks:

	- check_run_stats.ipynb gives a detailed view of check_run_stats (described
	  above in the "background" section).
	- analysis-template-with-qual.ipynb gives a detailed view of the variant
	  calling step (also described above).
	- prevalence_plotting.ipynb is used to view the prevalences of
	  user-configurable mutations. These analyses are very interactive and
	  don't lend themselves easily to automation, so there is no automated
	  equivalent like there is for the other two steps.

Analyzing an example dataset
----------------------------
To assist you, we've created a hypothetical dataset. The dataset is in 'An
example cross-sectional study' on the lefthand menu. The dataset includes a
project_resources folder, a species_resources folder, a sample sheet, and a
fastq directory with demultiplexed illumina paired end reads. We've created a
walkthrough that guides you through editing your config file to analyze this
dataset.

Resource Requirements
=====================
If you use the default processor counts, the example dataset should complete in
approximately five minutes each for the wrangling and variant calling steps,
with checking run stats completing substantially faster.

More generally, resources required vary widely depending on the project.
Wrangling and variant calling require the most RAM and processing power, and
both of these steps can be parallelized across multiple processors. Each
processor requires some amount of RAM. The more processors (also known as CPUs
or threads) you ask for, the faster the job will run, the more RAM will be
required, and the higher the probability that the job will crash if RAM is
insufficient.

The resources required by wrangling and variant calling are primarily controlled
by the :code:`general_cpu_count` and :code:`freebayes_cpu_count` variables in
the config file. These variables control the number of processors to use in
parallel, for general steps of the pipeline and extremely memory-intensive steps
of variant calling, respectively.

As an example, wrangling with ~7,000 samples can take up to three days to
complete. For variant calling on a dataset like this, a single processor might
require up to 150 GB of RAM and two days. If there are five such steps and you
request one processor (freebayes_cpu_count set to 1), variant calling may take
ten days and 150 GB of RAM. If two processors are running simultaneously, you
might need only five days, but you would need to have 300 GB of RAM or more
available to avoid a crash.

Even on a large dataset, in addition to a few steps that require large amounts
of RAM, there are often hundreds of small steps that each require small amounts
of RAM, and these steps can be drastically sped up with parallelization.
Internally, MIPTools uses snakemake so that if variant calling or wrangling
crashes partway through, you can rerun it and MIPTools will pick up where it
left off. If a step crashes, snakemake will still attempt to complete all steps
that aren't dependent on the crashed step. Therefore, you might consider running
variant calling or wrangling once and requesting a large number of processors
(e.g. 15) so that most of the steps finish quickly. Then, if it crashes, you
might edit the config file to request fewer processors (e.g. 4 or even 2 or 
1) so that any remaining particularly tricky steps can be run with a lower
likelihood of crashing.
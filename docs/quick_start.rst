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
	  is a tab delimited file called allInfo.tsv.gz

	- stat-checking: This is for figuring out which samples and which targeted
	  regions (aka MIPs) performed well and which did not. There are several
	  graphical outputs and comma separated files produced at this stage. Among
	  the most important are barcode_counts.csv (how many times each targeted
	  region of the genome was seen in each sample) and repool.csv (which MIPs
	  need to be re-sequenced or repooled in each sample). umi_heatmap.html is
	  also useful for visualizing the performance of each MIP in each sample
	  (open this with a web browser).

	- variant calling: This is for estimating the number of times that each MIP
	  was associated with a mutation in each sample. The outputs are a VCF file
	  containing mutated genomic positions for each sample and several tables
	  that annotate mutations with associated amino acid changes for each sample.
	  For the tutorial, the main outputs of interest are three tables that can be
	  used to infer which samples contain each mutation:

	  - coverage_AA_table.csv: how many times the mutation was sequenced

	  - reference_AA_table.csv: how many times the reference allele was seen in each sample

	  - alternate_AA_table.csv: how many times the alternate (mutant) allele was seen in each sample

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

Analyzing the Data
==================
Now that we know what steps will be performed and how files are organized, we can look at a
hypothetical dataset. The tutorial dataset contains 56 samples sequenced with 57 targeted
genomic regions corresponding to known drug resistance mutations of the P. falciparum genome.
It assumes that MIPs have already been designed, and that samples have already been sequenced
using the MIPs using illumina paired end reads and demultiplexed.

| You can download the tutorial dataset from here:
| https://baileylab.brown.edu/MIPTools/download/test-data.tar.gz
| The dataset includes a project_resources folder, a species_resources folder, a sample sheet,
and a fastq directory with demultiplexed output.

| You can obtain a copy of our latest sif file from here:
| https://baileylab.brown.edu/MIPTools/download/miptools_dev.sif
| This includes all executable programs needed for analysis

Wrangling
---------
| For convenience, settings can be passed in to the wrangler via a yaml file. Later in the tutorial, we'll show some more advanced usage options available for troubleshooting and passing more customizable inputs. For now, you can obtain an example simple settings file for wrangling with this command:
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/wrangler_by_sample.yaml`
| After downloading, open the file for editing with a text editor and make sure to **follow the instructions in this file**, editing it to contain the correct path to the project resources, species resources, sample sheet, and sif files you downloaded above, as well as the
location where you'd like the output to be sent.

| We've also provided a bash script for converting the yaml settings into instructions for the wrangler. You can obtain the bash script for wrangling with this command (put it in the same folder as the settings yaml file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/wrangler_by_sample.sh`

| After changing directory to a folder that can run your data, you can execute the wrangler script with:
| :code:`bash wrangler_by_sample.sh`

Checking run stats
------------------
| After wrangling is finished, you can obtain a settings file for checking run stats with this command:
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/variant_calling.yaml`
| **Make sure to follow the instructions in this file.**

| You can obtain the script for checking run stats here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/check_run_stats.sh`

| And you can execute it like this:
| :code:`bash check_run_stats.sh`

Variant Calling
---------------
Variant calling uses the same settings file as check_run_stats.

| You can obtain the script for variant calling here (put it in the same folder as the settings file):
| :code:`wget https://github.com/bailey-lab/MIPTools/raw/master/user_scripts/variant_calling.sh`

| And you can execute it like this:
| :code:`bash variant_calling.sh`

Resource Requirements
=====================
Resources required vary widely depending on the project. Wrangling and variant calling require the
most RAM and processing power, and both of these steps can be parallelized across multiple processors.
The more processors (also known as CPUs or threads) you ask for, the faster the job will run, the more
RAM will be required, and the higher the probability that the job will crash. Internally, MIPTools uses
snakemake so that if a job crashes partway through, you can rerun it and MIPTools will pick up where it
left off. Therefore, you might consider running a job once, requesting a large number of processors (e.g.
15) so that most of the steps finish quickly, and then editing the settings file to request fewer
processors (e.g. 4 or even 2 or 1) if the job crashes so that any remaining particularly tricky steps can
be run with fewer processors with a lower likelihood of crashing.

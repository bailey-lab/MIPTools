===========
Quick Start
===========

Although :code:`miptools.sif` contains all programs needed, it does not include
the data to be analyzed or other resources to be used. Every time we run
Singularity we will **bind** needed directories to the container. There are
three resources directories which are required for most operations. In addition
to those, some apps need a :code:`data_dir` and :code:`analysis_dir`. The
:code:`-B` option is used for each binding:

.. code-block:: shell

	singularity some-command -B path_on_host:path_on_container path_to_container

The path on the left side of the colon specifies where on *your* computer the
directory is and the right side is the location in the container where the
directory should be bound (mounted) to. You should only change the left side of
the column according to the location of the resource you are providing, and
should *never* change the path on the right side. Each binding is specified
with a separate :code:`-B` option. See below for examples.

Directory Structure
===================

Three resource directories are required for most operations. These live outside
the container and must be **bound** to the container at run time with the
:code:`-B` option. In addition, a data directory and an analysis directory will
be used for most operations.

	- **base_resources:** Provided in the GitHub repository. It contains common
	  resources across projects. It should be bound to the container with
	  :code:`-B [path to base resources container]:/opt/resources`. This makes the
	  base_resources directory available to the container and it would be reached
	  at :code:`/opt/resources` path within the container. The
	  :code:`/opt/resources` part of this argument must not be altered. For
	  example, if my base resources are located in my computer at
	  :code:`/home/base`, I would bind it to the container with :code:`-B
	  /home/base:/opt/resources`.

	- **species_resources:** Contains resources shared by projects using the same
	  target species (Pf, human, etc.). Bind this to
	  :code:`/opt/species_resources` in the container. For example, if I am
	  working with *Plasmodium falciparum* sequences and I have the necessary
	  files in my computer at :code:`/home/pf3d/`, then the binding parameter is
	  :code:`-B /home/pf3d:/opt/species_resources`. The directory contains the
	  following contents:

	  	- *file_locations.tsv*: This file is required for all operations. It is a
	  	  tab separated text file showing where each required file will be
	  	  located in the container. Each line corresponds to one file. First
	  	  field states the species for the file, second field states what kind of
	  	  file it is and the last field is the absolute path to the file.

	  	  For example, the line 
	  	  *"pf fasta_genome /opt/species_resources/genomes/genome.fa"* would mean 
	  	  that the fasta genome file for the species 'pf' will be found at 
	  	  :code:`/opt/species_resources/genomes/genome.fa` within the container. 
	  	  This also means that there is a file at
	  	  :code:`/home/pf3d/genomes/genome.fa` in my computer, assuming I bound
	  	  :code:`/home/pf3d` to :code:`/opt/species_resources` in the container.

	  	- *fasta file*: This file is required for all operations. Genome
	  	  reference sequence in fasta format.

	  	- *bowtie2_genome*: This file is required for probe design operations
	  	  only. It is the reference genome indexed using bowtie2. If this is not
	  	  available, it can be generated using MIPTools.

  		- *bwa_genome*: This file is required for data analysis operations only.
  		  It is the reference genome indexed using bwa. If this is not available,
  		  it can be generated using MIPTools.

  		- *snps*: This is an optional file. However, it is extremely useful in
  		  probe designs to avoid probe arms landing on variant regions, etc. So
  		  it should always be used except in rare cases where such a file is not
  		  available for the target species. The format of the file is vcf.
  		  Individual genotypes are not necessary (a.k.a. sites only vcf). The
  		  only requirement is that the INFO field for each variant has a field
  		  showing the population allele frequency of alternate alleles. By
  		  default, AF field is used. The AF field lists the allele frequencies of
  		  each alternate allele, and does not list the frequency of the reference
  		  allele. Vcf files may have other INFO fields that include allele
  		  frequency information. If such a field is to be used, there are two
  		  settings in the design settings file (.rinfo file) that must be
  		  modified. *allele_frequency_name* field must be set to the INFO field
  		  name to be used; *af_start_index* may have to be set to a 1 (instead of
  		  default 0) depending on whether the reference allele frequency is
  		  provided in the new field. For example, if we want to use the 1000
  		  genomes vcf file, the allele frequencies are provided in the CAF field
  		  and they include the reference allele. We would have to change the
  		  *allele_frequency_name* field to *CAF* from the default *AF*; and set
  		  *af_start_index* to 1 because the first alternate allele's frequency is
  		  provided in the second place (following the reference allele).

  		- *refgene*: RefGen style gene/gene prediction table in GenePred format.
  		  These are available at http://genome.ucsc.edu under Tools/Table Browser
  		  for most species. The fields in the file are "bin, name, chrom, strand,
  		  txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds,
  		  score, name2, cdsStartStat, cdsEndStat, exonFrames". This file is
  		  required for probe design operations if genic information is to be
  		  used. For example, if probes need to be designed for exons of a gene,
  		  or a gene name is given as design target. If a gene name will be
  		  provided, it must match the **name2** column of the RefGen file. If you
  		  are creating this file manually, the only fields necessary are: chrom,
  		  strand, exonStarts, exonEnds and name2. All other fields can be set to
  		  an arbitrary value (none, for example) but not left empty. The order of
  		  columns must not be changed.

  		  .. note::

  		  	If you have gff3/gtf formatted files, they can be converted to
  		  	GenePred format using Jim Kent's programs `gff3ToGenePred
  		  	<http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred>`_
  		  	and `gtfToGenePred
  		  	<http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred>`_.

  		- *refgene_tabix*: RefGen file, sorted and indexed using tabix. File
  		  requirement is the same as the refgene file. tabix is available within
  		  the MIPTools container, so you don't have to install it yourself.


	- **project_resources:** Contains project specific files. Bind this to 
	  :code:`/opt/project_resources`.

	- **data_dir:** Contains data to be analyzed. Typically, nothing will be
	  written to this directory. Bind this directory to :code:`/opt/data`.

	- **analysis_dir:** Where analysis will be carried out and all output files
	  will be saved. Bind it to :code:`/opt/analysis` This is the only directory 
	  that needs write permission as the output will be saved here.

:code:`data_dir` and :code:`analysis_dir` will have different content for
different app operations. Also, one app's analysis directory may be the
next app's data directory in the pipeline.

Resource Requirements
=====================

Resources required vary widely depending on the project. Both designs and data
analysis can be parallelized, so the more CPUs you have the better. Plenty of
storage is also recommended. For designs on large target regions (>5kb), files
can take up 10 GB or more per region. Consider allocating > 5 GB RAM for a
large design region (multiply the RAM requirement by CPU number if
parallelizing). For a typical MIP data analysis involving ~1,000 MIPs and ~1,000
samples, consider using at least 20 CPUs and 20 GB RAM to get the analysis done
within 10-12 h. You should expect ~200 GB disk space used for such an analysis
as well, although most files can be removed after processing steps to reduce
long term disk usage.

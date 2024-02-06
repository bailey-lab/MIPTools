======================
MIPTools Documentation
======================

Welcome to the MIPTools User Guide!

MIPTools is a suite of computational tools that are used for molecular inversion
probe (MIP) design, data processing, and analysis. Throughout much of this
tutorial, we assume a user interested in using MIPs as a cost effective way to
amplify and sequence hundreds to thousands of targeted regions of the genomes
from hundreds to thousands of pooled barcoded samples. Our group primarily uses
these MIPs to assess relatedness and drug resistance status of Plasmodium
falciparum targets, but we have attempted to generalize this tool for other
questions and datasets. This toolset also assumes the use of unique molecular
identifiers (UMIs) to define unique MIP capture events.

A typical pipeline might look something like this:
 - First, a user might design MIP probes (using the probe design tool of this
 program) that have UMIs added to each MIP probe.
 - Second, a user might perform mip capturing reactions, PCR, sample barcoding,
 and illumina sequencing. The output data should be demultiplexed, resulting in
 two fastq files per sample. Bench techniques for these experiments are
 described elsewhere.
 - Third, the data is wrangled to generate an output file describing which
 genotypes (or haplotypes) are found at which abundances in each sample for each
 targeted region, using:
   - a sample sheet that describes the samples
   - a fastq folder of samples
   - a project resources folder that describes the probes
 - Finally, the haplotype data is analyzed using a variant caller (Freebayes is
 currently our best-supported tool) to produce a VCF file and some output tables
 with frequencies and prevalences of mutations of interest, using:
   - a sample sheet that describes the samples
   - a folder containing the wrangled haplotype data
   - a folder of indexed genomes for your species of interest
   - a project resources folder that describes the probes
   
.. toctree::
  :caption: Quick Start
  :maxdepth: 2
              
  installation
  quick_start

.. toctree::
  :caption: Guides
  :maxdepth: 2

  guides/probe-design
  guides/analysis-pipeline
  guides/hpcc

.. toctree::
  :caption: Man Pages
  :maxdepth: 1

  app-reference/download-app
  app-reference/download-superseded-app
  app-reference/demux-app
  app-reference/demux-qc-app
  app-reference/wrangler-app
  app-reference/jupyter-app

.. toctree::
  :caption: Links
  :maxdepth: 1
  :hidden:
  
  Source Code <https://github.com/bailey-lab/MIPTools>
  CHANGELOG
  license

Need Help?
----------

- Join our |Google Group|
- See our |issues page|

.. |Google Group| raw:: html

  <a href="https://groups.google.com/g/miptools" target="_blank">Google
  Group</a>
  
.. |issues page| raw:: html

  <a href="https://github.com/bailey-lab/MIPTools/issues"
  target="_blank">Github Issues Page</a>

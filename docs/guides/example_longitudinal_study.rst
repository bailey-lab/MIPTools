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
study. We also chose to include only MIPs covering Kelch13 and dhps.

The data for this file can be obtained here:
(Download link not built yet)

This data is organized into x main folders:
	- **pf_species_resources:** This folder includes an indexed copy of the
	  Pf3D7 falciparum genome, gene annotations, common SNPs, and a directory of
	  key files.

	- **k13_DHPS_project_resources:** This is a special project resources
	  folder that includes MIPs for K13 and DHPS only.

	- **metadata_files:** These include 2016_metadata.tsv, 2019_metadata.tsv,
	  and 2022_metadata.tsv, with collection site information for each sample.

	- **fastq_files:** These are the raw paired end sequencing reads associated
	  with each sample.

	- **sample_sheets:** This folder contains sample sheets for each year,
	  including 2016_sample_sheet.tsv, 2019_sample_sheet.tsv, and
	  2022_sample_sheet.tsv

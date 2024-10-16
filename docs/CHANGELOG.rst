=========
Changelog
=========

MIPTools (development version)
==============================

New Features
------------
- New features that have been introduced since the most recent release will
  appear here until they get rolled into the next release.

MIPTools 0.5.0
==============

New Features
------------
-  User can call a single script, run_miptools.sh, which automatically launches
   the other steps of miptools.
-  shell scripts and config files now have version numbers associated to reduce
   confusion around which sif files are compatible with which scripts.
   for interactively exploring mutation prevalences.
-  Added a choropleth map jupyter notebook, called prevalence_plotting.ipynb,
   for interactively exploring mutation prevalences.
-  Introduced a config.yaml file that eliminates the need for user to provide
   binding arguments to singularity
-  wrangling and variant calling are now both subdivided into snakemake steps so
   that partially completed jobs can be resumed and so that crashes report the
   exact step where a crash occurred. This also improves scalability to large
   datasets.
-  When running the ``wrangler`` app, if the number of UMIs detected for a MIP
   is above a certain threshold, we reduce the UMI count to a lower value. This
   is done in order to increase the speed of our pipeline. Above a certain UMI
   count, the information becomes redundant (:github:user:`arisp99`,
   :github:pull:`40`).
-  Add an additional argument to the ``wrangler`` app to control the population
   clustering fraction cutoff defined by :github:repo:`MIPWrangler
   <bailey-lab/MIPWrangler>` (:github:user:`arisp99`, :github:pull:`39`).
-  Add the capability to freeze software version numbers when building the
   container. Additionally, the version number for key software tools has been
   fixed (:github:user:`arisp99`, :github:pull:`32`).
-  Install :github:repo:`mipscripts <bailey-lab/mipscripts>`, which contains
   additional tools for analysis pipelines.
-  Perform additional argument parsing to ensure arguments are formatted
   correctly (:github:issue:`28`, :github:issue:`37`).
-  New ``download`` app supersedes the previous ``download`` app, which has
   been renamed to ``download_superseded``. The new app improves the method for
   downloading data from the Illumina BaseSpace Sequence Hub by using the
   official command line tool (:github:user:`arisp99`, :github:pull:`25`,
   :github:pull:`13`).

Bug Fixes
---------
-  Removed an alarming (but ultimately harmless) error message associated with
   snakemake 'shell' commands.
-  Removed the redundant mip_info.json files so they're no longer required by 
   the pipeline.
-  Upgraded pandas to a newer version and fixed ambiguous merge commands
-  Upgraded jupyter notebooks to a newer version
-  Upgraded conda packages.
-  Added cached wrangler and conda sif files to avoid automatically rebuilding
   these on MIPTools rebuild.
-  Upgrade C and C++ compiler versions (:github:issue:`43`).
-  Don't install conda and mamba packages using defaults as this can cause the
   install process to hang.
-  Upgrade ``libgfortran4`` to ``libgfortran5`` (:github:issue:`38`).
-  Let Freebayes run with only one CPU thread (:github:issue:`33`).
-  Fix error when app arguments have whitespace characters (:github:issue:`26`,
   :github:issue:`37`).
-  Fix missing file error when MIP arms file is created from the MIP
   info dictionary (:github:user:`aydemiro`, :github:pull:`23`).
-  Improve sample sheet preparation. Avoid errors when sample file
   columns are empty. Throw an error if there are invalid samples or
   input fields (:github:user:`aydemiro`, :github:pull:`22`).
-  Fix build failure due to dependency changes in the McCOILR R package
   (:github:issue:`7`).

Maintenance
-----------

-  Remove the ``msa2vcf`` program and other conversion tools
   (:github:issue:`35`).
-  Reduce size of image by deleting source code after installation of programs.
-  Remove sequence aligners (:github:issue:`35`).
-  Remove unused analysis settings files (:github:issue:`35`).
-  Install programs from GitHub instead of storing source code
   (:github:user:`arisp99`, :github:pull:`36`).
-  Update LICENSE year.
-  Store containers using an HTTP directory (:github:issue:`12`).
-  Remove duplicated files.
-  Improve bash errors.
-  Make strings human readable (:github:user:`arisp99`, :github:pull:`5`).

Documentation Overhaul
----------------------

-  Add guides on :ref:`probe design <guides/probe-design:Probe Design>`,
   :ref:`data analysis <guides/analysis-pipeline:Analysis Pipeline>`, and
   :ref:`HPCC use <guides/hpcc:HPCC Use>`.
-  Generate online documentation using
   `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`__ and
   `Github Pages <https://pages.github.com/>`__.
-  Improve app documentation.
-  Add doc-strings to python functions.
-  Improve clarity of README and add additional instructions on
   :ref:`downloading <installation:Quick Start>` or :ref:`building the
   container <installation:Install From Source>`.

MIPTools 0.4.0
==============================

-  Latest stable build.

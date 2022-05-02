=========
Changelog
=========

MIPTools (development version)
==============================

New Features
------------

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

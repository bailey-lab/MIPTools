=========
Changelog
=========

MIPTools (development version)
==============================

New Features
------------

-  Perform additional argument parsing to ensure arguments are formatted
   correctly (:github:issue:`28`).
-  New ``download`` app supersedes the previous ``download`` app, which has
   been renamed to ``download_superseded``. The new app improves the method for
   downloading data from the Illumina BaseSpace Sequence Hub by using the
   official command line tool (:github:user:`arisp99`, :github:pull:`25`,
   :github:pull:`13`).

Bug Fixes
---------

-  Let Freebayes run with only one CPU thread (:github:issue:`33`).
-  Fix error when app arguments have whitespace characters (:github:issue:`26`).
-  Fix missing file error when MIP arms file is created from the MIP
   info dictionary (:github:user:`aydemiro`, :github:pull:`23`).
-  Improve sample sheet preparation. Avoid errors when sample file
   columns are empty. Throw an error if there are invalid samples or
   input fields (:github:user:`aydemiro`, :github:pull:`22`).
-  Fix build failure due to dependency changes in the McCOILR R package
   (:github:issue:`7`).

Maintenance
-----------

-  Install programs from GitHub instead of storing source code
   (:github:user:`arisp99`, :github:pull:`36`).
-  Update LICENSE year.
-  Store containers using an HTTP directory (:github:issue:`12`).
-  Remove duplicated files.
-  Improve bash errors.
-  Make strings human readable (:github:user:`arisp99`, :github:pull:`5`).

Documentation Overhaul
----------------------

-  Generate online documentation using
   `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`__ and
   `Github Pages <https://pages.github.com/>`__.
-  Add better documentation for the ``jupyter`` app.
-  Add better documentation for the ``wrangler`` app.
-  Add better documentation for the ``download`` app.
-  Add better documentation for the ``demux_qc`` app.
-  Add doc-strings to python functions.
-  Improve clarity of README and add additional instructions on
   downloading or building the container.

MIPTools 0.4.0
==============================

-  Latest stable build.

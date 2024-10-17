.. _installation:

============
Installation
============
All the programs of MIPTools are packaged into a single file known as a
singularity image or sif file.

If you have a working copy of singularity, you can use our sif file without the
need to install anything else.


Obtaining MIPTools
==================
The MIPTools sif file, built and ready to use, can be downloaded from here:
https://baileylab.brown.edu/MIPTools/download/. This is the recommended approach
for obtaining MIPTools.

You can download the development version or any previous release:

.. code-block:: shell

    # Download the development version
    wget https://baileylab.brown.edu/MIPTools/download/miptools_dev.sif

    # Download the latest stable release
    wget https://baileylab.brown.edu/MIPTools/download/miptools_v0.5.0.sif

If your machine has singularity installed, you can ignore the rest of this
guide.

Sudo Privileges
===============

There are two scenarios in which you might need ``sudo`` privileges:
  - First, you will need ``sudo`` in order to install singularity if it's not
    already on your system.
  - Second, you will need ``sudo`` in order to build a sif file if you don't
    want to use our prebuilt sif files.

If you want to run the container on an environment where you don't have ``sudo``
privileges, either download a prebuilt image (recommended approach, see above)
or build the image on a machine where you *do* have ``sudo`` privileges and copy
the image file to the computer where you don't have ``sudo`` privileges.

Obtaining Singularity
=====================

Singularity is available for most Linux systems and is usually already installed
on academic high-performance clusters. Remember that you will need sudo
permissions in order to install singularity in a way that is compatible with our
pipelines. It is also possible to install and use on Mac OS using virtual
machines with a little bit of extra work. We have tested two flavors of
singularity. One is called apptainer and the other is singularity
community-edition (aka singularity-ce).

Apptainer
---------
Singularity can be installed from Apptainer using instructions
from here:
https://apptainer.org/docs/admin/main/installation.html

Specifically, if you're using this on Linux, we would recommend using pre-built
packages, which are available here under "assets":
https://github.com/apptainer/apptainer/releases

singularity-ce
--------------
Singularity-ce can be installed from here:
https://docs.sylabs.io/guides/4.2/admin-guide/installation.html

Specifically, if you're using this on Linux, we would recommend using pre-built
packages, which are available here under "assets":
https://github.com/sylabs/singularity/releases

.. _install-source:

Advanced: Building a sif file from source 
=========================================
The MIPTools sif file can also be built from source using the build.sh shell
script provided in the build folder of the
`GitHub repository <https://github.com/bailey-lab/MIPTools>`_. This could be
useful if, for example, you want to tweak the behavior of MIPTools by editing
the source code. This relies on a working copy of singularity and a machine with
``sudo`` privileges. You can install the most recent release using the
following:

Step 1, clone the repository (either a stable version or the dev version):

.. code-block:: shell

	# Clone stable version
	git clone --branch v0.5.0 https://github.com/bailey-lab/MIPTools.git

	# Clone dev version
	git clone https://github.com/bailey-lab/MIPTools.git

Step 2, make any changes to the source code.

Step 3, build the sif file:

.. code-block:: shell

    cd MIPTools/build
    bash build.sh miptools.sif

:code:`miptools.sif` is a single **portable** file which has all the programs
needed for MIP design, data analysis, and a lot more. Once the sif file is
built, you can copy or move the sif file anywhere (even onto a different
computer).

Demultiplexing
==============

If you plan to use MIPTools to demultiplex bcl files, you must download
:code:`bcl2fastq` separately. Currently, you can download it from `here
<https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html>`_.
You must download the file: :code:`bcl2fastq2 Conversion Software v2.20
Installer (Linux rpm)` and place it in the :code:`MIPTools/programs` directory.
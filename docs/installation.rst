============
Installation
============

Dependencies
============

A working copy of `Singularity <https://www.sylabs.io/docs/>`_ is required.
Singularity is best installed with **sudo**. While it is said to be possible to
install as an unprivileged user with some features missing, MIPTools has yet to
be tested on such an installation.

Singularity is available for most Linux systems, including high-performance
clusters. It is also possible to install and use on Mac OS using virtual
machines with a little bit of extra work.

Note that the :code:`snap` package is a rapid way to install the go language
required by Singularity (e.g., on Ubuntu/Debian: :code:`sudo snap install go
--classic`).

Quick Start
===========
The MIPTools container, built and ready to use, can be downloaded |download|.
You can download the development version or any previous release:

.. |download| raw:: html

  <a href="https://baileylab.brown.edu/MIPTools/download/"
  target="_blank">here</a>

.. code-block:: shell
	
	# Download the latest stable release
	wget https://baileylab.brown.edu/MIPTools/download/miptools_v0.4.0.sif

	# Download the development version
	wget https://baileylab.brown.edu/MIPTools/download/miptools_dev.sif

.. note::
	
	These prebuilt versions do not include the :code:`bcl2fastq` software due to
	its license. You must build the container yourself if you plan to use 
	MIPTools to demultiplex bcl files.

.. _install-source:

Install From Source 
===================
MIPTools can also be built from source using the definition file provided in
the `GitHub repository <https://github.com/bailey-lab/MIPTools>`_. You can
install the most recent release using the following:

.. code-block:: shell

	# Install stable version
	git clone --branch v0.4.0 https://github.com/bailey-lab/MIPTools.git


You can alternatively install the development version:

.. code-block:: shell

	# Install dev version
	git clone https://github.com/bailey-lab/MIPTools.git

Next, build the container, and you should be all set to get started using
MIPTools!

.. code-block:: shell

	cd MIPTools
	sudo singularity build miptools.sif MIPTools.def

:code:`miptools.sif` is a single **portable** file which has all the programs
needed for MIP design, data analysis, and a lot more.

Sudo Privileges
---------------

.. warning::

	You must have ``sudo`` privileges to *build* the image. You do not need
	``sudo`` to *use* the image.

If you want to run the container on an environment without ``sudo``, either
download a prebuilt image (see above) or build the container on your own
machine where you *do* have ``sudo`` privileges and copy the image file to the
computer without ``sudo``. Note that the Singularity program itself must have
been installed with ``sudo``.

Demultiplexing
--------------

If you plan to use MIPTools to demultiplex bcl files, you must download
:code:`bcl2fastq` separately. Currently, you can download it from `here
<https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html>`_.
You must download the file: :code:`bcl2fastq2 Conversion Software v2.20
Installer (Linux rpm)` and place it in the :code:`MIPTools/programs` directory.

CPU Usage
---------

The build process can take about 30-60 minutes to build, depending on the
number of CPU cores available. By default, the build process will use 20 CPU
cores. If the computer used for building the container has less then 20 CPU
cores available, change the :code:`CPU_COUNT=20` value at the top of the
:code:`MIPTools.def` file to a suitable number before building the container.
On the other hand, if the computer has additional CPU's, by all means, use them
by setting the same parameter to a higher value.

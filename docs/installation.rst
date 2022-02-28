============
Installation
============

Dependencies
============

A working copy of Singularity is required: https://www.sylabs.io/docs/.
Singularity is best installed with **sudo**. While it is said to be possible to
install with unprivileged user with some features missing, MIPTools hasn't been
tested on such an installation.

Singularity is available for most Linux systems including high-performance clusters. It is also possible to install
and use on Mac OS using virtual machines with a little bit of extra work.

Note that the :code:`snap` package is a rapid way to install the go language
required by Singularity (e.g. on Ubuntu/Debian: :code:`sudo snap install go
--classic`). Install system dependencies

Quick Start
===========
The MIPTools container, built and ready to use, can be
downloaded from the [Sylabs Cloud](https://cloud.sylabs.io/). You can download
either the development version or the most recent stable release:

.. code-block:: bash
	
	# Download the development version
	# The development version is updated every two weeks
	singularity pull library://apascha1/miptools/miptools:dev

	# Download the latest stable release
	singularity pull library://apascha1/miptools/miptools:v1.0.0


Note that these prebuilt versions do not include the :code:`bcl2fastq` software
due to its license. If you plan to use MIPTools to demultiplex bcl files, you
must build the container yourself.

.. _install-source:

Install From Source 
===================
MIPTools can also be built from source code using the definition file provided
in this [GitHub repository](https://github.com/bailey-lab/MIPTools).

The process can take about 10-30 minutes to build, depending on the number of
CPU cores available. By default, the build process will use 6 CPU cores. This
should pose no problems with most modern computers, but if the computer used
for building the container has less then 6 cpu cores available, change the
:code:`"CPU_COUNT=6"` value at the top of the :code:`MIPTools.def` file to a
suitable number before running the following code. On the other hand, if
you have access to more CPU power, by all means, use them by setting the
same parameter to a higher value.

You must have **sudo** privelege to _build_ the image. You do not need sudo to
_use_ the image. So if you want to run the container on an environment without
sudo, either download a prebuilt image (see above) or build the container on
your own machine where you _do_ have sudo privilege and copy the image file to
the computer without sudo. Note that the Singularity program itself must have
been installed with sudo.

If you plan to use MIPTools to demultiplex bcl files, you should download
:code:`bcl2fastq` separately. Currently, you can download it from
[here](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html),
but this may change in the future. You must download the file: :code:`bcl2fastq2
Conversion Software v2.20 Installer (Linux rpm)` and place it in the
:code:`MIPTools/programs` directory.

You can install the most recent release using the following:

.. code-block:: bash

	# Install stable version v1.0.0
	git clone --b v1.0.0 https://github.com/bailey-lab/MIPTools.git


You can alternatively install the development version:

.. code-block:: bash

	# Install dev version
	git clone https://github.com/bailey-lab/MIPTools.git

Next, simply build the container and you should be all set to get started using
MIPTools!

.. code-block:: bash

	cd MIPTools
	sudo singularity build miptools.sif MIPTools.def

:code:`miptools.sif` is a single **portable** file which has all the programs
needed for MIP design, data analysis, and a lot more. More information
about the extra programs and their uses will be added over time.

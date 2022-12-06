============
Installation
============

Dependencies
============

A working copy of `Singularity <https://www.sylabs.io/docs/>`_ is required.
Singularity is best installed with **sudo**. While it is said to be possible to
install with unprivileged user with some features missing, MIPTools hasn't been
tested on such an installation.

Singularity is available for most Linux systems including high-performance clusters. It is also possible to install
and use on Mac OS using virtual machines with a little bit of extra work.

Note that the :code:`snap` package is a rapid way to install the go language
required by Singularity (e.g. on Ubuntu/Debian: :code:`sudo snap install go
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
	its license. If you plan to use MIPTools to demultiplex bcl files, you must
	build the container yourself.

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

Next, simply build the container and you should be all set to get started using
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
machine where you *do* have ``sudo`` privilege and copy the image file to the
computer without ``sudo``. Note that the Singularity program itself must have
been installed with ``sudo``.

Software Versioning
-------------------

MIPTools installs several software tools together into the final built
container. Software packages are installed in the ``%post`` section on the
definition file, ``MIPTools.def`` (for more information of the definition file
consult the `Singularity documentation <https://sylabs.io/docs>`_). Programs in
MIPTools are installed in a variety of ways including via ``wget``,
``apt-get``, building source code for programs downloaded via ``git``, and even
via ``mamba``.

In order to ensure reproducible builds, the version number has been fixed for
many of the key programs MIPTools uses. The exceptions to this rule include
software installed via ``apt-get`` and ``mamba``. Software installed via
``mamba`` is defined in an ``environment.yml`` file in the root of the MIPTools
directory. This ``environment.yml`` file does not contain package versions as
in many cases dependency conflicts may arise. It is, however, possible to
specify the version number of installed packages by defining an
``environment_versioned.yml`` file in the root of the MIPTools directory.
During the build process if this file exists it will be used to install
``mamba`` packages. If no ``environment_versioned.yml`` file exists, it will be
generated during the build process and saved within the MIPTools container.
Users may then save this file to the root of the MIPTools directory to ensure
package versions of software installed with ``mamba`` do not change. To save
this file locally you may use ``singularity exec``:

.. code-block:: shell

	singularity exec <container> cat /opt/environment_versioned.yml > environment_versioned.yml


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

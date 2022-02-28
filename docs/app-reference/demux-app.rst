=====
demux
=====

Synopsis
========
.. code-block:: console
	
	singularity run [options] --app demux <container> [app_options]

Description
===========
Demultiplex data. Generate per-sample FASTQ files from the raw sequence data
consisting of BCL files.

.. warning::

	The demultiplexing software :code:`bcl2fastq` is not shipped with prebuilt
	versions of MIPTools. To demultiplex BCL files, you must `download
	<https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html>`_
	the software, place it in the :code:`programs` directory, and build MIPtools
	from :ref:`source <install-source>`.

Options
=======
.. code-block:: shell
	
	# Required
	-s    Path to the sample sheet for demultiplexing.

	# Optional
	-h    Print the help page.

Sample Sheet
------------

The sample sheet must be present in the directory mounted to
:code:`/opt/analysis`.

Examples
========

.. code-block:: shell

	# Set paths
	resource_dir=/bin/MIPTools/base_resources
	bcl_dir=/work/usr/download
	fastq_root_dir=/work/usr/
	container=/work/bin/MIPTools/miptools.sif

	# Run app
	singularity run \\
	  -B ${resource_dir}:/opt/resources \\
	  -B ${bcl_dir}:/opt/data \\
	  -B ${fastq_root_dir}:/opt/analysis \\
	  --app demux ${container} -s SampleSheet.csv
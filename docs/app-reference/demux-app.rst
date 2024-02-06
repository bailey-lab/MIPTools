.. _demux-app:

=====
demux
=====

Synopsis
========
.. code-block:: shell
	
	singularity run [run options...] --app demux <container> [app options...]

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
.. code-block:: none
	
	# Required
	-s    Path to the sample sheet for demultiplexing.
	# Optional
	-h    Print the help page.
	-p    The sequencing platform used. Can either be 'nextseq' or 'miseq'.
=======
Sample Sheet
------------

The sample sheet must be present in the directory mounted to
:code:`/opt/analysis`.

Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B bcl_dir:/opt/data \
	  -B fastq_root_dir:/opt/analysis \
	  --app demux miptools.sif -s SampleSheet.csv
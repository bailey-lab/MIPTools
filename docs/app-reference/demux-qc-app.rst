========
demux_qc
========

Synopsis
========
.. code-block:: console
	
	singularity run [options] --app demux_qc <container> [app_options]

Description
===========
Compute quality control checks on demultiplexed data. Prints the total number
of sequencing reads and how many were of undetermined indices, meaning that
they contained indices that were not present in the sample file. Additionally
prints how many of the undetermined index reads belong to possible primer
pairs, i.e., they are likely due to errors in the provided sample list and not
just faulty reads from the sequencer.

Options
=======
.. code-block:: shell
	
	# Required
	-s    The sequencing platform used. Can either be 'nextseq' or 'miseq'.

	# Optional
	-h    Print the help page.

Examples
========

.. code-block:: shell

	# Set paths
	resource_dir=/bin/MIPTools/base_resources
	fastq_root=/work/usr/download
	container=/work/bin/MIPTools/miptools.sif

	# Run app
	singularity run \\
	  -B ${resource_dir}:/opt/resources \\
	  -B ${fastq_root}:/opt/analysis \\
	  --app demux_qc ${container} -p nextseq

	singularity run \\
	  -B ${resource_dir}:/opt/resources \\
	  -B ${fastq_root}:/opt/analysis \\
	  --app demux_qc ${container} -p miseq

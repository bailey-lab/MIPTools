========
demux_qc
========

Synopsis
========
.. code-block:: shell
	
	singularity run [run options...] --app demux_qc <container> [app options...]

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
.. code-block:: none
	
	# Required
	-p    The sequencing platform used. Can either be 'nextseq' or 'miseq'.
	# Optional
	-h    Print the help page.

Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B downloaded:/opt/analysis \
	  --app demux_qc miptools.sif -p 'nextseq'

	singularity run \
	  -B base_resources:/opt/resources \
	  -B downloaded:/opt/analysis \
	  --app demux_qc miptools.sif -p 'miseq'

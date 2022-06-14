.. _jupyter-app:

=======
jupyter
=======

Synopsis
========
.. code-block:: shell
	
	singularity run [run options...] --app jupyter <container> [app options...]

Description
===========
Open an interactive Jupyter Notebook. The notebook can be used for
post-wrangler mapping and variant calling.

Options
=======
.. code-block:: none
	
	# Optional
	-d    The notebook directory.
	-h    Print the help page.
	-p    The port to be used to load the Jupyter Notebook.

Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B project_resources:/opt/project_resources \
	  -B pf_species_resources:/opt/species_resources \
	  -B wrangler_dir:/opt/data \
	  -B variant_dir:/opt/analysis \
	  --app jupyter miptools.sif

Accessing the Notebooks
=======================
When running this app, a series of instructions will be pasted to the terminal
on how to access the notebook. Accessing the server will depend on whether the
app was run via a remote server or via your local machine. In the case where
the app was run via a remote server, you must forward the port to your local
machine. At this point you may open a browser, paste the URL, and navigate to
the :code:`analysis` directory to access the notebooks.

.. _jupyter-app-faq:

FAQ
===

I can't connect to the remote server via ssh.
	Make sure that you can access the server you are trying to connect to. In
	some cases, you may need to use a VPN.

How can I use this app on a HPCC?
	If you are using a cluster computing system, the login node will likely be
	different from the compute node. In that case, the printed remote server will
	refer to the compute node and you will need to change this to the login node.
	For example, you may login in to the HPCC using: :code:`ssh
	oa@login.hpcc.edu` and submit a job to the compute node :code:`compute1`. In
	this case, the app may print :code:`ssh -N -f -L
	localhost:9913:128.148.254.107:9913 oa@compute1.hpcc.edu`. You would change
	this to :code:`ssh -N -f -L localhost:9913:128.148.254.107:9913
	oa@login.hpcc.edu`.

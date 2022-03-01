=======
jupyter
=======

Synopsis
========
.. code-block:: console
	
	singularity run [options] --app jupyter <container> [app_options]

Description
===========
Open an interactive Jupyter Notebook. The notebook can be used for
post-wrangler mapping and variant calling.

Options
=======
.. code-block:: shell
	
	# Optional
	-d    The port to be used to load the Jupyter Notebook.
	-h    Print the help page.
	-p    The notebook directory.

Examples
========

.. code-block:: shell

	# Set paths
	resource_dir=/bin/MIPTools/base_resources
	project_resources=/work/usr/DR1_project_resources
	species_resources=/work/usr/pf_species_resources
	wrangler_dir=/work/usr/wrangler
	variant_dir=/work/usr/variant
	container=/work/bin/MIPTools/miptools.sif

	# Run app
	singularity run \\
	  -B ${resource_dir}:/opt/resources \\
	  -B ${project_resources}:/opt/project_resources \\
	  -B ${species_resources}:/opt/species_resources \\
	  -B ${wrangler_dir}:/opt/data \\
		-B ${variant_dir}:/opt/analysis \\
	  --app jupyter ${container}

Accessing the Notebooks
=======================
When running this app, a series of instructions will be pasted to the terminal
on how to access the notebook. Accessing the server will depend on whether
the app was run via a remote server or via your local machine. In the case where
the app was run via a remote server, you must forward the port to your local
machine. At this point you may open a browser and paste the URL to access the
notebooks.

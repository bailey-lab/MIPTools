=================
Analysis Pipeline
=================

This guide provides a walkthrough of some of the key analysis steps in
analyzing molecular inversion probe data. We use a test data set composed of
FASTQ files and a sample list.

.. note:: 
	
	This guide does not cover the :ref:`download-app` or :ref:`demux-app` apps.
	The test data provided has already been demultiplexed.

Download Test Data
==================

The test data set can be downloaded |download| or via the use of the command
line:

.. |download| raw:: html

  <a href="https://baileylab.brown.edu/MIPTools/download/test-data/"
  target="_blank">here</a>

.. code-block:: shell
	
	# Download and untar directory
	wget -qO- https://baileylab.brown.edu/MIPTools/download/test-data.tar.gz | tar -xvz

The test data set contains 5 directories that contain the test data, species
resources, as well as project resources:

.. code-block:: shell

  tree -FL 1 test-data

  #> test-data
  #> ├── DR1_project_resources/
  #> ├── hg38_host/
  #> ├── pf_species_resources/
  #> ├── test_data/
  #> └── test_design_resources_pf/

Wrangle Data
============

After downloading the test data, the next step is to run the :ref:`wrangler app
<wrangler-app>`. We first create a directory that will store our wrangler
analysis and copy our sample list into said directory:

.. code-block:: shell
	
	mkdir test-data/wrangler
	cp test-data/test_data/sample_list.tsv test-data/wrangler/

We additionally define several parameters needed to wrangle data:

.. code-block:: shell

	experiment_id='test_run'
	sample_list='sample_list.tsv'
	probe_sets_used='DR1,VAR4'
	sample_sets_used='JJJ'
	cpu_number=10
	min_capture_length=30

Next, we can run the :ref:`wrangler app <wrangler-app>`. For additional
instructions on what each flag represents, consult the :ref:`man page
<wrangler-app>` for the app or the built in documentation with
:code:`singularity run --app wrangler miptools_dev.sif -h`.

.. code-block:: shell

  singularity run \
    -B test-data/DR1_project_resources:/opt/project_resources \
    -B test-data/test_data/fastq:/opt/data \
    -B test-data/wrangler:/opt/analysis \
    --app wrangler miptools_dev.sif \
    -e ${experiment_id} -l ${sample_list} -p ${probe_sets_used} \
    -s ${sample_sets_used} -c ${cpu_number} -m ${min_capture_length}

The :ref:`wrangler app <wrangler-app>` will save the main outputs as compressed
files in the :code:`wrangler` directory. There will additionally be a
:code:`nohup` file that contains errors and warning messages logged by the
:ref:`wrangler app <wrangler-app>`. This file should be empty if the all went
well. In our example run, the :code:`nohup` file was empty and the main
outputs were aggregated into the three files:

* :code:`run_test_run_wrangled_20220314.txt.gz`
* :code:`extractInfoByTarget.txt.gz`
* :code:`extractInfoSummary.txt.gz`

.. tip::

	After confirming the :ref:`wrangler app <wrangler-app>` successfully ran, we
	recommend you delete the :code:`wrangler/analysis` directory. This will
	remove many small files and save space in the future.

	.. code-block::

		rm -rf test-data/wrangler/analysis

Variant Calling
===============

To further process our data and call and analyze variants, we will leverage an
interactive `Jupyter notebook <https://jupyter.org/>`_ by calling the
:ref:`jupyter app <jupyter-app>`. Our main variant calling method uses the
`Freebayes software <https://arxiv.org/abs/1207.3907>`_, a Bayesian genetic
variant detector. While we have optimized the algorithm for calling on
molecular inversion probes (MIPs), we use an interactive environment for
calling and initial assessment to provide the user with greater
customizability.

Before running the :ref:`jupyter app <jupyter-app>`, we must define a new
directory in which we will run our variant calling pipeline:

.. code-block:: shell

	mkdir test-data/variant

Then we can start our Jupyter notebook:

.. code-block:: shell

  singularity run \
    -B test-data/DR1_project_resources:/opt/project_resources \
    -B test-data/pf_species_resources:/opt/species_resources \
    -B test-data/wrangler:/opt/data \
    -B test-data/variant:/opt/analysis \
    --app jupyter miptools_dev.sif

A series of instructions will be printed to the terminal on how to access the
notebook. Follow these instructions to run the Jupyter notebooks in a web
browser. For more information refer to the :ref:`FAQ of the jupyter app
<jupyter-app-faq>`. Next, navigate to the :code:`analysis` directory. The
:code:`analysis-template-with-qual` notebook contains a demonstration of 
processing data, variant calling, and additional data analysis.

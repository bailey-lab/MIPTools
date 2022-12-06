============
Probe Design
============

This guide provides a walkthrough of how to design molecular inversion probes
(MIPs) using MIPTools. We use a test data set to demonstrate each step in the
design process.

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

Test designs will be carried out on a few *Plasmodium falciparum* targets. The
species resources directory includes resources for the *P. falciparum*
genome and the host resources directory contain the resources for the human
genome, which will be used as the host species.

.. attention::

	The MIP design manual can be found `here
	<https://docs.google.com/document/d/1k3SpO8B5zz6OVTn1wgxivep2GbUlCENunSwNap0cV2c/>`_,
	for the time being, and an example walkthrough `here
	<https://docs.google.com/document/d/1sY8EIbiWy_cW9TNc7jbM__amc1i_ne8KTzqhz69QQNU/>`_.
	We plan on merging those documents here in the future...

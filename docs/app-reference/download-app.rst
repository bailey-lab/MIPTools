.. _download-app:

========
download
========

Synopsis
========
.. code-block:: shell
	
	singularity run [run options...] --app download <container> [app options...]

Description
===========
Download data from the Illumina BaseSpace Sequence Hub.

Options
=======
.. code-block:: none
	
	# Required
	-r    The run ID of the data to download.

	# Optional
	-h    Print the help page.

Download Destination
--------------------
Data will be downloaded to :code:`/opt/analysis`. A directory may be mounted
to this path to customize the download destination.

Authentication Credential File
------------------------------

.. note::
	
	Users must first authenticate their account in order to download data from
	the BaseSpace Sequence Hub.

An authentication token must be generated in order to download data from the
BaseSpace Sequence Hub. The steps to do so are outlined below:

#. Go to `<https://developer.basespace.illumina.com/>`_ and log in.
#. From the toolbar, select My Apps.
#. In the applications tab, select Create a New Application.
#. Fill out the Applications Details and then select Create Application.
#. On the application page, select the Credentials tab and copy the Access Token.

Once generated, the authentication token must be copied to
:code:`base_resources/access-token.txt`.

Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B downloaded:/opt/analysis \
	  --app download miptools.sif -r 12345

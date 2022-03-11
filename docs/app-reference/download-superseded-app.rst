===================
download_superseded
===================

Synopsis
========
.. code-block:: console
	
	singularity run [options] --app download_superseded <container> [app_options]

Description
===========
Download data from the Illumina BaseSpace Sequence Hub.

.. warning:: 
	
	This app has been superseded by the :ref:`download app <download-app>`, which
	uses the BaseSpace CLI for downloading data.

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
BaseSpace Sequence Hub. In order to do so, you consult the :ref:`authentication
credential file section of the download app <authenticate-label>`. Once
generated authentication token must be copied to
:code:`base_resources/access-token.txt`.

.. note::
	
	The :ref:`authentication credential file section of the download app
	<authenticate-label>` will generate a configuration file that indicates an API
	server to contact and an access token to authenticate against BaseSpace
	Sequence Hub. Only the **access token value** must be copied to
	:code:`base_resources/access-token.txt`.


Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B downloaded:/opt/analysis \
	  --app download_superseded miptools.sif -r 12345

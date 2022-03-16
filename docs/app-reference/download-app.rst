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
	-i    The run ID of the data to download.

	# Optional
	-o    The path to the output directory.
	-c    The path to the authentication credentials file.
	-h    Print the help page.

Defaults
--------
.. code-block:: shell

	-o    Default: '/opt/analysis'
	-c    Default: '/opt/resources/basespace.cfg'

.. _authenticate-label:

Authentication Credential File
------------------------------

.. note::
	
	Users must first authenticate their account in order to download data from
	the BaseSpace Sequence Hub.

The authentication credential file can be generated via the `BaseSpace CLI
<https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview#Authenticate>`_.
The :code:`bs auth` command generates a configuration file that indicates an
API server to contact and an access token to authenticate against BaseSpace
Sequence Hub. The default API server defaults to Virginia, USA. Other options
include Europe and the UK, among others. To override the default API server, 
use the :code:`--api-server` option.

.. code-block:: shell
	
	singularity exec [options] <container> bs auth --force

	# Specify UK API server
	singularity exec [options] <container> bs auth --force --api-server https://api.euw2.sh.basespace.illumina.com

These commands will generate a configuration file: :code:`basespace.cfg` in
:code:`${HOME}/.basespace`. In order for the download app to access the configuration file, you must copy the file into :code:`base_resources`.

.. code-block:: shell

	cp ${HOME}/.basespace/basespace.cfg base_resources


Examples
========

.. code-block:: shell

	singularity run \
	  -B base_resources:/opt/resources \
	  -B downloaded:/opt/analysis \
	  --app download miptools.sif -i 12345

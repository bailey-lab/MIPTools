========
HPCC Use
========

In order to improve the computational efficiency of MIPTools, it is possible to
use MIPTools via a high-performance computing cluster (HPCC). While we leave
the details of HPCC use to the reader, briefly, clusters use a management and
job scheduling system to organize user requests. Users may submit jobs to these
management systems, which will in turn schedule and execute submitted jobs.

.. note::
	This guide covers HPCCs configured to use `slurm
	<https://slurm.schedmd.com/overview.html>`_ cluster manager. For HPCCs
	configured with another job manager, please review the documentation for the
	job manager of interest.

To submit a job via slurm, users can run the following:

.. code-block:: shell
	
	sbatch <jobscript>

Each job script is a batch script with commands for the HPCC to execute. These
files also contain configuration commands used by slurm. The bash shebang and
these configuration commands make up what we call the header section of the
script. An example is shown below:

.. literalinclude:: example-jobscript.sh
	:language: shell
	:lines: 1-14
	:lineno-match:

The next section of the job script contains the bash code to be executed by the
HPCC. This can contain a command to run a MIPTools app or some other MITPools
command. For example, users could run the :ref:`wrangler app <wrangler-app>`:

.. literalinclude:: example-jobscript.sh
	:language: shell
	:lines: 16-
	:lineno-match:

This is just one of the many possible job scripts a user could use. Almost
every aspect of the MIPTools pipeline could be configured to run via a HPCC.
For more examples of job scripts, see the :github:repo:`slurm-scripts
repository <bailey-lab/slurm-scripts>`. The MIPTools folder has scripts to run
different MIPTools commands. Please feel free to add more scripts via a pull
request!

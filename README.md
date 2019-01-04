MIP
=========

## Overview

add 3-4 sentence to provide the grand overview or is that elsewhere on github???

## Installation
Clone the github repository:
```bash
git clone git@github.com:bailey-lab/MIPMaker.git
```
Note: This currently requires access to the baileylab private repository
### Dependencies
Requires a working copy of Singularity: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps

Note: Snap package install is a rapid way to install the go language required by Singularity (e.g. on Ubuntu/Debian: `sudo snap install go --classic`) 

### Build MIPMaker from the definition file 
This takes  about 10 minutes and may appear stalled. 

User must have **sudo** privilege to _build_ the container image. You do not need sudo to _use_ the image. So if you want to run the container on an environment without sudo, build the container on your own machine and copy the image file to the non-sudo host machine, where singularity is already installed (Note: the pipeline has not been tested on a host with locally/user installed Singularity.)
```bash
cd MIPMaker
sudo singularity build mipmaker.sif MIPMaker.def
```
The created `mipmaker.sif` is a single **portable** container which has all the programs needed for MIP data analysis and a lot more.
More information about its uses will be added over time.

##Usage

### General Overview and Data Structure
The computer on which the **container** is running on is called the **host**. So if you're running this on your laptop, your laptop is the host. If you're connected to a HPC cluster and running the container there, the HPC is the **host**.  

In general, 3 external input host directories are required and need to be mapped to specific directories in the container. 
- Base (general) resource directory is included in the git repository directory `MIPMaker/base_resources`and needs to map to `/opt/resources` in  container
- Project specific resource directory must map to `/opt/project_resources` in container
- Analysis/working directory  must map to `/opt/work` in container 
    - this directory must have executing user write permissions
    - must contain **sample_list.tsv** file (see test_data for an example).
The -B option is used to map host directories to container directories (e.g `-B /hostpath/dir:/containerpath/dir`)

#### Pre-wrangler usage
Set up a MIPWrangler workspace and scripts to run MIPWrangler. Although currently the container does not include MIPWrangler itself, this  first step is necessary to generate some files used in "post wrangler" analysis. User needs to specify the sequencing **sequencing_platform** used (miseq or nextseq) and **experiment_id** (we use sequencing date in the format "YYMMDD", but it can be any name for the current experiment.
```bash
singularity run --app wrangler \
    -B [[user host general base resources]]:/opt/resources \
    -B [[user host project specific resources]]:/opt/project_resources \
    -B [[host working analysis dir]]:/opt/work \
    mipmaker.sif sequencing_platform experiment_id
```
This should generate a new directory in your "analysis_dir" named "experiment_id" that contains a few files to use for MIPWrangler program.

#### Post-wrangler usage
**analysis_dir** must have 2 more files in addition to the sample_list.tsv for post-wrangler analysis.
- **wrangler_output**: MIPWrangler output file (can be gzipped or regular text file).
- **settings.txt**: analysis settings file

Run the following command, arguments (-B) to fit your **host** directory structure.
```bash
singularity run --app jupyter \
    -B [[on-host base resources dir]]:/opt/resources \
    -B [[on-host project_resources dir]]:/opt/project_resources \
    -B [[on-host work analysis dir ]]:/opt/work \
    mipmaker.sif
```

This starts a jupyter notebook on the host computer and generates a unique token everytime that must be pasted into a browser: 
```bash
[I 09:31:00.231 NotebookApp] Serving notebooks from local directory: /opt
[I 09:31:00.231 NotebookApp] The Jupyter Notebook is running at:
[I 09:31:00.231 NotebookApp] http://localhost:8888/?token=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
[I 09:31:00.231 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 09:31:00.232 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:8888/?token=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

If the host is a remote server, e.g. HPC, forward the notebook port to your computer. Skip this if the host is a computer where you have access to the GUI, specifically a web browser. Make sure the port number matches the one you see in the notebook address above (8888 for this example).

```bash
ssh -N -f -L localhost:8888:localhost:8888 username@serveraddress
```

Finally, open a web browser (Chrome, Chromium, Firefox works in my experience), and then copy and paste the specified link from your terminal (http://localhost:8888/?token=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)  


To excute the notebook script and perform analysis:
- Navigate to `/opt/resources/analysis.ipyn`
- Make a copy of this file and move to the analysis directory
- Execute this copy `/opt/work/analysis.ipyn`
- Follow the instructions contained in the notebook and modify as needed for specific data set. (You can step through each chunk of code or execute it all at once.)

### Test run dataset for post wrangler analysis
An example data set containing MIPWrangler output from P. falciparum MIPs is available for testing: 
- a test_data folder is included in  git repo `MIPmaker/test_data` which is the working analysis directory to map to `/opt/work`
- the plasmodium_resources (not in the repo) as project_resources must be download separately from dropbox to map to `/opt/project_resources` 
    - Dropbox Link: https://www.dropbox.com/sh/pe63kb9wnx6j4ks/AAAsV_xz06jse-XFahOiNoSAa?dl=0 

After copy and running `analysis.ipyn`, you can step through each chunk of code or execute it all at once.  To ensure that it is working properly, you can then compare your results to the results in the original analysis.pynb.  



mipster
=========
## Installation
Clone the github repository:
```bash
git clone git@github.com:bailey-lab/MIPMaker.git
```
### Dependencies
Requires a working copy of Singularity: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
### Build MIPMaker from the definition file 
```bash
# go to MIPMaker directory
cd MIPMaker
# build container "mipmaker.sif". Container name and path can be different.
# this takes  about 10 minutes.
# user must have root privilege to build the container, but not to use it. So if you want to run the container on an environment without sudo, build the container on your own machine and copy the container to the host machine.
sudo singularity build mipmaker.sif MIPMaker.def
#
```
### Usage
The computer on which the container is running is called the **host**. So if you're running this on your laptop, your laptop is the host. If you're connected to a HPC cluster and running the container there, the HPC is the **host**. 
In addition to base_resources included in the repository, each project requires a project_resources directory with project specific files in it; and an analysis directory where MIP analysis will be carried out. These 3 directories must be on the host computer. This README refers to these 3 directories as base_resources, project_resources, and analysis; these can be paths to any directory on the host where you have +w permission.
Analysis directory must have sample_list.tsv file (see test_data for an example).
#### Pre-wrangler
Set up a MIPWrangler workspace and scripts to run MIPWrangler. Although currently the container does not include MIPWrangler itself, this  first step is necessary to generate some files used in "post wrangler" analysis. User needs to specify the sequencing **sequencing_platform** used (miseq or nextseq) and **experiment_id** (we use sequencing date in the format "YYMMDD", but it can be any name for the current experiment.
```bash
singularity run --app wrangler -B base_resources(on-host):/opt/resources -B project_resources(on-host):/opt/project_resources -B analysis_dir(on-host):/opt/work mipmaker.sif sequencing_platform experiment_id
```
This should generate a new directory in your "analysis_dir" named "experiment_id" that contains a few files to use for MIPWrangler program.
#### Post-wrangler
Analysis directory must have 2 more files relevant in addition to the sample_list.tsv.
wrangler_output: MIPWrangler output file (can be gzipped or regular text file).
settings.txt: analysis settings file


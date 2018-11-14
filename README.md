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
1. First usage of the container is to set up a MIPWrangler workspace and scripts to run MIPWrangler,  although currently the container does not include MIPWrangler itself. (This will be added soon)
```bash
# "shell" into the container
singularity shell -B base_resources(on-host):/opt/resources -B project_resources(on-host):/opt/project_resources -B analysis_dir(on-host):/opt/work mipmaker.sif
```
This should take you into the container. You'll notice the prompt has changed to "Singularity mipmaker.sif:~>  " like below.
```bash
Singularity mipmaker.sif:~>  

```

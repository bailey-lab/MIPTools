MIP
=========
## Installation
Clone the github repository:
```bash
git clone git@github.com:bailey-lab/MIPMaker.git
```
### Dependencies
Requires a working copy of Singularity: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
### Build MIPMaker from the definition file 
This takes  about 10 minutes.  
User must have **sudo** privilege to _build_ the image. You do not need sudo to _use_ the image. So if you want to run the container on an environment without sudo, build the container on your own machine and copy the image file to the host machine.
```bash
cd MIPMaker
sudo singularity build mipmaker.sif MIPMaker.def
```
mipmaker.sif is a single **portable** file which has all the programs needed for MIP data analysis and a lot more.  
More information about its uses will be added over time.
### Directory Structure
Following directories are needed for MIP data analysis.
*  **base_resources:** Included in the repository, contains common resources across projects.
*  **species_resources:** Contains resources shared by projects using the same target species (Pf, human, etc.)
*  **project_resources:** Contains project specific files.
*  **data_dir:** Contains data to be analyzed (see below).
*  **analysis_dir:** Where analysis will be carried out and all output files will be saved.
### Usage
The computer on which the container is running on is called the **host**. So if you're running this on your laptop, your laptop is the host. If you're connected to a HPC cluster and running the container there, the HPC is the **host**.  

The directories  listed above must be on the host computer. They can be paths to any directory on the host where you have +w permission.  
#### Pre-wrangler usage
Set up a MIPWrangler workspace and scripts to run MIPWrangler. Although currently the container does not include MIPWrangler itself, this  first step is necessary to generate some files used in "post wrangler" analysis. User needs to specify the **sequencing_platform** used (miseq or nextseq) and **experiment_id** (we use sequencing date in the format "YYMMDD", but it can be any name for the current experiment. The **data** directory must have **sample_list.tsv** file (see test_data for an example).  

Run the following command, changing only the **host** part of the binding arguments (-B) to fit your directory structure.
```bash
singularity run --app wrangler \
    -B base_resources(on-host):/opt/resources \
    -B project_resources(on-host):/opt/project_resources \
    -B data_dir(on-host):/opt/data \
    -B analysis_dir(on-host):/opt/analysis \
    mipmaker.sif -p sequencing_platform -e experiment_id -l /opt/data/samaple_list.tsv
```
This should generate a new directory in your "data_dir" named "experiment_id" that contains a few files to use for MIPWrangler program.
#### Post-wrangler usage
**data_dir** must have 1 more file in addition to the sample_list.tsv for post-wrangler analysis.  
*  **wrangler_output**: MIPWrangler output file (can be gzipped or regular text file).
**analysis_dir** must contain the analysis settings file.
*  **settings.txt**: analysis settings file  

Run the following command, changing only the **host** part of the binding arguments (-B) to fit your directory structure. Two optional arguments can set the notebook server port (default 8888) and notebook directory where the server is started (default /opt, **this is relative to the container and NOT the host**, so do not change from the default if you're not absolutely sure it is called for).  

Note that we are binding an additional directory now (species_resources) that we did not need for the *wrangler* app.
```bash
singularity run --app jupyter \
    -B base_resources(on-host):/opt/resources \
    -B species_resources(on-host):/opt/species_resources \
    -B project_resources(on-host):/opt/project_resources \
    -B data_dir(on-host):/opt/data \
    -B analysis_dir(on-host):/opt/analysis \
    mipmaker.sif -p port_to_use -d notebook_directory
```
**port_to_use** and **notebook_directory** are optional arguments defaulting to **8888** and **/opt**, respectively.
This starts a jupyter notebook on the host computer: 
```bash
[I 09:31:00.231 NotebookApp] Serving notebooks from local directory: /opt
[I 09:31:00.231 NotebookApp] The Jupyter Notebook is running at:
[I 09:31:00.231 NotebookApp] http://localhost:8888/?token=600be64ca79ef74ebe4fbd3698bc8a0d049e01d4e28b30ec
[I 09:31:00.231 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 09:31:00.232 NotebookApp] 
    
    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:8888/?token=600be64ca79ef74ebe4fbd3698bc8a0d049e01d4e28b30ec
```
If the host is a remote server, e.g. HPC, forward the notebook port to your computer. Skip this if the host is a computer where you have access to the GUI, specifically a web browser. Make sure the port number matches the one you see in the notebook address above (8888 for this example).
```bash
ssh -N -f -L localhost:8888:localhost:8888 username@serveraddress
```

Finally, open a web browser (Chrome, Chromium, Firefox works in my experience), paste the address from your terminal (http://localhost:8888/?token=600be64ca79ef74ebe4fbd3698bc8a0d049e01d4e28b30ec in this example)  
Navigate to **resources/analysis.ipyn** and follow the instructions contained in the notebook.
### Test run
A test_data folder and test_analysis folder are included in the base_resources directory. Use them as your data_dir and analysis_dir, respectively. You'll need 2 additional directories for project_resources and species_resources to start the jupyter notebook. Make a copy of the analysis.ipynb and compare if your results are the same as the one in the provided copy.

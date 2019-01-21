MIPTools
=========
MIPTools is a suit of computational tools that are used for molecular inversion probe design, data processing and analysis.
## Installation
Clone the github repository:
```bash
git clone git@github.com:bailey-lab/MIPMaker.git
```
### Dependencies
Requires a working copy of Singularity: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps  
Singularity is best installed with **sudo**. While it is said to be possible to install with unprivilidged user with some features missing, MIPTools hasn't been tested on such an installation yet.
### Build MIPMaker from the definition file 
This can take about 10-30 minutes, depending on number of cpu cores available.  
User must have **sudo** privilege to _build_ the image. You do not need sudo to _use_ the image. So if you want to run the container on an environment without sudo, build the container on your own machine and copy the image file to the host machine. Note that Singularity program itself must have been installed with sudo.
```bash
cd MIPMaker
sudo singularity build miptools.sif MIPMaker.def
```
mipmaker.sif is a single **portable** file which has all the programs needed for MIP data analysis and a lot more.  
More information about its uses will be added over time.
### Directory Structure
Following directories are needed for MIP data analysis.
*  **base_resources:** Included in the repository, contains common resources across projects.
*  **species_resources:** Contains resources shared by projects using the same target species (Pf, human, etc.)
*  **project_resources:** Contains project specific files (probe sequences, sample information, etc.)
*  **data_dir:** Contains data to be analyzed (see below).
*  **analysis_dir:** Where analysis will be carried out and all output files will be saved.
### Usage
The computer on which the container is running on is called the **host**. So if you're running this on your laptop, your laptop is the host. If you're connected to a HPC cluster and running the container there, the HPC is the **host**.  

The directories  listed above must be on the host computer. They can be paths to any directory on the host where you have **read** permission for all directories and **write** permission for the analysis_dir.  

miptools.sif is a Singularity image file. It contains all parts of an (Ubuntu) operating system that are required to do everything MIP. Although it contains all programs needed, it does not include the data to be analyzed or the resources to be used. Every time we run Singularity we will **bind** needed directories to the container. The 3 resources directories are constant across operations. In addition to those, each app needs a data_dir and analysis_dir. **-B** option is used for each binding: 
```bash
singularity -B path_on_host:path_on_container
```
Path on the left side of the column specifies where on the host computer the directory is and the right side is the location in the container where the directory should be bound (mounted) to. Each binding is specified with a separate -B option.

### Apps
Each MIPtool is an **app** in the container. This is a typical Singularity command to run an app:  
```bash
singularity run --app appName singularityOptions miptools.sif appOptions
```

#### wrangler
Runs MIPWrangler.   

Run the following command, changing only the **host** part of the binding arguments (-B) to fit your directory structure.
```bash
singularity run --app wrangler \
    -B base_resources:/opt/resources \
    -B project_resources:/opt/project_resources \
    -B species_resources:/opt/species_resources \
    -B data_dir:/opt/data \
    -B analysis_dir:/opt/analysis \
    miptools.sif -p probe_sets_used -s sample_set -e experiment_id -l sample_list_file
```

You must have the sample_list file in the analysis_dir. For the test_data, go to the test_data dir and run:
```bash
singularity run --app wrangler \
    -B base_resources/:/opt/base_resources \
    -B DR1_project_resources/:/opt/project_resources \
    -B pf_resources/:/opt/species_resources \
    -B fastq/:/opt/data \
    -B wrangler_dir/:/opt/analysis \
    miptools.sif -e wrangler -l 170828_samples.tsv -p DR1,CSP -s SM -c 4
```
After the app finishes, it creates a file with all _clean haplotype_ sequences (and a lot more) in a subdirectory called analsis within the wrangler_dir. Create a link to the main output file to set up for the next step in analysis:  
```bash
ln -s wrangler_dir/analysis/serverResources/mip1/popClusInfo/allInfo.tab.txt.gz wrangler_dir/
```
#### mapper
Runs post-wrangler mapping and variant calling pipeline. For the test_data:    
```bash
singularity run --app jupyter \
    -B base_resources/:/opt/base_resources \
    -B DR1_project_resources/:/opt/project_resources \
    -B pf_resources/:/opt/species_resources \
    -B wrangler_dir/:/opt/data \
    -B mapper_dir/:/opt/analysis \
    miptools.sif
```

This will start a Jupyter notebook where the rest of the analysis takes place.

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
Navigate to **analysis/analysis.ipynb** and follow the instructions contained in the notebook.

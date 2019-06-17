MIPTools
=========
MIPTools is a suit of computational tools that are used for molecular inversion probe design, data processing and analysis.
## Installation
Clone the github repository:
```bash
git clone git@github.com:bailey-lab/MIPTools.git
```
### Dependencies
Requires a working copy of Singularity: https://www.sylabs.io/docs/  
Singularity is best installed with **sudo**. While it is said to be possible to install with unprivilidged user with some features missing, MIPTools hasn't been tested on such an installation.

Note: Snap package install is a rapid way to install the go language required by Singularity (e.g. on Ubuntu/Debian: `sudo snap install go --classic`)

### Build MIPTools from the definition file 
This can take about 10-30 minutes, depending on number of cpu cores available.  
User must have **sudo** privilege to _build_ the image. You do not need sudo to _use_ the image. So if you want to run the container on an environment without sudo, build the container on your own machine and copy the image file to the computer without sudo. Note that Singularity program itself must have been installed with sudo.
```bash
cd MIPTools
sudo singularity build miptools.sif MIPTools.def
```
miptools.sif is a single **portable** file which has all the programs needed for MIP design, data analysis and a lot more.  
More information about the extra programs and their uses will be added over time.

Although miptools.sif contains all programs needed, it does not include the data to be analyzed or the resources to be used. Every time we run Singularity we will **bind** needed directories to the container. There are 3 resources directories which are required for most operations. In addition to those, each app needs a data_dir and analysis_dir. **-B** option is used for each binding: 
```bash
singularity -B path_on_host:path_on_container
```
Path on the left side of the column specifies where on *your* computer the directory is and the right side is the location in the container where the directory should be bound (mounted) to. Each binding is specified with a separate -B option. See below for examples.

### Directory Structure
3 resource directories are required for most operations. These live outside the container and must be **bound** to the container at run time with `-B` option.
* **base_resources:** Included in the GitHub repository. It contains common resources across projects. It should be bound to the container with `-B [path to base resources dir outside of the container]:/opt/resources`. This makes the base_resources directory available to the container and it would be reached at `/opt/resources` path within the container. `/opt/resources` part of this argument must not be altered. For example, if my base resources are located in my computer at `/home/base`, I would bind it to the container with `-B /home/base:/opt/resources`.
* **species_resources:** Contains resources shared by projects using the same target species (Pf, human, etc.). Bind this to `/opt/species_resources` in the container.  

  * What should be in the species_resources?
* **project_resources:** Contains project specific files (probe sequences, sample information, etc.). Bind this to `/opt/project_resources`
* **data_dir:** Contains data to be analyzed. Typically, this nothing will be written to this directory. Bind this directory to `/opt/data`.
* **analysis_dir:** Where analysis will be carried out and all output files will be saved. Bind it to `/opt/analysis` This is the only directory that needs write permission as the output will be saved here.

Usually, one app's analysis directory will be the next app's data directory in the pipeline.

### Usage
Each MIPtool is an **app** in the container. This is a typical Singularity command to run an app:  
```bash
singularity run --app appName singularityOptions miptools.sif appOptions
```
#### download
This app should be used to download sequencing run data from basespace. An access token obtained from BaseSpace is required for downloading data from basespace. This accestoken should be placed in the acces_token.txt file within the base resource directory.

```bash
singularity run --app download \  
    -B [path to base resources]:/opt/resources \  
    -B [path to save sequence data to (bcl directory)]:/opt/analysis \  
    [path to miptools.sif] -r [run ID from BaseSpace]
```
The above command will save the sequencing data to a subfolder with the same name as the run ID within the specified bcl  directory. This directory containing the sequence data (bcl directory/run ID) is referred to as bcl directory below.

#### demux
Sample demultiplexing. This app generates per-sample fastq files from the sequence data downloaded with the **download** app.

```bash
singularity run --app demux \  
    -B [path to base resources]:/opt/resources \  
    -B [path to bcl directory (from download app)]:/opt/data \  
    -B [path to fastq directory where fastq files should be saved to]:/opt/analysis |  
    [path to miptools.sif] -s sample_list -p sequencing platform
```

*  *sample_list*: A file listing the samples used in the study, primers used etc. An example can be found in the test data set. It is important that sample names contain only alphanumeric vaules and dashes. No underscores, spaces or other special characters are allowed. It is also important not to change the field names of the sample list provided. This file must be present in the data directory (bcl directory).
*  *sequencing platform*: The sequencing platform used. This must be one of miseq or nextseq.

#### demux_qc
This app goes through the demultiplexing stats and prints the total number of sequencing reads and how many were of undetermined indices, meaning they contained indices that were not present in the sample file. In addition, it prints how many of the undetermined index reads belongs to possible primer pairs, i.e. they are likely due to errors in provided sample list and not just faulty reads from the sequencer.

```bash
singularity run --app demux_qc \  
    -B [path to base resources]:/opt/resources \  
    -B [path to fastq directory where fastq files were saved]:/opt/analysis \  
   [path to miptools.sif] -p sequencing platform
```

#### wrangler

This app runs MIPWrangler on the demultiplexed fastq files.

```bash
singularity run --app wrangler \  
    -B [path to base resources]:/opt/resources \  
    -B [path to project resources]:/opt/project_resources \  
    -B [path to fastq directory]:/opt/data \  
    -B [path to MIPWrangler output directory(wrangler dir)]:/opt/analysis \  
    [path to miptools.sif] -p probe_sets_used -s sample_sets -e experiment_id -l sample_list -c cpu_number
```

*  *probe_sets_used*: A comma separated list of probe sets used. MIP probes are grouped into sets for ease of use. It is possible to use many probe sets in an experiment. All probe related information is located in project_resources/mip_ids directories. More information on this will be added. Probe sets must be specified in the sample list as well. However, the sets provided here may be a subset of those present in the sample list file. Only the probes in specified here will be analysed.
*  *sample_set*: Another field in the sample list is "sample_set". This is used for grouping the samples for analysis. Sample sets that are meant to be analysed together should be given as a comma separated list here. If multiple sample sets are to be analysed together, they must have the same *probe_sets*.
*  *experiment_id*: A unique id given to the experiment. This is typically the date of the sequencing run (YYMMDD) but can be anything alphanumeric.
*  *sample_list*: same as the sample list used in demultiplexing. This file must be present in the *analysis_dir*
After the app finishes, it generates the MIPWrangler output and saves it to the wrangler directory. A number of files for quality control, such as proportion of reads failing each step of the pipeline, are also generated and saved in the wrangler directory specified.
*  *cpu_number*: number of processors to use.

#### mapper
Runs post-wrangler mapping and variant calling pipeline. 
```bash
singularity run --app jupyter \
    -B [path to base resources]:/opt/resources \  
    -B [path to project resources]:/opt/project_resources \  
    -B [path to species resources]:/opt/species_resources \  
    -B [path to MIPWrangler output directory(wrangler dir)]:/opt/data \  
    -B [path to analysis directory]:/opt/analysis \  
    [path to miptools.sif]
```
If multiple sequence runs will be combined and analysed, they need to be present in the data dir. 
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
I cloned the MIPTools repository to ~/git location in my computer, so the base_resources location is ~/git/MIPTools/base_resources and the test data is located in ~/git/MIPTools/test_data/fastq.

I built miptools.sif to my home directory so miptools.sif path is ~/miptools.sif
The project resources path on my computer is ~/resources/project_resources/DR1
The species resources path on my computer is ~/resources/species_resources/pf/Pf_3D7_v3-PDB9.3/
I created an analysis directory for MIPWrangler output: ~/fast_data/analysis/test_analysis and copied the sample list located in the base_resources directory to the analylsis directory. So the sample list path is ~/fast_data/analysis/test_analysis/sample_list.tsv

I created an analysis directory~/fast_data/analysis/test_mapper
For the above setup, the MIPWrangler command is this:
```
bash
singularity run --app wrangler -B ~/git/MIPTools/base_resources/:/opt/resources -B ~/resources/project_resources/DR1/:/opt/project_resources -B ~/resources/species_resources/pf/Pf_3D7_v3-PDB9.3/:/opt/species_resources -B ~/git/MIPTools/test_data/fastq/:/opt/data -B ~/fast_data/analysis/test_analysis/:/opt/analysis ~/miptools.sif -p DR1,VAR4 -s JJJ -l sample_list.tsv -e test_id -c 20
```

That creates the file ~/fast_data/analysis/test_analysis/test_id_JJJ_DR1_VAR4_20190412.txt.gz. Note that the file name contains the data and you will have a different file name.



Finally, open a web browser, paste the address from your terminal (http://localhost:8888/?token=600be64ca79ef74ebe4fbd3698bc8a0d049e01d4e28b30ec in this example)  
Navigate to **analysis/analysis.ipynb** and follow the instructions contained in the notebook.

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
* **species_resources:** Contains resources shared by projects using the same target species (Pf, human, etc.). Bind this to `/opt/species_resources` in the container. For example, if I am working with *Plasmodium falciparum* sequences and I have the necessary files in my computer at /home/pf3d/, then the binding parameter is -B /home/pf3d:/opt/species_resources.
   species_resources contents:
   * file_locations.tsv: a tab separated text file showing where each required file will be located in the container. Each line corresponds to one file. First field states the species for the file, second field states what kind of file it is and the last field is the absolute path to the file.  
   For example, *'pf      fasta_genome    /opt/species_resources/genomes/genome.fa'* shows the fasta genome file for the species 'pf' will be found at '/opt/species_resources/genomes/genome.fa' within the container. This also means that there is a file at /home/pf3d/genomes/genome.fa in my computer, because I bound /home/pf3d to /opt/species_resources in the container.
   * fasta file: Genome reference sequence in fasta format.
   * bowtie2_genome: Reference genome indexed using bowtie2. If this is not available, it can be generated using MIPTools.
   * bwa_genome: Reference genome indexed using bwa. If this is not available, it can be generated using MIPTools.
   * snps: A vcf file containing genomic variation. Individual genotypes are not necessary. The only requirement is that either the INFO fields AC and AN are present for number of allele counts and total allele counts, respectively. Tese can be number of samples containing the allele and total number of samples, instead of allele  numbers. Because the information is used to get an idea about the population frequency of each allele, frequency of samples carrying an allele provides a good approximation. If these numbers are provided with different field names than AC and AN, the field names used must be specified in design settings file.
   * refgene: RefGen/RefSeq style gene/gene prediction table. These are available ad http://genome.ucsc.edu under Tools/Table Browser for most species. The fields in the file are "bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames". 
   * refgene_tabix: RefGen file, sorted and indexed using tabix.
* **project_resources:** Contains project specific files (probe sequences, sample information, etc.). Bind this to `/opt/project_resources`
* **data_dir:** Contains data to be analyzed. Typically, this nothing will be written to this directory. Bind this directory to `/opt/data`.
* **analysis_dir:** Where analysis will be carried out and all output files will be saved. Bind it to `/opt/analysis` This is the only directory that needs write permission as the output will be saved here.

Usually, one app's analysis directory will be the next app's data directory in the pipeline.

## Usage for MIP design
### Region Prep
The first step in the probe design process is to create target regions. Targets can be provided in four ways. All files explained below must be placed in the *project_resources* directory if they will be used. 
**Important:** all target names must not include any specieal characters except "-" (dash). Letters, numbers and dashes are the only allowed characters, except for gene names, where the exact gene name  must be used. If any specieal character is used, these will be replaced with a dash. 
1) SNPs/Targets can be specified in a tab separated text file with headers: *Chrom, Start, End, Name, InsertionLength*. All field names are case sensitive. First 3 fields describe the genomic location of the target. *Name* is a unique name given to the target. *InsertionLength* specifies the maximum length change for insertions relative to the reference genome. This information is used to shorten the MIP captured region to accomodate the largest insertion in the target region. For example, let's say we have a capture size limit of 10 bp (this limit depends on the sequencing platform used and will be discussed later) and a target specified as *chr1 10 13 test-target 0*. If a MIP captures the coordinates between 5 and 15 of chr1, that MIP is considered capturing this target. If we know that this target is in a tandem repeat region and known to have insertions up to 6 bases, we'd specify it as  *chr1 10 13 test-target 6*. Now the same MIP capturing chr1:5-15 will be discarded because the captured region may be 16 bp for samples carrying the 6 bp insertion and the size limit is reduced to 4 (10 - 6). Instead, another MIP capturing chr1:9-14 will be needed.
2) Region coordinates can be specified in a tab separated text file with headers *Name, Chrom, Start, End, CaptureType*. These coordinates provide a target region for MIP design as opposed to specific targets for MIPs to capture. *CaptureType* must be one of "exons" or "whole". For example, *chr1 100 300 test-target whole* would use the chr1:100-300 sequence (flanked by a length of bases on each side specified later in the process), as template to design MIPs on. If CaptureType was specified as "exons", then only the exons overlapping chr1:100-300 would be used as template.
3) Gene names can be specified in a tab separated text file with headers: *Gene, GeneID, CaptureType*. Gene field can be any name, but it makes sense to use actual gene names. GeneID must be the value that is present in the *name2* column of the refgene file. *CaptureType* must be one of "targets", "exons" or "whole". Specifying "targets" as the CaptureType here is useful to group the targets specifed in the other files under a more meaninful name. "exons" would use exons of the gene as template while "whole" would use introns as well.
4) Fasta sequences can be provided in a single (possibly multi-sequence) fasta file or multiple fasta files. The sequences provided here are aligned to the reference genome and the coordinates from this alignment is used as template, not the actual sequence provided in the fasta file.

"Shell" into the container using correct path bindings:
```bash
singularity shell \  
    -B [path to base resources]:/opt/resources \  
    -B [path to project resources]:/opt/project_resources \  
    -B [path to species resources]:/opt/species_resources \  
    -B [path to design directory]:/opt/analysis \  
    [path to miptools.sif]
```
After the above command you should get a command prompt different from your usual one, similar to:
```bash
Singularity miptools:~>
```
This means you are working in the container and all paths must refer to those in the container. For example, the species_resources directory is always at /opt/species_resources, regardless of where it is in your computer.  
If indexed genomes are not available (bowtie2, bwa), create them:
```bash
elucidator bioIndexGenome --genomeFnp [path to fasta file in the container] --verbose
```
If you created the indexed genomes now, update the bowtie_genome and bwa_genome entries in your `/opt/species_resources/file_locations.tsv` file.

```bash
python /opt/extras/region_prep.py  OPTIONS
```
```
-n A short name for the current design.
-s Target species name.
[-o ]: Host species name. Probes predicted to bind this genome will be discarded. default=None
[-p ]: Number of available processors. default=7
[--parallel-processes ]: Number of designs carried out in parallel. default=1
[--flank ]: Number of bases to flank target sites on each side. default=150
[--single-mip-threshold ]: Target regions smaller than this will have a single MIP designed. default=0
[--min-target-size ]: Length threshold to eliminate un-aligned targets. default=150
[--max-indel-size ]: Indel size limit between paralogs, allowing this size gapped alignment within paralogs. default=50

[--coordinates-file ]: Path to file containing target coordinates.
[--genes-file ]: Path to file containing target genes.
[--snps-file ]: Path to file containing SNP coordinates.
[--fasta-files]: Path(s) to fasta file(s) containing target sequences
[--fasta-capture-type ]: Capture type for targets specified in fasta files. default="whole"
[--genome-coverage ]: Minimum alignment length to reference genome. default=1000
[--genome-identity ]: Minimum sequence identity (0-100) for genomic alignments. default=100
[--local-coverage ]: Minimum alignment length to reference paralog. default=100
[--local-identity ]: Minimum sequence identity (0-100) for paralog alignments. default=100
[--design-dir ]: Path to location to output design files. default="/opt/analysis"
[--resource-dir ]: Path to location to output prep files. default="/opt/project_resources"
[--targets-rinfo-template ]: Path to rinfo template for 'targets' capture type.default= "/opt/resources/templates/rinfo_templates/template_rinfo.txt"
[--exons-rinfo-template ]: Path to rinfo template for 'exons' capture type.default= "/opt/resources/templates/rinfo_templates/template_rinfo.txt"
[--whole-rinfo-template ]: Path to rinfo template for 'whole' capture type.default= "/opt/resources/templates/rinfo_templates/template_rinfo.txt"
[--output-file OUTPUT_FILE]: Base name to save region prep results. default="/opt/project_resources/design_info"
```





### Usage for data analysis
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

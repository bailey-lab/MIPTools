changes
- update the miptools_analysis_no_jupyter and setup_run and create_profile to include a miptools analysis folder
- make a copy of the mip_functions file
- change the mip_functions file in all of the scripts mip_functions as mip
- update freebayes_call mip function to stop after generating the contig_dict_list and only return that as an output
- update the call_variants_freebayes.py script so the only return is the contig_dict_list and the output is a contig_dict_list yaml
- add a run_freebayes rule to the variant_calling snakemake file
- add a concatenate and fix vcf headers rule



# miptools_analysis_no_jupyter

This is a snakemake pipeline for running a modified analysis_template_with_qual
notebook (from miptools analysis). It consists of three parts:

 - setup_run: creates a singularity profile for running remaining steps
 - check_run_stats: checks which samples and mips worked and at what levels
 - variant_calling: calls variants and generates output tables

Currently each step requires that the previous steps have run - you can't check
run stats without setting up the singularity profile, and you can't run the
variant calling step without setting up the run and checking run stats.

## Installation

 - Install conda: https://github.com/conda-forge/miniforge#unix-like-platforms-mac-os--linux.
You'll need to follow the instructions to 'initialize' the conda environment at the end of the
installer, then sign out and back in again.
 - Create a conda environment and install snakemake there:
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

### Setup your environment:
 - Change directory to a folder where you want to run the analysis
 - Download the files of this repository into that folder

## Usage:
 - Edit the miptools_analysis_no_jupyter.yaml file to point to your files.
Use a text editor that outputs unix line endings (e.g. vscode, notepad++, gedit, micro, emacs, vim, vi, etc.)
 - If snakemake is not your active conda environment, activate snakemake with:
```bash
conda activate snakemake
```
 - To setup your singularity environment, use:
```bash
snakemake -s setup_run.smk --cores 4
```
 - To check run stats, use:
```bash
snakemake -s check_run_stats.smk --profile singularity_profile
```
 - To run variant calling, use:
```bash
snakemake -s variant_calling.smk --profile singularity_profile
```
 - to run all three steps at once, use:
 ```bash
 bash run_all_steps.sh
 ```

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

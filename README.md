# MIPTools

[![Build Singularity](https://github.com/bailey-lab/MIPTools/actions/workflows/build-container.yaml/badge.svg)](https://github.com/bailey-lab/MIPTools/actions/workflows/build-container.yaml)
[![Deploy Docs](https://github.com/bailey-lab/MIPTools/actions/workflows/deploy-docs.yaml/badge.svg)](https://github.com/bailey-lab/MIPTools/actions/workflows/deploy-docs.yaml)
![GitHub release (latest
SemVer)](https://img.shields.io/github/v/release/bailey-lab/MIPTools)
![License](https://img.shields.io/github/license/bailey-lab/MIPTools)

MIPTools is a suite of computational tools that are used for molecular
inversion probe design, data processing, and analysis.

## Installation

### Dependencies

A working copy of Singularity is required: https://www.sylabs.io/docs/.
Singularity is best installed with **sudo**. While it is said to be possible to
install with unprivileged user with some features missing, MIPTools hasn't been
tested on such an installation.

Singularity is available for most Linux systems. It is also possible to install
and use on Mac OS using virtual machines with a little bit of extra work.

Note that the `snap` package is a rapid way to install the go language required
by Singularity (e.g. on Ubuntu/Debian: `sudo snap install go --classic`).

### Obtaining MIPTools

#### Download prebuilt container

The MIPTools container, built and ready to use, can be
downloaded from the [Sylabs Cloud](https://cloud.sylabs.io/). You can download
either the development version or the most recent stable release:

```bash
# Download the development version
# The development version is updated every two weeks
singularity pull library://apascha1/miptools/miptools:dev

# Download the latest stable release
singularity pull library://apascha1/miptools/miptools:v1.0.0
```

Note that these prebuilt versions do not include the `bcl2fastq` software due
to its license. If you plan to use MIPTools to demultiplex bcl files, you must
build the container yourself.

#### Install from source

MIPTools can also be built from source code using the definition file provided
in this [GitHub repository](https://github.com/bailey-lab/MIPTools).

The process can take about 10-30 minutes to build, depending on the number of
CPU cores available. By default, the build process will use 6 CPU cores. This
should pose no problems with most modern computers, but if the computer used
for building the container has less then 6 cpu cores available, change the
`"CPU_COUNT=6"` value at the top of the `MIPTools.def` file to a suitable
number before running the following code. On the other hand, if you have access
to more CPU power, by all means, use them by setting the same parameter to a
higher value.

You must have **sudo** privelege to _build_ the image. You do not need sudo to
_use_ the image. So if you want to run the container on an environment without
sudo, either download a prebuilt image (see above) or build the container on
your own machine where you _do_ have sudo privilege and copy the image file to
the computer without sudo. Note that the Singularity program itself must have
been installed with sudo.

If you plan to use MIPTools to demultiplex bcl files, you should download
`bcl2fastq` separately. Currently, you can download it from
[here](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html),
but this may change in the future. You must download the file: **`bcl2fastq2 Conversion Software v2.20 Installer (Linux rpm)`** and place it in the
`MIPTools/programs` directory.

You can install the most recent release using the following:

```bash
# Install stable version v1.0.0
git clone --b v1.0.0 https://github.com/bailey-lab/MIPTools.git
```

You can alternatively install the development version:

```bash
# Install dev version
git clone https://github.com/bailey-lab/MIPTools.git
```

Next, simply build the container and you should be all set to get started using
MIPTools!

```bash
cd MIPTools
sudo singularity build miptools.sif MIPTools.def
```

`miptools.sif` is a single **portable** file which has all the programs needed
for MIP design, data analysis, and a lot more. More information about the extra
programs and their uses will be added over time.

### Using MIPTools

Although `miptools.sif` contains all programs needed, it does not include the
data to be analyzed or other resources to be used. Every time we run
Singularity we will **bind** needed directories to the container. There are
three resources directories which are required for most operations. In addition
to those, some apps need a `data_dir` and `analysis_dir`. The **`-B`** option is
used for each binding:

```bash
singularity some-command -B path_on_host:path_on_container path_to_container
```

The path on the left side of the colon specifies where on _your_ computer the
directory is and the right side is the location in the container where the
directory should be bound (mounted) to. You should only change the left side of
the column according to the location of the resource you are providing, and
should _never_ change the path on the right side. Each binding is specified
with a separate `-B` option. See below for examples.

### Directory Structure

Three resource directories are required for most operations. These live outside
the container and must be **bound** to the container at run time with the `-B`
option. In addition, a data directory and an analysis directory will be used for
most
operations.

<details><summary>Expand for details on the directory structure and container
binding. </summary> <p>

- **base_resources:** Provided in the GitHub repository. It contains common
  resources across projects. It should be bound to the container with `-B [path to base resources dir outside of the container]:/opt/resources`. This
  makes the base_resources directory available to the container and it would
  be reached at `/opt/resources` path within the container. The
  `/opt/resources` part of this argument must not be altered. For example, if
  my base resources are located in my computer at `/home/base`, I would bind
  it to the container with `-B /home/base:/opt/resources`.

- **species_resources:** Contains resources shared by projects using the same
  target species (Pf, human, etc.). Bind this to `/opt/species_resources` in
  the container. For example, if I am working with _Plasmodium falciparum_
  sequences and I have the necessary files in my computer at `/home/pf3d/`,
  then the binding parameter is `-B /home/pf3d:/opt/species_resources`.

  _Contents of species_resources directory:_

  - _file_locations.tsv:_ This file is required for all operations. It is a
    tab separated text file showing where each required file will be
    located in the container. Each line corresponds to one file. First
    field states the species for the file, second field states what kind of
    file it is and the last field is the absolute path to the file.

    For example, the line _"pf &nbsp; &nbsp; &nbsp; &nbsp; fasta_genome
    &nbsp; &nbsp; &nbsp; &nbsp; /opt/species_resources/genomes/genome.fa"_
    would mean that the fasta genome file for the species 'pf' will be
    found at '/opt/species_resources/genomes/genome.fa' within the
    container. This also means that there is a file at
    /home/pf3d/genomes/genome.fa in my computer, assuming I bound
    /home/pf3d to /opt/species_resources in the container.

  - _fasta file:_ This file is required for all operations. Genome
    reference sequence in fasta format.

  - _bowtie2_genome:_ This file is required for probe design operations
    only. It is the reference genome indexed using bowtie2. If this is not
    available, it can be generated using MIPTools.

  - _bwa_genome:_ This file is required for data analysis operations only.
    It is the reference genome indexed using bwa. If this is not available,
    it can be generated using MIPTools.

  - _snps:_ This is an optional file. However, it is extremely useful in
    probe designs to avoid probe arms landing on variant regions, etc. So
    it should always be used except in rare cases where such a file is not
    available for the target species. The format of the file is vcf.
    Individual genotypes are not necessary (a.k.a. sites only vcf). The
    only requirement is that the INFO field for each variant has a field
    showing the population allele frequency of alternate alleles. By
    default, AF field is used. The AF field lists the allele frequencies of
    each alternate allele, and does not list the frequency of the reference
    allele. Vcf files may have other INFO fields that include allele
    frequency information. If such a field is to be used, there are two
    settings in the design settings file (.rinfo file) that must be
    modified. _allele_frequency_name_ field must be set to the INFO field
    name to be used; _af_start_index_ may have to be set to a 1 (instead of
    default 0) depending on whether the reference allele frequency is
    provided in the new field. For example, if we want to use the 1000
    genomes vcf file, the allele frequencies are provided in the CAF field
    and they include the reference allele. We would have to change the
    _allele_frequency_name_ field to _CAF_ from the default _AF_; and set
    _af_start_index_ to 1 because the first alternate allele's frequency is
    provided in the second place (following the reference allele).

  - _refgene:_ RefGen style gene/gene prediction table in GenePred format.
    These are available at http://genome.ucsc.edu under Tools/Table Browser
    for most species. The fields in the file are "bin, name, chrom, strand,
    txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds,
    score, name2, cdsStartStat, cdsEndStat, exonFrames". This file is
    required for probe design operations if genic information is to be
    used. For example, if probes need to be designed for exons of a gene,
    or a gene name is given as design target. If a gene name will be
    provided, it must match the **name2** column of the RefGen file. If you
    are creating this file manually, the only fields necessary are: chrom,
    strand, exonStarts, exonEnds and name2. All other fields can be set to
    an arbitrary value (none, for example) but not left empty. The order of
    columns must not be changed.

    Note: If you have gff3/gtf formatted files, they can be converted to
    GenePred format using Jim Kent's programs
    [gff3ToGenePred](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred)
    and
    [gtfToGenePred](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred).

  - _refgene_tabix:_ RefGen file, sorted and indexed using tabix. File
    requirement is the same as the refgene file. tabix is available within
    the MIPTools container, so you don't have to install it yourself.

- **project_resources:** Contains project specific files (probe sequences,
  sample information, etc.). Bind this to `/opt/project_resources`

- **data_dir:** Contains data to be analyzed. Typically, nothing will be
  written to this directory. Bind this directory to `/opt/data`.

- **analysis_dir:** Where analysis will be carried out and all output files
  will be saved. Bind it to `/opt/analysis` This is the only directory that
  needs write permission as the output will be saved here.

`data_dir` and `analysis_dir` will have different content for different app
operations. Also, one app's analysis directory may be the next app's data
directory in the pipeline.

</p> </details>

## Resource requirements

Resources required vary widely depending on the project. Both designs and data
analysis can be parallelized, so the more CPUs you have the better. Plenty of
storage is also recommended. For designs on large target regions (>5kb), files
can take up 10 GB or more per region. Consider allocating > 5 GB RAM for a
large design region (multiply the RAM requirement by CPU number if
parallelizing). For a typical MIP data analysis involving ~1000 MIPs and ~1000
samples, consider using at least 20 CPUs and 20 GB RAM to get the analysis done
within 10-12 h. You should expect ~200 GB disk space used for such an analysis
as well, although most files can be removed after processing steps to reduce
long term disk usage.

## Further documentation

Further documentation for MIPTools is available
[here](https://drive.google.com/drive/folders/1Tmu7hdRYrdw-jqAN35lZpIjG2lBebuCK?usp=sharing)
for various use cases (MIP design, data analysis, etc.).

## Troubleshooting

Please send any questions/comments to: <miptools@googlegroups.com>.

Join our Google Group: <https://groups.google.com/g/miptools>

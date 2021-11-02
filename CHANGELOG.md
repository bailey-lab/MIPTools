# MIPTools Changelog

## MIPTools (development version)

### App Changes

- New `bs_download` app replaces the `download` app. The new app improves the
  method to download data from the Illumina BaseSpace Sequence Hub by using
  command line tools (@arisp99, #13).
- The `demux` app no longer requires both a run directory and a output
  directory. These directories have been combined so that fastq files are
  output to the run directory.

### Documentation Overhaul

- Add doc-strings to python functions.
- Improve clarity of README and add additional instructions on downloading or
  building the container.

### Bug Fixes

- Fix build failure due to dependency changes (#7).

### Maintenance

- Automatically build and deploy container using Github Actions (@arisp99, #11).
- Remove duplicated files.
- Improve bash errors.
- Make strings human readable (@arisp99, #5).

## MIPTools 1.0.0

- First major release.

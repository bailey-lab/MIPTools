# MIPTools Changelog

## MIPTools (development version)

### New Features

- New `download` app supersedes the previous `download` app, which has been
  renamed to `download_superseded`. The new app improves the method for
  downloading data from the Illumina BaseSpace Sequence Hub by using the
  official command line tool (@arisp99, #25, #13).

### Bug Fixes

- Fix missing file error when MIP arms file is created from the MIP info
  dictionary (@aydemiro, #23).
- Improve sample sheet preparation. Avoid errors when sample file columns are
  empty. Throw an error if there are invalid samples or input fields
  (@aydemiro, #22).
- Fix build failure due to dependency changes in the McCOILR R package (#7).

### Maintenance

- Update LICENSE year.
- Store copies of container on Sylabs Cloud (#12).
- Remove duplicated files.
- Improve bash errors.
- Make strings human readable (@arisp99, #5).

### Documentation Overhaul

- Generate online documentation using
  [Sphinx](https://www.sphinx-doc.org/en/master/index.html) and [Github
  Pages](https://pages.github.com/).
- Add better documentation for the `jupyter` app.
- Add better documentation for the `wrangler` app.
- Add better documentation for the `download` app.
- Add better documentation for the `demux_qc` app.
- Add doc-strings to python functions.
- Improve clarity of README and add additional instructions on downloading or
  building the container.

## MIPTools 0.4.0

- Latest stable build.

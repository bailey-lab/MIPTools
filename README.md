# MIPTools Documentation

[![Deploy Docs](https://github.com/bailey-lab/MIPTools/actions/workflows/deploy-docs.yaml/badge.svg)](https://github.com/bailey-lab/MIPTools/actions/workflows/deploy-docs.yaml)
![GitHub release (latest
SemVer)](https://img.shields.io/github/v/release/bailey-lab/MIPTools)
![License](https://img.shields.io/github/license/bailey-lab/MIPTools)

This folder contains the source code for the MIPTools documentation. Our
documentation can either be written in
[rst](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
or [MyST](https://myst-parser.readthedocs.io/en/latest/index.html) (a extension
of [markdown](https://www.markdownguide.org/)), both simple markup languages.
We use the documentation tool,
[sphinx](https://www.sphinx-doc.org/en/master/index.html), to build our
documentation, and host our final product on [Github
Pages](https://bailey-lab.github.io/MIPTools/).

## Contributing

Please feel free to contribute by submitting a pull request! To build and
preview your changes, you may use the `Makefile` provided.

```shell
# Build html docs
make html

# Preview docs
open _build/html/index.html
```

You may also use [sphinx autobuild](https://pypi.org/project/sphinx-autobuild/)
to view changes live.

```shell
sphinx-autobuild . _build/html
```

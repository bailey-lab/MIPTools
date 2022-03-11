#!/bin/bash

# Build container
sudo singularity build --force MIPTools/miptools_dev.sif MIPTools/MIPTools.def

# Upload to website
mv -f MIPTools/miptools_dev.sif /work/bailey_share/software/MIPTools/containers/

#!/bin/bash

# Build container
sudo singularity build --force miptools_dev.sif MIPTools.def

# Upload to website
mv -f miptools_dev.sif /work/bailey_share/software/MIPTools/download/

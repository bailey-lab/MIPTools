# be sure to update the version in MIPTools.def before running this script

eval $(grep -h "MIPTOOLS_VERSION=" MIPTools.def)

sudo singularity build miptools_$MIPTOOLS_VERSION.sif MIPTools.def
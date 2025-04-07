# this script reads the version of the miptools file in MIPTools.def and 
# creates a sif file with the same name.  It is a handy way to ensure the naming convention
# remains constant for newly created sif files

eval $(grep -h "MIPTOOLS_VERSION=" MIPTools.def)

sudo singularity build miptools_$MIPTOOLS_VERSION.sif MIPTools.def
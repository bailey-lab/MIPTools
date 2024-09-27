if [ -z $1 ]; then
    sif_output='miptools_dev.sif'
else
    sif_output=$1
fi

build_conda() {
    if [ ! -f conda.sif ]; then
        sudo singularity build conda.sif conda.def
    fi
}

build_wrangler() {
    if [ ! -f wrangler.sif ]; then
        sudo singularity build wrangler.sif wrangler.def
    fi
}

build_conda &
build_wrangler &

wait

sudo singularity build $sif_output MIPTools.def

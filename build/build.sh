#!/usr/bin/env bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$script_dir" || exit 1

# set the version here at build time, this will be applied to the run_miptools.sh script and the config.toml file
# inside the sif container, and will be used to tag the sif file itself
VERSION=dev

mkdir -p sif_files

if [ ! -e sif_files/miptools_base_$VERSION.sif ]; then
	sudo singularity build \
		--tmpdir=/var/tmp \
		--arch amd64 \
		sif_files/miptools_base_$VERSION.sif MIPTools_base.def
fi

sudo singularity build \
	--arch amd64 \
	--build-arg VERSION=$VERSION \
	--tmpdir=/var/tmp \
	--mksquashfs-args "-no-compression" \
	sif_files/miptools_$VERSION.sif MIPTools.def


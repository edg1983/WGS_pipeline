#!/usr/bin/bash
nextflow="/well/gel/HICF2/software/nextflow/nextflow_v20.10"
singularity_dir="/well/gel/HICF2/software/singularity"

mkdir -p $singularity_dir/cache_dir
mkdir -p $singularity_dir/tmp_dir

export SINGULARITY_CACHEDIR=${singularity_dir}/cache_dir
export SINGULARITY_TMPDIR=${singularity_dir}/tmp_dir

$nextflow WGS_analysis.nf "$@" --singularity_basedir ${singularity_dir}
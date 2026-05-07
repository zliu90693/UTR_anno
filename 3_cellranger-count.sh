#!/bin/bash

# Depends on conda environment: fastq2matrix (see .env/fastq2matrix.yml)

set -euo pipefail

########################################
# Parameters
########################################

project_name=""
localcores=20
transcriptome_name=""
create_bam=false

while getopts "p:l:t:c:" opt; do
    case $opt in
        p) project_name=$OPTARG ;;
        l) localcores=$OPTARG ;;
        t) transcriptome_name=$OPTARG ;;
        c) create_bam=$OPTARG ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires a parameter value." >&2
            exit 1
            ;;
    esac
done

# parameter valid?

if [ -z "$project_name" ]; then
    echo "Error: The project name must be specified using -p."
    exit 1
fi

if ! [[ "$localcores" =~ ^[0-9]+$ ]]; then
    echo "Error: -l must be a positive integer."
    exit 1
fi

if [ -z "$transcriptome_name" ]; then
    echo "Error: The transcriptome name must be specified using -t."
    exit 1
fi

if [[ "$create_bam" != "true" && "$create_bam" != "false" ]]; then
    echo "Error: -c must be true or false"
    exit 1
fi

########################################
# Functions
########################################

log() {
    echo
    echo "==================== $1 ===================="
    echo
}

########################################
# Pipeline Begin
########################################

mkdir -p "./${project_name}/cellranger-count-out"
cd "./${project_name}/cellranger-count-out"

get_sample_prefixes() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        echo "Error: '$dir' is not a directory." >&2
        return 1
    fi

    # 遍历目录下所有 .fastq.gz 文件，提取前缀，排序去重，用逗号拼接
    for f in "$dir"/*.fastq.gz; do
        file="${f##*/}"          # 去掉路径，只保留文件名
        echo "${file%%_*}"       # 取第一个 '_' 之前的部分
    done | sort -u | paste -sd,
}

for dir in ../fastq/*; do
    sample_prefixes=$(get_sample_prefixes $dir)
    GSM=$(basename "$dir")
    log "start processing ${GSM}"
    cellranger count --localcores "$localcores" \
        --id="$GSM" \
        --fastqs="../fastq/${GSM}" \
        --sample="$sample_prefixes" \
        --transcriptome="../ref/${transcriptome_name}" \
        --create-bam "$create_bam"
    log "${GSM} processing completed" # single GSM about 4.5 hour
done

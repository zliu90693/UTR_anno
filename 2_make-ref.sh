#!/bin/bash

# Depends on conda environment: fastq2matrix (see .env/fastq2matrix.xml)

set -euo pipefail

########################################
# Parameters
########################################

project_name=""
gtf_name=""
filtered_gtf_name=""
ref_fasta=""
mkref_output_name=""
cellranger_nthreads=20

while getopts "p:g:f:r:m:t:" opt; do
    case $opt in
        p) project_name=$OPTARG ;;
        g) gtf_name=$OPTARG ;;
        f) filtered_gtf_name=$OPTARG ;;
        r) ref_fasta=$OPTARG ;;
        m) mkref_output_name=$OPTARG ;;
        t) cellranger_nthreads=$OPTARG ;;
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

if [ -z "$gtf_name" ]; then
    echo "Error: The gtf file name must be specified using -g."
    exit 1
fi

if [ -z "$filtered_gtf_name" ]; then
    echo "Error: The filtered gtf file name must be specified using -f."
    exit 1
fi

if [ -z "$ref_fasta" ]; then
    echo "Error: The reference fasta file name must be specified using -r."
    exit 1
fi

if [ -z "$mkref_output_name" ]; then
    echo "Error: The mkref output name name must be specified using -m."
    exit 1
fi

if ! [[ "$cellranger_nthreads" =~ ^[0-9]+$ ]]; then
    echo "Error: -t must be a positive integer."
    exit 1
fi

# path exists?

if [ ! -d "${project_name}/ref" ]; then
    echo "Error: Lack ref directory"
    exit 1
fi

if [ ! -f "${project_name}/ref/${gtf_name}" ]; then
    echo "Error: Lack gtf file"
    exit 1
fi

if [ ! -f "${project_name}/ref/${ref_fasta}" ]; then
    echo "Error: Lack reference fasta file"
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

check_pattern() {
    local label="$1"
    local pattern="$2"
    local file="$3"
    local count
    count=$(grep -c "$pattern" "$file" 2>/dev/null || true)
    echo "==================== [CHECK] ===================="
    echo "Stage: $label"
    echo "target: $pattern"
    echo "file: $(basename $file)"
    echo "${count} fits in total"
}

check_biotypes() {
    local label="$1"
    local file="$2"
    echo "==================== [CHECK] ===================="
    echo "Stage: $label"
    echo "All gene_biotype: "
    grep -o 'gene_biotype "[^"]*"' "$file" | sort | uniq -c | sort -rn
}

########################################
# Pipeline Begin
########################################

cd "./${project_name}/ref"
mkdir -p .tmp

# refer: /home/liuzhiyu/Projects/caste/Caste_code_management/1_LZU_fastq2matrix/LZU_code/index.qmd

# > tmp.gtf
log "STEP 1: remove transcript_id \"\""
check_pattern "Before step 1" 'transcript_id "";' "$gtf_name"
sed 's/transcript_id ""; //' "$gtf_name" > ".tmp/tmp.gtf"
check_pattern "After step 1" 'transcript_id "";' ".tmp/tmp.gtf"
log "STEP 1 ended"

# > tmp1.gtf
log "STEP 2: remove composite fields"
check_pattern "Before step 2" 'note ".*; product "' ".tmp/tmp.gtf"
check_pattern "Before step 2" 'model_evidence ".*; transcript_biotype "' ".tmp/tmp.gtf"
check_pattern "Before step 2" ' $' ".tmp/tmp.gtf"
sed 's/ $//' ".tmp/tmp.gtf" | \
    perl -pe 's|note ".*?; product "|product "|' | \
    perl -pe 's|model_evidence ".*?; transcript_biotype "|transcript_biotype "|' \
    > ".tmp/tmp1.gtf"
check_pattern "After step 2" 'note ".*; product "' ".tmp/tmp1.gtf"
check_pattern "After step 2" 'model_evidence ".*; transcript_biotype "' ".tmp/tmp1.gtf"
check_pattern "After step 2" ' $' ".tmp/tmp1.gtf"
log "STEP 2 ended"

# > tmp2.gtf
log "STEP 3: Remove the quote in the middle of a value"
check_pattern "Before step 3" '" "' ".tmp/tmp1.gtf"
sed 's/" "/,/g' ".tmp/tmp1.gtf" > ".tmp/tmp2.gtf"
check_pattern "After step 3" '" "' ".tmp/tmp2.gtf"
log "STEP 3 ended"

# > filtered_gtf_name
log "STEP 4: Filt biotype using mkgtf"
check_biotypes "Before step 4" ".tmp/tmp2.gtf"
cellranger mkgtf \
    ".tmp/tmp2.gtf" \
    "$filtered_gtf_name" \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA
check_biotypes "After step 4" "$filtered_gtf_name"
log "STEP 4 ended"


log "STEP 5: Building reference index"
cellranger mkref --nthreads $cellranger_nthreads \
    --genome="$mkref_output_name" \
    --fasta="$ref_fasta" \
    --genes="$filtered_gtf_name"
log "STEP 5 ended"

rm -f .tmp/*
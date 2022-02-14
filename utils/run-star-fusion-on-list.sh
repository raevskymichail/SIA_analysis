#!/bin/bash

set -e

## Accept a list selected samples to find according .fastq.gz files in `EXP_DIR` to run STAR-Fusion

###########################
# RUN ON SINGLE-END FASTQ #
###########################

# ./run-star-fusion-on-list.sh Oncobox_BT_samples_for_star_fusion.txt

# SAMPLES_LIST="/home/mraevsky/star-fusion/samples_fusion_candidates_top_1_perc.txt"
export SAMPLES_LIST="/home/mraevsky/BT_samples/Oncobox_BT_samples_for_star_fusion.txt"
# export SAMPLES_LIST=$1
export EXP_DIR="/media/large/volume3/RNAseq"
export N_JOBS=1

export STAR_FUSION_DIR="/home/msorokin/STAR-Fusion"
export CTAT_RESOURCE_LIB="${STAR_FUSION_DIR}/CTAT_RESOURCE_LIB/STAR_Fusion_v1.9/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir"
export STAR_PATH="${STAR_FUSION_DIR}/STAR_versions/STAR-2.7.2b/source/STAR"
# export STAR_FUSION_OUTDIR="${STAR_FUSION_DIR}/OUTDIR"
export STAR_FUSION_OUTDIR="/home/mraevsky/BT_samples/star_fusion_results"
mkdir -p $STAR_FUSION_OUTDIR

run-star-fusion() {
    # SAMPLE_READS="/media/large/volume3/RNAseq/AB_5RNAseq_09062019_FA212/Luc-45_S5_R1_001.fastq.gz"
    SAMPLE_READS=$1
    SAMPLE_NAME="$(basename -- "${SAMPLE_READS%.fastq.gz}")"
    OUTDIR="${STAR_FUSION_OUTDIR}/${SAMPLE_NAME}"
    echo "Run STAR-Fusion on $SAMPLE_READS"

    $STAR_FUSION_DIR/STAR-Fusion --left_fq "$SAMPLE_READS" \
        -O "$OUTDIR" \
        --genome_lib_dir $CTAT_RESOURCE_LIB \
        --STAR_PATH $STAR_PATH \
        --FusionInspector inspect \
        --verbose_level 2
}

export -f run-star-fusion

cat "$SAMPLES_LIST" |
    xargs -I{} find $EXP_DIR -type f -iname *{}* 2>&1 |
    grep ".fastq.gz" |
    parallel -k --jobs $N_JOBS run-star-fusion >log.txt

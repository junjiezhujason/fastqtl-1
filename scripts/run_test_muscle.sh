#!/bin/bash
# load libraries
TDIR=/home/jjzhu/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/boost_1_58_0/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/cnpy/lib/

# input files
TISS=Muscle_Skeletal
IDIR=/share/PI/sabatti/controlled_access_gtex_data/our_analysis
VCF=${IDIR}/genotype_data/gtex.v6.allchr.impute.info04.maf01.hwep1e6.constrvarids.vcf.gz
BED=${IDIR}/phenotype_data/${TISS}.normalized.expr.bed.gz
COV=${IDIR}/covariate_data/${TISS}_Analysis.covariates.txt
SAMP=${IDIR}/sample_list/intersect_ids_${TISS}.txt

ODIR=/scratch/PI/sabatti/controlled_access_data/fastqtl_tmp/${TISS}
PFX=${ODIR}/test

mkdir -p ${PFX}_chunk_999_mtx

PROG=/home/jjzhu/source_code/fastqtl-1/bin/fastQTL
$PROG  --vcf $VCF  --bed $BED --cov $COV \
       --include-samples $SAMP \
       --region 22:17000000-18000000 \
       --window 1e6 \
       --ma-sample-threshold 10 \
       --maf-threshold 0.01 \
       --out ${PFX}_chunk_999.txt.gz \
       --log ${PFX}_chunk_999.log \
       --mtx ${PFX}_chunk_999_mtx

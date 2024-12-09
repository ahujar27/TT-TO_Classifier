#!/bin/bash

TMP_DIR = $1
ref = $STORAGE1/reference/genomes/GRCh37_r75/Homo_sapiens.GRCh37.75.gtf
CPUS = $3
OUTDIR = $4
SAMPNAME = $5
BAMPATH = $6

/opt/subread-2.0.8-Linux-x86_64/bin/featureCounts --tmpDir $TMP_DIR \
         -p \
         -t exon \
         -a $ref \
         -g gene_id \
         -T $CPUS \
         -s 0 \
         -M \
         -O \
         -d 40 \
         -D 1000 \
         -o $OUTDIR/$SAMPNAME.cnt \
         $BAMPATH
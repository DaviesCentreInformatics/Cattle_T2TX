#!/bin/bash -l
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=3-00:00:00
#SBATCH --mem=1000GB

prefix=X
BASE=tuliwu_chrX_94identity_400bp
container=$BASE/lib/CARP_02_w_py.sif

OUTDIR=${BASE}/output
K_INPUT=${BASE}/data/${prefix}.fa
THREADS=24

I_INPUT=${OUTDIR}/${prefix}.gff
I_OUTPUT=${OUTDIR}/${prefix}_krishna.json

G_OUTPUT=${OUTDIR}/${prefix}_krishna.igor.gff

S_INPUT=${OUTDIR}/${prefix}.mfa
S_OUTPUT=${OUTDIR}/${prefix}_krishna.igor.gff

# Load singulairty
module load Singularity
# Toward working directory
cd $BASE


mkdir -p $OUTDIR
cd $OUTDIR

#############################################################################
# 1.2.1 Use krishna to do pairwise alignment between human genome sequences #
#############################################################################
# 1.2.1 output: $OUTDIR/krishna-[date]-[digits].log, $OUTDIR/$prefix.gff, 
singularity exec $container matrix -threads=$THREADS -krishnaflags="-tmp=./ -threads=$THREADS -log -filtid=0.94 -filtlen=400 -target=$K_INPUT" $K_INPUT

###########################################################################
# 1.2.2 Use igor to report repeat feature family groupings in JSON format #
###########################################################################
# 1.2.2 output: $OUTDIR/$prefix_krishna.json
# Unable this when used bundle
singularity exec $container igor -in $I_INPUT -out $I_OUTPUT

#########################################################################
# 1.2.3 Use seqer to generate consensus sequences from genome intervals #
#########################################################################
# 1.2.3 output: $OUTDIR/$prefix_krishna.igor.gff, $OUTDIR/X.mfa, $OUTDIR/$prefix_krishna.igor.gff, $OUTDIR/consensus/
singularity exec $container gffer < $I_OUTPUT > $G_OUTPUT
cat $K_INPUT > $S_INPUT
singularity exec $container seqer -aligner=muscle -dir=consensus -fasta=true -maxFam=100 -subsample=true -minLen=0.95 -threads=24 -ref=$S_INPUT $S_OUTPUT

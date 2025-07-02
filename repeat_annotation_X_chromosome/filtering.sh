#!/bin/bash -l
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --qos=long
#SBATCH --time=7-00:00:00
#SBATCH --mem=100GB

BASE=tuliwu_chrX_94identity_400bp
container=$BASE/lib/CARP_02_w_py.sif

OUTDIR=${BASE}/output
LIBDIR=${BASE}/lib

REPBASE=${LIBDIR}/29.04_RepBase_All.fa
OUR_KNOWN=${LIBDIR}/our_known_reps_20130520.fasta

SPROT_LIB=${LIBDIR}/sprot/uniprot_sprot.fasta
BLAST_INPUT=${OUTDIR}/notKnown.fa
GB_TE_LIB=${LIBDIR}/BlastDB/GB_TE.new
RETROVIRUS_LIB=${LIBDIR}/BlastDB/all_retrovirus

# Load singulairty
module load Singularity
# Toward working directory
cd $BASE

###########################################################
# 1.5.1 Annotate consensus sequences with repeat families #
###########################################################
# Identifying repeats in RepBase by censor
# 1.5.1 output: $OUTDIR/censor.[digits].log, $OUTDIR/ConsensusSequences.fa (from find), $OUTDIR/ConsensusSequences.fa.map, $OUTDIR/ConsensusSequences.fa.aln, $OUTDIR/ConsensusSequences.fa.found, $OUTDIR/ConsensusSequences.fa.idx, $OUTDIR/ConsensusSequences.fa.masked
cd $OUTDIR
echo "====================== Identifying repeats in RepBase... ======================"
find $OUTDIR/consensus -maxdepth 1 -name '[!.]*.fq' -print0 | xargs -r0 cat > $OUTDIR/ConsensusSequences.fa
singularity exec $container censor -bprm cpus=24 -lib $REPBASE -lib $OUR_KNOWN ConsensusSequences.fa



######################################
# 1.5.2 Classify consensus sequences #
######################################
# 1.5.2 output: $OUTDIR/ClassifyConsensusSequences.log, $LIBDIR/ClassifyConsensusSequences.class, $OUTDIR/known.txt, $OUTDIR/partial.txt, $OUTDIR/check.txt, $OUTDIR/notKnown.fa, $OUTDIR/notknown.fa.gff
cd $LIBDIR
echo "======================Processing and classifying consensus sequences...======================"
singularity exec $container javac ClassifyConsensusSequences.java
singularity exec $container java ClassifyConsensusSequences > $OUTDIR/ClassifyConsensusSequences.log 2>&1


##########################
# 1.5.3 Filter sequences #
##########################

# 1.5.3.1 Blast sprot protein
# 1.5.3.1 output: $OUTDIR/sprot/notKnown.fa.spwb, $OUTDIR/notKnown.fa.spwb.gff
cd $LIBDIR/sprot
echo "======================Blastx with Sprot protein======================"
singularity exec $container blastx $SPROT_LIB $BLAST_INPUT -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > $OUTDIR/notKnown.fa.spwb
singularity exec $container python $LIBDIR/wublastx2gff.py $OUTDIR/notKnown.fa.spwb > $OUTDIR/notKnown.fa.spwb.gff


# 1.5.3.2 Blast GB_TE
# 1.5.3.2 output: $OUTDIR/notKnown.fa.tewb, $OUTDIR/notKnown.fa.tewb.gff
cd $LIBDIR/BlastDB
echo "======================Blastx with GB_TE======================"
singularity exec $container blastx $GB_TE_LIB $BLAST_INPUT -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > $OUTDIR/notKnown.fa.tewb
singularity exec $container python $LIBDIR/wublastx2gff.py $OUTDIR/notKnown.fa.tewb > $OUTDIR/notKnown.fa.tewb.gff


# 1.5.3.3 Blast retrovirus
# 1.5.3.3 output: $OUTDIR/notKnown.fa.ervwb, $OUTDIR/notKnown.fa.ervwb.gff
echo "======================tBlastx with retrovirus======================"
singularity exec $container tblastx $RETROVIRUS_LIB $BLAST_INPUT -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > $OUTDIR/notKnown.fa.ervwb
singularity exec $container python $LIBDIR/wublastx2gff.py $OUTDIR/notKnown.fa.ervwb > $OUTDIR/notKnown.fa.ervwb.gff


##########################################################
# 1.5.4 Get protein information from consensus sequences #
##########################################################
# 1.5.4 output: $OUTDIR/GetProteins.log, $LIBDIR/GetProteins.class, $OUTDIR/proteins.txt (a list of families that have been identified as proteins and the proteins they match); $OUTDIR/notKnownNotProtein.fa (a fasta file of the families that were not classified).
cd $LIBDIR
echo "======================Get protein======================"
singularity exec $container javac GetProteins.java
singularity exec $container java GetProteins > $OUTDIR/GetProteins.log 2>&1

#################################################
# 1.5.5 Check for simple sequence repeats (SSR) #
#################################################
# 1.5.5 output: $OUTDIR/phobos.log, $OUTDIR/notKnownNotProtein.phobos
cd $OUTDIR
echo "======================Check SSR======================"
singularity exec $container phobos-linux-gcc4.1.2 -r 7 --outputFormat 0 --printRepeatSeqMode 0 $OUTDIR/notKnownNotProtein.fa 2>&1 | tee $OUTDIR/notKnownNotProtein.phobos > $OUTDIR/phobos.log

#####################################################################
# 1.5.6 Identify the sequences that are SSRs from the phobos output #
#####################################################################
# 1.5.6 output: $OUTDIR/IdentifySSRs.log, $LIBDIR/IdentifySSRs.class, $OUTDIR/SSR.txt
cd $LIBDIR
echo "======================Identify the sequences that are SSRs from the phobos output======================"
singularity exec $container javac IdentifySSRs.java
singularity exec $container java IdentifySSRs > $OUTDIR/IdentifySSRs.log 2>&1

cd $OUTDIR
echo "=== GB_TE and retrovirus data are downloaded by using EDirect. Convert the data format before generating the annotated repeat library... ==="
singularity exec $container perl -pe "s/^>/>gi|GBTE|sp|/g" $LIBDIR/GB_TE.2022-07-14.fa > $OUTDIR/GB_TE.fa
singularity exec $container sed -i 's/ /| /' $OUTDIR/GB_TE.fa

###########################################
# 1.5.7 Generate annotated repeat library #
###########################################
# 1.5.7 output: $OUTDIR/GenerateAnnotatedLibrary.log, $OUTDIR/library/Denovo_TE_Library.fasta
echo "================== Finally, generating annotated repeat library... =================="
mkdir -p $OUTDIR/library
cd $LIBDIR
singularity exec $container javac GenerateAnnotatedLibrary.java
singularity exec $container java GenerateAnnotatedLibrary > $OUTDIR/GenerateAnnotatedLibrary.log 2>&1
echo "================== CARP has completed the annotation. Check your results at output/library/Denovo_TE_Library.fasta! =================="

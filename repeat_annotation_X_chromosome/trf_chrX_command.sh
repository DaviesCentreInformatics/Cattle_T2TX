#!/bin/bash
match=2 (matching weight)
mismatch=7 (mismatching penalty)
indel=7 (indel penalty)
PM=80 (match probability)
PI=5 (indel probability)
minscore=200 (minimum alignment score to report)
maxperiod=2000 (maximum period size to report)
inFile= 38-50Mb, chr X (separated applied)
cd $BASE
trf $inFile $match $mismatch $indel $PM $PI $minscore $maxperiod -d -ngs > trf_chrX_2000.txt

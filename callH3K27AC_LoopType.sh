#!/bin/bash

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_H3K27AC_HICHIP/doc
chipf=()

for type in MEF; do 
	chipf+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

sm=MEF
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--fin $DIN/${sm}_H3K27AC.txt \
--chipf ${chipf[@]} \
--bedout $DIN/${sm}_H3K27AC_10KB_LoopType.txt"
$cmd
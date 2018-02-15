#!/bin/bash

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_H3K27AC_HICHIP/doc
vchip=()

sm_grp=(ESC)
for type in ${sm_grp[@]}; do 
	vchip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

sm=ESC
enh_type=Super
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/${sm}_H3K27AC.txt \
--vchip ${vchip[@]} \
--bedout $DIN/${sm}_H3K27AC_10KB_${enh_type}_LoopType.txt \
--echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
eval $cmd

# sm=MEF
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/${sm}_H3K27AC.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/${sm}_H3K27AC_10KB_LoopType.txt"
# eval $cmd

# sm=CONSTANT
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/Spec_H3K27AC_${sm}.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/Spec_H3K27AC_${sm}_LoopType.txt"
# eval $cmd




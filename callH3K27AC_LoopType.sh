#!/bin/bash

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_H3K27AC_HICHIP/doc
vchip=()

for type in ESC; do 
	vchip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

sm=ESC
enh_type=Normal
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/${sm}_H3K27AC.txt \
--vchip ${vchip[@]} \
--bedout $DIN/${sm}_H3K27AC_10KB_${enh_type}_LoopType.txt \
--echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
echo $cmd

sm=ESC
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/${sm}_H3K27AC.txt \
--vchip ${vchip[@]} \
--bedout $DIN/${sm}_H3K27AC_10KB_LoopType.txt"
echo $cmd



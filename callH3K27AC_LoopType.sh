#!/bin/bash

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_H3K27AC_HICHIP/doc/
vchip=()
sm=ESC

sm_grp=($sm)
for type in ${sm_grp[@]}; do 
	vchip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

enh_type=Super
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/Spec_H3K27AC_${sm}.txt \
--vchip ${vchip[@]} \
--bedout $DIN/Spec_H3K27AC_${sm}_${enh_type}_LoopType.txt \
--echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
eval $cmd

# enh_type=Super
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/Spec_H3K27AC_${sm}.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/Spec_H3K27AC_${sm}_${enh_type}_LoopType.txt \
# --echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
# eval $cmd

# sm=ESC
# enh_type=Super
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/Spec_H3K27AC_${sm}.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/Spec_H3K27AC_${sm}_${enh_type}_LoopType.txt \
# --echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
# eval $cmd

# sm=ESC
# enh_type=Normal
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/${sm}_H3K27AC.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/${sm}_H3K27AC_100KB_${enh_type}_LoopType.txt \
# --echip ~/athena/CHIP/PubData/Whyte/doc/Whyte-${enh_type}EnhGencodeEnh.bed"
# eval $cmd

# sm=ESC
# cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
# --hichip $DIN/${sm}_H3K27AC.txt \
# --vchip ${vchip[@]} \
# --bedout $DIN/${sm}_H3K27AC_100KB_LoopType.txt"
# eval $cmd

sm=CONSTANT
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/Spec_H3K27AC_${sm}.txt \
--vchip ${vchip[@]} \
--bedout $DIN/Spec_H3K27AC_${sm}_LoopType.txt"
eval $cmd

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_KLF4_HICHIP/doc
vchip=()
echip=()
sm=DAY6

sm_grp=($sm)
for type in ${sm_grp[@]}; do
	vchip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/KLF4_${type}/KLF4_${type}_common.bed) 
	echip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip $DIN/${sm}_KLF4.txt \
--vchip ${vchip[@]} \
--bedout $DIN/${sm}_KLF4_50K_H3K27AC_LoopType.txt \
--echip ${echip[@]}"
eval $cmd

SDR=~/athena/ComScripts/RPack/GRFLoop
DIN=~/athena/Andreas_H3K27AC_HICHIP/CMP
vchip=()
sm=ESC

sm_grp=($sm)
for type in ${sm_grp[@]}; do 
	vchip+=(~/athena/CHIP/CHIP_seq/MergeRep/mergepeak/ReProgram/H3K27AC_${type}/H3K27AC_${type}_common.bed)
done

enh_type=TE
cmd="Rscript --no-save --no-restore $SDR/H3K27AC_LoopType.R \
--hichip ${DIN/CMP/doc}/Spec_H3K27AC_${sm}.txt \
--vchip ${vchip[@]} \
--bedout $DIN/Spec_H3K27AC_${sm}_${enh_type}_LoopType.txt \
--echip $DIN/Whyte${enh_type}mm10_Gencode.bed"
eval $cmd



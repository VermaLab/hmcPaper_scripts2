bam=$1
bname="${bam%.*}"

bam clipOverlap --in $bam --out ${bname}_clip.bam
samtools index ${bname}_clip.bam


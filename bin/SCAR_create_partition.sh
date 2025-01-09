#!/bin/bash
# Author: Maria Dalby
# Modified by: Nicolas Alcaraz <nicolas.alcaraz@cpr.ku.dk>
# downloaded from https://github.com/grothlab/SCARseq_Pipeline/blob/main/SCAR_create_partition.sh

BAM_FILE=""
PAIRED_END="true"
INPUT_PAIRED_END="true"
INPUT_BAM_FILE=""
OUT_BASE_DIR="MAPPING_mm10"
CHROM_SIZES="/maps/projects/dan1/people/ngl887/scripts/scar_example/data/mm10.base.chrom.sizes"
STRAND_LIBRARY="reverse"
WIN=1000
RADIUS=30
DRADIUS=30
ZRADIUS=1


usage() {
        echo
        echo ${0##/*}" Usage: SCAR_seq_process.sh [OPTIONS] <bam_file>"
        echo "Script that 1) Splits SCAR-seq <bam_file> by strand."
        echo "            2) Smooths the signal."
        echo "            3) Computes the breakpoints and the RFD scores."
        echo "requires partition_smooth.pl script and bdg2bw program to be in the path"
        echo
        echo "OPTIONS:"
        echo "[-b STRING]   BAM file to analyze (required)"
        echo "[-p STRING]   If SCAR BAM file is paired-end (default: true)"
        echo "[-e STRING]   stranded Input BAM (optional, required for Input Correction)"
        echo "[-i STRING]   If stranded Input BAM is paired-end (optional, default: true)"
        echo "[-s STRING]   If strand library is forward or reverse (default: reverse)"
        echo "[-o STRING]   Directory where all output will be written to"
        echo "[-g STRING]   Chromosome sizes file"
        echo "[-s STRING]   Strande library (reverse or forward, default: reverse)"
        echo "[-w INTEGER]  Window size for smoothing (default 1000)"
        echo "[-r INTEGER]  Radius for smoothing (default 30)"
        echo "[-d INTEGER]  D-radius for smoothing (default 30)"
        echo "[-z INTEGER]  Z-raduus for smoothing (default 1)"
        exit
}

## show usage if '-h' or  '--help' is the first argument or no argument is given
case $1 in
        ""|"-h"|"--help") usage ;;
esac

echo "Running \"`basename $0`\" with parameters \"$*\""

## read the paramters
while getopts b:p:e:i:s:o:g:w:r:d:z: opt; do
        case ${opt} in
                b) BAM_FILE=${OPTARG};;
                p) PAIRED_END=${OPTARG};;
                e) INPUT_BAM_FILE=${OPTARG};;
                i) INPUT_PAIRED_END=${OPTARG};;
                s) STRAND_LIBRARY=${OPTARG};;
                o) OUT_BASE_DIR=${OPTARG};;
                g) CHROM_SIZES=${OPTARG};;
                w) WIN=${OPTARG};;
                r) RADIUS=${OPTARG};;
                d) DRADIUS=${OPTARG};;
                z) ZRADIUS=${OPTARG};;
                *) echo "unrecognized option ${opt}"; usage;;
        esac
done


if [ ! -d ${OUT_BASE_DIR} ]; then
  mkdir ${OUT_BASE_DIR}
  chmod -R g+ws ${OUT_BASE_DIR}
fi

OUT_BAM_DIR="${OUT_BASE_DIR}/bam_files"
OUT_BW_DIR="${OUT_BASE_DIR}/bw_files"
OUT_RFD_DIR="${OUT_BASE_DIR}/rfd_files"

if [ ! -d ${OUT_BAM_DIR} ]; then mkdir ${OUT_BAM_DIR}; chmod -R g+ws ${OUT_BAM_DIR}; fi
if [ ! -d ${OUT_BW_DIR} ]; then mkdir ${OUT_BW_DIR}; chmod -R g+ws ${OUT_BW_DIR}; fi
if [ ! -d ${OUT_RFD_DIR} ]; then mkdir ${OUT_RFD_DIR}; chmod -R g+ws ${OUT_RFD_DIR}; fi


PREFIX=$(basename ${BAM_FILE} .bam)


TEMP_ROOT=${OUT_BASE_DIR}/tmp
if [ ! -d ${TEMP_ROOT} ]; then mkdir ${TEMP_ROOT}; chmod -R g+ws ${TMP_ROOT}; fi

TMP_DIR=$(mktemp -d -p ${TEMP_ROOT})




OUT_F_BAM=${TMP_DIR}/${PREFIX}_F.bam
OUT_R_BAM=${TMP_DIR}/${PREFIX}_R.bam

OUT_F_BG=${OUT_BW_DIR}/${PREFIX}_F.bdg
OUT_R_BG=${OUT_BW_DIR}/${PREFIX}_R.bdg
OUT_F_BW=${OUT_BW_DIR}/${PREFIX}_F.bw
OUT_R_BW=${OUT_BW_DIR}/${PREFIX}_R.bw



PREFIX=$(basename ${BAM_FILE} .bam)_SE
if [[ ${PAIRED_END} == 'true' ]]; then
        PREFIX=$(basename ${BAM_FILE} .bam)_PE
fi
INPUT_PREFIX=$(basename ${INPUT_BAM_FILE} .bam)_SE
if [[ ${INPUT_PAIRED_END} == 'true' ]]; then
  INPUT_PREFIX=$(basename ${INPUT_BAM_FILE} .bam)_PE
fi


INPUT_F_BAM=${TMP_DIR}/${INPUT_PREFIX}_F.bam
INPUT_R_BAM=${TMP_DIR}/${INPUT_PREFIX}_R.bam


INPUT_F_BG=${OUT_BW_DIR}/${INPUT_PREFIX}_F.bdg
INPUT_R_BG=${OUT_BW_DIR}/${INPUT_PREFIX}_R.bdg
INPUT_F_BW=${OUT_BW_DIR}/${INPUT_PREFIX}_F.bw
INPUT_R_BW=${OUT_BW_DIR}/${INPUT_PREFIX}_R.bw



if [[ ! -f ${OUT_F_BAM} ]] || [[ ! -f ${OUT_R_BAM} ]]; then
  echo "Splitting bam file in forward and reverse strands..."
  if [[ ${STRAND_LIBRARY} == 'reverse' ]]
  then
    echo "STRAND LIBRARY detected as reversed"
    samtools view -F 20 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_R_BAM}
    samtools view -f 16 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_F_BAM}
  elif [[ ${STRAND_LIBRARY} == 'forward' ]]
  then
    samtools view -F 20 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_F_BAM}
    samtools view -f 16 -h ${BAM_FILE} | samtools view -Sb -h > ${OUT_R_BAM}
  fi
  samtools index ${OUT_F_BAM}
  samtools index ${OUT_R_BAM}
fi

if [[ -f ${INPUT_BAM_FILE} ]]; then
  if [[ ! -f ${INPUT_F_BAM} ]] || [[ ! -f ${INPUT_R_BAM} ]]; then
    if [[ ${STRAND_LIBRARY} == 'reverse' ]]
    then
      echo "STRAND LIBRARY detected as reversed"
      samtools view -F 20 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_R_BAM}
      samtools view -f 16 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_F_BAM}
    elif [[ ${STRAND_LIBRARY} == 'forward' ]]
    then
      samtools view -F 20 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_F_BAM}
      samtools view -f 16 -h ${INPUT_BAM_FILE} | samtools view -Sb -h > ${INPUT_R_BAM}
    fi
    samtools index ${INPUT_F_BAM}
    samtools index ${INPUT_R_BAM}
  fi
fi



if [[ ! -f ${OUT_F_BW} ]] || [[ ! -f ${OUT_R_BW} ]]; then
  echo "Creating bedgraph files...."
  if [ ${PAIRED_END} == 'false' ]; then
  	genomeCoverageBed -ibam ${OUT_F_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${OUT_F_BG}
  	genomeCoverageBed -ibam ${OUT_R_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${OUT_R_BG}
  else
	  genomeCoverageBed -ibam ${OUT_F_BAM} -bg -pc -g ${CHROM_SIZES} > ${OUT_F_BG}
    genomeCoverageBed -ibam ${OUT_R_BAM} -bg -pc -g ${CHROM_SIZES} > ${OUT_R_BG}
  fi
  echo "Creating bigwig files.... ${OUT_F_BW}"
  bdg2bw ${OUT_F_BG} ${CHROM_SIZES} ${OUT_F_BW}
  bdg2bw ${OUT_R_BG} ${CHROM_SIZES} ${OUT_R_BW}
fi


if [[ -f ${INPUT_BAM_FILE} ]]; then
  if [[ ! -f ${INPUT_F_BW} ]] || [[ ! -f ${INPUT_R_BW} ]]; then
    echo "Creating bedgraph files...."
    if [ ${INPUT_PAIRED_END} == 'false' ]; then
      genomeCoverageBed -ibam ${INPUT_F_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${INPUT_F_BG}
      genomeCoverageBed -ibam ${INPUT_R_BAM} -bg -fs 150 -g ${CHROM_SIZES} > ${INPUT_R_BG}
    else
      genomeCoverageBed -ibam ${INPUT_F_BAM} -bg -pc -g ${CHROM_SIZES} > ${INPUT_F_BG}
      genomeCoverageBed -ibam ${INPUT_R_BAM} -bg -pc -g ${CHROM_SIZES} > ${INPUT_R_BG}
    fi
    echo "Creating bigwig files.... ${INPUT_F_BW}"
    bdg2bw ${INPUT_F_BG} ${CHROM_SIZES} ${INPUT_F_BW}
    bdg2bw ${INPUT_R_BG} ${CHROM_SIZES} ${INPUT_R_BW}
  fi
fi


### Smooth and partion
echo "Smoothing...."
bedtools makewindows -i srcwinnum -g ${CHROM_SIZES} -w ${WIN} -s ${WIN} > ${TMP_DIR}/windows.bed

echo "Splitting windows...."
## Split windows on chromosome
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	awk -v c=${CHR} '$1==c' ${TMP_DIR}/windows.bed > ${TMP_DIR}/${CHR}_windows.bed &
done
wait

echo "Aggregating counts..."
## Aggregate counts in windows
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	bigWigAverageOverBed ${OUT_F_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab &
	bigWigAverageOverBed ${OUT_R_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_R.tab &
done
wait

## CPM normalise strands individually (with pseudocount 1)
cat ${TMP_DIR}/*_windows_F.tab ${TMP_DIR}/*_windows_R.tab > ${TMP_DIR}/tmp_windows.bed
CPM=$(awk '{SUM += ($4+1)} END {print SUM/10^6}' ${TMP_DIR}/tmp_windows.bed)
echo "$CPM"
echo "CPM factors: $CPM"

if [[ -f ${INPUT_BAM_FILE} ]]; then
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	  bigWigAverageOverBed ${INPUT_F_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_input_F.tab &
    bigWigAverageOverBed ${INPUT_R_BW} ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_input_R.tab &
  done
  wait

  cat ${TMP_DIR}/*_windows_input_F.tab ${TMP_DIR}/*_windows_input_R.tab > ${TMP_DIR}/tmp_windows_input.bed
  CPM2=$(awk '{SUM += ($4+1)} END {print SUM/10^6}' ${TMP_DIR}/tmp_windows_input.bed)
  echo "CPM input factors: $CPM2"
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	  echo "Normalizing for chromosome ${CHR}..."
	  awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_F.tab > ${TMP_DIR}/${CHR}_windows_F_CPM.tab
	  awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_R.tab > ${TMP_DIR}/${CHR}_windows_R_CPM.tab
	  awk -v c=${CPM2} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_input_F.tab > ${TMP_DIR}/${CHR}_windows_input_F_CPM.tab
      awk -v c=${CPM2} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_input_R.tab > ${TMP_DIR}/${CHR}_windows_input_R_CPM.tab

      paste ${TMP_DIR}/${CHR}_windows_F_CPM.tab ${TMP_DIR}/${CHR}_windows_input_F_CPM.tab | \

      awk -v c=${CPM} 'BEGIN{OFS="\t"}{if (($4-$10) > 1/c) print $1, $2, $3, $4-$10, $5, $6; else print $1, $2, $3, 1/c, $5, $6}' - > ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab

      paste ${TMP_DIR}/${CHR}_windows_R_CPM.tab ${TMP_DIR}/${CHR}_windows_input_R_CPM.tab | \
            awk -v c=${CPM} 'BEGIN{OFS="\t"}{if (($4-$10) > 1/c) print $1, $2, $3, $4-$10, $5, $6; else print $1, $2, $3, 1/c, $5, $6}' - > ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab
  done
  wait
else
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
	  echo "Normalizing for chromosome ${CHR}..."
	  awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_F.tab > ${TMP_DIR}/${CHR}_windows_F_CPM.tab
	  awk -v c=${CPM} 'BEGIN{OFS="\t"}{print $1, $2, $3, ($4+1)/c, $5, $6}' ${TMP_DIR}/${CHR}_windows_R.tab > ${TMP_DIR}/${CHR}_windows_R_CPM.tab
  done
  wait
fi


for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
        echo "doing for ${CHR}"
        partition_smooth.pl ${TMP_DIR}/${CHR}_windows_F_CPM.tab \
          ${TMP_DIR}/${CHR}_windows_R_CPM.tab \
          ${RADIUS} ${DRADIUS} ${ZRADIUS} > ${TMP_DIR}/${CHR}_windows_RFD.txt &
done
wait



if [[ -f ${INPUT_BAM_FILE} ]]; then
  echo "Calculating RFD scores for Input"
  ## Calculate RFD, derivative and boundary scores
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
          echo "doing for ${CHR}"
          partition_smooth.pl ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab \
            ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab \
            ${RADIUS} ${DRADIUS} ${ZRADIUS} > ${TMP_DIR}/${CHR}_windows_RFD_minusinput.txt &
  done
  wait
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
          echo "doing for ${CHR}"
          partition_smooth.pl ${TMP_DIR}/${CHR}_windows_input_F_CPM.tab \
            ${TMP_DIR}/${CHR}_windows_input_R_CPM.tab \
            ${RADIUS} ${DRADIUS} ${ZRADIUS} > ${TMP_DIR}/${CHR}_windows_RFD_input.txt &
  done
  wait
fi



OUT_NORMAL=${OUT_RFD_DIR}/${PREFIX}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
if [ -f ${OUT_NORMAL} ]; then
  rm ${OUT_NORMAL}
fi

touch ${OUT_NORMAL}
for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
  paste ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab \
    ${TMP_DIR}/${CHR}_windows_R.tab ${TMP_DIR}/${CHR}_windows_F_CPM.tab \
      ${TMP_DIR}/${CHR}_windows_R_CPM.tab ${TMP_DIR}/${CHR}_windows_RFD.txt | \
       awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT_NORMAL}
done
wait

OUT_BASE=$(basename ${OUT_NORMAL} .txt)
OUT_RFD_BW=${OUT_BW_DIR}/${OUT_BASE}_RFD.bw
OUT_CPM_F_BW=${OUT_BW_DIR}/${OUT_BASE}_CPM_F.bw
OUT_CPM_R_BW=${OUT_BW_DIR}/${OUT_BASE}_CPM_R.bw
OUT_RFD_BDG=${TMP_DIR}/${OUT_BASE}_RFD.bdg
OUT_CPM_F_BDG=${TMP_DIR}/${OUT_BASE}_CPM_F.bdg
OUT_CPM_R_BDG=${TMP_DIR}/${OUT_BASE}_CPM_R.bdg

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$9}' ${OUT_NORMAL} | \
  sort -k1,1 -k2,2n - > ${OUT_RFD_BDG}
bedGraphToBigWig ${OUT_RFD_BDG} ${CHROM_SIZES} ${OUT_RFD_BW}

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$6}' ${OUT_NORMAL} | \
  sort -k1,1 -k2,2n - > ${OUT_CPM_F_BDG}
bedGraphToBigWig ${OUT_CPM_F_BDG} ${CHROM_SIZES} ${OUT_CPM_F_BW}

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' ${OUT_NORMAL} | \
  sort -k1,1 -k2,2n - > ${OUT_CPM_R_BDG}
bedGraphToBigWig ${OUT_CPM_R_BDG} ${CHROM_SIZES} ${OUT_CPM_R_BW}

gzip ${OUT_NORMAL}

if [[ -f ${INPUT_BAM_FILE} ]]; theth.pl n
  ## Gather output
  OUT_MINUSINPUT=${OUT_RFD_DIR}/${PREFIX}_SCARminusinput_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
  if [ -f ${OUT_MINUSINPUT} ]; then
    rm ${OUT_MINUSINPUT}
  fi

  touch ${OUT_MINUSINPUT}
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
    paste ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_F.tab \
    ${TMP_DIR}/${CHR}_windows_R.tab ${TMP_DIR}/${CHR}_windows_F_CPM_minusinput.tab \
    ${TMP_DIR}/${CHR}_windows_R_CPM_minusinput.tab ${TMP_DIR}/${CHR}_windows_RFD_minusinput.txt | \
      awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT_MINUSINPUT}
  done
  wait

  OUT_MINPUT_BASE=$(basename ${OUT_MINUSINPUT} .txt)
  OUT_RFD_BW=${OUT_BW_DIR}/${OUT_MINPUT_BASE}_RFD.bw
  OUT_CPM_F_BW=${OUT_BW_DIR}/${OUT_MINPUT_BASE}_CPM_F.bw
  OUT_CPM_R_BW=${OUT_BW_DIR}/${OUT_MINPUT_BASE}_CPM_R.bw
  OUT_RFD_BDG=${TMP_DIR}/${OUT_MINPUT_BASE}_RFD.bdg
  OUT_CPM_F_BDG=${TMP_DIR}/${OUT_MINPUT_BASE}_CPM_F.bdg
  OUT_CPM_R_BDG=${TMP_DIR}/${OUT_MINPUT_BASE}_CPM_R.bdg

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$9}' ${OUT_MINUSINPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_RFD_BDG}
  bedGraphToBigWig ${OUT_RFD_BDG} ${CHROM_SIZES} ${OUT_RFD_BW}

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$6}' ${OUT_MINUSINPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_CPM_F_BDG}
  bedGraphToBigWig ${OUT_CPM_F_BDG} ${CHROM_SIZES} ${OUT_CPM_F_BW}

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' ${OUT_MINUSINPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_CPM_R_BDG}
  bedGraphToBigWig ${OUT_CPM_R_BDG} ${CHROM_SIZES} ${OUT_CPM_R_BW}



  OUT_INPUT=${OUT_RFD_DIR}/${INPUT_PREFIX}_smooth_results_w${WIN}_s${RADIUS}_d${DRADIUS}_z${ZRADIUS}.txt
  if [ -f ${OUT_INPUT} ]; then
    rm ${OUT_INPUT}
  fi

  touch ${OUT_INPUT}
  for CHR in $( cut -f 1 ${CHROM_SIZES} | sort | uniq ); do
    paste ${TMP_DIR}/${CHR}_windows.bed ${TMP_DIR}/${CHR}_windows_input_F.tab \
      ${TMP_DIR}/${CHR}_windows_input_R.tab ${TMP_DIR}/${CHR}_windows_input_F_CPM.tab \
      ${TMP_DIR}/${CHR}_windows_input_R_CPM.tab ${TMP_DIR}/${CHR}_windows_RFD_input.txt | \
        awk '$8>0 || $14>0' | cut -f -3,8,14,20,26,30- >> ${OUT_INPUT}
  done
  wait

  OUT_INPUT_BASE=$(basename ${OUT_INPUT} .txt)
  OUT_RFD_BW=${OUT_BW_DIR}/${OUT_INPUT_BASE}_RFD.bw
  OUT_CPM_F_BW=${OUT_BW_DIR}/${OUT_INPUT_BASE}_CPM_F.bw
  OUT_CPM_R_BW=${OUT_BW_DIR}/${OUT_INPUT_BASE}_CPM_R.bw
  OUT_RFD_BDG=${TMP_DIR}/${OUT_INPUT_BASE}_RFD.bdg
  OUT_CPM_F_BDG=${TMP_DIR}/${OUT_INPUT_BASE}_CPM_F.bdg
  OUT_CPM_R_BDG=${TMP_DIR}/${OUT_INPUT_BASE}_CPM_R.bdg

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$9}' ${OUT_INPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_RFD_BDG}
  bedGraphToBigWig ${OUT_RFD_BDG} ${CHROM_SIZES} ${OUT_RFD_BW}

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$6}' ${OUT_INPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_CPM_F_BDG}
  bedGraphToBigWig ${OUT_CPM_F_BDG} ${CHROM_SIZES} ${OUT_CPM_F_BW}

  awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$7}' ${OUT_INPUT} | \
    sort -k1,1 -k2,2n - > ${OUT_CPM_R_BDG}
  bedGraphToBigWig ${OUT_CPM_R_BDG} ${CHROM_SIZES} ${OUT_CPM_R_BW}

  gzip ${OUT_MINUSINPUT} ${OUT_INPUT}
fi
rm -rf ${TMP_DIR}

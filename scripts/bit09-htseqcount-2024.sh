#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest, 2024
# Usage		: bash bit09-htseqcount.sh /home/pacoh/hisat2_filtered/ 
#             /data/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf
#             /home/pacoh/htseqcount/ CDS gene_id PE NO intersection-strict 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this HTSeqCount script requires 9 parameters:\n
		1. Path of the folder with the mapping files to run HTSeqCount on.\n
        2. Path of the GTF file.\n
		3. Path of the output folder.\n
		4. Feature type to use from GTF file (e.g. CDS, exon).\n
		5. Attribute to use as label for counting (e.g. gene_id, gene_name, transcript_id).\n
        6. Specify library layout either PE (paired end) or SE (single end).\n
        7. Specify strandedness either YES (HTSeq default) or NO  (or REVERSE).\n
        8. Mode either intersection-strict or union (HTSeq default).\n
		9. Number of threads to use.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 9 ]; then
	usage
fi
################################################################################
# INPUT FOLDER
################################################################################
inputFolder=$1;
# Remove trailing slash if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi
################################################################################
# OUTPUT FOLDER (CREATE IF NOT EXISTS)
################################################################################
outFolder=$3;
# Remove trailing slash if this is last char
len=${#outFolder};
lastPos=$(expr $len - 1);
lastChar=${outFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	outFolder=${outFolder:0:$lastPos};
fi
mkdir -p ${outFolder}
################################################################################
# RUN HTSEQCOUNT
################################################################################
# GTF annotation arguments
pathGTF=$2;
feature=$4;
attribute=$5;
# Library layout argument
layout=$6;
if [[ $layout == "PE" ]]; then
    layout="-r pos";
elif [[ $layout == "SE" ]]; then
    layout="";
else
    echo -e "Invalid layout specified. Please specify either PE or SE.\n";
    exit 1;
fi
# Strandedness argument
strandedness=$7;
if [[ $strandedness == "YES" || $strandedness == "Yes" || $strandedness == "yes" ]]; then
    strandedness="yes";
elif [[ $strandedness == "NO" || $strandedness == "No" || $strandedness == "no" ]]; then
    strandedness="no";
elif [[ $strandedness == "REVERSE" || $strandedness == "Reverse" || $strandedness == "reverse" ]]; then
    strandedness="reverse";
else
    echo -e "Invalid strandedness specified. Please specify either YES or NO (or REVERSE).\n";
    exit 1;
fi
# Mode argument
mode=$8;
if [[ $mode == "intersection-strict" ]]; then
    mode="intersection-strict";
elif [[ $mode == "union" || $mode == "Union" || $mode == "UNION" ]]; then
    mode="union";
else
    echo -e "Invalid mode specified. Please specify either intersection-strict or union.\n";
    exit 1;
fi
# Threads
threads=$9;
# Part of input/file name
tmpPart='"{}"';
tmpPart2='${IN}';
tmpPart3='$(basename ${IN%.*})';
# Compose command
htseqCommand="find ${inputFolder} -name '*_sorted.bam'";
htseqCommand="$htseqCommand | xargs --max-procs=${threads} -I {} sh -c 'IN=$tmpPart;";
htseqCommand="$htseqCommand htseq-count -f bam -m ${mode}";
htseqCommand="$htseqCommand ${layout} --stranded=${strandedness}";
htseqCommand="$htseqCommand -a 10 -t ${feature} -i ${attribute} ${tmpPart2} ${pathGTF}";
htseqCommand="$htseqCommand > ${outFolder}/counts_${tmpPart3}.txt';";
# Show command
echo -e "$htseqCommand";
# Execute
output=$(eval $htseqCommand);
# Show output
echo -e "$output\n";
################################################################################
#!/bin/bash
################################################################################
# Author	: Paco Hulpiau - Howest
# Usage		: ./bit09-filtering-samtools.sh /home/user/tophat2/ 20
#             /home/user/tophat2_filtered/ 4
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="This filtering script requires 4 parameters:\n
		1. Path of the folder with mapping files to run samtools on.\n
        2. MAPping Quality value to filter (use default 20).\n
        3. Path of the output folder.\n
		4. Number of threads to use (max. 8 on BIT server!).\n\n
        Run samtools (on TopHat2 output).";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 4 ]; then
	usage
fi
################################################################################
# INPUT FOLDER CONTAINING .fastq.gz files
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
# RUN SAMTOOLS
################################################################################
# MAPping Quality and number of threads to use
mapq=$2;
threads=$4;
# Part of input/file name found with "find" command
tmpPart='"{}"';
tmpPart2='${IN}';
tmpPart3='$(basename ${IN%.*})';
################################################################################
# SAMTOOLS VIEW:
# samtools view -bq 20 accepted_hits.bam > accepted_hits_filtered.bam
# --> filter alignment input to BAM output format (-b option) 
#     with minimum required MAPping Quality (MAPQ) of 20 (-q option)
echo -e "\n### SAMTOOLS VIEW ###";
viewCommand="find ${inputFolder} -maxdepth 1 -type d -not -path ${inputFolder}";
viewCommand="$viewCommand | xargs --max-procs=${threads} -I {} sh -c 'IN=$tmpPart;";
viewCommand="$viewCommand samtools view -bq ${mapq} ${tmpPart2}/accepted_hits.bam";
viewCommand="$viewCommand > ${outFolder}/${tmpPart3}_filtered.bam';";
# Show command
echo -e "$viewCommand";
# Execute command
outputViewCommand=$(eval $viewCommand);
echo -e "$outputViewCommand\n";
################################################################################
# SAMTOOLS SORT:
# samtools sort accepted_hits_filtered.bam accepted_hits_filtered_sorted
echo -e "\n### SAMTOOLS SORT ###";
sortCommand="find ${outFolder} -name '*_filtered.bam'";
sortCommand="$sortCommand | xargs --max-procs=${threads} -I {} sh -c 'IN=$tmpPart;";
sortCommand="$sortCommand samtools sort ${tmpPart2} ${outFolder}/${tmpPart3}_sorted';";
echo -e "$sortCommand";
outputSortCommand=$(eval $sortCommand);
echo -e "$outputSortCommand\n";
################################################################################
# SAMTOOLS INDEX:
# samtools index accepted_hits_filtered_sorted.bam
echo -e "\n### SAMTOOLS INDEX ###";
indexCommand="find ${outFolder} -name '*_filtered_sorted.bam'";
indexCommand="$indexCommand | xargs --max-procs=${threads} -I {} sh -c 'IN=$tmpPart;";
indexCommand="$indexCommand samtools index ${tmpPart2}';";
echo -e "$indexCommand";
outputIndexCommand=$(eval $indexCommand);
echo -e "$outputIndexCommand\n";
################################################################################
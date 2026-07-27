#!/bin/bash
################################################################################
# Author	: Saskia Alexander - Howest
# Usage		: ./bit09-search-trimmomatic-adapters.sh <sequence>
################################################################################
# VALIDATE INPUT
################################################################################
function usage(){
	errorString="Running this script requires 1 parameter:\n
		1. The adapter or primer sequence to search for.\n\n
        Searches this sequence and its shorter subsequences (fewer bases on\n
        both sides) against the Trimmomatic combined adapters file and reports\n
        the longest matches in FASTA format.\n";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 1 ]; then
	usage
fi
################################################################################
# SETTINGS AND INPUT
################################################################################
# Trimmomatic combined adapters file to search in
adapterFile="/opt/Trimmomatic-0.39/adapters/adapters.fa";
# Shortest subsequence length that is still reported as a match
minLength=10;
# Query sequence exactly as given
query="$1";
queryLength=${#query};
################################################################################
# LOAD ADAPTER ENTRIES (read the adapters file once into memory)
################################################################################
declare -a entryHeader entrySeq;
header="";
while IFS= read -r line || [ -n "$line" ]; do
	case "$line" in
		\>*) header="$line" ;;
		"" ) : ;;
		*  ) entryHeader+=("$header");
		     entrySeq+=("$line") ;;
	esac
done < ${adapterFile};
nEntries=${#entrySeq[@]};
################################################################################
# SEARCH
# Try the full query first, then progressively shorter subsequences, and stop
# at the longest length that produces at least one match. Each adapter entry
# is reported only once.
################################################################################
declare -A seen;
matchCount=0;
matchOutput="";
# Outer loop: subsequence length k, from the full query length down to minLength
for (( k=queryLength; k>=minLength; k-- )); do
	# Inner loop: start position i, sliding the window of length k across the query
	for (( i=0; i+k<=queryLength; i++ )); do
		# Cut out the subsequence of length k starting at position i
		subSeq=${query:i:k};
		# Compare this subsequence against every adapter entry
		for (( j=0; j<nEntries; j++ )); do
			# Keep the entry if its sequence contains the subsequence
			if [[ "${entrySeq[j]}" == *"$subSeq"* ]]; then
				# Unique key per adapter entry (header + sequence) to avoid duplicates
				matchKey="${entryHeader[j]}|${entrySeq[j]}";
				# Report each adapter entry only the first time it is matched
				if [ -z "${seen[$matchKey]}" ]; then
					seen[$matchKey]=1;
					matchCount=$(expr $matchCount + 1);
					matchOutput="${matchOutput}Match length (bp): ${k}\n";
                    matchOutput="${matchOutput}Matched subsequence: ${subSeq}\n";
					matchOutput="${matchOutput}Matched adapter: ${entryHeader[j]}\n";
					matchOutput="${matchOutput}                 ${entrySeq[j]}\n\n";
				fi
			fi
		done
	done
	# Stop at the first (longest) length that produces at least one match
	if [ $matchCount -gt 0 ]; then
		break;
	fi
done
################################################################################
# OUTPUT
################################################################################
echo "Total matches found: ${matchCount}";
echo "";
echo -e "${matchOutput}";
################################################################################
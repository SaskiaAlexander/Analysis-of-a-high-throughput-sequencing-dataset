#!/bin/bash

#######################################################
# Author	: Paco Hulpiau & Mathias Verbeke - Howest #
# Usage		: ./bit09-merge-srr.sh -h                 #
#######################################################

# Define a function to display usage information for the script
usage() {
    # Print usage instructions for the script
    echo "Usage: $0 -i <input_folder> -o <output_folder> -n <merge_number> [-p] [-h]"
    # Explain the purpose of each argument
    echo "  -i: Input folder containing FASTQ files (required)"
    echo "  -o: Output folder for merged files (required)"
    echo "  -n: Number of files to merge in each group (required)"
    echo "  -p: Flag to indicate paired-end data (optional)"
    echo "  -h: Display usage information"
    # Exit the script after displaying usage
    exit 1
}

# Initialize variables to store input folder, output folder, merge number, and pairing flag
input_folder=""
output_folder=""
merge_number=""
is_paired=false

# Parse command-line arguments using the getopts function
while getopts ":i:o:n:ph" opt; do
    case ${opt} in
        # Handle the -i flag for input folder
        i )
            input_folder=$OPTARG
            ;;
        # Handle the -o flag for output folder
        o )
            output_folder=$OPTARG
            ;;
        # Handle the -n flag for the number of files to merge
        n )
            merge_number=$OPTARG
            ;;
        # Handle the -p flag to indicate paired-end data
        p )
            is_paired=true
            ;;
        # Handle the -h flag to display usage information
        h )
            usage
            ;;
        # Handle invalid options
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        # Handle missing arguments for options
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Validate that all required arguments are provided
if [[ -z "$input_folder" ]] || [[ -z "$output_folder" ]] || [[ -z "$merge_number" ]]; then
    # Print an error message if required arguments are missing
    echo "Error: Missing required arguments" 1>&2
    usage
fi

# Check if the input folder exists
if [[ ! -d "$input_folder" ]]; then
    # Print an error message if the input folder does not exist
    echo "Error: Input folder does not exist" 1>&2
    usage
fi

# Check if the input folder and output folder are the same
if [[ "$input_folder" == "$output_folder" ]]; then
    # Print an error message if the input and output folders are the same
    echo "Error: Input and output folders cannot be the same" 1>&2
    usage
fi

# Create the output folder if it does not already exist
if [[ ! -d "$output_folder" ]]; then
    mkdir -p "$output_folder"
fi

# Define a function to merge single-end FASTQ files
merge_single_end() {
    # Extract arguments into local variables for readability
    local input_folder=$1
    local output_folder=$2
    local merge_number=$3

    # find all .fastq.gz files in the input folder
    local files=($(find "$input_folder" -maxdepth 1 -type f -name "*.fastq.gz" | sort))

    # Check if any FASTQ files were found
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No FASTQ files found in $input_folder" 1>&2
        exit 1
    fi

    # Check if the number of files is divisible by the merge number
    if [[ $(( ${#files[@]} % merge_number )) -ne 0 ]]; then
        echo "Error: Number of files in $input_folder is not divisible by $merge_number" 1>&2
        exit 1
    fi

    # Loop through the files in groups of the specified merge number
    for (( i=0; i<${#files[@]}; i+=merge_number )); do
        # Extract the group of files to merge
        local group=("${files[@]:i:merge_number}")

        # Define the output file name based on the first file in the group
        local output_file="$output_folder/$(basename "${group[0]}" .fastq.gz)_merged.fastq.gz"

        # Merge the group of files into a single output file
        zcat "${group[@]}" | gzip > "$output_file"

        # Print a message indicating the files that were merged
        echo "Merged files: ${group[@]} into $output_file"
    done
}

# Define a function to merge paired-end FASTQ files
merge_paired_end() {
    # Extract arguments into local variables for readability
    local input_folder=$1
    local output_folder=$2
    local merge_number=$3

    # find all .fastq.gz files in the input folder
    local files=($(find "$input_folder" -maxdepth 1 -type f -name "*.fastq.gz" | sort))

    # Check if any FASTQ files were found
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Error: No FASTQ files found in $input_folder" 1>&2
        exit 1
    fi

    # Find all _1 and _2 FASTQ files in the input folder
    local files_1=($(find "$input_folder" -maxdepth 1 -type f -name "*_1.fastq.gz" | sort))
    local files_2=($(find "$input_folder" -maxdepth 1 -type f -name "*_2.fastq.gz" | sort))

    # Ensure that the number of _1 and _2 files are equal
    if [[ ${#files_1[@]} -ne ${#files_2[@]} ]]; then
        echo "Error: Unequal number of _1 and _2 files" 1>&2
        exit 1
    fi

    # Check if the number of _1 and _2 files is divisible by the merge number
    if [[ $(( ${#files_1[@]} % merge_number )) -ne 0 ]]; then
        echo "Error: Number of _1 and _2 files is not divisible by $merge_number" 1>&2
        exit 1
    fi

    # Loop through the files in groups of the specified merge number
    for (( i=0; i<${#files_1[@]}; i+=merge_number )); do
        # Extract the group of files to merge
        local group_1=("${files_1[@]:i:merge_number}")
        local group_2=("${files_2[@]:i:merge_number}")

        # Define the output file names based on the first file in the group
        local base=$(basename "${group_1[0]}" _1.fastq.gz)
        local output_file_1="$output_folder/${base}_merged_1.fastq.gz"
        local output_file_2="$output_folder/${base}_merged_2.fastq.gz"

        # Merge the group of files into a single output file for each pair
        zcat "${group_1[@]}" | gzip > "$output_file_1"
        zcat "${group_2[@]}" | gzip > "$output_file_2"

        # Print a message indicating the files that were merged
        echo "Merged files: ${group_1[@]} into $output_file_1"
        echo "Merged files: ${group_2[@]} into $output_file_2"
    done
}

# Check if the data is paired-end or single-end and call the appropriate merge function
if [[ "$is_paired" == true ]]; then
    merge_paired_end "$input_folder" "$output_folder" "$merge_number"
else
    merge_single_end "$input_folder" "$output_folder" "$merge_number"
fi

# Exit the script with a success status
echo "Merge completed successfully"
exit 0
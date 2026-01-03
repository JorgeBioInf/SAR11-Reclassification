#!/bin/bash
# Executes prodigal for all genomes in a given directory
# Returns annotated genomes named as "anotated_{genome_name}.fa"

# Check arguments
if [ $# -lt 2 ] || [ $# -gt 2 ] ; then
        echo "Error: 2 arguments are needed ($# given)"
        echo "Use: prodigal_analysis.sh genomes_dir output_dirname"
        exit
fi

# Check if directories exist
if [ ! -d $1 ] ; then
        echo "Error: directory $1 does not exist"
        exit
fi


if [ ! -d $2 ] ; then
	echo "Creating directory $2 ..."
	mkdir $2
else
	echo "Directory $2 already exists! Please, select other name"
	exit
fi


GENOMES=$(find $1 -type f -name "*.fa")

for GENOME in $GENOMES ; do
	NAME=$(basename $GENOME)
	echo "Processing $NAME ..."
	prodigal -i $GENOME -d "$2/anotated_$NAME"
done

echo "Prodigal analysis finished!"
echo "Results can be found in $2"

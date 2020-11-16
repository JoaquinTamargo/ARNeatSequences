#! /bin/bash

## The next message will appear if no imputs are given.

if [ $# -eq 0 || $# -gt 1 ]
then
        echo "Usage: bash ARNeatSequences.sh <params>"
        echo ""
        echo "params: Input file including the parameters"
        exit
fi

## In any other case, the script will print the experiment name and number of
## samples specified in the parameters file.

PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }')

echo "Working directory: $WD"

EXP=$(grep experiment_name $PARAMS | awk '{ print $2 }')

echo "Experiment name: $EXP"

NSAM=$(grep number_samples: $PARAMS | awk '{ print $2 }')

echo "Number of samples: $NSAM"


## Once these parameters have been read, workspace is set:

cd $WD
mkdir $EXP
cd $EXP
mkdir genome annotation samples results
cd samples

i=1
while [ $i -le $NSAM ]
do
	mkdir sample_$i
	((i++))
done

## Everythings is looking good! Now, let's copy fastq, genome and annotation files
## and move them to the directory previously set; i.e. fastq to samples, genome
## to genome and annotation to annotation file.

FASTQD=$(grep fastq_directory: $PARAMS | awk '{ print $2 }')

echo "FASTQ directory: $FASTQD"

GEND=$(grep genome_directory: $PARAMS | awk '{ print $2 }')

echo "Genome directory: $GEND"

ANND=$(grep annotation_directory: $PARAMS | awk '{ print $2}')

echo "Annotation directory: $ANND"

mv $GEND/*.fa $EXP/genome
mv $ANND/*.gtf $EXP/annotation

i=1
while [ $i -le $NSAM ]
then
	echo $i
	((i++))
done




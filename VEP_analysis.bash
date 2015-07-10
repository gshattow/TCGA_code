#!/bin/bash   

OPTIND=1 

usage()
{
cat << EOF
usage: $0 options

This script run the test1 or test2 over a machine.

OPTIONS:
   -h		Show this message
   -d		directory
   -i		infile
   -c		chromosome
   -s		sift
   -l		catalogue
   -r		n_runs
EOF
}

DIR=./
INFILE=output_
CHROMOSOME=''
SIFTSCORE=0.05
CATALOGUE=HGNC
NRUNS=100
while getopts “h:d:i:c:s:l:r:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         d)
         	 DIR=$OPTARG
        	 ;;
         i)
             INFILE=$OPTARG
             ;;
         c)
             CHROMOSOME=$OPTARG
             ;;
         s)
             SIFTSCORE=$OPTARG
             ;;
         l)
             CATALOGUE=$OPTARG
             ;;
         r)
         	 NRUNS=$OPTARG
         	 ;;
         ?)
             usage
             exit
             ;;
     esac
done


time for i in 1;
do
	datafile=$DIR$INFILE$CHROMOSOME
	outfile=$INFILE$CHROMOSOME'_SIFT'
	echo "moving "$datafile" to "$outfile
	grep -m 1 '^\#[A-Z]' $datafile > $outfile
	grep SIFT $datafile >> $outfile
	
	
# create params.py file
	echo "chrom = '"$CHROMOSOME"'" > params.py				
	echo "datafile = '"$INFILE$CHROMOSOME"_SIFT'" >> params.py	
	echo "sift = "$SIFTSCORE >> params.py
	echo "catalogue = '"$CATALOGUE"'" >> params.py
	echo "n_runs = "$NRUNS"" >> params.py

	python vep_parser.py

	python vep_results.py

	python vep_figures.py
	
done
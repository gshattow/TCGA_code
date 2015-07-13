#!/bin/bash   

OPTIND=1 

usage()
{
cat << EOF
usage: $0 options

This script runs the VEP analysis suite.

OPTIONS:
   -h		Show this message
   -d		directory (e.g. ../VEP_output/)
   -i		infile prefix (e.g. VEP_output_)
   -c		chromosome (e.g. 3)
   -s		maximum sift score (e.g. 0.05)
   -l		catalogue of gene names (e.g. HGNC)
   -r		n_runs for individual pairings (e.g. 100)
EOF
}

DIR=../VEP_output/
INFILE=VEP_output_
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
	echo $DIR, $INFILE
	datafile=$DIR$INFILE$CHROMOSOME
	outfile=$INFILE$CHROMOSOME'_SIFT'
	echo "moving mutations with SIFT scores < "$SIFTSCORE
	echo "from "$datafile" to ./"$outfile
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
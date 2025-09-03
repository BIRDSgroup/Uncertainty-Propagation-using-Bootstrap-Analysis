#!/bin/bash
#$ -cwd
#$ -o outputfile.txt
#$ -e errorfile.txt

tissue=$1
B=$2
cores=50

echo $tissue $B $n
mkdir $tissue

#SEED GENERATION
echo "Seed Generation + Gene Expression Extraction"
mkdir seeds
Rscript seed.R $tissue $B

#ORIGINAL DEGREE COMPUTATION
echo "Computing Original Degree"
python network.py $tissue

#BOOTSTRAPPING
echo 'Bootstrapping'
mkdir degree
mkdir pagerank

counter=0
for b in `seq 1 $B`
do
    echo $b
    time python compute.py $tissue $((b-1)) &
    ((counter++)); ((counter % cores == 0)) && wait
done

##Merge centralities into one file
echo "merging"
time python merge.py $tissue $B

#DELETING INTERMEDIATE FILES
rm -r seeds
rm -r degree
rm -r pagerank

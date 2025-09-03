#!/bin/bash
#$ -cwd
#$ -o outputfile.txt
#$ -e errorfile.txt

for centrality in 'degree' 'pagerank'; 
do
	echo $centrality
	python "code 2.1.py" $centrality
	python "code 2.2.py" $centrality
done

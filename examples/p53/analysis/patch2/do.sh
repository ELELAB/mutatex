export BASEDIR=../../

export MUTLIST=${BASEDIR}/mutation_list.txt
export RESULTS=${BASEDIR}/results/interface_ddgs/3kmd_model0_Repair/A-D
export PDB=3kmd_195-204_221-227.pdb

ddg2heatmap -p $PDB -d $RESULTS -l $MUTLIST -s 50 -F Arial -n -3 -x 5 -c plasma -b labels.csv
ddg2distribution -p $PDB -l $MUTLIST -d $RESULTS -T stem -F Arial -s 500 -n -3 -x 22 -c viridis -b labels.csv


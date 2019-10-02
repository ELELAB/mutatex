export BASEDIR=../../

export MUTLIST=${BASEDIR}/mutation_list.txt
export RESULTS=${BASEDIR}/results/interface_ddgs/3kmd_model0_checked_Repair/A-B/
export PDB=./3kmd_173-184_238-247.pdb

ddg2heatmap -p $PDB -d $RESULTS -l $MUTLIST -s 50 -F Arial -n -3 -x 5 -c plasma -b labels.csv 
ddg2distribution -p $PDB -l $MUTLIST -d $RESULTS -T violin -F Arial -s 500 -n -3 -x 22 -c viridis -b labels.csv -v


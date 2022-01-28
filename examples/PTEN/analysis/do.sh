export MUTLIST=../mutation_list.txt
export RESULTS=../results/mutation_ddgs/1D5R_noheatm_model0_checked_Repair/
export PDB=./1D5R_clean_shorter.pdb

ddg2heatmap -p $PDB -d $RESULTS  -l $MUTLIST -s 50 -F Arial -n -3 -x 5 -c plasma -b labels_1D5R_clean_shorter.csv
ddg2distribution -p $PDB -l $MUTLIST -d $RESULTS -T box -F Arial -s 50 -n -3 -x 22 -b labels_1D5R_clean_shorter.csv

export PDB=../1D5R_noheatm.pdb

ddg2distribution -p $PDB -l $MUTLIST -d $RESULTS -T scatter -F Arial -s 500 -n -3 -x 22 -c viridis -u 10 -b labels_1D5R_clean.csv

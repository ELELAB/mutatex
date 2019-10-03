SYSTEMS=(CRK GRB2 SHP2)

for s in ${SYSTEMS[@]}; do
	
	PDB=./${s}_small.pdb
	MUTLIST=../${s}/mutation_list.txt
	RESULTS=../${s}/results/interface_ddgs/*Repair/A-*

	ddg2heatmap -p $PDB -d $RESULTS -l $MUTLIST -s 50 -F Arial -n -3 -x 5 -c plasma -o heatmap_${s}.pdf -b labels_${s}.csv
	ddg2distribution -p $PDB -l $MUTLIST -d $RESULTS -T box -F Arial -s 500 -n -3 -x 10 -o ${s}.pdf -b labels_${s}.csv

	for r in $(cut -d , labels_${s}.csv -f1 | grep -v Residue_name | tr '\n' ' '); do
		ddg2histo -l $MUTLIST -d $RESULTS -p ${s}_small.pdb -r ${r} -n -3 -x 15 -b labels_${s}.csv -F Arial -o histo_${s}.pdf
	done
done


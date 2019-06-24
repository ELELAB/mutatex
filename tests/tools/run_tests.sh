#!/bin/bash

export BASEDIR="../foldxsuite5/basic"
export MUTLIST=$BASEDIR"/mutation_list.txt"
export PDB=$BASEDIR"/2klz_1.pdb"
export RESULTS=$BASEDIR"/results/mutation_ddgs/2klz_1_model0_checked_Repair/"

export BASEOPT="-p $PDB"

export tools=(../../ddg2density ../../ddg2dg ../../ddg2distribution ../../ddg2histo ../../ddg2labels ../../ddg2logo ../../ddg2matrix ../../ddg2pdb ../../ddg2xlsx)
export options=("-l $MUTLIST -d $RESULTS"
"-d $RESULTS"
"-l $MUTLIST -d $RESULTS -T box"
"-l $MUTLIST -d $RESULTS -r QA11"
""
"-l $MUTLIST -d $RESULTS"
"-l $MUTLIST -d $RESULTS"
"-l $MUTLIST -d $RESULTS" 
"-l $MUTLIST -d $RESULTS")

for t in "${!tools[@]}"; do
	echo -n "now running " ${tools[$t]} $BASEOPT ${options[$t]}
	${tools[$t]} $BASEOPT ${options[$t]}
	if [[ $? -ne 0 ]]; then
		echo " ... FAILED"
	else
		echo " ... PASSED"
	fi
done




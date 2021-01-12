#!/bin/bash

export FOLDX_BINARY=/usr/local/foldx4_2020/foldx
export FOLDX_ROTABASE=/usr/local/foldx4_2020/rotabase.txt
export FOLDX_VERSION=suite4
export NP=8
export tests=(basic basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_ptms basic_skip_repair basic_targz basic_multimodel_skip_repair interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers interaction_skip_repair selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers selfmutate_skip_repair basic_poslist basic_multimers_poslist interaction_multimers_poslist selfmutate_interaction selfmutate_poslist)

runnables=()

if [ $# -gt 0 ]; then
	runnables=$@
else
	runnables=(${tests[@]})
fi

if [ "$0" = "$BASH_SOURCE" ]; then
	for t in ${runnables[@]}; do
		if [ -d $t ]; then
			cd $t;
			bash do.sh $NP;
			cd - > /dev/null
		else
			echo "target $t not found"
		fi
	done
fi

echo "All done."

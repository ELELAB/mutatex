#!/bin/bash

unset FOLDX_ROTABASE
export FOLDX_BINARY=/usr/local/foldx5_2024/foldx
export FOLDX_VERSION=suite5
export NP=8
export tests=(basic basic_verbose basic_custom_rotabase basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_ptms basic_targz basic_skip_repair basic_multimodel_skip_repair interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers interaction_skip_repair selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers basic_multimers_poslist basic_poslist interaction_multimers_poslist selfmutate_interaction selfmutate_poslist selfmutate_skip_repair)

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

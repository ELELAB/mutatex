#!/bin/bash

export FOLDX_BINARY=/usr/local/foldx4_2020/foldx
export FOLDX_ROTABASE=/usr/local/foldx4_2020/rotabase.txt
export FOLDX_VERSION=suite4
export NP=4
export tests=(basic basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_ptms basic_targz interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers basic_poslist basic_multimers_poslist interaction_multimers_poslist)

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

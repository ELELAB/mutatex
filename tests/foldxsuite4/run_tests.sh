#!/bin/bash

export FOLDX_BINARY=/usr/local/foldx_2017/foldx
export FOLDX_ROTABASE=/usr/local/foldx_2017/rotabase.txt
export NP=16
export tests=(basic basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_targz interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers)

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

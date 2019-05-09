#!/bin/bash

unset FOLDX_ROTABASE
export FOLDX_BINARY=/usr/local/foldx5/foldx
export FOLDX_VERSION=suite5
export NP=4
export tests=(basic basic_custom_rotabase basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_targz interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers)

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

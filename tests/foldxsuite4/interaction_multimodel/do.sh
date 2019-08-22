cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex 2n9x_1_2.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 --binding-energy -a &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

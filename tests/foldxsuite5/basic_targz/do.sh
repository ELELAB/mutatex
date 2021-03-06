cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex 2klz_1.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 --compress &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

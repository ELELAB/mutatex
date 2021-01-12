cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex repaired_3wim0.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 --binding-energy --skip-repair &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

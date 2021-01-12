cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex 2klz_repaired_models.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 -a --skip-repair &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

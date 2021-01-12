cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex repaired_2klz_1.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 --self-mutate --skip-repair &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

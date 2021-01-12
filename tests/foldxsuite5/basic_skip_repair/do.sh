cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex repaired_2klz_1.pdb -m mutation_list.txt --foldx-version=$FOLDX_VERSION --np $NP --nruns=2 --skip-repair &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

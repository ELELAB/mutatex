cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
python mutatex.py 2klz_1.pdb --foldx-version=$FOLDX_VERSION -m mutation_list.txt --np $NP --nruns=2 --clean none &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

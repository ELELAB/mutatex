cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
mutatex 2klz_1.pdb -m mutation_list.txt --foldx-version=$FOLDX_VERSION --np $NP --nruns=2 --clean deep &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

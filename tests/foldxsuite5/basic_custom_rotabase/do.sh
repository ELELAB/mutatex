cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
python mutatex.py 2klz_1.pdb -m mutation_list.txt --foldx-version=$FOLDX_VERSION --np $NP --nruns=2 --rotabase ./rotabase.txt &> mutatex.log 

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

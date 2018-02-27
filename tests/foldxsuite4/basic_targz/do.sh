cwd=$(pwd)

echo -n "Now running $(basename $cwd)... "
python mutatex.py 2KLZ_1.pdb -m mutation_list.txt --np $NP --nruns=2 --compress &> mutatex.log

if [[ $? -eq 0 ]]; then
	echo "PASSED"
else
	echo "FAILED"
fi

export FOLDX_BINARY=/usr/local/foldx4/foldx4
export FOLDX_ROTABASE=/usr/local/foldx4/rotabase.txt

python mutatex.py 2KLZ_1_2.pdb -m mutation_list.txt --np 4 --foldx-version=4 -u

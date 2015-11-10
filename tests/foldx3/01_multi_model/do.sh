export FOLDX_BINARY=/usr/local/foldx3b6/foldx3b6
export FOLDX_ROTABASE=/usr/local/foldx3b6/rotabase.txt

python mutatex.py 2KLZ_1_2.pdb -m mutation_list.txt --np 4 -u

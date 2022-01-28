#chains A (spike) and D (ace2) of 7T9L PDB ID
#water and NAG removed 

mutatex 7T9Lad_noHOH.pdb --np 4 \
	-s \
	-q poslist.txt \
	-m mutation_list.txt \
	-x /usr/local/foldx5_2021/foldx \
	--foldx-version=suite5 \
	-B \
	-C deep \
	-v \
	-L \
	-c \
	-l


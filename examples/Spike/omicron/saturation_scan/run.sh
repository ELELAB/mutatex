###GENERAL_INFO###

#chains A (spike) and D (ace2) of 7t9l PDB ID
#water and NAG removed 
#cryo-EM structure


#saturation scan
mutatex 7T9Lad_noHOH.pdb \
	--np 8 \
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


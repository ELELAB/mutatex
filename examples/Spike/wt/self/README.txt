###GENERAL_INFO###

#chains A (ace2) and B (spike) of 6LZG PDB ID
#crystallographic water and NAG removed 


###REQUIRED_FILES###

# 6LZGab_noHOH.pdb
# mutation_list.txt
# repair_runfile_template.txt
# mutate_runfile_template.txt
# interface_runfile_template.txt


###COMMANDLINE###

source /usr/local/envs/mutatex/bin/activate

#self scan
tsp -N 4 mutatex 6LZGab_noHOH.pdb --np 4 -s -q poslist.txt -m mutation_list.txt -x /usr/local/foldx5_2021/foldx --foldx-version=suite5 -B -C deep -v -L -c -l &> log &

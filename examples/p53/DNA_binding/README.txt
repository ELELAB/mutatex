###GENERAL_INFO###

#chain B  of 3KZ8 (p53)
#chain C,D  DNA two strands
#crystallographic water removed 


###REQUIRED_FILES###

# 3KZ8b_91-289_DNA.pdb
# mutation_list.txt
# repair_runfile_template.txt
# mutate_runfile_template.txt
# interface_DNA_runfile_template.txt


###COMMANDLINE###

source /usr/local/envs/mutatex/bin/activate

#saturation
tsp -N 8 mutatex 3KZ8b_91-289_DNA.pdb --np 8 -m mutation_list.txt -x /usr/local/foldx5_2022/foldx --foldx-version=suite5 -B -C deep -v -L -c -l 

#plotting

ddg2heatmap -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -x 5 -c viridis -s 40 -f 8

ddg2distribution -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -T scatter -u 20 

ddg2logo -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -x 10 -t 0.8 | -T 0.8

ddg2density -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt

ddg2dg -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/

ddg2excel -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -F csv

pdb2labels -p 3KZ8b_91-289_DNA_model0_checked.pdb 

ddg2histo -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -r AB276

ddg2summary -p 3KZ8b_91-289_DNA_model0_checked.pdb -d results/interface_ddgs/final_averages/B-CD/ -l mutation_list.txt -L mutations_of_interest.txt

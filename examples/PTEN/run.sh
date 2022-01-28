mutatex 1D5R_noheatm.pdb \
	-p 4 \
	-m mutation_list.txt \
	-x /usr/local/foldx5/foldx \
	-f suite5 \
	-R repair_runfile_template.txt \
	-M mutate_runfile_template.txt \

# manual cleaning for this run
cd mutations
find . -name runfile.txt | xargs rm
find . -name individual_list.txt | xargs rm
find . -name rotabase.txt | xargs rm
find . -name '*.pdb' | xargs rm
cd ..
tar cvzf mutations.tar.gz mutations
rm -rf mutations
	

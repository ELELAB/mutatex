. /home/teo/envs/mutatex/bin/activate
mutatex 3tl0.pdb \
	-p 4 \
	-m mutation_list.txt \
	-x /usr/local/foldx4/foldx \
	-f suite4 \
	--rotabase /usr/local/foldx4/rotabase.txt \
	-R repair_runfile_template.txt \
	-M mutate_runfile_template.txt \
	-I interface_runfile_template.txt \
	-B

rm -rf mutations

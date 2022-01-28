. /home/teo/envs/mutatex/bin/activate
mutatex 1ju5.pdb \
	-p 4 \
	-m mutation_list.txt \
	-x /usr/local/foldx4/foldx \
	--rotabase /usr/local/foldx4/rotabase.txt \
	-f suite4 \
	-R repair_runfile_template.txt \
	-M mutate_runfile_template.txt \
	-I interface_runfile_template.txt \
	-B

rm -rf mutations

. /home/teo/envs/mutatex/bin/activate
mutatex 1jyr.pdb \
	-p 4 \
	-m mutation_list.txt \
	-x /usr/local/foldx5/foldx \
	-f suite5 \
	-R ../../../templates/foldxsuite4/repair_runfile_template.txt \
	-M ../../../templates/foldxsuite4/mutate_runfile_template.txt \
	-I ../../../templates/foldxsuite4/interface_runfile_template.txt \
	-B

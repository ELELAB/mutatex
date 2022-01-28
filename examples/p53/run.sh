# notice that this was run with a previous version of mutatex, 
# the --skip-check option is not available anymore. Just remove
# it to rerun the scan with a more modern version

mutatex 3kmd.pdb \
	-p 4 \
	-m mutation_list.txt \
	-x /usr/local/foldx4/foldx \
	-f suite4 \
	--rotabase /usr/local/foldx4/rotabase.txt \
	-R repair_runfile_template.txt \
	-M mutate_runfile_template.txt \
	-I interface_runfile_template.txt \
	--clean-deep \
	--binding-energy \
	--skip-check \
	--foldx-log


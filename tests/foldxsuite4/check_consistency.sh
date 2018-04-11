#!/bin/bash
#cat 2KLZ_1_2.pdb | grep CA | grep ATOM | tr -s ' ' | cut -d ' ' -f1-6 | sort | uniq | wc -l

declare -A expected_lengths=( ["_2klz"]="46" ["_2n9x"]="137" ["_3i3c"]="170" ["_3wim"]="129" )
declare -A expected_lengths_nomulti=( ["_3i3c"]="227" )

export tests=(basic basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers basic_targz interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers selfmutate selfmutate_multimers selfmutate_multimodel selfmutate_nomultimers)
export test_mutations=(basic basic_deepclean basic_multimers basic_multimodel basic_noclean basic_nomultimers interaction interaction_deepclean interaction_multimers interaction_multimodel interaction_noclean interaction_nomultimers)

function check_nres {
	ret="true"
	pdb=$(ls $1/*.pdb | grep -v checked)
        pdbid="_"$(basename $pdb | cut -c1-4) 
	if [[ $1 = *"nomultimers"* ]]; then
		explen=${expected_lengths_nomulti[$pdbid]}
	else
		explen=${expected_lengths[$pdbid]}
	fi

	for d in $(ls $1/mutations); do
		nfiles=$(ls $1/mutations/$d/ | wc -l)
		#echo $pdbid
		#echo $pdbid A"${expected_lengths[$pdbid]}"A A"$nfiles"A
		if [[ $nfiles -ne  $explen ]]; then
			ret="false"
			echo $ret
			return
		fi
	done
	echo $ret
}

function parse_foldx_mutation_average {
	tail -n +10 $1 | cut -f2,3 | xargs printf "%.5f "
}

function parse_mutatex_average {
	tail -n +2 $1 | awk '// {print $2" "$1}' | tr '\n' ' '
}


echo "Test output size"
for t in ${test_mutations[@]}; do
	nout=$(check_nres $t)
	if [ $nout = "true" ]; then
		echo $t PASSED
	else
		echo $t FAILED
	fi
done

#export test_mutations=(basic)
echo -e "\n\nTest DDG values"
for t in ${test_mutations[@]}; do
	passed="true"
    mut_dirs_str=$(ls $t/*.pdb | grep checked | xargs -I% basename -s '.pdb' % | awk '{print $1"_Repair"}' | tr '\n' ' ')
    mut_dirs=(${mut_dirs_str})
    for mdir in ${mut_dirs[@]}; do
    	all_mutfiles=$(ls $t"/mutations/"$mdir"/") 
    	for m in $all_mutfiles; do
    		mut_fname=$(basename $m)
    		foldx_fname="Average_"$mdir".fxout"
    		foldx_values_str=$(parse_foldx_mutation_average $t"/mutations/"$mdir"/"$m"/"$foldx_fname)
    		mutatex_values_str=$(parse_mutatex_average $t"/results/mutation_ddgs/"$mdir"/"$mut_fname)
    		all_values=$foldx_values_str" "$mutatex_values_str
    		check_outcome=$(echo $all_values | awk '
    		{
 	        	for (i=1;i<=NF/2;i++) {
                	#print $i" "$(i+NF/2)" "$i-$(i+NF/2)
                	if ($i-$(i+NF/2)>0.0001) {
                        print "false"
                        exit 0
                	}
        		}
        	print "true"
        	exit 0
        	}')
    		if [ $check_outcome = 'true' ]; then
    			continue
    		elif [ $check_outcome = 'false' ]; then
    			passed="false"
    			break
    		else
    			passed="wrong"
    			break
    		fi
        done
        if [ $passed = 'false' ]; then
        	break
        fi
    done
    if [ $passed = 'true' ]; then
    	echo $t PASSED
    elif [ $passed = 'false' ]; then
    	echo $t FAILED
    else
    	echo $t testing failed!
    fi
done
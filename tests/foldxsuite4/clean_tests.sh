source run_tests.sh;

delenda=()

if [ $# -gt 0 ]; then
	delenda=$@
else
	delenda=(${tests[@]})
fi

for arg in ${delenda[@]}; do
	if [ -d $arg ]; then
		echo "cleaning $arg"
		cd $arg;
		rm -rf repair mutations selfmutations results *checked*.pdb mutatex.log mutations.tar.gz;
		cd - >/dev/null;
	else
		echo "target $arg doesn't exist"
	fi
done

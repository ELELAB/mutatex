### INFO ###
# The script calculates the cutoff change in free energy as ±1.2 kcal/mol and returns the residues overcoming the threshold.

### COMMANDLINE ####
more selfmutation_energies.dat | awk '{if ($2 > 1.2 || $4 > 1.2 || $5 > 1.2 || $2 < -1.2 || $4 < -1.2 || $5 < -1.2) print $0}' 

#      avg	   std	     min	max

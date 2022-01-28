### INFO ###
# The script calculates the cutoff change in free energy as ±1.2 kcal/mol and returns the residues overcoming the threshold.

### COMMANDLINE ####
more selfmutation_energies.dat | awk '{if ($2 > 1.2 || $4 > 1.2 || $5 > 1.2 || $2 < -1.2 || $4 < -1.2 || $5 < -1.2) print $0}' 

#      avg	   std	     min	max
EA406  -0.432840   0.641091  -1.300240   0.307589
CA432  -1.351256   1.660617  -3.411040   0.009135

ddg2heatmap -p ../wt/saturation_scan/6LZGab_noHOH.pdb -d ../wt/saturation_scan/results/interface_ddgs/final_averages/A-B/ -l ../wt/saturation_scan/mutation_list.txt -q poslist_WT.txt -s 50 -F Arial -n -3 -x 5 -c plasma -b labels_6LZG_spike.csv
ddg2distribution -p ../wt/saturation_scan/6LZGab_noHOH.pdb -d ../wt/saturation_scan/results/interface_ddgs/6LZGab_noHOH_model0_checked_Repair/A-B/ -l ../wt/saturation_scan/mutation_list.txt -q poslist_WT.txt -T box -F Arial -s 50 -n -3 -x 16 -b labels_6LZG_spike.csv
ddg2summary -p ../wt/saturation_scan/6LZGab_noHOH.pdb -d ../wt/saturation_scan/results/interface_ddgs/final_averages/A-B/ -l ../wt/saturation_scan/mutation_list.txt -L summary_list_WT.txt
rename 's/^/wt_/' *.pdf summary.txt


ddg2heatmap -p ../omicron/saturation_scan/7T9Lad_noHOH.pdb -d ../omicron/saturation_scan/results/interface_ddgs/final_averages/A-D/ -l ../omicron/saturation_scan/mutation_list.txt -q poslist_OMICRON.txt -s 50 -F Arial -n -3 -x 5 -c plasma -b labels_7T9L_spike.csv
ddg2distribution -p ../omicron/saturation_scan/7T9Lad_noHOH.pdb -d ../omicron/saturation_scan/results/interface_ddgs/final_averages/A-D/ -l ../omicron/saturation_scan/mutation_list.txt -q poslist_OMICRON.txt -T box -F Arial -s 50 -n -3 -x 16 -b labels_7T9L_spike.csv
ddg2summary -p ../omicron/saturation_scan/7T9Lad_noHOH.pdb -d ../omicron/saturation_scan/results/interface_ddgs/final_averages/A-D/ -l ../omicron/saturation_scan/mutation_list.txt -L summary_list_OMICRON.txt
rename 's/^/omicron_/' *.pdf summary.txt


#!/bin/bash
#
#==============================================================
# Try to submit many jobs with different arguments
#==============================================================


for i in {1..20};
do
	for j in {1..10};
	do	
		for k in {1..10};
		do
			if grep -Fxq "${i}_${j}_${k}" missing_files.txt
			then
			  echo "ww_slmuq-"$i"-"$j"-"$k"-"
			  condor_submit htc_job_grid_ww_slmuq.sub arguments="${i} ${j} ${k}"
              sleep 0.1		
			else 
			  echo "Already have ww_slmuq-"$i"-"$j"-"$k"-"
			fi
		done
	done
done

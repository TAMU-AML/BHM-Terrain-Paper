source path.sh
num_parallel=6
num_jobs=0
turbine=1
n_turbines=66
while [[ $turbine -le $n_turbines ]]; do
  if [[ $num_jobs -lt $num_parallel ]]; then
   Rscript "$SOURCE_DIR/BHM_holdout_training.R" --test_turb "$turbine" > $LOG_DIR/bhm_turb_${turbine}_all_vars.log 2>&1 &
   num_jobs=$((num_jobs+1))
   turbine=$((turbine+1))
  else
    wait
    num_jobs=0
  fi
done
wait

Rscript $SOURCE_DIR/binning_prediction_all_turbines.R
echo "All jobs completed... Generating results now..."

Rscript $SOURCE_DIR/Table4.R
Rscript $SOURCE_DIR/Figure9.R
echo "Saved Figure 9 in the results folder with filename: Figure9a.pdf and Figure9b.pdf"
echo "Saved Table 4 in the results folder with filename: Table4.txt"

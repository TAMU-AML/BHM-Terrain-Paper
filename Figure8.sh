source path.sh
num_parallel=6
terrain_combinations=("1,2" "1,3" "2,3" "1" "2" "3")
echo "Submitting "$num_parallel" jobs in parallel"
num_jobs=0
for terrain_vars in "${terrain_combinations[@]}"; do
  if [[ $num_jobs -lt $num_parallel ]]; then
    echo "Submitting job for terrain vars: $terrain_vars"
    Rscript "$SOURCE_DIR/BHM-2dim.R" "None" "$terrain_vars" > $LOG_DIR/bhm_all_turb_"$terrain_vars"_vars.log 2>&1 &
    num_jobs=$((num_jobs+1))
    sleep 5
  else
    wait
    num_jobs=0
    echo "Submitting job for terrain vars: $terrain_vars"
    Rscript "$SOURCE_DIR/BHM-2dim.R" "None" "$terrain_vars" > $LOG_DIR/bhm_all_turb_"$terrain_vars"_vars.log 2>&1 &
    num_jobs=$((num_jobs+1))
    sleep 5
  fi
done
wait

echo "Generating all the subplots for Figure 8"
Rscript $SOURCE_DIR/Figure8.R
echo "Saved Figure 8 in the results folder with filename: Figure8a.pdf to Figure8f.pdf"
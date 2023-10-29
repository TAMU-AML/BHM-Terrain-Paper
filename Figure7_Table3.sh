source path.sh

echo "Running bhm model for all turbines using all terrain variables"
Rscript "$SOURCE_DIR/BHM-2dim.R" "None" "all" > $LOG_DIR/bhm_all_turb_all_vars.log 2>&1
if [ $? -ne 0 ]; then
 echo "Error while running the bhm model. Check the log file $LOG_DIR/bhm_all_turb_all_vars.log for details"
 exit 1
fi
echo "bhm model training completed successfully; generating results..."
Rscript $SOURCE_DIR/Figure7_Table3.R
echo "Saved Figure 7 in the results folder with filename: Figure7a.pdf and Figure7b.pdf"
echo "Saved Table 3 in the results folder with filename: Table3.txt"
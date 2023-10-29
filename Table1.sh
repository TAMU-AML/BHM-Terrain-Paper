source path.sh

echo "Running script for Table 1"
Rscript $SOURCE_DIR/Table1.R --method "knn" > $LOG_DIR/Table1.log 2>&1
echo "Saved output in the results folder with filename: Table1.txt"
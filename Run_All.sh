echo "=========================================================="
echo -e "Running all the code for the paper:\n \
\"A Bayesian hierarchical model to understand the \
effect of terrain on wind turbine power curves\""
echo ""
echo "Repo link: https://github.com/TAMU-AML/BHM-Terrain-Paper"
echo "=========================================================="

echo "Generating Table 1"
bash Table1.sh
echo ""
echo "=========================================================="

echo "Generating Figure 2"
bash Figure2.sh
echo ""
echo "=========================================================="

echo "Generating Figure 3"
bash Figure3.sh
echo ""
echo "=========================================================="

echo "Generating Figure 4"
bash Figure4.sh
echo ""
echo "=========================================================="

echo "Generating Figure 5"
bash Figure5.sh
echo ""
echo "=========================================================="

echo "Generating Figure 7 and Table3"
bash Figure7_Table3.sh
echo ""
echo "=========================================================="

echo "Generating Figure 8"
bash Figure8.sh
echo ""
echo "=========================================================="

echo "Generating Table 4 and Figure 9"
bash Table4_Figure9.sh
echo ""
echo "=========================================================="

echo "Finished running all the code"

echo "=========================================================="


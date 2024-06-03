srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_0.qasm &> sabre_log/54QBT_35CYC_QSE_0.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_1.qasm &> sabre_log/54QBT_35CYC_QSE_1.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_2.qasm &> sabre_log/54QBT_35CYC_QSE_2.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_3.qasm &> sabre_log/54QBT_35CYC_QSE_3.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_35CYC_QSE_4.qasm &> sabre_log/54QBT_35CYC_QSE_4.out


srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_0.qasm &> sabre_log/54QBT_45CYC_QSE_0.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_1.qasm &> sabre_log/54QBT_45CYC_QSE_1.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_2.qasm &> sabre_log/54QBT_45CYC_QSE_2.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_3.qasm &> sabre_log/54QBT_45CYC_QSE_3.out
srun python3 -u run_mlqls.py --sabre --dt sycamore --qf ../QUEKO-benchmark/BNTF/54QBT_45CYC_QSE_4.qasm &> sabre_log/54QBT_45CYC_QSE_4.out

srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/adder_n28/adder_n28_transpiled.qasm &> sabre_log/36_adder_n28.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/knn_n31/knn_n31_transpiled.qasm &> sabre_log/36_knn_31.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 9 --qf ../QASMBench/large/knn_n67/knn_n67_transpiled.qasm &> sabre_log/81_knn_67.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 19 --qf ../QASMBench/large/knn_n341/knn_341_transpiled.qasm &> sabre_log/361_knn_341.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 10 --qf ../QASMBench/large/adder_n64/adder_n64_transpiled.qasm &> sabre_log/100_adder_n64.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/bv_n30/bv_n30_transpiled.qasm &> sabre_log/36_bv_30.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/ising_n34/ising_n34.qasm &> sabre_log/36_isling_34.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 10 --qf ../QASMBench/large/ising_n98/ising_n98.qasm &> sabre_log/100_isling_34.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/qft_n29/qft_n29_transpiled.qasm &> sabre_log/36_qft_n29.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 8 --qf ../QASMBench/large/qft_n63/qft_n63_transpiled.qasm &> sabre_log/64_qft_n63.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 7 --qf ../QASMBench/large/ghz_n40/ghz_n40_transpiled.qasm &> sabre_log/7_ghz_40.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/QV_n32/32.qasm &> sabre_log/36_QV_n32.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 7 --qf ../QASMBench/large/swap_test_n41/swap_test_n41_transpiled.qasm &> sabre_log/swap_test_n41.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 10 --qf ../QASMBench/large/swap_test_n83/swap_test_n83_transpiled.qasm &> sabre_log/swap_test_n83.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 11 --qf ../QASMBench/large/adder_n118/adder_n118_transpiled.qasm &> sabre_log/121_adder_n118.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 9 --qf ../QASMBench/large/bv_n70/bv_n70_transpiled.qasm &> sabre_log/81_bv_70.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 12 --qf ../QASMBench/large/bv_n140/bv_n140_transpiled.qasm &> sabre_log/144_bv_140.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 17 --qf ../QASMBench/large/bv_n280/bv_n280_transpiled.qasm &> sabre_log/289_bv_280.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 9 --qf ../QASMBench/large/cat_n65/cat_n65_transpiled.qasm &> sabre_log/81_cat_65.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 12 --qf ../QASMBench/large/cat_n130/cat_n130_transpiled.qasm &> sabre_log/144_cat_130.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/dnn_n33/dnn_n33_transpiled.qasm &> sabre_log/36_dnn_33.out
srun python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../QASMBench/large/wstate_n36/wstate_n36_transpiled.qasm &> sabre_log/36_wstate_33.out

# srun python3 -u run_mlqls.py --sabre --dt grid --d 7 --qf ../QASMBench/large/square_root_n45/square_root_n45_transpiled.qasm &> sabre_log/48_square_root_45.out
# srun python3 -u run_mlqls.py --sabre --dt grid --d 8 --qf ../QASMBench/large/square_root_n60/square_root_n60_transpiled.qasm &> sabre_log/48_square_root_60.out

srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/adder_n28/adder_n28_transpiled.qasm &> sabre_log/eagle_adder_n28.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/knn_n31/knn_n31_transpiled.qasm &> sabre_log/eagle_knn_31.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/knn_n67/knn_n67_transpiled.qasm &> sabre_log/eagle_knn_67.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/adder_n64/adder_n64_transpiled.qasm &> sabre_log/eagle_adder_n64.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/bv_n30/bv_n30_transpiled.qasm &> sabre_log/eagle_bv_30.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/ising_n34/ising_n34.qasm &> sabre_log/eagle_isling_34.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/ising_n98/ising_n98.qasm &> sabre_log/eagle_isling_98.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/qft_n29/qft_n29_transpiled.qasm &> sabre_log/eagle_qft_n29.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/qft_n63/qft_n63_transpiled.qasm &> sabre_log/eagle_qft_n63.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/ghz_n40/ghz_n40_transpiled.qasm &> sabre_log/eagle_ghz_40.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/QV_n32/32.qasm &> sabre_log/eagle_QV_n32.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/swap_test_n41/swap_test_n41_transpiled.qasm &> sabre_log/eagle_swap_test_n41.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/swap_test_n83/swap_test_n83_transpiled.qasm &> sabre_log/eagle_swap_test_n83.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/adder_n118/adder_n118_transpiled.qasm &> sabre_log/eagle_adder_n118.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/bv_n70/bv_n70_transpiled.qasm &> sabre_log/eagle_bv_70.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/cat_n65/cat_n65_transpiled.qasm &> sabre_log/eagle_cat_65.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/dnn_n33/dnn_n33_transpiled.qasm &> sabre_log/eagle_dnn_33.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/wstate_n36/wstate_n36_transpiled.qasm &> sabre_log/eagle_wstate_33.out

# srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/square_root_n45/square_root_n45_transpiled.qasm &> sabre_log/eagle_square_root_45.out
# srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../QASMBench/large/square_root_n60/square_root_n60_transpiled.qasm &> sabre_log/eagle_square_root_60.out

srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --all_commute &> sabre_log/eagle_qaoa_24_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_36_0.qasm --all_commute &> sabre_log/eagle_qaoa_36_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm --all_commute &> sabre_log/eagle_qaoa_40_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm --all_commute &> sabre_log/eagle_qaoa_50_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm --all_commute &> sabre_log/eagle_qaoa_60_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm --all_commute &> sabre_log/eagle_qaoa_70_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm --all_commute &> sabre_log/eagle_qaoa_80_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_90_0.qasm --all_commute &> sabre_log/eagle_qaoa_90_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_100_0.qasm --all_commute &> sabre_log/eagle_qaoa_100_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_110_0.qasm --all_commute &> sabre_log/eagle_qaoa_110_0.out
srun python3 -u run_mlqls.py --sabre --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_120_0.qasm --all_commute &> sabre_log/eagle_qaoa_120_0.out

srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 5 --qf ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --all_commute &> sabre_log/qaoa_24_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 6 --qf ../quantum_cir_benchmark/qaoa/qaoa_36_0.qasm --all_commute &> sabre_log/qaoa_36_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 7 --qf ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm --all_commute &> sabre_log/qaoa_40_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 8 --qf ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm --all_commute &> sabre_log/qaoa_50_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 8 --qf ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm --all_commute &> sabre_log/qaoa_60_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 9 --qf ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm --all_commute &> sabre_log/qaoa_70_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 9 --qf ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm --all_commute &> sabre_log/qaoa_80_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 10 --qf ../quantum_cir_benchmark/qaoa/qaoa_90_0.qasm --all_commute &> sabre_log/qaoa_90_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 10 --qf ../quantum_cir_benchmark/qaoa/qaoa_100_0.qasm --all_commute &> sabre_log/qaoa_100_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 11 --qf ../quantum_cir_benchmark/qaoa/qaoa_110_0.qasm --all_commute &> sabre_log/qaoa_110_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 11 --qf ../quantum_cir_benchmark/qaoa/qaoa_120_0.qasm --all_commute &> sabre_log/qaoa_120_0.out

srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 15 --qf ../quantum_cir_benchmark/qaoa/qaoa_200_0.qasm --all_commute &> sabre_log/qaoa_200_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 16 --qf ../quantum_cir_benchmark/qaoa/qaoa_250_0.qasm --all_commute &> sabre_log/qaoa_250_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 18 --qf ../quantum_cir_benchmark/qaoa/qaoa_300_0.qasm --all_commute &> sabre_log/qaoa_300_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 19 --qf ../quantum_cir_benchmark/qaoa/qaoa_350_0.qasm --all_commute &> sabre_log/qaoa_350_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 20 --qf ../quantum_cir_benchmark/qaoa/qaoa_400_0.qasm --all_commute &> sabre_log/qaoa_400_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 22 --qf ../quantum_cir_benchmark/qaoa/qaoa_450_0.qasm --all_commute &> sabre_log/qaoa_450_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 23 --qf ../quantum_cir_benchmark/qaoa/qaoa_500_0.qasm --all_commute &> sabre_log/qaoa_500_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 24 --qf ../quantum_cir_benchmark/qaoa/qaoa_550_0.qasm --all_commute &> sabre_log/qaoa_550_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 25 --qf ../quantum_cir_benchmark/qaoa/qaoa_600_0.qasm --all_commute &> sabre_log/qaoa_600_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 26 --qf ../quantum_cir_benchmark/qaoa/qaoa_650_0.qasm --all_commute &> sabre_log/qaoa_650_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 27 --qf ../quantum_cir_benchmark/qaoa/qaoa_700_0.qasm --all_commute &> sabre_log/qaoa_700_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 28 --qf ../quantum_cir_benchmark/qaoa/qaoa_750_0.qasm --all_commute &> sabre_log/qaoa_750_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 29 --qf ../quantum_cir_benchmark/qaoa/qaoa_800_0.qasm --all_commute &> sabre_log/qaoa_800_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 30 --qf ../quantum_cir_benchmark/qaoa/qaoa_850_0.qasm --all_commute &> sabre_log/qaoa_850_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 30 --qf ../quantum_cir_benchmark/qaoa/qaoa_900_0.qasm --all_commute &> sabre_log/qaoa_900_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 31 --qf ../quantum_cir_benchmark/qaoa/qaoa_950_0.qasm --all_commute &> sabre_log/qaoa_950_0.out
srun --mem=2G python3 -u run_mlqls.py --sabre --dt grid --d 32 --qf ../quantum_cir_benchmark/qaoa/qaoa_1000_0.qasm --all_commute &> sabre_log/qaoa_1000_0.out


# srun python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/adder_n28/adder_n28_transpiled.qasm &> molsq_log/36_adder_n28.out
# srun python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/knn_n31/knn_n31_transpiled.qasm &> molsq_log/36_knn_31.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 9 --qf ../QASMBench/large/knn_n67/knn_n67_transpiled.qasm &> molsq_log/81_knn_67.out

# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 10 --qf ../QASMBench/large/adder_n64/adder_n64_transpiled.qasm &> molsq_log/100_adder_n64.out
# srun python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/bv_n30/bv_n30_transpiled.qasm &> molsq_log/36_bv_30.out
# srun --mem=10G python3 -u run_mlqls.py --dt grid --d 9 --qf ../QASMBench/large/bv_n70/bv_n70_transpiled.qasm &> molsq_log/81_bv_70.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 7 --qf ../QASMBench/large/swap_test_n41/swap_test_n41_transpiled.qasm &> molsq_log/swap_test_n41.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 10 --qf ../QASMBench/large/swap_test_n83/swap_test_n83_transpiled.qasm &> molsq_log/swap_test_n83.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/qft_n29/qft_n29_transpiled.qasm &> molsq_log/36_qft_n29.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 11 --qf ../QASMBench/large/adder_n118/adder_n118_transpiled.qasm &> molsq_log/121_adder_n118.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 12 --qf ../QASMBench/large/bv_n140/bv_n140_transpiled.qasm &> molsq_log/144_bv_140.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 17 --qf ../QASMBench/large/bv_n280/bv_n280_transpiled.qasm &> molsq_log/289_bv_280.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/dnn_n33/dnn_n33_transpiled.qasm &> molsq_log/36_dnn_33.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 8 --qf ../QASMBench/large/qft_n63/qft_n63_transpiled.qasm &> molsq_log/eagle_qft_n63.out

srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/knn_n31/knn_n31_transpiled.qasm &> molsq_log/eagle_knn_31.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/knn_n67/knn_n67_transpiled.qasm &> molsq_log/eagle_knn_67.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/bv_n30/bv_n30_transpiled.qasm &> molsq_log/eagle_bv_30.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/swap_test_n41/swap_test_n41_transpiled.qasm &> molsq_log/eagle_swap_test_n41.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/swap_test_n83/swap_test_n83_transpiled.qasm &> molsq_log/eagle_swap_test_n83.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/qft_n29/qft_n29_transpiled.qasm &> molsq_log/eagle_qft_n29.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/adder_n28/adder_n28_transpiled.qasm &> molsq_log/eagle_adder_n28.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/adder_n64/adder_n64_transpiled.qasm &> molsq_log/eagle_adder_n64.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/bv_n70/bv_n70_transpiled.qasm &> molsq_log/eagle_bv_70.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/dnn_n33/dnn_n33_transpiled.qasm &> molsq_log/eagle_dnn_33.out
srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/adder_n118/adder_n118_transpiled.qasm &> molsq_log/eagle_adder_n118.out

# srun --mem=30G python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/qft_n63/qft_n63_transpiled.qasm &> molsq_log/eagle_qft_n63.out

# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/QV_n32/32.qasm &> molsq_log/36_QV_n32.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 10 --qf ../QASMBench/large/qft_n63/qft_n63_transpiled.qasm &> molsq_log/100_qft_n63.out
# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/QV_n32/32.qasm &> molsq_log/eagle_QV_n32.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 19 --qf ../QASMBench/large/knn_n341/knn_341_transpiled.qasm &> molsq_log/361_knn_341.out


# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/ghz_n40/ghz_n40_transpiled.qasm &> molsq_log/eagle_ghz_40.out
# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/wstate_n36/wstate_n36_transpiled.qasm &> molsq_log/eagle_wstate_33.out
# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/cat_n65/cat_n65_transpiled.qasm &> molsq_log/eagle_cat_65.out
# srun python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/ising_n34/ising_n34.qasm &> molsq_log/36_isling_34.out
# srun python3 -u run_mlqls.py --dt grid --d 10 --qf ../QASMBench/large/ising_n98/ising_n98.qasm &> molsq_log/100_isling_98.out
# srun python3 -u run_mlqls.py --dt grid --d 7 --qf ../QASMBench/large/ghz_n40/ghz_n40_transpiled.qasm &> molsq_log/7_ghz_40.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 6 --qf ../QASMBench/large/wstate_n36/wstate_n36_transpiled.qasm &> molsq_log/36_wstate_33.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 9 --qf ../QASMBench/large/cat_n65/cat_n65_transpiled.qasm &> molsq_log/81_cat_65.out
# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/ising_n34/ising_n34.qasm &> molsq_log/eagle_isling_34.out
# srun python3 -u run_mlqls.py --dt eagle --qf ../QASMBench/large/ising_n98/ising_n98.qasm &> molsq_log/eagle_isling_98.out
# srun --mem=30G python3 -u run_mlqls.py --dt grid --d 12 --qf ../QASMBench/large/cat_n130/cat_n130_transpiled.qasm &> molsq_log/144_cat_130.out


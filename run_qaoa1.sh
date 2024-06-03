srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm --all_commute &> molsq_log/eagle_qaoa_80_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_90_0.qasm --all_commute &> molsq_log/eagle_qaoa_90_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_100_0.qasm --all_commute &> molsq_log/eagle_qaoa_100_0.out
srun --mem=10G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_110_0.qasm --all_commute &> molsq_log/eagle_qaoa_110_0.out
srun --mem=10G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_120_0.qasm --all_commute &> molsq_log/eagle_qaoa_120_0.out

srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --all_commute &> molsq_log/eagle_qaoa_24_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_36_0.qasm --all_commute &> molsq_log/eagle_qaoa_36_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm --all_commute &> molsq_log/eagle_qaoa_40_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm --all_commute &> molsq_log/eagle_qaoa_50_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm --all_commute &> molsq_log/eagle_qaoa_60_0.out
srun --mem=5G python3 -u run_mlqls.py --dt eagle --qf ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm --all_commute &> molsq_log/eagle_qaoa_70_0.out
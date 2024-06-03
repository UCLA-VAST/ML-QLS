srun --mem=2G python3 -u run_mlqls.py --dt grid --d 5 --qf ../quantum_cir_benchmark/qaoa/qaoa_24_0.qasm --all_commute &> molsq_log/qaoa_24_0.out
srun --mem=2G python3 -u run_mlqls.py --dt grid --d 6 --qf ../quantum_cir_benchmark/qaoa/qaoa_36_0.qasm --all_commute &> molsq_log/qaoa_36_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 7 --qf ../quantum_cir_benchmark/qaoa/qaoa_40_0.qasm --all_commute &> molsq_log/qaoa_40_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 8 --qf ../quantum_cir_benchmark/qaoa/qaoa_50_0.qasm --all_commute &> molsq_log/qaoa_50_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 8 --qf ../quantum_cir_benchmark/qaoa/qaoa_60_0.qasm --all_commute &> molsq_log/qaoa_60_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 9 --qf ../quantum_cir_benchmark/qaoa/qaoa_70_0.qasm --all_commute &> molsq_log/qaoa_70_0.out

srun --mem=5G python3 -u run_mlqls.py --dt grid --d 9 --qf ../quantum_cir_benchmark/qaoa/qaoa_80_0.qasm --all_commute &> molsq_log/qaoa_80_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 10 --qf ../quantum_cir_benchmark/qaoa/qaoa_90_0.qasm --all_commute &> molsq_log/qaoa_90_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 10 --qf ../quantum_cir_benchmark/qaoa/qaoa_100_0.qasm --all_commute &> molsq_log/qaoa_100_0.out
srun --mem=5G python3 -u run_mlqls.py --dt grid --d 11 --qf ../quantum_cir_benchmark/qaoa/qaoa_110_0.qasm --all_commute &> molsq_log/qaoa_110_0.out
srun --mem=10G python3 -u run_mlqls.py --dt grid --d 11 --qf ../quantum_cir_benchmark/qaoa/qaoa_120_0.qasm --all_commute &> molsq_log/qaoa_120_0.out




srun --mem=20G python3 -u run_mlqls.py --dt grid --d 15 --qf ../quantum_cir_benchmark/qaoa/qaoa_200_0.qasm --all_commute &> molsq_log/qaoa_200_0.out
srun --mem=20G python3 -u run_mlqls.py --dt grid --d 16 --qf ../quantum_cir_benchmark/qaoa/qaoa_250_0.qasm --all_commute &> molsq_log/qaoa_250_0.out
srun --mem=20G python3 -u run_mlqls.py --dt grid --d 18 --qf ../quantum_cir_benchmark/qaoa/qaoa_300_0.qasm --all_commute &> molsq_log/qaoa_300_0.out
srun --mem=40G python3 -u run_mlqls.py --dt grid --d 19 --qf ../quantum_cir_benchmark/qaoa/qaoa_350_0.qasm --all_commute &> molsq_log/qaoa_350_0.out
srun --mem=40G python3 -u run_mlqls.py --dt grid --d 20 --qf ../quantum_cir_benchmark/qaoa/qaoa_400_0.qasm --all_commute &> molsq_log/qaoa_400_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 22 --qf ../quantum_cir_benchmark/qaoa/qaoa_450_0.qasm --all_commute &> molsq_log/qaoa_450_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 23 --qf ../quantum_cir_benchmark/qaoa/qaoa_500_0.qasm --all_commute &> molsq_log/qaoa_500_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 24 --qf ../quantum_cir_benchmark/qaoa/qaoa_550_0.qasm --all_commute &> molsq_log/qaoa_550_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 25 --qf ../quantum_cir_benchmark/qaoa/qaoa_600_0.qasm --all_commute &> molsq_log/qaoa_600_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 26 --qf ../quantum_cir_benchmark/qaoa/qaoa_650_0.qasm --all_commute &> molsq_log/qaoa_650_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 27 --qf ../quantum_cir_benchmark/qaoa/qaoa_700_0.qasm --all_commute &> molsq_log/qaoa_700_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 28 --qf ../quantum_cir_benchmark/qaoa/qaoa_750_0.qasm --all_commute &> molsq_log/qaoa_750_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 29 --qf ../quantum_cir_benchmark/qaoa/qaoa_800_0.qasm --all_commute &> molsq_log/qaoa_800_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 30 --qf ../quantum_cir_benchmark/qaoa/qaoa_850_0.qasm --all_commute &> molsq_log/qaoa_850_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 30 --qf ../quantum_cir_benchmark/qaoa/qaoa_900_0.qasm --all_commute &> molsq_log/qaoa_900_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 31 --qf ../quantum_cir_benchmark/qaoa/qaoa_950_0.qasm --all_commute &> molsq_log/qaoa_950_0.out
# srun --mem=40G python3 -u run_mlqls.py --dt grid --d 32 --qf ../quantum_cir_benchmark/qaoa/qaoa_1000_0.qasm --all_commute &> molsq_log/qaoa_1000_0.out


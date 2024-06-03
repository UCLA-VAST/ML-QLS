#!/usr/bin/env python3
import argparse
import timeit
import os
import sys
sys.path.append(os.getcwd())
sys.path.append(os.getcwd()+'/lib')
from src.pyolsq.apiPy import createCircuit, createDevice
from src.pyolsq.input import input_qasm
from src.pyolsq.device import get_device_by_name, get_nnGrid, get_heavyHex
from olsqPy import Device, mOLSQ, Circuit, OLSQ2
from run_heuristic import run_sabre, run_qmap

def run_olsq_tbolsq(circuit: Circuit, device: Device, connection, all_commute):
    # start = timeit.default_timer()

    # device.printDevice()
    lsqc_solver = mOLSQ(circuit, device)
    if all_commute:
        lsqc_solver.enableAllCommute() # for QAOA circuits
    lsqc_solver.run()

    # stop = timeit.default_timer()
    # print('Time: ', stop - start)  

if __name__ == "__main__":
    # Initialize parser
    parser = argparse.ArgumentParser()
    # Adding optional argument
    parser.add_argument("--dt", dest='device_type', type=str,
        help="grid or heavy-hexagon or ring")
    parser.add_argument("--d", dest='device', type=int,
        help="device (x-by-x grid)")
    parser.add_argument("--qf", dest="qasm", type=str,
        help="Input file name")
    parser.add_argument("--all_commute", action='store_true', default=False,
        help="gates are all commmutable")
    parser.add_argument("--sabre", action='store_true', default=False,
        help="use sabre")
    
    # Read arguments from command line
    args = parser.parse_args()


    circuit_info = open(args.qasm, "r").read()
    circuit_name = args.qasm.split("/")[-1]
    circuit_name = "result/"+circuit_name.split(".")[0]
    circuit = createCircuit(circuit_name, circuit_info)
    # circuit.printCircuit()
    if args.device_type == "grid":
        connection, device = get_nnGrid(args.device)
    else:
        connection, device = get_device_by_name(args.device_type)
    
    qubit_num, gate_list, _ = input_qasm(circuit_info)
    if args.sabre:
        run_sabre(gate_list, connection, device.nQubit())
    else:
        run_olsq_tbolsq(circuit, device, connection, args.all_commute)
    # run_qmap(circuit_info, connection, device.nQubit())

    # circuit.printCircuitLayout()
    

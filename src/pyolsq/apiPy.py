# import os
# import sys
# sys.path.append(os.getcwd())
# sys.path.append(os.getcwd()+'/lib')
# print(sys.path)
from src.pyolsq.input import input_qasm
# from src.pyolsq.run_h_compiler import run_sabre, test_sabre
from olsqPy import Device, Circuit, mOLSQ, addGate, setEdge, setDependency, setInitialMapping

def createCircuit(name, program, is_qasm= True, gate_duration: dict = None):
    """Translate input program to Circuit
    Args:
        name: circuit name
        program: a qasm string, or a list of the three things in IR.
        input_mode: (optional) can be "IR" if the input has ben
            translated to mOLSQ IR; can be "benchmark" to use one of
            the benchmarks.  Default mode assumes qasm input.
    """
    
    if not is_qasm:
        count_program_qubit = program[0]
        list_gate_qubits = program[1]
        list_gate_name = program[2]
    else:
        program = input_qasm(program)
        count_program_qubit = program[0]
        list_gate_qubits = program[1]
        list_gate_name = program[2]

    circuit = Circuit(name, program[0], len(program[1]))
    if gate_duration != None:
        for g, n in zip(list_gate_qubits, list_gate_name):
            addGate(circuit, n, g, gate_duration[n])
    else:
        for g, n in zip(list_gate_qubits, list_gate_name):
            addGate(circuit, n, g)
    return circuit

def createDevice(name: str, nqubits: int = None, connection: list = None, xy_coordinate: list = None):
    """Pass in parameters from the given device.  If in TB mode,
        swap_duration is set to 1 without modifying the device.

    Args:
        device: a qcdevice object for mOLSQ
    """
    device = Device(name, nqubits, len(connection))
    setEdge(device, connection)
    if xy_coordinate != None:
        setXYCoord(device, xy_coordinate)
    return device


# def useSabre(nqubits, connection: list, program, is_qasm):
#     swap_num, depth = test_sabre(is_qasm, program, connection, nqubits)
#     print("[Info] Run heuristic compiler SABRE to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
#     return swap_num

# def useSabre(for_swap: bool, for_mapping: bool, olsq: mOLSQ, circuit: Circuit, nqubits, connection: list, program, is_qasm):
#     swap_num, depth, initial_mapping, _ = run_sabre(is_qasm, program, connection, nqubits)
#     print("[Info] Run heuristic compiler SABRE to get upper bound for SWAP: {}, depth: {}".format(swap_num, depth))
#     if for_swap:
#         olsq.setSabreForSwap(True, swap_num)
#     if for_mapping:
#         olsq.useCircuitInitalMapping()
#         setInitialMapping(circuit, initial_mapping)

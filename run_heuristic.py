from qiskit.transpiler import CouplingMap, Layout
from qiskit import QuantumCircuit
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import SabreLayout, SabrePreLayout
from qiskit.converters import *
import timeit

def run_sabre(circuit_info, coupling, count_physical_qubit):
    # read qasm
    list_gate = circuit_info
    qc = QuantumCircuit(count_physical_qubit)
    two_qubit_gate = 0
    single_qubit_gate = 0
    for gate in list_gate:
        if len(gate) == 2:
            qc.cx(gate[0], gate[1])
            two_qubit_gate += 1
        elif len(gate) == 1:
            qc.h(gate[0])
            single_qubit_gate += 1
        else:
            raise TypeError("Currently only support one and two-qubit gate.")
    print("gate num: {}".format(len(list_gate)))
    print("2Q gate num: {}".format(two_qubit_gate))
    print("1Q gate num: {}".format(single_qubit_gate))
    
    # qc.draw(scale=0.7, filename = "cir_for_tket.png", output='mpl', style='color')
    device = CouplingMap(couplinglist = coupling, description="sabre_test")
    
    # initialize sabre
    # ["basic", "lookahead", "decay"]
    pm = PassManager(
    [
        # SabrePreLayout(coupling_map=device),
        SabreLayout(device, seed = 0, layout_trials = 4),
    ]
    )

    start = timeit.default_timer()
    
    sabre_cir = pm.run(qc)

    stop = timeit.default_timer()
    print('Time: ', stop - start)  

    sabre_layout = sabre_cir.layout
    count_swap = 0
    other_gate = 0
    qubit_push_back_depth = [ 0 for i in range(count_physical_qubit)]
    swap_time = []
    circuit_depth = 0
    for gate in sabre_cir.data:
        gate_time = 0
        for q in gate.qubits:
            gate_time = max(qubit_push_back_depth[q._index], gate_time)
        if gate[0].name == 'swap':
            count_swap += 1
            swap_time.append(gate_time)
        else:
            other_gate += 1
        gate_time += 1
        circuit_depth = max(circuit_depth, gate_time)
        for q in gate.qubits:
            qubit_push_back_depth[q._index] = gate_time 
    # print(sabre_layout.initial_layout)
    # for i in sabre_layout.initial_layout:
    #     l.append(i)
    # print(l)
    # import matplotlib.pyplot as plt
    # plt.hist(swap_time, bins=range(min(swap_time), max(swap_time) + 2, 2))
    # print(swap_time)
    # plt.hist(swap_time)
    # plt.savefig('swap_time.png')
    # results = {"count_swap: ": count_swap,
    #             "depth:": sabre_cir.depth()}
    print("count_swap: ", count_swap)
    print("depth:", sabre_cir.depth())
    
    assert(len(list_gate) == other_gate)
    # return count_swap, sabre_cir.depth()



def run_qmap(circuit_info, coupling, count_physical_qubit):
    from mqt import qmap
    # read qasm
    # list_gate = circuit_info
    # qc = QuantumCircuit(count_physical_qubit)
    # print(qubit_num)
    # print("gate num: {}".format(len(list_gate)))
    # print(list_gate)
    # for gate in list_gate:
    #     if len(gate) == 2:
    #         qc.cx(gate[0], gate[1])
    #     elif len(gate) == 1:
    #         qc.h(gate[0])
    #     else:
    #         raise TypeError("Currently only support one and two-qubit gate.")
    qc = QuantumCircuit.from_qasm_str(circuit_info)
    
    # qc.draw(scale=0.7, filename = "cir_for_tket.png", output='mpl', style='color')
    device = qmap.Architecture(
        count_physical_qubit, {
            (t[0], t[1]) for t in coupling
        }
    )

    qc_mapped, res = qmap.compile(qc, device, method="heuristic", post_mapping_optimizations=False)
    print("Additional SWAPs: %d" % res.output.swaps)
    print("Detph: %d" % qc_mapped.depth())
    print("Runtime:          %f" % res.time)
    return 
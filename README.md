# mOLSQ: Multilevel Quantum Layout Synthesis

## Build and Installation Instructions

```
# Clone mOLSQ repository
git clone git@github.com:WanHsuanLin/mOLSQ.git
cd mOLSQ

# Build
cmake . -Bbuild
cd build
make
```

- Required Dependencies: 
  - GMP: Please install GMP. You may need to compile it with `-fPIC` depending on your system.
  - Bitwuzla (https://github.com/bitwuzla/bitwuzla/tree/4eda0536800576cb2531ab9ce13292da8f21f0eb): Please install Bitwuzla.
  - pblib (https://github.com/master-keying/pblib): Please intall pblib in the `include/`. You may need to compile it with `-fPIC` depending on your system.

## Create a device from the Input Coupling Graph

To perform QLS, we need to know the connections between the qubits, which is information about the physical device.
We are going to use the `createDevice` function.

```
from olsqPy import Device
from src.pyolsq.apiPy import createDevice
device = createDevice(name="dev", nqubits=5, connection=[(0, 1), (1, 2), (1, 3), (3, 4)])
device.printDevice()
```

We use a minimalist class `qcdevice` to store the properties of the device that we need, which can be constructed with these arguments.
(The last three are only for fidelity optimization.)
- `name`
- `nqubits`: the number of physical qubits
- `connection`: a list of physical qubit pairs corresponding to edges in the coupling graph
- `swap_duration`: number of cycles a SWAP gate takes.
   Usually it is either one, or three meaning three CX gates.
- `fmeas`: a list of measurement fidelity
- `fsingle`: a list of single-qubit gate fidelity
- `ftwo`: a list of two-qubit gate fidelity, indices aligned with `connection`

If `name` starts with `"default_"`, a hard-coded device stored in `olsq/devices/` would be loaded.
Other arguments can still be specified, in which case the original device properties would be replaced by the input.
```
# use a hard-coded device in olsq/devices/ called ourense
# which actually has the same properties as the device we constructed above
lsqc_solver.setdevice( qcdevice("default_ourense") )
```

## Create a Circuit from the Input Program

Apart from the device, we need the quantum program/circuit to execute, which can be constructed with the `createCircuit` function.

mOLSQ has an intermediate representation (IR) of quantum programs. (For details, refer to [a later part](#olsq-ir) of this tutorial.)
In general, there are four ways to set the program: 
1. Use OLSQ IR
2. Use a string in QASM format
```
from olsqPy import Circuit
from src.pyolsq.apiPy import createCircuit
circuit_name = "toffoli"
circuit_str = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nh q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[1];\nt q[2];\n" \
              "cx q[0], q[1];\nh q[2];\nt q[0];\ntdg q[1];\ncx q[0], q[1];\n"

# input the quantum program as a QASM string
lsqc_solver.setprogram(circuit_str)
circuit = createCircuit(circuit_name, circuit_str, is_qasm = true)
circuit.printCircuit()
```


## Initialization and solve

```
from olsqPy import Device, Circuit, mOLSQ

# initialize your circuit
# initialize your device

lsqc_solver = mOLSQ(circuit, device)
lsqc_solver.enableAllCommute() # for circuits whose gates are all commutable
# lsqc_solver.disAllCommute() # for circuits whose gates are not commutable (default)

lsqc_solver.run()
```

There are two argument in the constructor of mOLSQ: `Circuit` and `Device`.
The former one stands for the input circuit, and the later one stands for the input device.
Idealy, the `run` method will return three solutions, `${circuit_name}_stage_0.qasm`, `${circuit_name}_stage_1.qasm`, and `${circuit_name}_stage_2.qasm`, from FastQLS (hueristic QLS algorithm), the first multilevel V cycle, and the second multilevel V cycle, respectively.

## OLSQ IR

OLSQ IR contains three things:
1. `count_program_qubit`: the number of qubits in the program.
2. `gates`: a list of tuples representing qubit(s) acted on by a gate, each tuple has one index if it is a single-qubit gate, two indices if it is a two-qubit gate.
3. `gate_spec`: list of type/name of each gate, which is not important to OLSQ, and only needed when generating output.

```
# For the following circuit
# q_0: ───────────────────■───
#                         │  
# q_1: ───────■───────────┼───
#      ┌───┐┌─┴─┐┌─────┐┌─┴─┐
# q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├─
#      └───┘└───┘└─────┘└───┘ 

# count_program_qubit = 3
# gates = ((2,), (1,2), (2,), (0,1))
# gate_spec = ("h", "cx", "tdg", "cx")
```

## Example: run_mlqls.py

run_mlqls.py is an example program to use mOLSQ2 to perform layout synthesis. The output will be stored in under `result/`
```
# compile an qaoa circuit on a 5-by-5 grid quantum device
python3 run_mlqls.py --dt grid --d 4 -qf benchmark/qaoa/qaoa_16_0.qasm
# The output files (Final IR output file and the intermediate qasm file) of running the command are in example/.

# compile an qaoa circuit on sycamore quantum device
python3 run_mlqls.py --dt sycamore --qf benchmark/qaoa/qaoa_16_0.qasm
```
- `--dt $(str)`: Type of the quantum device: ourense, sycamore, rochester, tokyo, aspen-4, eagle, or grid. When using a grid architecure, add `--d $(int)` to specify the grid length.
- `--d $(int)`: Grid length of the grid architecture
- `--qf $(str)`: Input QASM file name

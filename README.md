# mOLSQ: Multilevel Quantum Layout Synthesis

## Build and Installation Instructions

```
# Clone ML-QLS repository
git clone git@github.com:UCLA-VAST/ML-QLS.git
cd ML-QLS

# Build
cmake . -Bbuild
cd build
make
```

- Required Dependencies: 
  - [GMP](https://github.com/alisw/GMP/tree/master): Please install GMP. You may need to compile it with `-fPIC` depending on your system.
  - [Bitwuzla](https://github.com/bitwuzla/bitwuzla/tree/4eda0536800576cb2531ab9ce13292da8f21f0eb): Please install Bitwuzla.
  - [pblib](https://github.com/master-keying/pblib): Please install pblib in the `include/`. You may need to compile it with `-fPIC` depending on your system.
  - [pybind11](https://pybind11.readthedocs.io/en/stable/)
  - [Qiskit v1.1.0](https://docs.quantum.ibm.com/start): Please install Qiskit.

## Create a device from the Input Coupling Graph

To perform QLS, we need to know the connections between the qubits, which is information about the physical device.
We are going to use the `createDevice` function to return a C++ object for a QC device.

```
from olsqPy import Device
from src.pyolsq.apiPy import createDevice
device = createDevice(name="dev", nqubits=5, connection=[(0, 1), (1, 2), (1, 3), (3, 4)])
device.printDevice()
```

We also provide three functions `get_nnGrid`, `get_heavyHex`, and `get_device_by_name`, to construct coupling graphs by 1. grid length for grid-based architecture, 2. qubit number for heavy-hexagon-based archtecture,and 3. name, including qx, ourense, rochester, and eagle from IBM, Sycamore from Google, Aspen-4 from Rigetti.
The function will return a list of qubit connections and a device stored in C++ object.

An example to use above functions:
```
from src.pyolsq.device import get_device_by_name, get_nnGrid, get_heavyHex

# construct a 5-by-5 grid architecture
grid_connection, grid_device = get_nnGrid(5)

# construct eagle architecture
eagle_connection, eagle_device = get_device_by_name("eagle")
```

## Create a Circuit from the Input Program

Apart from the device, we need the quantum program/circuit to execute, which can be constructed with the `createCircuit` function.
The function inputs are the name and the qasm string of the circuit.
```
from src.pyolsq.apiPy import createCircuit
circuit_name = "toffoli"
circuit_str = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nh q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[2];\n" \
              "cx q[1], q[2];\ntdg q[2];\ncx q[0], q[2];\nt q[1];\nt q[2];\n" \
              "cx q[0], q[1];\nh q[2];\nt q[0];\ntdg q[1];\ncx q[0], q[1];\n"
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
Idealy, the `run` method will return two solutions, `${circuit_name}_stage_0.qasm` and `${circuit_name}_stage_1.qasm` from sRefine (hueristic QLS algorithm) and the multilevel V cycle, respectively.

## Example: run_mlqls.py

run_mlqls.py is an example program to use ML-QLS to perform layout synthesis. The solutions will be stored under `result/`
```
python3 run_mlqls.py --dt <device_type> --qf <circuit_file>
```
where `<circuit_file>` is the path to circuit qasm file, and `<device_type>` specifies the function used to construct devices. 
We input the name if we want to construct a device by name and input `grid` if we want to constrcut a grid-based architecture using grid length. 
For the second case, we use `--d <N>`, where `<N>` is the grid length, to specify the coupling graph.
In addition, we use `--all_commute` to compile circuit with all commutable gates.
For example, 
```
# compile a QUEKO circuit on a 8-by-8 grid quantum device
python3 run_mlqls.py --dt grid --d 8 --qf benchmark/54QBT_05CYC_QSE_0.qasm

# compile a QUEKO circuit on sycamore quantum device
python3 run_mlqls.py --dt sycamore --qf benchmark/54QBT_05CYC_QSE_0.qasm
```

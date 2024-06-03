def input_qasm(circ_str: str):
    """Process an input qasm string
    A simple qasm parser that works on many qasm files that I know of.
    The current mOLSQ formulation only support single and two-qubitgates,
    so only those are processed.  This parser has some more assumptions:
    1) each instruction in a newline, 2) two-qubit gates: cx or zz,
    and 3) there is only one array of qubits.

    Args:
        circ_str: a string in qasm format.

    Returns:
        a list [qubit_num, gates, gate_spec]
            qubit_num: logical qubit number.
            gates: which qubit(s) each gate acts on.
            gate_spec: type of each gate.

    Example:
        for the following circuit 
            q_0: ───────────────────■───
                                    │  
            q_1: ───────■───────────┼───
                 ┌───┐┌─┴─┐┌─────┐┌─┴─┐
            q_2: ┤ H ├┤ X ├┤ TDG ├┤ X ├─
                 └───┘└───┘└─────┘└───┘ 
        gates = ((2,), (1,2), (2,), (0,1))
        gate_spec = ("h", "cx", "tdg", "cx")
    """

    gates = list() # which qubit(s) each gate acts on
    gate_spec = list() # type of each gate
    qubit_num = 0

    for qasmline in circ_str.splitlines():
        words = qasmline.split()
        if not isinstance(words, list): continue # empty lines
        if len(words) == 0: continue
        grammar = len(words)
        gate_name = qasmline.split('(')[0]
        # grammer=2 -> single-qubit_gate_name one_qubit
        #              |two-qubit_gate_name one_qubit,one_qubit
        # grammer=3 ->  two-qubit_gate_name one_qubit, one_qubit
        # grammer=1 or >3 incorrect, raise an error
        
        if words[0] in ["OPENQASM", "include", "creg", "measure", "//", "barrier"]:
            continue

        if words[0] == 'qreg':
            qubit_num = words[1].split('[')[1]
            # print(qubit_num)
            qubit_num = qubit_num.split(']')[0]
            qubit_num = int(qubit_num)
            continue
        
        if words[0] in ['cx', 'zz', 'u4', 'rzz(pi/4)'] or gate_name in ['rzz']:
            if grammar == 3:
                try:
                    qubit0 = words[1].split('[')[1]
                    # print(qubit_num)
                    qubit0 = qubit0.split(']')[0]
                    qubit0 = int(qubit0)
                    qubit1 = words[2].split('[')[1]
                    # print(qubit_num)
                    qubit1 = qubit1.split(']')[0]
                    qubit1 = int(qubit1)
                    # qubit0 = int(words[1][2:-2])
                    # qubit1 = int(words[2][2:-2])
                except:
                    raise ValueError(f"{qasmline} invalid two-qubit gate.")
            
            elif grammar == 2:
                qubits = words[1].split(',')
                try:
                    # print(qubits)
                    qubit0 = qubits[0].split('[')[1]
                    # print(qubit0)
                    qubit0 = qubit0.split(']')[0]
                    qubit0 = int(qubit0)
                    qubit1 = qubits[1].split('[')[1]
                    # print(qubit_num)
                    qubit1 = qubit1.split(']')[0]
                    qubit1 = int(qubit1)
                    # qubit0 = int(qubits[0][2:-1])
                    # qubit1 = int(qubits[1][2:-2])
                except:
                    raise ValueError(f"{qasmline} invalid two-qubit gate.")
            
            else:
                raise ValueError(f"{qasmline}: invalid two-qubit gate.")
            
            gates.append([qubit0, qubit1])
            gate_spec.append(words[0])
        
        elif grammar == 2: # not two-qubit gate, and has two words
            try:
                # qubit = int(words[1][2:-2])
                qubit = words[1].split('[')[1]
                # print(qubit_num)
                qubit = qubit.split(']')[0]
                qubit = int(qubit)
            except:
                raise ValueError(f"{qasmline} invalid single-qubit gate.")
            gates.append([qubit,])
            gate_spec.append(words[0])
        
        else:
            raise ValueError(f"{qasmline} invalid gate.")

    if qubit_num == 0:
        raise ValueError("Qubit number is not specified.")
    return [qubit_num, tuple(gates), tuple(gate_spec)]

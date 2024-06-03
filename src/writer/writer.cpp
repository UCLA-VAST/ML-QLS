/***********************************************************************
  File        [ molsq.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ planner ]
  Synopsis    [ planner class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2022 ]
***********************************************************************/

#include "writer/writer.hpp"

MOLSQ_NAMESPACE_CPP_START


void Writer::outputQASM(string const & fileName){
    vector< vector<unsigned_t> > vvTimeGate(_circuit.circuitDepth(), vector<unsigned_t>());
    Gate gate;
    unsigned_t i, t, qId, j;
    string line;
    for (i = 0; i < _circuit.nGate(); ++i){
        Gate & gate = _circuit.gate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i);
    }
        for (i = 0; i < _circuit.nSwapGate(); ++i){
        Gate & gate = _circuit.swapGate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i+_circuit.nGate());
    }
    line = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" + to_string(_device.nQubit()) + "];\ncreg c[" + to_string(_device.nQubit()) + "];\n";
    for (t = 0; t < _circuit.circuitDepth(); ++t){
        for (i = 0; i < vvTimeGate[t].size(); ++i){
            if (vvTimeGate[t][i] < _circuit.nGate()){
                gate = _circuit.gate(vvTimeGate[t][i]);
            }
            else{
                gate = _circuit.swapGate(vvTimeGate[t][i] - _circuit.nGate());
            }
            
            if (gate.nTargetQubit() == 1){
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "];\n";
            }
            else{
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "], q[" + to_string(gate.targetPhysicalQubit(1)) + "];\n";
            }
        }
    }
    line = line + "\n// measurement\n";
    for (i = 0; i < _circuit.nProgramQubit(); ++i){
        line = line + "measure q[" + to_string(_circuit.finalMapping(i)) + "]->c[" + to_string(i) + "];\n";
    }
    FILE* fout = fopen(fileName.c_str(), "w");
    fprintf(fout, "%s", line.c_str());
    fclose(fout);
}


string Writer::outputQASMStr(){
    vector< vector<unsigned_t> > vvTimeGate(_circuit.circuitDepth(), vector<unsigned_t>());
    Gate gate;
    unsigned_t i, t, qId, j;
    string line;
    for (i = 0; i < _circuit.nGate(); ++i){
        Gate & gate = _circuit.gate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i);
    }
        for (i = 0; i < _circuit.nSwapGate(); ++i){
        Gate & gate = _circuit.swapGate(i);
        vvTimeGate[gate.executionTime()].emplace_back(i+_circuit.nGate());
    }
    line = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[" + to_string(_device.nQubit()) + "];\ncreg c[" + to_string(_device.nQubit()) + "];\n";
    for (t = 0; t < _circuit.circuitDepth(); ++t){
        for (i = 0; i < vvTimeGate[t].size(); ++i){
            if (vvTimeGate[t][i] < _circuit.nGate()){
                gate = _circuit.gate(vvTimeGate[t][i]);
            }
            else{
                gate = _circuit.swapGate(vvTimeGate[t][i] - _circuit.nGate());
            }
            
            if (gate.nTargetQubit() == 1){
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "];\n";
            }
            else{
                line = line + gate.name() + " q[" + to_string(gate.targetPhysicalQubit(0)) + "], q[" + to_string(gate.targetPhysicalQubit(1)) + "];\n";
            }
        }
    }
    line = line + "\n// measurement\n";
    for (i = 0; i < _circuit.nProgramQubit(); ++i){
        line = line + "measure q[" + to_string(_circuit.finalMapping(i)) + "]->c[" + to_string(i) + "];\n";
    }
    return line.c_str();
}


MOLSQ_NAMESPACE_CPP_END
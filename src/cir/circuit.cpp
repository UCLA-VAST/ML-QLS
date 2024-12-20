/***********************************************************************
  File        [ cir.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ cir ]
  Synopsis    [ Circuit class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "cir/circuit.hpp"

MOLSQ_NAMESPACE_CPP_START

void Circuit::addGate( string const & gateName, vector<unsigned_t> const & vTargetQubit, unsigned_t duration){
    unsigned_t idx = _vGate.size();
    _vGate.emplace_back(idx, gateName, duration);
    for(unsigned_t q: vTargetQubit){
        assert(isValidQubitIdx(q));
        _vGate[idx].addTargetProgramQubit(q);
    }
}

void Circuit::addGate( Gate & gate){
    unsigned_t idx = _vGate.size();
    _vGate.emplace_back(idx, gate.name(), gate.duration());
    for(unsigned_t i = 0; i < gate.nTargetQubit(); ++i){
        _vGate[idx].addTargetProgramQubit(gate.targetProgramQubit(i));
    }
}

void Circuit::addSwapGate(unsigned_t swapIdx, vector<unsigned_t> const & vTargetQubit, unsigned_t duration){
    unsigned_t idx = _vSwapGate.size(), i;
    _vSwapGate.emplace_back(swapIdx, "swap", duration);
    for(i = 0; i < vTargetQubit.size(); ++i){
        // assert(isValidQubitIdx(vTargetQubit[i]));
        _vSwapGate[idx].addTargetProgramQubit(0);
        _vSwapGate[idx].setTargetPhysicalQubit(i, vTargetQubit[i]);
    }
}

void Circuit::addSwapGate(unsigned_t swapIdx, pair<unsigned_t, unsigned_t> const & vTargetQubit, unsigned_t duration){
    unsigned_t idx = _vSwapGate.size(), i;
    _vSwapGate.emplace_back(swapIdx, "swap", duration);
    _vSwapGate[idx].addTargetProgramQubit(0);
    _vSwapGate[idx].addTargetProgramQubit(0);
    _vSwapGate[idx].setTargetPhysicalQubit(0, vTargetQubit.first);
    _vSwapGate[idx].setTargetPhysicalQubit(1, vTargetQubit.second);
}

void Circuit::addSwapGate(unsigned_t swapIdx, unsigned_t q0, unsigned_t q1,  unsigned_t t, unsigned_t duration){
    unsigned_t idx = _vSwapGate.size(), i;
    _vSwapGate.emplace_back(swapIdx, "swap", duration);
    _vSwapGate[idx].addTargetProgramQubit(0);
    _vSwapGate[idx].addTargetProgramQubit(0);
    _vSwapGate[idx].setTargetPhysicalQubit(0, q0);
    _vSwapGate[idx].setTargetPhysicalQubit(1, q1);
    _vSwapGate[idx].setExecutionTime(t);
}

void Circuit::constructDependency(){
    _vpGateDependency.clear();
    vector<int_t> vQubitLastGate(_nProgramQubit, -1);
    Gate gate;
    unsigned_t qId, i, j;
    for (i = 0; i < _vGate.size(); ++i){
        Gate & gate = _vGate[i];
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
                qId = gate.targetProgramQubit(j);
            if (vQubitLastGate[qId] > -1){
                addDependency(vQubitLastGate[qId], i);
            }
            vQubitLastGate[qId] = i;
        }
    }
}

void Circuit::setInitialMapping(vector<unsigned_t> const &  vInitialMapping){
    // assert(vInitialMapping.size() == _nProgramQubit);
    for(unsigned_t i = 0; i < _nProgramQubit; i++){
        // assert(isValidQubitIdx(vInitialMapping[i]));
        _vInitialMapping[i] = vInitialMapping[i];
    }
}

void Circuit::setFinalMapping(vector<unsigned_t> const &  vFinalMapping){
    // assert(vInitialMapping.size() == _nProgramQubit);
    for(unsigned_t i = 0; i < _nProgramQubit; i++){
        // assert(isValidQubitIdx(vInitialMapping[i]));
        _vFinalMapping[i] = vFinalMapping[i];
    }
}

void Circuit::printCircuit(){
    unsigned_t i;
    fprintf(stdout, "[Info] Circuit Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Gate list\n");
    for (Gate &gate: _vGate){
        fprintf(stdout, "        - Gate %d, name: %s, duration: %d, target qubit:", gate.idx(), gate.name().c_str(), gate.duration());
        for ( i = 0; i < gate.nTargetQubit();  ++i )
            fprintf(stdout, " %d", gate.targetProgramQubit(i));
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}


void Circuit::printCircuitLayout(){
    unsigned_t qId, i;
    Gate gate;
    fprintf(stdout, "[Info] Compiled Circuit Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       - #swap            : %lu\n", _vSwapGate.size());
    fprintf(stdout, "       - depth            : %lu\n", _circuitDepth);
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Initial mapping\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "        - Program qubit %d is mapped to physical qubit: %d\n", qId, _vInitialMapping[qId]);
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Final mapping\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "        - Program qubit %d is mapped to physical qubit: %d\n", qId, _vFinalMapping[qId]);
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Gate list\n");
    for (Gate &gate: _vGate){
        fprintf(stdout, "        - Gate %d, name: %s, duration: %d, time: %d, target program qubit:", gate.idx(), gate.name().c_str(), gate.duration(), gate.executionTime());
        for ( i = 0; i < gate.nTargetQubit();  ++i ){
            fprintf(stdout, " %d", gate.targetProgramQubit(i));
        }
        fprintf(stdout, ", target physical qubit: ");
        for ( i = 0; i < gate.nTargetQubit();  ++i ){
            fprintf(stdout, "%d ", gate.targetPhysicalQubit(i));
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - SWAP Gate list\n");
    for (Gate &gate: _vSwapGate){
        fprintf(stdout, "        - SWAP Gate %d, duration: %d, time: %d, target qubit: %d %d\n", gate.idx(), gate.duration(), gate.executionTime(), gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1));
    }
    fprintf(stdout, "\n");
}

void Circuit::printCircuitInitialMapping(){
    unsigned_t qId, i;
    fprintf(stdout, "[Info] Circuit Initial Mapping Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       ------------------------------------------\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "       - Program qubit %d is mapped to physical qubit: %d\n", qId, _vInitialMapping[qId]);
    }
}

void Circuit::printCircuitFinalMapping(){
    unsigned_t qId, i;
    fprintf(stdout, "[Info] Circuit Final Mapping Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - name             : %s\n", _name.c_str());
    fprintf(stdout, "       - #Program Qubit   : %d\n", _nProgramQubit);
    fprintf(stdout, "       ------------------------------------------\n");
    for (qId = 0; qId < _nProgramQubit; ++qId){
        fprintf(stdout, "       - Program qubit %d is mapped to physical qubit: %d\n", qId, _vFinalMapping[qId]);
    }
}


void Circuit::printDependency(){
    unsigned_t i;
    fprintf(stdout, "[Info] Circuit Dependency Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Dependency list\n");
    for ( i = 0; i < _vpGateDependency.size();  ++i ){
        fprintf(stdout, "        - Gate %d depends on gate %d\n", _vGate[_vpGateDependency[i].second].idx(), _vGate[_vpGateDependency[i].first].idx());
    }
}

void Circuit::printQubitRegion(){
    fprintf(stdout, "[Info] Qubit Region Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    unsigned_t i;
    for ( i = 0; i < _vsQubitRegion.size();  ++i ){
        fprintf(stdout, "        - Qubit %d: ", i);
        for(int_t q : _vsQubitRegion[i]){
            fprintf(stdout, "%d ", q);
        }
        fprintf(stdout, "\n");
    }
}

MOLSQ_NAMESPACE_CPP_END

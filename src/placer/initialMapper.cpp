/***********************************************************************
  File        [ olsq.cpp ]
  System      [ nOLSQ multilevel  quantum layout synthesis tool]
  Package     [ clusterer ]
  Synopsis    [ initialMapper class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "placer/initialMapper.hpp"

MOLSQ_NAMESPACE_CPP_START

bool InitialMapper::run(Circuit& cir, Device& device, unsigned_t solNum){
    _timer.start(TimeUsage::FULL);
    _solNum = solNum;
    _pCircuit = &cir;
    _pDevice = &device;
    _timer.setTimeout(_initialMapperParam.timeout);
    _hasSol = false;
    // if circuit has already gone through pure heuristic stage, use the current sol as the best disdance.
    if(_pCircuit->circuitDepth() > 0){
        _bestDisSum = 0;
        for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
            Gate & gate = _pCircuit->gate(i);
            if(gate.nTargetQubit() > 1){
                _bestDisSum += _pDevice->getDistance(_pCircuit->initialMapping(gate.targetProgramQubit(0)), _pCircuit->initialMapping(gate.targetProgramQubit(1)));
            }
        }
    }
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Find initial mapping which allows most gate execution                        \n");
    bool result = runSMT();
    // cerr << "result: " << result << endl;
    sort(_vvQubitMapingsolution.begin(), _vvQubitMapingsolution.end());
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Set circuit mapping (cost: %.2f) \n", _vvQubitMapingsolution[0].first);
    for (unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _pCircuit->setInitialMapping(i, _vvQubitMapingsolution[0].second[i]);
    }
    if(result){
        // cerr << "mapGateBasedOnInitialMapping" << endl;
        if(_verbose > 0)
            fprintf(stdout, "[Info] InitialMapper: Map gates based on initial mapping\n");
        mapGateBasedOnInitialMapping();
    }
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Finish\n");
    return result;
}

bool InitialMapper::runSMT(){
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Formulation generation                        \n");
    generateFormulation();
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Optimization                        \n");
    bool result = optimize();
    return result;
}

void InitialMapper::generateFormulation(){
    _smt.reset(_initialMapperParam.timeout);
    if(_verbose > 0)
        fprintf(stdout, "[Info]          constructing variables                       \n");
    constructVariable();
    if(_verbose > 0)
        fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints();
    if(_initialMapperParam.is_all_commute){
        for (unsigned_t i = 0; i < _pCircuit->nGate();  ++i ){
            if(_pCircuit->gate(i).nTargetQubit() == 2){
                addValidTwoQubitGateConstraints(i);
            }
        }
    }
}

void InitialMapper::constructVariable(){
    _smt.piSort = mk_bv_sort(ceil(log2(_pDevice->nQubit() + 1)));
    _smt.vPi.reserve(_pCircuit->nProgramQubit());
    unsigned_t i;
    string s;
    Sort sortbool = mk_bool_sort();
    if(_initialMapperParam.is_all_commute){
        _smt.vGate.reserve(_pCircuit->nGate());
        for (i = 0; i < _pCircuit->nGate(); ++i){
            if(_pCircuit->gate(i).nTargetQubit() == 2){
                s = "exe_g" + to_string(i);
                _smt.vGate.emplace_back(mk_const(sortbool, s.c_str()));        
            }
        }
    }
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        s = "map_q" + to_string(i);
        _smt.vPi.emplace_back(mk_const(_smt.piSort, s.c_str()));
    }
}

void InitialMapper::addInjectiveMappingConstraints(){
    unsigned_t i, j;
    const Term zero = mk_bv_zero(_smt.piSort);
    const Term nqubit = mk_bv_value_uint64(_smt.piSort, _pDevice->nQubit());
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vPi[i]}));
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {_smt.vPi[i], nqubit}));
        for (j = 0; j < i; ++j){
             _smt.pSolver->assert_formula(mk_term(Kind::DISTINCT, {_smt.vPi[i], _smt.vPi[j]}));
        }
    }
}

void InitialMapper::addValidTwoQubitGateConstraints(unsigned_t gateId){
    unsigned_t i, t, j;
    Term pro1EqPhy1, pro2EqPhy2, pro1EqPhy2, pro2EqPhy1, cond1, cond2, q1Bv, q2Bv;
    Gate & gate = _pCircuit->gate(gateId);
    Edge & edge = _pDevice->edge(0);
    q1Bv = mk_bv_value_uint64(_smt.piSort, edge.qubitId1());
    q2Bv = mk_bv_value_uint64(_smt.piSort, edge.qubitId2());
    pro1EqPhy1 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(0)], q1Bv});
    pro2EqPhy2 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(1)], q2Bv});
    pro1EqPhy2 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(0)], q2Bv});
    pro2EqPhy1 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(1)], q1Bv});
    cond1 = mk_term(Kind::AND, {pro1EqPhy1, pro2EqPhy2});
    cond2 = mk_term(Kind::AND, {pro2EqPhy1, pro1EqPhy2});
    Term clause = mk_term(Kind::OR, {cond1, cond2});
    
    for ( i = 1; i < _pDevice->nEdge();  ++i ){
        Edge & edge = _pDevice->edge(i);
        q1Bv = mk_bv_value_uint64(_smt.piSort, edge.qubitId1());
        q2Bv = mk_bv_value_uint64(_smt.piSort, edge.qubitId2());
        pro1EqPhy1 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(0)], q1Bv});
        pro2EqPhy2 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(1)], q2Bv});
        pro1EqPhy2 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(0)], q2Bv});
        pro2EqPhy1 = mk_term(Kind::EQUAL, {_smt.vPi[gate.targetProgramQubit(1)], q1Bv});
        cond1 = mk_term(Kind::AND, {pro1EqPhy1, pro2EqPhy2});
        cond2 = mk_term(Kind::AND, {pro2EqPhy1, pro1EqPhy2});
        clause = mk_term(Kind::OR, {clause, cond1, cond2});
    }                       
    if(_initialMapperParam.is_all_commute){
        clause = mk_term(Kind::IFF, {mk_term(Kind::NOT, {_smt.vGate[gateId]}), clause});
    }
     _smt.pSolver->assert_formula(clause);

}

void InitialMapper::addGateCountConstraints(unsigned_t bound){
    PB2CNF pb2cnf;
    vector<int_t> vSigmaLit;
    unsigned_t firstFreshVariable = 1 + _smt.vGate.size(), newFreashVariable;
    for (unsigned_t i = 1; i < firstFreshVariable; ++i){
        vSigmaLit.emplace_back(i);
    }
    vector< vector<int_t> > formula;
    newFreashVariable = pb2cnf.encodeAtMostK(vSigmaLit, bound, formula, firstFreshVariable) + 1;
    map<unsigned_t, Term> mAncillary;
    vector<Term> vOrs;
    unsigned_t var;
    Sort sortbool = mk_bool_sort();
    string s;
    // unsigned_t length = 0;
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        // length += clause.size();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                if (lit < 0){ 
                    vOrs.emplace_back(mk_term(Kind::NOT, {_smt.vGate[var - 1]}));
                }
                else{
                    vOrs.emplace_back(_smt.vGate[var - 1]);
                }
            }
            else{
                if(mAncillary.find(var) == mAncillary.end()){
                    s = "anc_" + to_string(var);
                    mAncillary[var] = mk_const(sortbool, s.c_str());
                }
                if (lit < 0){
                    vOrs.emplace_back(mk_term(Kind::NOT, {mAncillary[var]}));
                }
                else{
                    vOrs.emplace_back(mAncillary[var]);
                }
            }
        }    
        Term cnf = vOrs[0];
        for (unsigned_t i = 1; i < vOrs.size(); ++i){
            cnf = mk_term(Kind::OR, {vOrs[i], cnf});
        }
         _smt.pSolver->assert_formula(cnf);
    }
}

bool InitialMapper::checkModel(){
    Result status = _smt.pSolver->check_sat();
    if (status == Result::SAT){
        return true;
    }
    else{
        return false;
    }
}

bool InitialMapper::optimize(){
    bool success;
    unsigned_t gateCount = 0, failCount = 0;
    int_t upperBound = _pCircuit->nGate(), lowerBound = 0, bound = (lowerBound+upperBound)/2;
    if(_initialMapperParam.is_all_commute){
        // binary search for largest gate count
        while(lowerBound < upperBound && lowerBound >= 0){
            _smt.pSolver->push(1);
            addGateCountConstraints(bound);
            if(_verbose > 0)
                fprintf(stdout, "[Info] InitialMapper: Trying bound: %d, Current best distance sum: %d \r", bound, _bestDisSum);
            success = checkModel();
            if(success){
                _hasSol = true;
                extractModel();
                upperBound = bound;
            }
            else{
                lowerBound = bound + 1;
            }
            _smt.pSolver->pop(1);
            bound = (lowerBound+upperBound)/2;
        }
        if(upperBound > 0){
            addGateCountConstraints(upperBound);
        }
    }
    else{
        // iteratively add gate into the model
        vector<unsigned_t> vGateOrder(_pCircuit->nGate());
        for (unsigned_t i = 0; i < _pCircuit->nGate();  ++i ){
            vGateOrder[i] = i;
        }
        random_shuffle( vGateOrder.begin(), vGateOrder.end() );
        for (unsigned_t i = 0; i < _pCircuit->nGate();  ++i ){
            if(_pCircuit->gate(vGateOrder[i]).nTargetQubit() == 2){
                _smt.pSolver->push(1);
                addValidTwoQubitGateConstraints(vGateOrder[i]);
                success = checkModel();
                if(success){
                    _hasSol = true;
                    extractModel();
                    ++gateCount;
                }
                else{
                    ++failCount;
                    _smt.pSolver->pop(1);
                }
            }
        }
    }
    if(failCount > 0){
        sampleSolution();
    }
    if(_initialMapperParam.is_all_commute){
        if(_verbose > 0)
            fprintf(stdout, "[Info] Successfully-mapped gate count: %d, Unsuccessfully-mapped gate count: %d\n", _pCircuit->nGate() - upperBound, upperBound);
        return (upperBound == 0);
    }
    else{
        if(_verbose > 0)
            fprintf(stdout, "[Info] Successfully-mapped gate count: %d, Unsuccessfully-mapped gate count: %d\n", gateCount, failCount);
        return (failCount == 0);
    }
    return false;
}


void InitialMapper::sampleSolution(){
    _smt.pSolver->push(1);
    unsigned_t solCnt = 0, solCntBound = 100;
    bool success = true;
    
    while(success && solCnt < solCntBound){
        success = checkModel();
        if(success){
            extractModel(1);
        }
        ++solCnt;
    }
    _smt.pSolver->pop(1);
    return;
}

void InitialMapper::blockSolution(vector<unsigned_t> & vMapping){
    unsigned_t t, i;
    Term qBv = mk_bv_value_uint64(_smt.piSort, vMapping[0]);
    Term clause = mk_term(Kind::EQUAL, {_smt.vPi[0], qBv});
    // expr clause = !(_smt.vvPi[0][0] == (int_t)(m.eval(_smt.vvPi[0][0]).get_numeral_int64()));   
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        qBv = mk_bv_value_uint64(_smt.piSort, vMapping[i]);
        clause = mk_term(Kind::AND, {clause,
                        mk_term(Kind::EQUAL, {_smt.vPi[i], qBv})});
    }
     _smt.pSolver->assert_formula(mk_term(Kind::NOT, {clause}));
}


void InitialMapper::extractModel(bool blockSol){
    unsigned_t circuitDepth = 0, i;
    double_t disSum = 0;
    string s;

    // collect qubit mapping
    vector<unsigned_t> vMapping(_pCircuit->nProgramQubit());
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        vMapping[i] = stoi(_smt.pSolver->get_value(_smt.vPi[i]).value<string>(10));
    }
    if(blockSol){
        blockSolution(vMapping);
    }

    for (unsigned_t i = 0; i < _pCircuit->nGate();  ++i ){
        Gate & gate = _pCircuit->gate(i);
        if(gate.nTargetQubit() == 2){
            // cerr << "add distance: " << _pDevice->getDistance(vMapping[gate.targetProgramQubit(0)], vMapping[gate.targetProgramQubit(1)]) << endl;
            disSum += (_pDevice->getDistance(vMapping[gate.targetProgramQubit(0)], vMapping[gate.targetProgramQubit(1)]));
        }
    }

    // if(disSum < _bestDisSum){
    //     _bestDisSum = disSum;
    //     for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
    //         _pCircuit->setInitialMapping(i, vMapping[i]);
    //     }
        // if (_verbose > 0){
        //     fprintf(stdout, "[Info] Extract Qubit Mapping (cost: %d):\n", disSum);
        //         for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //             fprintf(stdout, "%d->%d ", i, vMapping[i]);
        //     }
        //     fprintf(stdout, "\n");
        // }
    // }
    if(_vvQubitMapingsolution.size() < _solNum){
        _vvQubitMapingsolution.emplace_back(disSum, vMapping);
    }
    else{
        unsigned_t maxIdx = 0;
        for(i = 1; i < _vvQubitMapingsolution.size(); ++i){
            if(_vvQubitMapingsolution[maxIdx].first < _vvQubitMapingsolution[i].first){
                maxIdx = i;
            }
        }
        if(_vvQubitMapingsolution[maxIdx].first > disSum){
            // cerr << "update sol with dis " << disSum << endl;
            _vvQubitMapingsolution.erase(_vvQubitMapingsolution.begin()+maxIdx); 
            _vvQubitMapingsolution.emplace_back(disSum, vMapping);
            // cerr << "_vvQubitMapingsolution.size() " << _vvQubitMapingsolution.size() << endl;
        }
    }
}

void InitialMapper::mapGateBasedOnInitialMapping(){
    for (unsigned_t i = 0; i < _pCircuit->nGate();  ++i ){
        Gate& gate = _pCircuit->gate(i);
        gate.setExecutionTime(0);
        for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
            gate.setTargetPhysicalQubit(j, _pCircuit->initialMapping(gate.targetProgramQubit(j)));
        }
    }
    for (unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _pCircuit->setFinalMapping(i, _pCircuit->initialMapping(i));
        _pCircuit->addQubitRegion(i, _pCircuit->initialMapping(i));
    }
    _pCircuit->setCircuitDepth(1);
}

void InitialMapper::runPostprocessing(){
    // Placer placer;
    Placer_dev placer;
    placer.run(*_pCircuit, *_pDevice, _hasSol);
    _bestDisSum = placer.getOptimalCost();
    return;
}

void InitialMapper::setCircuitIntialMappingBySolutionIdx(unsigned_t idx, Circuit & cir){
    if(_verbose > 0)
        fprintf(stdout, "[Info] InitialMapper: Set circuit mapping (cost: %.2f) \n", _vvQubitMapingsolution[idx].first);
    for (unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        cir.setInitialMapping(i, _vvQubitMapingsolution[idx].second[i]);
    }
}

MOLSQ_NAMESPACE_CPP_END
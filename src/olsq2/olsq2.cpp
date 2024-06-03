/***********************************************************************
  File        [ olsq.cpp ]
  System      [ OLSQ2: optimal quantum layout synthesis tool]
  Package     [ olsq ]
  Synopsis    [ OLSQ2 class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "olsq2/olsq2.hpp"

MOLSQ_NAMESPACE_CPP_START

void OLSQ2::setDependency(vector<pair<unsigned_t, unsigned_t> >& vDependencies){
    for (pair<unsigned_t, unsigned_t>& p : vDependencies){
        assert(isValidGateIdx(p.first));
        assert(isValidGateIdx(p.second));
        addDependency(p.first, p.second);
    }
}

bool OLSQ2::run(){
    _timer.start(TimeUsage::FULL);
    _hasSol = false;
    _timer.setTimeout(_olsqParam.timeout/1000);
    fprintf(stdout, "[Info] OLSQ2 Layout Synthesis                        \n");
    if(!_olsqParam.is_given_dependency && !_olsqParam.is_all_commute){
        constructDependency();
        fprintf(stdout, "[Info] Constructing dependency                        \n");
        if(_verbose == 2){
            printDependency();
        }
    }
    if(_olsqParam.is_all_commute && !_olsqParam.is_transition){
        constructCollisionList();
    }
    bool solve = runSMT();
    if (_olsqParam.is_transition && !_olsqParam.is_multilevel){
        asapScheduling();
    }
    if(_verbose > 0){
        _pCircuit->printCircuitLayout();
    }
    return solve;
}

void OLSQ2::dump(){
    // for formulation generation
    _olsqParam.min_depth = 13;
    increaseDepthBound();
    constructDependency();
    if (!_olsqParam.is_transition && _olsqParam.use_window_range_for_gate){
        constructGateTimeWindow();
    }
    fprintf(stdout, "[Info] Generating formulation                        \n");
    generateFormulation();
    addDepthConstraints();
    unsigned_t swap_bound = 19;
    addSwapCountConstraints(swap_bound);
    string fileName = to_string(_device.nQubit()) + "_" + to_string(_pCircuit->nGate()) + "_" + to_string(_olsqParam.min_depth) + "_" + to_string(swap_bound) + "_timewindow_wo_card.txt";
    std::ofstream fout(fileName);
    _smt.pSolver->print_formula(fout);
}

bool OLSQ2::runSMT(){
    fprintf(stdout, "[Info] OLSQ2 Layout Synthesis                        \n");
    if(!_olsqParam.is_transition){
        if (!_olsqParam.is_given_depth){
            _olsqParam.min_depth = extract_longest_chain();
            increaseDepthBound();
            fprintf(stdout, "[Info] Longest chain = %d\n", _olsqParam.min_depth);
        }
    }
    increaseDepthBound();
    bool solve = false;
    _iter = 0;
    if (!_olsqParam.is_transition && _olsqParam.use_window_range_for_gate){
        constructGateTimeWindow();
    }
    while (!solve && !_timer.isTimeout()){
        fprintf(stdout, "[Info] Iter %d: Solving with depth range (%d, %d)            \n", _iter, _olsqParam.min_depth, _olsqParam.max_depth);
        fprintf(stdout, "[Info] Iter %d: Generating formulation                        \n", _iter);
        generateFormulation();
        fprintf(stdout, "[Info] Iter %d: Optimizing model                             \n", _iter);
        solve = optimize();
        if(!solve){
            increaseDepthBound();
        }
        ++_iter;
    }
    return solve;
}

void OLSQ2::generateFormulation(){
    _smt.reset(_olsqParam.timeout, _iter);
    fprintf(stdout, "[Info]          constructing variables                       \n");
    constructVariable();
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints();
    fprintf(stdout, "[Info]          constructing dependency constraint           \n");
    addDependencyConstraints();
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addValidTwoQubitGateConstraints();
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraints();
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraints();
}

void OLSQ2::constructVariable(){
    unsigned_t bit_length_pi, bit_length_time;
    bit_length_pi = ceil(log2(_device.nQubit() + 1));
    bit_length_time = _olsqParam.max_depth_bit;
    _smt.vvPi.reserve(_olsqParam.max_depth);
    _smt.vTg.reserve(_pCircuit->nGate());
    _smt.vvSigma.reserve(_olsqParam.max_depth);

    string s;
    unsigned_t i, j;

    // Create a bit-vector sort of size 1.
    Sort sortbvpi = mk_bv_sort( bit_length_pi);
    Sort sortbvtime = mk_bv_sort( bit_length_time);
    Sort sortbool = mk_bool_sort();

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvPi.emplace_back(vector<Term> ());
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            s = "map_t" + to_string(i) + "_q" + to_string(j);
            _smt.vvPi[i].emplace_back(mk_const(sortbvpi, s.c_str()));
        }
    }

    for (i = 0; i < _pCircuit->nGate(); ++i){
        s = "time_" + to_string(i);
        _smt.vTg.emplace_back(mk_const(sortbvtime, s.c_str()));
    }

    for (i = 0; i < _olsqParam.max_depth; ++i){
        _smt.vvSigma.emplace_back(vector<Term> ());
        for (j = 0; j < _device.nEdge(); ++j){
            s = "ifswap_t" + to_string(i) + "_e" + to_string(j);
            _smt.vvSigma[i].emplace_back(mk_const(sortbool, s.c_str()));
        }
    }

}

void OLSQ2::addInjectiveMappingConstraints(unsigned_t boundOffset){
    unsigned_t i, j, k;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    Sort sortbvpi = mk_bv_sort( ceil(log2(_device.nQubit() + 1)));
    //
    // if(boundOffset == 0){
        // vector<unsigned_t> vMapping = {6, 15, 13, 0, 10, 7, 1, 12, 5, 11, 9, 14, 4, 3, 2, 8};
        // vector<unsigned_t> vMapping = {14, 2, 17,4, 22, 5, 18, 19, 0, 11, 9, 12, 1, 24, 16, 13, 8, 7, 6, 3, 10, 15, 21, 23};  // a* 13 swap
        // vector<unsigned_t> vMapping = {0, 13, 7, 22, 21, 14, 3, 15, 19, 8, 6, 2, 11, 5, 17, 4, 1, 18, 24, 23, 9, 12, 10, 16};  // a* 28 swap
        
    //     for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
    //         Term mapping = mk_bv_value_uint64(sortbvpi, vMapping[i]);
    //         _smt.pSolver->assert_formula(mk_term(Kind::EQUAL, {_smt.vvPi[0][i], mapping}));
    //     }
    // }
    //
    Term zero = mk_bv_zero(sortbvpi);
    Term nqubit = mk_bv_value_uint64(sortbvpi, _device.nQubit());
    for (i = begin; i < end; ++i){
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vvPi[i][j]}));
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {_smt.vvPi[i][j], nqubit}));
            for (k = 0; k < j; ++k){
                 _smt.pSolver->assert_formula(mk_term(Kind::DISTINCT, {_smt.vvPi[i][j], _smt.vvPi[i][k]}));
            }
        }
    }
}

void OLSQ2::addValidTwoQubitGateConstraints(unsigned_t boundOffset){
    unsigned_t i, t, j;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    Sort sortbvtime = mk_bv_sort( _olsqParam.max_depth_bit);
    Sort sortbvpi = mk_bv_sort( ceil(log2(_device.nQubit() + 1)));
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for ( i = 0; i < _pCircuit->nGate();  ++i ){
        Gate & gate = _pCircuit->gate(i);
        if ((gate).nTargetQubit() == 2){
            if(_olsqParam.use_window_range_for_gate){
                begin = _vpGateTimeWindow[i].first;
                end = _vpGateTimeWindow[i].second + 1;
                if(boundOffset > 0){
                    begin = end;
                    end += boundOffset;
                }
                // cerr << "construct valid two qubit constraint for gate "<< i <<" from " << begin << " to " << end << endl;
            }
            // cerr << "construct valid two qubit constraint for gate "<< i << endl;
            for (t = begin; t < end; ++t){
                Term bvt = mk_bv_value_uint64(sortbvtime, t);
                Term clause = mk_term(Kind::NOT, {mk_term(Kind:: EQUAL, {_smt.vTg[i], bvt})});
                Term pro1EqPhy1, pro2EqPhy2, pro1EqPhy2, pro2EqPhy1, cond1, cond2;
                for ( j = 0; j < _device.nEdge();  ++j ){
                    Edge & edge = _device.edge(j);
                    Term q1Bv = mk_bv_value_uint64(sortbvpi, edge.qubitId1());
                    Term q2Bv = mk_bv_value_uint64(sortbvpi, edge.qubitId2());
                    pro1EqPhy1 = mk_term(Kind:: EQUAL, {_smt.vvPi[t][gate.targetProgramQubit(0)], q1Bv});
                    pro2EqPhy2 = mk_term(Kind:: EQUAL, {_smt.vvPi[t][gate.targetProgramQubit(1)], q2Bv});
                    pro1EqPhy2 = mk_term(Kind:: EQUAL, {_smt.vvPi[t][gate.targetProgramQubit(0)], q2Bv});
                    pro2EqPhy1 = mk_term(Kind:: EQUAL, {_smt.vvPi[t][gate.targetProgramQubit(1)], q1Bv});
                    cond1 = mk_term(Kind::AND, {pro1EqPhy1, pro2EqPhy2});
                    cond2 = mk_term(Kind::AND, {pro2EqPhy1, pro1EqPhy2});
                    clause = mk_term( Kind::OR, {clause, cond1, cond2});
                }                       
                 _smt.pSolver->assert_formula(clause);
            }
        }
    }    
}

void OLSQ2::addDependencyConstraints(){
    unsigned_t i;
    Gate g;
    Sort sortbvtime = mk_bv_sort( _olsqParam.max_depth_bit);
    Term zero = mk_bv_zero(sortbvtime);
    for (i = 0; i < _pCircuit->nGate(); ++i){
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vTg[i]}));
    }
    if (_olsqParam.is_transition){
        if(!_olsqParam.is_all_commute){
            for ( i = 0; i < _vpGateDependency.size();  ++i ){
                 _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {_smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]}));
            }
        }
    }
    else if (_olsqParam.is_all_commute){
        for ( i = 0; i < _vGateCollision.size();  ++i ){
             _smt.pSolver->assert_formula(mk_term(Kind::DISTINCT, {_smt.vTg[_vGateCollision[i].first], _smt.vTg[_vGateCollision[i].second]}));
        }
    }
    else{
        for ( i = 0; i < _vpGateDependency.size();  ++i ){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {_smt.vTg[_vpGateDependency[i].first], _smt.vTg[_vpGateDependency[i].second]}));
        }
    }
}

void OLSQ2::addSwapConstraints(unsigned_t boundOffset){
    // No swap for t<s
    unsigned_t i, j, t, e, tt, q1, q2;
    Qubit qubit;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    unsigned_t bound = (_olsqParam.swap_duration < _olsqParam.max_depth) ? _olsqParam.swap_duration: _olsqParam.max_depth;
    if(_olsqParam.is_transition){
        --bound;
    }
    // unsigned_t begin = _olsqParam.swap_duration -1, end = _olsqParam.min_depth;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    else{
        for (i = 0; i < bound; ++i){
            for (j = 0; j < _device.nEdge(); ++j){
                 _smt.pSolver->assert_formula(mk_term(Kind::NOT, {_smt.vvSigma[i][j]}));
            }
        }
    }
    // cout << "begin: " << begin << ", end: " << end << endl;
    // swap gates can not overlap with swap in space
    for (t = begin; t < end; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            q1 = _device.edge(e).qubitId1();
            q2 = _device.edge(e).qubitId2();
            for (unsigned_t ee : _device.qubit(q1).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    if (ee < e){
                         _smt.pSolver->assert_formula(
                            mk_term(Kind::OR, 
                                {mk_term(Kind::NOT, {_smt.vvSigma[t][e]}),
                                mk_term(Kind::NOT, {_smt.vvSigma[tt][ee]})}));
                    }
                }
            }
            for (unsigned_t ee : _device.qubit(q2).vSpanEdge){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t + 1; ++tt){
                    if (ee < e){
                         _smt.pSolver->assert_formula(
                            mk_term(Kind::OR, 
                                {mk_term(Kind::NOT, {_smt.vvSigma[t][e]}),
                                mk_term(Kind::NOT, {_smt.vvSigma[tt][ee]})}));
                    }
                }
            }
        }
    }
    if (!_olsqParam.is_transition){
        begin = 0;
        end = _olsqParam.min_depth;
        Sort sortbvtime = mk_bv_sort( _olsqParam.max_depth_bit);
        Sort sortbvpi = mk_bv_sort( ceil(log2(_device.nQubit() + 1)));
        // swap gates can not overlap with swap in time
        for (t = begin; t < end; ++t){
            for (e = 0; e < _device.nEdge(); ++e){
                for (tt = t - _olsqParam.swap_duration + 1; tt < t; ++tt){
                     _smt.pSolver->assert_formula(
                            mk_term(Kind::OR, 
                                {mk_term(Kind::NOT, {_smt.vvSigma[t][e]}),
                                mk_term(Kind::NOT, {_smt.vvSigma[tt][e]})}));
                }
            }
        }
        // swap gates can not ovelap with other gates
        // the time interval should be modified
        for (i = 0; i < _pCircuit->nGate(); ++i){
            if(_olsqParam.use_window_range_for_gate){
                begin = _vpGateTimeWindow[i].first;
                end = _vpGateTimeWindow[i].second + 1;
                if(boundOffset > 0){
                    begin = end;
                    end += boundOffset;
                }
            }
            for (t = begin; t < end; ++t){
                for (e = 0; e < _device.nEdge(); ++e){
                    Gate& gate = _pCircuit->gate(i);
                    for (tt = t; tt < t + _olsqParam.swap_duration; ++tt){
                        Term bvt = mk_bv_value_uint64(sortbvtime, tt);
                        Term q1Bv = mk_bv_value_uint64(sortbvpi, _device.edge(e).qubitId1());
                        Term q2Bv = mk_bv_value_uint64(sortbvpi, _device.edge(e).qubitId2());
                        Term clause = mk_term(Kind:: EQUAL, {_smt.vTg[i], bvt});
                        Term clauseOr = mk_term(Kind::OR, 
                                                    {mk_term(Kind:: EQUAL, {_smt.vvPi[tt][gate.targetProgramQubit(0)], q1Bv}), //pro1EqPhy1,
                                                    mk_term(Kind:: EQUAL, {_smt.vvPi[tt][gate.targetProgramQubit(0)], q2Bv})});
                        Term clauseAnd;
                        Term cond;
                        cond = mk_term(Kind::NOT, {_smt.vvSigma[t][e]});
                        if (gate.nTargetQubit() == 1){
                            clauseAnd =mk_term(Kind::AND, {clause, clauseOr});
                            clauseAnd = mk_term(Kind::NOT, {clauseAnd});
                             _smt.pSolver->assert_formula(
                                mk_term(Kind::OR, {clauseAnd, cond}));
                        }
                        else if(gate.nTargetQubit() == 2){
                            clauseOr = mk_term( Kind::OR, 
                                                    {mk_term(Kind:: EQUAL, {_smt.vvPi[tt][gate.targetProgramQubit(1)], q1Bv}), // pro2EqPhy1,
                                                    mk_term(Kind:: EQUAL, {_smt.vvPi[tt][gate.targetProgramQubit(1)], q2Bv}), // pro2EqPhy2,
                                                    clauseOr});
                            clauseAnd =mk_term(Kind::AND, {clause, clauseOr});
                            clauseAnd = mk_term(Kind::NOT, {clauseAnd});
                             _smt.pSolver->assert_formula(
                                mk_term(Kind::OR, {clauseAnd, cond}));
                        }
                        else{
                            assert(false);
                        }
                    } 
                } 
            }
        }
    }
}

void OLSQ2::addTransformationConstraints(unsigned_t boundOffset){
    unsigned_t i, j, t, e, nSpanEdge;
    unsigned_t end = _olsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    // Mapping Not Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Not Transformations by SWAP Gates.                       \n");
    Sort sortbvpi = mk_bv_sort( ceil(log2(_device.nQubit() + 1)));
    for (t = begin; t < end; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (j = 0; j < _device.nQubit(); ++j){
                Qubit& qubit = _device.qubit(j);
                nSpanEdge = qubit.vSpanEdge.size();
                Term clause = _smt.vvSigma[t][qubit.vSpanEdge[0]];
                for (e = 1; e < nSpanEdge; ++e){
                    clause = mk_term(Kind::OR, {clause, _smt.vvSigma[t][qubit.vSpanEdge[e]]});
                }
                Term qBv = mk_bv_value_uint64(sortbvpi, j);
                clause = mk_term(Kind::AND,
                            {mk_term(Kind::NOT, {clause}), 
                            mk_term(Kind:: EQUAL, {_smt.vvPi[t][i], qBv})});
                 _smt.pSolver->assert_formula(
                                mk_term(Kind::OR,
                                    {mk_term(Kind::NOT, {clause}), 
                                    mk_term(Kind:: EQUAL, {_smt.vvPi[t+1][i], qBv})}));
            }
        }
    }
    // Mapping Transformations by SWAP Gates.
    // fprintf(stdout, "[Info]          Mapping Transformations by SWAP Gates.                       \n");
    for (t = begin; t < end; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            for (e = 0; e < _device.nEdge(); ++e){
                Term q1Bv = mk_bv_value_uint64(sortbvpi, _device.edge(e).qubitId1());
                Term q2Bv = mk_bv_value_uint64(sortbvpi, _device.edge(e).qubitId2());
                Term cond = mk_term(Kind::NOT,
                                        {mk_term(Kind::AND,
                                            {_smt.vvSigma[t][e], 
                                            mk_term(Kind:: EQUAL, {_smt.vvPi[t][i], q1Bv})})});
                 _smt.pSolver->assert_formula(mk_term(Kind::OR,
                                                {cond, 
                                                mk_term(Kind:: EQUAL, {_smt.vvPi[t+1][i], q2Bv})}));
                cond = mk_term(Kind::NOT,
                        {mk_term(Kind::AND,
                        {_smt.vvSigma[t][e], 
                        mk_term(Kind:: EQUAL, {_smt.vvPi[t][i], q2Bv})})});
                 _smt.pSolver->assert_formula(mk_term(Kind::OR,
                                                {cond, 
                                                mk_term(Kind:: EQUAL, {_smt.vvPi[t+1][i], q1Bv})}));
            }
        }
    }
}

void OLSQ2::addDepthConstraints(){
    Sort sortbvtime = mk_bv_sort( _olsqParam.max_depth_bit);
    if(!_olsqParam.use_window_range_for_gate){
        Term depthBv = mk_bv_value_uint64(sortbvtime, _olsqParam.min_depth);
        for (Term& t: _smt.vTg){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {t, depthBv}));
        }
    }
    else{
        for(unsigned_t i = 0; i < _vpGateTimeWindow.size(); ++i){
        // for (const Term t: _smt.vTg){
            Term depthBv = mk_bv_value_uint64(sortbvtime, _vpGateTimeWindow[i].second);
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE,  {_smt.vTg[i], depthBv}));
        }
    }
}

void OLSQ2::addSwapCountConstraints(unsigned_t bound){
    PB2CNF pb2cnf;
    vector<int_t> vSigmaLit;
    unsigned_t firstFreshVariable = 1 + _olsqParam.min_depth*_device.nEdge(), newFreashVariable;
    for (unsigned_t i = 1; i < firstFreshVariable; ++i){
        vSigmaLit.emplace_back(i);
    }
    vector< vector<int_t> > formula;
    newFreashVariable = pb2cnf.encodeAtMostK(vSigmaLit, bound, formula, firstFreshVariable) + 1;
    map<unsigned_t, Term> mAncillary;
    vector<Term> vOrs;
    unsigned_t sigmaT, sigmaE, var;
    Sort sortbool = mk_bool_sort();
    string s;
    // unsigned_t length = 0;
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        // length += clause.size();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                sigmaT = (var - 1) / _device.nEdge();
                sigmaE = (var - 1) % _device.nEdge();
                if (lit < 0){ 
                    vOrs.emplace_back(mk_term(Kind::NOT, {_smt.vvSigma[sigmaT][sigmaE]}));
                }
                else{
                    vOrs.emplace_back(_smt.vvSigma[sigmaT][sigmaE]);
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

bool OLSQ2::checkModel(){
    Result status= _smt.pSolver->check_sat();
    // cout << status << endl;
    // cout << _smt.smtSolver.reason_unknown() << endl;
    // getchar();
    // return true;
    if (status == Result::SAT){
        return true;
    }
    else{
        return false;
    }
}

bool OLSQ2::optimize(){
    bool success_optimize;
    if(!_olsqParam.is_optimize_swap){
        success_optimize = optimizeDepth();
    }
    else{
        success_optimize = optimizeSwap();
    }
    return success_optimize;
}

bool OLSQ2::optimizeDepth(){
    bool success, find_min_depth = false, has_jump = false;
    unsigned_t step;
    unsigned_t i;
    if (_olsqParam.is_transition){
        step = 1;
    }
    else{
        step = (_olsqParam.min_depth > 100) ? 10 : 1;
    }
    while(!find_min_depth && _olsqParam.min_depth < _olsqParam.max_depth ){
        fprintf(stdout, "[Info]          trying to optimize for depth bound %d            \n", _olsqParam.min_depth);
        _timer.start(TimeUsage::PARTIAL);
        if(!_hasSol){
            _smt.resetTimeState(_olsqParam.timeout);
        }
        _smt.pSolver->push(1);
        addDepthConstraints();
        success = checkModel();
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing depth", TimeUsage::FULL);
        if (success){
            extractModel();
            // getchar();
            _smt.pSolver->pop(1);
            for (i = 1; i < step && !find_min_depth && has_jump; ++i){
                fprintf(stdout, "[Info]          trying to optimize for depth bound %d            \n", _olsqParam.min_depth);
                _timer.start(TimeUsage::PARTIAL);
                --_olsqParam.min_depth;
                _smt.pSolver->push(1);
                addDepthConstraints();
                success = checkModel();
                fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
                _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
                _timer.showUsage("optimizing depth", TimeUsage::FULL);
                if (success) 
                    extractModel();
                else
                    find_min_depth = true;
                _smt.pSolver->pop(1);
            }
            find_min_depth = true;
        }
        else{
            _smt.pSolver->pop(1);
            if(_olsqParam.min_depth + step < _olsqParam.max_depth){
                updateSMT(step);
                if(step > 1)
                    has_jump = true;
            }
            _olsqParam.min_depth += step;
        }
    }
    if (find_min_depth){
        return true;
    }
    else{
        return false;
    }
}

bool OLSQ2::optimizeSwap(){
    bool success_optimize;
    success_optimize = optimizeDepth();
    // getchar();
    if (!success_optimize){
        return false;
    }
    if(_pCircuit->nSwapGate() == 0){
        return true;
    }

    unsigned_t lower_swap_bound = 0;
    unsigned_t upper_swap_bound = (_olsqParam.is_use_heuristic_bound_for_swap) ? _olsqParam.heuristic_swap_bound : _pCircuit->nGate();
    upper_swap_bound = (_pCircuit->nSwapGate() < upper_swap_bound) ? _pCircuit->nSwapGate() : upper_swap_bound;
    bool reduce_swap = true;
    bool firstRun = true;
    unsigned_t step = 2; 
    while (reduce_swap && (!_timer.isTimeout() || !_hasSol)){
        // cout << "enter loop" << endl;
        _smt.pSolver->push(1);
        addDepthConstraints();
        reduce_swap = optimizeSwapForDepth(lower_swap_bound, upper_swap_bound, firstRun);
        _smt.pSolver->pop(1);
        upper_swap_bound = _pCircuit->nSwapGate() - 1;
        firstRun = false;
        // getchar();
        if(reduce_swap){
            fprintf(stdout, "[Info] Successfully reduce SWAP count. Go to next run.            \n");
            fprintf(stdout, "[Info] Solving with depth %d            \n", _olsqParam.min_depth + step);
            fprintf(stdout, "[Info] Generating formulation                        \n");
            if(_olsqParam.min_depth + step < _olsqParam.max_depth){
                updateSMT(step);
                _olsqParam.min_depth += step;
            }
            else{
                _olsqParam.min_depth += step;
                increaseDepthBound();
                generateFormulation();
            }
        }
    }
    return true;
}

bool OLSQ2::optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun){
    unsigned_t swap_bound = upper_swap_bound;
    bool find_min_swap = false;
    bool success;
    while (!find_min_swap && lower_swap_bound <= swap_bound && swap_bound <= upper_swap_bound && (!_timer.isTimeout() || !_hasSol)){
        fprintf(stdout, "[Info]          trying to optimize for swap bound %d            \n", swap_bound);
        if(!_hasSol){
            _smt.resetTimeState(_olsqParam.timeout);
        }
        else{
            _smt.resetTimeState(_olsqParam.timeout_per_run);
        }
        _timer.start(TimeUsage::PARTIAL);
        _smt.pSolver->push(1);
        addSwapCountConstraints(swap_bound);

        success = checkModel();
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing swap", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing swap", TimeUsage::FULL);
        if (success){
            _hasSol = true;
            extractModel();
            if(swap_bound > _pCircuit->nSwapGate()){
                swap_bound = _pCircuit->nSwapGate();
            }
            --swap_bound;
        }
        else{
            if (swap_bound < upper_swap_bound){
                find_min_swap = true;
            }
            else if (!firstRun && swap_bound == upper_swap_bound){
                break;
            }
            else if (firstRun && swap_bound == upper_swap_bound){
                // cout << "line 518: increase swap bound" << endl;
                ++upper_swap_bound;
                lower_swap_bound = upper_swap_bound;
                swap_bound = upper_swap_bound;
            }
        }
        _smt.pSolver->pop(1);
    }
    return find_min_swap;
}

void OLSQ2::extractModel(){
    if (_verbose > 0){
        fprintf(stdout, "[Info] Extract Model Info                              \n");
    }
    unsigned_t circuitDepth = 0, i, gateTime, q, j, swapId, qId, t, e, tt;
    string s;

    // collect gate execution time
    vector<int_t> vQubitFirstGateTime(_pCircuit->nProgramQubit(), -1);
    vector<int_t> vQubitLastGateTime(_pCircuit->nProgramQubit(), -1);
    for ( i = 0; (i < _pCircuit->nGate());  ++i){
        Gate &gate = _pCircuit->gate(i);
        gateTime = stoi(_smt.pSolver->get_value( _smt.vTg[i]).value<string>(10)); 
        circuitDepth = (circuitDepth < gateTime) ? gateTime : circuitDepth;
        gate.setExecutionTime(gateTime);
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            if(vQubitFirstGateTime[gate.targetProgramQubit(j)] == -1){
                vQubitFirstGateTime[gate.targetProgramQubit(j)] = gateTime;
            }
            vQubitLastGateTime[gate.targetProgramQubit(j)] = (vQubitLastGateTime[gate.targetProgramQubit(j)] < (int_t)gateTime) ? (int_t)gateTime: vQubitLastGateTime[gate.targetProgramQubit(j)];
        }
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            gate.setTargetPhysicalQubit(j, stoi(_smt.pSolver->get_value(_smt.vvPi[gateTime][gate.targetProgramQubit(j)]).value<string>(10)));
        }
        if (_verbose > 0){
            fprintf(stdout, "        - Gate %d, name: %s, duration: %d, time: %d, target program qubit: ", i, gate.name().c_str(), gate.duration(), gate.executionTime());
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                fprintf(stdout, "%d ", gate.targetProgramQubit(j));
            }
            fprintf(stdout, ", target physical qubit: ");
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                fprintf(stdout, "%d ", gate.targetPhysicalQubit(j));
            }
            fprintf(stdout, "\n");
        }
    }

    // collect qubit mapping
    vector<vector<int_t> > vvPhy2Pro(circuitDepth,  vector<int_t>(_device.nQubit(), -1));
    for (t = 0; t < circuitDepth; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            vvPhy2Pro[t][stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10))] = i;
        }
    }

    // // get SWAP gate
    _pCircuit->clearSwap();
    vector<unsigned_t> swapTargetQubit(2,0);
    bool cond1, cond2;
    int_t tmp, bound = circuitDepth - _olsqParam.swap_duration;
    vector<vector<bool> > vvTimeSwap(circuitDepth, vector<bool>(_device.nEdge(), 0));
    // cerr << "get swap value"<< endl;
    for (t = _olsqParam.swap_duration -1; t < circuitDepth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            if(_smt.pSolver->get_value( _smt.vvSigma[t][e]).is_true()){
                vvTimeSwap[t][e] = true;    
            }
        }
    }
    // cerr << endl;
    if(!_olsqParam.is_transition){
        for (e = 0; e < _device.nEdge(); ++e){
            for (t = _olsqParam.swap_duration -1; t < bound; ++t){
                if(vvTimeSwap[t][e] && vvTimeSwap[t + _olsqParam.swap_duration][e]){
                    // only cancel two consecutive swap
                    vvTimeSwap[t][e] = 0;
                    vvTimeSwap[t + _olsqParam.swap_duration][e] = 0;
                }
            }
        }
    }
    for (t = _olsqParam.swap_duration -1; t < circuitDepth; ++t){
        for (e = 0; e < _device.nEdge(); ++e){
            if (vvTimeSwap[t][e]){
                cond1 = !((vvPhy2Pro[t][_device.edge(e).qubitId1()] == -1) && (vvPhy2Pro[t][_device.edge(e).qubitId2()] == -1));
                cond2 = cond1;
                // check if it is before the first gates of these qubits
                if(vvPhy2Pro[t][_device.edge(e).qubitId1()] != -1){
                    cond1 = (vQubitFirstGateTime[vvPhy2Pro[t][_device.edge(e).qubitId1()]] <= t);
                }
                if(vvPhy2Pro[t][_device.edge(e).qubitId2()] != -1){
                    cond1 = (vQubitFirstGateTime[vvPhy2Pro[t][_device.edge(e).qubitId2()]] <= t) || cond1;
                }
                // check if it is after the last gates of these qubits
                if(vvPhy2Pro[t][_device.edge(e).qubitId1()] != -1){
                    cond2 = (vQubitLastGateTime[vvPhy2Pro[t][_device.edge(e).qubitId1()]] > t);
                }
                if(vvPhy2Pro[t][_device.edge(e).qubitId2()] != -1){
                    cond2 = (vQubitLastGateTime[vvPhy2Pro[t][_device.edge(e).qubitId2()]] > t) || cond2;
                }
                // cerr << "cond1: " << cond1 << ", cond2: " << cond2 << endl;
                if(cond1 && cond2){
                    swapTargetQubit[0] = _device.edge(e).qubitId1();
                    swapTargetQubit[1] = _device.edge(e).qubitId2();
                    swapId = _pCircuit->nSwapGate();
                    _pCircuit->addSwapGate(swapId, swapTargetQubit, _olsqParam.swap_duration);
                    Gate & gate = _pCircuit->swapGate(swapId);
                    gate.setExecutionTime(t);
                    if (_verbose > 0){
                        fprintf(stdout, "        - SWAP Gate %d, duration: %d, time: %d, target qubit: %d %d\n", swapId + _pCircuit->nGate(), gate.duration(), gate.executionTime(), gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1));
                    }
                    if(vvPhy2Pro[t][swapTargetQubit[0]] != -1){
                        vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] = (vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] > t) ? t: vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[0]]];
                        vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] = (vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]] <= t) ? t + 1: vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[0]]];
                    }
                    if(vvPhy2Pro[t][swapTargetQubit[1]] != -1){
                        vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] = (vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] > t) ? t: vQubitFirstGateTime[vvPhy2Pro[t][swapTargetQubit[1]]];
                        vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] = (vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]] <= t) ? t + 1: vQubitLastGateTime[vvPhy2Pro[t][swapTargetQubit[1]]];
                    }
                }
                // else{
                //     fprintf(stdout, "        - SWAP Gate, time: %d, target qubit: %d %d\n", t, _device.edge(e).qubitId1(), _device.edge(e).qubitId2());
                // }
            }
        }
    }
    // collect qubit mapping
    if (_verbose > 0){
        fprintf(stdout, "        - Qubit mapping: \n");
        for (t = 0; t <= circuitDepth; ++t){
        fprintf(stdout, "        - Time %d: ",t);
            for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
                fprintf(stdout, "%d->%d ", i, stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10)));
            }
        fprintf(stdout, "\n");
        }
    }
    for(i = 0; i < vQubitFirstGateTime.size(); ++i){
        if(vQubitFirstGateTime[i] == -1){
            vQubitFirstGateTime[i] = 0;
        }
        if(vQubitLastGateTime[i] == -1){
            vQubitLastGateTime[i] = 0;
        }
    }
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _pCircuit->setInitialMapping(i,  stoi(_smt.pSolver->get_value( _smt.vvPi[vQubitFirstGateTime[i]][i]).value<string>(10)));
        _pCircuit->setFinalMapping(i,  stoi(_smt.pSolver->get_value( _smt.vvPi[vQubitFirstGateTime[i]][i]).value<string>(10)));
    }

    if(_olsqParam.is_collect_qubit_region){
        _pCircuit->resetQubitRegion();
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            _pCircuit->addQubitRegion(i, _pCircuit->initialMapping(i));
            _pCircuit->addQubitRegion(i, _pCircuit->finalMapping(i));
            for (t = vQubitFirstGateTime[i]+1; t < vQubitLastGateTime[i]; ++t){
                _pCircuit->addQubitRegion(i,  stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10)));
            }
        }
    }
    _pCircuit->setCircuitDepth(circuitDepth+1);
}

void OLSQ2::asapScheduling(){
    vector<int_t> vPushForwardDepth(_device.nQubit(), -1);
    int_t gateExecutionTime;
    unsigned_t block, i, j, q0, q1, qId, maxTime = 0;
    Gate gate;
    unordered_set<unsigned_t> sGateId;
    for (block = 0; block < _pCircuit->circuitDepth(); ++block){
        for (i = 0; i < _pCircuit->nGate(); ++i){
            Gate & gate = _pCircuit->gate(i);
            if (gate.executionTime() == block && sGateId.count(gate.idx()) == 0){
                gateExecutionTime = vPushForwardDepth[gate.targetPhysicalQubit(0)];
                for (j = 1; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    gateExecutionTime = (gateExecutionTime < vPushForwardDepth[qId]) ? vPushForwardDepth[qId] : gateExecutionTime;
                }
                ++gateExecutionTime;
                for (j = 0; j < gate.nTargetQubit(); ++j){
                    qId = gate.targetPhysicalQubit(j);
                    vPushForwardDepth[qId] = gateExecutionTime;
                }
                maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                gate.setExecutionTime(gateExecutionTime);
                sGateId.insert(gate.idx());
            }
        }
        if (block < _pCircuit->circuitDepth() - 1){
            for (j = 0; j < _pCircuit->nSwapGate(); ++j){
                Gate & gate = _pCircuit->swapGate(j);
                if (gate.executionTime() == block && sGateId.count(gate.idx()+ _pCircuit->nGate()) == 0){
                    q0 = gate.targetPhysicalQubit(0);
                    q1 = gate.targetPhysicalQubit(1);
                    gateExecutionTime = (vPushForwardDepth[q0] < vPushForwardDepth[q1]) ? vPushForwardDepth[q1] : vPushForwardDepth[q0];
                    ++gateExecutionTime;
                    vPushForwardDepth[q0] = gateExecutionTime;
                    vPushForwardDepth[q1] = gateExecutionTime;
                    maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                    gate.setExecutionTime(gateExecutionTime);
                    sGateId.insert(gate.idx()+ _pCircuit->nGate());
                }
            }
        }
    }
    ++maxTime;
    _pCircuit->setCircuitDepth(maxTime);
}



void OLSQ2::increaseDepthBound(){
    while(_olsqParam.min_depth >= _olsqParam.max_depth){
        _olsqParam.max_depth_bit += _olsqParam.max_depth_expand_factor; 
        _olsqParam.max_depth = _olsqParam.max_depth << _olsqParam.max_depth_expand_factor;
    }
}

unsigned_t OLSQ2::extract_longest_chain(){
    vector<unsigned_t> vPushForwardDepth(_pCircuit->nProgramQubit(),0);
    Gate gate;
    unsigned_t max, qId, i, j;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        if (gate.nTargetQubit() == 1){
            ++vPushForwardDepth[gate.targetProgramQubit(0)];
        }
        else{
            max = 0;
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                qId = gate.targetProgramQubit(j);
                max = (max > vPushForwardDepth[qId]) ? max : vPushForwardDepth[qId];
            }
            ++max;
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                vPushForwardDepth[gate.targetProgramQubit(j)] = max;
            }
        }
    }
    return *std::max_element(vPushForwardDepth.begin(), vPushForwardDepth.end());
}

void OLSQ2::constructDependency(){
    _vpGateDependency.clear();
    vector<int_t> vQubitLastGate(_pCircuit->nProgramQubit(), -1);
    unsigned_t qId, i, j;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            qId = gate.targetProgramQubit(j);
            if (vQubitLastGate[qId] > -1){
                addDependency(vQubitLastGate[qId], i);
            }
            vQubitLastGate[qId] = i;
        }
    }
}

void OLSQ2::constructCollisionList(){
    _vGateCollision.clear();
    unsigned_t qId, i, j;
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        for (j = 0; j < _pCircuit->nGate(); ++j){
            Gate & gate1 = _pCircuit->gate(j);
            if(gate.targetProgramQubit(0) == gate1.targetProgramQubit(0)){
                _vpGateDependency.emplace_back(i, j);
            }
            else if(gate.nTargetQubit() == 2){
                if(gate.targetProgramQubit(1) == gate1.targetProgramQubit(0)){
                    _vpGateDependency.emplace_back(i, j);
                }
                else if(gate1.nTargetQubit() == 2){
                    if((gate.targetProgramQubit(1) == gate1.targetProgramQubit(0))||(gate.targetProgramQubit(0) == gate1.targetProgramQubit(1))){
                        _vpGateDependency.emplace_back(i, j);
                    }
                }
            }
        }
    }
}

void OLSQ2::printDependency(){
    pair<unsigned_t, unsigned_t> p;
    unsigned_t i;
    fprintf(stdout, "[Info] Circuit Dependency Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Dependency list\n");
    for ( i = 0; i < _vpGateDependency.size();  ++i ){
        fprintf(stdout, "        - Gate %d depends on gate %d\n", _pCircuit->gate(_vpGateDependency[i].second).idx(), _pCircuit->gate(_vpGateDependency[i].first).idx());
    }
}

void OLSQ2::updateSMT(unsigned_t d){
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints(d);
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addValidTwoQubitGateConstraints(d); 
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraints(d);
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraints(d);
    if(_olsqParam.use_window_range_for_gate){
        updateGateTimeWindow(d);
    }
}

void OLSQ2::constructGateTimeWindow(){
    vector<int_t> vQubitArriveTime(_pCircuit->nProgramQubit(), 0);
    vector<int_t> vQubitRequiredTime(_pCircuit->nProgramQubit(), _olsqParam.min_depth-1);
    unsigned_t qId, j, t;
    int_t i;
    for ( i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        t = 0;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            qId = gate.targetProgramQubit(j);
            t = (vQubitArriveTime[qId] < t) ? t : vQubitArriveTime[qId];
        }
        _vpGateTimeWindow.emplace_back(make_pair(t, 0));
        ++t;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            vQubitArriveTime[gate.targetProgramQubit(j)] = t;
        }
        // cerr << "qubit arrive time: " << endl;
        // for (auto qat : vQubitArriveTime){
        //     cerr << qat << " ";
        // }
        // cerr << endl;
    }
    // cerr << "_olsqParam.min_depth: " << _olsqParam.min_depth << endl;
    for ( i = _pCircuit->nGate() - 1; i >= 0; --i){
        Gate & gate = _pCircuit->gate(i);
        t = _olsqParam.min_depth - 1;
        // cerr << "before t: " << t << endl;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            qId = gate.targetProgramQubit(j);
            t = (vQubitRequiredTime[qId] > t) ? t : vQubitRequiredTime[qId];
            // cerr << "update t: " << t << endl;
            
        }
        // cerr << "t: " << t << endl;
        _vpGateTimeWindow[i].second = t;
        --t;
        for ( j = 0; j < gate.nTargetQubit();  ++j ){
            vQubitRequiredTime[gate.targetProgramQubit(j)] = t;
        }
        
    }
    printGateTimeWindow();
}

void OLSQ2::updateGateTimeWindow(unsigned_t d){
    for(int_t i = 0; i < _vpGateTimeWindow.size(); ++i){
        _vpGateTimeWindow[i].second += d;
    }
    // printGateTimeWindow();
}

void OLSQ2::printGateTimeWindow(){
    fprintf(stdout, "[Info] Gate Time Window Info                              \n");
    fprintf(stdout, "       ------------------------------------------\n");
    fprintf(stdout, "       - Min depth: %d                           \n", _olsqParam.min_depth);
    fprintf(stdout, "       - Gate time window:\n");
    for (unsigned_t i = 0; i < _vpGateTimeWindow.size();  ++i ){
        fprintf(stdout, "        - Gate %d: [%d, %d]\n", i, _vpGateTimeWindow[i].first, _vpGateTimeWindow[i].second);
    }
}


MOLSQ_NAMESPACE_CPP_END
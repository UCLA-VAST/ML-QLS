/***********************************************************************
  File        [ olsq.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ sTBOLSQ2 ]
  Synopsis    [ sTBOLSQ2 class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 24, Feb., 2022 ]
***********************************************************************/

#include "stbolsq2/stbolsq2.hpp"

MOLSQ_NAMESPACE_CPP_START


bool sTBOLSQ2::run(unsigned_t bound, bool considerDis){
    assert(_pCircuit->nProgramQubit() <= _pDevice->nQubit());
    _timer.start(TimeUsage::FULL);
    // if(_stbolsqParam.is_refinement){
    //     _pvpGateDependency = _pCircuit->pvpGateDependency();
    // }
    _hasSol = false;
    _bestDisSum = MAX_DOUBLE;
    _stbolsqParam.considerDis = considerDis;
    _timer.setTimeout(_stbolsqParam.timeout);
    if(_verbose == 2){
        _pCircuit->printDependency();
    }
    // cout << "ENTER RUNSMT" << endl;

    runSMT(bound);
    // cout << "finish scheduling" << endl;
    // _pCircuit->printCircuitLayout();
    _timer.showUsage("sTBOLSQ2: Total runtime", TimeUsage::FULL);
    return _hasSol;
}

void sTBOLSQ2::initQubitRegion(){
    _vsQubitRegion.clear();
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vsQubitRegion.emplace_back(_pCircuit->sQubitRegion(i));
    }
}

void sTBOLSQ2::runSMT(unsigned_t bound){
    fprintf(stdout, "[Info] sTBOLSQ2                        \n");
    _iter = 0;
    increaseDepthBound();
    fprintf(stdout, "[Info]          relabeling qubit index                       \n");
    initQubitRegion();
    computeQubitIdx();
    computeOverlapRegion();
    bool solve = false;
    while (_iter < 2 && !solve && !_timer.isTimeout()){
        fprintf(stdout, "[Info] Iter %d: Solving with depth %d            \n", _iter, _stbolsqParam.min_depth);
        fprintf(stdout, "[Info] Iter %d: Generating formulation                        \n", _iter);
        generateFormulation();
        fprintf(stdout, "[Info] Iter %d \n", _iter);
        fprintf(stdout, "[Info]         Optimizing change pi cnt                             \n");
        solve = optimize(bound);
        ++_iter;
        // fprintf(stdout, "solve: %d, _iter: %d\n", solve, _iter);
        if(!solve || _iter < 2){
            if(!solve && _iter == 1){
                _iter = 0;
                increaseDepthBound();
            }
            else{
                bound = _changeCntBound;
                if(_pCircuit->nGate() > 40){
                    break;
                }
            }
            expandQubitRegion();
            // if(_iter > 2){
            //     expandQubitRegion();
            // }
        }
    }
}

void sTBOLSQ2::dumpSMT(){
    _stbolsqParam.max_depth = 5;
    generateFormulation();
    Sort sortbvtime = mk_bv_sort( ceil(log(_stbolsqParam.max_depth+1)));
    Term depthBv = mk_bv_value_uint64( sortbvtime, _stbolsqParam.max_depth);
    for (Term & t: _smt.vTg){
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {t, depthBv}));
    }
    fprintf(stdout, "[Info]          constructing change pi constraint           \n");
    addChangePiConstraint();
    fprintf(stdout, "[Info]          constructing change pi count constraint           \n");
    addChangeCountConstraints(30);
    string filename = to_string(_pDevice->nQubit())+"_"+to_string(_pCircuit->nProgramQubit())+"_depth_"+to_string(_stbolsqParam.max_depth)+".txt";
    std::ofstream fout(filename);
    _smt.pSolver->print_formula(fout);
}

void sTBOLSQ2::generateFormulation(){
    _smt.reset(_stbolsqParam.timeout_per_run, _iter);
    fprintf(stdout, "[Info]          constructing variables                       \n");
    constructVariable();
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints();
    cerr << "_stbolsqParam.fixInitialMapping: " << _stbolsqParam.fixInitialMapping << endl;
    if(_stbolsqParam.fixInitialMapping){
        fprintf(stdout, "[Info]          constructing initial mapping constraint    \n");
        addInitialMappingConstraints();
    }
    fprintf(stdout, "[Info]          constructing mapping variable region    \n");
    addMappingRegionConstraints();
    fprintf(stdout, "[Info]          constructing dependency constraint           \n");
    addDependencyConstraints();
    fprintf(stdout, "[Info]          constructing qubit moving constraint \n");
    addQubitMovingDistanceConstraints();
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addRestrictedValidTwoQubitGateConstraints();
    fprintf(stdout, "[Info]          constructing change pi constraint           \n");
    addChangePiConstraint();
}

void sTBOLSQ2::computeQubitIdx(){
    _vvQubitRegion.clear();
    _vvQubitRegion.resize(_pCircuit->nProgramQubit());
    unsigned_t idx;
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        // fprintf(stdout, "[Info]          Relabel region for qubit %d: ", i);
        idx = 0;
        for(unsigned_t j : _vsQubitRegion[i]){
            // fprintf(stdout, " %d->%d ", idx, j);
            _vvQubitRegion[i].emplace_back(j);
            ++idx;
        }
        // fprintf(stdout, "\n");
    }
}

void sTBOLSQ2::computeOverlapRegion(){
    _vmPhy2ProIdxChoiceIdx.clear();
    _vmPhy2ProIdxChoiceIdx.resize(_pDevice->nQubit());
    for(unsigned_t i = 0; i < _vvQubitRegion.size(); ++i){
        for(unsigned_t j = 0; j < _vvQubitRegion[i].size(); ++j){
            _vmPhy2ProIdxChoiceIdx[_vvQubitRegion[i][j]][i] = j;
        }    
    }
}

void sTBOLSQ2::expandQubitRegion(){
    unsigned_t i, j;
    for(i = 0; i < _vvQubitRegion.size(); ++i){
        if(_vvQubitRegion[i].size() < _pDevice->nQubit()){
            for(j = 0; j < _vvQubitRegion[i].size(); ++j){
                Qubit & qubit = _pDevice->qubit(_vvQubitRegion[i][j]);
                for(unsigned_t e : qubit.vSpanEdge){
                    Edge & edge = _pDevice->edge(e);
                    if(_vsQubitRegion[i].count(edge.qubitId1()) == 0){
                        _vsQubitRegion[i].insert(edge.qubitId1());
                        _vmPhy2ProIdxChoiceIdx[edge.qubitId1()][i] = _vvQubitRegion[i].size();
                        _vvQubitRegion[i].emplace_back(edge.qubitId1());
                    }
                    else if(_vsQubitRegion[i].count(edge.qubitId2()) == 0){
                        _vsQubitRegion[i].insert(edge.qubitId2());
                        _vmPhy2ProIdxChoiceIdx[edge.qubitId2()][i] = _vvQubitRegion[i].size();
                        _vvQubitRegion[i].emplace_back(edge.qubitId2());
                    }
                }
            }    
        }
    }
}


void sTBOLSQ2::constructVariable(){
    unsigned_t bit_length_pi, bit_length_time;
    bit_length_time = _stbolsqParam.max_depth_bit;
    // fprintf(stdout, "[Info]          bit_length_pi: %d, bit_length_time: %d\n", bit_length_pi, bit_length_time);
    // _smt.vvPi.resize(_stbolsqParam.max_depth, vector<expr> (_pCircuit->nProgramQubit()));
    _smt.vvPi.reserve(_stbolsqParam.max_depth);
    _smt.vTg.reserve(_pCircuit->nGate());
    _smt.vvChangeCnt.reserve(_stbolsqParam.max_depth);

    string s;
    unsigned_t i, j;

    // Create a bit-vector sort of size 1.
    Sort sortbvtime = mk_bv_sort(bit_length_time);
    Sort sortbool = mk_bool_sort();

    for (i = 0; i < _vvQubitRegion.size(); ++i){
        bit_length_pi = ceil(log2(_vvQubitRegion[i].size() + 1));
        _smt.vPiSort.emplace_back(mk_bv_sort( bit_length_pi));
    }

    for (i = 0; i < _stbolsqParam.max_depth; ++i){
        _smt.vvPi.emplace_back(vector<Term> ());
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            s = "map_t" + to_string(i) + "_q" + to_string(j);
            // _smt.vvPi[i].emplace_back(_smt.c.bv_const(s.c_str(), bit_length_pi));
            _smt.vvPi[i].emplace_back(mk_const(_smt.vPiSort[j], s.c_str()));
        }
    }

    for (i = 0; i < _pCircuit->nGate(); ++i){
        s = "time_" + to_string(i);
        // _smt.vTg.emplace_back(_smt.c.bv_const(s.c_str(), bit_length_time));
        _smt.vTg.emplace_back(mk_const(sortbvtime, s.c_str()));
    }

    for (i = 0; i < _stbolsqParam.max_depth - 1; ++i){
        _smt.vvChangeCnt.emplace_back(vector<Term> ());
        for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
            s = "ifchange_t" + to_string(i) + "_q" + to_string(j);
            // _smt.vvSigma[i].emplace_back(_smt.c.bool_const(s.c_str()));
            _smt.vvChangeCnt[i].emplace_back(mk_const(sortbool, s.c_str()));
        }
    }
}


void sTBOLSQ2::addMappingRegionConstraints(unsigned_t boundOffset){
    unsigned_t end = _stbolsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    // add mapping region constraint
    Term qBv, zero;
    for (unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        qBv = mk_bv_value_uint64( _smt.vPiSort[i], _vvQubitRegion[i].size());
        zero = mk_bv_zero(_smt.vPiSort[i]);
        for(unsigned_t t = begin; t < end; ++t){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vvPi[t][i]}));
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {_smt.vvPi[t][i], qBv}));
        }
    }
}

void sTBOLSQ2::addInitialMappingConstraints(){
    Term qBv;
    // cerr << "_pvInitialMapping->size(): " << _pvInitialMapping->size() << endl;
    for (unsigned_t i = 0; i < _pvInitialMapping->size(); ++i){
        // cerr << "at time 0, qubit " << i << " is mapped to " << (*_pvInitialMapping)[i] << "(" << _vmPhy2ProIdxChoiceIdx[(*_pvInitialMapping)[i]][i] << ")" << endl;
        qBv = mk_bv_value_uint64( _smt.vPiSort[i], _vmPhy2ProIdxChoiceIdx[(*_pvInitialMapping)[i]][i]);
         _smt.pSolver->assert_formula(mk_term(Kind::EQUAL, {_smt.vvPi[0][i], qBv}));
    }
}

void sTBOLSQ2::addInjectiveMappingConstraints(unsigned_t boundOffset){
    Term cond1, cond2, q1Bv, q2Bv;
    unsigned_t end = _stbolsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (unsigned_t i = 0; i < _vmPhy2ProIdxChoiceIdx.size(); ++i){
        for (auto iterj = _vmPhy2ProIdxChoiceIdx[i].begin(); iterj != _vmPhy2ProIdxChoiceIdx[i].end(); iterj++){
            q1Bv = mk_bv_value_uint64( _smt.vPiSort[iterj -> first], iterj -> second);
            for(auto iterk = next(iterj, 1); iterk != _vmPhy2ProIdxChoiceIdx[i].end(); ++iterk){
                q2Bv = mk_bv_value_uint64( _smt.vPiSort[iterk->first], iterk->second);
                for(unsigned_t t = begin; t < end; ++t){
                    cond1 = mk_term(Kind::EQUAL, {_smt.vvPi[t][iterj->first], q1Bv});
                    cond2 = mk_term(Kind::EQUAL, {_smt.vvPi[t][iterk->first], q2Bv});
                     _smt.pSolver->assert_formula(mk_term(Kind::NOT, 
                                {mk_term(Kind::AND, {cond1, cond2})}));
                }
            }
        }
    }
}


void sTBOLSQ2::addRestrictedValidTwoQubitGateConstraints(unsigned_t boundOffset){
    // cerr << "In addRestrictedValidTwoQubitGateConstraintsZ3" << endl;
    unsigned_t i, t, j, k, q0, q1;
    unsigned_t end = _stbolsqParam.min_depth, begin = 0;
    Sort sortbvtime = mk_bv_sort( _stbolsqParam.max_depth_bit);
    Term q1Bv, q2Bv, pro1EqPhy, pro2EqPhy;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for ( i = 0; i < _pCircuit->nGate();  ++i ){
        Gate & gate = _pCircuit->gate(i);
        if ((gate).nTargetQubit() == 2){
            for (t = begin; t < end; ++t){
                Term bvt = mk_bv_value_uint64( sortbvtime, t);
                Term clause = mk_term(Kind::NOT, {mk_term(Kind::EQUAL, {_smt.vTg[i], bvt})});
                q0 = gate.targetProgramQubit(0);
                q1 = gate.targetProgramQubit(1);
                for (j = 0; j < _vvQubitRegion[q0].size(); ++j ){
                    for (k = 0; k < _vvQubitRegion[q1].size(); ++k ){
                        if(_pDevice->isAdjacent(_vvQubitRegion[q0][j], _vvQubitRegion[q1][k])){
                            q1Bv = mk_bv_value_uint64( _smt.vPiSort[q0], j);
                            q2Bv = mk_bv_value_uint64( _smt.vPiSort[q1], k);
                            pro1EqPhy = mk_term(Kind::EQUAL, {_smt.vvPi[t][q0], q1Bv});
                            pro2EqPhy = mk_term(Kind::EQUAL, {_smt.vvPi[t][q1], q2Bv});
                            clause = mk_term(Kind::OR, {clause, 
                                        mk_term(Kind::AND, {pro1EqPhy, pro2EqPhy})});
                        }
                    }    
                }
                 _smt.pSolver->assert_formula(clause);
            }
        }
    }   
}



void sTBOLSQ2::addDependencyConstraints(){
    unsigned_t i;
    Gate g;
    Sort sortbvtime = mk_bv_sort( _stbolsqParam.max_depth_bit);
    Term zero = mk_bv_zero(sortbvtime);
    for (i = 0; i < _pCircuit->nGate(); ++i){
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vTg[i]}));
    }
    if(!_stbolsqParam.is_all_commute){
        for ( i = 0; i < _pvpGateDependency->size();  ++i ){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {_smt.vTg[(*_pvpGateDependency)[i].first], _smt.vTg[(*_pvpGateDependency)[i].second]}));
        }
    }
}

void sTBOLSQ2::addDepthConstraints(){
    // cout << " add depth constraints" << endl;
    Sort sortbvtime = mk_bv_sort( _stbolsqParam.max_depth_bit);
    Term depthBv = mk_bv_value_uint64( sortbvtime, _stbolsqParam.min_depth);
    for (Term & t: _smt.vTg){
        // _smt.smtSolver.add(ult(t, _stbolsqParam.min_depth));
         _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {t, depthBv}));
    }
    // cout << "finish add depth constraints" << endl;
}

void sTBOLSQ2::addQubitMovingDistanceConstraints(unsigned_t boundOffset){
    // we only use distance-1 location. may extend to other distance in the future
    unsigned_t t, q, i, j;
    unsigned_t end = _stbolsqParam.min_depth - 1, begin = 0;
    Term qBv, clause, cond;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (t = begin; t < end; ++t){ 
        for (q = 0; q < _pCircuit->nProgramQubit(); ++q){
            for (i = 0; i < _vvQubitRegion[q].size(); ++i){
                qBv = mk_bv_value_uint64( _smt.vPiSort[q], i);
                clause = mk_term(Kind::EQUAL, {_smt.vvPi[t+1][q], qBv});
                cond = mk_term(Kind::NOT, {mk_term(Kind::EQUAL, {_smt.vvPi[t][q], qBv})});
                for (j = 0; j < _vvQubitRegion[q].size(); ++j){
                    if(_pDevice->isAdjacent(_vvQubitRegion[q][i], _vvQubitRegion[q][j])){
                        qBv = mk_bv_value_uint64( _smt.vPiSort[q], j);
                        clause = mk_term(Kind::OR, {clause,
                                     mk_term(Kind::EQUAL, {_smt.vvPi[t+1][q], qBv})});
                    }
                }
                 _smt.pSolver->assert_formula(mk_term(Kind::OR, {cond, clause}));
            }
        }
    }
}

void sTBOLSQ2::addChangePiConstraint(unsigned_t boundOffset){
    unsigned_t t, q;
    unsigned_t end = _stbolsqParam.min_depth - 1, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (t = begin; t < end; ++t){ 
        for (q = 0; q < _pCircuit->nProgramQubit(); ++q){
            // _smt.smtSolver.add(_smt.vvChangeCnt[t][q] == (_smt.vvPi[t][q] != _smt.vvPi[t+1][q]));
             _smt.pSolver->assert_formula(
                mk_term(Kind::IFF, {_smt.vvChangeCnt[t][q],
                    mk_term(Kind::NOT, 
                    {mk_term(Kind::EQUAL, {_smt.vvPi[t][q], _smt.vvPi[t+1][q]})})}));
        }
    }
}

void sTBOLSQ2::addChangeCountConstraints(unsigned_t bound){
    PB2CNF pb2cnf;
    vector<int_t> vChangeLit;
    unsigned_t firstFreshVariable = 1 + (_stbolsqParam.min_depth-1)*_pCircuit->nProgramQubit(), newFreashVariable;
    for (unsigned_t i = 1; i < firstFreshVariable; ++i){
        vChangeLit.emplace_back(i);
    }
    vector< vector<int_t> > formula;
    newFreashVariable = pb2cnf.encodeAtMostK(vChangeLit, bound, formula, firstFreshVariable) + 1;
    map<unsigned_t, Term> mAncillary;
    vector<Term> vOrs;
    unsigned_t sigmaT, sigmaE, var;
    Sort sortbool = mk_bool_sort();
    string s;
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                sigmaT = (var - 1) / _pCircuit->nProgramQubit();
                sigmaE = (var - 1) % _pCircuit->nProgramQubit();
                if (lit < 0){ 
                    vOrs.emplace_back(mk_term(Kind::NOT, {_smt.vvChangeCnt[sigmaT][sigmaE]}));
                }
                else{
                    vOrs.emplace_back(_smt.vvChangeCnt[sigmaT][sigmaE]);
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



bool sTBOLSQ2::checkModel(){
    Result status = _smt.pSolver->check_sat();
    if (status == Result::SAT){
        return true;
    }
    else{
        return false;
    }
}

bool sTBOLSQ2::optimize(int_t change_bound){
    bool success = 1;
    if(_iter == 0){
        success = optimizeDepth();
        // getchar();
        if (!success){
            return false;
        }
    }

    unsigned_t step = 1;
    // optimize change cnt
    _smt.pSolver->push(1);
    addDepthConstraints();
    optimizeChangeCnt(change_bound);
    _smt.pSolver->pop(1);
    // increase depth and sample solution again
    // while (!_timer.isTimeout() && success && (int)_bestDisSum > 0){
    //     fprintf(stdout, "[Info] Successfully reduce Change count. Go to next run.            \n");
    //     fprintf(stdout, "[Info] Solving with depth %d            \n", _stbolsqParam.min_depth + step);
    //     fprintf(stdout, "[Info] Generating formulation                        \n");
    //     if(_stbolsqParam.min_depth + step < _stbolsqParam.max_depth){
    //         updateSMT(step);
    //         _stbolsqParam.min_depth += step;
    //     }
    //     else{
    //         _stbolsqParam.min_depth += step;
    //         increaseDepthBound();
    //         generateFormulation();
    //     }
    //     bitwuzla->push(1);
    //     addDepthConstraints();
    //     success = sampleSolution();
    //     pSolver->pop(1);
    // }
    return true;
}


bool sTBOLSQ2::optimizeDepth(){
    bool success, find_min_depth = false;
    unsigned_t i;
    bool tmp;
    while(!find_min_depth && _stbolsqParam.min_depth < _stbolsqParam.max_depth ){
        fprintf(stdout, "[Info]          Trying to optimize for depth bound %d            \n", _stbolsqParam.min_depth);
        _timer.start(TimeUsage::PARTIAL);
        if(!_hasSol){
            _smt.resetTimeState(_stbolsqParam.timeout_per_run);
        }
        _smt.pSolver->push(1);
        addDepthConstraints();
        success = checkModel();
        fprintf(stdout, "[Info]          Optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing depth", TimeUsage::FULL);
        if (success){
            tmp = extractModel();
            _smt.pSolver->pop(1);
            find_min_depth = true;
        }
        else{
            _smt.pSolver->pop(1);
            if(_stbolsqParam.min_depth + 1 < _stbolsqParam.max_depth){
                updateSMT(1);
            }
            ++_stbolsqParam.min_depth;
            _stbolsqParam.timeout_per_run *= 2;
        }
    }
    return find_min_depth;
}

void sTBOLSQ2::optimizeChangeCnt(int_t change_bound){
    // optimize the pi count change
    unsigned_t step = 1;
    unsigned_t init_change_cnt = change_bound, reduce_change_cnt = 0, upper_change_bound = change_bound;
    bool isIncrease = false, success, reduceDis = false;
    if (change_bound > 19){
        step = change_bound * 0.1;
    }
    _stbolsqParam.timeout_per_run = _stbolsqParam.timeout_per_run_base;
    unsigned_t changeBoundLimit = min(_pCircuit->nProgramQubit() * _stbolsqParam.min_depth, _changeCntBound);
    while ((!_timer.isTimeout() || !_hasSol) && change_bound >= 0 && (int)_bestDisSum > 0 && change_bound < changeBoundLimit){
        // cerr << "_hasSol: " << _hasSol << ", change_bound: " << change_bound << endl;
        if(!_hasSol){
            _smt.resetTimeState(_stbolsqParam.timeout_per_run_base);
        }
        fprintf(stdout, "[Info]          Trying to optimize for change bound %d            \n", change_bound);
        _timer.start(TimeUsage::PARTIAL);
        _smt.pSolver->push(1);
        addChangeCountConstraints(change_bound);
        success = checkModel();
        fprintf(stdout, "[Info]          Optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("sTBOLSQ2: optimizing change", TimeUsage::PARTIAL);
        _timer.showUsage("sTBOLSQ2: optimizing change", TimeUsage::FULL);
        if (success){
            _hasSol = true;
            reduceDis = extractModel();
            _smt.pSolver->pop(1);
            if(isIncrease){
                if(_stbolsqParam.considerDis){
                    reduceDis = sampleSolution();
                }
                return;
            }
            if((int)_bestDisSum * 2 == _changeCntBound){
                _smt.resetTimeState(_stbolsqParam.timeout_short);
            }
            change_bound -= step;
        }
        else{
            _smt.pSolver->pop(1);
            if (change_bound < init_change_cnt){
                // invoke solver several times to get better dis
                if(_stbolsqParam.considerDis){
                    reduceDis = sampleSolution();
                }
                return;
            }
            else if (change_bound == upper_change_bound){
                isIncrease = true;
                if(step > 5){
                    upper_change_bound = _changeCntBound;
                    change_bound = _changeCntBound;
                }
                else{
                    upper_change_bound += step;
                    change_bound += step;
                }
                step *= 3;
                _stbolsqParam.timeout_per_run += _stbolsqParam.timeout_per_run_increase_step;
            }
        }
    }
    return;
}

bool sTBOLSQ2::sampleSolution(unsigned_t solCntBound){
    unsigned_t change_bound = (unsigned_t)_bestDisSum;
    if(_changeCntBound == 0){
        change_bound = 4; // allow two qubit to swap?
    }
    else if(change_bound < 30){
        change_bound *= 2;
    }
    fprintf(stdout, "[Info] Start sample solutions with change bound %d                         \n", change_bound);
    _smt.pSolver->push(1);
    unsigned_t solCnt = 0;
    bool success = true, reduceDis, tmp;
    // if(change_bound < 100){
    //     change_bound *= 1.1; // may be more qubit change can lead to better dis result
    // }
    // else{
    //     change_bound += 5;
    // }
    addChangeCountConstraints(change_bound);
    // while(success && solCnt < solCntBound && !_timer.isTimeout()){
    while(success && solCnt < solCntBound){
        _smt.resetTimeState(_stbolsqParam.timeout_short);
        success = checkModel();
        if(success){
            tmp = extractModel(1);
        }
        if(tmp){
            reduceDis = 1;
        }
        ++solCnt;
    }
    _smt.pSolver->pop(1);
    return reduceDis;
}

void sTBOLSQ2::blockSolution(){
    unsigned_t t, i;
    Term qBv = mk_bv_value_uint64( _smt.vPiSort[0], stoi(_smt.pSolver->get_value( _smt.vvPi[0][0]).value<string>(10)));
    Term clause = mk_term(Kind::EQUAL, {_smt.vvPi[0][0], qBv});
    // expr clause = !(_smt.vvPi[0][0] == (int_t)(m.eval(_smt.vvPi[0][0]).get_numeral_int64()));   
    for (t = 0; t < _stbolsqParam.min_depth; ++t){
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            qBv = mk_bv_value_uint64( _smt.vPiSort[i], stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10)));
            clause = mk_term(Kind::AND, {clause,
                         mk_term(Kind::EQUAL, {_smt.vvPi[t][i], qBv})});
        }
    }
     _smt.pSolver->assert_formula(mk_term(Kind::NOT, {clause}));
}

bool sTBOLSQ2::extractModel(bool blockSol){
    unsigned_t i, gateTime, j, t, circuitDepth = 0, e, q1, q2, p1, p2;
    unsigned_t totalChangeCnt = 0, totalSwapCnt = 0;
    double_t totalDis = 0;
    string s;
    // fprintf(stdout, "[Info] Device Qubit Distance Info                              \n");
    // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){ 
    //     for(unsigned_t j = i; j < _pDevice->nQubit(); ++j){ 
    //         fprintf(stdout, "        - (%d,%d): dis %.2f\n", i, j, _pDevice->getDistance(i, j));
    //     }
    // }
    // collect circuit depth
    for ( i = 0; (i < _pCircuit->nGate());  ++i){
        Gate &gate = _pCircuit->gate(i);
        gateTime = stoi(_smt.pSolver->get_value( _smt.vTg[i]).value<string>(10)); 
        circuitDepth = (circuitDepth < gateTime) ? gateTime : circuitDepth;
        if(gate.nTargetQubit() > 1){
            q1 = gate.targetProgramQubit(0);
            q2 = gate.targetProgramQubit(1);
            p1 = _vvQubitRegion[q1][stoi(_smt.pSolver->get_value( _smt.vvPi[gateTime][q1]).value<string>(10))];
            p2 = _vvQubitRegion[q2][stoi(_smt.pSolver->get_value( _smt.vvPi[gateTime][q2]).value<string>(10))];
            // cerr << "_pDevice->getDistance(" << p1 << ", " << p2 << "): " << _pDevice->getDistance(p1, p2) << endl;
            // cerr << "totalDis: " << totalDis << endl;
            totalDis += _pDevice->getDistance(p1, p2);
        }
    }
    vector<bool> vIsMappingChange(circuitDepth+1, 0);
    vIsMappingChange[0] = 1;
    if(circuitDepth > 0){
        for (i = 0; i < circuitDepth; ++i){
            unordered_set<pair<int_t, int_t>> sChangePair;
            pair<int_t, int_t> changePair;
            for (j = 0; j < _pCircuit->nProgramQubit(); ++j){
                if (_smt.pSolver->get_value( _smt.vvChangeCnt[i][j]).is_true()){
                    vIsMappingChange[i+1] = 1;
                    q1 = _vvQubitRegion[j][stoi(_smt.pSolver->get_value( _smt.vvPi[i][j]).value<string>(10))];
                    q2 = _vvQubitRegion[j][stoi(_smt.pSolver->get_value( _smt.vvPi[i+1][j]).value<string>(10))];
                    if(q1 < q2){
                        changePair.first = q1;
                        changePair.second = q2;
                    }
                    else{
                        changePair.first = q2;
                        changePair.second = q1;
                    }
                    ++totalChangeCnt;
                    if (sChangePair.find(changePair) == sChangePair.end()){
                        if(_stbolsqParam.considerDis){
                            // cerr << "q1: " << _vvQubitRegion[i][q1] << ", q2: " << _vvQubitRegion[i][q2] << endl;
                            totalDis += _pDevice->getDistance(changePair.first, changePair.second);
                            // cerr << "totalDis: " << totalDis << ", dis for " << changePair.first << " " << changePair.second << " : " << _pDevice->getDistance(changePair.first, changePair.second) << endl;
                            ++ totalDis;
                        }
                        sChangePair.insert(changePair);
                    }
                }
            }
        }
    }
    else{
        totalDis = 0;
    }
    int_t p;
    if(!_stbolsqParam.considerDis || (totalDis < _bestDisSum)){
        fprintf(stdout, "[Info] Extract Model Info                              \n");
        vector<unsigned_t> vRemapTime(circuitDepth+1, 0);
        unsigned_t curTime = 0;
        _changeCntBound = totalChangeCnt;
        _vvQubitMapping.clear();
        _pCircuit->resetQubitRegion();
        _pCircuit->setCircuitDepth(circuitDepth+1);
        vector<unsigned_t> vQubitMapping(_pCircuit->nProgramQubit());
        for (t = 0; t <= circuitDepth; ++t){
            if(vIsMappingChange[t]){
                _vvQubitMapping.emplace_back(vector<unsigned_t>());
                for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
                    p = _vvQubitRegion[i][stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10))];
                    _pCircuit->addQubitRegion(i, p);
                    _vvQubitMapping.back().emplace_back(p);
                }
                curTime = t;
            }
            vRemapTime[t] = curTime;
        }
        _pCircuit->setCircuitDepth(_vvQubitMapping.size());
         // collect gate execution time
        for ( i = 0; (i < _pCircuit->nGate());  ++i){
            Gate &gate = _pCircuit->gate(i);
            gateTime = stoi(_smt.pSolver->get_value( _smt.vTg[i]).value<string>(10)); 
            gate.setExecutionTime(vRemapTime[gateTime]);
            for ( j = 0; j < gate.nTargetQubit();  ++j ){
                gate.setTargetPhysicalQubit(j, _vvQubitRegion[gate.targetProgramQubit(j)][stoi(_smt.pSolver->get_value( _smt.vvPi[gateTime][gate.targetProgramQubit(j)]).value<string>(10))]);
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
        // cerr << "line 741: hi" << endl;
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            // cerr << "i: " << i << endl; 
            _pCircuit->setInitialMapping(i, _vvQubitMapping[0][i]);
        }
        // cerr << "line 745: hi" << endl;
        for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
            // cerr << "i: " << i << endl;
            // cerr << "_vvQubitMapping.back()[i]: " << _vvQubitMapping.back()[i] << endl;
            _pCircuit->setFinalMapping(i,  _vvQubitMapping.back()[i]);
            // cerr << "_vvQubitMapping.back()[i]: " << _vvQubitMapping.back()[i] << endl;
        }
        // cerr << "qubit mapping " << endl;
        fprintf(stdout, "        - Qubit mapping: \n");
        for (t = 0; t < _vvQubitMapping.size(); ++t){
            fprintf(stdout, "        - Time %d: ",t);
            for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
                fprintf(stdout, "%d->%d ", i, _vvQubitMapping[t][i]);
            }
            fprintf(stdout, "\n");
        }
        if(_stbolsqParam.considerDis){
            _bestDisSum = totalDis;
            fprintf(stdout, "        - Total Dis: %.4f\n", _bestDisSum);
        }
        
        if(blockSol)
            blockSolution();
        return true;
    }

    if(blockSol)
        blockSolution();
    return false;
}

void sTBOLSQ2::increaseDepthBound(){
    while(_stbolsqParam.min_depth >= _stbolsqParam.max_depth){
        _stbolsqParam.max_depth_bit += _stbolsqParam.max_depth_expand_factor; 
        _stbolsqParam.max_depth = _stbolsqParam.max_depth << _stbolsqParam.max_depth_expand_factor;
    }
}

void sTBOLSQ2::updateSMT(unsigned_t d){
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints(d);
    fprintf(stdout, "[Info]          constructing mapping variable region    \n");
    addMappingRegionConstraints(d);
    fprintf(stdout, "[Info]          constructing qubit moving constraint \n");
    addQubitMovingDistanceConstraints(d);
    fprintf(stdout, "[Info]          constructing valid two-qubit gate constraint \n");
    addRestrictedValidTwoQubitGateConstraints(d);
    fprintf(stdout, "[Info]          constructing change pi constraint           \n");
    addChangePiConstraint(d);
}


MOLSQ_NAMESPACE_CPP_END
/***********************************************************************
  File        [ router.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ rotuer ]
  Synopsis    [ router class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "rolsq2/rolsq2.hpp"

MOLSQ_NAMESPACE_CPP_START

void rOLSQ2::run(unsigned_t bound){
    _timer.start(TimeUsage::FULL);
    _timer.setTimeout(_rolsqParam.timeout);
    _hasSol = false;
    computeQubitRegion();
    runSMT(bound);
}

void rOLSQ2::computeQubitRegion(){
    _vvQubitRegion.clear();
    _vvQubitRegion.resize(_pvInitialMapping->size());
    _vsQubitRegion.clear();
    _vsQubitRegion.resize(_pvInitialMapping->size());
    _sAllSwapRegion.clear();
    _vsProSwapRegion.clear();
    _vsProSwapRegion.resize(_pvInitialMapping->size());
    for(unsigned_t i = 0; i < _pvInitialMapping->size(); ++i){
        bfsSearch(i);
    }

    for(unsigned_t i = 1; i < _pvInitialMapping->size(); ++i){
        std::set_union(_vsProSwapRegion[i-1].begin(), _vsProSwapRegion[i-1].end(),
                _vsProSwapRegion[i].begin(), _vsProSwapRegion[i].end(),
                inserter(_sAllSwapRegion, _sAllSwapRegion.begin()));
    }
}

void rOLSQ2::computeQubitIdx(){
    _vvQubitRegion.clear();
    _vvQubitRegion.resize(_pvInitialMapping->size());
    unsigned_t idx;
    for(unsigned_t i = 0; i < _pvInitialMapping->size(); ++i){
        fprintf(stdout, "[Info]          Relabel region for qubit %d: ", i);
        idx = 0;
        for(unsigned_t j : _vsQubitRegion[i]){
            fprintf(stdout, " %d->%d ", idx, j);
            _vvQubitRegion[i].emplace_back(j);
            ++idx;
        }
        fprintf(stdout, "\n");
    }
}

void rOLSQ2::bfsSearch(unsigned_t qId){
    unsigned_t i, j, q, nIdx;
    vector<Node> vNode;
    priority_queue<int_t, std::vector<int>, std::greater<int> > priorityQ;
    vector<int_t> vBacktraceNode;
    unsigned_t cost, nSpanEdge, dis = 0;
    bool findQubit = false;
    set<unsigned_t> sQubitRegion;
    sQubitRegion.clear();
    vNode.emplace_back(-1, 0, (*_pvInitialMapping)[qId], -1, 0);
    sQubitRegion.insert((*_pvInitialMapping)[qId]);
    priorityQ.push(0);
    while(priorityQ.size() > 0 && (!findQubit || vNode[priorityQ.top()].dis <= dis)){
        nIdx = priorityQ.top();
        priorityQ.pop();
        Qubit& qubit = _device.qubit(vNode.at(nIdx).qIdx);
        // cerr << "expand node " << vNode.at(nIdx).idx << " phy q: " << vNode.at(nIdx).qIdx << " parent idx: " << vNode.at(nIdx).parentIdx << ", dis: " << vNode.at(nIdx).dis << endl;
        if(vNode.at(nIdx).qIdx == (*_pvFinalMapping)[qId]){
            if(!findQubit){
                dis = vNode.at(nIdx).dis;
                vBacktraceNode.emplace_back(vNode.at(nIdx).idx);
            }
            else if (dis == vNode.at(nIdx).dis){
                vBacktraceNode.emplace_back(vNode.at(nIdx).idx);
            }            
            findQubit = true;
        }
        nSpanEdge = qubit.vSpanEdge.size();
        for (j = 0; j < nSpanEdge; ++j){
            q = _device.edge(qubit.vSpanEdge[j]).qubitId1();
            if (qubit.idx == q){
                q = _device.edge(qubit.vSpanEdge[j]).qubitId2();
            }
            // cerr << "node idx " << vNode.at(nIdx).idx << endl;
            vNode.emplace_back(vNode.at(nIdx).idx, vNode.size(), q, qubit.vSpanEdge[j], vNode.at(nIdx).dis+1);
            // cerr << "add node " << vNode.back().idx << " qubit " << vNode.back().qIdx << " with distance " << vNode.back().dis << " parent idx " << vNode.back().parentIdx << endl;
            priorityQ.push(vNode.back().idx);
        }
        // cerr << "top node dis: " << vNode[priorityQ.top()].dis << endl;
    }
    ++dis;
    if(_rolsqParam.min_depth < dis){
        _rolsqParam.min_depth = dis;
    }
    // backtrack 
    int_t pIdx;
    for(int_t nIdx : vBacktraceNode){
        pIdx = nIdx;
        // cerr << "find node " << vNode.at(nIdx).idx << " qubit " << vNode.at(nIdx).qIdx << " with distance " << vNode.at(nIdx).dis << " parent idx " << vNode.at(nIdx).parentIdx << endl;
        while(pIdx > 0){
            // cerr << "pIdx: " << pIdx << endl;
            sQubitRegion.insert(vNode[pIdx].qIdx);
            // cerr << "add edge " << vNode[pIdx].parentEIdx << endl;
            pIdx = vNode[pIdx].parentIdx;
        }
    }
    // expend the qubit region by 1
    for(auto it = sQubitRegion.begin(); it != sQubitRegion.end(); ++ it){
        Qubit& qubit = _device.qubit(*it);
        _vsQubitRegion[qId].insert(*it);
        for(unsigned_t e : qubit.vSpanEdge){
            if(_vsQubitRegion[qId].count(_device.edge(e).qubitId1()) == 0){
                _vsQubitRegion[qId].insert(_device.edge(e).qubitId1());
                Qubit& qubit1 = _device.qubit(_device.edge(e).qubitId1());
                for(unsigned_t e1 : qubit1.vSpanEdge){
                    if(qubit1.idx != _device.edge(e1).qubitId1()){
                        if(_vsQubitRegion[qId].count(_device.edge(e1).qubitId1())){
                            _vsProSwapRegion[qId].insert(e1);
                        }
                    }
                    else if(qubit1.idx != _device.edge(e).qubitId2()){
                        if(_vsQubitRegion[qId].count(_device.edge(e1).qubitId2())){
                            _vsProSwapRegion[qId].insert(e1);
                        }
                    }
                }
            }
            else if(_vsQubitRegion[qId].count(_device.edge(e).qubitId2()) == 0){
                _vsQubitRegion[qId].insert(_device.edge(e).qubitId2());
                Qubit& qubit1 = _device.qubit(_device.edge(e).qubitId2());
                for(unsigned_t e1 : qubit1.vSpanEdge){
                    if(qubit1.idx != _device.edge(e1).qubitId1()){
                        if(_vsQubitRegion[qId].count(_device.edge(e1).qubitId1())){
                            _vsProSwapRegion[qId].insert(e1);
                        }
                    }
                    else if(qubit1.idx != _device.edge(e1).qubitId2()){
                        if(_vsQubitRegion[qId].count(_device.edge(e1).qubitId2())){
                            _vsProSwapRegion[qId].insert(e1);
                        }
                    }
                }
            }
        }
    }
}


void rOLSQ2::computeOverlapRegion(){
    _vmPhy2ProIdxChoiceIdx.clear();
    _vmPhy2ProIdxChoiceIdx.resize(_device.nQubit());
    for(unsigned_t i = 0; i < _vvQubitRegion.size(); ++i){
        for(unsigned_t j = 0; j < _vvQubitRegion[i].size(); ++j){
            _vmPhy2ProIdxChoiceIdx[_vvQubitRegion[i][j]][i] = j;
        }    
    }
    // cerr << "for debug: " << endl;
    // for(unsigned_t i = 0; i < _vmPhy2ProIdxChoiceIdx.size(); ++i){
    //     cerr << "for phy qubit " << i << ":" << endl;
    //     for(auto & j : _vmPhy2ProIdxChoiceIdx[i]){
    //         cerr << "idx for pro qubit " << j.first << " is " << j.second << " (correspondding in vvQubitRegion: " << _vvQubitRegion[j.first][j.second] << ")" << endl;
    //         assert(i == _vvQubitRegion[j.first][_vmPhy2ProIdxChoiceIdx[i][j.first]]);
    //     }    
    // }
}

void rOLSQ2::runSMT(unsigned_t bound){
    fprintf(stdout, "[Info] rOLSQ2                        \n");
    bool solve = false;
    _iter = 0;
    while (!solve){
        fprintf(stdout, "[Info] Iter %d: Generating formulation                        \n", _iter);
        generateFormulation();
        fprintf(stdout, "[Info] Iter %d: Optimizting model                             \n", _iter);
        solve = optimize(bound);
        if(!solve){
            increaseDepthBound(); 
            expandQubitRegion();
        }
        ++_iter;
        // assert(_iter < 2);
    }
}

void rOLSQ2::generateFormulation(){
    _smt.reset(_rolsqParam.timeout, _iter);
    fprintf(stdout, "[Info]          relabeling qubit index                       \n");
    computeQubitIdx();
    computeOverlapRegion();
    fprintf(stdout, "[Info]          constructing variables                       \n");
    constructVariable();
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInitialMappingConstraints();
    fprintf(stdout, "[Info]          constructing mapping variable region    \n");
    addInjectiveMappingConstraints();
    addMappingRegionConstraints();
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraints();
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraints();
}

void rOLSQ2::constructVariable(){
    unsigned_t bit_length_pi;
    _smt.vvPi.reserve(_rolsqParam.max_depth);
    _smt.vvSigma.reserve(_rolsqParam.max_depth);

    string s;
    unsigned_t i, j;

    // Create a bit-vector sort of size 1.
    const Sort sortbool = mk_bool_sort();

    for (i = 0; i < _vvQubitRegion.size(); ++i){
        bit_length_pi = ceil(log2(_vvQubitRegion[i].size() + 1));
        _smt.vPiSort.emplace_back(mk_bv_sort( bit_length_pi));
    }

    for (i = 0; i < _rolsqParam.max_depth; ++i){
        _smt.vvPi.emplace_back(vector<Term>());
        for (j = 0; j < _pvInitialMapping->size(); ++j){
            s = "map_t" + to_string(i) + "_q" + to_string(j);
            // _smt.vvPi[i].emplace_back(_smt.c.bv_const(s.c_str(), bit_length_pi));
            _smt.vvPi[i].emplace_back(mk_const(_smt.vPiSort[j], s.c_str()));
        }
    }

    for (i = 0; i < _rolsqParam.max_depth; ++i){
        _smt.vvSigma.emplace_back(vector<Term>());
        for (j = 0; j < _device.nEdge(); ++j){
            s = "ifswap_t" + to_string(i) + "_e" + to_string(j);
            // _smt.vvSigma[i].emplace_back(_smt.c.bool_const(s.c_str()));
            _smt.vvSigma[i].emplace_back(mk_const(sortbool, s.c_str()));
        }
    }
}

void rOLSQ2::addInitialMappingConstraints(){
    Term qBv;
    for (unsigned_t i = 0; i < _pvInitialMapping->size(); ++i){
        qBv = mk_bv_value_uint64(_smt.vPiSort[i], _vmPhy2ProIdxChoiceIdx[(*_pvInitialMapping)[i]][i]);
        // cerr << "at time 0, qubit " << i << " is mapped to " << (*_pvInitialMapping)[i] << "(" << _vmPhy2ProIdxChoiceIdx[(*_pvInitialMapping)[i]][i] << ")" << endl;
         _smt.pSolver->assert_formula(mk_term(Kind::EQUAL, {_smt.vvPi[0][i], qBv}));
    }
}

void rOLSQ2::addInjectiveMappingConstraints(unsigned_t boundOffset){
    Term cond1, cond2, q1Bv, q2Bv;
    unsigned_t end = _rolsqParam.min_depth, begin = 1;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (unsigned_t i = 0; i < _vmPhy2ProIdxChoiceIdx.size(); ++i){
        for (auto iterj = _vmPhy2ProIdxChoiceIdx[i].begin(); iterj != _vmPhy2ProIdxChoiceIdx[i].end(); iterj++){
            q1Bv = mk_bv_value_uint64(_smt.vPiSort[iterj -> first], iterj -> second);
            for(auto iterk = next(iterj, 1); iterk != _vmPhy2ProIdxChoiceIdx[i].end(); ++iterk){
                q2Bv = mk_bv_value_uint64(_smt.vPiSort[iterk->first], iterk->second);
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

void rOLSQ2::addMappingRegionConstraints(unsigned_t boundOffset){
    Term qBv, zero;
    unsigned_t end = _rolsqParam.min_depth, begin = 1;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    // add mapping region constraint
    for (unsigned_t j = 0; j < _vvQubitRegion.size(); ++j){
        qBv = mk_bv_value_uint64(_smt.vPiSort[j], _vvQubitRegion[j].size());
        // cerr << "qubit " << j << "'s choice is less than " << _vvQubitRegion[j].size() << endl;
        zero = mk_bv_zero( _smt.vPiSort[j]);
        for (unsigned_t t = begin; t < end; ++t){
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULE, {zero, _smt.vvPi[t][j]}));
             _smt.pSolver->assert_formula(mk_term(Kind::BV_ULT, {_smt.vvPi[t][j], qBv}));
        }
    }
}

void rOLSQ2::addSwapConstraints(unsigned_t boundOffset){
    // No swap for t<s
    unsigned_t end = _rolsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    unsigned_t i, j, t, e, tt, q1, q2;
    // swap gates can not overlap with swap in space
    for (t = begin; t < end; ++t){
        for (unsigned_t e : _sAllSwapRegion){
            q1 = _device.edge(e).qubitId1();
            q2 = _device.edge(e).qubitId2();
            for (unsigned_t ee : _device.qubit(q1).vSpanEdge){
                if ((ee < e) && (_sAllSwapRegion.find(ee) != _sAllSwapRegion.end())){
                     _smt.pSolver->assert_formula(
                            mk_term(Kind::OR, 
                                {mk_term(Kind::NOT, {_smt.vvSigma[t][e]}),
                                mk_term(Kind::NOT, {_smt.vvSigma[t][ee]})}));
                }
            }
            for (unsigned_t ee : _device.qubit(q2).vSpanEdge){
                if ((ee < e) && (_sAllSwapRegion.find(ee) != _sAllSwapRegion.end())){
                     _smt.pSolver->assert_formula(
                        mk_term(Kind::OR, 
                            {mk_term(Kind::NOT, {_smt.vvSigma[t][e]}),
                            mk_term(Kind::NOT, {_smt.vvSigma[t][ee]})}));
                }
            }
        }
    }
}

void rOLSQ2::addTransformationConstraints(unsigned_t boundOffset){
    unsigned_t i, j, t, e, nSpanEdge, idx1, idx2, idx;
    Term q1Bv, q2Bv, cond, clause, qBv;
    unsigned_t end = _rolsqParam.min_depth, begin = 0;
    if(boundOffset > 0){
        begin = end;
        end += boundOffset;
    }
    for (t = begin; t < end; ++t){
        for (i = 0; i < _vsProSwapRegion.size(); ++i){
            // mapping change by swap
            for (int_t e : _vsProSwapRegion[i]){
                idx1 = _vmPhy2ProIdxChoiceIdx[_device.edge(e).qubitId1()][i];
                idx2 = _vmPhy2ProIdxChoiceIdx[_device.edge(e).qubitId2()][i];
                q1Bv = mk_bv_value_uint64(_smt.vPiSort[i], idx1);
                q2Bv = mk_bv_value_uint64(_smt.vPiSort[i], idx2);
                cond = mk_term(Kind::NOT,
                                        {mk_term(Kind::AND,
                                            {_smt.vvSigma[t][e], 
                                            mk_term(Kind::EQUAL, {_smt.vvPi[t][i], q1Bv})})});
                 _smt.pSolver->assert_formula(mk_term(Kind::OR,
                                                {cond, 
                                                mk_term(Kind::EQUAL, {_smt.vvPi[t+1][i], q2Bv})}));
                cond = mk_term(Kind::NOT,
                        {mk_term(Kind::AND,
                        {_smt.vvSigma[t][e], 
                        mk_term(Kind::EQUAL, {_smt.vvPi[t][i], q2Bv})})});
                 _smt.pSolver->assert_formula(mk_term(Kind::OR,
                                                {cond, 
                                                mk_term(Kind::EQUAL, {_smt.vvPi[t+1][i], q1Bv})}));
            }
             // mapping not change by swap
            for (int_t q : _vvQubitRegion[i]){
                Qubit& qubit = _device.qubit(q);
                nSpanEdge = qubit.vSpanEdge.size();
                clause = _smt.vvSigma[t][qubit.vSpanEdge[0]];
                for (e = 1; e < nSpanEdge; ++e){
                    if(_vsProSwapRegion[i].find(qubit.vSpanEdge[e]) != _vsProSwapRegion[i].end()){
                        clause = mk_term(Kind::OR, {clause, _smt.vvSigma[t][qubit.vSpanEdge[e]]});
                    }
                }
                idx = _vmPhy2ProIdxChoiceIdx[q][i];
                qBv = mk_bv_value_uint64(_smt.vPiSort[i], idx);
                clause = mk_term(Kind::AND,
                            {mk_term(Kind::NOT, {clause}), 
                            mk_term(Kind::EQUAL, {_smt.vvPi[t][i], qBv})});
                 _smt.pSolver->assert_formula(
                                mk_term(Kind::OR,
                                    {mk_term(Kind::NOT, {clause}), 
                                    mk_term(Kind::EQUAL, {_smt.vvPi[t+1][i], qBv})}));
            }
        }
    }
}

void rOLSQ2::addFinalMappingConstraints(){
    Term qBv;
    unsigned_t d = _rolsqParam.min_depth;
    for (unsigned_t i = 0; i < _pvFinalMapping->size(); ++i){
        qBv = mk_bv_value_uint64(_smt.vPiSort[i], _vmPhy2ProIdxChoiceIdx[(*_pvFinalMapping)[i]][i]);
         _smt.pSolver->assert_formula(mk_term(Kind::EQUAL, {_smt.vvPi[d][i], qBv}));
    }
}

void rOLSQ2::addSwapCountConstraints(unsigned_t bound){
    vector<Term> vSwap;
    unsigned_t t, e;
    // cout << "add swap constraints" << endl;
    for (t = 0; t < _rolsqParam.min_depth; ++t){
        for (auto e : _sAllSwapRegion){
            vSwap.push_back(_smt.vvSigma[t][e]);
        }
    }

    PB2CNF pb2cnf;
    vector<int_t> vSigmaLit;
    unsigned_t firstFreshVariable = 1 + _rolsqParam.min_depth*_sAllSwapRegion.size(), newFreashVariable;
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
    for(vector<int_t>& clause : formula){
        vOrs.clear();
        for(int& lit : clause){
            var = abs(lit);
            if(var < firstFreshVariable){
                sigmaT = (var - 1) / _sAllSwapRegion.size();
                sigmaE = (var - 1) % _sAllSwapRegion.size();
                // cerr << "var: " << var << " sigmaT: " << sigmaT << ", sigmaE: " << sigmaE << endl;
                auto iter = _sAllSwapRegion.begin();
                for(unsigned_t i = 0; i < sigmaE; ++i){
                    ++iter;
                }
                if (lit < 0){ 
                    // vOrs.emplace_back(mk_term(Kind::NOT, _smt.vvSigma[sigmaT][(*iter)]));
                    vOrs.emplace_back(mk_term(Kind::NOT, {vSwap[var-1]}));
                }
                else{
                    // vOrs.emplace_back(_smt.vvSigma[sigmaT][(*iter)]);
                    vOrs.emplace_back(vSwap[var-1]);
                }
                // cerr << "finish" << endl;
            }
            else{
                // cerr << "add ancilla" << endl;
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
                // cerr << "finish" << endl;
            }
        }    
        // cerr << "add or list" << endl;
        Term cnf = vOrs[0];
        for (unsigned_t i = 1; i < vOrs.size(); ++i){
            cnf = mk_term(Kind::OR, {vOrs[i], cnf});
        }
         _smt.pSolver->assert_formula(cnf);
        //  _smt.pSolver->assert_formula(bitwuzla_mk_term(_smt.pSolver, Kind::OR, vOrs.size(), vOrs));
    }
    // cerr << "finish add swap count constraint" << endl;
}

bool rOLSQ2::checkModel(){
    Result status = _smt.pSolver->check_sat();
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

bool rOLSQ2::optimizeDepth(){
    bool success, find_min_depth = false;
    unsigned_t i;
    while(!find_min_depth && _rolsqParam.min_depth < _rolsqParam.max_depth ){
        cerr << "_rolsqParam.min_depth: " << _rolsqParam.min_depth << ", _rolsqParam.max_depth: " << _rolsqParam.max_depth << endl;
        fprintf(stdout, "[Info]          Trying to optimize for depth bound %d            \n", _rolsqParam.min_depth);
        _timer.start(TimeUsage::PARTIAL);
        if(!_hasSol){
            _smt.resetTimeState(_rolsqParam.timeout);
        }
        _smt.pSolver->push(1);
        addFinalMappingConstraints();
        success = checkModel();
        fprintf(stdout, "[Info]          Optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("optimizing depth", TimeUsage::PARTIAL);
        _timer.showUsage("optimizing depth", TimeUsage::FULL);
        if (success){
            extractModel();
            _smt.pSolver->pop(1);
            find_min_depth = true;
        }
        else{
            _smt.pSolver->pop(1);
            if(_rolsqParam.min_depth + 1 < _rolsqParam.max_depth){
                updateSMT(1);
            }
            ++_rolsqParam.min_depth;
        }
    }
    if (find_min_depth){
        return true;
    }
    else{
        return false;
    }
}

bool rOLSQ2::optimize(unsigned_t bound){
    bool success_optimize;
    success_optimize = optimizeDepth();
    // getchar();
    if (!success_optimize){
        return false;
    }
    if(_vSwap.size() == 0){
        return true;
    }
    unsigned_t step = 2, lower_swap_bound = 0, upper_swap_bound = min(bound, (unsigned_t)_vSwap.size()-1);
    cout << "swap lower bound = " << lower_swap_bound << " , swap upper bound = " << upper_swap_bound << endl;
    bool reduce_swap = true;
    bool firstRun = true;
    while (reduce_swap && (!_timer.isTimeout() || !_hasSol)){
        // cout << "enter loop" << endl;
        addFinalMappingConstraints();
        reduce_swap = optimizeSwapForDepth(lower_swap_bound, upper_swap_bound, firstRun);
        upper_swap_bound = _vSwap.size() - 1;
        firstRun = false;
        // getchar();
        if(reduce_swap){
            fprintf(stdout, "[Info] Successfully reduce SWAP count. Go to the next run.            \n");
            fprintf(stdout, "[Info] Solving with depth %d            \n", _rolsqParam.min_depth);
            _timer.start(TimeUsage::PARTIAL);
            fprintf(stdout, "[Info] Generating formulation                        \n");
            if(_rolsqParam.min_depth + step < _rolsqParam.max_depth){
                updateSMT(step);
                _rolsqParam.min_depth += step;
            }
            else if(_rolsqParam.max_depth - _rolsqParam.min_depth > 1){
                updateSMT(_rolsqParam.max_depth - _rolsqParam.min_depth - 1);
                _rolsqParam.min_depth += step;
            }
            else{
                _rolsqParam.min_depth += step;
                increaseDepthBound();
                generateFormulation();
            }
            _timer.showUsage("Generating formulation", TimeUsage::PARTIAL);
            _timer.start(TimeUsage::PARTIAL);
        }
    }
    return true;
}

bool rOLSQ2::optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun){
    unsigned_t swap_bound = upper_swap_bound, swap_step;
    bool find_min_swap = false;
    bool success;
    while (!find_min_swap && lower_swap_bound <= swap_bound && swap_bound <= upper_swap_bound && (!_timer.isTimeout() || !_hasSol)){
        if(!_hasSol){
            _smt.resetTimeState(_rolsqParam.timeout);
        }
        fprintf(stdout, "[Info]          trying to optimize for swap bound %d            \n", swap_bound);
        _timer.start(TimeUsage::PARTIAL);
        _smt.pSolver->push(1);
        addSwapCountConstraints(swap_bound);
        success = checkModel();
        fprintf(stdout, "[Info]          optimization results: %s                         \n", success ? "success" : "fail");
        _timer.showUsage("rOLSQ2: optimizing swap", TimeUsage::PARTIAL);
        _timer.showUsage("rOLSQ2: optimizing swap", TimeUsage::FULL);
        if (success){
            extractModel();
            _hasSol = true;
            assert(_vSwap.size() <= swap_bound);
            swap_bound = _vSwap.size();
            --swap_bound;
            _smt.pSolver->pop(1);
        }
        else{
            _smt.pSolver->pop(1);
            if(swap_bound < upper_swap_bound){
                find_min_swap = true;
            }
            else if (!firstRun){
                break;
            }
            else if (firstRun && swap_bound == upper_swap_bound){
                // cout << "line 518: increase swap bound" << endl;
                lower_swap_bound = upper_swap_bound + 1;
                upper_swap_bound = upper_swap_bound + 5; // make a big jump when unsat
                swap_bound = upper_swap_bound;
            }
        }
    }
    return find_min_swap;
}

void rOLSQ2::extractModel(){
        fprintf(stdout, "[Info] Extract Model Info                              \n");
    // model m = _smt.smtSolver.get_model();
    unsigned_t i, gateTime, q, j, swapId = 0, qId, t, e;
    // collect gate execution time
    // get SWAP gate
    _vSwap.clear();
    string s;
    vector<unsigned_t> swapTargetQubit(2,0);
    for (t = 0; t < _rolsqParam.min_depth; ++t){
        for (unsigned_t e: _sAllSwapRegion){
            // cout << "t = " << t << ", e = " << e << " " << m.eval(_smt.vvSigma[t][e]).bool_value() << endl;
            // if (m.eval(_smt.vvSigma[t][e]).bool_value() == 1){
            if (_smt.pSolver->get_value( _smt.vvSigma[t][e]).is_true()){
                swapTargetQubit[0] = _device.edge(e).qubitId1();
                swapTargetQubit[1] = _device.edge(e).qubitId2();
                _vSwap.emplace_back(swapTargetQubit[0], swapTargetQubit[1]);
                if (_verbose > 0){
                    fprintf(stdout, "        - SWAP Gate %d, time: %d, target qubit: %d %d\n", swapId, t, swapTargetQubit[0], swapTargetQubit[1]);
                }
                ++swapId;
            }
        }
    }
    // set initial and final mapping
    if (_verbose > 0){
        fprintf(stdout, "        - Qubit mapping: \n");
        for (t = 0; t <= _rolsqParam.min_depth; ++t){
        fprintf(stdout, "        - Time %d: ",t);
            for (i = 0; i < _pvInitialMapping->size(); ++i){
                qId = stoi(_smt.pSolver->get_value( _smt.vvPi[t][i]).value<string>(10));
                assert(qId < _vvQubitRegion[i].size());
                fprintf(stdout, "%d->%d(%d) ", i, _vvQubitRegion[i][qId], qId);
            }
        fprintf(stdout, "\n");
        }
    }
}


void rOLSQ2::expandQubitRegion(){
    unsigned_t i, j;
    for(i = 0; i < _vvQubitRegion.size(); ++i){
        for(j = 0; j < _vvQubitRegion[i].size(); ++j){
            Qubit & qubit = _device.qubit(_vvQubitRegion[i][j]);
            for(unsigned_t e : qubit.vSpanEdge){
                Edge & edge = _device.edge(e);
                _vsQubitRegion[i].insert(edge.qubitId1());
                _vsQubitRegion[i].insert(edge.qubitId2());
                _vsProSwapRegion[i].insert(e);
                _sAllSwapRegion.insert(e);
            }
        }    
    }
}

void rOLSQ2::increaseDepthBound(){
    while(_rolsqParam.min_depth >= _rolsqParam.max_depth){
        _rolsqParam.max_depth_bit += _rolsqParam.max_depth_expand_factor; 
        _rolsqParam.max_depth = _rolsqParam.max_depth << _rolsqParam.max_depth_expand_factor;
    }
}

void rOLSQ2::updateSMT(unsigned_t d){
    fprintf(stdout, "[Info]          constructing injective mapping constraint    \n");
    addInjectiveMappingConstraints(d);
    fprintf(stdout, "[Info]          constructing mapping variable region    \n");
    addMappingRegionConstraints(d);
    fprintf(stdout, "[Info]          constructing swap overlapping constraint     \n");
    addSwapConstraints(d);
    fprintf(stdout, "[Info]          constructing mapping transformation constraint\n");
    addTransformationConstraints(d);
}


MOLSQ_NAMESPACE_CPP_END
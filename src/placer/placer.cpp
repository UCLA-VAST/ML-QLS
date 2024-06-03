/***********************************************************************
  File        [ placer.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ placer ]
  Synopsis    [ placer class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 5, Sep., 2023 ]
***********************************************************************/
#include "placer/placer.hpp"

MOLSQ_NAMESPACE_CPP_START

void Placer::run(Circuit& cir, Device& device, bool givenInitialMapping){
    _pCircuit = &cir;
    _pDevice = &device;
    // _pvvGateBlock = &vvGateBlock;
    _pvvGateBlock = nullptr;
    _timer.start(TimeUsage::FULL);
    _isGivenInitialMapping = givenInitialMapping;
    unsigned_t i;
    _optSol = Solution();
    // construct QGG and rQGG for each subcircuit
    // generate cost function
    _vvWeightedEdge.clear();
    _vvWeightedEdge.reserve(1);
    // preComputeWeight();
    // cir.printCircuitInitialMapping();
    if(!_isGivenInitialMapping && _pCircuit->sQubitRegion(0).size() > 0){
        fprintf(stdout, "[Info] SA-Based Placer: Find initial mapping\n");
        decideQubitInitialMapping();
    }
    ////
    // for (i = 0; i < vvGateBlock.size(); ++i){
    //     for( unsigned_t j = 0; j < vvGateBlock[i].size(); ++j){
    //         Gate & gate = _pCircuit->gate(vvGateBlock[i][j]);
    //         _vvWeightedEdge[0].emplace_back(1, make_pair(gate.targetProgramQubit(0), gate.targetProgramQubit(1)));
    //     }
    // }
    _vvWeightedEdge.emplace_back(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >());
    for (i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        if(gate.nTargetQubit() > 1)
            _vvWeightedEdge[0].emplace_back(1, make_pair(gate.targetProgramQubit(0), gate.targetProgramQubit(1)));
    }
    ////

    // for (i = 0; i < _pvvGateBlock->size(); ++i){
    //     fprintf(stdout, "       - Generate QGG %d\n", i);
    //     _vvWeightedEdge.emplace_back(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >());
    //     // cerr << "finish emplace back vector" << endl;
    //     constructWeightedEdge(i, _vvWeightedEdge[i]); 
    // }
    // cerr << "_vvWeightedEdge->size(): " << _vvWeightedEdge.size() << endl;
    fprintf(stdout, "[INFO] SA-Based Placer: Start SA\n");
    runSA();
    fprintf(stdout, "[INFO] SA-Based Placer: Refinement \n");
    runRefinement();
    fprintf(stdout, "[INFO] SA-Based Placer: Final cost: %4f\n", _optSol.cost);
    fprintf(stdout, "[INFO] SA-Based Placer: Set initial mapping\n");
    setMapping();
}

void Placer::constructWeightedEdge(unsigned_t curI, vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge){
    int_t i = 0, j, q0, q1;
    vector<unsigned_t> vQubitLevel;
    unordered_map<pair<unsigned_t , unsigned_t>, double_t, pair_hash> umEdge2Weight;
    pair<unsigned_t , unsigned_t> edge;
    while(i < _qggParam.max_gate_level && curI < _pvvGateBlock->size()){
        for(j = 0; j < (*_pvvGateBlock)[curI].size(); ++j){
            Gate & gate = _pCircuit->gate((*_pvvGateBlock)[curI][j]);
            if(gate.nTargetQubit()>1){
                q0 = gate.targetProgramQubit(0);
                q1 = gate.targetProgramQubit(1);
                edge = (q0 < q1) ? make_pair(q0,q1) : make_pair(q1,q0); 
                if (umEdge2Weight.find(edge) == umEdge2Weight.end()){
                    umEdge2Weight[edge] = _vPrecomputeWeight[i];
                }
                else{
                    umEdge2Weight[edge] = umEdge2Weight[edge] + _vPrecomputeWeight[i];
                }
                // cout << "umEdge2Weight: " << umEdge2Weight[edge] << endl;
            }
        }
        ++i;
        ++curI;
    }
    // cerr << "before clear" << endl;
    vpWeightEdge.clear(); 
    // cerr << "before reserve " << umEdge2Weight.size() << endl;
    vpWeightEdge.reserve(umEdge2Weight.size());
    // cerr << "after reserve" << endl;
    // unsigned_t s = 0;
    for (auto & keyItemPair : umEdge2Weight){
        // cerr << "hi " << keyItemPair.second << " " << keyItemPair.first.first << " " << keyItemPair.first.second << endl;
        vpWeightEdge.emplace_back(keyItemPair.second, keyItemPair.first);
        // cerr << "hi " << vpWeightEdge[s].first << " " <<  vpWeightEdge[s].second.first << " " <<  vpWeightEdge[s].second.second << endl;
        // ++s;
    }
    // cerr << "after insert" << endl;
    sort(vpWeightEdge.rbegin(), vpWeightEdge.rend());
    // cerr << "after sort" << endl;
    printEdgeWeight(vpWeightEdge);
}

void Placer::preComputeWeight(){
    _vPrecomputeWeight.clear();
    _vPrecomputeWeight.reserve(_qggParam.max_gate_level+1);
    _vPrecomputeWeight.emplace_back(1);
    // cout << "_qggParam.max_gate_level: " << _qggParam.max_gate_level << endl;
    // cout << "Grapher precompute weight: " << endl;
    // cout << "level " << 0 << ", weight: " << _vPrecomputeWeight[0] << endl;
    for( unsigned_t i = 1; i <= _qggParam.max_gate_level; ++i ){
        _vPrecomputeWeight.emplace_back(_vPrecomputeWeight[i-1]*_qggParam.gate_decay_factor);
        // cout << "level " << i << ", weight: " << _vPrecomputeWeight[i] << endl;
        // getchar();
    }
}

void Placer::printEdgeWeight(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge){
    fprintf(stdout, "===========================================\n");
    fprintf(stdout, "[INFO] Print Graph Info\n");
    for (pair<double_t, pair<unsigned_t, unsigned_t> >& pWeightEdge : vpWeightEdge){
        fprintf(stdout, "       - Edge(%d, %d), weight: %4f\n", pWeightEdge.second.first, pWeightEdge.second.second, pWeightEdge.first);
    }
    fprintf(stdout, "===========================================\n");
}

void Placer::runSA() {
    unsigned_t iter = 0, noImproveIter = 0;
    unsigned_t numProcess = 0;
    double_t newCost = 0.0, lastCost = 0;
    bool firstTime = true;
    while(!_timer.isTimeout() && noImproveIter < 3){
        initSol();
        if (firstTime){
            fprintf(stdout, "[INFO] SA-Based Placer: Init cost %4f\n", _optSol.cost);
            firstTime = false;
        }
        
        // cerr << "finish generate initial solution" << endl;
        int_t rejectNum = 0;
        ++numProcess;
        _saParam.n = 0;
        // cout << "Print qubit mapping: " << endl; 
        // cout << "pro -> phy: " << endl; 
        // for (unsigned_t i = 0; i < _curSol.vvProQ2PhyQ[0].size(); ++i){
        //     cout << i << " -> " << _curSol.vvProQ2PhyQ[0][i] << endl;
        // }
        while(_saParam.t > _saParam.t_frozen){
            ++_saParam.n;
            ++iter;
            _saParam.delta_cost_cnt = 0;
            _saParam.delta_sum = 0;
            rejectNum = 0;
            for(int_t i = 0; i < _saParam.l; ++i){
                getMovement();
                // if(_pDevice->nQubit() > 20){
                //     cerr << "makeMovement ";
                // }
                makeMovement(); // make movement and calculate cost difference
                // if(_pDevice->nQubit() > 20){
                //     cerr << "finishMovement ";
                // }
                _saParam.delta_cost_cnt++;
                _saParam.delta_sum += fabs(_saParam.delta);

                if(_saParam.delta <= 0){
                    _curSol.cost += _saParam.delta;
                    if (_curSol.cost <= _optSol.cost){
                        // cout << "update new cost: " << _curSol.cost << endl;
                        // cout << "Print qubit mapping: " << endl; 
                        // cout << "pro -> phy: " << endl; 
                        // for (unsigned_t i = 0; i < _curSol.vvProQ2PhyQ[0].size(); ++i){
                        //     cout << i << " -> " << _curSol.vvProQ2PhyQ[0][i] << endl;
                        // }
                        updateOptSol();
                    }
                }
                else{   // delta > 0
                    if(acceptWorseSol()){
                        _curSol.cost += _saParam.delta;
                    }
                    else{
                        // undo current movement
                        recover();
                        ++rejectNum;
                    }
                }
            } 
            updateT();
            
            if(_saParam.n > _saParam.iter_limit) break;
        }
        if(_saParam.find_opt_sol) break;
        fprintf(stdout, "[INFO] SA-Based Placer: Iter %d\n", iter);
        fprintf(stdout, "                        Curret cost: %4f\n", _optSol.cost);
        fprintf(stdout, "                        Refinement \n");
        runRefinement();
        fprintf(stdout, "                        Final cost: %4f\n", _optSol.cost);
        if(fabs(lastCost - _optSol.cost) < 0.001){
            ++noImproveIter;
        }
        else{
            noImproveIter = 0;
        }
        lastCost = _optSol.cost;
    }

    return;
}

void Placer::getCost(bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    sol.cost = 0; 
    unsigned_t i = 0, verticalCost = 0;
    _vCostForEachTime.clear();
    _vCostForEachQubit.clear();
    _vCostForEachTime.resize(_vvWeightedEdge.size(), 0);
    _vCostForEachQubit.resize(_pCircuit->nProgramQubit(), 0);
    // calculate horizontal distance cost (between qubits)
    for (i = 0; i < _vvWeightedEdge.size(); ++i){
        calCostForTime(i, isCurSol);
        sol.cost += _vCostForEachTime[i];
        // cerr << "_vCostForEachTime[" << i << "]: " << _vCostForEachTime[i] << endl; 
    }
    // cerr << "Finish _vCostForEachTime " << endl;
    // calculate vertical distance cost (between time slots)
    for (i = 0; i < _pCircuit->nProgramQubit(); ++i){
        calCostForQubit(i, isCurSol);
        verticalCost += _vCostForEachQubit[i];
        // cerr << "_vCostForEachQubit[" << i << "]: " << _vCostForEachQubit[i] << endl; 
    }
    sol.cost += (verticalCost * (double_t)_saParam.verticalCostWeight);
    // cerr << "total cost: " << _curSol.cost << endl;
}

void Placer::calCostForTime(unsigned_t i, bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    _vCostForEachTime[i] = 0;
    unsigned_t qId0, qId1;
    for (pair<double_t, pair<unsigned_t, unsigned_t> > & weightedEdge : _vvWeightedEdge[i]){
        qId0 = sol.vvProQ2PhyQ[i][weightedEdge.second.first];
        qId1 = sol.vvProQ2PhyQ[i][weightedEdge.second.second];
        // if(qId1 >= _pDevice->nQubit()){
        //    cerr << "isCurSol? " << isCurSol << endl; 
        //    cerr << "cal cost for sum of distance for pro " << weightedEdge.second.first << " and " << weightedEdge.second.second << endl; 
        //    cerr << "cal cost for sum of distance for phy " << qId0 << " and " << qId1 << endl; 
        //    cerr << "print qubit mapping swap" << endl;
        //     cerr << "phyQ2proQ: (" << sol.vvPhyQ2ProQ[i].size() <<  ")";
        //     for (unsigned_t j = 0; j < _pDevice->nQubit(); ++j){
        //         cerr << j << " -> " << sol.vvPhyQ2ProQ[i][j] << endl;
        //     }
        //     cerr << "proQ2phyQ: (" << sol.vvProQ2PhyQ[i].size() <<  ")";
        //     for (unsigned_t j = 0; j < sol.vvProQ2PhyQ[i].size(); ++j){
        //         cerr << j << " -> " << sol.vvProQ2PhyQ[i][j] << endl;
        //     }
        // }
        // cerr << "cal cost for sum of distance for pro " << weightedEdge.second.first << " and " << weightedEdge.second.second << endl;
        assert(qId0 >= 0);
        assert(qId1 >= 0);
        assert(qId0 < _pDevice->nQubit());
        assert(qId1 < _pDevice->nQubit());
        _vCostForEachTime[i] += (weightedEdge.first * ((double_t)_pDevice->getDistance(qId0, qId1)));
    }
}

// input i: id of program qubit
void Placer::calCostForQubit(unsigned_t i, bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    _vCostForEachQubit[i] = 0;
    unsigned_t qId0, qId1;
    for (unsigned_t j = 0; j < _vvWeightedEdge.size()-1; ++j){
        qId0 = sol.vvProQ2PhyQ[j][i];
        qId1 = sol.vvProQ2PhyQ[j+1][i];
        if (qId0 != qId1){
            assert(qId0 >= 0);
            assert(qId1 >= 0);
            assert(qId0 < _pDevice->nQubit());
            assert(qId1 < _pDevice->nQubit());
            _vCostForEachQubit[i] += (_pDevice->getDistance(qId0, qId1));
        }
    }
}

void Placer::getMovement(){
    bool type = (((double_t)(rand() % _saParam.random_range) / _saParam.random_range) <= _saParam.prob_move_within_region);
    unsigned_t m = rand() % _vvWeightedEdge.size();
    unsigned_t qId = rand() % _pCircuit->nProgramQubit(); // randomly choose a program qubit to move
    unsigned_t neighborPhyQId = _curSol.vvProQ2PhyQ[m][qId]; 
    assert(neighborPhyQId >= 0);
    assert(neighborPhyQId < _pDevice->nQubit());
    get<0>(_saParam.movement) = m;
    get<1>(_saParam.movement) = _curSol.vvProQ2PhyQ[m][qId];
    unordered_set<int_t>& sQubitRegion = _pCircuit->sQubitRegion(qId);
    if(sQubitRegion.size() == 0 || !_restrictMovement){
        type = false;
    }
    int_t qId2 = _curSol.vvPhyQ2ProQ[m][neighborPhyQId];
    // cerr << "neighborPhyQId == _curSol.vvProQ2PhyQ[m][qId]: " << (neighborPhyQId == _curSol.vvProQ2PhyQ[m][qId]) << endl;
    // cerr <<  "!(qId2 == -1 || _pCircuit->sQubitRegion(qId2).find(_curSol.vvProQ2PhyQ[m][qId]) != _pCircuit->sQubitRegion(qId2).end()): " << !(qId2 == -1 || _pCircuit->sQubitRegion(qId2).find(_curSol.vvProQ2PhyQ[m][qId]) != _pCircuit->sQubitRegion(qId2).end()) << endl;
    // while( neighborPhyQId == _curSol.vvProQ2PhyQ[m][qId] || !(qId2 == -1 || _pCircuit->sQubitRegion(qId2).find(_curSol.vvProQ2PhyQ[m][qId]) != _pCircuit->sQubitRegion(qId2).end())){
    if(type){
        while( neighborPhyQId == _curSol.vvProQ2PhyQ[m][qId] ){
            auto it = sQubitRegion.cbegin();
            unsigned_t random = rand() % sQubitRegion.size();
            // cerr << "hi sQubitRegion.size(): " << sQubitRegion.size() << ", random: " << random << endl;
            advance(it, random);
            neighborPhyQId = *it;
            qId2 = _curSol.vvPhyQ2ProQ[m][neighborPhyQId];
        }
    }
    else{
        while( neighborPhyQId == _curSol.vvProQ2PhyQ[m][qId] ){
            neighborPhyQId = rand() % _pDevice->nQubit();
        }
    }
    assert(neighborPhyQId >= 0);
    assert(neighborPhyQId < _pDevice->nQubit());
    get<2>(_saParam.movement) = neighborPhyQId;
}

void Placer::decideQubitInitialMapping(){
    my_graph graph(_pCircuit->nProgramQubit() + _pDevice->nQubit());
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        unordered_set<int_t>& sQubitRegion = _pCircuit->sQubitRegion(i);
        for(unsigned_t j : sQubitRegion){
            add_edge(i, _pCircuit->nProgramQubit() + j, EdgeProperty(1), graph);    
        }
    }
    vector<V> mate(num_vertices(graph));
    vector<unsigned_t> vFinerQ2CoarserQ(_pDevice->nQubit(), _pDevice->nQubit());
    // use boost library to construct maximum card matching 
    // see: https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/maximum_matching.html
    edmonds_maximum_cardinality_matching(graph, &mate[0]);
    
    unsigned_t pairCount = 0;
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            _pCircuit->setInitialMapping(v, mate[v]-_pCircuit->nProgramQubit());
           ++pairCount;
        }
    }
    _firstInit = false;
    assert(pairCount == _pCircuit->nProgramQubit());
}

void Placer::initSol(){  
    /* Random generate a placement solution */
    // call get cost here
    _curSol.vvProQ2PhyQ.clear();
    _curSol.vvProQ2PhyQ.resize(_vvWeightedEdge.size(), vector<unsigned_t>(_pCircuit->nProgramQubit(), 0));
    _curSol.vvPhyQ2ProQ.clear();
    _curSol.vvPhyQ2ProQ.resize(_vvWeightedEdge.size(), vector<int_t>(_pDevice->nQubit(), -1));
    unsigned_t phyQ;
    if(!_isGivenInitialMapping && _pCircuit->sQubitRegion(0).size() == 0){
        vector<unsigned_t> vPhyQubit(_pDevice->nQubit());
        for (unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
            vPhyQubit[i] = i;
        }
        random_shuffle( vPhyQubit.begin(), vPhyQubit.end() );
        for (unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
            _pCircuit->setInitialMapping(i, vPhyQubit[i]);
        }
    }
    for (unsigned_t j = 0; j < _pCircuit->nProgramQubit(); ++j){
        phyQ = _pCircuit->initialMapping(j);
        _curSol.vvProQ2PhyQ[0][j] = phyQ;
        _curSol.vvPhyQ2ProQ[0][phyQ] = j; 
    }
    for (unsigned_t i = 1; i < _vvWeightedEdge.size(); ++i){
        for (unsigned_t j = 0; j < _pCircuit->nProgramQubit(); ++j){
            _curSol.vvProQ2PhyQ[i][j] = _curSol.vvProQ2PhyQ[0][j]; 
            _curSol.vvPhyQ2ProQ[i][_curSol.vvProQ2PhyQ[0][j]] = _curSol.vvPhyQ2ProQ[0][_curSol.vvProQ2PhyQ[0][j]]; 
        }
    }

    // cerr << "print qubit mapping swap" << endl;
    // cerr << "proQ2phyQ: ";
    // for (unsigned_t i = 0; i < _curSol.vvProQ2PhyQ[0].size(); ++i){
    //     cerr << i << " -> " << _curSol.vvProQ2PhyQ[0][i] << endl;
    // }

    getCost();
    // cerr << "first init cost: " << _curSol.cost << endl;
    if(_curSol.cost < _optSol.cost){
        updateOptSol();
    }
    // cerr << "_pCircuit->nProgramQubit(): " << _pCircuit->nProgramQubit() << endl;
    // cerr << "_curSol.vvProQ2PhyQ[0].size(): " << _curSol.vvProQ2PhyQ[0].size() << endl;
    /* Initial Perturbation */
    initPerturb();
    // cerr << "initPerturb" << endl;
    /* Initialize */
    _saParam.uphill_sum = 0;
    _saParam.uphill_avg_cnt = 0;
    _saParam.delta_sum = 0;
    _saParam.delta_cost_cnt = 0;
}

void Placer::initPerturb(){
    double_t uphillsum = 0;
    double_t sumCost = 0;
    int_t uphillcnt = 0;
    for(unsigned_t i = 0; i < _saParam.init_perturb_num; ++i){
        getMovement();
        makeMovement(); // make movement and calculate cost difference
        _curSol.cost += _saParam.delta;
        if(_curSol.cost < _optSol.cost){
            updateOptSol();
        }
        if(_saParam.delta > 0){
            uphillsum += _saParam.delta;
            ++uphillcnt;
        }
        sumCost += _curSol.cost;
    }
    // cerr << "uphillsum: " << uphillsum << endl;
    // cerr << "uphillcnt: " << uphillcnt << endl;
    _saParam.t1 = ((double_t)uphillsum / (double_t)uphillcnt) / ((-1)*log(_saParam.p));
    _saParam.t = _saParam.t1;
}

void Placer::updateT(){
    // update T based on fast SA
    if(_saParam.n <= _saParam.k){
        _saParam.t = (_saParam.t1 * abs(_saParam.delta_sum / _saParam.delta_cost_cnt) / _saParam.n / _saParam.c);
    }
    else if(_saParam.n > _saParam.k){
        _saParam.t = (_saParam.t1 * abs(_saParam.delta_sum / _saParam.delta_cost_cnt) / _saParam.n);
    }
    else{
        assert(false);
    }
}

bool Placer::acceptWorseSol(){
    // calculate the probability of accepting a worse solution
    return ((double_t)(rand() % _saParam.random_range) / _saParam.random_range) <= exp(-((double_t)_saParam.delta)/(_saParam.t));
}

void Placer::updateOptSol(){
    // update optimal solution to current solution
    double_t previousCost = _optSol.cost;
    _optSol = Solution(_curSol);
    // cout << "previousCost: " << previousCost << endl;
    // cout << "_optSol.cost " << _optSol.cost << endl;
    assert(previousCost >= _optSol.cost && _optSol.cost >= 0);
}

void Placer::recover(){
    // recover the current action, including cost array and qubit mapping
    unsigned_t m = get<0>(_saParam.movement);
    unsigned_t phyId0 = get<1>(_saParam.movement);
    unsigned_t phyId1 = get<2>(_saParam.movement);
    swapQubit(m, phyId0, phyId1);
    // recover cost
    _vCostForEachTime[m] = _oriCostTime;
    if (_curSol.vvPhyQ2ProQ[m][phyId0] > -1){
        _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId0]] = _oriCostQubit1;
    }
    if (_curSol.vvPhyQ2ProQ[m][phyId1] > -1){
        _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId1]] = _oriCostQubit0;
    }
}

void Placer::makeMovement(){
    // make movement, record previous cost, calculate cost difference
    unsigned_t m = get<0>(_saParam.movement);
    unsigned_t phyId0 = get<1>(_saParam.movement);
    unsigned_t phyId1 = get<2>(_saParam.movement);
    // record previous cost
    unsigned_t newCostForQubit = 0;
    // make the movement
    swapQubit(m, phyId0, phyId1);
    // update cost
    _oriCostQubit0 = 0;
    _oriCostQubit1 = 0;
    if (_curSol.vvPhyQ2ProQ[m][phyId0] > -1){
        // cout << "cal cost for qubit 0, id :" << _curSol.vvPhyQ2ProQ[m][phyId0] << endl;
        _oriCostQubit0 = _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId0]];
        calCostForQubit(_curSol.vvPhyQ2ProQ[m][phyId0]);
        newCostForQubit += _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId0]];
    }
    if (_curSol.vvPhyQ2ProQ[m][phyId1] > -1){
        // cout << "cal cost for qubit 0, id :" << _curSol.vvPhyQ2ProQ[m][phyId1] << endl;
        _oriCostQubit1 = _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId1]];
        calCostForQubit(_curSol.vvPhyQ2ProQ[m][phyId1]);
        newCostForQubit += _vCostForEachQubit[_curSol.vvPhyQ2ProQ[m][phyId1]];
    }
    _oriCostTime = _vCostForEachTime[m];
    calCostForTime(m);
    // calculate cost difference
    // cerr << "_saParam.delta: " << _saParam.delta << endl;
     // calculate cost difference
    _saParam.delta = ((double_t)newCostForQubit - (double_t)_oriCostQubit0 - (double_t)_oriCostQubit1)*_saParam.verticalCostWeight + (double_t)(_vCostForEachTime[m] - _oriCostTime);
    // cout << "swap phy qubit " << phyId0 << ", " << phyId1 << endl;
    // cout << "_oriCostTime: " << _oriCostTime << endl;
    // cout << "_vCostForEachTime[m]: " << _vCostForEachTime[m] << endl;
    // cout << "_saParam.delta: " << _saParam.delta << endl;
    // getchar();
}

void Placer::swapQubit(unsigned_t m, unsigned_t phyId0, unsigned_t phyId1, bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    // swap the mapping of phyId0 and phyId1
    int_t proId0 = sol.vvPhyQ2ProQ[m][phyId0];
    int_t proId1 = sol.vvPhyQ2ProQ[m][phyId1];
    // cerr << "print qubit mapping" << endl;
    // cerr << "phyQ2proQ: ";
    // for (unsigned_t i = 0; i < _curSol.vvPhyQ2ProQ[m].size(); ++i){
    //     cerr << i << " -> " << _curSol.vvPhyQ2ProQ[m][i] << endl;
    // }
    // cerr << "proQ2phyQ: ";
    // for (unsigned_t i = 0; i < _curSol.vvProQ2PhyQ[m].size(); ++i){
    //     cerr << i << " -> " << _curSol.vvProQ2PhyQ[m][i] << endl;
    // }
    // cerr << "swap physical qubits " << phyId0 << "(pro " << proId0 << ") and " << phyId1 << "(pro " << proId1 << ") in slice " << m << endl;
    assert(proId0 < 0 || phyId0 == sol.vvProQ2PhyQ[m][proId0]);
    assert(proId1 < 0 || phyId1 == sol.vvProQ2PhyQ[m][proId1]);
    if (0 <= proId0){
        sol.vvProQ2PhyQ[m][proId0] = phyId1;
    }
    if (0 <= proId1){
        sol.vvProQ2PhyQ[m][proId1] = phyId0;
    }
    sol.vvPhyQ2ProQ[m][phyId1] = proId0;
    sol.vvPhyQ2ProQ[m][phyId0] = proId1;
    // cerr << "print qubit mapping after swap" << endl;
    // cerr << "phyQ2proQ: ";
    // for (unsigned_t i = 0; i < _curSol.vvPhyQ2ProQ[m].size(); ++i){
    //     cerr << i << " -> " << _curSol.vvPhyQ2ProQ[m][i] << endl;
    // }
    // cerr << "proQ2phyQ: ";
    // for (unsigned_t i = 0; i < _curSol.vvProQ2PhyQ[m].size(); ++i){
    //     cerr << i << " -> " << _curSol.vvProQ2PhyQ[m][i] << endl;
    // }
}

void Placer::setMapping(){
    // set solution based on SA results
    // unsigned_t qId;
    // for(unsigned_t i = 0; i < _optSol.vvProQ2PhyQ.size(); ++i){
    //     fprintf(stdout, "           - Time %d\n", i);
        // cerr << "_optSol.vvProQ2PhyQ[i].size(): " << _optSol.vvProQ2PhyQ[i].size() << endl;
    for(unsigned_t i = 0; i < _optSol.vvProQ2PhyQ[0].size(); ++i){
        _pCircuit->setInitialMapping(i, _optSol.vvProQ2PhyQ[0][i]);
    }
    if(_verbose > 0){
        for(unsigned_t i = 0; i < _optSol.vvProQ2PhyQ[0].size(); ++i){
            fprintf(stdout, "              - Program Qubit %d is mapped to physical qubit %d\n", i, _optSol.vvProQ2PhyQ[0][i]);
        }
    }
}

void Placer::runRefinement(){
    unsigned_t m, i, j;
    int_t proId0, proId1;
    double_t costBeforeRefine = _optSol.cost, newCostForQubit;
    getCost(false);
    // cout << "costBeforeRefine: " << costBeforeRefine << ", _optSol.cost: " << _optSol.cost << endl;
    assert(fabs(costBeforeRefine - _optSol.cost) < 0.001);
    for (m = 0; m < _vvWeightedEdge.size(); ++m){
        for (i = 0; i < _pDevice->nQubit(); ++i){
            for (j = i; j < _pDevice->nQubit(); ++j){
                if(_optSol.vvPhyQ2ProQ[m][i] == -1 && _optSol.vvPhyQ2ProQ[m][j] == -1){
                    continue;
                }
                // record previous cost
                newCostForQubit = 0;
                _oriCostQubit0 = 0;
                _oriCostQubit1 = 0;
                _oriCostTime = _vCostForEachTime[m];
                // make the movement
                // swap qubit
                swapQubit(m, i, j, false);
                // update cost
                calCostForTime(m, false);
                if (_optSol.vvPhyQ2ProQ[m][i] > -1){
                    _oriCostQubit0 = _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][i]];
                    calCostForQubit(_optSol.vvPhyQ2ProQ[m][i], false);
                    newCostForQubit += _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][i]];
                    // cerr << "_oriCostQubit0: " << _oriCostQubit0 << ", _vCostForEachQubit[" << _optSol.vvPhyQ2ProQ[m][i] << "]: " << _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][i]] << endl;
                }
                if (_optSol.vvPhyQ2ProQ[m][j] > -1){
                    _oriCostQubit1 = _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][j]];
                    calCostForQubit(_optSol.vvPhyQ2ProQ[m][j], false);
                    newCostForQubit += _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][j]];
                    // cerr << "_oriCostQubit1: " << _oriCostQubit1 << ", _vCostForEachQubit[" << _optSol.vvPhyQ2ProQ[m][j] << "]: " << _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][j]] << endl;
                }
                // newCostForQubit += _vCostForEachQubit[i];
                // _oriCostQubit0 = _vCostForEachQubit[j];
                // calCostForQubit(j, false);
                // newCostForQubit += _vCostForEachQubit[j];
                _saParam.delta = ((double_t)newCostForQubit - (double_t)_oriCostQubit0 - (double_t)_oriCostQubit1)*_saParam.verticalCostWeight + (double_t)(_vCostForEachTime[m] - _oriCostTime);

                if (_saParam.delta > 0){
                    swapQubit(m, i, j, false);
                    // recover cost
                    _vCostForEachTime[m] = _oriCostTime;
                    if (_optSol.vvPhyQ2ProQ[m][i] > -1){
                        _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][i]] = _oriCostQubit1;
                        // cerr << "_oriCostQubit0: " << _oriCostQubit0 << ", _vCostForEachQubit[" << _optSol.vvPhyQ2ProQ[m][i] << "]: " << _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][i]] << endl;
                    }
                    if (_optSol.vvPhyQ2ProQ[m][j] > -1){
                        _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][j]] = _oriCostQubit0;
                        // cerr << "_oriCostQubit1: " << _oriCostQubit1 << ", _vCostForEachQubit[" << _optSol.vvPhyQ2ProQ[m][j] << "]: " << _vCostForEachQubit[_optSol.vvPhyQ2ProQ[m][j]] << endl;
                    }
                    // _vCostForEachQubit[i] = _oriCostQubit0;
                    // _vCostForEachQubit[j] = _oriCostQubit1;
                }
                else{
                    _optSol.cost += _saParam.delta;
                    // cerr << "_saParam.delta: " << _saParam.delta << endl;
                    // cerr << "_optSol.cost: " << _optSol.cost << endl;
                    // getCost(false);
                    // cerr << "after getCost, _optSol.cost: " << _optSol.cost << endl;
                }
            }
        }
        // costBeforeRefine = _optSol.cost;
        // getCost(false);
        // cerr << "costBeforeRefine: " << costBeforeRefine << endl;
        // cerr << "_optSol.cost: " << _optSol.cost << endl;       
    }

    assert(costBeforeRefine > _optSol.cost || fabs(costBeforeRefine - _optSol.cost) < 0.001 );
}

MOLSQ_NAMESPACE_CPP_END

/***********************************************************************
  File        [ placer.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ placer ]
  Synopsis    [ placer class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 5, Sep., 2023 ]
***********************************************************************/
#include "placer/placer_dev.hpp"

MOLSQ_NAMESPACE_CPP_START

void Placer_dev::run(Circuit& cir, Device& device, bool givenInitialMapping, bool allowOverlapping){
    if(_verbose > 0)
        fprintf(stdout, "[INFO] SA-Based Placer_dev: Start\n");
    _pCircuit = &cir;
    _pDevice = &device;
    _timer.start(TimeUsage::FULL);
    _isGivenInitialMapping = givenInitialMapping;
    // if(!_isGivenInitialMapping && _pCircuit->sQubitRegion(0).size() > 0){
    //     fprintf(stdout, "[Info] SA-Based Placer_dev: Find initial mapping\n");
    //     decideQubitInitialMapping();
    //     return;
    // }
    if(!_saParam.is_all_commute){
        constructGateWeight();
    }
    if(_verbose > 0)
        fprintf(stdout, "[INFO] SA-Based Placer_dev: Start SA: Run 1\n");
    _randomSol = false;
    _saParam.allow_overlapping = false;
    _optSol = Solution();
    runSA();
    // if the optSol has overlapping, recalculate the cost using a super large weight
    // so that the solution without overlapping in other stage can be selected.
    // Solution optSol = _optSol;
    if(_verbose > 0)
        fprintf(stdout, "                            Optimal cost: %4f\n", _optSol.cost);
    // fprintf(stdout, "[INFO] SA-Based Placer_dev: Start SA: Run 2\n");
    // _timer.restart(300);
    // _saParam.allow_overlapping = allowOverlapping;
    // _optSol = Solution();
    // runSA();
    // if(_saParam.allow_overlapping){
    //     _saParam.overlapping_weight = 10*_pCircuit->nGate();
    //     getCost(false);
    // }
    // fprintf(stdout, "                            Optimal cost of run 2: %4f\n", _optSol.cost);
    // if(((optSol.cost - _optSol.cost) > 0.00001) && checkInjectivity()){
    //     optSol = _optSol;
    // }
    // if(!_timer.isTimeout()){
    //     fprintf(stdout, "[INFO] SA-Based Placer_dev: Start SA: Run 3\n");
    //     _randomSol = true;
    //     _isGivenInitialMapping = false;
    //     _saParam.allow_overlapping = true;
    //     _optSol = Solution();
    //     runSA();
    //     _saParam.overlapping_weight = 10*_pCircuit->nGate();
    //     getCost(false);
    //     fprintf(stdout, "                            Optimal cost of run 3: %4f\n", _optSol.cost);
    //     if(((optSol.cost - _optSol.cost) > 0.00001) && checkInjectivity()){
    //         optSol = _optSol;
    //     }
    // }
    // if(!_timer.isTimeout()){
    //     fprintf(stdout, "[INFO] SA-Based Placer_dev: Start SA: Run 4\n");
    //     _isGivenInitialMapping = false;
    //     _saParam.allow_overlapping = false;
    //     _optSol = Solution();
    //     runSA();
    //     fprintf(stdout, "                           Optimal cost of run 4: %4f\n", _optSol.cost);
    // }
    // if((_optSol.cost - optSol.cost) > 0.00001){
    //     _optSol = optSol;
    // }
    // fprintf(stdout, "[INFO] SA-Based Placer_dev: Final cost: %4f\n", _optSol.cost);
    if(_verbose > 0)
        fprintf(stdout, "[INFO] SA-Based Placer_dev: Set initial mapping\n");
    setMapping();
    if(_verbose > 0)
        fprintf(stdout, "[INFO] SA-Based Placer_dev: Finish\n");
}

void Placer_dev::runSA() {
    unsigned_t iter = 0, noImproveIter = 0;
    unsigned_t numProcess = 0;
    double_t newCost = 0.0, lastCost = 0;
    bool firstTime = true;
    if(_pCircuit->nProgramQubit() > 85){
        _saParam.overlapping_weight = 0.001;
    }
    else{
        _saParam.overlapping_weight = 1;
    }
    _saParam.t = 100000.0;
    _saParam.uphill_avg_cnt = 0;;
    _saParam.uphill_sum = 0;
    _saParam.delta_cost_cnt = 0;
    _saParam.delta_sum = 0;
    _saParam.delta = 0;
    while(!_timer.isTimeout() && noImproveIter < 3){
        initSol();
        // if(_saParam.allow_overlapping && noImproveIter < 1){
        //     updateOverlappingWeight();
        //     getCost(false);
        // }
        if (firstTime){
            if(_verbose > 0)
                fprintf(stdout, "[INFO] SA-Based Placer_dev: Init cost %4f\n", _optSol.cost);
            firstTime = false;
        }
        // else{
        //     _saParam.allow_overlapping = false;
        // }
        
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
                    if (_optSol.cost - _curSol.cost > 1e-9){
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
            if(_saParam.allow_overlapping && noImproveIter < 1){
                // cerr << "_saParam.overlapping_weight: " << _saParam.overlapping_weight << endl;
                updateOverlappingWeight();
                getCost(false);
            }
            if(_saParam.n > _saParam.iter_limit) break;
        }
        if(_saParam.find_opt_sol) break;
        if(_verbose > 1){
            fprintf(stdout, "[INFO] SA-Based Placer_dev: Iter %d\n", iter);
            fprintf(stdout, "                        Curret cost: %4f\n", _optSol.cost);
            fprintf(stdout, "                        Refinement \n");
            fprintf(stdout, "                        Final cost: %4f\n", _optSol.cost);
        }
        runRefinement();
        if(fabs(lastCost - _optSol.cost) < 0.001){
            ++noImproveIter;
        }
        else{
            noImproveIter = 0;
        }
        lastCost = _optSol.cost;
        if(noImproveIter > 1){
            // cerr << "_saParam.overlapping_weight: " << _saParam.overlapping_weight << endl;
            _saParam.overlapping_weight = 2*_pCircuit->nGate();
            getCost(false);
        }
    }

    return;
}


void Placer_dev::constructGateWeight(){
    // _vGateWeight.clear();
    // _vGateWeight.resize(_pCircuit->nGate(), 0);
    vector<unsigned_t> vQubitLevel(_pCircuit->nProgramQubit(), 0);
    unsigned_t q0, q1, gLevel;
    // for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
    //     Gate & gate = _pCircuit->gate(i);
    //     if(gate.nTargetQubit() > 1){
    //         q0 = gate.targetProgramQubit(0);
    //         q1 = gate.targetProgramQubit(1);
    //         gLevel = max(vQubitLevel[q0], vQubitLevel[q1]);
    //         ++vQubitLevel[q0];
    //         ++vQubitLevel[q1];
    //         _vGateWeight[i] = pow(_saParam.gate_weight_base, gLevel/_saParam.decay_level);
    //     }
    // }

    vector<int_t> vQubitLastQubit(_pCircuit->nProgramQubit(), -1);
    vector<vector<double_t>> vQubitWeight(_pCircuit->nProgramQubit(), vector<double_t>(_pCircuit->nProgramQubit(), 0));
    double_t weight;
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        if(gate.nTargetQubit() > 1){
            if (gate.targetProgramQubit(0) < gate.targetProgramQubit(1)){
                q0 = gate.targetProgramQubit(0);
                q1 = gate.targetProgramQubit(1);
            }
            else{
                q1 = gate.targetProgramQubit(0);
                q0 = gate.targetProgramQubit(1);
            }
            gLevel = max(vQubitLevel[q0], vQubitLevel[q1]);
            ++vQubitLevel[q0];
            ++vQubitLevel[q1];
            weight = pow(_saParam.gate_weight_base, gLevel/_saParam.decay_level);
            vQubitWeight[q0][q1] += weight;
            if(_considerOneHopCost){
                if(vQubitLastQubit[q0] > -1){
                    if(vQubitLastQubit[q0] < q1){
                        vQubitWeight[vQubitLastQubit[q0]][q1] += weight * _saParam.weight_for_neighbor_qubit;
                    }
                    else{
                        vQubitWeight[q1][vQubitLastQubit[q0]] += weight * _saParam.weight_for_neighbor_qubit;

                    }
                }
                if(vQubitLastQubit[q1] > -1){
                    if(vQubitLastQubit[q1] < q0){
                        vQubitWeight[vQubitLastQubit[q1]][q0] += weight * _saParam.weight_for_neighbor_qubit;
                    }
                    else{
                        vQubitWeight[q0][vQubitLastQubit[q1]] += weight * _saParam.weight_for_neighbor_qubit;
                    }
                }
            }
            vQubitLastQubit[q0] = q1;
            vQubitLastQubit[q1] = q0;
        }
    }
    _vQubitWeight.clear();
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        for(unsigned_t j = i+1; j < _pCircuit->nProgramQubit(); ++j){
            if(vQubitWeight[i][j] > 0.0000001){
                _vQubitWeight.emplace_back(make_pair(i,j), vQubitWeight[i][j]);
            }
        }    
    }
}

void Placer_dev::getCost(bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    // calculate horizontal distance cost (between qubits)
    calCostForTime(isCurSol);
    sol.cost = _costForTime;
}

void Placer_dev::calCostForTime(bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    _costForTime = 0;
    unsigned_t qId0, qId1, i;
    if(_saParam.is_all_commute){
        for (i = 0; i < _pCircuit->nGate(); ++i){
            Gate & gate = _pCircuit->gate(i);
            if(gate.nTargetQubit() == 2){
                qId0 = sol.vProQ2PhyQ[gate.targetProgramQubit(0)];
                qId1 = sol.vProQ2PhyQ[gate.targetProgramQubit(1)];
                _costForTime += (_pDevice->getDistance(qId0, qId1));
            }
        }
    }
    else{
        for(i = 0; i < _vQubitWeight.size(); ++i){
            qId0 = sol.vProQ2PhyQ[_vQubitWeight[i].first.first];
            qId1 = sol.vProQ2PhyQ[_vQubitWeight[i].first.second];
            if(!_saParam.is_all_commute){
                _costForTime += ((_pDevice->getDistance(qId0, qId1)) * _vQubitWeight[i].second);
            }
            else{
                _costForTime += (_pDevice->getDistance(qId0, qId1));
            }
        }
    }
}

void Placer_dev::decideQubitInitialMapping(){
    my_graph graph(_pCircuit->nProgramQubit() + _pDevice->nQubit());
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        unordered_set<int_t>& sQubitRegion = _pCircuit->sQubitRegion(i);
        for(unsigned_t j : sQubitRegion){
            add_edge(i, _pCircuit->nProgramQubit() + j, EdgeProperty(1), graph);    
        }
    }
    vector<V> mate(num_vertices(graph));
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

void Placer_dev::initSol(){  
    /* Random generate a placement solution */
    // call get cost here
    _curSol.vProQ2PhyQ.clear();
    _curSol.vProQ2PhyQ.resize(_pCircuit->nProgramQubit(), 0);
    _curSol.vvPhyQ2ProQ.clear();
    _curSol.vvPhyQ2ProQ.resize(_pDevice->nQubit(), vector<int_t>());
    unsigned_t phyQ;
    if((!_isGivenInitialMapping) || _randomSol){
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
        _curSol.vProQ2PhyQ[j] = phyQ;
        _curSol.vvPhyQ2ProQ[phyQ].emplace_back(j); 
    }


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

void Placer_dev::initPerturb(){
    double_t uphillsum = 0;
    double_t sumCost = 0;
    int_t uphillcnt = 0;
    // cerr << "in initPerturb" << endl;
    for(unsigned_t i = 0; i < _saParam.init_perturb_num; ++i){
        // if ( i >= 40 && i < 50){
        //     cerr << "i: " << i << endl;
        //     cerr << "get movement" << endl;
        // }
        getMovement();
        // if ( i >= 40 && i < 50){
        //     cerr << "make movement" << endl;
        // }
        makeMovement(); // make movement and calculate cost difference
        // if ( i >= 40 && i < 50){
        //     cerr << "finish movement" << endl;
        // }
        // cerr << "new qubit mapping" << endl;
        // cerr << "proQ2phyQ: ";
        // for (unsigned_t i = 0; i < _curSol.vProQ2PhyQ.size(); ++i){
        //     cerr << i << " -> " << _curSol.vProQ2PhyQ[i] << endl;
        // }
        // cerr << "phyQ2pro: ";
        // for (unsigned_t i = 0; i < _curSol.vvPhyQ2ProQ.size(); ++i){
        //     for (unsigned_t j = 0; j < _curSol.vvPhyQ2ProQ[i].size(); ++i){
        //         cerr << i << " -> " << _curSol.vvPhyQ2ProQ[i][j] << endl;
        //     }
        // }
        _curSol.cost += _saParam.delta;
        if(_optSol.cost - _curSol.cost > 0.00001 ){
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

void Placer_dev::updateT(){
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

void Placer_dev::updateOverlappingWeight(){
    _saParam.overlapping_weight = 10.0*(exp(1.0/(_saParam.t)) - 0.9);
}

bool Placer_dev::acceptWorseSol(){
    // calculate the probability of accepting a worse solution
    return ((double_t)(rand() % _saParam.random_range) / _saParam.random_range) <= exp(-((double_t)_saParam.delta)/(_saParam.t));
}

void Placer_dev::updateOptSol(){
    // update optimal solution to current solution
    double_t previousCost = _optSol.cost;
    _optSol = Solution(_curSol);
    // cout << "previousCost: " << previousCost << endl;
    // cout << "_optSol.cost " << _optSol.cost << endl;
    assert(previousCost >= _optSol.cost && _optSol.cost >= 0);
}

void Placer_dev::recover(){
    // recover the current action, including cost array and qubit mapping
    bool type = get<0>(_saParam.movement);
    // make the movement
    // swap two qubits 
    if(type){
        moveProQtoPhyQ(get<1>(_saParam.movement), _oldPhyQforTargetProQ);
    }
    else{
        swapQubit(get<1>(_saParam.movement), get<2>(_saParam.movement));
    }
    _costForTime = _oriCostTime;
}


void Placer_dev::getMovement(){
    bool moveWithinRegion = (((double_t)(rand() % _saParam.random_range) / _saParam.random_range) <= _saParam.prob_move_within_region);
    unsigned_t qId = rand() % _pCircuit->nProgramQubit(); // randomly choose a program qubit to move
    unsigned_t neighborPhyQId = _curSol.vProQ2PhyQ[qId]; 
    assert(neighborPhyQId >= 0);
    assert(neighborPhyQId < _pDevice->nQubit());
    unordered_set<int_t>& sQubitRegion = _pCircuit->sQubitRegion(qId);
    if(sQubitRegion.size() <= 1 || !_restrictMovement){
        moveWithinRegion = false;
    }
    if(moveWithinRegion){
        while( neighborPhyQId == _curSol.vProQ2PhyQ[qId] ){
            auto it = sQubitRegion.cbegin();
            unsigned_t random = rand() % sQubitRegion.size();
            // cerr << "hi sQubitRegion.size(): " << sQubitRegion.size() << ", random: " << random << endl;
            advance(it, random);
            neighborPhyQId = *it;
        }
    }
    else{
        while( neighborPhyQId == _curSol.vProQ2PhyQ[qId] ){
            neighborPhyQId = rand() % _pDevice->nQubit();
        }
    }
    assert(neighborPhyQId >= 0);
    assert(neighborPhyQId < _pDevice->nQubit());

    get<1>(_saParam.movement) = qId;
    if(_curSol.vvPhyQ2ProQ[neighborPhyQId].size() > 0){
        if(_saParam.allow_overlapping){
            get<0>(_saParam.movement) = rand() % 2;
            if(get<0>(_saParam.movement) == 0){
                unsigned_t randIdx = rand() % _curSol.vvPhyQ2ProQ[neighborPhyQId].size();
                get<2>(_saParam.movement) = _curSol.vvPhyQ2ProQ[neighborPhyQId][randIdx];    
            }
            else{
                get<0>(_saParam.movement) = 1;
                get<2>(_saParam.movement) = neighborPhyQId;    
            }
        }
        // nooverlapping is allowed, _curSol.vvPhyQ2ProQ[neighborPhyQId].size() == 1. neighborPhyQId is occupied by pro qubit. swap two qubits.
        else{
            get<0>(_saParam.movement) = 0;
            get<2>(_saParam.movement) = _curSol.vvPhyQ2ProQ[neighborPhyQId][0];    
        }
    }
    // no pro qubit is on neighborPhyQId
    else{
        get<0>(_saParam.movement) = 1;
        get<2>(_saParam.movement) = neighborPhyQId;
    }
    // if(get<0>(_saParam.movement) == 1){
    //     cerr << "move program qubit " << get<1>(_saParam.movement) << " to phy qubit " << get<2>(_saParam.movement) << endl;
    // }
    // else{
    //     cerr << "swap program qubit " << get<1>(_saParam.movement) << " and " << get<2>(_saParam.movement) << endl;
    // }
}

void Placer_dev::makeMovement(){
    // make movement, record previous cost, calculate cost difference
    bool type = get<0>(_saParam.movement);
    // make the movement
    // swap two qubits 
    if(type){
        moveProQtoPhyQ(get<1>(_saParam.movement), get<2>(_saParam.movement));
    }
    else{
        swapQubit(get<1>(_saParam.movement), get<2>(_saParam.movement));
    }
    // record previous cost
    _oriCostTime = _costForTime;
    // update cost
    calCostForTime();
    // calculate cost difference
    _saParam.delta = (double_t)_costForTime - (double_t)_oriCostTime;
    // cout << "_oriCostTime: " << _oriCostTime << endl;
    // cout << "_costForTime: " << _costForTime << endl;
    // cout << "_saParam.delta: " << _saParam.delta << endl;
    // getchar();
}


// move a program qubit to a new physical qubit
void Placer_dev::moveProQtoPhyQ(unsigned_t proId0, unsigned_t phyId1, bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    // move proId0 to phyId1
    _oldPhyQforTargetProQ = sol.vProQ2PhyQ[proId0]; // record the old position for recovering

    sol.vProQ2PhyQ[proId0] = phyId1;

    auto it = find(sol.vvPhyQ2ProQ[_oldPhyQforTargetProQ].begin(), sol.vvPhyQ2ProQ[_oldPhyQforTargetProQ].end(), proId0);
    assert(it != sol.vvPhyQ2ProQ[_oldPhyQforTargetProQ].end());
    sol.vvPhyQ2ProQ[_oldPhyQforTargetProQ].erase(it);

    sol.vvPhyQ2ProQ[phyId1].emplace_back(proId0);
}

// SWAP two program qubits
void Placer_dev::swapQubit(unsigned_t proId0, unsigned_t proId1, bool isCurSol){
    Solution & sol = isCurSol ? _curSol : _optSol;
    // swap the mapping of proId0 and proId1
    unsigned_t phyId0 = sol.vProQ2PhyQ[proId0];
    unsigned_t phyId1 = sol.vProQ2PhyQ[proId1];


    auto it = find(sol.vvPhyQ2ProQ[phyId0].begin(), sol.vvPhyQ2ProQ[phyId0].end(), proId0);
    assert(it != sol.vvPhyQ2ProQ[phyId0].end());
    sol.vvPhyQ2ProQ[phyId0].erase(it);
    sol.vvPhyQ2ProQ[phyId0].emplace_back(proId1);

    it = find(sol.vvPhyQ2ProQ[phyId1].begin(), sol.vvPhyQ2ProQ[phyId1].end(), proId1);
    assert(it != sol.vvPhyQ2ProQ[proId1].end());
    sol.vvPhyQ2ProQ[phyId1].erase(it);
    sol.vvPhyQ2ProQ[phyId1].emplace_back(proId0);

    sol.vProQ2PhyQ[proId0] = phyId1;
    sol.vProQ2PhyQ[proId1] = phyId0;
}

void Placer_dev::setMapping(){
    // set solution based on SA results
    for(unsigned_t i = 0; i < _optSol.vProQ2PhyQ.size(); ++i){
        _pCircuit->setInitialMapping(i, _optSol.vProQ2PhyQ[i]);
    }
    if(_verbose > 1){
        for(unsigned_t i = 0; i < _optSol.vProQ2PhyQ.size(); ++i){
            fprintf(stdout, "              - Program Qubit %d is mapped to physical qubit %d\n", i, _optSol.vProQ2PhyQ[i]);
        }
    }
}

void Placer_dev::runRefinement(){
    unsigned_t i, j;
    int_t proId0, proId1;
    double_t costBeforeRefine = _optSol.cost, newCostForQubit;
    getCost(false);
    // cout << "costBeforeRefine: " << costBeforeRefine << ", _optSol.cost: " << _optSol.cost << endl;
    assert(fabs(costBeforeRefine - _optSol.cost) < 0.001);
    unsigned_t type, pro0, pro1;
    for (i = 0; i < _pDevice->nQubit(); ++i){
        for (j = i; j < _pDevice->nQubit(); ++j){
            if(_optSol.vvPhyQ2ProQ[i].size() == 0 && _optSol.vvPhyQ2ProQ[j].size() == 0){
                continue;
            }
            // record previous cost
            newCostForQubit = 0;
            _oriCostTime = _costForTime;
            // make the movement
            // assume no overlapping in the final solution
            if(_optSol.vvPhyQ2ProQ[i].size() == 0){
                type = 0;
                pro0 = _optSol.vvPhyQ2ProQ[j][0];
                // cerr << "move pro " << pro0 << " to phy " << i << endl;
                moveProQtoPhyQ(_optSol.vvPhyQ2ProQ[j][0], i, false);
            }
            else if(_optSol.vvPhyQ2ProQ[j].size() == 0){
                type = 1;
                pro0 = _optSol.vvPhyQ2ProQ[i][0];
                // cerr << "move pro " << pro0 << " to phy " << j << endl;
                moveProQtoPhyQ(_optSol.vvPhyQ2ProQ[i][0], j, false);
            }
            // if physical qubit is occupied, swap qubits
            else{
                type = 2;
                pro0 = _optSol.vvPhyQ2ProQ[i][0];
                pro1 = _optSol.vvPhyQ2ProQ[j][0];
                swapQubit(_optSol.vvPhyQ2ProQ[i][0], _optSol.vvPhyQ2ProQ[j][0], false);
            }
            // update cost
            calCostForTime(false);
            _saParam.delta = (double_t)_costForTime - (double_t)_oriCostTime;

            if (_saParam.delta > 0){
                if(type == 0){
                    moveProQtoPhyQ(pro0, j, false);
                }
                else if(type == 1){
                    moveProQtoPhyQ(pro0, i, false);
                }
                // if physical qubit is occupied, swap qubits
                else{
                    swapQubit(pro0, pro1, false);
                }

                // recover cost
                _costForTime = _oriCostTime;
            }
            else{
                _optSol.cost += _saParam.delta;
                // cerr << "_saParam.delta: " << _saParam.delta << endl;
                // cerr << "_optSol.cost: " << _optSol.cost << endl;
                // getCost(false);
                // cerr << "after getCost, _optSol.cost: " << _optSol.cost << endl;
            }
        }
        // costBeforeRefine = _optSol.cost;
        // getCost(false);
        // cerr << "costBeforeRefine: " << costBeforeRefine << endl;
        // cerr << "_optSol.cost: " << _optSol.cost << endl;       
    }

    assert(costBeforeRefine > _optSol.cost || fabs(costBeforeRefine - _optSol.cost) < 0.001 );
}

bool Placer_dev::checkInjectivity(){
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        if(_optSol.vvPhyQ2ProQ[i].size() > 1){
            return false;
        }
    }
    return true;
}

MOLSQ_NAMESPACE_CPP_END

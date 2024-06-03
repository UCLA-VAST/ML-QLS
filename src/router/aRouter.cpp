/***********************************************************************
  File        [ aRouter.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ rotuer ]
  Synopsis    [ aRouter class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "router/aRouter.hpp"

MOLSQ_NAMESPACE_CPP_START

void aRouter::setQubitRegion(Circuit& cir){
    _vsQubitRegion.clear(); 
    _vsQubitRegion.resize(cir.nProgramQubit());
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        unordered_set<int_t> & sQubitRegion = cir.sQubitRegion(i);
        for(int_t j : sQubitRegion){
            _vsQubitRegion[i].insert(j);
        }
    }
}

bool aRouter::run(Circuit& cir, Device& device){
    if(cir.nProgramQubit() < 60 || (_aRouterParam.is_all_commute && cir.nProgramQubit() < 200)){
        _aRouterParam.node_limit = 800000;
    }
    if(!_aRouterParam.is_all_commute){
        _aRouterParam.gate_per_astar_base = cir.nGate();
        _aRouterParam.gate_dis_per_astar = 10000;
    }
    // _aRouterParam.node_limit = 100;
    // if(_aRouterParam.restrict_region){
    //     fprintf(stdout, "[Info] Qubit Region Info                              \n");
    //     fprintf(stdout, "       ------------------------------------------\n");
    //     unsigned_t i;
    //     for ( i = 0; i < _vsQubitRegion.size();  ++i ){
    //         fprintf(stdout, "        - Qubit %d: ", i);
    //         for(int_t q : _vsQubitRegion[i]){
    //             fprintf(stdout, "%d ", q);
    //         }
    //         fprintf(stdout, "\n");
    //     }
    // }
    _aRouterParam.gate_per_astar = _aRouterParam.gate_per_astar_base;
    fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter: Start\n");
    // fprintf(stdout, "[INFO] g_cost limit: %.4f\n",_aRouterParam.g_cost_limit);
    _pCircuit = &cir;
    _pCircuit->clearSwap();
    _pCircuit->setCircuitDepth(0);
    _pDevice = &device;
    _pvpGateDependency = cir.pvpGateDependency();
    _executedGateCount = 0;
    _vExecutedGate.clear();
    _vExecutedGate.resize(_pCircuit->nGate(), false);
    _vQubitHasGate.clear();
    _vQubitHasGate.resize(_pDevice->nQubit(), false);
    _vCurMapping.clear();
    _vCurMapping.resize(cir.nProgramQubit(), 0);
    _vCurSolMappingPro2Phy.clear();
    _vCurSolMappingPro2Phy.resize(cir.nProgramQubit(), 0);
    _vCurSolMappingPhy2Pro.clear();
    _vCurSolMappingPhy2Pro.resize(_pDevice->nQubit(), -1);
    if(_aRouterParam.restrict_region && _vsQubitRegion.size() != cir.nProgramQubit()){
        setQubitRegion(cir);
    }
    cir.resetQubitRegion();
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurMapping[i] = _pCircuit->initialMapping(i);
        // _pCircuit->addQubitRegion(i, _pCircuit->initialMapping(i));
        _pCircuit->addQubitRegion(i, _pCircuit->initialMapping(i));
        // cerr << "_pCircuit->initialMapping(i): " << _pCircuit->initialMapping(i) << endl;
    }
    _timer.start(TimeUsage::FULL);
    constructDependencyInfo();
    // start swap insertion
    #ifdef DEBUG
        cerr << "start astar" << endl;
    #endif
    // delete gates that can be executed under initial mapping
    firstPeeling();
    bool reduceKRun = false;
    fprintf(stdout, "[INFO] #Gate: %d\n", _pCircuit->nGate());
    while (_executedGateCount < _pCircuit->nGate()){
        // _vCurTargetGates: set to store current gate to route 
        // cerr << "_executedGateCount: " << _executedGateCount << endl;
        _timer.restart(_aRouterParam.timeout);
        fprintf(stdout, "[INFO] #Routed Gate: %d\n", _executedGateCount);
        // fflush(stdout);
        
        // for(unsigned_t i = 0; i < _vExecutedGate.size(); ++i){
        //     cerr << i << " : " << _vExecutedGate[i] << ", ";
        // }
        collectGatesForSearch();
        fprintf(stdout, "[INFO] Plan to route %d gate in this A* run\n", _sCurTargetGates.size());
        // if(_sCurTargetGates.size() < 1){
        //     _pDevice->printDevice();
        //     fprintf(stdout, "[Info] Qubit Region Info                              \n");
        //     fprintf(stdout, "       ------------------------------------------\n");
        //     unsigned_t i;
        //     for ( i = 0; i < _vsQubitRegion.size();  ++i ){
        //         fprintf(stdout, "        - Qubit %d: ", i);
        //         for(int_t q : _vsQubitRegion[i]){
        //             fprintf(stdout, "%d ", q);
        //         }
        //         fprintf(stdout, "\n");
        //     }
        //     Gate & gate = _pCircuit->gate(*_sCurTargetGates.begin());
        //     cerr << " " << *_sCurTargetGates.begin() << ": " << gate.targetProgramQubit(0) << ", " << gate.targetProgramQubit(1) << " ("<< _vCurMapping[gate.targetProgramQubit(0)] << "," << _vCurMapping[gate.targetProgramQubit(1)] << ")" << endl;
        //     getchar();
        // }
        // run astar
        aSolution* pFinalSol = run_astar();

        if(pFinalSol == nullptr){
            // recover circuit
            fprintf(stdout, "[INFO] Fail to find a better solution with SWAP less than %f\n", _aRouterParam.g_cost_limit);
            fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter: Finish\n");
            cleanAllNode();     
            return false;
            // for(unsigned_t i : _sCurTargetGates){
            //     for(unsigned_t j : _vsGateChild[i]){
            //         _vsGateUncollectedParent[j].insert(i);
            //     }
            // }
            // // recollect gates
            // if(_aRouterParam.gate_per_astar > 100){
            //     _aRouterParam.gate_per_astar /= 2;
            // }
            // else if(_aRouterParam.gate_per_astar > 20){
            //     _aRouterParam.gate_per_astar -= 10;
            // }
            // else if(_aRouterParam.gate_per_astar > 10){
            //     _aRouterParam.gate_per_astar -= 5;
            // }
            // else if(_aRouterParam.gate_per_astar > 1){
            //     --_aRouterParam.gate_per_astar;
            // }
            // if(_aRouterParam.gate_per_astar == 1){
            //     _aRouterParam.timeout = 86400;
            // }
            // if(!_aRouterParam.is_all_commute && _aRouterParam.dependency_graph_depth_limit > 0){
            //     _aRouterParam.dependency_graph_depth_limit -= 1;
            // }
            // _aRouterParam.node_limit += 10000;
            // cerr << "update _aRouterParam.gate_per_astar: " << _aRouterParam.gate_per_astar << endl;
        }
        else{
            // fprintf(stdout, "node num: %d\n", _sAllSolution.size());
            // update circuit
            #ifdef DEBUG
                cerr << "update circuit" << endl;
            #endif
            _aRouterParam.gate_per_astar = _aRouterParam.gate_per_astar_base;
            _aRouterParam.dependency_graph_depth_limit = _aRouterParam.dependency_graph_depth_limit_base;
            updateCircuit(pFinalSol);
        }
        cleanAllNode();       
    }
    fprintf(stdout, "[INFO] #Routed Gate: %d\n", _executedGateCount);
    _pCircuit->setCircuitDepth(_pCircuit->circuitDepth()+1);
    fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter: Finish\n");
    return true;
}

void aRouter::firstPeeling(){
    bool keepFindNewGate = true;
    bool canExecute = false;
    #ifdef DEBUG
        cerr << "In first peeling: remove gates that can be directly executed under the initial mapping" << endl;
    #endif
    while(keepFindNewGate){
        keepFindNewGate = false;
        for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
            if(_vsGateUncollectedParent[i].size() == 0){
                Gate & gate = _pCircuit->gate(i);
                canExecute = (gate.nTargetQubit() == 1) || _pDevice->isAdjacent(_pCircuit->initialMapping(gate.targetProgramQubit(0)), _pCircuit->initialMapping(gate.targetProgramQubit(1)));
                if(canExecute && !_vExecutedGate[i]){
                    gate.setExecutionTime(0);
                    for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
                        gate.setTargetPhysicalQubit(j, _pCircuit->initialMapping(gate.targetProgramQubit(j)));
                        _vQubitHasGate[_pCircuit->initialMapping(gate.targetProgramQubit(j))] = true;
                    }
                    #ifdef DEBUG
                        cerr << "gate " << i << " can be directly executed on program qubit " << gate.targetProgramQubit(0) << "("<< _pCircuit->initialMapping(gate.targetProgramQubit(0)) << ")";
                        if(gate.nTargetQubit() == 2){
                            cerr << ", " << gate.targetProgramQubit(1) << "("<< _pCircuit->initialMapping(gate.targetProgramQubit(1)) << ")" << endl; 
                        }
                        else{
                            cerr << endl;
                        }
                    #endif
                    _vExecutedGate[i] = true;
                    ++_executedGateCount;
                    keepFindNewGate = true;
                    for(unsigned_t j : _vsGateChild[i]){
                        _vsGateUncollectedParent[j].erase(i);
                    }
                }
            }
        }
    }
}

void aRouter::constructDependencyInfo(){
    _vsGateUncollectedParent.clear();
    _vsGateUncollectedParent.resize(_pCircuit->nGate(), unordered_set<unsigned_t>());
    _vsGateUnexecutedParent.clear();
    _vsGateUnexecutedParent.resize(_pCircuit->nGate(), unordered_set<unsigned_t>());
    _vsGateChild.clear();
    _vsGateChild.resize(_pCircuit->nGate(), unordered_set<unsigned_t>());
    if(_aRouterParam.is_all_commute){
        return;
    }
    for(pair<unsigned_t, unsigned_t>& pDependency : (*_pvpGateDependency)){
        _vsGateUncollectedParent[pDependency.second].insert(pDependency.first);
        _vsGateUnexecutedParent[pDependency.second].insert(pDependency.first);
        _vsGateChild[pDependency.first].insert(pDependency.second);
    }
    #ifdef DEBUG
        cout << "Print _vsGateUncollectedParent and _vsGateChild for gates" << endl;
        for (unsigned_t i = 0; i < _vsGateUncollectedParent.size(); ++i){
            cout << "Print _vsGateUncollectedParent for gate " << i << ": ";
            for (unsigned_t gId : _vsGateUncollectedParent[i]){
                cout << gId << " ";
            }
            cout << endl;
            cout << "Print _vsGateChild for gate " << i << ": ";
            for (unsigned_t gId : _vsGateChild[i]){
                cout << gId << " ";
            }
            cout << endl;
        }
    #endif
}

void aRouter::collectGatesForSearch(){
    // add gate to _vCurTargetGates

    _curDisSum = 0;
    _sCurTargetGates.clear();
    unsigned_t p0, p1, twoQubitGateCount = 0;
    vector<pair<double_t, unsigned_t>> vGateOrder;
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        if(!_vExecutedGate[i] && (_vsGateUncollectedParent[i].size() == 0) && (_sCurTargetGates.count(i) == 0)){
            Gate & gate = _pCircuit->gate(i);
            if(gate.nTargetQubit() > 1){
                p0 = _vCurMapping[gate.targetProgramQubit(0)];
                p1 = _vCurMapping[gate.targetProgramQubit(1)];
                vGateOrder.emplace_back(_pDevice->getDistance(p0, p1), i);
                // _curDisSum += _pDevice->getDistance(p0, p1);
            }
            else{
                vGateOrder.emplace_back(0, i);
            }
        }
    }
    sort(vGateOrder.begin(), vGateOrder.end());
    if(_aRouterParam.is_all_commute){
        for(auto & pairWeightGate : vGateOrder){
            _sCurTargetGates.emplace(pairWeightGate.second);
            Gate & gate = _pCircuit->gate(pairWeightGate.second);
            if(gate.nTargetQubit() > 1){
                p0 = _vCurMapping[gate.targetProgramQubit(0)];
                p1 = _vCurMapping[gate.targetProgramQubit(1)];
                _curDisSum += _pDevice->getDistance(p0, p1);
                ++twoQubitGateCount;
            }
            if(twoQubitGateCount > _aRouterParam.gate_per_astar || _curDisSum > _aRouterParam.gate_dis_per_astar){
                break;
            }
        }
    }
    else{
        unsigned_t beginCount = 0, threshold = 5;
        for(auto & pairWeightGate : vGateOrder){
            if(pairWeightGate.first < threshold || beginCount < 3){
                Gate & gate = _pCircuit->gate(pairWeightGate.second);
                collectCurTargetGatesFromGate(pairWeightGate.second, 0);
                ++beginCount;
                if(_sCurTargetGates.size() > _aRouterParam.gate_per_astar || _curDisSum > _aRouterParam.gate_dis_per_astar){
                    break;
                }
            }
        }
    }
    #ifdef DEBUG
        cerr << "_curDisSum: " << _curDisSum << endl;
        cerr << "_vCurTargetGates in the current run: ";
        for(unsigned_t i : _sCurTargetGates)
            cerr << i << " ";
        cerr << endl;
    #endif
}

void aRouter::collectCurTargetGatesFromGate(unsigned_t g, unsigned_t curLevel ){
    if(_sCurTargetGates.size() > _aRouterParam.gate_per_astar  || _curDisSum > _aRouterParam.gate_dis_per_astar || curLevel > _aRouterParam.dependency_graph_depth_limit){
        return;
    }
    _sCurTargetGates.emplace(g);
    Gate & gate = _pCircuit->gate(g);
    if(gate.nTargetQubit() > 1){
        unsigned_t p0 = _vCurMapping[gate.targetProgramQubit(0)];
        unsigned_t p1 = _vCurMapping[gate.targetProgramQubit(1)];
        _curDisSum += _pDevice->getDistance(p0, p1);
    }
    for(unsigned_t i : _vsGateChild[g]){
        _vsGateUncollectedParent[i].erase(g);
        if(!_vExecutedGate[i] && _vsGateUncollectedParent[i].size() == 0 && (_sCurTargetGates.count(i) == 0)){
            // collectCurTargetGatesFromGate(i, curLevel+1);
            if(_pCircuit->gate(i).nTargetQubit() > 1){
                collectCurTargetGatesFromGate(i, curLevel+1);
            }
            else{
                collectCurTargetGatesFromGate(i, curLevel);
            }
        }
    }
}


aRouter::aSolution* aRouter::run_astar() {
    // Early return
    // aSolution& initSol = _vAllSolution[0];
    aSolution* pInitSol = new aSolution();
    _sAllSolution.insert(pInitSol);
    #ifdef DEBUG
        cerr << "start initialization" << endl;
    #endif
    initialize(pInitSol);
    if (pInitSol->sUnmappedGate.size() == 0 && pInitSol->sReadyMappedGate.size() == 0){
        return pInitSol;
    }
    // Initialize priority queue (frontier)
    _priorityQ.clear();
    _priorityQ.push(pInitSol);

    // Declare Node* pCurNode for the extracted min node
    aSolution* pCurNode = nullptr;

    while ( !_priorityQ.empty() ) {
        pCurNode = _priorityQ.top();
        if(_aRouterParam.gate_per_astar > 1){
            if(_timer.isTimeout()){
                fprintf(stdout, "[INFO] Timeout in aRouter\n");
                // if(_aRouterParam.is_all_commute)
                //     return nullptr;
                while(pCurNode->vExecutedGate.size() == 0 && pCurNode->pParent!=nullptr){
                    _priorityQ.pop();
                    pCurNode = _priorityQ.top();
                }
                return pCurNode;
            }
            else if(_sAllSolution.size() > _aRouterParam.node_limit){
                // if(checkMemUsage()){
                if(1){
                    fprintf(stdout, "[INFO] Too many nodes (>%d)\n", _aRouterParam.node_limit);
                    // if(_aRouterParam.is_all_commute)
                        // return nullptr;
                    while(pCurNode->vExecutedGate.size() == 0 && pCurNode->pParent!=nullptr){
                        _priorityQ.pop();
                        pCurNode = _priorityQ.top();
                    }
                    return pCurNode;
                }
                else{
                    _aRouterParam.node_limit += 10000;
                }
            }
        }
        // Extract min node
        _priorityQ.pop();
        assert(pCurNode->done == false);
        pCurNode->done = true; // mark the extracted node as DONE

        #ifdef DEBUG
            cerr << "====================================" << endl;
            cerr << "explore node" << endl;
            printSolution(pCurNode);
            getchar();
        #endif
        
        // Termination condition
        if (pCurNode->sUnmappedGate.size() == 0 && pCurNode->sReadyMappedGate.size() == 0) {
            return pCurNode;
        }
        updateCurMappingForSol(pCurNode);
        // construct adjacent nodes
        collectNeighbor(pCurNode);
        // Relax neighbors
        #ifdef DEBUG
            cerr << endl << "In astar: neighbor size: " << pCurNode->vAdjs.size() << endl;
        #endif
        for (aSolution* pSol : pCurNode->vAdjs) {
            // Skip adjNode which is "done"
            // assert(_sAllSolution.find(pSol) != _sAllSolution.end());
            if (pSol->done)
                continue;
            // update _vCurSolMappingPro2Phy and _vCurSolMappingPhy2Pro based on swapEdge
            swapQubit(pSol->swapEdge);
            // Initialize node if necessary

            // Calculate the actual cost from pSrcNode to pCurNode (i.e. g(n))
            const double_t tentative_g = pCurNode->g_cost + 1;
            // Relax if necessary
            if (((tentative_g + (double_t)_pCircuit->nSwapGate() + 1e-9) < _aRouterParam.g_cost_limit) && ((tentative_g + 1e-9) < pSol->g_cost)) {
                initNode(pSol);
                
                // Set g cost (i.e. current actual cost from src to this adj node)
                pSol->g_cost = tentative_g;

                const double heuristic_cost = calHeuristicCost(pSol); 
                // const double heuristic_cost = ((double)(Point3D::Mdist(adjPoint, tarPoint))) * _routerParam.heuristic_factor;
                pSol->f_cost = tentative_g + heuristic_cost; // g + h (A*)

                // Update pqueue for pAdjNode
                _priorityQ.push(pSol); // insert to pqueue if the node is first visited
            }
            else{
                _sAllSolution.erase(pSol);
            }
            // recover _vCurSolMappingPro2Phy and _vCurSolMappingPhy2Pro based on swapEdge
            swapQubit(pSol->swapEdge);
        }
    }

    return nullptr;
}

void aRouter::initialize(aSolution* pSol){
    // 1. map gates that can be directly executed under the initial mapping 
    // 2. add gates that can requrie SWAP gate to _sGateReadyToMap
    // 3. add gates that does not belong to 1 and 2 to _sGateToMap
    // 4. add gates to _vqQubit2UnmapGate
    // _vpAllSolution.clear();
    // _vpAllSolution.emplace_back(pInitSol);
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        _vCurSolMappingPhy2Pro[i] = -1;
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy[i] = _vCurMapping[i];
        _vCurSolMappingPhy2Pro[_vCurMapping[i]] = i;
    }
    for(unsigned_t i : _sCurTargetGates){
        pSol->sUnmappedGate.insert(i);
    }
    initNode(pSol);
    pSol->g_cost = 0;
    pSol->f_cost = calHeuristicCost(pSol);
    
    #ifdef DEBUG
        cerr << "initialize" << endl;
        printSolution(pSol);
    #endif
}

void aRouter::swapQubit(int_t eId){
    if(eId > -1){
        // cerr << "before swap " << eId << "(" << _pDevice->edge(eId).qubitId1() << "," << _pDevice->edge(eId).qubitId2() << ")" << endl;
        // cerr << "print _vCurSolMappingPro2Phy" << endl;
        // for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPro2Phy[i] << endl;
        // }
        // cerr << "print _vCurSolMappingPhy2Pro" << endl;
        // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPhy2Pro[i] << endl;
        // }
        Edge & edge = _pDevice->edge(eId);

        int_t proId0 = _vCurSolMappingPhy2Pro[edge.qubitId1()];
        int_t proId1 = _vCurSolMappingPhy2Pro[edge.qubitId2()];
        if (0 <= proId0){
            _vCurSolMappingPro2Phy[proId0] = edge.qubitId2();
        }
        if (0 <= proId1){
            _vCurSolMappingPro2Phy[proId1] = edge.qubitId1();
        }
        _vCurSolMappingPhy2Pro[edge.qubitId2()] = proId0;
        _vCurSolMappingPhy2Pro[edge.qubitId1()] = proId1;
        // cerr << "after swap" << endl;
        // cerr << "print _vCurSolMappingPro2Phy" << endl;
        // for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPro2Phy[i] << endl;
        // }
        // cerr << "print _vCurSolMappingPhy2Pro" << endl;
        // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPhy2Pro[i] << endl;
        // }
    }
}

void aRouter::initNode(aSolution* pSol){
    unsigned_t q0, q1;
    vector<unsigned> vGateToErase;
    bool keepFindNewGate = true;
    while (keepFindNewGate){
        #ifdef DEBUG
            cerr << "keep Find New Gate" << endl;
            // printSolution(pSol);
        #endif
        keepFindNewGate = false;
        vGateToErase.clear();
        for(unsigned_t gateId : pSol->sReadyMappedGate){
            Gate & gate = _pCircuit->gate(gateId);
            if(gate.nTargetQubit() > 1){
                q0 = gate.targetProgramQubit(0);
                q1 = gate.targetProgramQubit(1);
                // cerr << "gate " << gateId << endl;
                // cerr << "q0 " << q0 << endl;
                // cerr << "q1 " << q1 << endl;
                // cerr << "_vCurSolMappingPro2Phy[q0] " << _vCurSolMappingPro2Phy[q0] << endl;
                // cerr << "_vCurSolMappingPro2Phy[q1] " << _vCurSolMappingPro2Phy[q1] << endl;

                if( _pDevice->isAdjacent(_vCurSolMappingPro2Phy[q0], _vCurSolMappingPro2Phy[q1])){
                    keepFindNewGate = true;
                    vGateToErase.emplace_back(gateId);
                    #ifdef DEBUG
                        cerr << "In sReadyMappedGate: " << gateId << " can be executed on qubit " << q0 << "("<< _vCurSolMappingPro2Phy[q0] << ") " << q1 << "("<< _vCurSolMappingPro2Phy[q1] << ")" << endl;
                    #endif
                }
            }
            else{
                keepFindNewGate = true;
                vGateToErase.emplace_back(gateId);
            }
        }
        for(unsigned_t gateId : vGateToErase){
            pSol->sReadyMappedGate.erase(gateId);
            pSol->vExecutedGate.emplace_back(gateId);
        }
        vGateToErase.clear();
        for(unsigned_t gateId : pSol->sUnmappedGate){
            Gate & gate = _pCircuit->gate(gateId);
            bool canExecute = true;
            for(unsigned_t pId : _vsGateUnexecutedParent[gateId]){
                if(_pCircuit->gate(pId).nTargetQubit() == 1){
                    // if its parent gate is a single qubit gate, it will not added to sReadyMappedGate or sUnmappedGate.
                    if(!_vExecutedGate[pId]){
                        unsigned_t tmp = pId;
                        while(_pCircuit->gate(tmp).nTargetQubit() == 1){
                            if(_vsGateUnexecutedParent[tmp].size() > 0){
                                tmp = (*_vsGateUnexecutedParent[tmp].begin());
                            }
                            else{
                                break;
                            }
                        }
                        if(pSol->sReadyMappedGate.find(tmp) != pSol->sReadyMappedGate.end() || pSol->sUnmappedGate.find(tmp) != pSol->sUnmappedGate.end()){
                            canExecute = false;
                            break;
                        }
                    }
                }
                else if(pSol->sReadyMappedGate.find(pId) != pSol->sReadyMappedGate.end() || pSol->sUnmappedGate.find(pId) != pSol->sUnmappedGate.end()){
                    canExecute = false;
                    break;
                }
            }
            if(canExecute){
                if(gate.nTargetQubit() > 1){
                    q0 = gate.targetProgramQubit(0);
                    q1 = gate.targetProgramQubit(1);
                    if( _pDevice->isAdjacent(_vCurSolMappingPro2Phy[q0], _vCurSolMappingPro2Phy[q1])){
                        keepFindNewGate = true;
                        #ifdef DEBUG
                            cerr << "In sUnmappedGate: " << gateId << " can be executed " << q0 << "("<< _vCurSolMappingPro2Phy[q0] << ") " << q1 << "("<< _vCurSolMappingPro2Phy[q1] << ")" << endl;
                        #endif
                        pSol->vExecutedGate.emplace_back(gateId);
                        // for(unsigned_t i : _vsGateChild[gateId]){
                        //     _vsGateUncollectedParent[i].erase(gateId);
                            // if(!_vExecutedGate[i] && _vsGateUncollectedParent[i].size() == 0 && (find(_sCurTargetGates.begin(), _sCurTargetGates.end(), i) == _sCurTargetGates.end())){
                            //     collectCurTargetGatesFromGate(i);
                            // }
                        // }
                    }
                    else{
                        pSol->sReadyMappedGate.emplace(gateId);
                        #ifdef DEBUG
                            cerr << "In sUnmappedGate: " << gateId << " is ready to be executed " << q0 << "("<< _vCurSolMappingPro2Phy[q0] << ") " << q1 << "("<< _vCurSolMappingPro2Phy[q1] << ")" << endl;
                        #endif
                    }
                }
                else{
                    keepFindNewGate = true;
                    pSol->vExecutedGate.emplace_back(gateId);
                    // for(unsigned_t i : _vsGateChild[gateId]){
                    //     _vsGateUncollectedParent[i].erase(gateId);
                        // if(!_vExecutedGate[i] && _vsGateUncollectedParent[i].size() == 0 && (find(_sCurTargetGates.begin(), _sCurTargetGates.end(), i) == _sCurTargetGates.end())){
                        //     collectCurTargetGatesFromGate(i);
                        // }
                    // }
                }
                vGateToErase.emplace_back(gateId);
            }
        }
        for(unsigned_t gateId : vGateToErase){
            pSol->sUnmappedGate.erase(gateId);
        }
    }
}

double_t aRouter::calHeuristicCost(aSolution* pSol){
    // return 0;
    double_t cost = 0;
    double_t costUnmap = 0;
    unsigned_t p0, p1, count = 0;
    for (unsigned_t gId : pSol->sReadyMappedGate){
        Gate& gate = _pCircuit->gate(gId);
        if(gate.nTargetQubit() > 1){
            p0 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(0)];
            p1 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)];
            if(_pDevice->getDistance(p0, p1) > 0){
                // if(_aRouterParam.is_all_commute){
                //     cost += (double_t)pow(_pDevice->getDistance(p0, p1), 2);
                // }
                // else{
                //     cost += (double_t)(_pDevice->getDistance(p0, p1));
                // }
                // cost += (double_t)pow(_pDevice->getDistance(p0, p1), 2);
                cost += (double_t)(_pDevice->getDistance(p0, p1));
                ++count;
            }
        }
    }
    // return cost;
    for (unsigned_t gId : pSol->sUnmappedGate){
        Gate& gate = _pCircuit->gate(gId);
        if(gate.nTargetQubit() > 1){
            p0 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(0)];
            p1 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)];
            if(_pDevice->getDistance(p0, p1) > 0){
                // costUnmap += ((double_t)(_pDevice->getDistance(p0, p1))*(1.0/(double_t)_vsGateUnexecutedParent[gId].size()));
                costUnmap += ((double_t)(_pDevice->getDistance(p0, p1)));
                // cost += (double_t)_pCircuit->nProgramQubit();
                ++count;
            }
        }
    }
    if(pSol->sUnmappedGate.size() > 0){
         costUnmap /= (double_t)(pSol->sUnmappedGate.size()*_pCircuit->nProgramQubit()*2);
    }
    // for (unsigned_t gId : pSol->sUnmappedGate){
    //     Gate& gate = _pCircuit->gate(gId);
    //     bool calculateWeight = true;
    //     for(unsigned_t j : _vsGateUnexecutedParent[gId]){
    //         if(pSol->sUnmappedGate.count(j)){
    //             calculateWeight = false;
    //             break;
    //         }
    //     }
    //     if(gate.nTargetQubit() > 1 && calculateWeight){
    //         p0 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(0)];
    //         p1 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)];
    //         if(_pDevice->getDistance(p0, p1) > 0){
    //             costUnmap += ((double_t)(_pDevice->getDistance(p0, p1)));
    //             ++count;
    //         }
    //     }
    // }
    // if(count > 0){
    //     costUnmap /= (double_t)(count*_pCircuit->nProgramQubit()*2);
    // }
    // if(_aRouterParam.is_all_commute){
    //     return cost;
    // }
    // if(isnan(cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + costUnmap/2 + count)){
    //     cerr << "cost: " << cost << endl;
    //     cerr << "pSol->sReadyMappedGate.size(): " << pSol->sReadyMappedGate.size() << endl;
    //     cerr << "costUnmap: " << costUnmap << endl;
    //     cerr << "count: " << count << endl;
    //     cerr << "h_cost_1: " << cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) << endl;
    //     cerr << "h_cost: " << cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + costUnmap + ((double_t)(pSol->sUnmappedGate.size()+pSol->sReadyMappedGate.size())) << endl;
    //     getchar();
    // }
    if(_aRouterParam.use_admissible_heuristic && pSol->sReadyMappedGate.size() > 0){
        cost = cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + costUnmap/2 + count;
    }
    // cost = cost/(double_t)pSol->sReadyMappedGate.size() + costUnmap/(double_t)pSol->sUnmappedGate.size();
    else if(pSol->sReadyMappedGate.size() > 0){
        cost = cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + costUnmap + ((double_t)(pSol->sUnmappedGate.size()+pSol->sReadyMappedGate.size()));///(double_t)_sCurTargetGates.size());
    }
        // cost = cost / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + costUnmap + ((double_t)(count));///(double_t)_sCurTargetGates.size());
    return cost;
}

void aRouter::updateCurMappingForSol(aSolution* pSol){
    vector<aSolution*> vpSol= {pSol};
    while(pSol->pParent != nullptr){
        vpSol.emplace_back(pSol->pParent);
        pSol = pSol->pParent;
    }
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        _vCurSolMappingPhy2Pro[i] = -1;
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy[i] = _vCurMapping[i];
        _vCurSolMappingPhy2Pro[_vCurMapping[i]] = i;
    }
    
    for(int_t i = vpSol.size() - 1; i >= 0; --i){
        swapQubit(vpSol[i]->swapEdge);
    }
}

void aRouter::collectNeighbor(aSolution* pSol){
    // 
    unsigned_t q0, q1, p0, p1, p2, disBound;
    // aSolution& sol = _vAllSolution[idx];
    // 211287
    unordered_set<unsigned_t> sEdgeToHaveSwap;
    sEdgeToHaveSwap.clear();
    bool expandNode;
    for(unsigned_t gId : pSol->sReadyMappedGate){
        Gate & gate = _pCircuit->gate(gId);
        if(gate.nTargetQubit() > 1){
            q0 = gate.targetProgramQubit(0);
            q1 = gate.targetProgramQubit(1);
            p0 = _vCurSolMappingPro2Phy[q0];
            p1 = _vCurSolMappingPro2Phy[q1];
            disBound =  _pDevice->getDistance(p1,p0);
            // change p0 to p2
            unordered_set<int_t>& qubitRegion0 = _vsQubitRegion[q0]; //_pCircuit->sQubitRegion(q0);
            unordered_set<int_t>& qubitRegion1 = _vsQubitRegion[q1];
            for (unsigned_t eId : _pDevice->qubit(p0).vSpanEdge){
                // !  enable routing region
                p2 = (_pDevice->edge(eId).qubitId1() == p0) ? _pDevice->edge(eId).qubitId2() : _pDevice->edge(eId).qubitId1();
                if(_aRouterParam.restrict_region){
                    expandNode = false;
                    // p0 is in q0's qubit refion
                    if(qubitRegion0.count(p0) > 0){
                        if(qubitRegion0.count(p2) > 0){
                            expandNode = true;
                        }
                        else{
                            expandNode = ((double_t)(rand() % 1000 / 1000.0) < _aRouterParam.expand_prob);
                        }
                    }
                    else{
                        expandNode = true;
                    }
                }
                else{
                    expandNode = true;
                }
                expandNode = (expandNode && (_pDevice->getDistance(p1, p2) <= disBound));
                if(expandNode){
                    sEdgeToHaveSwap.insert(eId);
                }
            }
            // change p1 to p2
            for (unsigned_t eId : _pDevice->qubit(p1).vSpanEdge){
                // ! not enable routing region
                p2 = (_pDevice->edge(eId).qubitId1() == p1) ? _pDevice->edge(eId).qubitId2() : _pDevice->edge(eId).qubitId1();
                if(_aRouterParam.restrict_region){
                    expandNode = false;
                    if(qubitRegion1.count(p1) > 0){
                        if(qubitRegion1.count(p2) > 0){
                            expandNode = true;
                        }
                        else{
                            expandNode = ((double_t)(rand() % 1000 / 1000.0) < _aRouterParam.expand_prob);
                        }
                    }
                    else{
                        expandNode = true;
                    }
                }
                else{
                    expandNode = true;
                }
                expandNode = (expandNode && (_pDevice->getDistance(p0, p2) <= disBound));
                if(expandNode){
                    sEdgeToHaveSwap.insert(eId);
                }
            }
        }
        
    }
    unsigned_t swapIdx;
    for (unsigned_t eId : sEdgeToHaveSwap){
        aSolution * pAdjSol = new aSolution(pSol->sUnmappedGate, pSol->sReadyMappedGate);
        _sAllSolution.insert(pAdjSol);
        // _vpAllSolution.emplace_back(pAdjSol);
        pSol->vAdjs.emplace_back(pAdjSol);

        pAdjSol->done = false;
        // insert swap gate
        pAdjSol->swapEdge = eId;
        pAdjSol->pParent = pSol;
        // #ifdef DEBUG
        //     cerr << "construct neighbor with swap edge " << eId << "(" << p0 << "," << p1 << ")" << endl;
        //     printSolution(pAdjSol);
        // #endif
        // ++solIdx;
    }
    #ifdef DEBUG
        cerr << endl << "neighbor size: " << pSol->vAdjs.size() << endl;
    #endif
}


void aRouter::updateCircuit(aSolution* pSol){
    // backtrack solutions
    vector<aSolution*> vpSol= {pSol};
    bool printIsTrue = 0; //(pSol->sReadyMappedGate.size() >  0);//(pSol->sReadyMappedGate.size() >  0 && _pCircuit->nProgramQubit() < 25);
    for(unsigned_t i : pSol->sUnmappedGate){
        for(unsigned_t j : _vsGateChild[i]){
            _vsGateUncollectedParent[j].insert(i);
        }
    }
    while(pSol->pParent != nullptr){
        vpSol.emplace_back(pSol->pParent);
        pSol = pSol->pParent;
    }
    // printIsTrue = vpSol.size() < 10;
    // mark gates in _sCurTargetGates as executed
    // for(unsigned_t i : _sCurTargetGates){
    //     _vExecutedGate[i] = true;
    //     for(unsigned_t j : _vsGateChild[i]){
    //         _vsGateUnexecutedParent[j].erase(i);
    //     }
    //     ++_executedGateCount;
    // }
    // collect a set of gates that are ready to be executed 
    unordered_set<unsigned_t> sGateReadyToMap;
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        if(!_vExecutedGate[i] && _vsGateUnexecutedParent[i].size() == 0){
            sGateReadyToMap.insert(i);
        }
    }
    // check each layer
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        _vCurSolMappingPhy2Pro[i] = -1;
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy[i] = _vCurMapping[i];
        _vCurSolMappingPhy2Pro[_vCurMapping[i]] = i;
    }
    if(printIsTrue){
        _pDevice->printDevice();
        cout << "vpSol.size(): " << vpSol.size() << endl;
        for(int_t i = vpSol.size() - 1; i >= 0; --i){
            swapQubit(vpSol[i]->swapEdge);
            aSolution & sol = *vpSol[i];
            cerr << "Print solution" << endl;
            cerr << "g cost: " << sol.g_cost << ", f cost: " << sol.f_cost << endl;
            cerr << "#neighbor: " << sol.vAdjs.size() << endl;
            cerr << "vExecutedGate: (" << sol.vExecutedGate.size() << ")";
            for(unsigned_t i : sol.vExecutedGate){
                cerr << " " << i;
            }
            cerr << endl;
            cerr << "sReadyMappedGate: (" << sol.sReadyMappedGate.size() << ")";
            for(unsigned_t i : sol.sReadyMappedGate){
                Gate & gate = _pCircuit->gate(i);
                cerr << " " << i << "("<< _vCurSolMappingPro2Phy[gate.targetProgramQubit(0)] << "," << _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)] << ")";
            }
            // getchar();
            cerr << endl;
            cerr << "sUnmappedGate: (" << sol.sUnmappedGate.size() << ")";
            // getchar();
            for(unsigned_t i : sol.sUnmappedGate){
                cerr << " " << i;
            }
            // getchar();
            cerr << endl;
            if(sol.swapEdge > -1)
                cerr << "swapEdge: " << sol.swapEdge << "("<< _pDevice->edge(sol.swapEdge).qubitId1() << "," << _pDevice->edge(sol.swapEdge).qubitId2() << ")" << endl;
            cerr << "print its neighbor:" << endl;
            for(auto i : sol.vAdjs){
                printSolution(i);
            }
            getchar();
        }
    }
    unsigned_t circuitDepthBias = _pCircuit->circuitDepth(); 
    bool canExecute, keepFindNewGate;
    #ifdef DEBUG
        cerr << "start glue solution" << endl;
    #endif 
    for(int_t i = vpSol.size() - 1; i >= 0; --i){
        swapQubit(vpSol[i]->swapEdge);
        // add swap gate
        if(vpSol[i]->swapEdge >= 0){
            Edge & edge = _pDevice->edge(vpSol[i]->swapEdge);
            int_t proQ0 = _vCurSolMappingPhy2Pro[edge.qubitId1()], proQ1 = _vCurSolMappingPhy2Pro[edge.qubitId2()];
            if(_vQubitHasGate[edge.qubitId1()] || _vQubitHasGate[edge.qubitId2()]) {
                // if(printIsTrue){
                //     cerr << "add a swap gate _pCircuit->nSwapGate()" << endl;
                // }
                _pCircuit->addSwapGate(_pCircuit->nSwapGate(), edge.qubitId1(), edge.qubitId2(), circuitDepthBias);
                _vQubitHasGate[edge.qubitId1()] = true;
                _vQubitHasGate[edge.qubitId2()] = true;
            }
            else{
                if(proQ0 > -1){
                    _pCircuit->setInitialMapping(proQ0, edge.qubitId1());
                }
                if(proQ1 > -1){
                    _pCircuit->setInitialMapping(proQ1, edge.qubitId2());
                }
            }
            if(proQ0 > -1){
                _pCircuit->addQubitRegion(proQ0, edge.qubitId1());
                // _vsQubitRegion[proQ0].insert(edge.qubitId1());
            }
            if(proQ1 > -1){
                _pCircuit->addQubitRegion(proQ1, edge.qubitId2());
                // _vsQubitRegion[proQ1].insert(edge.qubitId2());
            }
            // if((edge.qubitId1() == 9 && edge.qubitId2() == 15) || (edge.qubitId1() == 15 && edge.qubitId2() == 9)){
            //     cerr << "_vQubitHasGate[edge.qubitId1()]: " << _vQubitHasGate[edge.qubitId1()] << endl;
            //     cerr << "_vQubitHasGate[edge.qubitId2()]: " << _vQubitHasGate[edge.qubitId2()] << endl;
            //     _pCircuit->printCircuitLayout();
            //     getchar();
            // }
        }
        ++circuitDepthBias;
        // set gate execution time and target qubits
        #ifdef DEBUG
            cerr << endl;
            cerr << "i = " << i << endl;
            printSolution(vpSol[i]);
            cerr << "circuitDepthBias : " << circuitDepthBias << endl;
        #endif
        for(unsigned_t j : vpSol[i]->vExecutedGate){
            if(!_vExecutedGate[j]){
                Gate& gate = _pCircuit->gate(j);
                gate.setExecutionTime(circuitDepthBias);
                for(unsigned_t k = 0; k < gate.nTargetQubit(); ++k){
                    gate.setTargetPhysicalQubit(k, _vCurSolMappingPro2Phy[gate.targetProgramQubit(k)]);
                    _vQubitHasGate[_vCurSolMappingPro2Phy[gate.targetProgramQubit(k)]] = true;
                }
                _vExecutedGate[j] = true;
                for(unsigned_t k : _vsGateChild[j]){
                    _vsGateUnexecutedParent[k].erase(j);
                    _vsGateUncollectedParent[k].erase(j);
                    if(_vsGateUnexecutedParent[k].size() == 0 && !_vExecutedGate[k]){
                        sGateReadyToMap.insert(k);
                        keepFindNewGate = true;
                    }
                }
                #ifdef DEBUG
                    cerr << "gate " << j << " is executed on program qubit " << gate.targetProgramQubit(0) << "("<< gate.targetPhysicalQubit(0) << ")";
                    if(gate.nTargetQubit() == 2){
                        cerr << ", " << gate.targetProgramQubit(1) << "("<< gate.targetPhysicalQubit(1) << ") on time " << circuitDepthBias << endl; 
                    }
                    else{
                        cerr << " on time " << circuitDepthBias << endl;
                    }
                #endif
                ++_executedGateCount;
            }
        }
        // check if other gates in the circuit can be executed. update executedGateCount
        keepFindNewGate = true;
        unordered_set<unsigned_t> sNewGateReadyToMap;
        unordered_set<unsigned_t> sGateIsExectued;
        while(keepFindNewGate){
            keepFindNewGate = false;
            sNewGateReadyToMap.clear();
            sGateIsExectued.clear();
            for(unsigned_t j : sGateReadyToMap){
                Gate& gate = _pCircuit->gate(j);
                canExecute = (((gate.nTargetQubit() == 1) || _pDevice->isAdjacent(_vCurSolMappingPro2Phy[gate.targetProgramQubit(0)], _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)])) && !_vExecutedGate[j]);
                if(canExecute){
                    sGateIsExectued.insert(j);
                    _vExecutedGate[j] = true;
                    ++_executedGateCount;
                    Gate& gate = _pCircuit->gate(j);
                    gate.setExecutionTime(circuitDepthBias);
                    for(unsigned_t k = 0; k < gate.nTargetQubit(); ++k){
                        gate.setTargetPhysicalQubit(k, _vCurSolMappingPro2Phy[gate.targetProgramQubit(k)]);
                        _vQubitHasGate[_vCurSolMappingPro2Phy[gate.targetProgramQubit(k)]] = true;
                    }
                    for(unsigned_t k : _vsGateChild[j]){
                        _vsGateUnexecutedParent[k].erase(j);
                        _vsGateUncollectedParent[k].erase(j);
                        if(_vsGateUnexecutedParent[k].size() == 0 && !_vExecutedGate[k]){
                            sNewGateReadyToMap.insert(k);
                            keepFindNewGate = true;
                        }
                    }
                    #ifdef DEBUG
                        cerr << "gate " << j << " can be directly executed on program qubit " << gate.targetProgramQubit(0) << "("<< gate.targetPhysicalQubit(0) << ")";
                        if(gate.nTargetQubit() == 2){
                            cerr << ", " << gate.targetProgramQubit(1) << "("<< gate.targetPhysicalQubit(1) << ") on time " << circuitDepthBias << endl; 
                        }
                        else{
                            cerr << " on time " << circuitDepthBias << endl;
                        }
                    #endif
                }
            }
            for(unsigned_t j : sGateIsExectued){
                sGateReadyToMap.erase(j);
            }
            for(unsigned_t j : sNewGateReadyToMap){
                sGateReadyToMap.insert(j);
            }
        }
        
        // getchar();
    }
    _pCircuit->setCircuitDepth(circuitDepthBias);
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurMapping[i] = _vCurSolMappingPro2Phy[i];
    }
    // if all gates are routed, set final mapping
    if(_executedGateCount == _pCircuit->nGate()){
        for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); i++){
            _pCircuit->setFinalMapping(i, _vCurSolMappingPro2Phy[i]);
        }
    }
    // cerr << "unexecuted gate: ";
    // for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
    //     if(!_vExecutedGate[i]){
    //         cerr << i << " ";
    //     }
    // }
    // cerr << endl;
    // getchar();
}

void aRouter::printSolution(aSolution* pSol){
    aSolution & sol = *pSol;
    cerr << "Print solution" << endl;
    cerr << "g cost: " << sol.g_cost << ", f cost: " << sol.f_cost << endl;
    // cerr << "Current mapping: (" << sol.vCurMapping.size() << ")" << endl;
    // // getchar();
    // for(unsigned_t i = 0; i < sol.vCurMapping.size(); ++i){
    //     cerr << "  pro " << i << " to phy " << sol.vCurMapping[i] << endl;
    // }
    // getchar();
    cerr << "vExecutedGate: (" << sol.vExecutedGate.size() << ")";
    // getchar();
    for(unsigned_t i : sol.vExecutedGate){
        cerr << " " << i;
    }
    cerr << endl;
    cerr << "sReadyMappedGate: (" << sol.sReadyMappedGate.size() << ")";
    // getchar();
    for(unsigned_t i : sol.sReadyMappedGate){
        cerr << " " << i;
    }
    // getchar();
    cerr << endl;
    cerr << "sUnmappedGate: (" << sol.sUnmappedGate.size() << ")";
    // getchar();
    for(unsigned_t i : sol.sUnmappedGate){
        cerr << " " << i;
    }
    // getchar();
    cerr << endl;
    cerr << "swapEdge: " << sol.swapEdge << endl;
}

void aRouter::cleanAllNode(){
    _priorityQ.clear();
    for_each(_sAllSolution.begin(), _sAllSolution.end(), [](aSolution* obj){ delete obj; });
    _sAllSolution.clear();
}

bool aRouter::checkMemUsage(){
    // double vm_usage     = 0.0; 
    // double resident_set = 0.0;
    // the two fields we want
    unsigned long vsize; // virtual memory usage in byte
    long rss;
        {
            std::string ignore;
            std::ifstream ifs("/proc/self/stat", std::ios_base::in);
            ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                    >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                    >> ignore >> ignore >> vsize >> rss;
        }
    // long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    //     vm_usage = vsize / 1024.0;
    //     resident_set = rss * page_size_kb;
    cerr << "vsize: " << vsize << endl;
    return (_aRouterParam.mem_limit < vsize);
}

MOLSQ_NAMESPACE_CPP_END

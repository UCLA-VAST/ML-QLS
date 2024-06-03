/***********************************************************************
  File        [ aRouter_dev.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ rotuer ]
  Synopsis    [ aRouter_dev class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#include "router/aRouter_dev.hpp"

MOLSQ_NAMESPACE_CPP_START

void aRouter_dev::setQubitRegion(Circuit& cir){
    _vsQubitRegion.clear(); 
    _vsQubitRegion.resize(cir.nProgramQubit());
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        unordered_set<int_t> & sQubitRegion = cir.sQubitRegion(i);
        for(int_t j : sQubitRegion){
            _vsQubitRegion[i].insert(j);
        }
    }
}

bool aRouter_dev::run(Circuit& cir, Device& device){
    // _aRouter_devParam.node_limit = 350000 - (cir.nGate()*30);
    // if (_aRouter_devParam.node_limit < 3000){
    //     _aRouter_devParam.node_limit = 3000;
    // }
    // if(_aRouter_devParam.is_all_commute){
    //     _aRouter_devParam.node_limit = 600000 - (cir.nGate()*30);
    //     _aRouter_devParam.timeout = 300;
    //     if (_aRouter_devParam.node_limit < 3000){
    //         _aRouter_devParam.node_limit = 3000;
    //     }
    //     _aRouter_devParam.node_num_after_trim = 1000;
    // }
    
    if(_verbose > 0)
        fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter_dev: Start\n");
    _pCircuit = &cir;
    _pCircuit->clearSwap();
    _pCircuit->setCircuitDepth(0);
    _pDevice = &device;
    _executedGateCount = 0;
    _vCurMapping.clear();
    _vCurMapping.resize(cir.nProgramQubit(), 0);
    _vCurSolMappingPro2Phy.clear();
    _vCurSolMappingPro2Phy.resize(cir.nProgramQubit(), 0);
    _vCurSolMappingPhy2Pro.clear();
    _vCurSolMappingPhy2Pro.resize(_pDevice->nQubit(), -1);
    _mvSwapEdge.clear();
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurMapping.at(i) = _pCircuit->initialMapping(i);
        _pCircuit->addQubitRegion(i, _pCircuit->initialMapping(i));
        
    }
    _timer.start(TimeUsage::FULL);
    // if(_aRouter_devParam.restrict_region && _vsQubitRegion.size() != cir.nProgramQubit()){
    //     setQubitRegion(cir);
    // }
    cir.resetQubitRegion();
    constructDependencyInfo();
    // start swap insertion
    #ifdef DEBUG
        cerr << "start astar" << endl;
    #endif
    if(_verbose > 0){
        fprintf(stdout, "[INFO] #Gate: %d\n", _pCircuit->nGate());
        // fprintf(stdout, "[INFO] #Routed Gate: %d\r", _executedGateCount);
    }
    _timer.restart(_aRouter_devParam.timeout);
    // run astar
    shared_ptr<aSolution> pFinalSol = run_astar();
    if(pFinalSol == nullptr){
        if(_verbose > 0){
            fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter: Finish\n");
        }
        cleanAllNode();     
        return false;
    }
    _cost = pFinalSol->g_cost;
    updateCircuit(pFinalSol);
    cleanAllNode();       
    _pCircuit->setCircuitDepth(_pCircuit->circuitDepth());
    if(_verbose > 0){
        // fprintf(stdout, "[INFO] #Routed Gate: %d\n", _executedGateCount);
        fprintf(stdout, "[INFO] A*-based SWAP insertion aRouter_dev: Finish\n");
    }
    return true;
}

void aRouter_dev::constructDependencyInfo(){
    // _vGateWeight.clear();
    // _vGateWeight.resize(_pCircuit->nGate(), 1);
    _vsGateParent.clear();
    _vsGateParent.resize(_pCircuit->nGate(), unordered_set<unsigned_t>());
    _vsGateChild.clear();
    _vsGateChild.resize(_pCircuit->nGate(), unordered_set<unsigned_t>());
    if(_aRouter_devParam.is_all_commute){
        // cerr << "all commute" << endl;
        return;
    }
    for(pair<unsigned_t, unsigned_t>& pDependency : (*_pCircuit->pvpGateDependency())){
        _vsGateParent.at(pDependency.second).insert(pDependency.first);
        _vsGateChild.at(pDependency.first).insert(pDependency.second);
    }
    // double_t maxWeight = 0.0, weightBase = 10 / _pCircuit->nGate();
    // for(int_t i = _pCircuit->nGate() - 1; i >= 0; --i){
    //     if(_pCircuit->gate(i).nTargetQubit() > 1){
    //         maxWeight = 0.0;
    //         for(unsigned_t j : _vsGateParent[i]){
    //             maxWeight = min(maxWeight, _vGateWeight[j]);
    //         }
    //         _vGateWeight[i] = maxWeight + weightBase;
    //     }
    // }
    #ifdef DEBUG
        cout << "Print _vsGateParent and _vsGateChild for gates" << endl;
        for (unsigned_t i = 0; i < _vsGateParent.size(); ++i){
            cout << "Print _vsGateParent for gate " << i << ": ";
            for (unsigned_t gId : _vsGateParent.at(i)){
                cout << gId << " ";
            }
            cout << endl;
            cout << "Print _vsGateChild for gate " << i << ": ";
            for (unsigned_t gId : _vsGateChild.at(i)){
                cout << gId << " ";
            }
            cout << endl;
        }
    #endif
}

shared_ptr<aRouter_dev::aSolution> aRouter_dev::run_astar() {
    // Early return
    // aSolution& initSol = _vAllSolution[0];
    // bool toPrint = _pCircuit->nGate() < 30;
    _vAllSolution.clear();
    _vAllSolution.emplace_back(std::make_shared<aSolution>());
    #ifdef DEBUG
        cerr << "start initialization" << endl;
    #endif
    initialize(_vAllSolution.at(0));
    if (checkFinalState(_vAllSolution.at(0))){
        return _vAllSolution.at(0);
    }
    // Initialize priority queue (frontier)
    PQueue_t priorityQ;
    priorityQ.push(_vAllSolution.at(0));
    // printSolution(_vAllSolution.at(0));
    // getchar();
    // Declare Node* pCurNode for the extracted min node
    // Start Dijstra
    // bool toPrint = false;
    unsigned_t trimStateCount = 0;
    // bool toPrint = (_pCircuit->nGate() < 40);
    while ( !priorityQ.empty() ) {
        // Extract min node
        shared_ptr<aSolution> pCurNode = priorityQ.top();
        priorityQ.pop();
        assert(pCurNode->done == false);
        pCurNode->done = true; // mark the extracted node as DONE

        
        // Termination condition
        if (checkFinalState(pCurNode)){
            // if(pCurNode->g_cost < 20){
            //     printSolution(pCurNode);
            // }
            return pCurNode;
        }
        updateCurMappingForSol(pCurNode);
        #ifdef DEBUG
            cerr << "====================================" << endl;
            cerr << "explore node" << endl;
            printSolution(pCurNode);
            getchar();
        #endif
        // if (toPrint){
        //     cerr << "====================================" << endl;
        //     cerr << "node size" <<  _vAllSolution.size() << endl;
        //     cerr << "explore node" << endl;
        //     printSolution(pCurNode);
        //     cerr << endl;
        //     getchar();
        // }
        // if(toPrint){
        //     cerr << endl;
        //     cerr << "print _vCurSolMappingPro2Phy" << endl;
        //     for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //         cerr << i << ": " << _vCurSolMappingPro2Phy.at(i) << endl;
        //     }
        //     cerr << "print _vCurSolMappingPhy2Pro" << endl;
        //     for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        //         cerr << i << ": " << _vCurSolMappingPhy2Pro.at(i) << endl;
        //     }
        //     printSolution(pCurNode);
        //     getchar();
        // }
        // construct adjacent nodes
        // if(toPrint){
        //     cerr << "begin collect neighbor" << endl;
        // }
        collectNeighbor(pCurNode);
        // if(toPrint){
        //     cerr << "finiish collect neighbor" << endl;
        // }
        // Relax neighbors
        #ifdef DEBUG
            cerr << endl << "In astar: neighbor size: " << pCurNode->vAdjs.size() << endl;
        #endif
        for (shared_ptr<aSolution> pSol : pCurNode->vAdjs) {
            // Skip adjNode which is "done"
            if (pSol->done)
                continue;
            // update _vCurSolMappingPro2Phy and _vCurSolMappingPhy2Pro based on swapEdge
            swapQubit(pSol->swapEdge); // qubit is swapped back in updateHeuristicCost
            // Initialize node
            initNode(pSol);
            pSol->g_cost = pCurNode->g_cost + 1;            
            pSol->f_cost = pSol->g_cost + calHeuristicCost(pSol);
            // Update pqueue for pAdjNode
            priorityQ.push(pSol); 
            // if(toPrint){
            //     printSolution(pCurNode);
            // }
            // recover _vCurSolMappingPro2Phy and _vCurSolMappingPhy2Pro based on swapEdge
            swapQubit(pSol->swapEdge);
            #ifdef DEBUG
                printSolution(pSol);
            #endif
        }
        // if(toPrint){
        //     cerr << "update neighbor cost" << endl;
        // }
        // condition to trim states
        if(_aRouter_devParam.kill_if_timeout){
            if(_timer.isTimeout()){
                return nullptr;
            }
        }
        if(toTrimState()){
            // toPrint = true;
            
            ++trimStateCount;
            unordered_set<shared_ptr<aSolution>> sNewSolution;
            for(unsigned_t i = 0; i < _aRouter_devParam.node_num_after_trim && !priorityQ.empty(); ++i){
                shared_ptr<aSolution> pSol = priorityQ.top();
                sNewSolution.insert(pSol);
                priorityQ.pop();
            }
            trimState(sNewSolution);
            while(!priorityQ.empty()){
                shared_ptr<aSolution> pSol = priorityQ.top();
                // pSol->pParent = nullptr;
                // pSol->vAdjs.clear();
                priorityQ.pop();
            }
            // _aRouter_devParam.weight *= 1.1;
            // _aRouter_devParam.weight = min(_aRouter_devParam.weight, 100.0);
            bool isTimeout = _timer.isTimeout();
            if(isTimeout){
                fprintf(stdout, "[INFO] Timeout in aRouter_dev. Increase weight for unexecuted gates\n");
                _timer.restart(_aRouter_devParam.timeout);
                _aRouter_devParam.weight += 0.1;
            }
            for(shared_ptr<aSolution> pSol : _vAllSolution){
                if(isTimeout){
                    pSol->f_cost = pSol->g_cost + calHeuristicCost(pSol);
                }
                priorityQ.push(pSol);
            }
            // cerr << "f_cost: " <<  priorityQ.top()->f_cost << endl;
            // cerr << "g_cost: " <<  priorityQ.top()->g_cost << endl;
            // cerr << "unexecutedGateCount: " <<  priorityQ.top()->unexecutedGateCount << endl;
            
            if(_verbose > 1){
                fprintf(stdout, "[INFO] #Routed Gate: %d, SWAP count: %0f\r", _pCircuit->nGate()-priorityQ.top()->unexecutedGateCount, priorityQ.top()->g_cost);
                fprintf(stdout, "[INFO] g_cost: %.4f\n", priorityQ.top()->g_cost);
                fprintf(stdout, "[INFO] f_cost: %.4f\n", priorityQ.top()->f_cost);
            }
        }
    }
    return nullptr;
}

bool aRouter_dev::toTrimState(){
    if((int_t)_vAllSolution.size() > _aRouter_devParam.node_limit){
        // fprintf(stdout, "[INFO] Too many nodes (>%d)\n", _aRouter_devParam.node_limit);
        return true;
    }
    return false;
}

void aRouter_dev::trimState(unordered_set<shared_ptr<aSolution>> & sNewSolution){
    // save n state
    unordered_map<shared_ptr<aSolution>, vector<unsigned_t>> mvSwapEdge;
    // printSolution(*sNewSolution.begin());
    for(shared_ptr<aSolution> pNewSol : sNewSolution){
        // cerr << "new sol " << pNewSol << " (g cost: " << pNewSol->g_cost << ")";
        shared_ptr<aSolution> pSol = pNewSol;
        vector<shared_ptr<aSolution>> vpSol = {pSol};
        if(_mvSwapEdge.find(pNewSol) != _mvSwapEdge.end()){
            mvSwapEdge[pNewSol] = _mvSwapEdge[pNewSol];
            continue;
        }
        while(_mvSwapEdge.find(pSol->pParent) == _mvSwapEdge.end() && pSol->pParent != nullptr){
            vpSol.emplace_back(pSol->pParent);
            pSol = pSol->pParent;
        }
        // cerr << "keep pNewSol " << pNewSol << endl;
        // cerr << "reach parent state " << pSol->pParent << endl;
        mvSwapEdge[pNewSol] = vector<unsigned_t>();
        if(pSol->pParent == nullptr){
            // cerr << " has no parent swaps" << endl;
            for(int_t i = vpSol.size() - 2; i >= 0; --i){
                mvSwapEdge[pNewSol].emplace_back(vpSol.at(i)->swapEdge);
            }
        }
        else{
            // cerr << " has parent swap list " << pSol->pParent << " with size " << _mvSwapEdge[pSol->pParent].size() << endl;
            for(int_t i : _mvSwapEdge[pSol->pParent]){
                mvSwapEdge[pNewSol].emplace_back(i);
            }
            for(int_t i = vpSol.size() - 1; i >= 0; --i){
                mvSwapEdge[pNewSol].emplace_back(vpSol.at(i)->swapEdge);
            }
        }
    }
    // cerr << "_mvSwapEdge before tirm state" << endl;
    // for(auto m : _mvSwapEdge){
    //     cerr << "pSol: " << m.first << " , swap size(): " << m.second.size() << endl;
    //     // for(int_t i : _mvSwapEdge[pSol]){
    //     //     cerr << "add mvedge " <<  i << " " << _pDevice->edge(i).qubitId1() << " " << _pDevice->edge(i).qubitId2() << endl;
    //     // }   
    // }
    _mvSwapEdge = mvSwapEdge;
    // cerr << "_mvSwapEdge after tirm state" << endl;
    // for(auto m : _mvSwapEdge){
    //     cerr << "pSol: " << m.first << " , swap size(): " << m.second.size() << endl;
    //     // for(int_t i : _mvSwapEdge[pSol]){
    //     //     cerr << "add mvedge " <<  i << " " << _pDevice->edge(i).qubitId1() << " " << _pDevice->edge(i).qubitId2() << endl;
    //     // }   
    // }

    // unsigned_t count = 0;
    for(shared_ptr<aSolution> pSol : _vAllSolution){
        pSol->pParent = nullptr;
        pSol->vAdjs.clear();
    }
    _vAllSolution.clear();
    for(shared_ptr<aSolution> pNewSol : sNewSolution){
        pNewSol->pParent = nullptr;
    }
    // _vAllSolution = move(vNewSolution);      
    _vAllSolution.reserve(sNewSolution.size());
    for (auto it = sNewSolution.begin(); it != sNewSolution.end(); ) {
        _vAllSolution.emplace_back(std::move(sNewSolution.extract(it++).value()));
    }
    #ifdef DEBUG
        cerr << "print _mvSwapEdge " << _mvSwapEdge.size() << endl;
        for(auto m : _mvSwapEdge){
            cerr << "m.first: " << m.first << endl;
            for(auto i : m.second){
                cerr << i << " ";
            }
            cerr << endl;
        }
    #endif
    // fprintf(stdout, "[INFO] Finish trimming states\n");
    _timer.restart(_aRouter_devParam.timeout);
}

bool aRouter_dev::checkFinalState(shared_ptr<aSolution> pSol){
    return (pSol->unexecutedGateCount == 0);
}

void aRouter_dev::initialize(shared_ptr<aSolution> pSol){
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        _vCurSolMappingPhy2Pro.at(i) = -1;
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy.at(i) = _vCurMapping.at(i);
        _vCurSolMappingPhy2Pro.at(_vCurMapping.at(i)) = i;
    }
    _vvQubitGate.clear();
    _vvQubitGate.resize(_pCircuit->nProgramQubit(), vector<unsigned_t>());
    pSol->vGateState.clear();
    pSol->vGateState.resize(_pCircuit->nGate(), make_pair(false, false));
    pSol->unexecutedGateCount = 0;
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        Gate & gate = _pCircuit->gate(i);
        if(gate.nTargetQubit() > 1){
            _vvQubitGate[gate.targetProgramQubit(0)].emplace_back(i);
            _vvQubitGate[gate.targetProgramQubit(1)].emplace_back(i);
            ++pSol->unexecutedGateCount;
        }
        if(_vsGateParent.at(i).size() == 0){
            pSol->vGateState.at(i).first = true;
            pSol->vGateState.at(i).second = false;
            pSol->sReadyMappedGate.insert(i);
        }
    }
    // _aRouter_devParam.weight = pSol->unexecutedGateCount
    initNode(pSol);
    pSol->g_cost = 0;
    pSol->f_cost = calHeuristicCost(pSol); 
    #ifdef DEBUG
        cerr << "initialize" << endl;
        printSolution(pSol);
    #endif
}

void aRouter_dev::swapQubit(int_t eId){
    if(eId > -1){
        // cerr << "before swap " << eId << "(" << _pDevice->edge(eId).qubitId1() << "," << _pDevice->edge(eId).qubitId2() << ")" << endl;
        // cerr << "print _vCurSolMappingPro2Phy" << endl;
        // for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPro2Phy.at(i) << endl;
        // }
        // cerr << "print _vCurSolMappingPhy2Pro" << endl;
        // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPhy2Pro.at(i) << endl;
        // }
        Edge & edge = _pDevice->edge(eId);

        int_t proId0 = _vCurSolMappingPhy2Pro.at(edge.qubitId1());
        int_t proId1 = _vCurSolMappingPhy2Pro.at(edge.qubitId2());
        if (0 <= proId0){
            _vCurSolMappingPro2Phy.at(proId0) = edge.qubitId2();
        }
        if (0 <= proId1){
            _vCurSolMappingPro2Phy.at(proId1) = edge.qubitId1();
        }
        _vCurSolMappingPhy2Pro.at(edge.qubitId2()) = proId0;
        _vCurSolMappingPhy2Pro.at(edge.qubitId1()) = proId1;
        // cerr << "after swap" << endl;
        // cerr << "print _vCurSolMappingPro2Phy" << endl;
        // for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPro2Phy.at(i) << endl;
        // }
        // cerr << "print _vCurSolMappingPhy2Pro" << endl;
        // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        //     cerr << i << ": " << _vCurSolMappingPhy2Pro.at(i) << endl;
        // }
    }
}


void aRouter_dev::initNode(shared_ptr<aSolution> pSol){
    vector<unsigned_t> vGateToErase;
    unordered_set<unsigned_t> sGateToAdd;
    bool canExecute, isReady, hasOneLevelDependency;
    unordered_set<unsigned_t> sAffectedNode;
    vector<unsigned_t> vGateExecuted;
    // cerr << endl;
    // cerr << "print _vCurSolMappingPro2Phy" << endl;
    // for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
    //     cerr << i << ": " << _vCurSolMappingPro2Phy.at(i) << endl;
    // }
    // cerr << "print _vCurSolMappingPhy2Pro" << endl;
    // for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
    //     cerr << i << ": " << _vCurSolMappingPhy2Pro.at(i) << endl;
    // }
    for(unsigned_t i : pSol->sReadyMappedGate){
        Gate& gate = _pCircuit->gate(i);
        canExecute = (((gate.nTargetQubit() == 1) || _pDevice->isAdjacent(_vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0)), _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1)))));
        if(canExecute){
            pSol->vGateState.at(i).second = true;
            if(gate.nTargetQubit() == 2){
                --pSol->unexecutedGateCount;
            }
            for(unsigned_t j : _vsGateChild.at(i)){
                sAffectedNode.emplace(j);
            }
            // cerr << "gate " << i << " can be executed" << endl;
            vGateExecuted.emplace_back(i);
        }
        // else if(gate.nTargetQubit() == 2){
        //     pSol->g_cost += _pDevice->getDistance(_vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0)), _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1))); // !
        // }
    }
    for(unsigned_t i : vGateExecuted){
        pSol->sReadyMappedGate.erase(i);
    }
    unordered_set<unsigned_t> sNewAffectedNode;
    while(sAffectedNode.size() > 0){
        vGateExecuted.clear();
        sNewAffectedNode.clear();
        // cerr << "sAffectedNode.size(): " << sAffectedNode.size() << endl;
        for(unsigned_t i : sAffectedNode){
            if(pSol->vGateState.at(i).first && !pSol->vGateState.at(i).second){
                Gate& gate = _pCircuit->gate(i);
                canExecute = (((gate.nTargetQubit() == 1) || _pDevice->isAdjacent(_vCurSolMappingPro2Phy[gate.targetProgramQubit(0)], _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)])));
                if(canExecute){
                    pSol->vGateState.at(i).second = true;
                    if(gate.nTargetQubit() == 2){
                        --pSol->unexecutedGateCount;
                    }
                    for(unsigned_t j : _vsGateChild.at(i)){
                        sNewAffectedNode.emplace(j);
                    }
                    vGateExecuted.emplace_back(i);
                }   
                // else if(gate.nTargetQubit() == 2){
                //     pSol->g_cost += _pDevice->getDistance(_vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0)), _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1))); // !
                // }
                // cerr << "gate " << i << " canexecute: " << canExecute << endl;
            }
            // gate in G_unmmapped
            else if(!pSol->vGateState.at(i).first && pSol->vGateState.at(i).second){
                isReady = true;
                for(unsigned_t j : _vsGateParent.at(i)){
                    // if the parent of the gate is not in 11: !(xx)
                    if(!pSol->vGateState.at(j).first || !pSol->vGateState.at(j).second){
                        isReady = false;
                        break;
                    }
                }
                if(isReady){
                    pSol->vGateState.at(i).first = true;
                    pSol->vGateState.at(i).second = false;
                    sNewAffectedNode.emplace(i);
                    for(unsigned_t j : _vsGateChild.at(i)){
                        sNewAffectedNode.emplace(j);
                    }
                    pSol->sReadyMappedGate.insert(i);
                    --pSol->oneLevelDependencyGateCount;
                }
                // cerr << "gate " << i << " isReady: " << isReady << endl;
            }
            // gate not in G_unmapped
            else if(!pSol->vGateState.at(i).first && !pSol->vGateState.at(i).second){
                hasOneLevelDependency = true;
                for(unsigned_t j : _vsGateParent.at(i)){
                    // if the parent of the gate is in 00 or 01: !x!x || !xx = 
                    if(!pSol->vGateState.at(j).first){
                        hasOneLevelDependency = false;
                        break;
                    }
                }
                if(hasOneLevelDependency){
                    pSol->vGateState.at(i).second = true;
                    sNewAffectedNode.emplace(i);
                    for(unsigned_t j : _vsGateChild.at(i)){
                        sNewAffectedNode.emplace(j);
                    }
                    ++pSol->oneLevelDependencyGateCount;
                }
                // cerr << "gate " << i << " hasOneLevelDependency: " << hasOneLevelDependency << endl;
            }
        }
        sAffectedNode.clear();
        // cerr << "sNewAffectedNode.size(): " << sNewAffectedNode.size() << endl;
        for(unsigned_t i : sNewAffectedNode){
            sAffectedNode.emplace(i);
        }
        for(unsigned_t i : vGateExecuted){
            pSol->sReadyMappedGate.erase(i);
        }
    }

    #ifdef DEBUG
        // ! for check bug
        unsigned_t count = 0;
        for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
            if(pSol->vGateState.at(i).first && pSol->vGateState.at(i).second){
                ++count;
            }
        }
        cerr << "count: " << count << endl;
        cerr << "pSol->unexecutedGateCount: " << pSol->unexecutedGateCount << endl;
        printSolution(pSol);
        assert((_pCircuit->nGate()-count) == pSol->unexecutedGateCount);
    #endif
}

double_t aRouter_dev::calHeuristicCost(shared_ptr<aSolution> pSol){
    // return (static_cast<double_t>(pSol->unexecutedGateCount));
    double_t cost = 0, costDis = 0;
    double_t costUnmap = 0, costRelated = 0;
    unsigned_t p0, p1;
    unsigned_t count = 0, p2, p3;
    vector<bool> vGateIsCount(_pCircuit->nGate(), 0);
    for(unsigned_t i : pSol->sReadyMappedGate){
        Gate& gate = _pCircuit->gate(i);
        p0 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(0)];
        p1 = _vCurSolMappingPro2Phy[gate.targetProgramQubit(1)];
        costDis += _pDevice->getDistance(p0, p1);
        for(unsigned_t j : _vsGateChild[i]){
            Gate& gateChild = _pCircuit->gate(j);
            if(!pSol->vGateState.at(j).first && pSol->vGateState.at(j).second & !vGateIsCount[j]){
                vGateIsCount[j] = 1;
                if(gateChild.nTargetQubit() > 1){
                    p2 = _vCurSolMappingPro2Phy.at(gateChild.targetProgramQubit(0));
                    p3 = _vCurSolMappingPro2Phy.at(gateChild.targetProgramQubit(1));
                    costUnmap += _pDevice->getDistance(p2, p3);
                    costRelated += _pDevice->getDistance(p0, p2);
                    costRelated += _pDevice->getDistance(p0, p3);
                    costRelated += _pDevice->getDistance(p1, p2);
                    costRelated += _pDevice->getDistance(p1, p3);
                    // if(gate.targetProgramQubit(0) != gateChild.targetProgramQubit(0)){
                    // }
                    // else if(gate.targetProgramQubit(0) != gateChild.targetProgramQubit(1)){
                    //     costUnmap += _pDevice->getDistance(p0, p3);
                    // }
                    // else if(gate.targetProgramQubit(1) != gateChild.targetProgramQubit(0)){
                    //     costUnmap += _pDevice->getDistance(p1, p2);
                    // }
                    // else{
                    //     costUnmap += _pDevice->getDistance(p1, p3);
                    // }
                    ++count;
                }
            }
        }
    }
    if(count > 0){
        costUnmap /= 2;
        costRelated /= 3;
        costUnmap += costRelated;
        costUnmap /= (count * _pCircuit->nProgramQubit());
    }
    if(pSol->sReadyMappedGate.size() > 0){
        cost = costUnmap + costDis / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + _aRouter_devParam.weight * (static_cast<double_t>(pSol->unexecutedGateCount));
        // cost /= 2;
        // cost = costUnmap + costDis / (double_t)(pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit()) + (static_cast<double_t>(pSol->unexecutedGateCount)/_pCircuit->nGate());
        // cost = costUnmap + costDis + (static_cast<double_t>(pSol->unexecutedGateCount) * _pCircuit->nProgramQubit());
    }
    // if(cost > 10000000){
    //     // for debug purpose
    //     // for(unsigned_t i = 0 ; i < vGateIsCount.size(); ++i){
    //     //     if(vGateIsCount[i]){
    //     //         Gate& gate = _pCircuit->gate(i);
    //     //         p0 = _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0));
    //     //         p1 = _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1));
    //     //         cerr << "gate distance for " << i << " : " << p0 << ", " << p1 << " dis: " << _pDevice->getDistance(p0, p1) << endl;
    //     //     }
    //     // }
    //     cerr << "cost: " << cost << endl;
    //     cerr << "costDis: " << costDis << endl;
    //     cerr << "pSol->sReadyMappedGate: " << pSol->sReadyMappedGate.size() << endl;
    //     cerr << "pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit(): " << pSol->sReadyMappedGate.size()*_pCircuit->nProgramQubit() << endl;
    //     // cerr << "costUnmap: " << costUnmap << endl;
    //     // cerr << "count: " << count << endl;
    //     cerr << "pSol->unexecutedGateCount: " << (static_cast<double_t>(pSol->unexecutedGateCount)) << endl;
    //     // cerr << "_aRouter_devParam.weight * (static_cast<double_t>(pSol->unexecutedGateCount)): " << _aRouter_devParam.weight * (static_cast<double_t>(pSol->unexecutedGateCount)) << endl;
    // }
    return cost;
}

void aRouter_dev::updateCurMappingForSol(shared_ptr<aSolution> pSol){
    // cerr << "in updateCurMappingForSol, for sol " << pSol << endl;
    for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
        _vCurSolMappingPhy2Pro.at(i) = -1;
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy.at(i) = _vCurMapping.at(i);
        _vCurSolMappingPhy2Pro.at(_vCurMapping.at(i)) = i;
    }
    vector<shared_ptr<aSolution>> vpSol;
    if(_mvSwapEdge.find(pSol) != _mvSwapEdge.end()){
        for(int_t i : _mvSwapEdge[pSol]){
            // cerr << i << " ";
            swapQubit(i);
        }
    }
    else{
        vpSol.emplace_back(pSol);
        while(_mvSwapEdge.find(pSol->pParent) == _mvSwapEdge.end() && pSol->pParent != nullptr){
            vpSol.emplace_back(pSol->pParent);
            pSol = pSol->pParent;
        }
        // cerr << "add swap: ";
        if(_mvSwapEdge.size() > 0){
            if(_mvSwapEdge.find(pSol->pParent) != _mvSwapEdge.end()){
                pSol = pSol->pParent;
            }
            for(int_t i : _mvSwapEdge[pSol]){
                // cerr << i << " ";
                swapQubit(i);
            }
        }
    }
    for(int_t i = vpSol.size() - 1; i >= 0; --i){
        swapQubit(vpSol.at(i)->swapEdge);
        // cerr << vpSol.at(i)->swapEdge << " ";
    }
    #ifdef DEBUG
        cerr << endl;
        cerr << "print _vCurSolMappingPro2Phy" << endl;
        for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
            cerr << i << ": " << _vCurSolMappingPro2Phy.at(i) << endl;
        }
        cerr << "print _vCurSolMappingPhy2Pro" << endl;
        for(unsigned_t i = 0; i < _pDevice->nQubit(); ++i){
            cerr << i << ": " << _vCurSolMappingPhy2Pro.at(i) << endl;
        }
    #endif
}

void aRouter_dev::collectNeighbor(shared_ptr<aSolution> pSol){
    // 
    unsigned_t q0, q1, p0, p1, p2, disBound, disBoundBias = 0;
    // aSolution& sol = _vAllSolution[idx];
    unordered_set<unsigned_t> sEdgeToHaveSwap;
    sEdgeToHaveSwap.clear();
    bool expandNode;
    while(sEdgeToHaveSwap.empty()){
        for(unsigned_t i : pSol->sReadyMappedGate){
            // gate in G_ready
            Gate & gate = _pCircuit->gate(i);
            if(gate.nTargetQubit() > 1){
                q0 = gate.targetProgramQubit(0);
                q1 = gate.targetProgramQubit(1);
                p0 = _vCurSolMappingPro2Phy.at(q0);
                p1 = _vCurSolMappingPro2Phy.at(q1);
                disBound =  _pDevice->getDistance(p1,p0) + disBoundBias;
                // cerr << "q0: " << q0 << endl;
                // cerr << "q1: " << q1 << endl;
                // cerr << "p0: " << p0 << endl;
                // cerr << "p1: " << p1 << endl;
                // cerr << "disBound: " << disBound << endl;
                // change p0 to p2
                unordered_set<int_t>& qubitRegion0 = _vsQubitRegion[q0]; //_pCircuit->sQubitRegion(q0);
                unordered_set<int_t>& qubitRegion1 = _vsQubitRegion[q1];
                // cerr << "for p0" << endl;
                for (unsigned_t eId : _pDevice->qubit(p0).vSpanEdge){
                    //  enable routing region
                    p2 = (_pDevice->edge(eId).qubitId1() == p0) ? _pDevice->edge(eId).qubitId2() : _pDevice->edge(eId).qubitId1();
                    
                    if(_aRouter_devParam.restrict_region){
                        expandNode = false;
                        // p0 is in q0's qubit refion
                        if(qubitRegion0.count(p0) > 0){
                            if(qubitRegion0.count(p2) > 0){
                                expandNode = true;
                            }
                            else{
                                expandNode = (((double_t)(rand() % 1000) / 1000.0) < _aRouter_devParam.expand_prob);
                            }
                        }
                        else{
                            expandNode = true;
                        }
                    }
                    else{
                        expandNode = true;
                    }
                    expandNode = (expandNode && (disBound - _pDevice->getDistance(p1, p2) >= 1e-9));
                    // cerr << "try p2 " << p2  << ", expand=" << expandNode << endl;
                    if(expandNode){
                        sEdgeToHaveSwap.insert(eId);
                    }
                }
                // change p1 to p2
                // cerr << "for p1" << endl;
                for (unsigned_t eId : _pDevice->qubit(p1).vSpanEdge){
                    // ! not enable routing region
                    p2 = (_pDevice->edge(eId).qubitId1() == p1) ? _pDevice->edge(eId).qubitId2() : _pDevice->edge(eId).qubitId1();
                    if(_aRouter_devParam.restrict_region){
                        expandNode = false;
                        if(qubitRegion1.count(p1) > 0){
                            if(qubitRegion1.count(p2) > 0){
                                expandNode = true;
                            }
                            else{
                                expandNode = (((double_t)(rand() % 1000) / 1000.0) < _aRouter_devParam.expand_prob);
                            }
                        }
                        else{
                            expandNode = true;
                        }
                    }
                    else{
                        expandNode = true;
                    }
                    expandNode = (expandNode && (disBound - _pDevice->getDistance(p0, p2) >= 1e-9));
                    // cerr << "try p2 " << p2  << ", expand=" << expandNode << endl;
                    if(expandNode){
                        sEdgeToHaveSwap.insert(eId);
                    }
                }
            }
        }
        disBoundBias += 1;
    }
    for (unsigned_t eId : sEdgeToHaveSwap){
        // _vpAllSolution.emplace_back(pAdjSol);
        _vAllSolution.emplace_back(make_shared<aSolution>(pSol->sReadyMappedGate, pSol->vGateState, pSol->unexecutedGateCount, eId));
        shared_ptr<aSolution> pAdjSol = _vAllSolution.back();
        pSol->vAdjs.emplace_back(pAdjSol);
        pAdjSol->g_cost = pSol->g_cost + _pDevice->getDistance(_pDevice->edge(eId).qubitId1(), _pDevice->edge(eId).qubitId2()) + 1;
        pAdjSol->done = false;
        // insert swap gate
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

void aRouter_dev::updateCircuit(shared_ptr<aSolution> pSol){
    // cerr << "final sol: " << pSol << endl;
    // backtrack pSol
    // printSolution(pSol);
    vector<shared_ptr<aSolution>> vpSol;
    if(_mvSwapEdge.find(pSol) == _mvSwapEdge.end()){
        vpSol.emplace_back(pSol);
        while(_mvSwapEdge.find(pSol->pParent) == _mvSwapEdge.end() && pSol->pParent != nullptr){
            vpSol.emplace_back(pSol->pParent);
            pSol = pSol->pParent;
            // printSolution(pSol);
            // getchar();
        }
        pSol = pSol->pParent;
    }
    
    // collect the circuit initial mapping
    _vCurMapping.clear();
    _vCurMapping.resize(_pCircuit->nProgramQubit(), 0);
    _vCurSolMappingPro2Phy.clear();
    _vCurSolMappingPro2Phy.resize(_pCircuit->nProgramQubit(), 0);
    _vCurSolMappingPhy2Pro.clear();
    _vCurSolMappingPhy2Pro.resize(_pDevice->nQubit(), -1);
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurMapping.at(i) = _pCircuit->initialMapping(i);
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurSolMappingPro2Phy.at(i) = _vCurMapping.at(i);
        _vCurSolMappingPhy2Pro.at(_vCurMapping.at(i)) = i;
    }
    vector<bool> vExecutedGate(_pCircuit->nGate(), 0);
    vector<unsigned_t> vQubitLastTime(_pDevice->nQubit(), 0);
    // collect a set of gates that are ready to be executed 
    _sGateReadyToMap.clear();
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        if(_vsGateParent.at(i).size() == 0){
            _sGateReadyToMap.insert(i);
        }
    }
    unsigned_t t;
    executeGateUnderCurrentMapping(vExecutedGate, vQubitLastTime);
    // check each layer for _mvSwapEdge
    if(pSol != nullptr){
        for(int_t i : _mvSwapEdge[pSol]){
            // cerr << i << endl;
            swapQubit(i);
            // add swap gate
            if(i >= 0){
                Edge & edge = _pDevice->edge(i);
                t = vQubitLastTime.at(edge.qubitId1()) < vQubitLastTime.at(edge.qubitId2()) ? vQubitLastTime.at(edge.qubitId2()) : vQubitLastTime.at(edge.qubitId1()); 
                // cerr << "time " << t << " a swap on edge " << i << ", " << edge.qubitId1() << " " << edge.qubitId2() << "(" << proQ0 << ", " << proQ1 << ")"<< endl;
                addSwapGate(i, t);
                if(t > 0){
                    ++t;
                    vQubitLastTime.at(edge.qubitId1()) = t;
                    vQubitLastTime.at(edge.qubitId2()) = t;
                }
            }
            // collect executable gates under current mapping
            executeGateUnderCurrentMapping(vExecutedGate, vQubitLastTime);
        }
    }
    // check each layer for the final solution
    for(int_t i = vpSol.size() - 1; i >= 0; --i){
        // cerr << vpSol.at(i)->swapEdge << endl;
        swapQubit(vpSol.at(i)->swapEdge);
        // add swap gate
        if(vpSol.at(i)->swapEdge >= 0){
            Edge & edge = _pDevice->edge(vpSol.at(i)->swapEdge);
            t = vQubitLastTime.at(edge.qubitId1()) < vQubitLastTime.at(edge.qubitId2()) ? vQubitLastTime.at(edge.qubitId2()) : vQubitLastTime.at(edge.qubitId1()); 
            // cerr << "time " << t << " a swap on edge " << vpSol.at(i)->swapEdge << ", " << edge.qubitId1() << " " << edge.qubitId2() << "(" << proQ0 << ", " << proQ1 << ")"<< endl;
            addSwapGate(vpSol.at(i)->swapEdge, t);
            if(t > 0){
                ++t;
                vQubitLastTime.at(edge.qubitId1()) = t;
                vQubitLastTime.at(edge.qubitId2()) = t;
            }
        }
        // collect executable gates under current mapping
        executeGateUnderCurrentMapping(vExecutedGate, vQubitLastTime);
    }
    // set circuit depth
    unsigned_t depth = 1;
    // for(unsigned_t i=0; i<_pCircuit->nProgramQubit(); i++) {
    for(unsigned_t d : vQubitLastTime) {
        if(depth < d){
            depth = d;
        }
    }
    _pCircuit->setCircuitDepth(depth);
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
        _vCurMapping.at(i) = _vCurSolMappingPro2Phy.at(i);
    }
    for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); i++){
        _pCircuit->setFinalMapping(i, _vCurSolMappingPro2Phy.at(i));
        _pCircuit->addQubitRegion(i, _vCurSolMappingPro2Phy.at(i));
    }

    // fprintf(stderr, "[INFO] get a solution with %d SWAP and depth %d\n", _pCircuit->nSwapGate(), _pCircuit->circuitDepth());
}
void aRouter_dev::addSwapGate(unsigned_t edgeId, unsigned_t t ){
    Edge & edge = _pDevice->edge(edgeId);
    // cerr << "edge: " << edgeId << ", " << edge.qubitId1() << " " << edge.qubitId2() << endl;
    if(t > 0){
        _pCircuit->addSwapGate(_pCircuit->nSwapGate(), edge.qubitId1(), edge.qubitId2(), t);
    }
    else{
        // cerr << "change initial mapping for ";
        int_t proQ0 = -1, proQ1 = -1;
        for(unsigned_t i = 0; i < _pCircuit->nProgramQubit(); ++i){
            if(_pCircuit->initialMapping(i) == static_cast<int>(edge.qubitId1())){
                proQ0 = i;
            }
            else if(_pCircuit->initialMapping(i) == static_cast<int>(edge.qubitId2())){
                proQ1 = i;
            }
        }
        if(proQ0 > -1){
            // cerr << proQ0;
            _pCircuit->setInitialMapping(proQ0, edge.qubitId2());
            _pCircuit->addQubitRegion(proQ0, edge.qubitId2());
        }
        if(proQ1 > -1){
            // cerr << " " << proQ1 << endl;
            _pCircuit->setInitialMapping(proQ1, edge.qubitId1());
            _pCircuit->addQubitRegion(proQ1, edge.qubitId1());
        }
        // cerr << endl;
    }
}

//void aRouter_dev::executeGateUnderCurrentMapping(vector<bool> & vExecutedGate, vector<unsigned_t> & vQubitLastTime){
void aRouter_dev::executeGateUnderCurrentMapping(vector<bool> & vExecutedGate, vector<unsigned_t> & vQubitLastTime){
    // return;
    unsigned_t q0, q1;
    vector<unsigned_t> vExecutedGateInTheRun;
    bool keepFindNewGate = true, canExecute;
    unordered_set<unsigned_t> sNewGateReadyToMap;
    while(keepFindNewGate){
        keepFindNewGate = false;
        sNewGateReadyToMap.clear();
        vExecutedGateInTheRun.clear();
        for(unsigned_t i : _sGateReadyToMap){
            Gate& gate = _pCircuit->gate(i);
            canExecute = false;
            if(!vExecutedGate.at(i)){
                if(gate.nTargetQubit() == 1){
                    canExecute = true;
                }
                else if(_pDevice->isAdjacent(_vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0)), _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1)))){
                    canExecute = true;
                }
            }
            if(canExecute){
                vExecutedGateInTheRun.emplace_back(i);
                vExecutedGate.at(i) = true;
                ++_executedGateCount;
                for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
                    gate.setTargetPhysicalQubit(j, _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(j)));
                    _pCircuit->addQubitRegion(gate.targetProgramQubit(j), _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(j)));
                }
                q0 = _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0));
                if( gate.nTargetQubit() == 1){
                    gate.setExecutionTime(vQubitLastTime.at(q0));
                    ++vQubitLastTime.at(q0);
                }
                else{
                    q1 = _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1));
                    unsigned_t gateTime = (vQubitLastTime.at(q0) < vQubitLastTime.at(q1)) ? vQubitLastTime.at(q1): vQubitLastTime.at(q0);
                    gate.setExecutionTime(gateTime);
                    assert(q0 < _pDevice->nQubit() && 0 <= q0);
                    assert(q1 < _pDevice->nQubit() && 0 <= q1);
                    ++gateTime;
                    vQubitLastTime.at(q0) = gateTime;
                    vQubitLastTime.at(q1) = gateTime;
                    // cerr << &_vQubitLastTime.at(q0) << " " << &_vQubitLastTime.at(q1) << endl;
                }
                for(unsigned_t j : _vsGateChild.at(i)){
                    _vsGateParent.at(j).erase(i);
                    if(_vsGateParent.at(j).size() == 0 && !vExecutedGate.at(j)){
                        sNewGateReadyToMap.insert(j);
                        keepFindNewGate = true;
                    }
                }
                #ifdef DEBUG
                    cerr << "gate " << i << " can be directly executed on program qubit " << gate.targetProgramQubit(0) << "("<< gate.targetPhysicalQubit(0) << ")";
                    if(gate.nTargetQubit() == 2){
                        cerr << ", " << gate.targetProgramQubit(1) << "("<< gate.targetPhysicalQubit(1) << ") on time " << gate.executionTime() << endl; 
                    }
                    else{
                        cerr << " on time " << gate.executionTime() << endl;
                    }
                #endif
            }
        }
        for(unsigned_t i : vExecutedGateInTheRun){
            _sGateReadyToMap.erase(i);
        }
    for(unsigned_t i : sNewGateReadyToMap){
            _sGateReadyToMap.insert(i);
        }
    }
    // cerr << "finish executeGateUnderCurrentMapping" << endl;
}

void aRouter_dev::printSolution(shared_ptr<aSolution> pSol){
    aSolution & sol = *pSol;
    cerr << "Print solution " << pSol << endl;
    cerr << "neighbor size: " << pSol->vAdjs.size() << endl;
    cerr << "g cost: " << sol.g_cost << ", f cost: " << sol.f_cost << endl;
    if(sol.swapEdge > -1){
        cerr << "swapEdge: " << sol.swapEdge << " (" << _pDevice->edge(sol.swapEdge).qubitId1() << ", " << _pDevice->edge(sol.swapEdge).qubitId2() << ")" << endl;
    }
    cerr << "pParent: " << sol.pParent << endl;
    unordered_set<unsigned_t> sGateIsExectued;
    unordered_set<unsigned_t> sGateIsReady;
    unordered_set<unsigned_t> sGateIsNotReady;
    for(unsigned_t i = 0; i < _pCircuit->nGate(); ++i){
        // gate in G_ready
        if(pSol->vGateState.at(i).first && !pSol->vGateState.at(i).second){
            sGateIsReady.insert(i);
        }
        // gate in G_unmmapped
        else if(!pSol->vGateState.at(i).first && pSol->vGateState.at(i).second){
            sGateIsNotReady.insert(i);
        }
        // gate not in G_unmapped
        else if(pSol->vGateState.at(i).first && pSol->vGateState.at(i).second){
            sGateIsExectued.insert(i);
        }
    }
    cerr << "unexecutedGateCount: " << pSol->unexecutedGateCount << endl;
    cerr << "sGateIsExectued: (" << sGateIsExectued.size() << "):";
    // getchar();
    for(unsigned_t i : sGateIsExectued){
        cerr << " " << i;        
    }
    cerr << endl;
    assert(pSol->unexecutedGateCount == static_cast<int_t>(_pCircuit->nGate() - sGateIsExectued.size()));
    cerr << "sGateIsReady: (" << sGateIsReady.size() << "): "  << pSol->sReadyMappedGate.size() << ": ";
    // getchar();
    for(unsigned_t i : sGateIsReady){
        Gate & gate = _pCircuit->gate(i);
        cerr << " " << i << " (" << gate.targetProgramQubit(0) << ", " << gate.targetProgramQubit(1) << ")" << "(" << _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(0)) << ", " << _vCurSolMappingPro2Phy.at(gate.targetProgramQubit(1)) << ")";
    }
    // getchar();
    cerr << endl;
    cerr << "sGateIsNotReady: (" << sGateIsNotReady.size() << ")";
    // getchar();
    for(unsigned_t i : sGateIsNotReady){
        cerr << " " << i;
    }
    // getchar();
    cerr << endl;
}

void aRouter_dev::cleanAllNode(){
    _vAllSolution.clear();
}

bool aRouter_dev::checkMemUsage(){
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
    return false;
}

MOLSQ_NAMESPACE_CPP_END

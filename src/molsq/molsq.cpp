/***********************************************************************
  File        [ molsq.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ planner ]
  Synopsis    [ planner class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2022 ]
***********************************************************************/

#include "molsq/molsq.hpp"

MOLSQ_NAMESPACE_CPP_START


void mOLSQ::run(){
    if(_device.nQubit() < _circuit.nProgramQubit()){
        fprintf(stdout, "[ERROR] #physical qubit (%d) is less than #program qubit (%d)\n", _device.nQubit(), _circuit.nProgramQubit());
        return;
    }
    _timer.start(TimeUsage::FULL);
    Writer writer(_circuit, _device);
    // preprocessing: construct gate dependency graph
    preprocessing();
    if(!_molsqParam.is_all_commute){
        _molsqParam.initialMapper_threshold_qubit *= 2;
    }
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====      Stage 0: Heuristic Algorithm     ====\n");
    fprintf(stdout, "===============================================\n");
    Circuit reverseCir("reverseCir", _circuit.nProgramQubit(), _circuit.nGate());
    constructReverseCircuit(_circuit, reverseCir);
    flatNonSMT(reverseCir);
    // runVCyclePureNonSMT();
    fprintf(stdout, "[INFO] get a solution with %d SWAP\n", _circuit.nSwapGate());
    writer.outputQASM(_circuit.name()+"_stage0.qasm");
    // _circuit.printCircuitLayout();
    // _circuit.printCircuitLayout();
    verify(_device, _circuit);
    fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\n", _circuit.nSwapGate(), _circuit.circuitDepth());
    _timer.showUsage("Heuristic algorithm compilation time", TimeUsage::FULL);
    // return;
    if(_circuit.nSwapGate() == 0){
        _timer.showUsage("Total compilation time", TimeUsage::FULL);
        return;
    }
    // if(!_molsqParam.is_all_commute){
    //     return;
    // }
    // check if use multilevel approach
    // if no, directly solve by olsq2
    if(!_molsqParam.use_multilevel){
        OLSQ2 olsq2(_circuit, _device);
        if(_molsqParam.is_all_commute){
            olsq2.enableAllCommute();
        }
        olsq2.setHeuristicForSwap(_circuit.nSwapGate());
        olsq2.setOptimizeForSwap();
        olsq2.run();
        writer.outputQASM(_circuit.name()+"_olsq2.qasm");
    }
    // if yes, use multilevel approach
    else{
        fprintf(stdout, "===============================================\n");
        fprintf(stdout, "====     V Cycle 1: Non-SMT Refinement     ====\n");
        fprintf(stdout, "===============================================\n");
        runVCyclePureNonSMT(reverseCir);
        if(_molsqParam.use_reverse_circuit){
            updateCircuitWithReverse(_circuit, reverseCir);
        }
        asapScheduling(_circuit);
        verify(_device, _circuit);
        fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\n", _circuit.nSwapGate(), _circuit.circuitDepth());
        writer.outputQASM(_circuit.name()+"_stage1.qasm");
        _timer.showUsage("Non-SMT refinement compilation time", TimeUsage::FULL);
        if(_circuit.nSwapGate() == 0){
            _timer.showUsage("Total compilation time", TimeUsage::FULL);
            return;
        }
        return;
        fprintf(stdout, "===============================================\n");
        fprintf(stdout, "====       V Cycle 2: SMT Refinement       ====\n");
        fprintf(stdout, "===============================================\n");
        runVCycleMix(reverseCir);  
        // runVCyclePureNonSMT();      
        // asapScheduling(_circuit);
        verify(_device, _circuit);
        fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\n", _circuit.nSwapGate(), _circuit.circuitDepth());
        _timer.showUsage("SMT refinement compilation time", TimeUsage::FULL);
        writer.outputQASM(_circuit.name()+"_stage2.qasm");
    }
    // _circuit.printCircuitLayout();
    // verify(_device, _circuit);
    _timer.showUsage("Total compilation time", TimeUsage::FULL);
}

void mOLSQ::preprocessing(){
    _circuit.constructDependency();
    _device.calAPSP();
}

bool mOLSQ::useMultilevel(){
    // TODO
    return true;
}

void mOLSQ::runCoarestLevel(Circuit & cir, Device & device){
    // for abalation study
    // if(1){
    OLSQ2 olsq2(cir, device);
    if(_molsqParam.is_all_commute){
        olsq2.enableAllCommute();
    }
    olsq2.setMutilevelMode();
    bool solve = olsq2.run();
    if(!solve){
        InitialMapper initialMapper;
        bool result = initialMapper.run(cir, device);
        if(result){
            return;
        }
        Placer_dev placer;
        if(_molsqParam.placer_one_hop_cost){
            placer.enableOneHopCost();
        }
        else{
            placer.disableOneHopCost();
        }
        placer.enableAllCommute();
        placer.run(cir, device, false);
        aRouter_dev arouter;
        arouter.enableAllCommute();
        aRouterThreePass(arouter, device, cir, true);
    }
    return;

    // if(cir.nProgramQubit() > 16 || cir.nGate() > 20 || device.nQubit() > 30){
    //     InitialMapper initialMapper(cir, device);
    //     bool result = initialMapper.run();
    //     if(result){
    //         return;
    //     }
    //     Placer_dev placer;
    //     if(_molsqParam.placer_one_hop_cost){
    //         placer.enableOneHopCost();
    //     }
    //     else{
    //         placer.disableOneHopCost();
    //     }
    //     placer.enableAllCommute();
    //     // if(_molsqParam.is_all_commute){
    //     //     placer.enableAllCommute();
    //     // }
    //     // else{
    //     //     placer.disableAllCommute();
    //     // }
    //     placer.run(cir, device, false);
    //     aRouter_dev arouter;
    //     arouter.enableAllCommute();
    //     aRouterThreePass(arouter, device, cir, true);
    // }
    // else{
    //     cir.printCircuit();
    //     device.printDevice();
    //     OLSQ2 olsq2(cir, device);
    //     // if(_molsqParam.is_all_commute){
    //     olsq2.enableAllCommute();
    //     // }
    //     olsq2.setMutilevelMode();
    //     bool solve = olsq2.run();
    // }
}

void mOLSQ::runVCyclePureNonSMT(Circuit & reverseCir){
    // calculate how many level we need
    int_t proLevel = ceil(log2(_circuit.nProgramQubit() + 1)) - 2;
    // int_t proLevel = 0;
    int_t phyLevel = proLevel;
    _timer.setTimeout(3600 * 8);
    // a vector of devices and circuits
    vector<Device> vCoarseDevice;
    vector<Circuit> vCoarseCir;
    _vmvCoarsePhyQ2FinerPhyQ.clear();
    _vmvCoarseProQ2FinerProQ.clear();

    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====          Stage 1: Clutsering          ====\n");
    fprintf(stdout, "===============================================\n");
    for(unsigned_t i = 0; i <= proLevel; ++i){
        vCoarseDevice.emplace_back(_device.name()+"_level_"+to_string(i));
        _vmvCoarsePhyQ2FinerPhyQ.emplace_back();
        vCoarseCir.emplace_back(_circuit.name()+"_level_"+to_string(i));
        _vmvCoarseProQ2FinerProQ.emplace_back();
    }
    Clusterer3 clusterer(0);
    // Clusterer2 clusterer(0);

    // obtain circuit mapping to do clustering
    // InitialMapper initialMapper(_circuit, _device);
    // if(_molsqParam.is_all_commute){
    //     initialMapper.enableAllCommute();
    // }
    // else{
    //     initialMapper.disableAllCommute();
    // }
    // bool result = initialMapper.run();
    // if(result){
    //     // asapScheduling(_circuit);
    //     // _circuit.printCircuitLayout();
    //     // verify(_device, _circuit);
    //     // _timer.showUsage("Total compilation time", TimeUsage::FULL);
    //     _circuit.clearSwap();
    //     asapScheduling(_circuit);
    //     return;
    // }
    if(_molsqParam.use_reverse_circuit){
        clusterer.cluster(reverseCir, _device, vCoarseCir, vCoarseDevice, _vmvCoarseProQ2FinerProQ, _vmvCoarsePhyQ2FinerPhyQ, _molsqParam.is_all_commute);
    }
    else{
        clusterer.cluster(_circuit, _device, vCoarseCir, vCoarseDevice, _vmvCoarseProQ2FinerProQ, _vmvCoarsePhyQ2FinerPhyQ, _molsqParam.is_all_commute);
    }
    
    // for(unsigned_t i = 0; i <= proLevel; ++i){
    //     if(!_molsqParam.is_all_commute){
    //         vCoarseCir[i].constructDependency();
    //     }
    // }

    // return;
    // begin with proper level
    // cerr << "proLevel: " << proLevel << endl;
    // cerr << "vCoarseCir.size(): " << vCoarseCir.size() << endl;
    while(proLevel > 0 && phyLevel > 0){
        if(vCoarseCir.at(proLevel).nProgramQubit() < _molsqParam.coarse_level_qubit && vCoarseCir.at(proLevel).nGate() < _molsqParam.coarse_level_gate){
            break;
        }
        --proLevel;
        --phyLevel;
        // cerr << "proLevel: " << proLevel << endl;
    }
    // cerr << "proLevel: " << proLevel << endl;
    while(proLevel > 0 && (vCoarseCir.at(proLevel).nProgramQubit() < 3 ||  vCoarseCir.at(proLevel).nGate() < 10)){
        --proLevel;
        --phyLevel;
    }
    // cerr << "line 193" << endl;
    // compile coarsest level circuit by OLSQ2
    fprintf(stdout, "[Info] %d Levels in V Cycle\n", proLevel);
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====   Stage 2: Coarsest level solving     ====\n");
    fprintf(stdout, "===============================================\n");
    // runCoarestLevel(vCoarseCir.at(proLevel), vCoarseDevice.at(phyLevel));
    while(proLevel >= 0 && vCoarseCir.at(proLevel).nProgramQubit() < 30){
        runCoarestLevel(vCoarseCir.at(proLevel), vCoarseDevice.at(phyLevel));
        --phyLevel;
        --proLevel;
        if(proLevel >= 0 && vCoarseCir.at(proLevel+1).nSwapGate() > 3){
            break;
        }
    }
    ++phyLevel;
    ++proLevel;
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====         Stage 3: Refinement           ====\n");
    fprintf(stdout, "===============================================\n");
    nonsmtRefinment(phyLevel, proLevel, vCoarseCir, vCoarseDevice, reverseCir);
    return;
}

void mOLSQ::runVCycleMix(Circuit & reverseCir){
    // calculate how many level we need
    int_t proLevel = ceil(log2(_circuit.nProgramQubit() + 1)) - 2;
    int_t phyLevel = proLevel;

    // a vector of devices and circuits
    vector<Device> vCoarseDevice;
    vector<Circuit> vCoarseCir;
    _vmvCoarsePhyQ2FinerPhyQ.clear();
    _vmvCoarseProQ2FinerProQ.clear();

    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====          Stage 1: Clutsering          ====\n");
    fprintf(stdout, "===============================================\n");
    for(unsigned_t i = 0; i <= proLevel; ++i){
        vCoarseDevice.emplace_back(_device.name()+"_level_"+to_string(i));
        _vmvCoarsePhyQ2FinerPhyQ.emplace_back();
        vCoarseCir.emplace_back(_circuit.name()+"_level_"+to_string(i));
        _vmvCoarseProQ2FinerProQ.emplace_back();
    }
    Clusterer3 clusterer(0);

    // obtain circuit mapping to do clustering
    // InitialMapper initialMapper(_circuit, _device);
    // if(_molsqParam.is_all_commute){
    //     initialMapper.enableAllCommute();
    // }
    // else{
    //     initialMapper.disableAllCommute();
    // }
    // bool result = initialMapper.run();
    // if(result){
    //     _circuit.clearSwap();
    //     asapScheduling(_circuit);
    //     return;
    // }
    clusterer.cluster(_circuit, _device, vCoarseCir, vCoarseDevice, _vmvCoarseProQ2FinerProQ, _vmvCoarsePhyQ2FinerPhyQ, _molsqParam.is_all_commute);
    
    for(unsigned_t i = 0; i <= proLevel; ++i){
        if(!_molsqParam.is_all_commute){
            vCoarseCir[i].constructDependency();
        }
    }
    // begin with proper level
    while(proLevel > 0 && phyLevel > 0){
        if(vCoarseCir.at(proLevel).nProgramQubit() < _molsqParam.coarse_level_qubit && vCoarseCir.at(proLevel).nGate() < _molsqParam.coarse_level_gate){
            break;
        }
        --proLevel;
        --phyLevel;
        // cerr << "proLevel: " << proLevel << endl;
    }
    while(vCoarseCir.at(proLevel).nProgramQubit() == 2 ||  vCoarseCir.at(proLevel).nGate() < 10){
        --proLevel;
        --phyLevel;
    }
    // compile coarsest level circuit by OLSQ2
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====   Stage 2: Coarsest level solving     ====\n");
    fprintf(stdout, "===============================================\n");
    runCoarestLevel(vCoarseCir.at(proLevel), vCoarseDevice.at(phyLevel));
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====         Stage 3: Refinement           ====\n");
    fprintf(stdout, "===============================================\n");
    // do one more level compilation to refine the coarest level device clustering
    sTBOLSQ2 stbolsq2;
    if(_molsqParam.is_all_commute){
        stbolsq2.enableAllCommute();
    }
    unsigned_t depth = vCoarseCir.at(proLevel).circuitDepth();
    unsigned_t bound = vCoarseCir.at(proLevel).nSwapGate() * 2;

    // compile the whole circuit in next level by stboslq2
    // iteratively declustering 
    stbolsq2.setChangeCntBound(bound);
    stbolsq2.setDepth(depth);
    // solve by sTBOLSQ2
    assert(phyLevel == proLevel);
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====      Stage 3.1: SMT Refinement        ====\n");
    fprintf(stdout, "===============================================\n");
    smtRefinment(stbolsq2, phyLevel, proLevel, vCoarseCir, vCoarseDevice, bound);
    // map _device and _circuit
    if(proLevel == 0 && _circuit.nProgramQubit() < _molsqParam.smt_refinement_qubit_cnt && _circuit.nGate() < _molsqParam.smt_refinement_gate_cnt){
        finestLevelRefinement(stbolsq2, vCoarseCir[0], vCoarseDevice[0]);
        if(stbolsq2.bestDisSum() < _circuit.nSwapGate()){
            _circuit.setFinalMapping(_vvQubitMapping.back());
            _circuit.setCircuitDepth(_vvQubitMapping.size()+1);
            unsigned_t nSwapGate = stbolsq2.bestDisSum();
            _circuit.setInitialMapping(_vvQubitMapping[0]);
            aRouter_dev arouter;
            arouter.enableRestrictRegion();
            arouter.setQubitRegion(_circuit);
            if(_molsqParam.is_all_commute){
                arouter.enableAllCommute();
            }
            else{
                arouter.disableAllCommute();
            }
            // arouter.setGCostLimit(nSwapGate);
            arouter.setQubitRegion(_circuit);
            bool result = arouter.run(_circuit, _device);
            if(!result){
                glueMapping(_vvQubitMapping.size()-1);
            }
            _circuit.printCircuitLayout();
            asapScheduling(_circuit);
            return;
        }
    }
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====    Stage 3.2: Non-SMT Refinement      ====\n");
    fprintf(stdout, "===============================================\n");
    nonsmtRefinment(phyLevel, proLevel, vCoarseCir, vCoarseDevice, reverseCir);
    return;
}

void mOLSQ::collectQubitRegion(Device& device, Circuit& finerCir, Circuit& coarserCir, bool declusterCir, unsigned_t phyLevel, unsigned_t proLevel){
    // device.printDevice();
    // coarserCir.printCircuit();
    // coarserCir.printQubitRegion();
    // cerr << "in collectQubitRegion" << endl;
    if(declusterCir){
        finerCir.resetQubitRegion();
    }
    for(unsigned_t i = 0; i < coarserCir.nProgramQubit(); ++i){
        // collect circuit qubit in the current device
        unordered_set<int_t> sQubitRegion(coarserCir.sQubitRegion(i));
        assert(sQubitRegion.size());
        // cerr << "generate region for coarse qubit " << i << endl;
        if(_vmvCoarseProQ2FinerProQ.at(proLevel)[i].size() > 2){
            bfsSearch(device, sQubitRegion);
        }
        else{
            bfsSearch(device, sQubitRegion);
        }
        // cerr << "finish bfs" << endl;
        // project the current region to the next level region
        if(!declusterCir){
            coarserCir.resetQubitRegion(i);
        }
        if(declusterCir){
            for(unsigned_t j : _vmvCoarseProQ2FinerProQ.at(proLevel)[i]){
                // cerr << "generate region for finer qubit: " << j << ":";
                for (unordered_set<int_t>::iterator it = sQubitRegion.begin(); it != sQubitRegion.end(); ++it) {
                    for(unsigned_t k : _vmvCoarsePhyQ2FinerPhyQ.at(phyLevel)[*it]){
                        // cerr << " " << k << endl;
                        finerCir.addQubitRegion(j, k);
                    }
                    // cerr << endl;
                }
            }
        }
        else{
            for (unordered_set<int_t>::iterator it = sQubitRegion.begin(); it != sQubitRegion.end(); ++it) {
                for(unsigned_t j : _vmvCoarsePhyQ2FinerPhyQ.at(proLevel)[*it]){
                    coarserCir.addQubitRegion(i, j);
                }
            }
        }
    }
}

void mOLSQ::smtRefinment(sTBOLSQ2& stbolsq2, int_t & phyLevel, int_t & proLevel, vector<Circuit>& vCoarseCir, vector<Device>& vCoarseDevice, unsigned_t bound, bool oneRun){
    while((phyLevel > 0) && (vCoarseCir[proLevel-1].nProgramQubit() < _molsqParam.smt_refinement_qubit_cnt) && (vCoarseCir[proLevel-1].nGate() < _molsqParam.smt_refinement_gate_cnt)){
        fprintf(stdout, "[INFO] Refine level-%d,%d prolem\n", proLevel, phyLevel);
        // cerr << "in normal SMT refinement" << endl;
        // if(_verbose > 0){
        //     vCoarseCir.at(proLevel).printQubitRegion();
        // }
        collectQubitRegion(vCoarseDevice.at(phyLevel), vCoarseCir[proLevel-1], vCoarseCir.at(proLevel), true, phyLevel, proLevel);
        --proLevel;
        --phyLevel;
        if(_verbose > 0){
            vCoarseCir.at(proLevel).printQubitRegion();
        }
        // return;    
        // solve by sTBOLSQ2
        stbolsq2.setCircuit(&vCoarseCir.at(proLevel));
        stbolsq2.setDevice(&vCoarseDevice.at(phyLevel));
        stbolsq2.run(bound, 1);
        bound = min(stbolsq2.changeCntBound() * 2, _circuit.nSwapGate());
        if(oneRun){
            ++phyLevel;
            ++proLevel;
            return;
        }
    }
}


void mOLSQ::finestLevelRefinement(sTBOLSQ2 & stbolsq2, Circuit& cir, Device& device){
    stbolsq2.setDevice(&_device);
    collectQubitRegion(device, _circuit, cir, true, 0, 0);
    if(_verbose > 0){
        _circuit.printQubitRegion();
    }
    // return;    
    // solve by sTBOLSQ2
    stbolsq2.setCircuit(&_circuit);
    stbolsq2.run(stbolsq2.changeCntBound(), 1);
    _vvQubitMapping.clear();
    _vvQubitMapping = stbolsq2.vvQubitMapping();
}

void mOLSQ::nonsmtRefinment(int_t & phyLevel, int_t & proLevel, vector<Circuit>& vCoarseCir, vector<Device>& vCoarseDevice, Circuit & reverseCir){
    // Placer placer;
    Placer_dev placer;
    placer.setRestrictMovement(true);
    if(_molsqParam.placer_one_hop_cost){
        placer.enableOneHopCost();
    }
    else{
        placer.disableOneHopCost();
    }
    if(!_molsqParam.is_all_commute){
        placer.disableAllCommute();
    }
    else{
        placer.enableAllCommute();
    }
    while(phyLevel >= 0){
        Circuit & coarserCir = vCoarseCir.at(proLevel);
        Device & coarserDevice = vCoarseDevice.at(phyLevel);
        --proLevel;
        --phyLevel;
        Circuit & finerCir = (proLevel >= 0)? vCoarseCir.at(proLevel) : ((_molsqParam.use_reverse_circuit)? reverseCir : _circuit);
        Device & finerDevice = (phyLevel >= 0)? vCoarseDevice.at(phyLevel) : _device;
        aRouter_dev arouter;
        if(proLevel >= 0){
            fprintf(stdout, "[INFO] Refine level-%d,%d problem\n", proLevel, phyLevel);
        }
        else{
            fprintf(stdout, "[INFO] Refine the finest-level problem\n", proLevel, phyLevel);
            
        }
        if(!_molsqParam.is_all_commute){
            arouter.disableAllCommute();
        }
        else{
            arouter.enableAllCommute();
        }
        // cerr << "in non SMT refinement" << endl;
        // coarserCir.printQubitRegion();
        collectQubitRegion(coarserDevice, finerCir, coarserCir, true, phyLevel+1, proLevel+1);
        // finerCir.printQubitRegion();
        // cerr << "finish construct qubit region" << endl;
        if(_verbose > 0){
            finerCir.printQubitRegion();
        }
        unsigned_t phyQ;
        my_graph graph(finerCir.nProgramQubit() + finerDevice.nQubit());
        for(unsigned_t i = 0; i < coarserCir.nProgramQubit(); ++i){
            phyQ = coarserCir.initialMapping(i);
            for(unsigned_t j : _vmvCoarseProQ2FinerProQ[proLevel+1][i]){
                for(unsigned_t k : _vmvCoarsePhyQ2FinerPhyQ[phyLevel+1][phyQ]){
                    add_edge(j, finerCir.nProgramQubit() + k, EdgeProperty(1), graph);
                    Qubit & qubit = finerDevice.qubit(k);
                    for(unsigned_t e : qubit.vSpanEdge){
                        Edge & edge = finerDevice.edge(e);
                        if(edge.qubitId1() == k){
                            add_edge(j, finerCir.nProgramQubit() + edge.qubitId2(), EdgeProperty(1), graph);
                        }
                        else{
                            add_edge(j, finerCir.nProgramQubit() + edge.qubitId1(), EdgeProperty(1), graph);
                        }
                    }
                }
            }
        }
        setQubitMappingViaMatching(graph, finerCir, finerDevice.nQubit());
        placer.run(finerCir, finerDevice, true);
        
        
        // arouter.enableRestrictRegion();
        // arouter.setQubitRegion(finerCir);
        Circuit cir(finerCir);
        // if(proLevel < 0){
        //     arouter.setGCostLimit(_circuit.nSwapGate());
        // }
        aRouterThreePass(arouter, finerDevice, finerCir, false);
        fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\r", finerCir.nSwapGate(), finerCir.circuitDepth());

        // return;
        arouter.enableKillIfTimeout();
        
        for(unsigned_t c = 1; c < _molsqParam.trials + 2; ++c){
            // find a new matching
            for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
                if(rand() % 2){
                    remove_edge(i, cir.nProgramQubit() + cir.initialMapping(i), graph);
                }
            }
            setQubitMappingViaMatching(graph, cir, finerDevice.nQubit());
            placer.setProbMoveWithinRegion(0.6 - c * 0.1);
            arouter.setExpandProb(0.15 * c);
            placer.run(cir, finerDevice, true);
            aRouterThreePass(arouter, finerDevice, cir, false);
            if((cir.nSwapGate() < finerCir.nSwapGate()) || (cir.nSwapGate() == finerCir.nSwapGate() && cir.circuitDepth() < finerCir.circuitDepth())){
                finerCir = cir;
            }
            // bool hasResult = arouter.run(cir, finerDevice);
            // if (hasResult){
            //     fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\r", cir.nSwapGate(), cir.circuitDepth());
            //     if(cir.nSwapGate() < finerCir.nSwapGate()){
            //         for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
            //             finerCir.setInitialMapping(i, cir.initialMapping(i));
            //             finerCir.setFinalMapping(i, cir.finalMapping(i));
            //         }
            //         for(unsigned_t i = 0; i < finerCir.nGate(); ++i){
            //             Gate & gate = finerCir.gate(i);
            //             Gate & gateForward = cir.gate(i);
            //             for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
            //                 gate.setTargetPhysicalQubit(j, gateForward.targetPhysicalQubit(j));
            //             }
            //             gate.setExecutionTime(gateForward.executionTime());
            //         }
            //         finerCir.clearSwap();
            //         for(unsigned_t i = 0; i < cir.nSwapGate(); ++i){
            //             Gate & gate = cir.swapGate(i);
            //             finerCir.addSwapGate(i, gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1), gate.executionTime());
            //         }
            //         finerCir.setCircuitDepth(cir.circuitDepth());
            //         finerCir.resetQubitRegion();
            //         for(unsigned_t i = 0; i < finerCir.nProgramQubit(); ++i){
            //             unordered_set<int_t> sQubitRegion = cir.sQubitRegion(i);
            //             for(unsigned_t j : sQubitRegion){
            //                 finerCir.addQubitRegion(i, j);
            //             }
            //         }
            //         // arouter.setGCostLimit(finerCir.nSwapGate());
            //     }
            // }
            // else{
            //     break;
            // }
        }
    }
}

void mOLSQ::setQubitMappingViaMatching(my_graph & graph, Circuit & cir, unsigned_t nQubit){
    vector<V> mate(num_vertices(graph));
    // use boost library to construct maximum card matching 
    // see: https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/maximum_matching.html
    edmonds_maximum_cardinality_matching(graph, &mate[0]);
    unsigned_t pairCount = 0;
    vector<bool> vProIsMapped(cir.nProgramQubit(), 0);
    vector<bool> vPhyIsMapped(nQubit, 0);
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            cir.setInitialMapping(v, mate[v]-cir.nProgramQubit());
            vProIsMapped[v] = 1;
            vPhyIsMapped[mate[v]-cir.nProgramQubit()] = 1;
            ++pairCount;
        }
    }
    if(pairCount < cir.nProgramQubit()){
        for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
            if(!vProIsMapped[i]){
                for(unsigned_t j = 0; j < nQubit; ++j){
                    if(!vPhyIsMapped[j]){
                        cir.setInitialMapping(i,j);
                        vProIsMapped[i] = 1;
                        vPhyIsMapped[j] = 1;
                        break;
                    }
                }
            }
        }
    }
}

void mOLSQ::flatNonSMT(Circuit & reverseCir){
    // Placer placer;
    InitialMapper initialMapper;
    bool result = false, hasRunInitialMapper = false;
    if(_circuit.nProgramQubit() < 200){
        result = initialMapper.run(_circuit, _device);
        hasRunInitialMapper == true;
    }
    if(result){
        _circuit.clearSwap();
        asapScheduling(_circuit);
        return;
    }
    Placer_dev placer;
    aRouter_dev arouter;
    if(_molsqParam.is_all_commute){
        placer.enableAllCommute();
        arouter.enableAllCommute();
    }
    else{
        arouter.disableAllCommute();
        // vector<unsigned_t> vPhyQubit = {64, 72, 74, 66, 68, 56, 30, 38, 24, 6, 4, 21, 3, 28, 9, 11,
        //                                 75, 76, 41, 42, 44, 25, 40, 59, 23, 53, 62, 26, 36, 27, 18, 
        //                                 14, 22, 54, 63, 73, 65, 67, 55, 39, 47, 33, 15, 5, 12, 13, 
        //                                 19, 10, 2, 57, 58, 50, 51, 43, 34, 48, 49, 31, 52, 61, 35, 
        //                                 37, 46, 29, 32, 20, 45}; // knn 67
        // vector<unsigned_t> vPhyQubit =  {19, 18, 25, 24, 26, 27, 21, 28, 22, 15, 23, 20, 16, 29, 
        //                                     10, 17, 11, 9, 14, 8, 4, 13, 7, 12, 6, 1, 2, 3, 0, 34, 
        //                                     31, 35, 33, 5, 32, 30};//qft
        // vector<unsigned_t> vPhyQubit = {16, 29, 12, 4, 9, 14, 3, 13, 7, 17, 8, 2, 35, 21, 10, 22, 1, 20, 24, 6,
        //                                 23, 26, 19, 25, 32, 27, 31, 34, 33, 15, 5, 28, 0, 11, 30, 18}; // bv
        // vector<unsigned_t> vPhyQubit(_device.nQubit());
        // for (unsigned_t i = 0; i < _device.nQubit(); ++i){
        //     vPhyQubit[i] = i;
        // }
        // for(unsigned_t i = 0; i < _circuit.nProgramQubit(); ++i){
        //     _circuit.setInitialMapping(i, vPhyQubit[i]);
        // }
        // aRouterThreePass(arouter, _device, _circuit, true);
        // return;
        placer.disableAllCommute();
    } 
    for(unsigned_t i = 0; i < _circuit.nProgramQubit(); ++i){
        reverseCir.setInitialMapping(i, _circuit.initialMapping(i));
    }
    Circuit cir(_circuit);
    placer.enableOneHopCost();
    placer.run(_circuit, _device, hasRunInitialMapper);
    _bestSAcost = placer.getOptimalCost();
    aRouterThreePass(arouter, _device, _circuit, true);
    // _circuit.printCircuitLayout();
    // return;
    if(_timer.isTimeout()){
        _molsqParam.trials = 3;
        return;
    }
    placer.disableOneHopCost();
    placer.run(cir, _device, hasRunInitialMapper);
    aRouterThreePass(arouter, _device, cir, true);
    if(_timer.isTimeout()){
        _molsqParam.trials = 3;
        return;
    }
    if(cir.nSwapGate() < _circuit.nSwapGate()){
        updateCircuit(_circuit, cir);
        _molsqParam.placer_one_hop_cost = true;
    }
    else{
        placer.disableOneHopCost();
        _molsqParam.placer_one_hop_cost = false;
    }
    _timer.setTimeout(3600);
    if(_timer.isTimeout()){
        _molsqParam.trials = 3;
        return;
    }
    placer.run(reverseCir, _device, hasRunInitialMapper);
    aRouterThreePass(arouter, _device, reverseCir, true);
    // return;
    initialMapper.setTimeout(100);
    if(reverseCir.nSwapGate() < _circuit.nSwapGate()){
        Circuit rCir(reverseCir);
        updateCircuitWithReverse(_circuit, reverseCir);
        for(unsigned_t c = 1; c < _molsqParam.trials; ++c){
            if(_timer.isTimeout()){
                _molsqParam.trials = c;
            }
            if(_circuit.nProgramQubit() < 100){
                initialMapper.run(rCir, _device);
                placer.run(rCir, _device, true);
            }
            else{
                placer.run(rCir, _device, false);
            }
            // arouter.setGCostLimit(_circuit.nSwapGate());
            // bool reduce = arouter.run(cir, _device);
            aRouterThreePass(arouter, _device, rCir, true);
            if(rCir.nSwapGate() < _circuit.nSwapGate() || (rCir.nSwapGate() == _circuit.nSwapGate() && rCir.circuitDepth() < _circuit.circuitDepth())){
                updateCircuitWithReverse(_circuit, rCir);
            }
        }
        _molsqParam.use_reverse_circuit = true;
        for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
            reverseCir.setInitialMapping(i, _circuit.finalMapping(i));
            reverseCir.setFinalMapping(i, _circuit.initialMapping(i));
        }
    }
    else{
        for(unsigned_t c = 1; c < _molsqParam.trials; ++c){
            if(_timer.isTimeout()){
                _molsqParam.trials = c;
            }
            if(_circuit.nProgramQubit() < 100){
                initialMapper.run(cir, _device);
                placer.run(cir, _device, true);
            }
            else{
                placer.run(cir, _device, false);
            }
            // arouter.setGCostLimit(_circuit.nSwapGate());
            // bool reduce = arouter.run(cir, _device);
            aRouterThreePass(arouter, _device, cir, true);
            if(cir.nSwapGate() < _circuit.nSwapGate() || (cir.nSwapGate() == _circuit.nSwapGate() && cir.circuitDepth() < _circuit.circuitDepth())){
                updateCircuit(_circuit, cir);
            }
        }
    }
}

void mOLSQ::updateCircuitWithReverse(Circuit & cir, Circuit & reverseCir){
    fprintf(stdout, "[INFO] Update circuit by reverseCir\n");
    unsigned_t reverseCircuitDepth = reverseCir.circuitDepth();
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        cir.setInitialMapping(i, reverseCir.finalMapping(i));
        cir.setFinalMapping(i, reverseCir.initialMapping(i));
    }
    int_t ir = cir.nGate() - 1;
    for(unsigned_t i = 0; i < cir.nGate(); ++i){
        Gate & gate = cir.gate(i);
        Gate & gateReverse = reverseCir.gate(ir);
        for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
            gate.setTargetPhysicalQubit(j, gateReverse.targetPhysicalQubit(j));
        }
        gate.setExecutionTime(reverseCircuitDepth - gateReverse.executionTime());
        --ir;
    }
    cir.clearSwap();
    ir = 0;
    for(int_t i = reverseCir.nSwapGate()-1; i >= 0; --i){
        Gate & gate = reverseCir.swapGate(i);
        cir.addSwapGate(ir, gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1), reverseCircuitDepth - gate.executionTime());
        ++ir;
    }
    cir.setCircuitDepth(reverseCir.circuitDepth() + 1);
    asapScheduling(cir);
}

void mOLSQ::updateCircuit(Circuit & cir1, Circuit & cir2){
    for(unsigned_t i = 0; i < cir1.nProgramQubit(); ++i){
        cir1.setInitialMapping(i, cir2.initialMapping(i));
        cir1.setFinalMapping(i, cir2.finalMapping(i));
    }
    for(unsigned_t i = 0; i < cir1.nGate(); ++i){
        Gate & gate = cir1.gate(i);
        Gate & gateForward = cir2.gate(i);
        for(unsigned_t j = 0; j < gate.nTargetQubit(); ++j){
            gate.setTargetPhysicalQubit(j, gateForward.targetPhysicalQubit(j));
        }
        gate.setExecutionTime(gateForward.executionTime());
    }
    cir1.clearSwap();
    for(unsigned_t i = 0; i < cir2.nSwapGate(); ++i){
        Gate & gate = cir2.swapGate(i);
        cir1.addSwapGate(i, gate.targetPhysicalQubit(0), gate.targetPhysicalQubit(1), gate.executionTime());
    }
    cir1.setCircuitDepth(cir2.circuitDepth());
}

void mOLSQ::constructReverseCircuit(Circuit & cir, Circuit& reverseCir){
    for(int_t i = cir.nGate()-1; i >=0 ; --i){
        // cerr << i << endl;
        Gate & gate = cir.gate(i);
        reverseCir.addGate(gate.name(), gate.vTargetProgramQubit(), 1);
    }
    if(!_molsqParam.is_all_commute){
        reverseCir.constructDependency();
    }
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        reverseCir.setInitialMapping(i, cir.finalMapping(i));
    }
}

void mOLSQ::aRouterThreePass(aRouter_dev& arouter, Device& device, Circuit& cir, bool disableRestrictRegion){
    arouter.run(cir, device);
    fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\r", cir.nSwapGate(), cir.circuitDepth());
    Circuit reverseCir("reverseCir", cir.nProgramQubit(), cir.nGate());
    constructReverseCircuit(cir, reverseCir);
    arouter.run(reverseCir, device);
    fprintf(stdout, "[INFO] get a solution for reverse circuit with %d SWAP and depth %d\r", reverseCir.nSwapGate(), reverseCir.circuitDepth());
    bool bestIsReverse = reverseCir.nSwapGate() < cir.nSwapGate();
    Circuit bestCir;
    if(bestIsReverse){
        bestCir = reverseCir;
    }
    else{
        bestCir = cir;
    }
    // while(!_timer.isTimeout()){
    while(1){
        for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
            cir.setInitialMapping(i, reverseCir.finalMapping(i));
        }
        arouter.run(cir, device);
        fprintf(stdout, "[INFO] get a solution with %d SWAP and depth %d\r", cir.nSwapGate(), cir.circuitDepth());
        for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
            reverseCir.setInitialMapping(i, cir.finalMapping(i));
        }
        arouter.run(reverseCir, device);
        fprintf(stdout, "[INFO] get a solution for reverse circuit with %d SWAP and depth %d\r", reverseCir.nSwapGate(), reverseCir.circuitDepth());
        if(cir.nSwapGate() < reverseCir.nSwapGate() && (cir.nSwapGate() < bestCir.nSwapGate()  || (cir.nSwapGate() == bestCir.nSwapGate() && cir.circuitDepth() < bestCir.circuitDepth()))){
            bestCir = cir;
            bestIsReverse = false;
        }
        else if(reverseCir.nSwapGate() < bestCir.nSwapGate() || (reverseCir.nSwapGate() == bestCir.nSwapGate() && reverseCir.circuitDepth() < bestCir.circuitDepth())){
                bestCir = reverseCir;
                bestIsReverse = true;
        }
        else{
            break;
        }
    }
    if(bestIsReverse){
        fprintf(stdout, "[INFO] Update circuit by reverseCir with %d SWAP and depth %d\r", bestCir.nSwapGate(), bestCir.circuitDepth());
        updateCircuitWithReverse(cir, bestCir);
    }
    else{
        fprintf(stdout, "[INFO] Update circuit by forwardCir with %d SWAP and depth %d\r", bestCir.nSwapGate(), bestCir.circuitDepth());
        cir = bestCir;
    }
}

void mOLSQ::bfsSearch(Device& device, unordered_set<int_t>& sQubitRegion, unsigned_t expansionTime){
    // collect qubit region from mapping
    map<int_t, int_t> mQubitPathLength;
    unordered_set<int_t>::iterator it = sQubitRegion.begin();
    ++it;
    for (; it != sQubitRegion.end(); ++it) {
        mQubitPathLength[*it] = -1;

    }
    // cerr << "in bfs: ";
    // cerr << "init qubit region: ";
    // for (int_t s : sQubitRegion){
    //     cerr << s << " ";
    // }
    // cerr << endl;
    unordered_set<int_t> sVisitedQubitRegion;
    unsigned_t i, j, q, nIdx;
    vector<Node> vNode;
    priority_queue<int_t, std::vector<int>, std::greater<int> > priorityQ;
    vector<int_t> vBacktraceNode;
    unsigned_t cost, nSpanEdge, maxDis = 0;
    vNode.emplace_back(-1, 0, (*sQubitRegion.begin()), -1, 0);
    priorityQ.push(0);
    while(priorityQ.size() > 0 && (sVisitedQubitRegion.size() != sQubitRegion.size() || vNode[priorityQ.top()].dis <= maxDis)){
        nIdx = priorityQ.top();
        priorityQ.pop();
        // cerr << "expand node " << vNode.at(nIdx).idx << " phy q: " << vNode.at(nIdx).qIdx << " parent idx: " << vNode.at(nIdx).parentIdx << ", dis: " << vNode.at(nIdx).dis << endl;
        // if(nIdx >= vNode.size()){
        //     cerr << "line 741: out of the bound access with nIdx " << nIdx << endl;
        // }
        // if(vNode.at(nIdx).qIdx >= device.nQubit()){
        //     cerr << "line 744: out of the bound access with qubit idx " << vNode.at(nIdx).qIdx << endl;
        // }
        Qubit& qubit = device.qubit(vNode.at(nIdx).qIdx);
        if(sQubitRegion.find(vNode.at(nIdx).qIdx) != sQubitRegion.end() && (mQubitPathLength[vNode.at(nIdx).qIdx] == -1 || mQubitPathLength[vNode.at(nIdx).qIdx] >= vNode.at(nIdx).dis)){
            // cerr << "find node " << vNode.at(nIdx).idx << " qubit " << vNode.at(nIdx).qIdx << " with distance " << vNode.at(nIdx).dis << " parent idx " << vNode.at(nIdx).parentIdx << endl;
            maxDis = (vNode.at(nIdx).dis > maxDis) ? vNode.at(nIdx).dis : maxDis;
            sVisitedQubitRegion.insert(vNode.at(nIdx).qIdx);
            vBacktraceNode.emplace_back(vNode.at(nIdx).idx);
        }
        // sTraversedQubit.insert(node.qIdx);
        nSpanEdge = qubit.vSpanEdge.size();
        for (j = 0; j < nSpanEdge; ++j){
            q = device.edge(qubit.vSpanEdge.at(j)).qubitId1();
            if (qubit.idx == q){
                q = device.edge(qubit.vSpanEdge.at(j)).qubitId2();
            }
                // if (sTraversedQubit.find(q) == sTraversedQubit.end()){
                // sTraversedQubit.insert(q);
            // cerr << "node idx " << vNode.at(nIdx).idx << endl;
            vNode.emplace_back(vNode.at(nIdx).idx, vNode.size(), q, qubit.vSpanEdge.at(j), vNode.at(nIdx).dis+1);
            // cerr << "add node " << vNode.back().idx << " qubit " << vNode.back().qIdx << " with distance " << vNode.back().dis << " parent idx " << vNode.back().parentIdx << endl;
            // queue.emplace_back(vNode.back().idx);
            priorityQ.push(vNode.back().idx);
        }
        // cerr << "top node dis: " << vNode[priorityQ.top()].dis << endl;
        // getchar();
    }
    // backtrack 
    // cerr << "start backtrace (" << vBacktraceNode.size() << ")" << endl;
    int_t pIdx;
    for(int_t nIdx : vBacktraceNode){
        pIdx = nIdx;
        // cerr << "find node " << vNode.at(nIdx).idx << " qubit " << vNode.at(nIdx).qIdx << " with distance " << vNode.at(nIdx).dis << " parent idx " << vNode.at(nIdx).parentIdx << endl;
        while(pIdx > 0){
            // cerr << "pIdx: " << pIdx << endl;
            // if(pIdx >= vNode.size()){
            //     cerr << "line 780: out of the bound access with pIdx " << pIdx << endl;
            // }
            sQubitRegion.insert(vNode.at(pIdx).qIdx);
            // cerr << "add edge " << vNode[pIdx].parentEIdx << endl;
            pIdx = vNode.at(pIdx).parentIdx;
        }
    }
    // cerr << "qubit region before expansion: ";
    // for (int_t s : sQubitRegion){
    //     cerr << s << " ";
    // }
    // cerr << endl;
    unordered_set<int_t> sExpandQubit;
    unordered_set<int_t> sExpandQubit2;
    // cerr << "sQubitRegion.size(): " << sQubitRegion.size() << endl;
    sExpandQubit.clear();
    sExpandQubit2.clear();
    for (it = sQubitRegion.begin(); it != sQubitRegion.end(); ++it) {
        q = *(it);
        // if(q >= device.nQubit()){
        //     cerr << "line 802: out of the bound access with qubit q " << q << endl;
        // }
        Qubit& qubit = device.qubit(q);
        nSpanEdge = qubit.vSpanEdge.size();
        for(j = 0; j < nSpanEdge; ++j){
            Edge& e = device.edge(qubit.vSpanEdge[j]);
            sExpandQubit.insert(e.qubitId1());
            sExpandQubit.insert(e.qubitId2());
        }
    }
    // cerr << "merge" << endl;
    sQubitRegion.merge(sExpandQubit);
    // cerr << "enter exansion " << expansionTime << endl;
    for(i = 1; i < expansionTime; ++i){
        if(i % 2 == 0){
            sExpandQubit.clear();
            for (it = sExpandQubit2.begin(); it != sExpandQubit2.end(); ++it) {
                q = *(it);
                // if(q >= device.nQubit()){
                //     cerr << "line 819: out of the bound access with qubit q " << q << endl;
                // }
                Qubit& qubit = device.qubit(q);
                nSpanEdge = qubit.vSpanEdge.size();
                for(j = 0; j < nSpanEdge; ++j){
                    Edge& e = device.edge(qubit.vSpanEdge.at(j));
                    sExpandQubit.insert(e.qubitId1());
                    sExpandQubit.insert(e.qubitId2());
                }
            }
            // cerr << "merge" << endl;
            // sQubitRegion.merge(sExpandQubit);
            // cerr << "finsh merge" << endl;
        }
        else{
            sExpandQubit2.clear();
            for (it = sExpandQubit.begin(); it != sExpandQubit.end(); ++it) {
                q = *(it);
                // if(q >= device.nQubit()){
                //     cerr << "line 836: out of the bound access with qubit q " << q << endl;
                // }
                Qubit& qubit = device.qubit(q);
                nSpanEdge = qubit.vSpanEdge.size();
                for(j = 0; j < nSpanEdge; ++j){
                    Edge& e = device.edge(qubit.vSpanEdge.at(j));
                    sExpandQubit2.insert(e.qubitId1());
                    sExpandQubit2.insert(e.qubitId2());
                }
            }
            // cerr << "merge" << endl;
            // sQubitRegion.merge(sExpandQubit2);
            // cerr << "finsh merge" << endl;
        }

    }
    // if(sQubitRegion.size() == 1){
    //     q = *(sQubitRegion.begin());
    //     Qubit& qubit = _pDevice->qubit(q);
    //     nSpanEdge = qubit.vSpanEdge.size();
    //     for(i = 0; i < nSpanEdge; ++i){
    //         sQubitRegion.insert(i);
    //         sSwapRegion.insert(qubit.vSpanEdge[i]);
    //     }
    // }
    // cerr << "final qubit region: ";
    // for (int_t s : sQubitRegion){
    //     cerr << s << " ";
    // }
    // cerr << endl;
}

void mOLSQ::insertSWAP(vector<vector<pair<unsigned_t, unsigned_t> > >& vvSwap){
    unsigned_t swapId = 0, level = 0;
    _circuit.clearSwap();
    for (vector<pair<unsigned_t, unsigned_t> >& vSwap : vvSwap){
        for (pair<unsigned_t, unsigned_t> & qubitPair : vSwap){
            _circuit.addSwapGate(swapId, qubitPair, 1);
            Gate & gate = _circuit.swapGate(swapId);
            gate.setExecutionTime(level);
            ++swapId;
        }   
        ++level;
    }                
}


void mOLSQ::asapScheduling(Circuit & cir){
    fprintf(stdout, "[INFO] ASAP Scheduling: Start                                               \n");
    vector<int_t> vPushForwardDepth(_device.nQubit(), -1);
    int_t gateExecutionTime;
    unsigned_t block, i, j, q0, q1, qId, maxTime = 0;
    Gate gate;
    unordered_set<unsigned_t> sGateId;
    for (block = 0; block < cir.circuitDepth(); ++block){
        // cerr << "add gate for block " << block << endl; 
        for (i = 0; i < cir.nGate(); ++i){
            Gate & gate = cir.gate(i);
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
                // cerr << "add gate " << i << " at time " << gateExecutionTime << endl; 
            }
        }
        if (block < cir.circuitDepth() - 1){
            for (j = 0; j < cir.nSwapGate(); ++j){
                Gate & gate = cir.swapGate(j);
                if (gate.executionTime() == block && sGateId.count(gate.idx()+ cir.nGate()) == 0){
                    q0 = gate.targetPhysicalQubit(0);
                    q1 = gate.targetPhysicalQubit(1);
                    gateExecutionTime = (vPushForwardDepth[q0] < vPushForwardDepth[q1]) ? vPushForwardDepth[q1] : vPushForwardDepth[q0];
                    ++gateExecutionTime;
                    vPushForwardDepth[q0] = gateExecutionTime;
                    vPushForwardDepth[q1] = gateExecutionTime;
                    maxTime = (maxTime < gateExecutionTime) ? gateExecutionTime : maxTime;
                    gate.setExecutionTime(gateExecutionTime);
                    sGateId.insert(gate.idx()+ cir.nGate());
                    // cerr << "add swap gate " << gate.idx() << " at time " << gateExecutionTime << endl; 
                }
            }
        }
    }
    ++maxTime;
    cir.setCircuitDepth(maxTime);
    fprintf(stdout, "[INFO] ASAP Scheduling: Finish                                               \n");
}


void mOLSQ::glueMapping(unsigned_t level){
    vector<pair<unsigned_t, unsigned_t> > vSwap;
    rOLSQ2 rolsq2(_device);
    // finest level
    vector<unsigned_t> vDepth;
    vector<vector<pair<unsigned_t, unsigned_t> > > vvSwap;
    for(unsigned_t i = 0; i < level; ++i){
        // iteratively glue two mapping
        rolsq2.setInitialMapping(_vvQubitMapping[i]);
        rolsq2.setFinalMapping(_vvQubitMapping[i+1]);
        // calculate the largest moving distance to serve as the depth bonud
        vDepth.emplace_back(estimateSwapDepth(_vvQubitMapping[i], _vvQubitMapping[i+1]));
        // cerr << "max swap count: " << vDepth[i] << endl;
        rolsq2.setMinDepth(2);
        rolsq2.run(vDepth[i]);
        vvSwap.emplace_back(rolsq2.vSwap());
    }
    insertSWAP(vvSwap);
}

unsigned_t mOLSQ::estimateSwapDepth(vector<unsigned_t>& vInitialMapping, vector<unsigned_t>& vFinalMapping){
    // unsigned_t maxDis = 0;
    // for(unsigned_t i = 0; i < vInitialMapping.size(); ++i){
    //     if(vInitialMapping[i] != vFinalMapping[i]){
    //         maxDis = max(maxDis, _device.getDistance(vInitialMapping[i], vFinalMapping[i]));
    //     }
    // }
    // return maxDis;
    unsigned_t totalDis = 0;
    unordered_set<pair<int_t, int_t>> sChangePair;
    pair<int_t, int_t> changePair;
    for(unsigned_t i = 0; i < vInitialMapping.size(); ++i){
        if(vInitialMapping[i] != vFinalMapping[i]){
            if(vInitialMapping[i] < vFinalMapping[i]){
                changePair.first = vInitialMapping[i];
                changePair.second = vFinalMapping[i];
            }
            else{
                changePair.first = vFinalMapping[i];
                changePair.second = vInitialMapping[i];
            }
            if (sChangePair.find(changePair) == sChangePair.end()){
                ++totalDis;
            }
        }
    }
    return totalDis;
}


void mOLSQ::verify(Device & device, Circuit & cir){
    fprintf(stdout, "===============================================\n");
    fprintf(stdout, "====            Verify Result              ====\n");
    fprintf(stdout, "===============================================\n");
    bool success = true;
    vector<unsigned_t> vCurMapping(cir.nProgramQubit(), 0);
    vector<bool> vPhyIsOccupied(device.nQubit(), 0);
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        if(vPhyIsOccupied[cir.initialMapping(i)]){
            fprintf(stdout, "[Error] injectivity violation: physical qubit %d is occuplied by more than one qubit\n", cir.initialMapping(i));    
            success = false;
        }
        vPhyIsOccupied[cir.initialMapping(i)] = 1;
        vCurMapping[i] = cir.initialMapping(i);
    }
    vector<vector<unsigned_t>> vvGateBlock(cir.circuitDepth(), vector<unsigned_t>());
    for(unsigned_t i = 0; i < cir.nGate(); ++i){
        Gate & gate = cir.gate(i);
        vvGateBlock[gate.executionTime()].emplace_back(i);
    }
    for(unsigned_t i = 0; i < cir.nSwapGate(); ++i){
        Gate & gate = cir.swapGate(i);
        vvGateBlock[gate.executionTime()].emplace_back(gate.idx() + cir.nGate());
    }
    unsigned_t p0, p1, q0, q1;
    unordered_set<unsigned_t> isBusyQubit;
    for(unsigned_t i = 0; i < vvGateBlock.size(); ++i){
        isBusyQubit.clear();
        for(unsigned_t j = 0; j < vvGateBlock[i].size(); ++j){
            if((vvGateBlock[i][j] < cir.nGate())){
                Gate & gate = cir.gate(vvGateBlock[i][j]);
                q0 = gate.targetProgramQubit(0);
                if(!(vCurMapping[q0] == gate.targetPhysicalQubit(0))){
                    fprintf(stdout, "[Error] inconsistency between gate and mapping: gate %d's target qubit %d is schedule on physical qubit %d its actual mapping is %d\n", vvGateBlock[i][j], q0,  gate.targetPhysicalQubit(0), vCurMapping[q0]);    
                    success = false;
                }
                else if(isBusyQubit.count(vCurMapping[q0]) > 0){
                    fprintf(stdout, "[Error] two gates act on the same qubit %d at time %d\n", vCurMapping[q0], i);    
                    success = false;
                }
                isBusyQubit.insert(vCurMapping[q0]);
                if(gate.nTargetQubit() > 1){
                    q1 = gate.targetProgramQubit(1);
                    if(!(vCurMapping[q1] == gate.targetPhysicalQubit(1))){
                        fprintf(stdout, "[Error] inconsistency between gate and mapping: gate %d's target qubit %d is schedule on physical qubit %d its actual mapping is %d\n", vvGateBlock[i][j], q1,  gate.targetPhysicalQubit(1), vCurMapping[q1]);    
                        success = false;
                    }
                    if(!device.isAdjacent(vCurMapping[q0], vCurMapping[q1])){
                        fprintf(stdout, "[Error] gate %d (%d,%d) acts on nonadjacent physical qubits %d(%d) and %d(%d)\n", vvGateBlock[i][j], gate.targetProgramQubit(0), gate.targetProgramQubit(1), gate.targetPhysicalQubit(0), vCurMapping[q0], gate.targetPhysicalQubit(1), vCurMapping[q1]);    
                        success = false;
                    }
                    if(isBusyQubit.count(vCurMapping[q1]) > 0){
                        fprintf(stdout, "[Error] two gates act on the same qubit %d at time %d\n", vCurMapping[q1], i);    
                        success = false;
                    }
                    isBusyQubit.insert(vCurMapping[q1]);
                }
            }
            else{
                Gate & gate = cir.swapGate(vvGateBlock[i][j] - cir.nGate());
                p0 = gate.targetPhysicalQubit(0);
                p1 = gate.targetPhysicalQubit(1);
                if(!device.isAdjacent(p0, p1)){
                    fprintf(stdout, "[Error] swap gate %d acts on nonadjacent physical qubits %d and %d\n", vvGateBlock[i][j] - cir.nGate(), p0, p1);    
                    success = false;
                }
                if(isBusyQubit.count(p0) > 0){
                    fprintf(stdout, "[Error] two gates act on the same qubit %d at time %d\n", p0, i);    
                    success = false;
                }
                if(isBusyQubit.count(p1) > 0){
                    fprintf(stdout, "[Error] two gates act on the same qubit %d at time %d\n", p1, i);    
                    success = false;
                }
                isBusyQubit.insert(p0);
                isBusyQubit.insert(p1);
                auto res1 = find(vCurMapping.begin(), vCurMapping.end(), p0);
                auto res2 = find(vCurMapping.begin(), vCurMapping.end(), p1);
                // cerr << "eid: " << eId << endl;
                // cerr << "swap on " << p0 << " " << p1 << " at time " << i << endl;
                if(res1 != vCurMapping.end() && res2 != vCurMapping.end()){
                    // cerr << "both physical qubits are mapped by program qubits" << endl; 
                    q0 = res1 - vCurMapping.begin();
                    q1 = res2 - vCurMapping.begin();
                    vCurMapping[q0] = p1;
                    vCurMapping[q1] = p0;
                }
                else if(res1 != vCurMapping.end()){
                    q0 = res1 - vCurMapping.begin();
                    // cerr << "only " << p0 << " is mapped by program qubit " << q0 << endl; 
                    vCurMapping[q0] = p1;

                }
                else if(res2 != vCurMapping.end()){
                    q1 = res2 - vCurMapping.begin();
                    // cerr << "only " << p1 << " is mapped by program qubits " << q1 << endl; 
                    vCurMapping[q1] = p0;
                }
                else{
                    fprintf(stdout, "[Error] swap on two unused qubit %d and %d at time %d\n", p0, p1, i);    
                    success = false;
                }   
            }
        }    
    }
    if(!_molsqParam.is_all_commute){
        vector<pair<unsigned_t, unsigned_t> >* pvpGateDependency = cir.pvpGateDependency();
        for(auto& pGateDependency : (*pvpGateDependency)){
            Gate& gate0 = cir.gate(pGateDependency.first);
            Gate& gate1 = cir.gate(pGateDependency.second);
            if(gate0.executionTime() >= gate1.executionTime()){
                fprintf(stdout, "[Error] dependency violation for gate %d and gate %d \n", pGateDependency.first,  pGateDependency.second);        
                success = false;
            }
        }

    }
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        if(! vCurMapping[i] == cir.finalMapping(i)){
            fprintf(stdout, "[Error] inconsistency for final mapping: program qubit %d in vCurmapping is %d but in cir.finalMapping is %d\n", i,  vCurMapping[i], cir.finalMapping(i));    
            success = false;
        }
    }
    if(success){
        fprintf(stdout, "PASS!!!\n");
    }
    else{
        fprintf(stdout, "FAIL...\n");
        _circuit.printCircuitLayout();
    }
}


MOLSQ_NAMESPACE_CPP_END

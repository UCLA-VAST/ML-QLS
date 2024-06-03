/***********************************************************************
  File        [ clusterer.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2022 ]
***********************************************************************/

#include "clusterer/clusterer2.hpp"

MOLSQ_NAMESPACE_CPP_START

void Clusterer2::cluster(Circuit & cir, Device & device, vector<Circuit>& vCir, vector<Device>& vDevice, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserProQ2FinerProQ, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserPhyQ2FinerPhyQ, bool allCommute){
    // cluster the first device
    _isAllCommute = allCommute;
    clusterFinestDevice(cir, device, vDevice[0], vmvCoarserPhyQ2FinerPhyQ[0]);
    
    // cluster circuit
    clusterCircuit(cir, vCir[0], vmvCoarserPhyQ2FinerPhyQ[0], vmvCoarserProQ2FinerProQ[0], device.nQubit());

    for(unsigned_t i = 1; i < vCir.size(); ++i){
        // cluster device
        clusterPhysicalQubit(vDevice[i-1], vDevice[i], vmvCoarserPhyQ2FinerPhyQ[i]);
        // cluster circuit
        clusterCircuit(vCir[i-1], vCir[i], vmvCoarserPhyQ2FinerPhyQ[i], vmvCoarserProQ2FinerProQ[i], vDevice[i-1].nQubit());
    }
    
}

void Clusterer2::clusterFinestDevice(Circuit & cir, Device & device, Device & newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    unordered_set<unsigned_t> sQubitRegion;
    collectFinestQubitRegion(cir, device, sQubitRegion);
    my_graph graph(sQubitRegion.size());
    consructBoostGraph(graph, device, sQubitRegion);
    vector<V> mate(num_vertices(graph));
    vector<unsigned_t> vFinerQ2CoarserQ(device.nQubit(), device.nQubit());
    // use boost library to construct maximum card matching 
    // see: https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/maximum_matching.html
    edmonds_maximum_cardinality_matching(graph, &mate[0]);
    if(_verbose > 1 ){
        printMatching(graph, mate);
    }
    // construct new device
    extractCoarserQubit(device, graph, mate, mvCoarserQ2FinerQ, vFinerQ2CoarserQ);
    constructNewDevice(device, newDevice, mvCoarserQ2FinerQ, vFinerQ2CoarserQ);
    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserQ2FinerQ);
    }
}

void Clusterer2::collectFinestQubitRegion(Circuit & cir, Device & device, unordered_set<unsigned_t> & sQubitRegion){
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        sQubitRegion.insert(cir.initialMapping(i));
    }
    unsigned_t expansionTime = 2;
    unordered_set<unsigned_t> sExpandQubit;
    unordered_set<unsigned_t> sExpandQubit2;
    // cerr << "sQubitRegion.size(): " << sQubitRegion.size() << endl;
    sExpandQubit.clear();
    sExpandQubit2.clear();
    unordered_set<unsigned_t>::iterator it;
    for (it = sQubitRegion.begin(); it != sQubitRegion.end(); ++it) {
        Qubit& qubit = device.qubit(*(it));
        for(unsigned_t j = 0; j < qubit.vSpanEdge.size(); ++j){
            Edge& e = device.edge(qubit.vSpanEdge[j]);
            sExpandQubit.insert(e.qubitId1());
            sExpandQubit.insert(e.qubitId2());
        }
    }
    sQubitRegion.merge(sExpandQubit);
    
    for(unsigned_t i = 1; i < expansionTime; ++i){
        if(i % 2 == 0){
            sExpandQubit.clear();
            for (it = sExpandQubit2.begin(); it != sExpandQubit2.end(); ++it) {
                Qubit& qubit = device.qubit(*(it));
                for(unsigned_t j = 0; j < qubit.vSpanEdge.size(); ++j){
                    Edge& e = device.edge(qubit.vSpanEdge[j]);
                    sExpandQubit.insert(e.qubitId1());
                    sExpandQubit.insert(e.qubitId2());
                }
            }
            sQubitRegion.merge(sExpandQubit);
        }
        else{
            sExpandQubit2.clear();
            for (it = sExpandQubit.begin(); it != sExpandQubit.end(); ++it) {
                Qubit& qubit = device.qubit(*(it));
                for(unsigned_t j = 0; j < qubit.vSpanEdge.size(); ++j){
                    Edge& e = device.edge(qubit.vSpanEdge[j]);
                    sExpandQubit2.insert(e.qubitId1());
                    sExpandQubit2.insert(e.qubitId2());
                }
            }
            sQubitRegion.merge(sExpandQubit);
        }
    }
}

void Clusterer2::consructBoostGraph(my_graph& graph, Device& device, unordered_set<unsigned_t>& sQubitRegion){
    unordered_set<unsigned_t>::iterator iter, iter1;
    for(iter = sQubitRegion.begin(); iter != sQubitRegion.end(); ++iter){
        iter1 = iter;
        ++iter1;
        for(; iter1 != sQubitRegion.end(); ++iter1){
            if(device.isAdjacent((*iter), (*iter1))){
                add_edge((*iter), (*iter1), EdgeProperty(1), graph);    
            }
        }
    }
}

void Clusterer2::clusterCircuit(Circuit & cir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, unsigned_t nPhyQ){
    vector<int> vPhy2Pro(nPhyQ, -1);
    // cerr << "print phy->pro" << endl;
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        vPhy2Pro[cir.initialMapping(i)] = i;
        // cerr << cir.initialMapping(i) << "->" << i << endl;
    }
    int_t qId = -1;
    bool hasInitialized;
    vector<unsigned_t> vFinerQ2CoarserQ(cir.nProgramQubit(), 0);
    vector<unsigned_t> vProId2PhyId;
    for(auto& pair_key_value : mvCoarserPhyQ2FinerPhyQ){
        hasInitialized = 0;
        for(unsigned_t i : pair_key_value.second){
            if(vPhy2Pro[i] >= 0){
                // cerr << "qubit " << vPhy2Pro[i] << "(" << i << ")" << " is in coarser qubit " << qId << endl;
                if(!hasInitialized){
                    ++qId;
                    mvCoarserProQ2FinerProQ[qId] = vector<unsigned_t>();
                    hasInitialized = 1;
                    vProId2PhyId.emplace_back(pair_key_value.first);
                }
                mvCoarserProQ2FinerProQ[qId].emplace_back(vPhy2Pro[i]);
                vFinerQ2CoarserQ[vPhy2Pro[i]] = qId;
            }
        }
    }
    constructNewCircuit(cir, newCir, vFinerQ2CoarserQ, mvCoarserProQ2FinerProQ);
    for(unsigned_t i = 0; i < newCir.nProgramQubit(); ++i){
        newCir.setInitialMapping(i, vProId2PhyId[i]);
    }

    if(_verbose > 2 ){
        fprintf(stdout, "[Info] Original Circuit:                        \n");
        cir.printCircuit();
    }
    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserProQ2FinerProQ);
    }
    if(_verbose > 0 ){
        fprintf(stdout, "[Info] Coarse Circuit:                        \n");
        newCir.printCircuit();
    }
}

// // ! for abalation study
void Clusterer2::clusterCircuitNaive(Circuit & cir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, unsigned_t nPhyQ){
    vector<int> vPhy2Pro(nPhyQ, -1);
    // cerr << "print phy->pro" << endl;
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        vPhy2Pro[cir.initialMapping(i)] = i;
        // cerr << cir.initialMapping(i) << "->" << i << endl;
    }
    vector<pair<double_t, pair<unsigned_t, unsigned_t> > > vpWeightedEdge;
    int_t i = 0, j, q0, q1;
    vector<unsigned_t> vQubitLevel;
    unordered_map<pair<unsigned_t , unsigned_t>, double_t, pair_hash> umEdge2Weight;
    pair<unsigned_t , unsigned_t> edge;
    for(unsigned_t i = 0; i < cir.nGate(); ++i){
        Gate & gate = cir.gate(i);
        if(gate.nTargetQubit()>1){
            q0 = gate.targetProgramQubit(0);
            q1 = gate.targetProgramQubit(1);
            edge = (q0 < q1) ? make_pair(q0,q1) : make_pair(q1,q0); 
            if (umEdge2Weight.find(edge) == umEdge2Weight.end()){
                umEdge2Weight[edge] = 1;
            }
            else{
                ++umEdge2Weight[edge];
            }
            // cout << "umEdge2Weight: " << umEdge2Weight[edge] << endl;
        }
    }
    vpWeightedEdge.clear(); 
    vpWeightedEdge.reserve(umEdge2Weight.size());
    for (auto & keyItemPair : umEdge2Weight){
        vpWeightedEdge.emplace_back(keyItemPair.second, keyItemPair.first);
    }
    sort(vpWeightedEdge.rbegin(), vpWeightedEdge.rend());

    int_t qId = 0;
    vector<unsigned_t> vFinerQ2CoarserQ(cir.nProgramQubit(), 0);
    vector<unsigned_t> vProId2PhyId;
    vector<bool> vIsQubitClustered(cir.nProgramQubit(), 0);
    unsigned_t clusterNum = 0;
    for(unsigned_t i = 0; i < vpWeightedEdge.size(); ++i){
        if(!vIsQubitClustered[vpWeightedEdge[i].second.first] && !vIsQubitClustered[vpWeightedEdge[i].second.second] && mvCoarserProQ2FinerProQ.size() < mvCoarserPhyQ2FinerPhyQ.size()){
            mvCoarserProQ2FinerProQ[qId] = vector<unsigned_t>();
            mvCoarserProQ2FinerProQ[qId].emplace_back(vpWeightedEdge[i].second.first);
            mvCoarserProQ2FinerProQ[qId].emplace_back(vpWeightedEdge[i].second.second);
            vFinerQ2CoarserQ[vpWeightedEdge[i].second.first] = qId;
            vFinerQ2CoarserQ[vpWeightedEdge[i].second.second] = qId;
            vIsQubitClustered[vpWeightedEdge[i].second.first] = true;
            vIsQubitClustered[vpWeightedEdge[i].second.second] = true;
            clusterNum += 2;
            ++qId;
        }
    }
    unsigned_t bound = 4;
    while(clusterNum < cir.nProgramQubit()){
        for(unsigned_t i = 0; i < vpWeightedEdge.size() && clusterNum < cir.nProgramQubit(); ++i){
            if(!vIsQubitClustered[vpWeightedEdge[i].second.first] && vIsQubitClustered[vpWeightedEdge[i].second.second]){
                if(mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[vpWeightedEdge[i].second.second]].size() < bound){
                    mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[vpWeightedEdge[i].second.second]].emplace_back(vpWeightedEdge[i].second.first);
                    vFinerQ2CoarserQ[vpWeightedEdge[i].second.first] = vFinerQ2CoarserQ[vpWeightedEdge[i].second.second];
                    vIsQubitClustered[vpWeightedEdge[i].second.first] = true;
                    ++clusterNum;
                }
            }
            else if(vIsQubitClustered[vpWeightedEdge[i].second.first] && !vIsQubitClustered[vpWeightedEdge[i].second.second]){
                if(mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[vpWeightedEdge[i].second.first]].size() < bound){
                    mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[vpWeightedEdge[i].second.first]].emplace_back(vpWeightedEdge[i].second.second);
                    vFinerQ2CoarserQ[vpWeightedEdge[i].second.second] = vFinerQ2CoarserQ[vpWeightedEdge[i].second.first];
                    vIsQubitClustered[vpWeightedEdge[i].second.second] = true;
                    ++clusterNum;
                }
            }
        }
        for(unsigned_t i = 0; i < cir.nProgramQubit() && clusterNum < cir.nProgramQubit(); ++i){
            if(!vIsQubitClustered[i]){
                for(unsigned_t j = i+1; j < cir.nProgramQubit() && clusterNum < cir.nProgramQubit(); ++i){
                    if(vIsQubitClustered[j] && mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[j]].size() < bound){
                        mvCoarserProQ2FinerProQ[vFinerQ2CoarserQ[j]].emplace_back(i);
                        vFinerQ2CoarserQ[i] = vFinerQ2CoarserQ[j];
                        vIsQubitClustered[i] = true;
                        ++clusterNum;
                    }
                }
            }
        }
        ++bound;
    }
    constructNewCircuit(cir, newCir, vFinerQ2CoarserQ, mvCoarserProQ2FinerProQ);
    // for(unsigned_t i = 0; i < newCir.nProgramQubit(); ++i){
    //     newCir.setInitialMapping(i, vProId2PhyId[i]);
    // }

    if(_verbose > 2 ){
        fprintf(stdout, "[Info] Original Circuit:                        \n");
        cir.printCircuit();
    }
    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserProQ2FinerProQ);
    }
    if(_verbose > 0 ){
        fprintf(stdout, "[Info] Coarse Circuit:                        \n");
        newCir.printCircuit();
    }
}

void Clusterer2::clusterPhysicalQubit(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    // construct boost graph
    my_graph graph(oriDevice.nQubit());
    consructBoostGraph(graph, oriDevice);
    vector<V> mate(num_vertices(graph));
    // use boost library to construct maximum card matching 
    // see: https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/maximum_matching.html
    edmonds_maximum_cardinality_matching(graph, &mate[0]);
    if(_verbose > 1 ){
        printMatching(graph, mate);
    }

    // construct new device
    vector<unsigned_t> vFinerQ2CoarserQ(oriDevice.nQubit(), oriDevice.nQubit());
    extractCoarserQubit(oriDevice, graph, mate, mvCoarserQ2FinerQ, vFinerQ2CoarserQ);

    // ! a special case for 6*6 grid
    // unsigned_t j = 0;
    // for(unsigned_t i = 0; i < oriDevice.nQubit(); i+=2){
    //     vFinalMapping[i] = j;
    //     vFinalMapping[i+1] = j;
    //     ++j
    //     mvCoarserQ2FinerQ[j] = {i , i + 1}
    // }

    constructNewDevice(oriDevice, newDevice, mvCoarserQ2FinerQ, vFinerQ2CoarserQ);
    if(_verbose > 2 ){
        fprintf(stdout, "[Info] Original Device:                        \n");
        oriDevice.printDevice();
    }

    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserQ2FinerQ);
    }
    if(_verbose > 0 ){
        fprintf(stdout, "[Info] Coarse Device:                        \n");
        newDevice.printDevice();
    }
}

void Clusterer2::extractCoarserQubit(Device& oriDevice, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ){
    // construct new qubit by matching result
    unsigned_t countFinerQubit = 0, q1, q2, qId = 0, i, j;
    vector<bool> vIsCluster(oriDevice.nQubit(), 0);
    // collect finer to coaerser mapping
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            mvCoarserQ2FinerQ[qId] = vector<unsigned_t>();
            vFinerQ2CoarserQ[v] = qId;
            vFinerQ2CoarserQ[mate[v]] = qId;
            mvCoarserQ2FinerQ[qId].emplace_back(v);
            mvCoarserQ2FinerQ[qId].emplace_back(mate[v]);
            vIsCluster[v] = 1;
            vIsCluster[mate[v]] = 1;
            countFinerQubit += 2;
           ++qId;
        }
    }
    // collect the uncluster phy qubit to its neighbor
    unsigned_t bound = 3;
    while(countFinerQubit < oriDevice.nQubit()){
        for(i = 0; i < oriDevice.nQubit(); ++i){
            if(!vIsCluster[i]){
                for(j = i + 1; j < oriDevice.nQubit(); ++j){
                    if(!vIsCluster[j]){
                        mvCoarserQ2FinerQ[qId] = vector<unsigned_t>();
                        vFinerQ2CoarserQ[i] = qId;
                        vFinerQ2CoarserQ[j] = qId;
                        mvCoarserQ2FinerQ[qId].emplace_back(i);
                        mvCoarserQ2FinerQ[qId].emplace_back(j);
                        vIsCluster[i] = 1;
                        vIsCluster[j] = 1;
                        countFinerQubit += 2;
                        ++qId;
                        break;
                    }
                }
            }
            if(!vIsCluster[i]){
                Qubit & qubit = oriDevice.qubit(i);
                for( unsigned_t j : qubit.vSpanEdge){
                    Edge & edge = oriDevice.edge(j);
                    q1 = (i == edge.qubitId1()) ? edge.qubitId2() : edge.qubitId1();
                    if(mvCoarserQ2FinerQ[vFinerQ2CoarserQ[q1]].size() < bound){
                        // cerr << "cluster qubit " << vFinerQ2CoarserQ[q1] << " add finer qubit " << i << " (because " << q1 << " )"<< endl;
                        mvCoarserQ2FinerQ[vFinerQ2CoarserQ[q1]].emplace_back(i);
                        vFinerQ2CoarserQ[i] = vFinerQ2CoarserQ[q1];
                        vIsCluster[i] = 1;
                        ++countFinerQubit;
                        break;
                    }
                }
            }
        }
        ++bound;
    }
    // check correctness
    assert((countFinerQubit == oriDevice.nQubit()) || (countFinerQubit == num_vertices(graph)));
    for (unsigned_t i = 0; i < mvCoarserQ2FinerQ.size(); ++i){
        for (unsigned_t j : mvCoarserQ2FinerQ[i]){
            // cerr << "i " << i << ", vFinerQ2CoarserQ[mvCoarserQ2FinerQ[i][j]: " << vFinerQ2CoarserQ[j] << ", mvCoarserQ2FinerQ[i][j]: " << j << endl;
            assert(i == vFinerQ2CoarserQ[j]);
        }
    }
}


void Clusterer2::constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    // construct new circuit by matching result
    newCir.setQubitNum(mvCoarserQ2FinerQ.size());
    vector<unsigned_t> vTargetQubit(2);
    vector<int_t> vQubitLastGate(newCir.nProgramQubit(), -1);
    vector<vector<bool>> vvQubitEntangle;
    if(_isAllCommute){
        vvQubitEntangle.resize(newCir.nProgramQubit(), vector<bool>(newCir.nProgramQubit(), false));
    }
    bool cond1, cond2, cond3 = 1;
    // construct new gates and add gate dependency
    for (unsigned_t i = 0; i < oriCir.nGate();  ++i ){
        Gate & gate = oriCir.gate(i);
        if(gate.nTargetQubit() == 2){
            vTargetQubit[0] = vFinerQ2CoarserQ[gate.targetProgramQubit(0)];
            vTargetQubit[1] = vFinerQ2CoarserQ[gate.targetProgramQubit(1)];
            // if two finer qubits act on two different coarser qubit, add gate
            cond1 = vTargetQubit[0] != vTargetQubit[1];
            // if this is not the consecutive gate acting on the same pair of qubits as the previous gates, e.g., (0, 1), (0, 1)
            cond2 = (vQubitLastGate[vTargetQubit[0]] != vQubitLastGate[vTargetQubit[1]]) || (vQubitLastGate[vTargetQubit[0]] == -1);
            if(_isAllCommute){
                cond3 = !vvQubitEntangle[vTargetQubit[0]][vTargetQubit[1]];
                vvQubitEntangle[vTargetQubit[0]][vTargetQubit[1]] = true;
                vvQubitEntangle[vTargetQubit[1]][vTargetQubit[0]] = true;
            }
            if(cond1 && ((!_isAllCommute && cond2) || (_isAllCommute && cond3))){
                // cerr << "gate.targetProgramQubit(0): " << gate.targetProgramQubit(0) << ", gate.targetProgramQubit(1): " << gate.targetProgramQubit(1) << endl;
                // cerr << "vTargetQubit[0]: " << vTargetQubit[0] << ", vTargetQubit[1]: " << vTargetQubit[1] << endl;
                newCir.addGate(to_string(i), vTargetQubit);
                if(vQubitLastGate[vTargetQubit[0]] != -1){
                    newCir.addDependency(vQubitLastGate[vTargetQubit[0]], newCir.nGate()-1);
                }
                if(vQubitLastGate[vTargetQubit[1]] != -1){
                    newCir.addDependency(vQubitLastGate[vTargetQubit[1]], newCir.nGate()-1);
                }
                vQubitLastGate[vTargetQubit[0]] = newCir.nGate() - 1;
                vQubitLastGate[vTargetQubit[1]] = vQubitLastGate[vTargetQubit[0]];
            }
        }
    }
    // newCir.printDependency();
}

void Clusterer2::constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ){
    // constructe new device
    unordered_set<pair<unsigned_t, unsigned_t>, pair_hash > sEdge;
    unsigned_t q1, q2;
    for(unsigned_t i = 0; i < oriDevice.nEdge(); ++i){
        Edge& e = oriDevice.edge(i);
        q1 = e.qubitId1();
        q2 = e.qubitId2();
        if(vFinerQ2CoarserQ[q1] < oriDevice.nQubit() && vFinerQ2CoarserQ[q2] < oriDevice.nQubit()){
            if(vFinerQ2CoarserQ[q1] < vFinerQ2CoarserQ[q2]){
                sEdge.emplace(vFinerQ2CoarserQ[q1], vFinerQ2CoarserQ[q2]);
            }
            else if(vFinerQ2CoarserQ[q1] > vFinerQ2CoarserQ[q2]){
                sEdge.emplace(vFinerQ2CoarserQ[q2], vFinerQ2CoarserQ[q1]);
            }
        }
    }
    newDevice.setQubit(mvCoarserQ2FinerQ.size());
    for(pair<unsigned_t, unsigned_t> pEdge : sEdge){
        newDevice.addEdge(pEdge.first, pEdge.second);
    }
}


void Clusterer2::consructBoostGraph(my_graph& graph, Device& device){
    for(unsigned_t i = 0; i < device.nEdge(); ++i){
        Edge& e = device.edge(i);
        add_edge(e.qubitId1(), e.qubitId2(), EdgeProperty(1), graph);    
    }
}

void Clusterer2::printMatching(my_graph& graph, vector<V>& mate, Weight sum){
    fprintf(stdout, "[Info] Found a matching:                        \n");
    fprintf(stdout, "[Info]         matching size: %lu                       \n", matching_size(graph, &mate[0]));
    if(sum > 0)
        fprintf(stdout, "[Info]         matching weight sum: %4f                       \n", (double_t)sum);
    fprintf(stdout, "[Info]         matching:                       \n");
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            fprintf(stdout, "                        {%lu, %lu}                       \n", v, mate[v]); 
        }
    }
    cout << endl;
}


void Clusterer2::printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    fprintf(stdout, "[Info] Coarser qubit -> finer qubit:                        \n");
    for (unsigned_t i = 0; i < mvCoarserQ2FinerQ.size(); ++i){
        fprintf(stdout, "         %d -> ", i);
        for (unsigned_t j : mvCoarserQ2FinerQ[i]){
            fprintf(stdout, "%d ", j);
        }
        fprintf(stdout, "\n");
    }
}

MOLSQ_NAMESPACE_CPP_END
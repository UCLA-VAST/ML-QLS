/***********************************************************************
  File        [ clusterer.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2022 ]
***********************************************************************/

#include "clusterer/clusterer.hpp"

MOLSQ_NAMESPACE_CPP_START

void Clusterer::clusterPorgramQubit(Circuit& oriCir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    // queko 15
    // vector<pair<double_t, pair<unsigned_t, unsigned_t> > > vpWeightEdge{ make_pair(1, make_pair(0,4)), make_pair(1, make_pair(1,2)),
    //                                                                      make_pair(1, make_pair(6,12)), make_pair(1, make_pair(8,13)), 
    //                                                                      make_pair(1, make_pair(9,15)), make_pair(1, make_pair(10, 11)), 
    //                                                                      make_pair(1, make_pair(3,14)), make_pair(1, make_pair(5,7)) };
    // queko 35
    // vector<pair<double_t, pair<unsigned_t, unsigned_t> > > vpWeightEdge{ make_pair(1, make_pair(9,10)), make_pair(1, make_pair(4,13)),
    //                                                                      make_pair(1, make_pair(8,14)), make_pair(1, make_pair(0,2)), 
    //                                                                      make_pair(1, make_pair(3,5)), make_pair(1, make_pair(12, 15)), 
    //                                                                      make_pair(1, make_pair(1,11)), make_pair(1, make_pair(6,7)) };
    
    /*
    vector<pair<double_t, pair<unsigned_t, unsigned_t> > > vpWeightEdge;
    constructQAG(oriCir, vpWeightEdge);
    printEdgeWeight(vpWeightEdge);
    // coarser two qubits based on vpWeightEdge
    // collect finer to coaerser mapping
    unsigned_t qId = 0, countFinerQubit = 0;
    vector<unsigned_t> vFinerQ2CoarserQ(oriCir.nProgramQubit(), oriCir.nProgramQubit());
    vector<bool> vIsCluster(oriCir.nProgramQubit(), 0);
    vector<unsigned_t> vFinerQ;
    for (pair<double_t, pair<unsigned_t, unsigned_t> > & pairWeightEdge : vpWeightEdge) {
        // if both endpoints haven't been cluster, then cluster them
        if(!vIsCluster[pairWeightEdge.second.first] && !vIsCluster[pairWeightEdge.second.second]){
            // cerr << "cluster qubit " << qId << ": " << pairWeightEdge.second.first << ", " << pairWeightEdge.second.second << endl;
            mvCoarserQ2FinerQ[qId] = vector<unsigned_t>();
            vIsCluster[pairWeightEdge.second.first] = 1;
            vIsCluster[pairWeightEdge.second.second] = 1;
            mvCoarserQ2FinerQ[qId].emplace_back(pairWeightEdge.second.first);
            mvCoarserQ2FinerQ[qId].emplace_back(pairWeightEdge.second.second);
            vFinerQ2CoarserQ[pairWeightEdge.second.first] = qId;
            vFinerQ2CoarserQ[pairWeightEdge.second.second] = qId;
            ++qId;
            countFinerQubit += 2;
        }
    }
    // cluster qubits that haven't been cluster
    unsigned_t bound = 3;
    for (pair<double_t, pair<unsigned_t, unsigned_t> > & pairWeightEdge : vpWeightEdge) {
        if(!vIsCluster[pairWeightEdge.second.first] && mvCoarserQ2FinerQ[vFinerQ2CoarserQ[pairWeightEdge.second.second]].size() < bound){
            // cerr << "cluster qubit " << vFinerQ2CoarserQ[pairWeightEdge.second.second] << " add finer qubit " << pairWeightEdge.second.first << " (because " << pairWeightEdge.second.second << " )"<< endl;
            mvCoarserQ2FinerQ[vFinerQ2CoarserQ[pairWeightEdge.second.second]].emplace_back(pairWeightEdge.second.first);
            vFinerQ2CoarserQ[pairWeightEdge.second.first] = vFinerQ2CoarserQ[pairWeightEdge.second.second];
            vIsCluster[pairWeightEdge.second.first] = 1;
            ++countFinerQubit;
        }
        if(!vIsCluster[pairWeightEdge.second.second] && mvCoarserQ2FinerQ[vFinerQ2CoarserQ[pairWeightEdge.second.first]].size() < bound){
            // cerr << "cluster qubit " << vFinerQ2CoarserQ[pairWeightEdge.second.first] << " add finer qubit " << pairWeightEdge.second.second << " (because " << pairWeightEdge.second.first << " )"<< endl;
            mvCoarserQ2FinerQ[vFinerQ2CoarserQ[pairWeightEdge.second.first]].emplace_back(pairWeightEdge.second.second);
            vFinerQ2CoarserQ[pairWeightEdge.second.second] = vFinerQ2CoarserQ[pairWeightEdge.second.first];
            vIsCluster[pairWeightEdge.second.second] = 1;
            ++countFinerQubit;
        }
    }
    for (unsigned_t i = 0; i < vIsCluster.size(); ++i){
        if(!vIsCluster[i]){
            for (unsigned_t j = 0; j < vIsCluster.size(); ++j){
                if(i != j){
                    if(!vIsCluster[j]){
                        mvCoarserQ2FinerQ[qId] = vector<unsigned_t>();
                        vIsCluster[i] = 1;
                        vIsCluster[j] = 1;
                        mvCoarserQ2FinerQ[qId].emplace_back(i);
                        mvCoarserQ2FinerQ[qId].emplace_back(j);
                        vFinerQ2CoarserQ[i] = qId;
                        vFinerQ2CoarserQ[j] = qId;
                        ++qId;
                        countFinerQubit += 2;
                        break;
                    }
                }
            }
        }
        if(!vIsCluster[i]){
            for (unsigned_t j = 0; j < vIsCluster.size(); ++j){
                if(i != j){
                   if(mvCoarserQ2FinerQ[vFinerQ2CoarserQ[j]].size() < bound){
                        mvCoarserQ2FinerQ[vFinerQ2CoarserQ[j]].emplace_back(i);
                        vFinerQ2CoarserQ[i] = vFinerQ2CoarserQ[j];
                        vIsCluster[i] = 1;
                        ++countFinerQubit;
                        break;
                    }
                }
            }
        }
    }
    cerr << "countFinerQubit: " << countFinerQubit << endl;
    // check correctness
    assert(countFinerQubit == oriCir.nProgramQubit());
    for (unsigned_t i = 0; i < mvCoarserQ2FinerQ.size(); ++i){
        for (unsigned_t j : mvCoarserQ2FinerQ[i]){
            // cerr << "i " << i << ", vFinerQ2CoarserQ[mvCoarserQ2FinerQ[i][j]: " << vFinerQ2CoarserQ[j] << ", mvCoarserQ2FinerQ[i][j]: " << j << endl;
            assert(i == vFinerQ2CoarserQ[j]);
        }
    }
    */

    /*
    // construct boost graph
    my_graph graph(oriCir.nProgramQubit());
    vector<vector<bool>> vAdjacent(oriCir.nProgramQubit(), vector<bool>(oriCir.nProgramQubit(), 0));
    consructBoostGraph(graph, oriCir, vAdjacent);
    vector<V> mate(num_vertices(graph));
    // use boost library to construct maximum card matching 
    // see: https://www.boost.org/doc/libs/1_76_0/libs/graph/doc/maximum_matching.html
    edmonds_maximum_cardinality_matching(graph, &mate[0]);
    if(_verbose > 1 ){
        printMatching(graph, mate);
    }
    // construct new circuit 
    vector<unsigned_t> vFinerQ2CoarserQ(oriCir.nProgramQubit(), oriCir.nProgramQubit());
    extractCoarserQubit(oriCir, graph, mate, mvCoarserQ2FinerQ, vFinerQ2CoarserQ, vAdjacent);
    */

    // construct optimal clutering for queko_54_15_1
    vector<unsigned_t> vFinerQ2CoarserQ(oriCir.nProgramQubit(), oriCir.nProgramQubit());
    if(oriCir.nProgramQubit() == 54){
        mvCoarserQ2FinerQ[0] = {36, 39, 48};
        mvCoarserQ2FinerQ[1] = {2, 13, 23};
        mvCoarserQ2FinerQ[2] = {34, 37};
        mvCoarserQ2FinerQ[3] = {20, 30, 44};
        mvCoarserQ2FinerQ[4] = {11, 14, 15};
        mvCoarserQ2FinerQ[5] = {40, 41};
        mvCoarserQ2FinerQ[6] = {3, 45};
        mvCoarserQ2FinerQ[7] = {25, 47};
        mvCoarserQ2FinerQ[8] = {9, 32};
        mvCoarserQ2FinerQ[9] = {35, 52, 53};
        mvCoarserQ2FinerQ[10] = {28, 42};
        mvCoarserQ2FinerQ[11] = {26, 31};
        mvCoarserQ2FinerQ[12] = {43, 19};
        mvCoarserQ2FinerQ[13] = {5, 21};
        mvCoarserQ2FinerQ[14] = {0, 4, 38};
        mvCoarserQ2FinerQ[15] = {22, 46};
        mvCoarserQ2FinerQ[16] = {12, 49};
        mvCoarserQ2FinerQ[17] = {17, 51};
        mvCoarserQ2FinerQ[18] = {1, 29};
        mvCoarserQ2FinerQ[19] = {6, 18};
        mvCoarserQ2FinerQ[20] = {7, 16};
        mvCoarserQ2FinerQ[21] = {33, 50};
        mvCoarserQ2FinerQ[22] = {8, 27};
        mvCoarserQ2FinerQ[23] = {10, 24};
    }
    else if(oriCir.nProgramQubit() == 24){
        mvCoarserQ2FinerQ[0] = {0, 1};
        mvCoarserQ2FinerQ[1] = {2, 3};
        mvCoarserQ2FinerQ[2] = {4, 5};
        mvCoarserQ2FinerQ[3] = {6, 7};
        mvCoarserQ2FinerQ[4] = {8, 9};
        mvCoarserQ2FinerQ[5] = {10, 11};
        mvCoarserQ2FinerQ[6] = {12, 13};
        mvCoarserQ2FinerQ[7] = {14, 15};
        mvCoarserQ2FinerQ[8] = {16, 17};
        mvCoarserQ2FinerQ[9] = {18, 19};
        mvCoarserQ2FinerQ[10] = {20, 21};
        mvCoarserQ2FinerQ[11] = {22, 23};
    }
    else {
        mvCoarserQ2FinerQ[0] = {0, 1};
        mvCoarserQ2FinerQ[1] = {2, 3};
        mvCoarserQ2FinerQ[2] = {4, 6};
        mvCoarserQ2FinerQ[3] = {5, 10};
        mvCoarserQ2FinerQ[4] = {7, 8};
        mvCoarserQ2FinerQ[5] = {9, 11};
    }

    for(auto& pair_key_value : mvCoarserQ2FinerQ){
        for(auto & i : pair_key_value.second){
            vFinerQ2CoarserQ[i] = pair_key_value.first;
        }
    }
    
    constructNewCircuit(oriCir, newCir, vFinerQ2CoarserQ, mvCoarserQ2FinerQ);
    if(_verbose > 2 ){
        fprintf(stdout, "[Info] Original Circuit:                        \n");
        oriCir.printCircuit();
    }
    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserQ2FinerQ);
    }
    if(_verbose > 0 ){
        fprintf(stdout, "[Info] Coarse Circuit:                        \n");
        newCir.printCircuit();
    }
}

void Clusterer::clusterPhysicalQubit(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    /*
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
    */ 

    // construct optimal clutering for queko_54_15_1
    vector<unsigned_t> vFinerQ2CoarserQ(oriDevice.nQubit(), oriDevice.nQubit());
    if(oriDevice.nQubit() == 54){
        mvCoarserQ2FinerQ[0] = {42, 49, 48};
        mvCoarserQ2FinerQ[1] = {24, 30, 36};
        mvCoarserQ2FinerQ[2] = {12 , 18};
        mvCoarserQ2FinerQ[3] = {0 ,6, 13};
        mvCoarserQ2FinerQ[4] = {25, 31, 37};
        mvCoarserQ2FinerQ[5] = {43, 50};
        mvCoarserQ2FinerQ[6] = {19, 26};
        mvCoarserQ2FinerQ[7] = {14, 20};
        mvCoarserQ2FinerQ[8] = {45, 51};
        mvCoarserQ2FinerQ[9] = {38, 39, 44};
        mvCoarserQ2FinerQ[10] = {46 ,52};
        mvCoarserQ2FinerQ[11] = {33, 40};
        mvCoarserQ2FinerQ[12] = {27, 32};
        mvCoarserQ2FinerQ[13] = {15, 21};
        mvCoarserQ2FinerQ[14] = {1, 2, 7};
        mvCoarserQ2FinerQ[15] = {3 , 8};
        mvCoarserQ2FinerQ[16] = {4, 9};
        mvCoarserQ2FinerQ[17] = {10, 16};
        mvCoarserQ2FinerQ[18] = {22, 28};
        mvCoarserQ2FinerQ[19] = {29, 34};
        mvCoarserQ2FinerQ[20] = {47, 53};
        mvCoarserQ2FinerQ[21] = {35, 41};
        mvCoarserQ2FinerQ[22] = {17, 23};
        mvCoarserQ2FinerQ[23] = {5, 11};
    }
    else if(oriDevice.nQubit() == 24){
        mvCoarserQ2FinerQ[0] = {0, 1};
        mvCoarserQ2FinerQ[1] = {2, 3};
        mvCoarserQ2FinerQ[2] = {4, 5};
        mvCoarserQ2FinerQ[3] = {6, 7};
        mvCoarserQ2FinerQ[4] = {8, 9};
        mvCoarserQ2FinerQ[5] = {10, 11};
        mvCoarserQ2FinerQ[6] = {12, 13};
        mvCoarserQ2FinerQ[7] = {14, 15};
        mvCoarserQ2FinerQ[8] = {16, 17};
        mvCoarserQ2FinerQ[9] = {18, 19};
        mvCoarserQ2FinerQ[10] = {20, 21};
        mvCoarserQ2FinerQ[11] = {22, 23};
    }
    else {
        mvCoarserQ2FinerQ[0] = {0, 1};
        mvCoarserQ2FinerQ[1] = {2, 3};
        mvCoarserQ2FinerQ[2] = {4, 6};
        mvCoarserQ2FinerQ[3] = {5, 10};
        mvCoarserQ2FinerQ[4] = {7, 8};
        mvCoarserQ2FinerQ[5] = {9, 11};
    }

    for(auto& pair_key_value : mvCoarserQ2FinerQ){
        for(auto & i : pair_key_value.second){
            vFinerQ2CoarserQ[i] = pair_key_value.first;
        }
    }

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

void Clusterer::constructQAG(Circuit& oriCir, vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge){
    vector<unsigned_t> vQubitLevel;
    unordered_map<pair<unsigned_t , unsigned_t>, double_t, pair_hash> umEdge2Weight;
    unsigned_t gateLevel, i, q0, q1;
    pair<unsigned_t , unsigned_t> edge;
    double_t weight;
    vQubitLevel.resize(oriCir.nProgramQubit(), 0);
    for(i = 0; i < oriCir.nGate(); ++i){
        if(oriCir.gate(i).nTargetQubit()>1){
            q0 = oriCir.gate(i).targetProgramQubit(0);
            q1 = oriCir.gate(i).targetProgramQubit(1);
            gateLevel = min(vQubitLevel[q0], vQubitLevel[q1]); 
            // cout << "q0: " << q0 << ", q1: " << q1 << ", gateLevel: " << gateLevel << endl;
            vQubitLevel[q0] = gateLevel+1;
            vQubitLevel[q1] = gateLevel+1;
            edge = (q0 < q1) ? make_pair(q0,q1) : make_pair(q1,q0); 
            weight = 1 - gateLevel * _programQParam.level_weight_decrease;
            if(weight - 0.1 < 0.001){
                weight = 0.1;
            }
            if (umEdge2Weight.find(edge) == umEdge2Weight.end()){
                umEdge2Weight[edge] = weight;
            }
            else{
                umEdge2Weight[edge] = umEdge2Weight[edge] + weight;
            }
        }
    }
    vpWeightEdge.clear(); 
    vpWeightEdge.reserve(umEdge2Weight.size());
    for (auto & keyItemPair : umEdge2Weight){
        vpWeightEdge.emplace_back(keyItemPair.second, keyItemPair.first);
    }
    sort(vpWeightEdge.rbegin(), vpWeightEdge.rend());
}

void Clusterer::printEdgeWeight(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge){
    fprintf(stdout, "===========================================\n");
    fprintf(stdout, "[INFO] Print Graph Info\n");
    for (pair<double_t, pair<unsigned_t, unsigned_t> >& pWeightEdge : vpWeightEdge){
        fprintf(stdout, "       - Edge(%d, %d), weight: %4f\n", pWeightEdge.second.first, pWeightEdge.second.second, pWeightEdge.first);
    }
    fprintf(stdout, "===========================================\n");
}


void Clusterer::extractCoarserQubit(Circuit& oriCir, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ, vector<vector<bool>>& vAdjacent){
    // construct new qubit by matching result
    unsigned_t countFinerQubit = 0, q1, q2, qId = 0, i, j;
    vector<bool> vIsCluster(oriCir.nProgramQubit(), 0);
    // collect finer to coaerser mapping
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            // cerr << "cluster qubit " << qId << ": " << v << ", " << mate[v] << endl;
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
    // collect the uncluster phy qubit to the other unnclustered qubit
    unsigned_t bound = 3;
    while(countFinerQubit < oriCir.nProgramQubit()){
        for(i = 0; i < oriCir.nProgramQubit(); ++i){
            if(!vIsCluster[i]){
                for(j = i + 1; j < oriCir.nProgramQubit(); ++j){
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
                for(j = 0; j < oriCir.nProgramQubit(); ++j){
                    if(vIsCluster[j] && mvCoarserQ2FinerQ[vFinerQ2CoarserQ[j]].size() < bound){
                        // cerr << "cluster qubit " << vFinerQ2CoarserQ[q1] << " add finer qubit " << i << " (because " << q1 << " )"<< endl;
                        mvCoarserQ2FinerQ[vFinerQ2CoarserQ[j]].emplace_back(i);
                        vFinerQ2CoarserQ[i] = vFinerQ2CoarserQ[j];
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
    assert(countFinerQubit == oriCir.nProgramQubit());
    for (unsigned_t i = 0; i < mvCoarserQ2FinerQ.size(); ++i){
        for (unsigned_t j : mvCoarserQ2FinerQ[i]){
            // cerr << "i " << i << ", vFinerQ2CoarserQ[mvCoarserQ2FinerQ[i][j]: " << vFinerQ2CoarserQ[j] << ", mvCoarserQ2FinerQ[i][j]: " << j << endl;
            assert(i == vFinerQ2CoarserQ[j]);
        }
    }
}


void Clusterer::extractCoarserQubit(Device& oriDevice, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ){
    // construct new qubit by matching result
    unsigned_t countFinerQubit = 0, q1, q2, qId = 0, i, j;
    vector<bool> vIsCluster(oriDevice.nQubit(), 0);
    // collect finer to coaerser mapping
    for (V v : boost::make_iterator_range(vertices(graph))) {
        if (mate[v] != graph.null_vertex() && v < mate[v]) {
            // cerr << "cluster qubit " << qId << ": " << v << ", " << mate[v] << endl;
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
    assert(countFinerQubit == oriDevice.nQubit());
    for (unsigned_t i = 0; i < mvCoarserQ2FinerQ.size(); ++i){
        for (unsigned_t j : mvCoarserQ2FinerQ[i]){
            // cerr << "i " << i << ", vFinerQ2CoarserQ[mvCoarserQ2FinerQ[i][j]: " << vFinerQ2CoarserQ[j] << ", mvCoarserQ2FinerQ[i][j]: " << j << endl;
            assert(i == vFinerQ2CoarserQ[j]);
        }
    }
}


void Clusterer::constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    // construct new circuit by matching result
    newCir.setQubitNum(mvCoarserQ2FinerQ.size());
    vector<unsigned_t> vTargetQubit(2);
    vector<int_t> vQubitLastGate(newCir.nProgramQubit(), -1);
    bool cond1, cond2;
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
             if(cond1 && cond2){
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


void Clusterer::constructNewDevice(Circuit& cir, Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ){
    vector<bool> vIsCluster(oriDevice.nQubit(), 0);
    vector<unsigned_t> vFinerQ2CoarserQ(oriDevice.nQubit(), oriDevice.nQubit());
    unsigned_t countFinerQubit = 0;
    mvCoarserPhyQ2FinerPhyQ.clear();

    // colllect new device clustering based on the circuit mapping
    for(auto & key_value : mvCoarserProQ2FinerProQ){
        mvCoarserPhyQ2FinerPhyQ[key_value.first] = vector<unsigned_t>();
        for(unsigned_t i : key_value.second){
            mvCoarserPhyQ2FinerPhyQ[key_value.first].emplace_back(cir.initialMapping(i));
            vFinerQ2CoarserQ[i] = key_value.first;
            ++countFinerQubit;
        }
    }

    // collect the uncluster phy qubit to its neighbor
    unsigned_t bound = 3, q1;
    while(countFinerQubit < oriDevice.nQubit()){
        for(unsigned_t i = 0; i < oriDevice.nQubit(); ++i){
            if(!vIsCluster[i]){
                Qubit & qubit = oriDevice.qubit(i);
                for( unsigned_t j : qubit.vSpanEdge){
                    Edge & edge = oriDevice.edge(j);
                    q1 = (i == edge.qubitId1()) ? edge.qubitId2() : edge.qubitId1();
                    if(mvCoarserPhyQ2FinerPhyQ[vFinerQ2CoarserQ[q1]].size() < bound){
                        // cerr << "cluster qubit " << vFinerQ2CoarserQ[q1] << " add finer qubit " << i << " (because " << q1 << " )"<< endl;
                        mvCoarserPhyQ2FinerPhyQ[vFinerQ2CoarserQ[q1]].emplace_back(i);
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

    // construct new device
    newDevice.clearEdge();
    constructNewDevice(oriDevice, newDevice, mvCoarserPhyQ2FinerPhyQ, vFinerQ2CoarserQ);
    if(_verbose > 2 ){
        fprintf(stdout, "[Info] Original Device:                        \n");
        oriDevice.printDevice();
    }

    if(_verbose > 1 ){
        printCoarserQ2FinerQ(mvCoarserPhyQ2FinerPhyQ);
    }
    if(_verbose > 0 ){
        fprintf(stdout, "[Info] Coarse Device:                        \n");
        newDevice.printDevice();
    }
}

void Clusterer::constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ){
    // constructe new device
    unordered_set<pair<unsigned_t, unsigned_t> > sEdge;
    unsigned_t q1, q2;
    for(unsigned_t i = 0; i < oriDevice.nEdge(); ++i){
        Edge& e = oriDevice.edge(i);
        q1 = e.qubitId1();
        q2 = e.qubitId2();
        if(vFinerQ2CoarserQ[q1] < vFinerQ2CoarserQ[q2]){
            sEdge.emplace(vFinerQ2CoarserQ[q1], vFinerQ2CoarserQ[q2]);
        }
        else if(vFinerQ2CoarserQ[q1] > vFinerQ2CoarserQ[q2]){
            sEdge.emplace(vFinerQ2CoarserQ[q2], vFinerQ2CoarserQ[q1]);
        }
    }
    newDevice.setQubit(mvCoarserQ2FinerQ.size());
    for(pair<unsigned_t, unsigned_t> pEdge : sEdge){
        newDevice.addEdge(pEdge.first, pEdge.second);
    }
}

void Clusterer::consructBoostGraph(my_graph& graph, Circuit & cir, vector<vector<bool>>& vAdjacent){
    unsigned_t q0, q1;
    for(unsigned_t i = 0; i < cir.nGate(); ++i){
        Gate & gate = cir.gate(i);
        if(gate.nTargetQubit() == 2){
            q0 = gate.targetProgramQubit(0);
            q1 = gate.targetProgramQubit(1);
            if(!vAdjacent[q0][q1]){
                add_edge(q0, q1, EdgeProperty(1), graph);        
                vAdjacent[q0][q1] = 1;
                vAdjacent[q1][q0] = 1;
            }
        }
    }
}

void Clusterer::consructBoostGraph(my_graph& graph, Device& device){
    for(unsigned_t i = 0; i < device.nEdge(); ++i){
        Edge& e = device.edge(i);
        add_edge(e.qubitId1(), e.qubitId2(), EdgeProperty(1), graph);    
    }
}

void Clusterer::printMatching(my_graph& graph, vector<V>& mate, Weight sum){
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


void Clusterer::printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
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
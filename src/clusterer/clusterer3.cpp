/***********************************************************************
  File        [ clusterer.cpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class implementation ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2022 ]
***********************************************************************/

#include "clusterer/clusterer3.hpp"

MOLSQ_NAMESPACE_CPP_START

void Clusterer3::cluster(Circuit & cir, Device & device, vector<Circuit>& vCir, vector<Device>& vDevice, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserProQ2FinerProQ, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserPhyQ2FinerPhyQ, bool allCommute){
    // cluster the first device
    // cir.printCircuitLayout();
    _isAllCommute = allCommute;
    // _isAllCommute = true;
    clusterOneLevel(cir, vCir[0], device, vDevice[0], vmvCoarserProQ2FinerProQ[0], vmvCoarserPhyQ2FinerPhyQ[0]);
    // return;
    for(unsigned_t i = 1; i < vCir.size() && vCir[i-1].nProgramQubit() > 10; ++i){
        clusterOneLevel(vCir[i-1], vCir[i], vDevice[i-1], vDevice[i], vmvCoarserProQ2FinerProQ[i], vmvCoarserPhyQ2FinerPhyQ[i]);
    }
    
}

void Clusterer3::clusterOneLevel(Circuit & cir, Circuit& newCir, Device & device, Device & newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ){
    // calculate affinity between qubit
    unordered_map<pair<unsigned_t, unsigned_t>, unsigned_t, pair_hash> umEdge2Weight;
    int_t proQId = 0, phyQId = 0, p0, p1, q0, q1;
    pair<unsigned_t, unsigned_t> pairEdge;
    vector<unordered_set<unsigned_t>> vsInteractQubit(cir.nProgramQubit());
    // calculate affinity between qubits caused by 2Q gates on them
    for(unsigned_t i = 0; i < cir.nGate(); ++i){
        Gate & gate = cir.gate(i);
        if(gate.nTargetQubit()>1){
            q0 = gate.targetProgramQubit(0);
            q1 = gate.targetProgramQubit(1);
            vsInteractQubit[q0].insert(q1);
            vsInteractQubit[q1].insert(q0);
            pairEdge = (q0 < q1) ? make_pair(q0,q1) : make_pair(q1,q0); 
            if (umEdge2Weight.find(pairEdge) == umEdge2Weight.end()){
                umEdge2Weight[pairEdge] = 10;
            }
            else{
                umEdge2Weight[pairEdge] += 10;
            }
        }
    }
    // add affinity between qubits if they interact with the same other qubit
    unordered_set<unsigned_t>::iterator it;
    unordered_set<unsigned_t>::iterator it2;
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        for(it = vsInteractQubit[i].begin(); it != vsInteractQubit[i].end(); ++it){
            q0 = *it;
            for(it2 = next(it, 1); it2 != vsInteractQubit[i].end(); ++it2){
                q1 = *it2;
                pairEdge = (q0 < q1) ? make_pair(q0,q1) : make_pair(q1,q0); 
                if (umEdge2Weight.find(pairEdge) == umEdge2Weight.end()){
                    umEdge2Weight[pairEdge] = 1;
                }
                else{
                    ++umEdge2Weight[pairEdge];
                }
            }    
        }
    }
    vector<pair<unsigned_t, pair<unsigned_t, unsigned_t>>> vpWeightEdge;
    vpWeightEdge.clear(); 
    vpWeightEdge.reserve(umEdge2Weight.size());
    for (auto & keyItemPair : umEdge2Weight){
        vpWeightEdge.emplace_back(keyItemPair.second, keyItemPair.first);
    }
    sort(vpWeightEdge.rbegin(), vpWeightEdge.rend());
    // cerr << "print qubit weight" << endl;
    // for(auto p : vpWeightEdge){
    //     cerr << "w: " << p.first << ", qubit: " << p.second.first << " " << p.second.second << endl;
    // }
    // cluster qubits
    vector<bool> vPhyIsClustered(device.nQubit(), 0);
    vector<bool> vProIsClustered(cir.nProgramQubit(), 0);
    vector<unsigned_t> vFinerProQ2CarserProQ(cir.nProgramQubit(), 0);
    vector<unsigned_t> vFinerPhyQ2CoarserPhyQ(device.nQubit(), 0);
    vector<int_t> vPhyQ2ProQ(device.nQubit(), -1);
    for(unsigned_t i = 0; i < cir.nProgramQubit(); ++i){
        vPhyQ2ProQ[cir.initialMapping(i)] = i;
    }
    unsigned_t nPhyIsClustered = 0, nProIsClustered = 0;
    vector<unsigned_t> vProId2PhyId;
    // cerr << "in pro qubit clustering 0" << endl;
    for(auto & pairWEdge : vpWeightEdge){
        q0 = pairWEdge.second.first;
        q1 = pairWEdge.second.second;
        p0 = cir.initialMapping(q0);
        p1 = cir.initialMapping(q1);
        if(!vProIsClustered[q0] && !vProIsClustered[q1]){
            if(device.isAdjacent(p0, p1)){
                // cerr << "form coarser qubit " << proQId << " by clustering pro qubit " << q0 << " " << q1 << ", and phy qubit " << p0 << " " << p1 << endl;
                mvCoarserProQ2FinerProQ[proQId] = {static_cast<unsigned_t>(q0), static_cast<unsigned_t>(q1)};
                mvCoarserPhyQ2FinerPhyQ[proQId] = {static_cast<unsigned_t>(p0), static_cast<unsigned_t>(p1)};
                vFinerProQ2CarserProQ[q0] = proQId;
                vFinerProQ2CarserProQ[q1] = proQId;
                vFinerPhyQ2CoarserPhyQ[p0] = proQId;
                vFinerPhyQ2CoarserPhyQ[p1] = proQId;
                nPhyIsClustered += 2;
                nProIsClustered += 2;
                vProIsClustered[q0] = 1;
                vProIsClustered[q1] = 1;
                vPhyIsClustered[p0] = 1;
                vPhyIsClustered[p1] = 1;
                vProId2PhyId.emplace_back(proQId);
                ++proQId;
            }
        }
    }
    // allow size three cluster
    // cerr << "in pro qubit clustering 1" << endl;
    size_t bound = 3;
    if(nProIsClustered < cir.nProgramQubit()){
        for(auto & pairWEdge : vpWeightEdge){
            q0 = pairWEdge.second.first;
            q1 = pairWEdge.second.second;
            p0 = cir.initialMapping(q0);
            p1 = cir.initialMapping(q1);
            if(!vProIsClustered[q0] && device.isAdjacent(p0, p1) && mvCoarserProQ2FinerProQ[vFinerProQ2CarserProQ[q1]].size() <= bound){
                mvCoarserProQ2FinerProQ[vFinerProQ2CarserProQ[q1]].emplace_back(q0);
                mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[p1]].emplace_back(p0);
                vFinerProQ2CarserProQ[q0] = vFinerProQ2CarserProQ[q1];
                vFinerPhyQ2CoarserPhyQ[p0] = vFinerPhyQ2CoarserPhyQ[p1];
                ++nPhyIsClustered;
                ++nProIsClustered;
                vProIsClustered[q0] = 1;
                vPhyIsClustered[p0] = 1;
                // cerr << " add pro qubit" << q0  << " and phy qubit " << p0 << " to coarser qubit " << vFinerProQ2CarserProQ[q1] << endl;
            }
            else if(!vProIsClustered[q1] && device.isAdjacent(p0, p1) && mvCoarserProQ2FinerProQ[vFinerProQ2CarserProQ[q0]].size() <= bound){
                mvCoarserProQ2FinerProQ[vFinerProQ2CarserProQ[q0]].emplace_back(q1);
                mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[p0]].emplace_back(p1);
                vFinerProQ2CarserProQ[q1] = vFinerProQ2CarserProQ[q0];
                vFinerPhyQ2CoarserPhyQ[p1] = vFinerPhyQ2CoarserPhyQ[p0];
                ++nPhyIsClustered;
                ++nProIsClustered;
                vProIsClustered[q1] = 1;
                vPhyIsClustered[p1] = 1;
                // cerr << " add pro qubit" << q1  << " and phy qubit " << p1 << " to coarser qubit " << vFinerProQ2CarserProQ[q0] << endl;
            }
            if(nProIsClustered == cir.nProgramQubit()){
                break;
            }
        }
    }
    int_t neighborQ, bestNeighborQ;
    phyQId = proQId;
    // cerr << "in pro qubit clustering 2" << endl;
    while(nProIsClustered < cir.nProgramQubit()){
        for(unsigned_t i = 0; i < cir.nProgramQubit() && nProIsClustered < cir.nProgramQubit(); ++i){
            if(!vProIsClustered[i]){
                p0 = cir.initialMapping(i);
                Qubit & qubit = device.qubit(p0);
                for(unsigned_t j : qubit.vSpanEdge){
                    Edge & edge = device.edge(j);
                    neighborQ = (edge.qubitId1() == p0) ? edge.qubitId2() : edge.qubitId1();
                    if(!vPhyIsClustered[neighborQ]){
                        mvCoarserPhyQ2FinerPhyQ[phyQId] = {static_cast<unsigned_t>(p0), static_cast<unsigned_t>(neighborQ)};
                        vFinerPhyQ2CoarserPhyQ[p0] = phyQId;
                        vFinerPhyQ2CoarserPhyQ[neighborQ] = phyQId;
                        mvCoarserProQ2FinerProQ[proQId] = {i};
                        // cerr << "form coarser qubit " << phyQId << " by clustering phy qubit " << p0 << " " << neighborQ << endl;
                        q1 = vPhyQ2ProQ[neighborQ];
                        vFinerProQ2CarserProQ[i] = proQId;
                        nPhyIsClustered += 2;
                        ++nProIsClustered;
                        if(q1 > -1){
                            mvCoarserProQ2FinerProQ[proQId].emplace_back(q1);
                            vFinerProQ2CarserProQ[q1] = proQId;
                            ++nProIsClustered;
                        }
                        // cerr << "form coarser qubit " << proQId << " by clustering pro qubit " << i << " " << q1 << endl;
                        vProId2PhyId.emplace_back(phyQId);
                        ++proQId;
                        ++phyQId;
                        vProIsClustered[i] = 1;
                        if(q1 > -1){
                            vProIsClustered[q1] = 1;
                        }
                        vPhyIsClustered[neighborQ] = 1;
                        vPhyIsClustered[p0] = 1;
                        break;
                    }
                }
            }
            if(!vProIsClustered[i]){
                p0 = cir.initialMapping(i);
                Qubit & qubit = device.qubit(p0);
                for(unsigned_t j : qubit.vSpanEdge){
                    Edge & edge = device.edge(j);
                    neighborQ = (edge.qubitId1() == p0) ? edge.qubitId2() : edge.qubitId1();
                    if(mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[neighborQ]].size() <= bound){
                        mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[neighborQ]].emplace_back(p0);
                        vFinerPhyQ2CoarserPhyQ[p0] = vFinerPhyQ2CoarserPhyQ[neighborQ];
                        // cerr << " add phy qubit " << p0 << " to coarser phy qubit " << vFinerPhyQ2CoarserPhyQ[neighborQ] << endl;
                        ++nPhyIsClustered;
                        for(unsigned_t k : mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[neighborQ]]){
                            q1 = vPhyQ2ProQ[k];
                            if(q1 > -1){
                                break;
                            }
                        }
                        if(q1 > -1){
                            mvCoarserProQ2FinerProQ[vFinerPhyQ2CoarserPhyQ[q1]].emplace_back(i);
                            vFinerProQ2CarserProQ[i] = vFinerProQ2CarserProQ[q1];
                            // cerr << " add pro qubit" << i  << " to coarser pro qubit " << vFinerProQ2CarserProQ[q1] << endl;
                        }
                        else{
                            mvCoarserProQ2FinerProQ[proQId] = {i};
                            vFinerProQ2CarserProQ[i] = proQId;
                            vProId2PhyId.emplace_back(vFinerPhyQ2CoarserPhyQ[neighborQ]);
                            // cerr << "form coarser qubit " << proQId << " by clustering pro qubit " << q0 << endl;
                            ++proQId;
                        }
                        ++nProIsClustered;
                        vProIsClustered[i] = 1;
                        vPhyIsClustered[p0] = 1;
                        break;
                    }
                }
            }
        }
        ++bound;
    }
    // cerr << "in pure phy qubit clustering" << endl;
    if(nPhyIsClustered < device.nQubit()){
        for(unsigned_t i = 0; i < device.nQubit() && nPhyIsClustered < device.nQubit(); ++i){
            if(!vPhyIsClustered[i]){
                Qubit & qubit = device.qubit(i);
                bestNeighborQ = -1;
                for(unsigned_t j : qubit.vSpanEdge){
                    Edge & edge = device.edge(j);
                    neighborQ = (edge.qubitId1() == i) ? edge.qubitId2() : edge.qubitId1();
                    if(!vPhyIsClustered[neighborQ]){
                        bestNeighborQ = neighborQ;
                        break;
                    }
                    else if(bestNeighborQ == -1){
                        bestNeighborQ = neighborQ;
                    }
                    else if(mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[bestNeighborQ]].size() > mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[neighborQ]].size()){
                        bestNeighborQ = neighborQ;
                    }
                }
                if(!vPhyIsClustered[bestNeighborQ]){
                    mvCoarserPhyQ2FinerPhyQ[phyQId] = {i, static_cast<unsigned_t>(bestNeighborQ)};
                    mvCoarserPhyQ2FinerPhyQ[phyQId] = {i, static_cast<unsigned_t>(bestNeighborQ)};
                    vFinerPhyQ2CoarserPhyQ[i] = phyQId;
                    vFinerPhyQ2CoarserPhyQ[bestNeighborQ] = phyQId;
                    vPhyIsClustered[i] = 1;
                    vPhyIsClustered[bestNeighborQ] = 1;
                    // cerr << "form coarser phy qubit " << phyQId << " by clustering phy qubit " << i << " " << bestNeighborQ << endl;
                    nPhyIsClustered += 2;
                    ++phyQId;
                }
                else{
                    mvCoarserPhyQ2FinerPhyQ[vFinerPhyQ2CoarserPhyQ[bestNeighborQ]].emplace_back(i);
                    vFinerPhyQ2CoarserPhyQ[i] = vFinerPhyQ2CoarserPhyQ[bestNeighborQ];
                    vPhyIsClustered[i] = 1;
                    // cerr << " add phy qubit " << i << " to coarser phy qubit " << vFinerPhyQ2CoarserPhyQ[neighborQ] << endl;
                    ++nPhyIsClustered;
                }
            }
        }
    }
    // cerr << "print device mvFinerPhyQ2CoarserPhyQ:" << endl;
    // printCoarserQ2FinerQ(mvCoarserPhyQ2FinerPhyQ);
    // cerr << "device before clustering: " << endl;
    // device.printDevice();
    constructNewDevice(device, newDevice, mvCoarserPhyQ2FinerPhyQ, vFinerPhyQ2CoarserPhyQ);
    // cerr << "device after clustering: " << endl;
    // newDevice.printDevice();
    // getchar();
    // cerr << "print circuit mvCoarserQ2FinerQ:" << endl;
    // printCoarserQ2FinerQ(mvCoarserProQ2FinerProQ);
    // cerr << "circuit before clustering: " << endl;
    // cir.printCircuitLayout();
    constructNewCircuit(cir, newCir, vFinerProQ2CarserProQ, mvCoarserProQ2FinerProQ);
    for(unsigned_t i = 0; i < newCir.nProgramQubit(); ++i){
        newCir.setInitialMapping(i, vProId2PhyId[i]);
    }
    // cerr << "circuit after clustering: " << endl;
    // newCir.printCircuitLayout();
    // getchar();
}

void Clusterer3::constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
    // construct new circuit by matching result
    newCir.setQubitNum(mvCoarserQ2FinerQ.size());
    vector<unsigned_t> vTargetQubit(2);
    vector<int_t> vQubitLastGate(newCir.nProgramQubit(), -1);
    vector<int_t> vCoarserGateId(oriCir.nGate(), -1);
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
                vCoarserGateId[i] = newCir.nGate();
                newCir.addGate(to_string(i), vTargetQubit);
                // if(vQubitLastGate[vTargetQubit[0]] != -1){
                //     newCir.addDependency(vQubitLastGate[vTargetQubit[0]], newCir.nGate()-1);
                // }
                // if(vQubitLastGate[vTargetQubit[1]] != -1){
                //     newCir.addDependency(vQubitLastGate[vTargetQubit[1]], newCir.nGate()-1);
                // }
                vQubitLastGate[vTargetQubit[0]] = newCir.nGate() - 1;
                vQubitLastGate[vTargetQubit[1]] = vQubitLastGate[vTargetQubit[0]];
            }
        }
    }
    for(pair<unsigned_t, unsigned_t>& pDependency : (*oriCir.pvpGateDependency())){
        if(vCoarserGateId[pDependency.first] > -1 && vCoarserGateId[pDependency.second] > -1){
            newCir.addDependency(vCoarserGateId[pDependency.first], vCoarserGateId[pDependency.second]);
        }
    }
    // newCir.printDependency();
}

void Clusterer3::constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ){
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
    // calculate qubit distance
    double_t dis;
    for(unsigned_t i = 0; i < newDevice.nQubit(); ++i){
        for(unsigned_t j = i+1; j < newDevice.nQubit(); ++j){
            dis = 0;
            for(unsigned_t p : mvCoarserQ2FinerQ[i]){
                for(unsigned_t pp : mvCoarserQ2FinerQ[j]){
                    dis += oriDevice.getDistance(p, pp);
                }
            }
            dis /= static_cast<double_t>(mvCoarserQ2FinerQ[i].size() * mvCoarserQ2FinerQ[j].size());
            newDevice.setDistance(i, j, dis);
            newDevice.setDistance(j, i, dis);
        }   
    }
    // fprintf(stdout, "[Info] Device Qubit Distance Info                              \n");
    // for(unsigned_t i = 0; i < newDevice.nQubit(); ++i){ 
    //     for(unsigned_t j = i; j < newDevice.nQubit(); ++j){ 
    //         fprintf(stdout, "        - (%d,%d): dis %.2f\n", i, j, newDevice.getDistance(i, j));
    //     }
    // }
}


void Clusterer3::printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ){
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
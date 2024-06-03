/***********************************************************************
  File        [ clusterer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2023 ]
***********************************************************************/
#ifndef CLUSTERER2_HPP
#define CLUSTERER2_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include <fstream>
#include <map>
#include <set>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

MOLSQ_NAMESPACE_HPP_START

using Weight = // boost::multiprecision::cpp_dec_float_50;
    boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<50>,
        boost::multiprecision::et_off >;

using EdgeProperty =
    boost::property<boost::edge_weight_t, Weight>;

using my_graph =
    boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::undirectedS,
        boost::no_property,
        EdgeProperty>;

using V = boost::graph_traits<my_graph>::vertex_descriptor;
using E = boost::graph_traits<my_graph>::edge_descriptor;



class Clusterer2{
    public:
        Clusterer2(unsigned_t v = 0): _verbose(v){};
        ~Clusterer2() {};

        void cluster(Circuit & cir, Device & device, vector<Circuit>& vCir, vector<Device>& vDevice, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserProQ2FinerProQ, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserPhyQ2FinerPhyQ, bool allCommute);


    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void clusterFinestDevice(Circuit & cir, Device & device, Device & newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);
        void collectFinestQubitRegion(Circuit & cir, Device & device, unordered_set<unsigned_t> & sQubitRegion);
        void consructBoostGraph(my_graph& graph, Device& device, unordered_set<unsigned_t>& sQubitRegion);
        void clusterCircuit(Circuit & cir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, unsigned_t nPhyQ);
        void clusterPhysicalQubit(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);

        void extractCoarserQubit(Device& oriDevice, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ);
        void constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ);
        void constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ);
        void consructBoostGraph(my_graph& graph, Device& device);
        void printMatching(my_graph& graph, vector<V>& mate, Weight sum = -1);
        void printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);

        void clusterCircuitNaive(Circuit & cir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, unsigned_t nPhyQ);
    ////////////////////////////
    // Private member
    ////////////////////////////
        unsigned_t _verbose;
        bool       _isAllCommute;
    

};

MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

/***********************************************************************
  File        [ clusterer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2023 ]
***********************************************************************/
#ifndef CLUSTERER_HPP
#define CLUSTERER_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
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

class Clusterer{
    public:
        Clusterer(unsigned_t v = 0): _verbose(v){};
        ~Clusterer() {};
        void clusterPorgramQubit(Circuit& oriCir, Circuit& newCir, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);
        void clusterPhysicalQubit(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);
        void constructNewDevice(Circuit& cir, Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ);

    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void constructQAG(Circuit& oriCir, vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge);
        void printEdgeWeight(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge);
        void extractCoarserQubit(Circuit& oriCir, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ, vector<vector<bool>>& vAdjacent);
        void extractCoarserQubit(Device& oriDevice, my_graph& graph, vector<V>& mate, map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ);
        void constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ);
        void constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ);
        void consructBoostGraph(my_graph& graph, Circuit& cir, vector<vector<bool>>& vAdjacent);
        void consructBoostGraph(my_graph& graph, Device& device);
        void printMatching(my_graph& graph, vector<V>& mate, Weight sum = -1);
        void printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);

    ////////////////////////////
    // Private struct
    ////////////////////////////
        struct ProgramQParam {
            double_t   level_weight_decrease = 0.05;  // reduce at most 5 SWAP count per subcircuit
            double_t   min_weight = 0.1;  // reduce at most 5 SWAP count per subcircuit
        } _programQParam;

    ////////////////////////////
    // Private member
    ////////////////////////////
        unsigned_t _verbose;
    

};

MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

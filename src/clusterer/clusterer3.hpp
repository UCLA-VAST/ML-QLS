/***********************************************************************
  File        [ clusterer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ clusterer ]
  Synopsis    [ clusterer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 11, May., 2023 ]
***********************************************************************/
#ifndef CLUSTERER3_HPP
#define CLUSTERER3_HPP

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
 

class Clusterer3{
    public:
        Clusterer3(unsigned_t v = 0): _verbose(v){};
        ~Clusterer3() {};

        void cluster(Circuit & cir, Device & device, vector<Circuit>& vCir, vector<Device>& vDevice, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserProQ2FinerProQ, vector<map<unsigned_t, vector<unsigned_t>>>& vmvCoarserPhyQ2FinerPhyQ, bool allCommute);

    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void clusterOneLevel(Circuit & cir, Circuit& newCir, Device & device, Device & newDevice, map<unsigned_t, vector<unsigned_t>>& mvCoarserProQ2FinerProQ, map<unsigned_t, vector<unsigned_t>>& mvCoarserPhyQ2FinerPhyQ);
        void constructNewCircuit(Circuit& oriCir, Circuit& newCir, vector<unsigned_t>& vFinerQ2CoarserQ, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ);
        void constructNewDevice(Device& oriDevice, Device& newDevice, map<unsigned_t, vector<unsigned_t> >& mvCoarserQ2FinerQ, vector<unsigned_t>& vFinerQ2CoarserQ);
        void printCoarserQ2FinerQ(map<unsigned_t, vector<unsigned_t>>& mvCoarserQ2FinerQ);
    ////////////////////////////
    // Private member
    ////////////////////////////
        unsigned_t _verbose;
        bool       _isAllCommute;
    

};

MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

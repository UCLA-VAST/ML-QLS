/***********************************************************************
  File        [ placer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ Placer ]
  Synopsis    [ placer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 5, Sep., 2023 ]
***********************************************************************/
#ifndef PLACER_HPP
#define PLACER_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <unordered_set>

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


class Placer{
    public:
        Placer(unsigned_t s = 300, bool constructInitialMappingForFirstCir = 1)
        : _pCircuit(nullptr), _pDevice(nullptr), _verbose(2), _timer(s), _restrictMovement(0){
        }
        ~Placer() {}

        void setTimeout(unsigned_t s){
            _timer.setTimeout(s);
        }
        void setRestrictMovement() { _restrictMovement = true; };
        void run(Circuit& cir, Device& device, bool givenInitialMapping = true);
        // void run(Circuit& cir, Device& device, vector<vector<unsigned_t>>&  vvGateBlock);


    ////////////////////////////
    // Struct RouterParam
    ////////////////////////////
    private:
        struct QGGParam {
            unsigned_t   max_gate_level                = 20;
            double_t     gate_decay_factor             = 0.8;
        } _qggParam;

        struct SAParam{
            /* SA arguments initial value */
            double_t t           = 100000.0;
            double_t t1          = 4.0;
            double_t t_frozen    = 0.000001;
            double_t p           = 0.987; // the initial probability to accept up-hill so-
            int_t l              = 200 ;
            // int_t l              = 10 ;
            // int_t l              = 5 ;
            int_t n              = 0;     
            int_t k              = 7;     // use for updating t
            int_t c              = 100;   // use for updating t
            int_t iter_limit     = 1000;
            // int_t iter_limit     = 50;

            /* Other define */
            int_t random_range   = 100000;
            int_t null_size      = 33;

            int_t init_perturb_num = 100;
            int_t uphill_avg_cnt   = 0;;
            double_t uphill_sum    = 0;
            int_t delta_cost_cnt   = 0;
            double_t delta_sum     = 0;
            double_t delta         = 0;

            // param for movement
            double_t prob_move_within_region = 0.7;
            bool allow_overlapping           = false; // if overlapping is allowed, move p0 on top of p1


            bool find_opt_sol = false;
            double_t verticalCostWeight = 1.5;
            tuple<unsigned_t, unsigned_t, unsigned_t> movement = {0, 0, 0};
        } _saParam;

        struct Solution{
            Solution():
                cost(MAX_DOUBLE)
                {
                    vvProQ2PhyQ.clear();
                    vvPhyQ2ProQ.clear();
                }
            vector<vector<unsigned_t> > vvProQ2PhyQ;
            vector<vector<int_t> >      vvPhyQ2ProQ; // -1 means the physical qubit is not occupied by any program qubit
            double_t cost; // update only when the new solution is accepted
        };

    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        TimeUsage                               _timer;
        Device*                                 _pDevice;
        Circuit*                                _pCircuit;
        vector<vector<unsigned_t>>*             _pvvGateBlock;
        unsigned_t                              _verbose;
        bool                                    _firstInit;
        bool                                    _isGivenInitialMapping;
        bool                                    _restrictMovement;
        Solution                                _optSol;
        Solution                                _curSol;
        vector<double_t>                        _vPrecomputeWeight;
        vector<double_t>                        _vCostForEachTime;  // for current sol
        vector<unsigned_t>                      _vCostForEachQubit; // for current sol
        double_t                                _oriCostTime;       // Used in recover function to recover the correct cost if the new sol is rejected
        unsigned_t                              _oriCostQubit0;     // Used in recover function to recover the correct cost if the new sol is rejected
        unsigned_t                              _oriCostQubit1;     // Used in recover function to recover the correct cost if the new sol is rejected
        vector<vector<pair<double_t, pair<unsigned_t, unsigned_t> > > > _vvWeightedEdge;

    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void constructWeightedEdge(unsigned_t curI, vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge);
        void runSA();  // https://dl.acm.org/doi/pdf/10.1145/1055137.1055161
        void setMapping();
        
        // construct weighted edge
        void preComputeWeight();
        void printEdgeWeight(vector<pair<double_t, pair<unsigned_t, unsigned_t> > >& vpWeightEdge); 

        // cal cost
        void getCost(bool isCurSol = true);
        void calCostForTime(unsigned_t i, bool isCurSol = true);
        void calCostForQubit(unsigned_t i, bool isCurSol = true);
    
       // SA
        void initPerturb();
        void initSol();
        void decideQubitInitialMapping();
        void getMovement();
        void updateOptSol();
        void makeMovement();
        void recover();
        void updateT();
        void swapQubit(unsigned_t m, unsigned_t phyId0, unsigned_t phyId1, bool isCurSol = true);
        bool acceptWorseSol();

        void runRefinement();


};
MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

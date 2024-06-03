/***********************************************************************
  File        [ placer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ Placer_dev ]
  Synopsis    [ placer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 5, Sep., 2023 ]
***********************************************************************/
#ifndef PLACER_DEV_HPP
#define PLACER_DEV_HPP

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


class Placer_dev{
    public:
        Placer_dev(unsigned_t s = 300, bool constructInitialMappingForFirstCir = 1)
        : _pCircuit(nullptr), _pDevice(nullptr), _verbose(0), _timer(s), _restrictMovement(0), _considerOneHopCost(0){
            srand(0);
        }
        ~Placer_dev() {
        }

        void setTimeout(unsigned_t s){
            _timer.setTimeout(s);
        }
        void setRestrictMovement(bool b) { _restrictMovement = b; };
        void setProbMoveWithinRegion(double_t p) { _saParam.prob_move_within_region = p; };
        void run(Circuit& cir, Device& device, bool givenInitialMapping = true, bool allowOverlapping = true);
        double_t getOptimalCost() { return _optSol.cost; };
        void enableAllCommute()                                 { _saParam.is_all_commute = true; }
        void disableAllCommute()                                { _saParam.is_all_commute = false; }
        void enableOneHopCost()                                 { _considerOneHopCost = true; }
        void disableOneHopCost()                                { _considerOneHopCost = false; }
        
        // void run(Circuit& cir, Device& device, vector<vector<unsigned_t>>&  vvGateBlock);


    ////////////////////////////
    // Struct SAParam and Solution
    ////////////////////////////
    private:

        struct SAParam{
            /* SA arguments initial value */
            double_t t           = 100000.0;
            double_t t1          = 4.0;
            double_t t_frozen    = 0.000001;
            double_t p           = 0.987; // the initial probability to accept up-hill so-
            int_t l              = 400 ;
            // int_t l              = 10 ;
            // int_t l              = 5 ;
            int_t n              = 0;     
            int_t k              = 7;     // use for updating t
            int_t c              = 100;   // use for updating t
            int_t iter_limit     = 2000;
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
            double_t prob_move_within_region = 0.5;
            double_t weight_for_neighbor_qubit = 0.75;
            double_t overlapping_weight = 0;
            bool allow_overlapping           = true; // if overlapping is allowed, move p0 on top of p1
            bool is_all_commute              = false; // if overlapping is allowed, move p0 on top of p1
            double_t gate_weight_base        = 0.95;
            unsigned_t decay_level           = 1;


            bool find_opt_sol = false;
            tuple<bool, unsigned_t, unsigned_t> movement = {0, 0, 0}; 
            // type 0->swap two program qubits), proid, phyid
            // type 1->move a program qubit to a new phy location, proid, proid
        } _saParam;

        struct Solution{
            Solution():
                cost(MAX_DOUBLE)
                {
                    vProQ2PhyQ.clear();
                    vvPhyQ2ProQ.clear();
                }
            vector<unsigned_t>       vProQ2PhyQ;
            vector<vector<int_t> >   vvPhyQ2ProQ; // -1 means the physical qubit is not occupied by any program qubit
            double_t cost; // update only when the new solution is accepted
        };

    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        TimeUsage                               _timer;
        Device*                                 _pDevice;
        Circuit*                                _pCircuit;
        unsigned_t                              _verbose;
        vector<double_t>                        _vGateWeight;
        bool                                    _firstInit;
        bool                                    _isGivenInitialMapping;
        bool                                    _restrictMovement;
        bool                                    _randomSol;
        bool                                    _considerOneHopCost;
        Solution                                _optSol;
        Solution                                _curSol;
        double_t                                _costForTime; // for current sol
        double_t                                _oriCostTime; // for old sol
        unsigned_t                              _oldPhyQforTargetProQ; 
        vector<pair<pair<unsigned_t, unsigned_t>, double_t>> _vQubitWeight;

    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void runSA();  // https://dl.acm.org/doi/pdf/10.1145/1055137.1055161
        void setMapping();
        
        // cal cost
        void constructGateWeight();
        void getCost(bool isCurSol = true);
        void calCostForTime(bool isCurSol = true);
        bool checkInjectivity();
    
       // SA
        void initPerturb();
        void initSol();
        void decideQubitInitialMapping();
        void getMovement();
        void updateOptSol();
        void makeMovement();
        void recover();
        void updateT();
        void updateOverlappingWeight();
        void swapQubit(unsigned_t proId0, unsigned_t proId1, bool isCurSol = true);
        void moveProQtoPhyQ(unsigned_t proId0, unsigned_t phyId1, bool isCurSol = true);
        bool acceptWorseSol();

        void runRefinement();


};
MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

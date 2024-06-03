/***********************************************************************
  File        [ aRouter.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ router ]
  Synopsis    [ aRouter class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 6, Sep., 2023 ]
***********************************************************************/
#ifndef AROUTER_HPP
#define AROUTER_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include "ds/pqueue.hpp"
#include <unistd.h>
#include <fstream>

// #define DEBUG

MOLSQ_NAMESPACE_HPP_START

class aRouter{
    public:
        aRouter()
        : _pCircuit(nullptr), _pDevice(nullptr), _verbose(0), _pvpGateDependency(nullptr), _curDisSum(0){
        }
        ~aRouter() {}
        
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }
        void enableAllCommute()                                 { _aRouterParam.is_all_commute = true; }
        void disableAllCommute()                                { _aRouterParam.is_all_commute = false; }
        void enableAdmissibleHeuristic()                        { _aRouterParam.use_admissible_heuristic = true; }
        void disableAdmissibleHeuristic()                       { _aRouterParam.use_admissible_heuristic = false; }
        void enableRestrictRegion()                             { _aRouterParam.restrict_region = true; }
        void disableRestrictRegion()                            { _aRouterParam.restrict_region = false; }
        void setGatePerSeacrh(unsigned_t i)                     { _aRouterParam.gate_per_astar = i; }
        void setGCostLimit(unsigned_t i)                        { _aRouterParam.g_cost_limit = (double_t)i; }
        void setExpandProb(double_t i)                          { _aRouterParam.expand_prob = i; }
        void setQubitRegion(Circuit& cir);
        bool run(Circuit& cir, Device& device);
    
    ////////////////////////////
    // Struct RouterParam
    ////////////////////////////
    private:
        struct aRouterParam {
            unsigned_t   timeout                             = 300;
            unsigned_t   node_limit                          = 450000;// 400000;
            unsigned_t   gate_per_astar                      = 40; //40
            unsigned_t   gate_per_astar_base                 = 40; //40
            unsigned_t   gate_dis_per_astar                  = 2000;
            unsigned_t   dependency_graph_depth_limit        = 200;
            unsigned_t   dependency_graph_depth_limit_base   = 200;
            bool         is_all_commute                      = false;
            bool         restrict_region                     = false;
            bool         use_admissible_heuristic            = false;
            double_t     g_cost_limit                        = MAX_DOUBLE;
            double_t     expand_prob                         = 0.05;
            // long int     mem_limit                           = 96636764160; // 90G
            long int     mem_limit                           = 40636764160; // 90G
        } _aRouterParam;
    
    ////////////////////////////
    // Struct aSolution
    ////////////////////////////

        struct aSolution {
            aSolution():
                swapEdge(-1), f_cost(MAX_DOUBLE), g_cost(MAX_DOUBLE), done(false), vAdjs(0), pParent(nullptr)
                {
                    sUnmappedGate.clear();
                    sReadyMappedGate.clear();
                    vExecutedGate.clear();
                    vAdjs.clear();
                }
            aSolution(unordered_set<unsigned_short_t> unmappedGate, unordered_set<unsigned_short_t> readyMappedGate):
                swapEdge(-1), sUnmappedGate(unmappedGate), sReadyMappedGate(readyMappedGate),  f_cost(MAX_DOUBLE), g_cost(MAX_DOUBLE), done(false), vAdjs(0), pParent(nullptr)
                {
                    vExecutedGate.clear();
                    vAdjs.clear();
                }
            ~aSolution() {}
            // size_t                    idx;
            int_t                              swapEdge;
            vector<unsigned_short_t>           vExecutedGate;
            unordered_set<unsigned_short_t>    sUnmappedGate;
            unordered_set<unsigned_short_t>    sReadyMappedGate;
            double_t               f_cost; // f = g + h
            double_t               g_cost;
            bool                   done;
            vector<aSolution*>     vAdjs;
            aSolution*             pParent;
        } ;

        struct aSolutionCmp {
            bool operator () (const aSolution* n1, const aSolution* n2) const {
                return (n1->f_cost > (n2->f_cost + 1e-9));
                // if (fabs(n1->f_cost - n2->f_cost) > 1e-9) return (n1->f_cost > (n2->f_cost + 1e-9));
                // else return (n1->g_cost > (n2->g_cost + 1e-9));
            }
        };

    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        TimeUsage                               _timer;
        Circuit*                                _pCircuit;
        Device*                                 _pDevice;
        vector<pair<unsigned_t, unsigned_t> > * _pvpGateDependency;
        unsigned_t                              _verbose;
        unsigned_t                              _executedGateCount;

        vector<bool>                            _vExecutedGate; // 1 indicates gate i is executed
        vector<unordered_set<unsigned_t> >      _vsGateUncollectedParent;
        vector<unordered_set<unsigned_t> >      _vsGateUnexecutedParent;
        vector<unordered_set<unsigned_t> >      _vsGateChild;
        unordered_set<aSolution*>               _sAllSolution;
        unordered_set<unsigned_t>               _sCurTargetGates;
        unsigned_t                              _curDisSum;
        vector<unsigned_t>                      _vCurMapping;
        vector<unsigned_t>                      _vCurSolMappingPro2Phy;
        vector<int_t>                           _vCurSolMappingPhy2Pro;
        vector<bool>                            _vQubitHasGate;

        vector<unordered_set<int_t> >           _vsQubitRegion;


        

        typedef PairingHeap<aSolution*, aSolutionCmp>                       PQueue_t;
        typedef PQueue_t::point_iterator                                    PQueueIter_t;
        
        PQueue_t                                _priorityQ; // used in run_astar
        // vector<pair<size_t, PQueueIter_t>>      _vPqIter; // ref, pqIter
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        void initialize(aSolution* pSol);
        aSolution* run_astar() ;
        void firstPeeling();
        void collectGatesForSearch();
        void collectCurTargetGatesFromGate(unsigned_t g, unsigned_t curLevel = 0);
        void constructDependencyInfo();
        void collectNeighbor(aSolution* pSol);
        void initNode(aSolution* pSol);
        double_t calHeuristicCost(aSolution* pSol);
        void updateCircuit(aSolution* pSol);
        void cleanAllNode();
        void printSolution(aSolution* pSol);
        void updateCurMappingForSol(aSolution* pSol);
        void swapQubit(int_t eId);
        bool checkMemUsage();
};
MOLSQ_NAMESPACE_HPP_END

#endif // HROUTER_HPP

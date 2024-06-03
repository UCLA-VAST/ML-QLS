/***********************************************************************
  File        [ aRouter_dev.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ router ]
  Synopsis    [ aRouter_dev class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 6, Sep., 2023 ]
***********************************************************************/
#ifndef AROUTER_DEV_HPP
#define AROUTER_DEV_HPP

#include <memory>
#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include "ds/pqueue.hpp"

// #define DEBUG

MOLSQ_NAMESPACE_HPP_START

class aRouter_dev{
    public:
        aRouter_dev()
        : _pCircuit(nullptr), _pDevice(nullptr), _verbose(0), _curDisSum(0){
            _vsGateParent.clear();
            _vsGateChild.clear();
            _vCurMapping.clear();
            _vCurSolMappingPro2Phy.clear();
            _vCurSolMappingPhy2Pro.clear();
            _sGateReadyToMap.clear();
            _mvSwapEdge.clear();
            srand(0);
        }
        ~aRouter_dev() {
        }
        
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }
        void enableAllCommute()                                 { _aRouter_devParam.is_all_commute = true; }
        void disableAllCommute()                                { _aRouter_devParam.is_all_commute = false; }
        void enableRestrictRegion()                             { _aRouter_devParam.restrict_region = true; }
        void disableRestrictRegion()                            { _aRouter_devParam.restrict_region = false; }
        void setExpandProb(double_t i)                          { _aRouter_devParam.expand_prob = i; }
        void setNodeLinit(int_t i)                              { _aRouter_devParam.node_limit = i; _aRouter_devParam.node_num_after_trim = _aRouter_devParam.node_limit * 0.75;}
        void setQubitRegion(Circuit& cir);
        void enableKillIfTimeout()                              { _aRouter_devParam.kill_if_timeout = true; }
        double_t cost()                                         { return _cost; }
        bool run(Circuit& cir, Device& device);
    
    ////////////////////////////
    // Struct RouterParam
    ////////////////////////////
    private:
        struct aRouter_devParam {
            unsigned_t   timeout                             = 180;
            int_t        node_limit                          = 200;//200;// 100000; // 70000;// 400000;
            unsigned_t   node_num_after_trim                 = 175;// 175;//50000;
            bool         is_all_commute                      = false;
            bool         restrict_region                     = false;
            bool         kill_if_timeout                     = false;
            double_t     expand_prob                         = 0.05;
            double_t     weight                              = 0.1; // 0.1;
        } _aRouter_devParam;
    
    ////////////////////////////
    // Struct aSolution
    ////////////////////////////

        struct aSolution {
            explicit aSolution():
                swapEdge(-1), unexecutedGateCount(0), f_cost(MAX_DOUBLE), g_cost(MAX_DOUBLE), done(false), vAdjs(0), pParent(nullptr)
                {
                    vAdjs.clear();
                    vGateState.clear();
                    sReadyMappedGate.clear();
                }
            explicit aSolution(unordered_set<unsigned_t> & readyMappedGate, vector<pair<bool, bool>> & gateState, unsigned_t count, int_t eId):
                swapEdge(eId), unexecutedGateCount(count), vGateState(gateState), sReadyMappedGate(readyMappedGate), f_cost(MAX_DOUBLE), g_cost(MAX_DOUBLE), done(false), vAdjs(0), pParent(nullptr)
                {   
                    vAdjs.clear();
                }
            aSolution(const aSolution&) = delete;
            aSolution& operator=(const aSolution&) = delete;
            ~aSolution(){
                // cerr << "delete node" << endl;
            };
            // size_t                    idx;
            int_t                        swapEdge;
            int_t                        unexecutedGateCount;
            unsigned_t                   oneLevelDependencyGateCount;
            vector<pair<bool, bool>>     vGateState;
            unordered_set<unsigned_t>    sReadyMappedGate;
            double_t               f_cost; // f = g + h
            double_t               g_cost;
            bool                   done;
            vector<shared_ptr<aSolution>>     vAdjs;
            shared_ptr<aSolution>             pParent;
        } ;

        struct aSolutionCmp {
            bool operator () (const shared_ptr<aSolution> n1, const shared_ptr<aSolution> n2) const {
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
        unsigned_t                              _verbose;
        unsigned_t                              _executedGateCount;
        double_t                                _cost;

        vector<unordered_set<unsigned_t> >      _vsGateParent;
        vector<unordered_set<unsigned_t> >      _vsGateChild;
        vector<unordered_set<int_t> >           _vsQubitRegion;
        vector<vector<unsigned_t> >             _vvQubitGate;
        unsigned_t                              _curDisSum;
        vector<double_t>                        _vGateWeight;
        vector<unsigned_t>                      _vCurMapping;
        vector<unsigned_t>                      _vCurSolMappingPro2Phy;
        vector<int_t>                           _vCurSolMappingPhy2Pro;

        unordered_set<unsigned_t>               _sGateReadyToMap;

        unordered_map<shared_ptr<aSolution>, vector<unsigned_t>> _mvSwapEdge;
        vector<shared_ptr<aSolution>>                 _vAllSolution;


        typedef PairingHeap<shared_ptr<aSolution>, aSolutionCmp>            PQueue_t;
        typedef PQueue_t::point_iterator                                    PQueueIter_t;
        
        // vector<pair<size_t, PQueueIter_t>>      _vPqIter; // ref, pqIter
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        bool toTrimState();
        void trimState(unordered_set<shared_ptr<aSolution>> & sNewSolution);
        bool checkFinalState(shared_ptr<aSolution> pSol);
        void initialize(shared_ptr<aSolution> pSol);
        shared_ptr<aSolution> run_astar() ;
        void constructDependencyInfo();
        void collectNeighbor(shared_ptr<aSolution> pSol);
        void initNode(shared_ptr<aSolution> pSol);
        double_t calHeuristicCost(shared_ptr<aSolution> pSol);
        void cleanAllNode();
        void printSolution(shared_ptr<aSolution> pSol);
        void updateCurMappingForSol(shared_ptr<aSolution> pSol);
        void swapQubit(int_t eId);
    
        void updateCircuit(shared_ptr<aSolution> pSol);
        void addSwapGate(unsigned_t edgeId, unsigned_t t );
        void executeGateUnderCurrentMapping(vector<bool> & vExecutedGate, vector<unsigned_t> & vQubitLastTime);
        bool checkMemUsage();
};
MOLSQ_NAMESPACE_HPP_END

#endif // HROUTER_HPP

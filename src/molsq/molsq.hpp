/***********************************************************************
  File        [ molsq.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ mOLSQ ]
  Synopsis    [ mOLSQ class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 24, Feb., 2023 ]
***********************************************************************/
#ifndef MOLSQ_HPP
#define MOLSQ_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <fstream>
#include <map>
#include "clusterer/clusterer2.hpp"
#include "clusterer/clusterer3.hpp"
#include "placer/initialMapper.hpp"
#include "olsq2/olsq2.hpp"
#include "rolsq2/rolsq2.hpp"
#include "stbolsq2/stbolsq2.hpp"
#include "router/aRouter.hpp"
#include "router/aRouter_dev.hpp" 
#include "writer/writer.hpp"
#include <queue>

MOLSQ_NAMESPACE_HPP_START

class mOLSQ{
    public:
        mOLSQ(Circuit& cir, Device& device, unsigned_t timeoutSTBOLSQ2 = 1200, unsigned_t timeoutOLSQ2 = 1200)
        : _circuit(cir), _device(device), _verbose(0), _timeoutSTBOLSQ2(timeoutSTBOLSQ2),
        _timeoutOLSQ2(timeoutOLSQ2){ srand (0); }
        ~mOLSQ() {}
        void run();
        void enableAllCommute() {_molsqParam.is_all_commute = true; }
        void disableAllCommute() {_molsqParam.is_all_commute = false; }


    ////////////////////////////
    // Private function
    ////////////////////////////
    private:
        void preprocessing();
        bool useMultilevel();
        void runVCyclePureNonSMT(Circuit & reverseCir);
        void runVCycleMix(Circuit & reverseCir);
        void runCoarestLevel(Circuit & cir, Device & device);
        void smtRefinment(sTBOLSQ2& stbolsq2, int_t & phyLevel, int_t & proLevel, vector<Circuit>& vCoarseCir, vector<Device>& vCoarseDevice, unsigned_t bound, bool oneRun = 0);
        void finestLevelRefinement(sTBOLSQ2 & stbolsq2, Circuit& cir, Device& device);
        void nonsmtRefinment(int_t & phyLevel, int_t & proLevel, vector<Circuit>& vCoarseCir, vector<Device>& vCoarseDevice, Circuit & reverseCir);
        void flatNonSMT(Circuit & reverseCir);
        void glueMapping(unsigned_t level);
        
        // declusting stage
        void collectQubitRegion(Device& device, Circuit& finerCir, Circuit& coarserCir, bool declusterCir, unsigned_t phyLevel, unsigned_t proLevel);
        void bfsSearch(Device& device, unordered_set<int_t>& sQubitRegion, unsigned_t expansionTime = 1);
        void insertSWAP(vector<vector<pair<unsigned_t, unsigned_t> > >& vvSwap);
        void asapScheduling(Circuit & cir);
        void aRouterThreePass(aRouter_dev& arouter, Device& device, Circuit& cir, bool disableRestrictRegion = false);
        void setQubitMappingViaMatching(my_graph & graph, Circuit & cir, unsigned_t nQubit);
        unsigned_t estimateSwapDepth(vector<unsigned_t>& vInitialMapping, vector<unsigned_t>& vFinalMapping);
        void constructReverseCircuit(Circuit & cir, Circuit& reverseCir);
        void updateCircuit(Circuit & cir1, Circuit & cir2);
        void updateCircuitWithReverse(Circuit & cir, Circuit & reverseCir);

        void verify(Device & device, Circuit & cir);

    ////////////////////////////
    // Struct OLSQ2Param
    ////////////////////////////
    private:
        struct mOLSQParam {
            unsigned_t   coarse_level_qubit                     = 16;
            unsigned_t   coarse_level_gate                      = 50;
            unsigned_t   smt_refinement_qubit_cnt               = 30;
            unsigned_t   smt_refinement_gate_cnt                = 100;
            unsigned_t   threshold_for_initial_mapper           = 200;
            unsigned_t   initialMapper_threshold_qubit          = 50;
            unsigned_t   initialMapper_threshold_gate           = 100;
            unsigned_t   trials                                 = 5;
            bool         use_multilevel                         = true;
            bool         is_all_commute                         = false;
            bool         placer_one_hop_cost                    = false;
            bool         use_reverse_circuit                    = false;
        } _molsqParam;

          struct Node {
            Node(int_t p, int_t idx, int_t qIdx, int_t eIdx, int_t dis): parentIdx(p), parentEIdx(eIdx), idx(idx), qIdx(qIdx), dis(dis){};
            ~Node() {};
            int_t parentIdx;
            int_t parentEIdx;
            int_t idx;
            int_t qIdx;
            int_t dis;
        };
    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        Circuit&                                _circuit;
        Device&                                 _device;
        unsigned_t                              _verbose;

        unsigned_t                              _timeoutSTBOLSQ2;
        unsigned_t                              _timeoutOLSQ2;
        
        unsigned_t                              _timeoutROLSQ2;

        double_t                                _bestSAcost;

        TimeUsage                               _timer;
        vector<map<unsigned_t, vector<unsigned_t>>> _vmvCoarseProQ2FinerProQ;
        vector<map<unsigned_t, vector<unsigned_t>>> _vmvCoarsePhyQ2FinerPhyQ;
        vector<vector<unsigned_t>>                  _vvQubitMapping;
};
MOLSQ_NAMESPACE_HPP_END

#endif // OLSQ_HPP

/***********************************************************************
  File        [ olsq.hpp ]
  System      [ OLSQ2: optimal quantum layout synthesis tool]
  Package     [ olsq ]
  Synopsis    [ OLSQ2 class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef OLSQ2_HPP
#define OLSQ2_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <algorithm>
#include <bitwuzla/cpp/bitwuzla.h>
#include <pb2cnf.h>
#include "misc/bitwuzlaTerminator.hpp"
#include <set>
#include <map>
#include <queue>

MOLSQ_NAMESPACE_HPP_START

using namespace bitwuzla;

class OLSQ2{
    public:
        OLSQ2(Circuit& cir, Device& device)
        : _pCircuit(&cir), _device(device), _swapIdx(cir.nGate()), _verbose(0){
            // cout << "+++++++++++++++++++++=" << &cir << endl;
            // cout << "+++++++++++++++++++++=" << _pCircuit << endl;
            _vpGateDependency.clear();
            _vpGateTimeWindow.clear();
        }
        ~OLSQ2() {}
        bool                          isValidGateIdx(unsigned_t idx)       const { return 0 <= idx && idx < _pCircuit->nGate();}
        bool                          isValidDependencyIdx(unsigned_t idx) const { return 0 <= idx && idx < _vpGateDependency.size();}
        pair<unsigned_t, unsigned_t>& dependency(unsigned_t idx)                 { assert(isValidDependencyIdx(idx)); return _vpGateDependency[idx]; }        
        
        void reset()                                            { _smt.reset(_olsqParam.timeout, 0);  
                                                                    if(_olsqParam.is_transition)
                                                                        initializeTransitionMode(); 
                                                                    else
                                                                        initializeNormalMode();
                                                                    if(_olsqParam.is_optimize_swap)
                                                                        setOptimizeForSwap();}
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }
        void setCircuit(Circuit & cir)                          { _pCircuit = &cir; }
        void setSwapDuration(unsigned_t d)                      { _olsqParam.swap_duration = d; }                 
        void setHeuristicForSwap(unsigned_t bound)              { _olsqParam.is_use_heuristic_bound_for_swap = true; _olsqParam.heuristic_swap_bound = bound; }                 
        void enableGateTimeWindow()                             { _olsqParam.use_window_range_for_gate = 1; }                 
        void enableCollectQubitRegion()                         { _olsqParam.is_collect_qubit_region = 1; }                 
        void setOptimizeForSwap(){ 
            _olsqParam.is_optimize_swap = true; 
        }                 
        void initializeTransitionMode(unsigned_t min_depth = 1){
            _olsqParam.is_transition = true;
            _olsqParam.is_given_depth = true;
            _olsqParam.min_depth = min_depth;
            _olsqParam.max_depth_bit = ceil(log2(min_depth + 1));
            _olsqParam.max_depth = pow(2, _olsqParam.max_depth_bit);
        }
        void initializeNormalMode(unsigned_t min_depth = 0){
            _olsqParam.is_transition = false;
            if (min_depth == 0){
                _olsqParam.is_given_depth = false;
            }
            else{
                _olsqParam.min_depth = min_depth;
                _olsqParam.max_depth_bit = ceil(log2(min_depth + 1));
                _olsqParam.max_depth = pow(2, _olsqParam.max_depth_bit);
            }
        }
        void setMutilevelMode()  { 
            _olsqParam.is_multilevel = true;    
            initializeTransitionMode();
            setOptimizeForSwap();
            enableCollectQubitRegion(); 
        }
        void enableAllCommute() {_olsqParam.is_all_commute = true; }
        void disableAllCommute() {_olsqParam.is_all_commute = false; }
        bool run();
        void dump();
        void setDependency(vector<pair<unsigned_t, unsigned_t> > & vDependencies);
        void printDependency();
        unsigned_t minDepth()     { return _olsqParam.min_depth; };


    
    ////////////////////////////
    // Struct OLSQ2Param
    ////////////////////////////
    private:
        struct OLSQ2Param {
            bool         is_transition                 = true;
            bool         is_optimize_swap              = false;
            bool         is_use_heuristic_bound_for_swap = false;
            bool         is_given_dependency           = false;
            bool         is_given_depth                = false;
            bool         is_collect_qubit_region       = false;
            bool         is_multilevel                 = false;
            bool         is_all_commute                = false;
            bool         use_window_range_for_gate     = false;
            unsigned_t   max_depth                     = 8;  //  always (power of 2), for bit length 
            unsigned_t   max_depth_bit                 = 3;  
            unsigned_t   min_depth                     = 1;
            unsigned_t   max_depth_expand_factor       = 1;
            unsigned_t   heuristic_swap_bound          = 0;
            unsigned_t   timeout                       = 600000; //ms
            unsigned_t   timeout_per_run               = 100000; //ms
            unsigned_t   swap_duration                 = 1;

        } _olsqParam;

    ////////////////////////////
    // Struct smt
    ////////////////////////////
    struct smt {
            smt(): options(), terminator(1000){
                options.set(Option::PRODUCE_MODELS, true);
                pSolver = unique_ptr<Bitwuzla>(new Bitwuzla(options));
                vvPi.clear();
                vTg.clear();
                vvSigma.clear();
            }
            ~smt(){
            }
            void resetTimeState(unsigned_t timeout){
                terminator = RuntimeTerminator(timeout);
                pSolver->configure_terminator(&terminator);
            }
            void reset(unsigned_t timeout, unsigned_t iter){
                if(iter == 0){
                    terminator = RuntimeTerminator(timeout);
                }
                pSolver.reset(new Bitwuzla(options));
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_PRODUCE_MODELS, 1);
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_INCREMENTAL, 1);
                pSolver->configure_terminator(&terminator);
                vvPi.clear();
                vTg.clear();
                vvSigma.clear();
            }
            unique_ptr<Bitwuzla>            pSolver;
            vector<vector<Term>>            vvPi;         // t->qId
            vector<Term>                    vTg;
            vector<vector<Term>>            vvSigma;      // t->qId
            RuntimeTerminator               terminator;
            Options                         options;
        } _smt;
    
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
        TimeUsage                               _timer;
        Circuit*                                _pCircuit;
        Device&                                 _device;
        vector<pair<unsigned_t, unsigned_t> >   _vpGateDependency;
        vector<pair<unsigned_t, unsigned_t> >   _vGateCollision;
        unsigned_t                              _swapIdx;
        unsigned_t                              _iter;
        unsigned_t                              _verbose;
        bool                                    _hasSol;
        vector<pair<unsigned_t, unsigned_t> >   _vpGateTimeWindow; // pair<start time, end time>
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        bool runSMT();

        void generateFormulation();
        void constructVariable();
        void addInjectiveMappingConstraints(unsigned_t boundOffset = 0);
        void addValidTwoQubitGateConstraints(unsigned_t boundOffset = 0); // if boundOffset > 0, it means we increase min_depth so we need to add constraints for partial t
        void addDependencyConstraints();
        void addSwapConstraints(unsigned_t boundOffset = 0);
        void addTransformationConstraints(unsigned_t boundOffset = 0);
        void addDepthConstraints();
        void addSwapCountConstraints(unsigned_t bound);
        bool checkModel();

        bool optimize();
        bool optimizeDepth();
        bool optimizeSwap();
        bool optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun);

        void extractModel();
        void asapScheduling();

        void increaseDepthBound();
        void constructDependency();
        void constructCollisionList();
        void addDependency(Gate& g1, Gate& g2){
            _vpGateDependency.emplace_back(make_pair(g1.idx(), g2.idx()));
        }
        void addDependency(unsigned_t g1, unsigned_t g2){
            _vpGateDependency.emplace_back(make_pair(g1, g2));
        }
        unsigned_t extract_longest_chain();
        void updateSMT(unsigned_t d);

        // construct overconstrained problem
        void constructGateTimeWindow();
        void updateGateTimeWindow(unsigned_t d);
        void printGateTimeWindow();
};
MOLSQ_NAMESPACE_HPP_END

#endif // MOLSQ_HPP

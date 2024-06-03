/***********************************************************************
  File        [ router.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ router ]
  Synopsis    [ router class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef rOLSQ2_HPP
#define rOLSQ2_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <algorithm>
#include <fstream>
#include <bitwuzla/cpp/bitwuzla.h>
#include "misc/bitwuzlaTerminator.hpp"
#include <queue>
#include <pb2cnf.h>
#include <set>

MOLSQ_NAMESPACE_HPP_START

using namespace bitwuzla;

class rOLSQ2{
    public:
        rOLSQ2(Device& device)
        : _device(device), _verbose(1){
            // cout << "+++++++++++++++++++++=" << &cir << endl;
            // cout << "+++++++++++++++++++++=" << _pCircuit << endl;
        }
        ~rOLSQ2() {}
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }                 
        void setMinDepth(unsigned_t d)                          { _rolsqParam.min_depth = d; increaseDepthBound(); }                 
        void setInitialMapping(vector<unsigned_t>& vInitialMapping )  { _pvInitialMapping = &vInitialMapping; }
        void setFinalMapping(vector<unsigned_t>& vFinalMapping )  { _pvFinalMapping = &vFinalMapping; }
        void run(unsigned_t bound);
        vector<pair<unsigned_t, unsigned_t> >& vSwap()          { return _vSwap; }
        unsigned_t minDepth()                                   { return _rolsqParam.min_depth; }                 
    
    ////////////////////////////
    // Struct rOLSQ2Param
    ////////////////////////////
    private:
        struct rOLSQ2Param {
            unsigned_t   max_depth                     = 8;  //  always (power of 2), for bit length 
            unsigned_t   max_depth_bit                 = 3;  
            unsigned_t   min_depth                     = 2;
            unsigned_t   max_depth_expand_factor       = 1;
            unsigned_t   depth_expand_factor           = 2;
            unsigned_t   timeout                       = 1800000; // ms
        } _rolsqParam;

    ////////////////////////////
    // Struct smtParam
    ////////////////////////////
    struct smt {
            smt(): options(), terminator(1000){
                options.set(Option::PRODUCE_MODELS, true);
                pSolver = unique_ptr<Bitwuzla>(new Bitwuzla(options));
                vvPi.clear();
                vvSigma.clear();
                vPiSort.clear();
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
                pSolver->configure_terminator(&terminator);
                vvPi.clear();
                vvSigma.clear();
                vPiSort.clear();
            }
            unique_ptr<Bitwuzla>           pSolver;
            vector<vector<Term >>          vvPi;         // t->qId
            vector<vector<Term >>          vvSigma;      // t->qId
            vector<Sort>                   vPiSort;      // t->qId
            RuntimeTerminator              terminator;
            Options                        options;
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
        TimeUsage                                         _timer;
        vector<unsigned_t> *                              _pvInitialMapping; 
        vector<unsigned_t> *                              _pvFinalMapping; 
        vector<vector<unsigned_t>>                        _vvQubitRegion; 
        vector<set<unsigned_t>>                           _vsQubitRegion; 
        Device&                                           _device;
        vector<pair<unsigned_t, unsigned_t> >             _vSwap;
        vector<map<unsigned_t, unsigned_t> >              _vmPhy2ProIdxChoiceIdx;
        vector<set<unsigned_t > >                         _vsProSwapRegion;
        set<unsigned_t >                                  _sAllSwapRegion;
        unsigned_t                                        _verbose;
        unsigned_t                                        _iter;
        bool                                              _hasSol;
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        void computeQubitRegion();
        void bfsSearch(unsigned_t qId);

        void computeOverlapRegion();
        void computeQubitIdx();

        void runSMT(unsigned_t bound);

        void generateFormulation();
        void constructVariable();
        void addInitialMappingConstraints();
        void addFinalMappingConstraints();
        void addInjectiveMappingConstraints(unsigned_t boundOffset = 0);
        void addMappingRegionConstraints(unsigned_t boundOffset = 0);
        void addSwapConstraints(unsigned_t boundOffset = 0);
        void addTransformationConstraints(unsigned_t boundOffset = 0);
        void addSwapCountConstraints(unsigned_t bound);
        bool checkModel();

        bool optimize(unsigned_t bound);
        bool optimizeDepth();
        bool optimizeSwapForDepth(unsigned_t lower_swap_bound, unsigned_t upper_swap_bound, bool firstRun);

        void extractModel();

        void increaseDepthBound();
        void updateSMT(unsigned_t d);
        void expandQubitRegion();
};
MOLSQ_NAMESPACE_HPP_END

// todo list: 1. add more constraints when depth increase, 2. graduately relax qubit mapping region

#endif // OLSQ_HPP
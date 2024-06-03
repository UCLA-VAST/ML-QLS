/***********************************************************************
  File        [ olsq.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ clusterer ]
  Synopsis    [ initialMapper class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef INITIALMAPPER_HPP
#define INITIALMAPPER_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include "misc/bitwuzlaTerminator.hpp"
#include <algorithm>
#include <bitwuzla/cpp/bitwuzla.h>
#include <set>
#include <map>
#include <pb2cnf.h>
#include <queue>
#include "placer/placer_dev.hpp"
// #include "placer/placer.hpp"

MOLSQ_NAMESPACE_HPP_START

using namespace bitwuzla;

class InitialMapper{
    public:
        InitialMapper()
        : _pCircuit(0), _pDevice(0), _verbose(0), _bestDisSum(MAX_UNSIGNED){
            srand(0);
        }
        ~InitialMapper() {}
        void reset()                                            { _smt.reset(_initialMapperParam.timeout, 0); }
        void setVerbose(unsigned_t verbose)                     { _verbose = verbose; }
        void setTimeout(unsigned_t t)                           { _initialMapperParam.timeout = t; reset(); }
        bool run(Circuit& cir, Device& device, unsigned_t solNum = 1);
        void enableAllCommute() {_initialMapperParam.is_all_commute = true; }
        void disableAllCommute() {_initialMapperParam.is_all_commute = false; }
        void setCircuitIntialMappingBySolutionIdx(unsigned_t idx, Circuit & cir);
        double_t getOptimalCost() { return _bestDisSum; }

    
    ////////////////////////////
    // Struct initialMapperParam
    ////////////////////////////
    private:
        struct initialMapperParam {
            unsigned_t   timeout                       = 1000000; //ms
            unsigned_t   timeout_per_run               = 100000; //ms
            bool         is_all_commute                = 0;
        } _initialMapperParam;

    ////////////////////////////
    // Struct smt
    ////////////////////////////
    struct smt {
            smt(): options(), terminator(1000){
                options.set(Option::SEED, (uint64_t)0);
                options.set(Option::PRODUCE_MODELS, true);
                pSolver = unique_ptr<Bitwuzla>(new Bitwuzla(options));
                vPi.clear();
            }
            ~smt(){
            }

            void reset(unsigned_t timeout, unsigned_t iter = 0){
                if(iter == 0){
                    terminator = RuntimeTerminator(timeout);
                }
                pSolver.reset(new Bitwuzla(options));
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_PRODUCE_MODELS, 1);
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_INCREMENTAL, 1);

                pSolver->configure_terminator(&terminator);
                vPi.clear();
            }
            unique_ptr<Bitwuzla>           pSolver;
            Sort                           piSort;
            vector<Term>                   vPi;         // t->qId
            vector<Term>                   vGate;         // t->qId
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
        TimeUsage                               _timer;
        Circuit*                                _pCircuit;
        Device*                                 _pDevice;
        double_t                                _bestDisSum;
        unsigned_t                              _verbose;
        unsigned_t                              _solNum;
        bool                                    _hasSol;
        vector<pair<double_t, vector<unsigned_t>>> _vvQubitMapingsolution;
    ////////////////////////////
    // Private functions
    ////////////////////////////
    private:
        bool runSMT();
        void runPostprocessing();

        void generateFormulation();
        void constructVariable();
        void addInjectiveMappingConstraints();
        void addValidTwoQubitGateConstraints(unsigned_t gateId);
        void addGateCountConstraints(unsigned_t bound);
        
        bool checkModel();

        bool optimize();

        void extractModel(bool blockSol = 0);
        void blockSolution(vector<unsigned_t> & vMapping);
        void sampleSolution();

        void mapGateBasedOnInitialMapping();

};
MOLSQ_NAMESPACE_HPP_END

#endif // INITIALMAPPER_HPP

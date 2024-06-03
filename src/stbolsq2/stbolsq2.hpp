/***********************************************************************
  File        [ olsq.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ sTBOLSQ2 ]
  Synopsis    [ sTBOLSQ2 class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 27, Dec., 2022 ]
***********************************************************************/
#ifndef STBOLSQ2_HPP
#define STBOLSQ2_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include "misc/timeUsage.hpp"
#include <fstream>
#include <set>
#include <map>
#include <algorithm>
#include <iterator>
#include <bitwuzla/cpp/bitwuzla.h>
#include "misc/bitwuzlaTerminator.hpp"
#include <pb2cnf.h>

MOLSQ_NAMESPACE_HPP_START

using namespace bitwuzla;

class sTBOLSQ2{
  public:
    sTBOLSQ2(unsigned_t verbose = 1): _pCircuit(0), _pDevice(0), _verbose(verbose), _bestDisSum(MAX_UNSIGNED){};
    sTBOLSQ2(Circuit* pCir, Device* pDevice, unsigned_t verbose = 1)
        : _pCircuit(pCir), _pDevice(pDevice), _verbose(verbose), _bestDisSum(MAX_UNSIGNED), _changeCntBound(0){
			_pvpGateDependency = _pCircuit->pvpGateDependency(); 
      _vvQubitMapping.clear();
    }
    ~sTBOLSQ2() {};
    bool run(unsigned_t bound, bool considerDis = true);
    void dumpSMT();
    void setCircuit(Circuit* pCir, bool constructDependency = false) {
      	_pCircuit = pCir; 
        _vvQubitMapping.clear(); };
    void setDevice(Device* pDevice) {
            _pDevice = pDevice; 
    };
    unsigned_t minDepth() { return _stbolsqParam.min_depth; };
    void setDepth(unsigned_t min_depth = 1){
        _stbolsqParam.min_depth = min_depth;
        _stbolsqParam.max_depth_bit = ceil(log2(min_depth + 1));
        _stbolsqParam.max_depth = pow(2, _stbolsqParam.max_depth_bit);
    }

    vector<vector<unsigned_t> > & vvQubitMapping() { return _vvQubitMapping; };
    vector<unsigned_t> & vFinalQubitMapping() { return _vvQubitMapping.back(); };
    vector<unsigned_t> & vQubitMapping(unsigned_t i) { return _vvQubitMapping[i]; };
    unsigned_t changeCntBound() { return _changeCntBound; };
    unsigned_t bestDisSum() { return _bestDisSum; };
    void setChangeCntBound(unsigned_t bound) { _changeCntBound = bound; };
    void enableFixInitialMapping() { _stbolsqParam.fixInitialMapping = 1; };
    void disableFixInitialMapping() { cerr << "_stbolsqParam.fixInitialMapping: " << _stbolsqParam.fixInitialMapping << endl; _stbolsqParam.fixInitialMapping = 0; };
    void setInitialMapping(vector<unsigned_t>& vInitialMapping )  { _pvInitialMapping = &vInitialMapping; }
    void resetLookaheadCir() { _vpLookaheadCir.clear(); };
    void addLookaheadCir(Circuit & cir) { _vpLookaheadCir.emplace_back(&cir); };
    void enableAllCommute() {_stbolsqParam.is_all_commute = true; }
    void disableAllCommute() {_stbolsqParam.is_all_commute = false; }

  ////////////////////////////
  // Struct OLSQ2Param
  ////////////////////////////
  private:
    struct sTBOLSQParam {
        unsigned_t   max_depth                     = 4;  //  always (power of 2), for bit length 
        unsigned_t   max_depth_bit                 = 2;  
        unsigned_t   min_depth                     = 1;
        unsigned_t   max_depth_expand_factor       = 1;
        unsigned_t   timeout                       = 3600000; //ms
        unsigned_t   timeout_short                 = 100000; //ms
        unsigned_t   timeout_per_run               = 3600000; //ms
        unsigned_t   timeout_per_run_base          = 1800000; //ms
        unsigned_t   timeout_per_run_increase_step = 100000; //ms
        bool         considerDis                   = 1;
        bool         restrictQubitMoving           = 1; // restrict that the qubit moving distance is one for each transition.
        bool         fixInitialMapping             = 0; 
        bool         is_all_commute                  = 0; 
        double_t     lookAheadCostDecay            = 0.5;
    } _stbolsqParam;

  struct smt {
            smt(): options(), terminator(1000){
                options.set(Option::PRODUCE_MODELS, true);
                pSolver = unique_ptr<Bitwuzla>(new Bitwuzla(options));
                vvPi.clear();
                vTg.clear();
                vvChangeCnt.clear();
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
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_PRODUCE_MODELS, 1);
                // bitwuzla_set_option(pSolver, BITWUZLA_OPT_INCREMENTAL, 1);

                pSolver->configure_terminator(&terminator);
                vvPi.clear();
                vTg.clear();
                vvChangeCnt.clear();
                vPiSort.clear();
            }
            unique_ptr<Bitwuzla>     pSolver;
            vector<vector<Term> >    vvPi;         // t->qId
            vector<Term>             vTg;
            vector<vector<Term> >    vvChangeCnt;      // t->qId
            vector<Sort>             vPiSort;      // t->qId
            RuntimeTerminator        terminator;
            Options                  options;
        } _smt;



  private:
    TimeUsage                               _timer;
    Circuit*                                _pCircuit;
    Device*                                 _pDevice; 
    unsigned_t                              _verbose;
    double_t                                _bestDisSum;
    unsigned_t                              _iter;
    unsigned_t                              _changeCntBound;
    bool                                    _hasSol;
    vector<pair<unsigned_t, unsigned_t> > * _pvpGateDependency;
    vector<vector<unsigned_t>>              _vvQubitRegion; 
    vector<map<unsigned_t, unsigned_t> >    _vmPhy2ProIdxChoiceIdx;
    vector<vector<unsigned_t> >             _vvQubitMapping;
    vector<unsigned_t> *                    _pvInitialMapping; 
    vector<Circuit*>                        _vpLookaheadCir; 
    vector<unordered_set<int_t>>            _vsQubitRegion;
    
  ////////////////////////////
  // Private functions
  ////////////////////////////
  private:
    void runSMT(unsigned_t bound);

    void generateFormulation();
    void constructVariable();
    void initQubitRegion();
    void computeQubitIdx();
    void computeOverlapRegion();
    void expandQubitRegion();
    void addInjectiveMappingConstraints(unsigned_t boundOffset = 0);
    void addInitialMappingConstraints();
    void addDependencyConstraints();
    void addDepthConstraints();

    void addQubitMovingDistanceConstraints(unsigned_t boundOffset = 0);
    // change pi related constriants
    void addChangePiConstraint(unsigned_t boundOffset = 0);
    void addChangeCountConstraints(unsigned_t change_bound);
	// restricted formulation
    void addMappingRegionConstraints(unsigned_t boundOffset = 0);
    void addRestrictedValidTwoQubitGateConstraints(unsigned_t boundOffset = 0);

    bool checkModel();

    bool optimize(int_t change_bound);
    bool optimizeDepth();
    void optimizeChangeCnt(int_t change_bound);
    void blockSolution();
    bool sampleSolution(unsigned_t solCntBound = 300);

    bool extractModel(bool blockSol = 0);
    void increaseDepthBound();
    void updateSMT(unsigned_t d);
};

MOLSQ_NAMESPACE_HPP_END

#endif
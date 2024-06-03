/***********************************************************************
  File        [ device.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis tool]
  Package     [ device ]
  Synopsis    [ device class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 22, Nov., 2022 ]
***********************************************************************/
#ifndef DEVICE_HPP
#define DEVICE_HPP

#include "misc/global.hpp"
#include <unordered_map>
#include <unordered_set>
#include <deque>

MOLSQ_NAMESPACE_HPP_START


////////////////////////////
// Struct Qubit
////////////////////////////
struct Qubit {
    Qubit():
        idx(0), vSpanEdge(0)
        {}
    Qubit(unsigned_t i):
        idx(i), vSpanEdge(0)
        {}
    ~Qubit() {}
    unsigned_t            idx; 
    vector<unsigned_t>    vSpanEdge;
};

////////////////////////////
// Struct Edge
////////////////////////////
struct Edge {
    Edge():
        idx(0), prEndPoints(0,0)
        {}
    Edge(unsigned_t i, unsigned_t p0, unsigned_t p1):
        idx(i), prEndPoints(p0,p1)
        {}
    ~Edge() {}
    unsigned_t  qubitId1()   const { return prEndPoints.first; }
    unsigned_t  qubitId2()   const { return prEndPoints.second; }
    unsigned_t            idx;
    pair<unsigned_t, unsigned_t>   prEndPoints;
};



class Device{
    public:
        Device(const string& name) :   
            _name(name), _nQubit(0)
            {
                _vEdge.clear();
                _vvAdjacentMatrix.clear();
            }
        Device(const string& name, unsigned_t nQubit, unsigned_t nEdge) :   
            _name(name), _nQubit(nQubit)
            {
                constructQubit();
                _vEdge.clear();
                _vvAdjacentMatrix.resize(nQubit, vector<bool>(nQubit, false));
            }
        ~Device() {}
        
        string             name()                               const { return _name; }
        unsigned_t         nQubit()                             const { return _nQubit; }
        unsigned_t         nEdge()                              const { return _vEdge.size(); }
        
        bool               isValidQubitIdx(unsigned_t idx)      const { return 0 <= idx && idx < _nQubit;}
        bool               isValidEdgeIdx(unsigned_t idx)       const { return 0 <= idx && idx < _vEdge.size();}


        Qubit&             qubit(unsigned_t idx)  { assert(isValidQubitIdx(idx)); return _vQubit[idx]; }
        Edge&              edge(unsigned_t idx)   { assert(isValidEdgeIdx(idx)); return _vEdge[idx]; }

        void setQubit(unsigned_t nQubit) { _nQubit = nQubit; constructQubit(); _vvAdjacentMatrix.clear(); _vvAdjacentMatrix.resize(nQubit, vector<bool>(nQubit, false)); };
        void setEdge(vector<pair<unsigned_t, unsigned_t> >& vEdge);
        void addEdge(unsigned_t q0, unsigned_t q1);

        void clearEdge(){ _vEdge.clear(); }

        void printDevice();

        void calAPSP();
        double_t getDistance(unsigned_t i, unsigned_t j){
            assert(isValidQubitIdx(j) && isValidQubitIdx(j));
            return _vvASAP[i][j];
        }
        void setDistance(unsigned_t i, unsigned_t j, double_t d){
            assert(isValidQubitIdx(j) && isValidQubitIdx(j));
            _vvASAP[i][j] = d;
        }
        bool isAdjacent(unsigned_t i, unsigned_t j) { assert(isValidQubitIdx(j) && isValidQubitIdx(j)); return _vvAdjacentMatrix[i][j]; }

    private:
        ////////////////////////////
        // Private member
        ////////////////////////////
        string          _name;
        unsigned_t      _nQubit;
        vector<Qubit>   _vQubit;
        vector<Edge>    _vEdge;

        vector<vector<double_t> >      _vvASAP;
        vector<vector<bool> >       _vvAdjacentMatrix;

        ////////////////////////////
        // Private functions
        ////////////////////////////
        void constructQubit();
};


MOLSQ_NAMESPACE_HPP_END

#endif //DEVICE_HPP
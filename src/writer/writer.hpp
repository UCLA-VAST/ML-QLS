/***********************************************************************
  File        [ writer.hpp ]
  System      [ mOLSQ: multilevel quantum layout synthesis toolool]
  Package     [ writer ]
  Synopsis    [ writer class header ]
  Author      [ ]
  
  Affiliation [ UCLA ]
  Date        [ 24, Feb., 2023 ]
***********************************************************************/
#ifndef WRITER_HPP
#define WRITER_HPP

#include "misc/global.hpp"
#include "cir/circuit.hpp"
#include "device/device.hpp"
#include <fstream>

MOLSQ_NAMESPACE_HPP_START

class Writer{
    public:
        Writer(Circuit& cir, Device& device)
        : _circuit(cir), _device(device) {}
        ~Writer() {}

        void outputQASM(string const & fileName);
        string outputQASMStr();

    ////////////////////////////
    // Private member
    ////////////////////////////
    private:
        Circuit&                                _circuit;
        Device&                                 _device;
};
MOLSQ_NAMESPACE_HPP_END

#endif // WRITER_HPP

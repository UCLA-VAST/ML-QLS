#include "misc/global.hpp"
#include "molsq/molsq.hpp"
#include "olsq2/olsq2.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>


MOLSQ_NAMESPACE_CPP_START

namespace py = pybind11;


PYBIND11_MODULE(olsqPy, m) {
    // optional module docstring
    m.doc() = "pybind11 olsq plugin";

    // bindings to Circuit class
    py::class_<Circuit>(m, "Circuit")
        .def(py::init<const string &, unsigned_t, unsigned_t>())
        .def("printCircuit", &Circuit::printCircuit, "Print circuit information")
        .def("printCircuitLayout", &Circuit::printCircuitLayout, "Print compiled circuit information")
        .def("printDependency", &Circuit::printDependency, "Print dependency information");

    m.def("addGate", [](Circuit& self, const string& gateName, pybind11::list pyTargetQubit, unsigned_t duration){
        // convert pyTargetQubit into vTargetQubit
        vector<unsigned_t> vTargetQubit;
        for (py::handle item : pyTargetQubit){
            vTargetQubit.emplace_back(item.cast<unsigned_t>());
        }
        return self.addGate(gateName, vTargetQubit, duration); 
    }, "Add a gate to circuit", py::arg("self"), py::arg("gateName"), py::arg("pyTargetQubit"), py::arg("duration") = 1);
    m.def("setInitialMapping", [](Circuit& self, pybind11::list pyInitialMapping){
        vector<unsigned_t> vInitialMapping;
        for (py::handle item : pyInitialMapping){
            vInitialMapping.emplace_back(item.cast<unsigned_t>());
        }
        return self.setInitialMapping(vInitialMapping); 
    }, "Set initial mapping");
    m.def("setDependency", [](Circuit& self, pybind11::list pyDependencies){
        // convert pyDependencies into vDependencies
        vector<pair<unsigned_t, unsigned_t> > vDependencies;
        for (py::handle item : pyDependencies){
            vDependencies.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        }
        return self.setDependency(vDependencies); 
    }, "Set dependency for circuit");

    // bindings to Device class
    py::class_<Device>(m, "Device")
        .def(py::init<const string &, unsigned_t, unsigned_t>())
        .def("nQubit", &Device::nQubit, "get number of physical qubits")
        .def("addEdge", &Device::addEdge, "Add an edge to device")
        .def("printDevice", &Device::printDevice, "Print device information");
    m.def("setEdge", [](Device& self, pybind11::list pyEdge){
        // convert pyEdge into vEdge
        vector<pair<unsigned_t, unsigned_t> > vEdge;
        for (py::handle item : pyEdge){
            vEdge.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        }
        return self.setEdge(vEdge); 
    }, "Set device edges");

    // bindings to mOLSQ class
    py::class_<mOLSQ>(m, "mOLSQ")
        .def(py::init<Circuit&, Device&, unsigned_t, unsigned_t>(), "initialize mOLSQ", py::arg("cir"), py::arg("device") ,py::arg("timeoutSTBOLSQ2") = 120, py::arg("timeoutOLSQ2") = 120)
        .def("run", &mOLSQ::run, "run multilevel OLSQ")
        .def("enableAllCommute", &mOLSQ::enableAllCommute, "enable all commute")
        .def("disableAllCommute", &mOLSQ::disableAllCommute, "disable all commute");

    // bindings to OLSQ2 class
    py::class_<OLSQ2>(m, "OLSQ2")
        .def(py::init<Circuit&, Device&>())
        .def("setSwapDuration", &OLSQ2::setSwapDuration, "use given initial mapping (default = 1)")
        .def("setHeuristicForSwap", &OLSQ2::setHeuristicForSwap, "set a given swap count as SWAP upper bound")
        .def("initializeTransitionMode", &OLSQ2::initializeTransitionMode, "use transition based mode (default)", py::arg("min_depth") = 1)
        .def("initializeNormalMode", &OLSQ2::initializeNormalMode, "use normal based mode", py::arg("min_depth") = 0)
        .def("setOptimizeForSwap", &OLSQ2::setOptimizeForSwap, "Set optimization object to swap count")
        .def("run", &OLSQ2::run, "run quantum layout synthesis")
        .def("dump", &OLSQ2::dump, "dump quantum layout synthesis formulation")
        .def("enableGateTimeWindow", &OLSQ2::enableGateTimeWindow, "dump quantum layout synthesis formulation")
        // .def("setDependency", [](OLSQ& self, pybind11::list pyDependencies){
        //     // convert pyDependencies into vDependencies
        //     vector<pair<unsigned_t, unsigned_t> > vDependencies;
        //     for (py::handle item : pyDependencies){
        //         vDependencies.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        //     }
        //     return self.setDependency(vDependencies); 
        // }, "Set dependency for circuit")
        .def("printDependency", &OLSQ2::printDependency, "Print dependency information");
    m.def("setDependency", [](OLSQ2& self, pybind11::list pyDependencies){
        // convert pyDependencies into vDependencies
        vector<pair<unsigned_t, unsigned_t> > vDependencies;
        for (py::handle item : pyDependencies){
            vDependencies.emplace_back(item.cast<pair<unsigned_t, unsigned_t> >());
        }
        return self.setDependency(vDependencies); 
    }, "Set dependency for circuit");
}

MOLSQ_NAMESPACE_CPP_END

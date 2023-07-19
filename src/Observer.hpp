//An Observer writes out information about a System to file

#ifndef OBSERVER_HPP
#define OBSERVER_HPP

#include <string>
#include <iomanip>
#include <experimental/filesystem>
#include "System.hpp"
#include <H5Cpp.h>
//#include "hdf5"
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace fs = std::experimental::filesystem;

class System;

class Observer
{
private:
    friend class System;

public:
    std::string output_dir = "./"; //where to dump output

    int network_freq=10;
    int thermo_freq=10;
    int do_h5md=1;

    /*** Methods ***/

    //constructor
    Observer(ParamDict &theParams);

    //destructor
    ~Observer();

    //Output
    //Style is that of LAMMPS dump files
    //void dump_nodes(System &theSys, std::string subdir); //write out positions
    //void dump_springs(System &theSys, std::string subdir); //write out topology
    //void dump_thermo(System &theSys, std::string subdir); //write out energy, pressure, etc.

    void open_h5md(System &theSys, std::string subdir);
    void dump_h5md(System &theSys, std::string subdir); //write all data to hdf5 file

};


#endif

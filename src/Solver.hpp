//A Solver advances a System in time according to some (for now, configuration-space) dynamics. 
//For now, it only supports overdamped Langevin dynamics solved via the Euler-Maruyama method.

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <string>
#include "CustomRandom.hpp"
#include "System.hpp"

class Solver
{
public:
    //TODO: allow for alternative dynamics + numerical solvers
    //std::string dynamics;
    //std::string method;
    double dt;
    double gamma; //friction
    double D; //D=kT/gamma (Einstein relation)
    int do_subtract_com; //if 1, subtract center of mass force from random forces
    int do_pin_node; //if 1, pin the first node (don't let its position change)

    gsl_rng *rg;

    /*** Methods ***/

    //constructor
    Solver(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Solver();

    //solve
    void update(System &theSys);

};

#endif
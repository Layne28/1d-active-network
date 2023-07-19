#include "Solver.hpp"

Solver::Solver(ParamDict &theParams, gsl_rng *&the_rg)
{
    dt = 0.005;
    gamma = 1.0;
    D = 1.0;
    do_subtract_com = 0;

    double kT=1.0; //temporary, used to compute D
    if(theParams.is_key("dt")) dt = std::stod(theParams.get_value("dt"));
    if(theParams.is_key("gamma")) gamma = std::stod(theParams.get_value("gamma"));
    if(theParams.is_key("kT")) kT = std::stod(theParams.get_value("kT")); 
    if(theParams.is_key("do_subtract_com")) do_subtract_com = std::stoi(theParams.get_value("do_subtract_com")); 
    //std::cout << dt << std::endl;

    //TODO: check that this kT is the same as that of the System being solved!
    D = kT/gamma;

    //Set RNG
    rg = the_rg;
}

Solver::~Solver() {}

void Solver::update(System &theSys)
{
    std::vector<arma::vec> forces(theSys.N);
    for(int i=0; i<theSys.N; i++)
    {
        arma::vec force = theSys.get_force(theSys.network[i]);
        forces[i] = force;
    }

    for(int i=0; i<theSys.N; i++)
    {
        arma::vec pos = theSys.network[i].get_pos();
        for(int k=0; k<3; k++)
        {
            //TODO: implement alternatives to Euler (RK2, etc.)
            double incr = forces[i](k)/gamma*dt + sqrt(2*D)*gsl_ran_gaussian(rg, sqrt(dt));
            pos(k) += incr;
            theSys.network[i].vel[k] = incr/dt;
        }

        theSys.network[i].set_pos(pos(0));
    }
    theSys.apply_pbc();
    theSys.zero_com();
    theSys.time++;
}

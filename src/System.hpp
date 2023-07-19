//A System consists of a SpringNetwork along with boundary conditions and control parameters (temperature, external field, etc.)

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include "Node.hpp"
#include "Spring.hpp"
#include "ParamDict.hpp"
#include "Observer.hpp"

class Observer;

class System
{
private:
    Observer *obs; //use to write output

public:
    //For now we assume NVT ensemble
    int N; //No. of particles (nodes)
    double rho; //density
    double kT; //temperature
    int dim; //# of spatial dimensions
    double dt; //timestep

    double a; //lattice constant
    int La; //box dimensions in units of a
    double Lx; //box dimensions
    int is_p_x; //periodic or not
    int can_break; //can the bonds break or not
    int repulsion_on; //turn on short-range WCA repulsion
    double repulsive_sigma; //range of repulsive force
    double repulsive_eps; //strength of repulsive force
    double drmax; //parameter for fene potential
    std::string potential_type = "harmonic";

    int time; //No. of timesteps taken (can be reset)

    std::vector<Node> network; //nodes and springs
    std::vector<std::vector<int>> image; //which periodic image each node is in (starts at (0,0,0))
    std::vector<double> noise_flucs;

    /*** Methods ***/

    //constructor
    //TODO: add default ParamDict
    System(ParamDict &theParams);

    //destructor
    ~System();

    //Make changes to System state
    void apply_pbc(); //apply periodic boundary conditions
    void zero_com(); //subtract center of mass from each particle position
    void set_obs(Observer &anObs);

    //Get (some of these could be made static)
    arma::vec get_com();
    double get_energy();
    std::vector<arma::vec> get_forces();
    arma::vec get_force(Node &n);
    arma::vec get_disp_vec(Node &n1, Node &n2);
    double get_dist(Node &n1, Node &n2);
    double get_strain(Spring &s);
    Observer get_obs();

    arma::vec get_wca_force(double r, arma::vec rvec, double sig, double eps);
    arma::vec get_lj_force(double r, arma::vec rvec, double sig, double eps);
    arma::vec get_harmonic_force(double r, arma::vec rvec, double K, double l0);
    arma::vec get_fene_force(double r, arma::vec rvec, double K, double l0, double drmax);

    double get_lj_potential(double r, double sig, double eps, double rc);
    double get_wca_potential(double r, double sig, double eps);
    double get_harmonic_potential(double r, double K, double l0);
    double get_fene_potential(double r, double K, double l0, double drmax);
};

#endif

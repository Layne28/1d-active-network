#include "ActiveSolver.hpp"

ActiveSolver::ActiveSolver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg) : Solver(theParams, the_rg)
{
    if(theParams.is_key("va")) va = std::stod(theParams.get_value("va"));
    //Modify parameters to make network and active noise generator consistent 
    theParams.add_entry("D", std::to_string(va*va));
    double dee_x = theSys.Lx/int(std::stoi(theParams.get_value("nx")));
    //Set precision of dx
    std::ostringstream out;
    out.precision(15);
    out << std::fixed << dee_x;
    theParams.add_entry("dx", out.str());
    anGen = new Generator(theParams, the_rg);

    //std::cout << theParams.get_value("dx") << std::endl;
    //std::cout << theParams.get_value("D") << std::endl;
    
    //Check that active noise gen. params are consistent with particle system
    if(fabs(anGen->Lx-theSys.Lx)>1e-10)
    {
        std::cout << "ERROR! System and Active Noise Generator do not have the same size (" << std::fixed << theSys.Lx << " vs " << anGen->Lx << ")." << std::endl;
        exit(1);
    }
}

ActiveSolver::~ActiveSolver()
{
    delete anGen;
}

void ActiveSolver::update(System &theSys)
{
    //Get conservative forces
    std::vector<arma::vec> potential_forces(theSys.N);
    for(int i=0; i<theSys.N; i++)
    {
        arma::vec force = theSys.get_force(theSys.network[i]);
        potential_forces[i] = force;
    }

    //Get active noise force on each particle
    std::vector<arma::vec> active_forces = get_active_noise_forces(theSys, *anGen); 
    for(int k=0; k<theSys.dim; k++) theSys.noise_flucs[k] = 0.0; 
    //std::cout << "D: " << anGen->D << std::endl; //make this into a unit test

    for(int i=0; i<theSys.N; i++) {
        if(i==0 && do_pin_node) continue; //skip updating position of 1st node
        arma::vec pos = theSys.network[i].get_pos();
        for(int k=0; k<theSys.dim; k++) {
            //Euler step
            theSys.noise_flucs[k] += active_forces[i](k)*active_forces[i](k);
            double incr = potential_forces[i](k)/gamma*dt 
                   + active_forces[i](k)*dt;
            if (theSys.kT>0) incr += sqrt(2*theSys.kT/gamma)*gsl_ran_gaussian(rg, sqrt(dt)); 
            pos(k) += incr; 

            theSys.network[i].pos[k] = pos(k);
            theSys.network[i].vel[k] = incr/dt;
        }
    }
    anGen->step(dt); //advance active noise in time
    theSys.apply_pbc();
    //theSys.zero_com(); //don't do this, messes up with pbc
    theSys.time++;
    for(int k=0; k<theSys.dim; k++) theSys.noise_flucs[k] /= theSys.N;
    //make this into a unit test
    //std::cout << "mean square active force, x, y, z: " << theSys.noise_flucs[0] << " " << theSys.noise_flucs[1] << " " << theSys.noise_flucs[2] << " " << std::endl;
}

void ActiveSolver::update_adaptive(System &theSys, double deet, int level)
{
    //Get conservative forces
    std::vector<arma::vec> potential_forces(theSys.N);
    potential_forces = theSys.get_forces();

    //Get active noise force on each particle
    std::vector<arma::vec> active_forces = get_active_noise_forces(theSys, *anGen); 
    for(int k=0; k<theSys.dim; k++) theSys.noise_flucs[k] = 0.0; 
    //std::cout << "D: " << anGen->D << std::endl; //make this into a unit test

    //Get thermal force on each particle
    std::vector<arma::vec> thermal_forces;
    thermal_forces = get_thermal_forces(theSys, deet);

    //Compute particle motion due to forces
    std::vector<arma::vec> incr(theSys.N);
    for(int i=0; i<theSys.N; i++){
        arma::vec v(theSys.dim,arma::fill::zeros);
        incr[i] = v;
    }
    for(int i=0; i<theSys.N; i++){
        incr[i] = potential_forces[i]/gamma*deet
                  + active_forces[i]*deet //TODO: check whether there should be a factor of 1/gamma here
                  + thermal_forces[i];
    }

    //Check whether the new position will result in a really large force
    for(int i=0; i<theSys.N; i++){
        if(i==0 && do_pin_node) continue; //skip updating position of 1st node
        theSys.network[i].old_pos = theSys.network[i].pos;
        theSys.network[i].pos += incr[i];
    }
    theSys.apply_pbc();
    std::vector<arma::vec> new_forces = theSys.get_forces();
    double max_force = 0;
    for(int i=0; i<theSys.N; i++){
        for(int k=0; k<theSys.dim; k++){
            if(fabs(new_forces[i][k])>max_force) max_force = new_forces[i][k];
        }
    }    

    //only decrease time step if force is above threshold
    //and timestep is not already tiny
    if(max_force > force_thresh && deet>1e-10){
        //if(level==0) std::cout << "Force too high. Decreasing time step by a factor of 4 (now =" << deet/4 << ")." << std::endl;
        //num_decrease_dt;
        //revert to old position
        for(int i=0; i<theSys.N; i++){
            theSys.network[i].pos = theSys.network[i].old_pos;
        }
        for(int k=0; k<4; k++){
            update_adaptive(theSys, deet/4, level+1);
        }
    }
    else{
        for(int i=0; i<theSys.N; i++){
            theSys.network[i].vel = incr[i]/deet;
        }
        anGen->step(deet); //advance active noise in time
        theSys.apply_pbc();
        for(int k=0; k<theSys.dim; k++) theSys.noise_flucs[k] /= theSys.N;
    }

    if(level==0) theSys.time++;
}

std::vector<arma::vec> ActiveSolver::get_active_noise_forces(System &theSys, Generator &gen)
{
    std::vector<arma::vec> active_forces(theSys.N);

    //Initialize to zero
    for(int i=0; i<theSys.N; i++)
    {
        arma::vec v(1,arma::fill::zeros);
        active_forces[i] = v;
    }

    //if active velocity is zero, don't bother doing calculation
    if (va<1e-10) {
        //std::cout << "Active speed is zero." << std::endl;
        return active_forces;
    }

    //Compute real-space noise
    arma::field<arma::cx_vec> xi = gen.get_xi_r(1);
    
    //Assign active force to each particle
    //based on location in noise grid
    //TODO: write tests for this
    for(int i=0; i<theSys.N; i++)
    {
        arma::vec pos = theSys.network[i].get_pos();
        int xind = int((pos(0)+0.5*theSys.Lx)/gen.dx);
        //std::cout << xind << " " << yind << " " << zind << std::endl;
        for(int j=0; j<1; j++) active_forces[i](j) = xi(xind)(j).real();
    }

    //Subtract out center of mass force (velocity)
    if (do_subtract_com==1) {
        arma::vec com_vel(theSys.dim, arma::fill::zeros);
        for(int i=0; i<theSys.N; i++) {
            com_vel += active_forces[i];
        }
        for(int i=0; i<theSys.N; i++) {
            active_forces[i] -= com_vel/theSys.N;
        }
    }

    return active_forces;
}

std::vector<arma::vec> ActiveSolver::get_thermal_forces(System &theSys, double deet){
    
    std::vector<arma::vec> thermal_forces(theSys.N);
    
    //Initialize to zero
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        thermal_forces[i] = v;
    }

    //if temperature is zero, don't bother doing calculation
    if (theSys.kT<1e-10) {
        //std::cout << "kT is zero." << std::endl;
        return thermal_forces; 
    }

    for (int i=0; i<theSys.N; i++) {
        for (int k=0; k<theSys.dim; k++) {
            thermal_forces[i][k] = sqrt(2*theSys.kT/gamma)*gsl_ran_gaussian(rg, sqrt(deet)); 
        }
    }

    //Subtract out center of mass force (velocity)
    if (do_subtract_com==1) {
        arma::vec com_vel(theSys.dim, arma::fill::zeros);
        for(int i=0; i<theSys.N; i++) {
            com_vel += thermal_forces[i];
        }
        for(int i=0; i<theSys.N; i++) {
            thermal_forces[i] -= com_vel/theSys.N;
        }
    }
    
    return thermal_forces;
}

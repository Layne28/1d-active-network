#include "System.hpp"

System::System(ParamDict &theParams)
{
    time = 0;

    //First assign default values to parameters
    kT = 0.0;
    rho = 4.0;
    dim = 1;
    dt = 0.005;
    a = 1.0;
    N = 3;
    La = 3;
    is_p_x = false;
    can_break = 0;
    obs = nullptr;
    repulsion_on = 1;
    repulsive_sigma = 0.3;
    repulsive_eps = 0.1;
    drmax = 0.5;

    //These are temporary, won't need to be accessed again outside of constructor
    std::string node_protocol = "uniform";
    std::string spring_protocol = "uniform";
    double l0 = 1.0;
    double K = 1.0;

    //Then assign from ParamDict if there
    if(theParams.is_key("kT")) kT = std::stod(theParams.get_value("kT"));
    if(theParams.is_key("rho")) rho = std::stod(theParams.get_value("rho"));
    if(theParams.is_key("dim")) dim = std::stoi(theParams.get_value("dim"));
    if(theParams.is_key("dt")) dt = std::stod(theParams.get_value("dt"));
    if(theParams.is_key("a")) a = std::stod(theParams.get_value("a"));
    if(theParams.is_key("N")) N = std::stod(theParams.get_value("N"));
    if(theParams.is_key("La")) La = std::stod(theParams.get_value("La"));
    if(theParams.is_key("kT")) kT = std::stod(theParams.get_value("kT"));
    if(theParams.is_key("is_p_x")) is_p_x = std::stoi(theParams.get_value("is_p_x"));
    if(theParams.is_key("node_protocol")) node_protocol = theParams.get_value("node_protocol");
    if(theParams.is_key("spring_protocol")) spring_protocol = theParams.get_value("spring_protocol");
    if(theParams.is_key("l0")) l0 = std::stod(theParams.get_value("l0"));
    if(theParams.is_key("K")) K = std::stod(theParams.get_value("K"));
    if(theParams.is_key("potential_type")) potential_type = theParams.get_value("potential_type");
    if(theParams.is_key("repulsion_on")) repulsion_on = std::stod(theParams.get_value("repulsion_on"));
    if(theParams.is_key("repulsive_sigma")) repulsive_sigma = std::stod(theParams.get_value("repulsive_sigma"));
    if(theParams.is_key("repulsive_eps")) repulsive_eps = std::stod(theParams.get_value("repulsive_eps"));
    if(theParams.is_key("drmax")) drmax = std::stod(theParams.get_value("drmax"));

    //Compute parameters that are inferred from other parameters
    //TODO: fix for non-cubic 
    Lx = a*La;

    /*
    double noise_arr[dim];
    for(int k=0; k<dim; k++) noise_arr[k] = 0.0;
    noise_flucs = noise_arr;
    */
    for(int k=0; k<dim; k++) noise_flucs.push_back(0.0);
    //Initialize nodes
    if (node_protocol=="zeros")
    {
        for (int i=0; i<N; i++)
        {
            Node n;
            network.push_back(n);
        }
    }
    else if (node_protocol=="uniform")
    {
        for(int i=0; i<N; i++)
        {
            Node n1(a*i-N*a/2);
            network.push_back(n1);
            network[i].old_pos = network[i].pos;
        }
    }
    else if (node_protocol=="random")
    {
        std::cout << "to be implemented" << std::endl;
    }
    else
    {
        throw std::runtime_error("Error: node initialization protocol not supported.");
    }
    this->apply_pbc();
    this->zero_com();

    //Initialize springs
    if (spring_protocol=="uniform")
    {
        //add springs to all nodes within rest length + epsilon
        double eps = 1e-2;
        for (int i=0; i<N-1; i++)
        {
            for (int j=i+1; j<N; j++)
            {
                if (get_dist(network[i],network[j])<(l0+eps))
                {
                    Spring::add_spring(network[i], network[j], K, l0);
                }
            }
        }
    }
    else
    {
        throw std::runtime_error("Error: spring initialization protocol not supported.");
    }

    //Initialize periodic image indices
    for (int i=0; i<N; i++)
    {
        std::vector<int> index(dim, 0);
        image.push_back(index);
    }

    //Update old position to equal position after pbc
    for(int i=0; i<N; i++){
        network[i].old_pos = network[i].pos;
    }

}

System::~System() {}

void System::apply_pbc()
{
    for (int i=0; i<N; i++)
    {
        //arma::vec pos = network[i].get_pos();
        if (is_p_x)
        {
            if (network[i].pos[0]<-0.5*Lx)
            {
                network[i].pos[0] += Lx;
                image[i][0] -= 1;
            }
            if (network[i].pos[0]>=0.5*Lx)
            {
                network[i].pos[0] -= Lx;
                image[i][0] += 1;
            }
        }
    }
}

void System::zero_com()
{
    arma::vec com(dim, arma::fill::zeros);
    for(int i=0; i<N; i++)
    {
        com += network[i].pos;
    }
    for(int i=0; i<N; i++)
    {
        network[i].pos -= com/N;
    }
}

arma::vec System::get_com()
{
    arma::vec com(dim, arma::fill::zeros);
    for(int i=0; i<N; i++)
    {
        com += network[i].pos;
    }
    return com;
}

void System::set_obs(Observer &anObs)
{
    obs = &anObs;
}

//TODO: redo implementation to avoid double-counting
double System::get_energy()
{
    double energy = 0;
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<network[i].get_num_springs(); j++)
        {
            double K = network[i].springs[j].get_stiffness();
            //std::cout << "K: " << K << std::endl;
            double l0 = network[i].springs[j].get_rest_length();
            //std::cout << "l0: " << l0 << std::endl;
            Node *n = network[i].springs[j].node1;
            if ((*n).is_equal(network[i])) n = network[i].springs[j].node2;
            double dist = get_dist(network[i],*n);
            if(potential_type=="harmonic"){
                energy += get_harmonic_potential(dist, K, l0);
            }
            else if(potential_type=="harmonic_wca"){
                energy += get_harmonic_potential(dist, K, l0) + get_wca_potential(dist, repulsive_sigma, repulsive_eps);
            }
            else if(potential_type=="fene"){
                energy += get_fene_potential(dist, K, l0, drmax);
            }
            else{
                std::cout << "WARNING: potential type not recognized!" << std::endl;
            }
            energy += 0.5*K*(dist-l0)*(dist-l0);
            if(repulsion_on) energy += get_wca_potential(dist,repulsive_sigma,repulsive_eps);
        }
    }
    energy *= 0.5; //correct for double-counting bonds
    return energy;
}

std::vector<arma::vec> System::get_forces(){

    std::vector<arma::vec> forces(N);
    for(int i=0; i<N; i++) {
        arma::vec force(dim);
        force = get_force(network[i]);
        forces[i] = force;
    }
    return forces;
}

arma::vec System::get_force(Node &n1)
{
    arma::vec force(1, arma::fill::zeros);
    for(int j=0; j<n1.get_num_springs(); j++)
    {
        double K = n1.springs[j].get_stiffness();
        double l0 = n1.springs[j].get_rest_length();
        Node *n2 = n1.springs[j].node2;
        if(n1.is_equal(*n2)) n2 = n1.springs[j].node1;
        //if(n1.get_id()==(*n2).get_id()) std::cout << "id1: " << n1.get_id() << " id2: " << (*n2).get_id() << std::endl;
        arma::vec disp = get_disp_vec(n1, *n2);
        double dist = get_dist(n1, *n2);
        //std::cout << n1 << std::endl;
        //std::cout << *n2 << std::endl;
        if (dist<1e-15) throw std::runtime_error("ERROR: attempting to divide by zero in force calculation!");

        if (potential_type=="harmonic"){
            force += get_harmonic_force(dist, disp, K, l0);//-K*(dist-l0)*disp/dist; //harmonic force
        }
        else if(potential_type=="harmonic_wca"){
            force += get_harmonic_force(dist, disp, K, l0) + get_wca_force(dist, disp, repulsive_sigma, repulsive_eps);
        }
        else if(potential_type=="fene"){
            force += ;
        }
        else{
            std::cout << "WARNING: potential type not recognized!" << std::endl;
        }
    }
    //std::cout << "force: " << force << std::endl;
    //std::cout << std::endl;
    return force;
}

arma::vec System::get_disp_vec(Node &n1, Node &n2)
{
    //arma::vec pos1 = n1.get_pos();
    //arma::vec pos2 = n2.get_pos();
    arma::vec disp = {0};
    for (int i=0; i<dim; i++) disp[i] = n1.pos[i]-n2.pos[i];

    if (is_p_x)
    {
        if (disp[0]<(-0.5*Lx)) disp[0] += Lx;
        if (disp[0]>=(0.5*Lx)) disp[0] -= Lx;
    }

    return disp;
}

double System::get_dist(Node &n1, Node &n2)
{
    arma::vec disp = get_disp_vec(n1, n2);
    double len = 0;
    for (int i=0; i<dim; i++) len += disp[i]*disp[i];
    len = sqrt(len);

    return len;
}

Observer System::get_obs()
{
    if(obs)
    {
        return *obs;
    }
    else
    {
        throw("Error: no observer set!");
    }
}

double System::get_lj_potential(double r, double sig, double eps, double rc) {

    double rat = sig/r;
    double r2 = rat*rat;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;

    double rat_cut = sig/rc;
    double r2_cut = rat_cut*rat_cut;
    double r6_cut = r2_cut*r2_cut*r2_cut;
    double r12_cut = r6_cut*r6_cut;

    return 4*eps*(r12-r6) - 4*eps*(r12_cut-r6_cut);
}

double System::get_wca_potential(double r, double sig, double eps) {

    double rc = sig*pow(2,1.0/6.0);
    if (r<rc) return get_lj_potential(r, sig, eps, rc);
    else return 0;
}

double System::get_harmonic_potential(double r, double K, double l0) {

    return 0.5*K*(r-l0)*(r-l0);
}

double System::get_fene_potential(double r, double K, double l0, double drmax) {

    return -0.5*K*drmax*drmax*log(1-((r-l0)/drmax)*((r-l0)/drmax));
}

arma::vec System::get_lj_force(double r, arma::vec rvec, double sig, double eps) {

    double rat = sig/r;
    double r2 = rat*rat;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;

    double r14 = r12*r2;
    double r8 = r6*r2;

    return 24*(eps/sig)*(2*r14-r8)*rvec;
}

arma::vec System::get_wca_force(double r, arma::vec rvec, double sig, double eps) {

    arma::vec force = arma::zeros(dim);

    double rc = sig*pow(2,1.0/6.0);
    if (r<rc) return get_lj_force(r, rvec, sig, eps);
    else return force;
}

arma::vec System::get_harmonic_force(double r, arma::vec rvec, double K, double l0) {

    return -K*(r-l0)*rvec/r;
}

arma::vec System::get_fene_force(double r, arma::vec rvec, double K, double l0, double drmax) {

    return -K*(r-l0)/(1-((r-l0)/drmax)*((r-l0)/drmax))*rvec/r;
}

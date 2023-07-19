#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams), obs(theParams), solver(sys, theParams, theGen)
{

    params = theParams;
    sys.set_obs(obs);

    if(theParams.is_key("equil_steps")) equil_steps = std::stoi(theParams.get_value("equil_steps"));
    if(theParams.is_key("production_steps")) production_steps = std::stoi(theParams.get_value("production_steps"));
    if(theParams.is_key("info_freq")) info_freq = std::stoi(theParams.get_value("info_freq"));
    if(theParams.is_key("experiment")) experiment = theParams.get_value("experiment");
}

LabBench::~LabBench() {}

void LabBench::run(int nstps, std::string subdir, int net_freq, int therm_freq)
{
    if(nstps==-1) nstps = this->production_steps;
    if(net_freq==-1) net_freq = this->obs.network_freq;
    if(therm_freq==-1) therm_freq = this->obs.thermo_freq;

    if(obs.do_h5md)
    {
        obs.open_h5md(sys, subdir);
    }

    for(int i=0; i<nstps; i++)
    {
        //Record data
        if (i%info_freq==0) std::cout << "step " << i << std::endl;

        if(obs.do_h5md)
        {
            if(i%net_freq==0)
            {
                obs.dump_h5md(sys, subdir);
            }
        }
/*
        else
        {
            if (i%net_freq==0)
            {
                obs.dump_nodes(sys, subdir);
                //only output topology after first step
                //if it changes
                if(i==0 || sys.can_break==1)
                {
                    obs.dump_springs(sys, subdir);
                }
            }
            if (i%therm_freq==0) obs.dump_thermo(sys, subdir);
        }
*/
        //Advance dynamics
        //solver.update(sys);
        solver.update_adaptive(sys, sys.dt, 0);
    }
}

void LabBench::do_experiment(std::string expt)
{
    if(expt=="standard")
    {
        std::cout << "Running standard experiment." << std::endl;
        this->run_standard_experiment();
    } 
    else if(expt=="swollen") 
    {
        std::cout << "Running swollen experiment." << std::endl;
        this->run_swollen_experiment();
    }
    else if(expt=="compression") 
    {
        std::cout << "Running compression experiment." << std::endl;
        this->run_compression_experiment();
    }
    else if(expt=="shear")
    {
        std::cout << "Running shear experiment." << std::endl;
        this->run_shear_experiment();
    }
    else if(expt=="sine_perturb")
    {
        std::cout << "Running sine perturbation experiment." << std::endl;
        this->run_sine_perturb_experiment();
    }
    else{
        std::cout << "This experiment has not been designed yet.\n" << std::endl;
        exit(0);
    }
}

void LabBench::run_standard_experiment()
{
    std::cout << "Equilibrating..." << std::endl;
    this->run(this->equil_steps, "/equil", this->obs.network_freq, this->obs.thermo_freq);

    std::cout << "Doing production run..." << std::endl;
    this->run(this->production_steps, "/prod", this->obs.network_freq, this->obs.thermo_freq);
}

void LabBench::run_swollen_experiment()
{
    double factor = 1.0;
    if(this->params.is_key("swell_factor"))
    {
        factor = std::stod(this->params.get_value("swell_factor"));
    }
    std::cout << "Swelling/shrinking network by multiplying natural bond length by a factor " << factor << std::endl;
    for(int i=0; i<this->sys.N; i++)
    {
        for(int j=0; j<this->sys.network[i].get_num_springs(); j++)
        {
            double l0 = this->sys.network[i].springs[j].get_rest_length();
            this->sys.network[i].springs[j].set_rest_length(l0*factor);
        }
    }

    this->run_standard_experiment();

}


void LabBench::run_compression_experiment()
{

}

void LabBench::run_shear_experiment()
{

}
void LabBench::run_sine_perturb_experiment()
{
/*
    double A = 0.1;
    if(this->params.is_key("amplitude"))
    {
        A = std::stod(this->params.get_value("amplitude"));
    }
    std::cout << "Running sine perturbation experiment with amplitude " << A << "..." << std::endl;
    for(int i=0; i<this->sys.N; i++)
    {
        arma::vec pos_curr = sys.network[i].get_pos();
        pos_curr(1) += A*sin(2*M_PI*pos_curr(0)/sys.Lx);
        sys.network[i].set_pos(pos_curr(0),pos_curr(1),pos_curr(2));
    }
    this->run(this->production_steps, "/", this->obs.network_freq, this->obs.thermo_freq);

*/
}

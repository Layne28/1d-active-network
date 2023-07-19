#include "Observer.hpp"

Observer::Observer(ParamDict &theParams)
{
    if(theParams.is_key("output_dir")) output_dir = theParams.get_value("output_dir") + "/";
    if(theParams.is_key("network_freq")) network_freq = std::stoi(theParams.get_value("network_freq"));
    if(theParams.is_key("thermo_freq")) thermo_freq = std::stoi(theParams.get_value("thermo_freq"));

    fs::create_directories(output_dir);
}

Observer::~Observer() {}

void Observer::open_h5md(System &theSys, std::string subdir)
{
    //Create an empty h5md file for storing the trajectory
    using namespace HighFive;
    fs::create_directories(output_dir + subdir);
    std::string name = output_dir + subdir + "/traj.h5";
    std::cout << name << std::endl;
    if(fs::exists(name))
    {
        std::cout << "Warning: file already exists. Overwriting..." << std::endl; 
        fs::remove(name);
    }
    File file(name, File::ReadWrite | File::Create | File::Truncate);

    Group h5md = file.createGroup("/h5md");
    Group particles = file.createGroup("/particles");
    Group observables = file.createGroup("/observables");
    Group parameters = file.createGroup("/parameters");

    //Subgroups of "observables"
    Group potential_energy = file.createGroup("/observables/potential_energy");
    Group noise_flucs = file.createGroup("/observables/noise_fluctuations");

    //Subgroups of "particles"
    Group all_particles = file.createGroup("particles/all");
    Group box = file.createGroup("/particles/all/box");
    Group position = file.createGroup("/particles/all/position");
    Group velocity = file.createGroup("/particles/all/velocity");
    Group image = file.createGroup("/particles/all/image");

    //Sugroups of "connectivity"
    if(!theSys.can_break) {
        //Get connectivity
        int tot_num_springs = 0;
        for(int i=0; i<theSys.N; i++) {
            tot_num_springs += theSys.network[i].get_num_springs();
        }
        std::vector<std::vector<int>> bond_list(tot_num_springs, std::vector<int>(2,0));
        std::vector<int> bonds_to;
        std::vector<int> bonds_from;

        int cnt = 0;
        for(int i=0; i<theSys.N; i++) {
            for(int j=0; j<theSys.network[i].get_num_springs(); j++) {
                bond_list[cnt][0] = theSys.network[i].springs[j].node1->get_id();
                bond_list[cnt][1] = theSys.network[i].springs[j].node2->get_id();
                cnt++;
            }
        }
        std::sort(bond_list.begin(), bond_list.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
            if (a[0] == b[0]) return a[1]<b[1];
            else return a[0]<b[0];
            //return ((a[0] < b[0]) && (a[1]<b[1]));
            //return a[0] < b[0];
        });
        bond_list.erase( unique( bond_list.begin(), bond_list.end() ), bond_list.end() );
        for(int i=0; i<bond_list.size(); i++) {
            //std::cout << bond_list[i][0] << " " << bond_list[i][1] << std::endl;
            bonds_from.push_back(bond_list[i][0]+1);
            bonds_to.push_back(bond_list[i][1]+1);
        }

        //write connectivity to file
        DataSet all_bonds_from = file.createDataSet<int>("/parameters/vmd_structure/bond_from", DataSpace::From(bonds_from));
        all_bonds_from.write(bonds_from);

        DataSet all_bonds_to = file.createDataSet<int>("/parameters/vmd_structure/bond_to", DataSpace::From(bonds_to));
        all_bonds_to.write(bonds_to);

        Reference myRef1 = Reference(file, all_particles);
        Attribute bonds_to_group = all_bonds_to.createAttribute<Reference>("particles_group", DataSpace::From(myRef1));
        bonds_to_group.write(myRef1);
        Reference myRef2 = Reference(file, all_particles);
        Attribute bonds_from_group = all_bonds_from.createAttribute<Reference>("particles_group", DataSpace::From(myRef2));
        bonds_from_group.write(myRef2);
    }
    else
    {
        std::cout << "i/o for changing connectivity not currently supported." << std::endl;
    }

    //Assets of "box"
    std::vector<std::string> boundary_types(theSys.dim);
    std::vector<double> edges{theSys.Lx, 0.0, 0.0};

    for(int d=0; d<theSys.dim; d++)
    {
        if(d==0)
        {
            if(theSys.is_p_x)
            {
                boundary_types[d] = "periodic";
            }
            else
            {
                boundary_types[d] = "none";
            }
        }
    }

    Attribute box_dimension = box.createAttribute<int>("dimension", DataSpace::From(theSys.dim));
    box_dimension.write(theSys.dim);

    Attribute box_boundary = box.createAttribute<std::string>("boundary", DataSpace::From(boundary_types));
    box_boundary.write(boundary_types);

    DataSet edge_values = file.createDataSet<double>("/particles/all/box/edges", DataSpace::From(edges));
    edge_values.write(edges);

    DataSet a = file.createDataSet<double>("/particles/all/box/lattice_constant", DataSpace::From(theSys.a));
    a.write(theSys.a);
}

void Observer::dump_h5md(System &theSys, std::string subdir)
{
    //Write trajectory data in the h5md format
    using namespace HighFive;
    try
    {
        std::string name = output_dir + subdir + "/traj.h5";
        if(!fs::exists(name)) {
            std::cout << "Error: file does not exist!" << std::endl; 
            exit(0);
        }
        File file(name, File::ReadWrite);

        //Create necessary data structures for storage
        std::vector<std::vector<std::vector<double>>> all_pos(1, std::vector<std::vector<double>>(theSys.N, std::vector<double>(theSys.dim,0.0)));
        std::vector<std::vector<std::vector<double>>> all_vel(1, std::vector<std::vector<double>>(theSys.N, std::vector<double>(theSys.dim,0.0)));
        std::vector<std::vector<std::vector<int>>> all_image(1, std::vector<std::vector<int>>(theSys.N, std::vector<int>(theSys.dim,0)));

        DataSpace part_val_space = DataSpace({1,theSys.N,theSys.dim},{DataSpace::UNLIMITED, theSys.N,theSys.dim});
        DataSpace part_vec_space = DataSpace({1,theSys.dim},{DataSpace::UNLIMITED, theSys.dim});
        DataSpace part_t_space = DataSpace({1},{DataSpace::UNLIMITED});
        DataSetCreateProps props_val;
        props_val.add(Chunking(std::vector<hsize_t>{1,theSys.N,theSys.dim}));
        DataSetCreateProps props_vec;
        props_vec.add(Chunking(std::vector<hsize_t>{1,theSys.dim}));
        DataSetCreateProps props_time;
        props_time.add(Chunking(std::vector<hsize_t>{1}));

        //Fill in positions and velocities
        for(int i=0; i<theSys.N; i++)
        {
            for(int j=0; j<theSys.dim; j++)
            {
                all_pos[0][i][j] = theSys.network[i].pos(j);
                all_vel[0][i][j] = theSys.network[i].vel(j);
                all_image[0][i][j] = theSys.image[i][j];
            }
        }

        //Update position
        if(!file.exist("/particles/all/position/step"))
        {
            //std::cout << "creating position data" << std::endl;
            DataSet step = file.createDataSet<int>("/particles/all/position/step", part_t_space, props_time);
            step.write(theSys.time);
            DataSet time = file.createDataSet<double>("/particles/all/position/time", part_t_space, props_time);
            time.write(theSys.time*theSys.dt);
            DataSet value = file.createDataSet<double>("/particles/all/position/value", part_val_space, props_val);
            value.select({0,0,0},{1,theSys.N,theSys.dim}).write(all_pos);
        }
        else
        {
            //Update step
            DataSet step = file.getDataSet("/particles/all/position/step");
            std::vector<long unsigned int> step_dim = step.getDimensions();
            std::vector<long unsigned int> step_dim_old = step_dim;
            step_dim[0] += 1;
            step.resize(step_dim);
            step.select(step_dim_old,{1}).write(theSys.time);

            //Update time
            DataSet time = file.getDataSet("/particles/all/position/time");
            std::vector<long unsigned int> time_dim = time.getDimensions();
            std::vector<long unsigned int> time_dim_old = time_dim;
            time_dim[0] += 1;
            time.resize(time_dim);
            time.select(time_dim_old,{1}).write(theSys.time*theSys.dt);

            //Update position values
            DataSet value = file.getDataSet("/particles/all/position/value");
            std::vector<long unsigned int> value_dim = value.getDimensions();
            std::vector<long unsigned int> value_dim_old = value_dim;
            value_dim[0] += 1;
            value.resize(value_dim);
            value.select({value_dim_old[0],0,0},{1,theSys.N,theSys.dim}).write(all_pos);
        }

        //Update velocity
        if(!file.exist("/particles/all/velocity/step"))
        {
            //std::cout << "creating velocity data" << std::endl;
            DataSet step = file.createDataSet<int>("/particles/all/velocity/step", part_t_space, props_time);
            step.write(theSys.time);
            DataSet time = file.createDataSet<double>("/particles/all/velocity/time", part_t_space, props_time);
            time.write(theSys.time*theSys.dt);
            DataSet value = file.createDataSet<double>("/particles/all/velocity/value", part_val_space, props_val);
            value.select({0,0,0},{1,theSys.N,theSys.dim}).write(all_vel);
        }
        else
        {
            //Update step
            DataSet step = file.getDataSet("/particles/all/velocity/step");
            std::vector<long unsigned int> step_dim = step.getDimensions();
            std::vector<long unsigned int> step_dim_old = step_dim;
            step_dim[0] += 1;
            step.resize(step_dim);
            step.select(step_dim_old,{1}).write(theSys.time);

            //Update time
            DataSet time = file.getDataSet("/particles/all/velocity/time");
            std::vector<long unsigned int> time_dim = time.getDimensions();
            std::vector<long unsigned int> time_dim_old = time_dim;
            time_dim[0] += 1;
            time.resize(time_dim);
            time.select(time_dim_old,{1}).write(theSys.time*theSys.dt);

            //Update velocity values
            DataSet value = file.getDataSet("/particles/all/velocity/value");
            std::vector<long unsigned int> value_dim = value.getDimensions();
            std::vector<long unsigned int> value_dim_old = value_dim;
            value_dim[0] += 1;
            value.resize(value_dim);
            value.select({value_dim_old[0],0,0},{1,theSys.N,theSys.dim}).write(all_vel);
        }

        //Update image
        if(!file.exist("/particles/all/image/step"))
        {
            DataSet step = file.createDataSet<int>("/particles/all/image/step", part_t_space, props_time);
            step.write(theSys.time);
            DataSet time = file.createDataSet<double>("/particles/all/image/time", part_t_space, props_time);
            time.write(theSys.time*theSys.dt);
            DataSet value = file.createDataSet<int>("/particles/all/image/value", part_val_space, props_val);
            value.select({0,0,0},{1,theSys.N,theSys.dim}).write(all_image);
        }
        else
        {
            //Update step
            DataSet step = file.getDataSet("/particles/all/image/step");
            std::vector<long unsigned int> step_dim = step.getDimensions();
            std::vector<long unsigned int> step_dim_old = step_dim;
            step_dim[0] += 1;
            step.resize(step_dim);
            step.select(step_dim_old,{1}).write(theSys.time);

            //Update time
            DataSet time = file.getDataSet("/particles/all/image/time");
            std::vector<long unsigned int> time_dim = time.getDimensions();
            std::vector<long unsigned int> time_dim_old = time_dim;
            time_dim[0] += 1;
            time.resize(time_dim);
            time.select(time_dim_old,{1}).write(theSys.time*theSys.dt);

            //Update image values
            DataSet value = file.getDataSet("/particles/all/image/value");
            std::vector<long unsigned int> value_dim = value.getDimensions();
            std::vector<long unsigned int> value_dim_old = value_dim;
            value_dim[0] += 1;
            value.resize(value_dim);
            value.select({value_dim_old[0],0,0},{1,theSys.N,theSys.dim}).write(all_image);
        }

        //Update energy
        if(!file.exist("/observables/potential_energy/step"))
        {
            //std::cout << "creating energy data" << std::endl;
            DataSet step = file.createDataSet<int>("/observables/potential_energy/step", part_t_space, props_time);
            step.write(theSys.time);
            DataSet time = file.createDataSet<double>("/observables/potential_energy/time", part_t_space, props_time);
            time.write(theSys.time*theSys.dt);
            DataSet value = file.createDataSet<double>("/observables/potential_energy/value", part_t_space, props_time);
            value.write(theSys.get_energy());
        }
        else
        {
            //Update step
            DataSet step = file.getDataSet("/observables/potential_energy/step");
            std::vector<long unsigned int> step_dim = step.getDimensions();
            std::vector<long unsigned int> step_dim_old = step_dim;
            step_dim[0] += 1;
            step.resize(step_dim);
            step.select(step_dim_old,{1}).write(theSys.time);

            //Update time
            DataSet time = file.getDataSet("/observables/potential_energy/time");
            std::vector<long unsigned int> time_dim = time.getDimensions();
            std::vector<long unsigned int> time_dim_old = time_dim;
            time_dim[0] += 1;
            time.resize(time_dim);
            time.select(time_dim_old,{1}).write(theSys.time*theSys.dt);

            //Update energy values
            DataSet value = file.getDataSet("/observables/potential_energy/value");
            std::vector<long unsigned int> value_dim = value.getDimensions();
            std::vector<long unsigned int> value_dim_old = value_dim;
            value_dim[0] += 1;
            value.resize(value_dim);
            value.select(value_dim_old,{1}).write(theSys.get_energy());
        }

        //Update noise fluctuations
        if(!file.exist("/observables/noise_fluctuations/step"))
        {
            //std::cout << "creating energy data" << std::endl;
            DataSet step = file.createDataSet<int>("/observables/noise_fluctuations/step", part_t_space, props_time);
            step.write(theSys.time);
            DataSet time = file.createDataSet<double>("/observables/noise_fluctuations/time", part_t_space, props_time);
            time.write(theSys.time*theSys.dt);
            DataSet value = file.createDataSet<double>("/observables/noise_fluctuations/value", part_vec_space, props_vec);
            value.select({0,0},{1,theSys.dim}).write(theSys.noise_flucs);
        }
        else
        {
            //Update step
            DataSet step = file.getDataSet("/observables/noise_fluctuations/step");
            std::vector<long unsigned int> step_dim = step.getDimensions();
            std::vector<long unsigned int> step_dim_old = step_dim;
            step_dim[0] += 1;
            step.resize(step_dim);
            step.select(step_dim_old,{1}).write(theSys.time);

            //Update time
            DataSet time = file.getDataSet("/observables/noise_fluctuations/time");
            std::vector<long unsigned int> time_dim = time.getDimensions();
            std::vector<long unsigned int> time_dim_old = time_dim;
            time_dim[0] += 1;
            time.resize(time_dim);
            time.select(time_dim_old,{1}).write(theSys.time*theSys.dt);

            //Update fluctuation values
            DataSet value = file.getDataSet("/observables/noise_fluctuations/value");
            std::vector<long unsigned int> value_dim = value.getDimensions();
            std::vector<long unsigned int> value_dim_old = value_dim;
            value_dim[0] += 1;
            value.resize(value_dim);
            value.select({value_dim_old[0],0},{1,theSys.dim}).write(theSys.noise_flucs);
        }

        if(theSys.can_break)
        {
            std::cout << "Warning: i/o for changing connectivity not yet supported. No changes to bond topology will be written." << std::endl;
        }

    }
    catch(Exception& err)
    {
        std::cerr << err.what() << std::endl;
    }
}


#include <experimental/filesystem>
#include "catch.hpp"
#include "../../src/Node.hpp"
#include "../../src/Spring.hpp"
#include "../../src/System.hpp"
#include "../../src/ParamDict.hpp"

namespace fs = std::experimental::filesystem;

TEST_CASE("System")
{
    ParamDict emptyParams;

    SECTION("Testing PBCs")
    {
        emptyParams.add_entry("N", "2");
        System smallSys(emptyParams);
        smallSys.network[0].set_pos(0.1, 0, 0);
        smallSys.network[1].set_pos(0.15, 0, 0);
        smallSys.Lx = 0.2;
        smallSys.apply_pbc();
        arma::vec disp = smallSys.get_disp_vec(smallSys.network[0], smallSys.network[1]);
        REQUIRE(disp(0)==Approx(-0.05));
        REQUIRE(disp(1)==0);
        REQUIRE(disp(2)==0);
        REQUIRE(smallSys.get_dist(smallSys.network[0],smallSys.network[1])==Approx(0.05));

        ParamDict myParams;
        myParams.read_params("test/test_fcc.conf");
        System system(myParams);
        for (int i=0; i<system.N; i++)
        {
            REQUIRE(system.network[i].get_num_springs()==12);
        }
        for (int i=0; i<system.N; i++)
        {
            system.network[i].incr_pos(0.2*system.Lx, 0.2*system.Ly, 0.2*system.Lz);
        }
        system.apply_pbc();
        for (int i=0; i<system.N; i++)
        {
            arma::vec pos = system.network[i].get_pos();
            REQUIRE(std::abs(pos(0))<=0.5*system.Lx);
            REQUIRE(std::abs(pos(1))<=0.5*system.Ly);
            REQUIRE(std::abs(pos(2))<=0.5*system.Lz);
        }
        
    }
    SECTION("Testing zero_com")
    {
        ParamDict myParams;
        myParams.add_entry("node_protocol", "fcc");
        myParams.add_entry("a", "1.41421356237");
        myParams.add_entry("rho", "1.41421356237");
        myParams.add_entry("is_p_x","1");
        myParams.add_entry("is_p_y","1");
        myParams.add_entry("is_p_z","1");

        System mySys(myParams);

        //Move all the nodes by some amount
        for(int i=0; i<mySys.N; i++)
        {
            mySys.network[i].pos(0) += 2.0;
            mySys.network[i].pos(1) += 2.0;
            mySys.network[i].pos(2) += 2.0;
        }
        mySys.zero_com();

        arma::vec com(3, arma::fill::zeros);
        for(int i=0; i<mySys.N; i++)
        {
            com += mySys.network[i].pos;
        }
        REQUIRE(abs(com(0))==Approx(0.0));
        REQUIRE(abs(com(1))==Approx(0.0));
        REQUIRE(abs(com(2))==Approx(0.0));

    }
    SECTION("Testing forces")
    {
        ParamDict myParams;
        myParams.add_entry("node_protocol", "fcc");
        myParams.add_entry("a", "1.41421356237");
        myParams.add_entry("rho", "1.41421356237");
        myParams.add_entry("is_p_x","1");
        myParams.add_entry("is_p_y","1");
        myParams.add_entry("is_p_z","1");

        System mySys(myParams);

        //Nodes should start off at their equilibrium positions
        REQUIRE(abs(mySys.get_energy())<1e-10);
        for(int i=0; i<mySys.N; i++)
        {
            arma::vec force = mySys.get_force(mySys.network[i]);
            for(int k=0; k<3; k++) REQUIRE(abs(force(k))<1e-10);
        }
    }
    SECTION("Testing image")
    {
        ParamDict myParams;
        myParams.add_entry("node_protocol", "fcc");
        myParams.add_entry("is_p_x","1");
        myParams.add_entry("is_p_y","1");
        myParams.add_entry("is_p_z","1");

        System mySys(myParams);

        //Make sure each node is in the 0th periodic box
        for(int i=0; i<mySys.N; i++)
        {
            for(int j=0; j<mySys.dim; j++)
            {
                REQUIRE(mySys.image[i][j]==0);
            }
        }
    }
}
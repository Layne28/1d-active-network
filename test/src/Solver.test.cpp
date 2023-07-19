#include "catch.hpp"
#include "../../src/System.hpp"
#include "../../src/ParamDict.hpp"
#include "../../src/Solver.hpp"
#include "../../src/CustomRandom.hpp"

TEST_CASE("Solver")
{
    ParamDict defaultParams;
    defaultParams.add_entry("node_protocol", "fcc");
    defaultParams.add_entry("a", "1.41421356237");
    defaultParams.add_entry("rho", "1.41421356237");
    defaultParams.add_entry("is_p_x","1");
    defaultParams.add_entry("is_p_y","1");
    defaultParams.add_entry("is_p_z","1");
    defaultParams.add_entry("kT","0.2");
    defaultParams.add_entry("gamma","10.0");

    SECTION("Testing RNG")
    {
        System system1(defaultParams);
        gsl_rng *myGen1 = CustomRandom::init_rng(1);        
        Solver mySolver1(defaultParams, myGen1);
        for(int i=0; i<100; i++) mySolver1.update(system1);

        //Same seed for RNG
        System system2(defaultParams);
        gsl_rng *myGen2 = CustomRandom::init_rng(1);        
        Solver mySolver2(defaultParams, myGen2);
        for(int i=0; i<100; i++) mySolver2.update(system2);

        //Different seed for RNG
        System system3(defaultParams);
        gsl_rng *myGen3 = CustomRandom::init_rng(2);        
        Solver mySolver3(defaultParams, myGen3);
        for(int i=0; i<100; i++) mySolver3.update(system3);

        for(int i=0; i<system1.N; i++)
        {
            for(int k=0; k<3; k++)
            {
                //Same seeds should produce same changes in node positions
                REQUIRE(system1.network[i].get_pos()(k)==system2.network[i].get_pos()(k));
                //Different seeds should produce different changes in node positions
                REQUIRE(system1.network[i].get_pos()(k)!=system3.network[i].get_pos()(k));
            }
        }
    }
}
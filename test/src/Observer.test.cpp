#include <experimental/filesystem>
#include "catch.hpp"
#include "../../src/Node.hpp"
#include "../../src/Spring.hpp"
#include "../../src/System.hpp"
#include "../../src/ParamDict.hpp"

namespace fs = std::experimental::filesystem;

TEST_CASE("Observer")
{

    ParamDict emptyParams;

    SECTION("Testing Observer in System")
    {
        std::string thedir = "myspecialtestdir";
        emptyParams.add_entry("output_dir",thedir);
        emptyParams.add_entry("network_freq","23");
        emptyParams.add_entry("thermo_freq","41");
        Observer myObs(emptyParams);
        System testSys(emptyParams);
        REQUIRE_THROWS(testSys.get_obs());
        REQUIRE(fs::exists(thedir));
        testSys.set_obs(myObs);     
        REQUIRE(testSys.get_obs().output_dir==thedir+"/");
        REQUIRE(testSys.get_obs().network_freq==23);
        REQUIRE(testSys.get_obs().thermo_freq==41);

        emptyParams.add_entry("network_freq","1");
        emptyParams.add_entry("thermo_freq","1");
        Observer realObs(emptyParams);
        testSys.set_obs(realObs);
        testSys.time = 0;
        testSys.get_obs().dump_nodes(testSys, "");
        testSys.get_obs().dump_springs(testSys, "");
        testSys.get_obs().dump_thermo(testSys, "");
        REQUIRE(fs::exists(thedir+"/trajectory/0.dump"));
        REQUIRE(fs::exists(thedir+"/topology/0.dump"));
        REQUIRE(fs::exists(thedir+"/thermo.dat"));
        fs::remove_all(thedir); //probably not safe
    }
}
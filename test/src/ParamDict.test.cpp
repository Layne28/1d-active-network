#include "catch.hpp"
#include "../../src/ParamDict.hpp"

TEST_CASE("ParamDict")
{
    ParamDict myParams;

    SECTION("Check for errors upon reading bad file")
    {
        REQUIRE_THROWS_AS(myParams.read_params("sample.txt"), std::runtime_error);
        REQUIRE_THROWS(myParams.read_params("nothere.conf"));
    }
    SECTION("Check that good file is parsed correctly")
    {
        myParams.read_params("test/test.conf"); //might have to change path to ensure this is found
        REQUIRE(myParams.get_value("param1")=="1.00");
        REQUIRE(myParams.get_value("param2")=="mystring");
        REQUIRE(myParams.get_size()==2);
    }
    
}   
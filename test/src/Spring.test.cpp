#include "catch.hpp"
#include "../../src/Node.hpp"
#include "../../src/Spring.hpp"
#include <iostream>

TEST_CASE("Spring")
{

    SECTION("Checking construction")
    {

        Node n1;
        Node n2;
        Node n3;

        Spring s12 = Spring::add_spring(n1, n2, 1.0, 2.0);
        Spring s13 = Spring::add_spring(n1, n3);
        Spring s23 = Spring::add_spring(n2, n3);

        REQUIRE(n1.get_num_springs()==2);
        REQUIRE(n2.get_num_springs()==2);
        REQUIRE(n3.get_num_springs()==2);

        REQUIRE(s12.get_stiffness()==1.0);
        REQUIRE(s23.get_stiffness()==0.0);
        REQUIRE(s12.get_rest_length()==2.0);
        REQUIRE(s23.get_rest_length()==0.0);
    }
    

    SECTION("Check set methods")
    {
        Node n1;
        Node n3;

        Spring s13 = Spring::add_spring(n1, n3);
        s13.set_rest_length(1.0);
        s13.set_stiffness(12.0);
        REQUIRE(s13.get_rest_length()==1.0);
        REQUIRE(s13.get_stiffness()==12.0);
    }

    SECTION("Check length and strain methods")
    {
        Node na;
        Node nb;
        nb.set_pos(0.0,0.0,0.0);
        nb.set_pos(0.0,0.0,2.2);
        Spring s = Spring::add_spring(na,nb);
        s.set_rest_length(2.0);
        REQUIRE(s.get_length()==2.2);
        REQUIRE(s.get_strain()==Approx(0.1));

        //Since Springs contain pointers to nodes,
        //they should update as nodes are updated.
        na.set_pos(0.234, -1.34, 5.6);
        nb.set_pos(0.0, -2.4, 0.001);
        REQUIRE(s.get_length()==Approx(5.703258454602946));
    }

    SECTION("Check for no duplicate springs")
    {
        Node n1;
        Node n2;

        Spring::add_spring(n1, n2);
        REQUIRE_THROWS_AS(Spring::add_spring(n1,n2), std::runtime_error);
    }

    SECTION("Check spring removal")
    {
        Node n2;
        Node n3;

        Spring::add_spring(n2, n3);
        Spring s = Spring::get_spring(n2, n3);
        REQUIRE(s.get_rest_length()==0.0);
        REQUIRE(s.get_stiffness()==0.0);
        REQUIRE(s.node1->is_equal(n2));
        REQUIRE(s.node2->is_equal(n3));

        Spring::remove_spring(n2, n3);
        REQUIRE(!Spring::is_spring(n2,n3));
        REQUIRE(!n2.has_connection(n3));
        REQUIRE(!n3.has_connection(n2));
        REQUIRE(s.node1->is_equal(n2)); //"Disembodied" Spring
    }
}
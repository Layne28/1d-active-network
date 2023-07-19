#include "catch.hpp"
#include "../../src/Node.hpp"
#include "../../src/Spring.hpp"

TEST_CASE("Node")
{
    Node n;

    REQUIRE(n.get_pos().size()==3);
    arma::vec pos = n.get_pos();
    for (int i=0; i<3; i++) REQUIRE(pos[i]==0.0);
    REQUIRE(n.get_num_springs()==0); //Nodes should have no springs when first constructed

    SECTION("Initialize node with position")
    {
        Node n2(1.0,2.0,3.0);
        pos = n2.get_pos();
        REQUIRE(pos[0]==1.0);
        REQUIRE(pos[1]==2.0);
        REQUIRE(pos[2]==3.0);
    }

    SECTION("Set node position")
    {
        n.set_pos(10.0,-2.5,0.334);
        REQUIRE(n.get_pos()[0]==10.0);
        REQUIRE(n.get_pos()[1]==-2.5);
        REQUIRE(n.get_pos()[2]==0.334);      
    }

    SECTION("Test id")
    {
        Node n2;
        REQUIRE(n.is_equal(n));
        REQUIRE(!n.is_equal(n2));
    }

    //This section implicitly tests the "Node::add_spring" function, which is called whenever a Spring is constructed.
    SECTION("Add springs")
    {
        Node n1;
        Node n2;

        SECTION("Check num springs")
        {
            Spring::add_spring(n1, n2);

            REQUIRE(n1.get_num_springs()==1);
            REQUIRE(n2.get_num_springs()==1);
        }
        SECTION("Check spring scope")
        {
            REQUIRE(n1.get_num_springs()==0);
            REQUIRE(n2.get_num_springs()==0);

            Spring::add_spring(n1, n2);
            REQUIRE(n1.get_num_springs()==1);
            REQUIRE(n2.get_num_springs()==1);
            REQUIRE(n1.springs[0].get_rest_length()==0.0);
            REQUIRE(n1.springs[0].get_stiffness()==0.0);
            REQUIRE(n2.springs[0].get_rest_length()==0.0);
            REQUIRE(n2.springs[0].get_stiffness()==0.0);
        }

        SECTION("Check connection")
        {
            REQUIRE(!n1.has_connection(n2));
            Spring::add_spring(n1, n2);
            REQUIRE(n1.has_connection(n2));
            REQUIRE(n2.has_connection(n1));
            REQUIRE(Spring::is_spring(n1, n2));
        }
    }
}
// A node object has a position (x, y, z) in 3D space and a vector of Springs connecting it to other Nodes.

#ifndef NODE_HPP
#define NODE_HPP

#include <array>
#include <vector>
#include <armadillo>
#include "Spring.hpp"

class Node
{
private:
    static int counter; //Warning -- not thread-safe!!
    int id;
    //arma::vec pos;

    friend class Spring;
    
public:
    std::vector<Spring> springs;
    arma::vec pos;
    arma::vec old_pos;
    arma::vec vel;

    //constructor
    Node(double x=0.0);

    //destructor
    ~Node();

    void set_pos(double x);
    void incr_pos(double dx);
    int get_num_springs();
    int get_id();
    arma::vec get_pos();

    bool has_connection(Node &n);
    bool is_equal(Node &n);

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, Node& n);
};

#endif

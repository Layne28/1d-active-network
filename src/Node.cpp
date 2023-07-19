#include <iostream>
#include "Node.hpp"

//TODO: write unit tests for object id/counter
int Node::counter = 0;

Node::Node(double x)
{
    pos.zeros(1);
    old_pos.zeros(1);
    vel.zeros(1);
    this->set_pos(x);
    id = Node::counter;
    Node::counter++;
}

Node::~Node()
{
    //TODO: Check whether springs' pointers need to be deleted here.
}

void Node::set_pos(double x)
{
    this->pos(0) = x;
}

void Node::incr_pos(double dx)
{
    this->pos(0) += dx;
}

int Node::get_num_springs()
{
    return this->springs.size();
}

int Node::get_id()
{
    return id;
}

arma::vec Node::get_pos()
{
    return pos;
}

bool Node::has_connection(Node &n)
{
    for (int i=0; i<(this->get_num_springs()); i++)
    {
        if (this->springs[i].node1->is_equal(n) || this->springs[i].node2->is_equal(n))
        {
            return true;
        }
    }
    return false;
}

bool Node::is_equal(Node &n)
{
    if (this->id==n.id)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::ostream& operator<<(std::ostream& os, Node& n) 
{
    arma::vec pos = n.get_pos();
    return os << "Node with id: " << n.get_id() << " and Pos: (" << pos(0) << "); No. of Springs: " << n.get_num_springs();
}

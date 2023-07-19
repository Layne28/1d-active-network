//A Spring object has a rest length, a stiffness, and pointers to the two Nodes it connects.

#ifndef SPRING_HPP
#define SPRING_HPP

class Node;

class Spring
{
private:
    friend class Node;

public:
    double K; //stiffness
    double l0; //rest length
    Node* node1;
    Node* node2;

    //constructor (make this private in the future?)
    Spring(Node& n1, Node& n2, double K_1=0.0, double l0_1=0.0);

    //destructor
    ~Spring();

    void set_stiffness(double K_1);
    void set_rest_length(double l0_1);
    double get_stiffness();
    double get_rest_length();

    //Move these two to System
    double get_length(); //need to make this respect pbc
    double get_strain(); //need to make this respect pbc

    static bool is_spring(Node &n1, Node &n2);
    static Spring add_spring(Node &n1, Node &n2, double K_1=0.0, double l0_1=0.0);
    static void remove_spring(Node &n1, Node &n2);
    static Spring get_spring(Node &n1, Node &n2);

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, Spring& s);
};

#endif
#include "Node.hpp"
#include <iostream>

Node::Node(void) {
    left = right = primAncestor = hybAncestor = NULL;
    index = 0;
    name = "";
    branchLength = 0;
    time = -1;
}

void Node::setTime(double t) {
    time = t;
}

double Node::getTime() {
    return time;
}

void Node::print(void) {
    std::cout << "Node " << index << "(" << this << ")" << std::endl;
    std::cout << " Lft: " << left << std::endl;
    std::cout << " Rht: " << right << std::endl;
    std::cout << " PrimAnc: " << primAncestor << std::endl;
    std::cout << " HybAnc: " << hybAncestor << std::endl;
    std::cout << " Name: \"" << name << "\"" << std::endl;
    std::cout << " Brlen: " << branchLength << std::endl;
}
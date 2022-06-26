#include "Node.hpp"
#include <iostream>

Node::Node(void) {
    left = right = majorAncestor = minorAncestor = NULL;
    index = 0;
    name = "";
    majorBranchLength = 0;
    minorBranchLength = 0;
    gamma = gammaLft = gammaRht = 0;
    bootSupport = -1;
    time = -1;
}

void Node::printInfo(void) {
    if(time != -1)
        std::cout << "\tTime: " << time << std::endl;

    if(getLft() != NULL) {
        std::cout << "\tLeft: " << left->name << " (" << left << ")";
        if(gammaLft != 0)
            std::cout << ", gamma=" << gammaLft << std::endl;
        else
            std::cout << std::endl;
    }


    if(getRht() != NULL) {
        std::cout << "\tRight: " << right->name << " (" << right << ")";
        if(gammaRht != 0)
            std::cout << ", gamma=" << gammaRht << std::endl;
        else
            std::cout << std::endl;
    }
        
    if(getMajorAnc() != NULL)
        std::cout << "\tMajor anc: " << majorAncestor->name << std::endl;
    if(getMinorAnc() != NULL)
        std::cout << "\tMinor anc: " << minorAncestor->name << std::endl;
    
    if(bootSupport != -1)
        std::cout << "\tBoot supp: " << bootSupport << std::endl;
}
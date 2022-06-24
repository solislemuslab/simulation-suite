#include "Node.hpp"
#include <iostream>

Node::Node(void) {
    left = right = majorAncestor = minorAncestor = NULL;
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

void Node::printInfo(void) {
        if(getLft() != NULL)
            std::cout << "\tLeft: " << getLft()->getName();
        if(gammaLft != 0)
            std::cout << ", gamma=" << gammaLft << std::endl << std::flush;
        else
            std::cout << std::endl << std::flush;


        if(getRht() != NULL)
            std::cout << "\tRight: " << getRht()->getName();
        if(gammaRht != 0)
            std::cout << ", gamma=" << gammaRht << std::endl << std::flush;
        else
            std::cout << std::endl << std::flush;
            
        if(getMajorAnc() != NULL)
            std::cout << "\tMajor anc: " << getMajorAnc()->getName() << std::endl << std::flush;
        if(getMinorAnc() != NULL)
            std::cout << "\tMinor anc: " << getMinorAnc()->getName() << std::endl << std::flush;
        
        std::cout << "\tBoot supp: " << bootSupport << std::endl << std::flush;
}
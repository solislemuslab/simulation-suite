#include "Node.hpp"
#include <iostream>

Node::Node(void) {
    left = right = majorAncestor = minorAncestor = NULL;
    index = 0;
    name = "";
    majorBranchLength = 0;
    minorBranchLength = 0;
    time = -1;
}

void Node::printInfo(void) {
    std::cout << "\tTime: " << time << std::endl;

    if(getLft() != NULL)
        std::cout << "\tLeft: " << getLft()->getName();
    if(gammaLft != 0)
        std::cout << ", gamma=" << gammaLft << std::endl;
    else
        std::cout << std::endl;


    if(getRht() != NULL)
        std::cout << "\tRight: " << getRht()->getName();
    if(gammaRht != 0)
        std::cout << ", gamma=" << gammaRht << std::endl;
    else
        std::cout << std::endl;
        
    if(getMajorAnc() != NULL)
        std::cout << "\tMajor anc: " << getMajorAnc()->getName() << std::endl;
    if(getMinorAnc() != NULL)
        std::cout << "\tMinor anc: " << getMinorAnc()->getName() << std::endl;
    
    std::cout << "\tBoot supp: " << bootSupport << std::endl;
}
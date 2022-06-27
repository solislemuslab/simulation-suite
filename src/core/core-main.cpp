#include "Network.hpp"
#include "MSEvents.hpp"
#include <iostream>

int main(int narg, char *argv[]) {
    // Check for Newick string equivalence
    std::string n1 = std::string("((A:0.6,(B:0.5,#H1:0.1::0.6):0.1):0.2,((C:0.4)#H1:0.3::0.4,D:0.7):0.1);");
    std::string n2 = std::string("((A:0.6,(B:0.5,(C:0.4)#H1:0.1::0.6):0.1):0.2,(#H1:0.3::0.4,D:0.7):0.1);");

    std::cout << "Building net1" << std::endl;
    Network net1(n1, "newick");
    std::cout << "Building net2" << std::endl;
    Network net2(n2, "newick");

    if(isomorphicNewick(n1, n2))
        std::cout << "Equivalent!" << std::endl;
    else
        std::cout << "Not equivalent :(" << std::endl;

    return 0;
}
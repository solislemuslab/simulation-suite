#include "Network.hpp"
#include "MSEvents.hpp"
#include <iostream>

int main(int narg, char *argv[]) {
    // Building a network from a Newick string
    std::string myNewick = "(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;";
    Network *myNetwork = new Network(myNewick, "newick");
    std::cout << myNetwork->getNewickRepresentation();

    return 0;
}
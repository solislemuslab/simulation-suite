#include "Network.hpp"
#include "MSEvents.hpp"
#include <iostream>
#include <algorithm>

int main(int narg, char *argv[]) {
    // Building a network from a Newick string
    std::string myNewick = "(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;";
    Network *myNetwork = new Network(myNewick, "newick");

    // Retrieving an ms command that corresponds to an equivalent network
    // WARNING: at the moment this function essentially destroys the underlying network,
    //          so myNetwork cannot be reused after this. (will be fixed later)
    std::string msCmd = myNetwork->getMSString();
    std::cout << "\n\nCommand in ms format: " << msCmd << std::endl;

    // Reading a network from a sequence of ms joins and splits
    Network *msNet = new Network(msCmd, "ms");

    // A more explicit ms example
    Network net("-ej 1.0 1 2 -es 1.5 3 0.25 -es 2.0 4 0.33 -ej 2.1 2 4 -ej 2.4 3 4 -ej 2.6 5 4", "ms");

    return 0;
}
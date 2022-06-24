#include "Network.hpp"

int main(int narg, char *argv[]) {
    // Read in network from Newick string
    std::string myNewick = "(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;";
    Network myNetwork(myNewick);

    // Get the ms parameters
    myNetwork.toms();

    return 0;
}
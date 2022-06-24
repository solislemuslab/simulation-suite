#include "Network.hpp"

int main(int narg, char *argv[]) {
    // Read in network from Newick string
    std::string myNewick = "(((A,(B)#H1:::0.9)e,(C,#H1:::0.1)f)g,D)h;";
    Network myNetwork(myNewick);

    // Get the ms parameters
    myNetwork.toms();

    return 0;
}
#include "Network.hpp"

int main(int narg, char *argv[]) {
    // Read in network from Newick string
    std::string myNewick = "((1,((2,(3,(4)Y#H1)g)e,(((Y#H1,5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;";
    Network myNetwork(myNewick);

    // Get the ms parameters
    myNetwork.toms();

    return 0;
}
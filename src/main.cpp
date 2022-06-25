#include "Network.hpp"
#include "MSEvents.hpp"
#include <iostream>
#include <algorithm>

int main(int narg, char *argv[]) {
    // Read in network from Newick string
    std::string myNewick = "(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;";
    Network myNetwork(myNewick);

    // Get the ms parameters
    std::vector<MSEvent*> myEvents = myNetwork.toms();
    std::sort(myEvents.begin(), myEvents.end(), [](MSEvent *a, MSEvent*b) {
        return a->getTime() < b->getTime();
    });

    std::cout << "ms arguments corresponding to \"" << myNewick << "\", ordered by time:" << std::endl;
    for(MSEvent *e : myEvents) {
        std::cout << "\t";
        e->print();
    }

    return 0;
}
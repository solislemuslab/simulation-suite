#include "Network.hpp"
#include "MSEvents.hpp"
#include <iostream>

int main(int narg, char *argv[]) {
    // Check for Newick string equivalence
    std::string n1 = std::string("((A:0.6,(B:0.5,#H1:0.1::0.6):0.1):0.2,((C:0.4)#H1:0.3::0.4,D:0.7):0.1);");
    std::string n2 = std::string("((A:0.6,(B:0.5,(C:0.4)#H1:0.1::0.6):0.1):0.2,(#H1:0.3::0.4,D:0.7):0.1);");
    std::string n3 = std::string("((1:0.1,((2:0.2,(3:0.3,(4:0.4)Y#H1:0.5)g:0.6)e:0.7,(((Y#H1:0.8,5:0.9)h:1.0,6:1.1)f:1.2)X#H2:1.3)c:1.4)a:1.5,((X#H2:1.6,7:1.7)d:1.8,8:1.9)b:2.0)r;");

    std::cout << "Building net1" << std::endl;
    Network net1(n1, "newick");
    std::cout << "Building net2" << std::endl;
    Network net2(n2, "newick");

    if(isomorphicNewick(n1, n2))
        std::cout << "Equivalent!" << std::endl;
    else
        std::cout << "Not equivalent :(" << std::endl;

    Network net3(n3, "newick");
    std::cout << "\nBase Newick representation:\n";
    std::cout << net3.getNewickRepresentation();

    // std::vector<std::string> newicks = net3.getRandomNewickRepresentations(5);
    // std::cout << "\n\nRandom equivalent Newick representations:\n";
    // for(std::string str : newicks)
    //     std::cout << str << std::endl;

    return 0;
}
#include "SimSuite.hpp"
#include "core/Network.hpp"
#include <iostream>

int main(int narg, char **argv) {
    // -----------------------------------------------------------------------
    // - Example content; should be moved to a different file at some point. -
    // -----------------------------------------------------------------------

    // Convert extended Newick to MS
    std::string msString = SimSuite::newickToMS("(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;");
    std::cout << msString << std::endl;
    
    // ----------------------------
    // - Testing section; ignore. -
    // ----------------------------

    std::string newick = std::string("((1:0.1,((2:0.2,(3:0.3,(4:0.4)Y#H1:0.5)g:0.6)e:0.7,(((Y#H1:0.8,5:0.9)h:1.0,6:1.1)f:1.2)X#H2:1.3)c:1.4)a:1.5,((X#H2:1.6,7:1.7)d:1.8,8:1.9)b:2.0)r;");
    Network net(newick, "newick");

    std::vector<std::string> randomNewicks = net.getRandomNewickRepresentations(10);

    // std::cout << randomNewicks[0] << std::endl;
    for(unsigned int i=0; i < randomNewicks.size()-1; i++) {
        // std::cout << "\t\t" << (isomorphicNewick(randomNewicks[i], randomNewicks[i+1]) ? "EQUIVALENT" : "NOT EQUIVALENT") << std::endl << randomNewicks[1] << std::endl;
        std::cout << randomNewicks[i] << std::endl << randomNewicks[i+1] << std::endl;
        std::cout << "\t" << (isomorphicNewick(randomNewicks[i], randomNewicks[i+1]) ? "EQUIVALENT" : "NOT EQUIVALENT") << std::endl;
    }
}

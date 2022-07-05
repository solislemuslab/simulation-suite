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

    // Some tests by GAB
    // Not sure if this should fail because it lacks the root name h.
    // Also, I've only seen it in use in hybrid-Lambda... do we need it at all?
    std::cout << "GAB: No brlens, no h for root" << std::endl;
    std::string netNoBrlensNoh = SimSuite::newickToMS("((A,B),C);");
    std::cout << netNoBrlensNoh << std::endl;

    // This should probably fail because it is not ultrametric and coalescent trees should be
    std::cout << "GAB: No brlens" << std::endl;
    std::string netNoBrlens = SimSuite::newickToMS("((A,B),C)h;");
    std::cout << netNoBrlens << std::endl;    

    // Positive simplest test
    std::cout << "GAB: Tinyest example possible, two branches of the same length" << std::endl;
    std::string tinyNet = SimSuite::newickToMS("(A:0.5,B:0.5)h;");
    std::cout << tinyNet << std::endl;    

    // A bit larger tree.
    std::cout << "GAB: Three-leaf example, ultrametric" << std::endl;
    std::string threeTaxTree = SimSuite::newickToMS("((A:0.5,B:0.5):0.5,C:1.0)h;");
    std::cout << threeTaxTree << std::endl;    
    
    // ----------------------------
    // - Testing section; ignore. -
    // ----------------------------

    //    std::string newick = std::string("((1:0.1,((2:0.2,(3:0.3,(4:0.4)Y#H1:0.5)g:0.6)e:0.7,(((Y#H1:0.8,5:0.9)h:1.0,6:1.1)f:1.2)X#H2:1.3)c:1.4)a:1.5,((X#H2:1.6,7:1.7)d:1.8,8:1.9)b:2.0)r;");
    //    Network net(newick, "newick");
    //
       //    std::vector<std::string> randomNewicks = net.getRandomNewickRepresentations(10);
    //
       //    // std::cout << randomNewicks[0] << std::endl;
       //    for(unsigned int i=0; i < randomNewicks.size()-1; i++) {
      //        // std::cout << "\t\t" << (isomorphicNewick(randomNewicks[i], randomNewicks[i+1]) ? "EQUIVALENT" : "NOT EQUIVALENT") << std::endl << randomNewicks[1] << std::endl;
      //        std::cout << randomNewicks[i] << std::endl << randomNewicks[i+1] << std::endl;
      //        std::cout << "\t" << (isomorphicNewick(randomNewicks[i], randomNewicks[i+1]) ? "EQUIVALENT" : "NOT EQUIVALENT") << std::endl;
      //    }
}

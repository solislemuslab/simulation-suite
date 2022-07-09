#include "SimSuite.hpp"
#include "core/Network.hpp"
#include <iostream>

int main(int narg, char **argv) {
    // -----------------------------------------------------------------------
    // - Example content; should be moved to a different file at some point. -
    // -----------------------------------------------------------------------

    // Convert extended Newick to MS
    std::string msString = SimSuite::newickToMS("(15:11.0,(1:10.0,((14:8.0,(((7:2.8,(10:1.6,(9:0.4,8:0.4):1.2):1.2):0.8,(11:2.8,(13:0.4,12:0.4):2.4):0.8):3.4,#H1:0.4::0.3):1.0):1.2,(((2:0.4,3:0.4):5.2,((4:3.6,5:3.6):1.2,6:4.8):0.8):1.0)#H1:2.6::0.7):0.8):1.0);");
    std::cout << msString << std::endl;
}

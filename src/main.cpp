#include "SimSuite.hpp"
#include <iostream>

int main(int narg, char **argv) {
    // Convert extended Newick to MS
    std::string msString = SimSuite::newickToMS("(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;");
    std::cout << msString << std::endl;
}
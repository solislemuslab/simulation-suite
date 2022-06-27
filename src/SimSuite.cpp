#include "core/Network.hpp"
#include "core/MSEvents.hpp"


std::string newickToMS(std::string newickStr) {
    Network *net = new Network(newickStr, "newick");
    std::string msString = net->getMSString();
    delete net;
    return msString;
}
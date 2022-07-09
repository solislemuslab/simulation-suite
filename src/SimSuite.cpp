#include "core/Network.hpp"
#include "core/MSEvents.hpp"

#include <iostream>
#include <string>
#include <iomanip>

namespace SimSuite {
    std::string newickToMS(std::string newickStr) {
        Network net(newickStr, "newick");
        std::string msString = net.getMSString();
        return msString;
    }

    std::string msToNewick(std::string msStr) {
        Network net(msStr, "ms");
        std::string newick = net.getNewickRepresentation();
        return newick;
    }
}

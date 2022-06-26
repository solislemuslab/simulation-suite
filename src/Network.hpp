#ifndef DEBUG_NEWICK_PARSER
#define DEBUG_NEWICK_PARSER
#endif

#ifndef Network_hpp
#define Network_hpp

#include "MSEvents.hpp"

#include <vector>
#include <set>
#include <string>
#include <sstream>
class Node;

class Network {
    public:
        void listNodes(void);
        ~Network(void);
        std::string getNewickRepresentation(void);
        int totalNodes(void) { return nodes.size(); }
        Network(std::string newickStr);
        Network(std::vector<MSEvent*> events);
        std::vector<MSEvent*> toms(void);

    protected:
        Network(void) { }
        Node* root;
        std::vector<Node*> nodes;
        void writeNetwork(Node *p, std::stringstream& ss);
        bool isHybridName(std::string);
        int hybridNameIndex(std::string val, std::vector<std::string> list);
        int activeNodesIdx(Node *p, std::vector<Node*>);
        void easyDbg(void) { }
        int getLength(std::string);
        bool blankName(Node *p);

    private:
        std::vector<std::string> parseNewick(std::string ns);
        void patchNetwork(void);
        void setTimes(void);
        void setTimeRecur(Node *p);
        double totalTime;
        void postmsPatchAndRename(void);
        std::string getInternalName(int);
        std::string getLeafName(int);
        int postmsPatchAndRenameRecur(Node*, int);
};

#endif
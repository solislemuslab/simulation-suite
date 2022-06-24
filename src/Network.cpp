#include "Node.hpp"
#include "Network.hpp"

#include <iostream>
#include <iomanip>
#include <set>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <typeinfo>

using namespace std;
std::string BLANK_NAME = "__$!&*#%*__";

int main(int argc, char *argv[]) {
    // set the parameters of the simulation
    int numSites = 100;
    std::string myNewick = "(((A,(B)#H1),(C,#H1)),D);";
    
    // instantiate the Network and print the ms -ej and -es parameters
    Network n(myNewick);
    n.toms();

    return 0;
}

Network::Network(std::string newickStr) {
    hybrids = 0;
    std::vector<std::string> tokens = parseNewickFirstPass(newickStr);

    // build up the Network from the parsed Newick string
    bool readingBranchLength = false;
    bool namingInternalNode = false;
    Node* p = NULL;
    for (int i=0; i<tokens.size(); i++) {
        std::string token = tokens[i];
        //std::cout << token << std::endl;
        if (token == "(") {
            namingInternalNode = false;
            // new node
            Node* newNode = new Node;
            nodes.push_back(newNode);
            if (p == NULL)
                root = newNode;
            else {
                newNode->setPrimAnc(p);
                if (p->getLft() == NULL)
                    p->setLft(newNode);
                else
                    p->setRht(newNode);
            }
            
            p = newNode;
        } else if (token == ")" || token == ",") {
            namingInternalNode = false;
            if(token == ")")
                namingInternalNode = true;

            // move down one node
            if (p->getPrimAnc() == NULL) {
                std::cout << "Error: We cannot find an expected ancestor at i=" << i << "; p == root gives: " << (p == root) << std::endl;
                exit(1);
            }
            p = p->getPrimAnc();
        } else if (token == ":") {
            namingInternalNode = false;
            // begin reading a branch length
            readingBranchLength = true;
        } else if (token == ";") {
            namingInternalNode = false;
            // finished!
            if (p != root) {
                std::cout << "Error: We expect to finish at the root node" << std::endl;
                exit(1);
            }
        } else {
            // we have a taxon name or a branch length
            if (readingBranchLength == true) {
                double x = std::stod(token);
                p->setBranchLength(x);
                readingBranchLength = false;
            } else if(namingInternalNode) {
                cout << "\tNaming internal node at i=" << i <<  "." << std::endl;
                p->setName(token);
            } else {
                Node* newNode = new Node;
                nodes.push_back(newNode);
                newNode->setPrimAnc(p);
                
                if (p->getLft() == NULL)
                    p->setLft(newNode);
                else
                    p->setRht(newNode);
                
                // Check if this is a hybrid node to add to the hybrids ticker
                char ch;
                int j = 0;
                do {
                    ch = token[j];
                    if(ch == '#') {
                        hybrids++;
                        break;
                    }
                    j += 1;
                } while(ch != '\0');

                newNode->setName(token);
                p = newNode;
            }
            namingInternalNode = false;
        }
    }

    // index the nodes
    int ndeIdx = 0;
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->getLft() == NULL)
            nodes[i]->setIndex(ndeIdx++);
    }
    
    for (int i=0; i<nodes.size(); i++) {
        if (nodes[i]->getLft() != NULL)
            nodes[i]->setIndex(ndeIdx++);
    }

    patchNewick();

    cout << std::endl << std::endl;
    for(int i = 0; i < nodes.size(); i++) {
        cout << "Node " << i << ": " << nodes[i]->getName() << std::endl << std::flush;
        if(nodes[i]->getLft() != NULL)
            cout << "\tLeft: " << nodes[i]->getLft()->getName() << std::endl << std::flush;
        if(nodes[i]->getRht() != NULL)
            cout << "\tRight: " << nodes[i]->getRht()->getName() << std::endl << std::flush;
            
        if(nodes[i]->getPrimAnc() != NULL)
            cout << "\tPrim anc: " << nodes[i]->getPrimAnc()->getName() << std::endl << std::flush;
        if(nodes[i]->getHybAnc() != NULL)
            cout << "\tHyb anc: " << nodes[i]->getHybAnc()->getName() << std::endl << std::flush;
    }
    cout << "Hybrids: " << hybrids << std::endl;
    easyDbg();
}

std::vector<std::string> Network::parseNewickFirstPass(std::string ns) {
    std::vector<std::string> tks;
    for(int i = 0; i < ns.size(); i++) {
        char c = ns[i];
        if(c == '(' || c == ')' || c == ',' || c == ':' || c == ';') {
            std::string tempStr;
            tempStr = c;
            tks.push_back(tempStr);
        } else {
            int j = i;
            std::string tempStr = "";
            while(!(c == '(' || c == ')' || c == ',' || c == ':' || c == ';')) {
                tempStr += c;
                j++;
                c = ns[j];
            }
            i = j - 1;
            tks.push_back(tempStr);
        }
    }

    #if defined(DEBUG_NEWICK_PARSER)
    std::cout << "The Newick string, broken into parts: " << std::endl;
    for(int i = 0; i < tks.size(); i++) {
        std::cout << " tks[" << i << "] = \"" << tks[i] << "\"" << std::endl;
    }
    #endif
    return tks;
}

// Patches up the hybrid edges on the network 
void Network::patchNewick() {
    cout << std::endl << std::endl;
    std::vector<std::string> hybridNames;
    std::vector<int> nodeIndices;
    std::vector<int> removeMe;

    int count = 0;
    for(int i = 0; i < nodes.size(); i++) {
        Node *p = nodes[i];
        // Check if this node is a hybrid
        if(isHybridName(p->getName())) {
            // It's a hybrid! Is it already in the list?
            int idx = hybridNameIndex(p->getName(), hybridNames);
            if(idx == -1) {
                hybridNames.push_back(p->getName());
                nodeIndices.push_back(i);
                count++;
            } else {
                // This is the copy of the node that we are killing
                Node *dead = nodes[nodeIndices[idx]];

                // Make the primAnc of the dead node's child into p
                // The dead node's child should always be a left child, but better safe than sorry...
                Node *deadChild = dead->getLft() != NULL ? dead->getLft() : dead->getRht();
                deadChild->setPrimAnc(p);

                // Now make the child of the dead node p's child
                // Make sure we aren't overriding a node if p already has a child
                if(p->getLft() != NULL)
                    p->setRht(deadChild);
                else
                    p->setLft(deadChild);

                // Move the edge running from dead node's anc->dead node to be dead node's anc->p
                if(dead->getPrimAnc()->getLft() == dead)
                    dead->getPrimAnc()->setLft(p);
                else
                    dead->getPrimAnc()->setRht(p);

                // Make p's hybrid ancestor the dead node's ancestor
                p->setHybAnc(dead->getPrimAnc());

                // Set the dead node to null
                removeMe.push_back(nodeIndices[idx]);
            }
        }
    }

    // Remove the dead nodes from our node list
    // Reverse sort the list so that we don't mess up indices as we go.
    std::sort(removeMe.begin(), removeMe.end(),
        [](int a, int b) {
            return(a > b);
        }
    );
    for(int idx : removeMe) {
        nodes.erase(std::next(nodes.begin(), idx));
    }

    cout << std::endl << std::endl;
}

int Network::hybridNameIndex(std::string val, std::vector<std::string> list) {
    // list has length <hybrids>
    for(int i = 0; i < list.size(); i++) {
        if(val.compare(list[i]) == 0)
            return i;
    }
    return -1;
}

int Network::getLength(std::string val) { return val.length(); }

bool Network::isHybridName(std::string val) {
    int len = getLength(val);
    for(int i = 0; i < len; i++) {
        if(val[i] == '#')
            return true;
    }
    return false;
}

std::string Network::getNewickRepresentation(void) {
    std::stringstream ss;
    
    if (root->getLft() != NULL && root->getRht() != NULL)
        writeNetwork(root, ss);
    else
        writeNetwork(root->getLft(), ss);
    
    std::string newick = ss.str();
    return newick;
}

void Network::writeNetwork(Node* p, std::stringstream& ss) {
    if (p != NULL) {
        if (p->getLft() == NULL) {
            // Use this line if you do not want to include branch lengths
            ss << p->getName();
            // Use this line if you want to include branch lengths
            // ss << p->getName() << ":" << std::fixed << std::setprecision(5) << p->getBranchLength();
        } else {
            ss << "(";
            writeNetwork (p->getLft(), ss);
            if(p->getRht() != NULL) {
                ss << ",";
                writeNetwork (p->getRht(), ss);
            }
            if (p->getPrimAnc() == NULL)
                ss << ")";
            else
                // Use this line if you do not want to include branch lengths
                ss << ")";
                // Use this line if you want to include branch lengths
                // ss << "):" << std::fixed << std::setprecision(5) << p->getBranchLength();
        }
    }
}

Network::~Network(void) {
    printf("Calling the Network destructor with %lu nodes.\n", nodes.size());
    for(int i = 0; i < nodes.size(); i++)
        delete nodes[i];
}

void Network::listNodes(void) {
    for(int i = 0; i < nodes.size(); i++)
        nodes[i]->print();
}

void Network::toms(void) {
    // First, gather & name all leaves while setting the names of non-leaves to null
    std::vector<Node*> activeNodes;
    int popnCounter = 1;
    for(int i = 0; i < nodes.size(); i++) {
        if(nodes[i] != NULL) {
            if(nodes[i]->getLft() == NULL && nodes[i]->getRht() == NULL) {
                activeNodes.push_back(nodes[i]);
                cout << "Setting " << nodes[i]->getName() << " to next val" << std::endl;
                nodes[i]->setName(popnCounter++);
            } else {
                cout << "Setting " << nodes[i]->getName() << " to BLANK_NAME" << std::endl;
                nodes[i]->setName(BLANK_NAME);
            }
        }
    }

    /* Go back in time from each of the leaves. Don't worry about times for now. Right now, if the node we're looking at
     * has two ancestors, this is a hybridization event and we need to name both of the ancestors.
     * If the node has one ancestor, this is a coalescence event. If this node's sister is not named, do not perform this
     * coalsecence yet. Wait for the node to be named.
     */
    // TODO: Definitely need to go through every event chronologically. Otherwise directionality is lost.
    // In order to do it chronologically: idk. have to think about this more.
    // Also: maybe we don't actually have to do it chronologically, and the direction can be preserved b/c we specifically
    //       assign hybAnc() in a way that can be used to maintain directionality.
    while(activeNodes.size() > 0) {
        // Run a for loop that tries to resolve every node in activeNodes
        std::vector<int> removeMe;
        std::vector<Node*> addMe;
        int i;
        for(i = 0; i < activeNodes.size(); i++) {
            Node *p = activeNodes[i];
            if(p == NULL) {
                cout << "ACTIVE NODE IS NULL" << std::endl << std::flush;
                exit(-1);
            }

            Node *pAnc = p->getPrimAnc();
            Node *hAnc = p->getHybAnc();
            // a. if only one ancestor (this will only ever take place when there is a primAnc and NOT a hybAnc)
            if(pAnc != NULL && hAnc == NULL) {
                // I. if the anc is NOT named: take the anc
                if(blankName(pAnc)) {
                    pAnc->setName(p->getName());

                    // Remove p and add pAnc
                    removeMe.push_back(i);
                    addMe.push_back(pAnc);
                }
                // II. if anc named: coalesce *INTO* the ancestor
                else {
                    cout << "-ej t " << p->getName() << " " << pAnc->getName() << std::endl;

                    // Remove p; pAnc is already named, so if it still has more potential things to do, it should be
                    // in activeNodes already. So we don't need to touch it.
                    removeMe.push_back(i);
                }
            }
            // b. if we have two ancestors
            else if(pAnc != NULL && hAnc != NULL) {
                // I. both ancs unnamed: split, give primAnc our name and give hybAnc a new name
                if(blankName(pAnc) && blankName(hAnc)) {
                    // Split
                    cout << "-es t " << p->getName() << " gamma" << std::endl;

                    // Give primAnc our name and add primAnc to activeNodes
                    pAnc->setName(p->getName());
                    addMe.push_back(pAnc);

                    // Give hybAnc its new name and add it to activeNodes
                    hAnc->setName(popnCounter++);
                    addMe.push_back(hAnc);

                    // Remove p from activeNodes
                    removeMe.push_back(i);
                }
                // II. if only one anc is named: split, give the blank anc a new name, and then join with the named anc
                else if(blankName(pAnc) != blankName(hAnc)) {
                    // Split
                    cout << "-es t " << p->getName() << " gamma" << std::endl;
                    
                    Node *namedAnc = blankName(pAnc) ? hAnc : pAnc;
                    Node *unnamedAnc = blankName(pAnc) ? pAnc : hAnc;

                    // Give the blank anc a new name and add it to the queue
                    unnamedAnc->setName(popnCounter++);
                    addMe.push_back(unnamedAnc);

                    // Join
                    cout << "-ej t " << p->getName() << " " << namedAnc->getName() << std::endl;

                    // Remove p
                    removeMe.push_back(i);
                }
                // III. if both ancs are named, split, then join both the new lineages with the existing ones
                else {
                    // Split
                    cout << "-es t " << p->getName() << " gamma" << std::endl;

                    // Join left
                    cout << "-ej t " << p->getName() << " " << pAnc->getName() << std::endl;

                    // Join right
                    cout << "-ej t " << popnCounter++ << " " << hAnc->getName() << std::endl;

                    // Remove p
                    removeMe.push_back(i);
                }
            }
            // c. if we don't have ancestors, we are root, so just remove us
            else {
                removeMe.push_back(i);
            }
        }


        // Reverse sort the list so that we don't mess up indices as we go.
        std::sort(removeMe.begin(), removeMe.end(),
            [](int a, int b) {
                return(a > b);
            }
        );
        for(int idx : removeMe) {
            activeNodes.erase(std::next(activeNodes.begin(), idx));
        }
        // Add the addMe nodes to activeNodes
        for(Node *node : addMe) {
            activeNodes.push_back(node);
        }
    }
}

bool Network::blankName(Node *p) {
    return p->getName().compare(BLANK_NAME) == 0;
}

int Network::activeNodesIdx(Node *p, std::vector<Node*> list) {
    for(int i = 0; i < list.size(); i++) {
        if(list[i] == p)
            return i;
    }
    return -1;
}
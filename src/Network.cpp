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

Network::Network(std::string newickStr) {
    hybrids = 0;
    std::vector<std::string> tokens = parseNewick(newickStr);

    // build up the Network from the parsed Newick string
    bool readingBranchLength = false;
    bool readingBootSupport = false;
    bool readingGamma = false;
    bool namingInternalNode = false;
    Node* p = NULL;
    for (int i=0; i<tokens.size(); i++) {
        std::string token = tokens[i];
        //std::cout << token << std::endl;
        if (token == "(") {
            readingBranchLength = false;
            readingBootSupport = false;
            readingGamma = false;
            namingInternalNode = false;
            // new node
            Node* newNode = new Node;
            nodes.push_back(newNode);
            if (p == NULL) {
                root = newNode;
            } else {
                newNode->setMajorAnc(p);
                if (p->getLft() == NULL)
                    p->setLft(newNode);
                else
                    p->setRht(newNode);
            }
            
            p = newNode;
        } else if (token == ")" || token == ",") {
            readingBranchLength = false;
            readingBootSupport = false;
            readingGamma = false;
            namingInternalNode = false;
            if(token == ")")
                namingInternalNode = true;

            // move down one node
            if (p->getMajorAnc() == NULL) {
                std::cout << "Error: We cannot find an expected ancestor at i=" << i << "; p == root gives: " << (p == root) << std::endl;
                exit(1);
            }
            p = p->getMajorAnc();
        } else if (token == ":") {
            namingInternalNode = false;

            // if the field for branch length, boot support, or inheritance probability are left blank, then we end up here.
            // so, we want to progress to the next step of reading :br_len:boot_supp:inher_prob
            if(readingBranchLength) {
                readingBranchLength = false;
                readingBootSupport = true;
            } else if(readingBootSupport) {
                readingBootSupport = false;
                readingGamma = true;
            } else if(readingGamma) {
                std::cout << "Error: Read a sequence of four colons (possibly with names/numbers in between some of them). This is not allowed by the format; quitting." << std::endl;
                exit(-1);
            } else {
                // begin reading a branch length
                readingBranchLength = true;
            }
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
                p->setMajorBranchLength(x);
                readingBranchLength = false;

                // Is there another ':' coming up? We expect them in batches of 3
                if(tokens[i+1] == ":") {
                    readingBootSupport = true;
                    i++;
                }
            } else if(readingBootSupport) {
                double x = std::stod(token);
                p->setBootSupport(x);
                readingBootSupport = false;

                // Is there another ':' coming up? We expect them in batches of 3
                if(tokens[i+1] == ":") {
                    readingGamma = true;
                    i++;
                }
            } else if(readingGamma) {
                double x = std::stod(token);
                p->setGamma(x);
                readingGamma = false;
            } else if(namingInternalNode) {
                p->setName(token);
            } else {
                Node* newNode = new Node;
                nodes.push_back(newNode);
                newNode->setMajorAnc(p);
                
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

    // Stitch the hybrid nodes together
    patchNetwork();

    // Add time information to the nodes
    setTimes();

    cout << std::endl << std::endl;
    listNodes();
    easyDbg();
}

void Network::setTimes(void) {
    // First, find all the leaves. We are storing branch lengths of INCOMING branches,
    // so it will be easiest to calculate node times going from the ground up.
    // This is also how ms expects times to be, which is the main reason this is being
    // implemented, so yay!
    std::vector<Node*> leaves;
    for(int i = 0; i < nodes.size(); i++) {
        if(nodes[i] == NULL)
            continue;
        
        if(nodes[i]->getLft() == NULL && nodes[i]->getRht() == NULL)
            leaves.push_back(nodes[i]);
    }

    // Now, go through and set the time on all nodes. Time is time from present.
    // We do this recursively because it is the easiest to implement
    for(Node *leaf : leaves) {
        leaf->setTime(0);
        setTimeRecur(leaf);
    }
}

void Network::setTimeRecur(Node *p) {
    Node *majAnc = p->getMajorAnc();
    Node *minAnc = p->getMinorAnc();

    // check if majAnc exists
    if(majAnc != NULL) {
        // if it exists and its time is not set, set its time and then continue recursively
        if(majAnc->getTime() == -1) {
            majAnc->setTime(p->getTime() + p->getMajorBranchLength());
            setTimeRecur(majAnc);
        }

        // if it exists and its time IS set, two options:
        //   1. kill this line of recursion here
        //   2. continue the recursion anyways and check for consistency (should only
        //      be used for debugging and testing)
        //
        // 2. is currently implemented.
        else {
            // return;
            if(majAnc->getTime() != p->getTime() + p->getMajorBranchLength()) {
                cout << "\t\tTIMES DO NOT MATCH - " << majAnc->getTime() << " != " << p->getTime() << " + " << p->getMajorBranchLength() << std::endl;
            }
            setTimeRecur(majAnc);
        }
    }
    if(minAnc != NULL) {
        // all the same stuff as with majAnc.
        // if it exists and its time is not set, set its time and then continue recursively
        if(minAnc->getTime() == -1) {
            minAnc->setTime(p->getTime() + p->getMinorBranchLength());
            setTimeRecur(minAnc);
        }

        // if it exists and its time IS set, two options:
        //   1. kill this line of recursion here
        //   2. continue the recursion anyways and check for consistency (should only
        //      be used for debugging and testing)
        //
        // 2. is currently implemented.
        else {
            // return;
            if(minAnc->getTime() != p->getTime() + p->getMinorBranchLength()) {
                cout << "\t\tTIMES DO NOT MATCH - " << minAnc->getTime() << " != " << p->getTime() << " + " << p->getMinorBranchLength() << std::endl;
            }
            setTimeRecur(minAnc);
        }
    }
}

std::vector<std::string> Network::parseNewick(std::string ns) {
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
void Network::patchNetwork() {
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

                // Make the MajorAnc of the dead node's child into p
                // The dead node's child should always be a left child, but better safe than sorry...
                Node *deadChild = dead->getLft() != NULL ? dead->getLft() : dead->getRht();
                deadChild->setMajorAnc(p);

                // Now make the child of the dead node p's child
                // Make sure we aren't overriding a node if p already has a child
                if(p->getLft() != NULL)
                    p->setRht(deadChild);
                else
                    p->setLft(deadChild);

                // We need to make sure that we set MajorAnc and MinorAnc properly!!!
                // MajorAnc will be the ancestor from which the node receives >50% of its genes, and MinorAnc
                // its complement. (In the case of a tie, it doesn't matter, so we don't do anything special...)
                if(dead->getGamma() > p->getGamma()) {
                    // Higher inheritance probability coming from the ancestor of dead
                    p->setMinorAnc(p->getMajorAnc());
                    p->setMinorBranchLength(p->getMajorBranchLength());
                    p->setMajorAnc(dead->getMajorAnc());
                    p->setMajorBranchLength(dead->getMajorBranchLength());

                    // Set the gammas for each of p's ancestors
                    if(p->getMajorAnc()->getLft() == p)
                        p->getMajorAnc()->setGammaLft(dead->getGamma());
                    else
                        p->getMajorAnc()->setGammaRht(dead->getGamma());
                    
                    // Now the minor ancestor
                    if(p->getMinorAnc()->getLft() == p)
                        p->getMinorAnc()->setGammaLft(p->getGamma());
                    else
                        p->getMinorAnc()->setGammaRht(p->getGamma());
                } else {
                    // Higher inheritance probability coming from the ancestor of p
                    // p's MajorAnc is already correct, so just set the MinorAnc
                    p->setMinorAnc(dead->getMajorAnc());
                    p->setMinorBranchLength(dead->getMajorBranchLength());

                    // Set the gammas for each of p's ancestors
                    if(p->getMajorAnc()->getLft() == p)
                        p->getMajorAnc()->setGammaLft(p->getGamma());
                    else
                        p->getMajorAnc()->setGammaRht(p->getGamma());
                    
                    // Now the minor ancestor
                    if(p->getMinorAnc()->getLft() == p)
                        p->getMinorAnc()->setGammaLft(dead->getGamma());
                    else
                        p->getMinorAnc()->setGammaRht(dead->getGamma());
                }

                // Set p as the proper child of its new ancestor
                if(dead->getMajorAnc()->getLft() == dead)
                    dead->getMajorAnc()->setLft(p);
                else
                    dead->getMajorAnc()->setRht(p);

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
            if (p->getMajorAnc() == NULL)
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
    for(int i = 0; i < nodes.size(); i++) {
        std::cout << "Node " << i << ": " << nodes[i]->getName() << std::endl << std::flush;
        nodes[i]->printInfo();
    }
    cout << "Hybrids: " << hybrids << std::endl;
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
    //       assign MinorAnc() in a way that can be used to maintain directionality.
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

            Node *majAnc = p->getMajorAnc();
            Node *minAnc = p->getMinorAnc();
            // a. if only one ancestor (this will only ever take place when there is a MajorAnc and NOT a MinorAnc)
            if(majAnc != NULL && minAnc == NULL) {
                // I. if the anc is NOT named: take the anc
                if(blankName(majAnc)) {
                    majAnc->setName(p->getName());

                    // Remove p and add majAnc
                    removeMe.push_back(i);
                    addMe.push_back(majAnc);
                }
                // II. if anc named: coalesce *INTO* the ancestor
                else {
                    cout << "-ej " << majAnc->getTime() << " " << p->getName() << " " << majAnc->getName() << std::endl;

                    // Remove p; majAnc is already named, so if it still has more potential things to do, it should be
                    // in activeNodes already. So we don't need to touch it.
                    removeMe.push_back(i);
                }
            }
            // b. if we have two ancestors
            else if(majAnc != NULL && minAnc != NULL) {
                // I. both ancs unnamed: split, give MajorAnc our name and give MinorAnc a new name
                if(blankName(majAnc) && blankName(minAnc)) {
                    // p is coalescing into majAnc, so gamma on this split comes from majAnc
                    double gamma = (majAnc->getLft() == p) ? majAnc->getGammaLft() : majAnc->getGammaRht();

                    // Split
                    cout << "-es " << p->getTime() << " " << p->getName() << " " << gamma << std::endl;

                    // Give MajorAnc our name and add MajorAnc to activeNodes
                    majAnc->setName(p->getName());
                    addMe.push_back(majAnc);

                    // Give MinorAnc its new name and add it to activeNodes
                    minAnc->setName(popnCounter++);
                    addMe.push_back(minAnc);

                    // Remove p from activeNodes
                    removeMe.push_back(i);
                }
                // II. if only one anc is named: split, give the blank anc a new name, and then join with the named anc
                else if(blankName(majAnc) != blankName(minAnc)) {
                    Node *namedAnc = blankName(majAnc) ? minAnc : majAnc;
                    Node *unnamedAnc = blankName(majAnc) ? majAnc : minAnc;

                    // p is coalescing into unnamedAnc, so gamma comes from unnamedAnc 
                    double gamma = (unnamedAnc->getLft() == p) ? unnamedAnc->getGammaLft() : unnamedAnc->getGammaRht();

                    // Split
                    cout << "-es " << p->getTime() << " " << p->getName() << " " << gamma << std::endl;

                    // Give the blank anc a new name and add it to the queue
                    unnamedAnc->setName(popnCounter++);
                    addMe.push_back(unnamedAnc);

                    // Join
                    cout << "-ej " << namedAnc->getTime() << " " << p->getName() << " " << namedAnc->getName() << std::endl;

                    // Remove p
                    removeMe.push_back(i);
                }
                // III. if both ancs are named, split, then join both the new lineages with the existing ones
                else {
                    // p is coalescing into majAnc, so gamma comes from majAnc
                    double gamma = (majAnc->getLft() == p) ? majAnc->getGammaLft() : majAnc->getGammaRht();

                    // Split
                    cout << "-es " << p->getTime() << " " << p->getName() << " " << gamma << std::endl;

                    // Join left
                    cout << "-ej " << majAnc->getTime() << " " << p->getName() << " " << majAnc->getName() << std::endl;

                    // Join right
                    cout << "-ej " << minAnc->getTime() << " " << popnCounter++ << " " << minAnc->getName() << std::endl;

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
#include "Node.hpp"
#include "Network.hpp"
#include "MSEvents.hpp"

#include <iostream>
#include <iomanip>
#include <set>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <typeinfo>

using namespace std;
std::string BLANK_NAME = "__$!&*#%*__";

// IDEA: When building the network from MSEvents:
//    1. build the network NOT accounting for the fact that we will have additional nodes
//    2. after the network is built, just do a traversal of the network and anytime we find
//       a node that has exactly 1 ancestor and exactly 1 child: remove it.

Network::Network(std::vector<MSEvent*> events) {
    // First, order the events backwards in time (from tips to root) so that they
    // can be read in order.
    sort(events.begin(), events.end(), [](MSEvent *a, MSEvent *b) {
        return a->getTime() < b->getTime();
    });

    // Start by figuring out how many leaves we have. This is equal to the largest number
    // in the MSEvent's minus the total number of MSSplitEvent's
    int max = 0;
    int splits = 0;
    cout << "join: " << join << ", split: " << split << std::endl;
    for(MSEvent *e : events) {
        cout << e->getEventType() << std::endl;
        // Add a split. Splits necessarily create a new, larger number, so we don't need to
        // check these events for a maximum number
        if(e->getEventType() == split)
            splits++;
        else {
            max = std::max(std::max(max, ((MSJoinEvent*)e)->getMajorTaxa()), ((MSJoinEvent*)e)->getMinorTaxa());
        }
    }
    int ntaxa = max - splits;
    cout << "ntaxa = " << max << "-" << splits << " = " << (max-splits) << std::endl;

    // activeNodes includes the most recent version of any given taxa. I.e. after a split event where taxa
    // 3 is turned into taxa 3 & taxa 7, the old node for taxa 3 will be removed, and both the new nodes added.
    // Also, if taxa 3 & taxa 4 are joined with command -ej <t> 3 4 (i.e. deleting taxa 3), the destroyed node
    // for taxa 3 will remain in activeNodes, but the old version of taxa 4 will be removed and the new one added.
    //
    // We have to do things in this way because otherwise we run into issues when a node is both itself a hybrid
    // and apart of a hybridization event.
    std::vector<Node*> activeNodes;
    
    // These are the intermediate nodes that we are required to create when split events occurred, but that we
    // ultimately don't want to leave left over in the network. They are removed from nodes at the end of this section
    std::vector<Node*> globalRemoveMe;

    // Populate the nodes list with the leaves
    for(int i=1; i<=ntaxa; i++) {
        Node *p = new Node();
        p->setName(i);
        p->setTime(0);
        nodes.push_back(p);
        activeNodes.push_back(p);
    }
    
    // Now, just process each event, and build the network backwards in time
    for(unsigned int i=0; i < events.size(); i++) {
        // This line for debugging
        std::cout << "activeNodes: {";
        for(Node *n : activeNodes)
            std::cout << n->getName() << ", ";
        std::cout << "}" << std::endl;
        //

        cout << "i=" << i << std::endl;
        // Because we are doing things in chronological order, the time that a join event occurs should be
        // greater than or equal to the times of *each* of the nodes involved in the event. So, all we have
        // to do is create a parent pointing to each of the nodes involved.
        if(events[i]->getEventType() == join) {
            std::cout << "JOIN" << std::endl;
            MSJoinEvent *e = (MSJoinEvent*)events[i];

            // Both taxa involved in the join event should already exist
            Node *fromTaxa = NULL;
            Node *toTaxa = NULL;
            int fromIdx;  // used for removing the node from activeNodes

            std::cout << "Looking for nodes with names \"" << std::to_string(e->getMajorTaxa()) << "\" and \"" << e->getMinorTaxa() << "\"..." << std::endl;
            for(unsigned int i=0; i<activeNodes.size(); i++) {
                if(fromTaxa == NULL && activeNodes[i]->getName().compare(std::to_string(e->getMajorTaxa())) == 0) {
                    fromTaxa = activeNodes[i];
                    fromIdx = i;
                } else if(toTaxa == NULL && activeNodes[i]->getName().compare(std::to_string(e->getMinorTaxa()))) {
                    toTaxa = activeNodes[i];
                }

                if(fromTaxa != NULL && toTaxa != NULL)
                    break;
            }

            // If we couldn't find one or both of the taxa involved, this is an error; quit
            if(fromTaxa == NULL || toTaxa == NULL) {
                std::cout << "ERROR: When finding both taxa in a join event, one or both taxa could not be found; quitting." << std::endl;
                exit(-1);
            }

            // Create a new parent node pointing back to both of these nodes
            Node *p = new Node();
            p->setTime(e->getTime());
            p->setName(fromTaxa->getName());

            // If either fromTaxa or toTaxa (or both) is (are) the product of splits, then we want to get rid of them right now
            // and form connections from p to the correct, corresponding children of this (these) node(s)
            //
            // What we need to do to transfer the information contained in each intermediate node:
            //    1. copy over the gamma value (stored in the intermediate node inside node->getGamma())
            //    2. re-route the ancestors of the intermediate node's child (should only ever be 1 child)
            if(fromTaxa->justSplit()) {
                // 1. copy over the gamma value
                // (we always set fromTaxa as the left child)
                p->setGammaLft(fromTaxa->getGamma());

                // 2. re-route the ancestors of the intermediate node's child (should only ever have 1 child)
                Node *child = fromTaxa->getLft();
                if(child == NULL || fromTaxa->getRht() != NULL) {
                    // A node which is the product of a split should only ever have 1 child, and that child should always be left
                    std::cout << "ERROR: A node which is the product of a split had either 2 children, or the left child was NULL; quitting." << std::endl;
                    exit(-1);
                }
                
                // We don't know which ancestor we are, so we have to figure that out before creating ancestor connections
                if(child->getMajorAnc() == fromTaxa) {
                    child->setMajorAnc(p);
                    child->setMajorBranchLength(p->getTime() - child->getTime());
                } else {
                    child->setMinorAnc(p);
                    child->setMinorBranchLength(p->getTime() - child->getTime());
                }

                globalRemoveMe.push_back(fromTaxa);
            } else {
                // If the node is not the product of a split, then creating ancestor connections is more straightforward.
                //
                // Make the children point to their new parent and
                // Set the branch lengths. Remember, branch lengths are INCOMING, not outgoing.
                fromTaxa->setMajorAnc(p);
                fromTaxa->setMajorBranchLength(p->getTime() - fromTaxa->getTime());
            }
            if(toTaxa->justSplit()) {
                // 1. copy over the gamma values
                // (we always set toTaxa as the right child)
                p->setGammaRht(toTaxa->getGamma());

                // 2. re-route the ancestors of the intermediate node's child (should only ever have 1 child)
                Node *child = toTaxa->getLft();
                if(child == NULL || toTaxa->getRht() != NULL) {
                    // A node which is the product of a split should only ever have 1 child, and that child should always be left
                    std::cout << "ERROR: A node which is the product of a split had either 2 children, or the left child was NULL; quitting." << std::endl;
                    exit(-1);
                }

                // We don't know which ancestor we are, so we have to figure that out before creating ancestor connections
                if(child->getMajorAnc() == toTaxa) {
                    child->setMajorAnc(p);
                    child->setMajorBranchLength(p->getTime() - child->getTime());
                } else {
                    child->setMinorAnc(p);
                    child->setMinorBranchLength(p->getTime() - child->getTime());
                }

                globalRemoveMe.push_back(toTaxa);
            } else {
                // If the node is not the product of a split, then creating ancestor connections is more straightforward.
                //
                // Make the children point to their new parent and
                // Set the branch lengths. Remember, branch lengths are INCOMING, not outgoing.
                toTaxa->setMajorAnc(p);
                toTaxa->setMajorBranchLength(p->getTime() - toTaxa->getTime());
            }

            // Now fromTaxa and toTaxa are the nodes that we actually want to connect to, so we create the connections
            // There shouldn't be an issue always setting fromTaxa as left and toTaxa as right...I believe.
            p->setLft(fromTaxa);
            p->setRht(toTaxa);

            // Remove fromTaxa from activeNodes
            activeNodes.erase(std::next(activeNodes.begin(), fromIdx));

            // Add p to activeNodes AND nodes
            activeNodes.push_back(p);
            nodes.push_back(p);
        } else if(events[i]->getEventType() == split) {
            std::cout << "SPLIT" << std::endl;
            MSSplitEvent *e = (MSSplitEvent*)events[i];
            // Split events are straightforward, but require the use of an intermediate node that will be removed later.
            // This is because ms's protocol creates a new node everytime a split occurs, we essentially end up with a
            // duplicate node in our network. There are many ways that this could potentially be dealt with, but this is
            // how we choose to deal with it.
            Node *p = NULL;
            int pIdx;

            std::cout << "Looking for node with name \"" << std::to_string(e->getTaxa()) << "\"..." << std::endl;
            for(unsigned int i=0; i<activeNodes.size(); i++) {
                if(activeNodes[i]->getName().compare(std::to_string(e->getTaxa())) == 0) {
                    p = activeNodes[i];
                    pIdx = i;
                    break;
                }
            }
            
            // If p is still NULL, this is an error; quit
            if(p == NULL) {
                std::cout << "ERROR: Node involved in split event not found when split event was reached; quitting." << std::endl;
                exit(-1);
            }

            // If p is a leaf then we actually need to create a new node
            // getLft will equal getRht exactly when both are NULL.
            if(p->getLft() == p->getRht()) {
                Node *newNode = new Node();
                newNode->setName(p->getName());
                newNode->setTime(e->getTime());
                newNode->setLft(p);
                p->setMajorAnc(newNode);
                p->setMajorBranchLength(newNode->getTime());  // p->getTime() is 0
                nodes.push_back(newNode);

                // Remove the leaf from activeNodes
                std::cout << "pIdx: " << pIdx << std::endl;
                activeNodes.erase(std::next(activeNodes.begin(), pIdx));

                // Now newNode is the node we want to operate on
                p = newNode;
                activeNodes.push_back(p);
                pIdx = activeNodes.size() - 1;
            } else {
                // This event corresponds to Node p splitting, so this event gives us p's time
                p->setTime(e->getTime());
            }

            // p might have just come from a split. In this case, its children do not have branch lengths yet.
            if(p->justSplit()) {
                // Hybrid ancestors only ever have 1 child, and we set the child to be the left child,
                // so that's the only node we need to check. (We shouldn't actually even have to check it,
                // it should always be there)
                if(p->getLft() != NULL) {
                    if(p == p->getLft()->getMajorAnc())
                        // p is the majorAnc of its left child
                        p->getLft()->setMajorBranchLength(p->getTime() - p->getLft()->getTime());
                    else
                        // p is the minorAnc of its left child
                        p->getLft()->setMinorBranchLength(p->getTime() - p->getLft()->getTime());
                } else {
                    std::cout << "ERROR: p should have a left child here. Unsure what is going on; quitting." << std::endl;
                    exit(-1);
                }
            }

            // Create both of the new nodes and mark them as having just come from a split
            // NOTE: We do NOT know the times on EITHER the major or the minor hybrid nodes. We WILL know
            //       those times once either: 1) they split (the split event gives us their time), or
            //       2) they are joined with something else. In this case, the nodes we make here are
            //       only placeholders and will just be removed anyways, so they don't need a time.
            Node *maj = new Node();
            maj->setName(p->getName());
            maj->setLft(p);
            maj->setGamma(e->getGamma());
            maj->setJustSplit(true);

            Node *min = new Node();
            min->setName((ntaxa++)+1);
            min->setLft(p);
            min->setGamma(1-e->getGamma());
            min->setJustSplit(true);

            // Link p to its ancestors
            p->setMajorAnc(maj);
            p->setMinorAnc(min);

            // We do not know the times on either ancestor, so we cannot set p's branch length yet!
            // Remove p from activeNodes
            activeNodes.erase(std::next(activeNodes.begin(), pIdx));
            
            // Add min and maj to activeNodes
            activeNodes.push_back(min);
            activeNodes.push_back(maj);

            // Add min and maj to nodes
            nodes.push_back(min);
            nodes.push_back(maj);
        } else {
            std::cout << "ERROR: MSEvent is neither a split or join event; quitting." << std::endl;
            exit(-1);
        }
    }

    ////////
    std::cout << "\tnodes[5] info.\n\t\tnodes[5]->justSplit(): " << nodes[5]->justSplit() << "\n\t\tnodes[5]->getGamma(): " << nodes[5]->getGamma() << std::endl;
    ////////

    // Now we are going to remove the placeholder nodes that we had to add for split events
    // We remove by index, so first we have to find the indices of all the nodes
    std::vector<int> removeIdcs;
    for(unsigned int i=0; i < nodes.size(); i++) {
        for(unsigned int j=0; j < globalRemoveMe.size(); j++) {
            if(nodes[i] == globalRemoveMe[j]) {
                // Add the idx to our list
                removeIdcs.push_back(i);

                // Remove the node from our search space
                globalRemoveMe.erase(std::next(globalRemoveMe.begin(), j));

                // Break out of the loop that is looking for this node
                break;
            }
        }
    }

    ////////
    for(int idx : removeIdcs)
        std::cout << "\t" << idx << ", ";
    std::cout << std::endl;
    ////////

    // Sort the list of indices in descending order so that we can removing them without any issues
    std::sort(removeIdcs.begin(), removeIdcs.end(),
        [](int a, int b) {
            return(a > b);
        }
    );

    // This line for debugging
    std::cout << "nodes: {";
    for(Node *n : nodes)
        std::cout << n->getName() << ", ";
    std::cout << "}" << std::endl;
    //

    nodes[5]->printInfo();

    // Now, just loop through and remove the nodes
    for(int idx : removeIdcs) {
        nodes.erase(std::next(nodes.begin(), idx));
    }

    // This line for debugging
    std::cout << "nodes: {";
    for(Node *n : nodes)
        std::cout << n->getName() << ", ";
    std::cout << "}" << std::endl;
    //

    // All of the nodes in the network have indistinguishable names right now, so let's fix that
    postmsParseNetworkRename();
}

void Network::postmsParseNetworkRename(void) {
    // Internal nodes are given lower case letter names: a, b, c, ..., z, za, zb, ..., zz, ..., zzz, ...
    // Leaves are given natural number names: 1, 2, 3, ...

    // First: find root.
    Node *root = nodes[0];
    while(root->getMajorAnc() != NULL)
        root = root->getMajorAnc();


}

Network::Network(std::string newickStr) {
    hybrids = 0;
    std::vector<std::string> tokens = parseNewick(newickStr);

    // build up the Network from the parsed Newick string
    bool readingBranchLength = false;
    bool readingBootSupport = false;
    bool readingGamma = false;
    bool namingInternalNode = false;
    Node* p = NULL;
    for(unsigned int i=0; i<tokens.size(); i++) {
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
    for(unsigned int i=0; i<nodes.size(); i++) {
        if (nodes[i]->getLft() == NULL)
            nodes[i]->setIndex(ndeIdx++);
    }
    
    for(unsigned int i=0; i<nodes.size(); i++) {
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
    for(unsigned int i = 0; i < nodes.size(); i++) {
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
    for(unsigned int i = 0; i < ns.size(); i++) {
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
    for(unsigned int i = 0; i < tks.size(); i++) {
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
    for(unsigned int i = 0; i < nodes.size(); i++) {
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
    for(unsigned int i = 0; i < list.size(); i++) {
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
    cout << "Calling the Network destructor with " << nodes.size() << " nodes.\n";
    for(unsigned int i = 0; i < nodes.size(); i++)
        delete nodes[i];
}

void Network::listNodes(void) {
    for(unsigned int i = 0; i < nodes.size(); i++) {
        std::cout << "Node " << i << ": " << std::endl << std::flush;
        nodes[i]->printInfo();
    }
    cout << "Hybrids: " << hybrids << std::endl;
}

std::vector<MSEvent*> Network::toms(void) {
    std::vector<MSEvent*> events;
    // First, gather & name all leaves while setting the names of non-leaves to null
    std::vector<Node*> activeNodes;
    int popnCounter = 1;
    for(unsigned int i = 0; i < nodes.size(); i++) {
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
        unsigned int i;
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
                    events.push_back(new MSJoinEvent(majAnc->getTime(), p->getName(), majAnc->getName()));

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
                    events.push_back(new MSSplitEvent(p->getTime(), p->getName(), gamma));

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
                    events.push_back(new MSSplitEvent(p->getTime(), p->getName(), gamma));

                    // Give the blank anc a new name and add it to the queue
                    unnamedAnc->setName(popnCounter++);
                    addMe.push_back(unnamedAnc);

                    // Join
                    events.push_back(new MSJoinEvent(namedAnc->getTime(), p->getName(), namedAnc->getName()));

                    // Remove p
                    removeMe.push_back(i);
                }
                // III. if both ancs are named, split, then join both the new lineages with the existing ones
                else {
                    // p is coalescing into majAnc, so gamma comes from majAnc
                    double gamma = (majAnc->getLft() == p) ? majAnc->getGammaLft() : majAnc->getGammaRht();

                    // Split
                    events.push_back(new MSSplitEvent(p->getTime(), p->getName(), gamma));

                    // Join left
                    events.push_back(new MSJoinEvent(majAnc->getTime(), p->getName(), majAnc->getName()));

                    // Join right
                    events.push_back(new MSJoinEvent(minAnc->getTime(), popnCounter++, minAnc->getName()));

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

    return events;
}

bool Network::blankName(Node *p) {
    return p->getName().compare(BLANK_NAME) == 0;
}

int Network::activeNodesIdx(Node *p, std::vector<Node*> list) {
    for(unsigned int i = 0; i < list.size(); i++) {
        if(list[i] == p)
            return i;
    }
    return -1;
}
#ifndef Node_hpp
#define Node_hpp

#include <string>
#include <iostream>

class Node {
    public:
                Node(void);
        Node*   getLft(void) { return left; }
        Node*   getRht(void) { return right; }

        Node*   getPrimAnc(void) { return primAncestor; }
        Node*   getHybAnc(void) { return hybAncestor; }

        int     getIndex(void) { return index; }
        std::string getName(void) { return name; }
        double getBranchLength(void) { return branchLength; }
        void setLft(Node* p) { left = p; }
        void setRht(Node* p) { right = p; }
        
        void setPrimAnc(Node* p) { primAncestor = p; }
        void setHybAnc(Node* p) { hybAncestor = p; }

        void setIndex(int x) { index = x; }
        void setName(std::string s) { name = s; }
        void setName(int i) { name = std::to_string(i); }
        void setBranchLength(double x) { branchLength = x; }
        void print(void);
        void setTime(double);
        double getTime();
        bool touched = false;
    
    protected:
        Node* left;
        Node* right;
        Node* primAncestor;
        Node* hybAncestor;
        int index;
        std::string name;
        double branchLength;
        double time;
};

#endif
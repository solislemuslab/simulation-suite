#ifndef Node_hpp
#define Node_hpp

#include <string>
#include <iostream>

class Node {
    public:
                Node(void);
        Node*   getLft(void) { return left; }
        Node*   getRht(void) { return right; }

        Node*   getMajorAnc(void) { return majorAncestor; }
        Node*   getMinorAnc(void) { return minorAncestor; }

        int     getIndex(void) { return index; }
        std::string getName(void) { return name; }
        double getBranchLength(void) { return branchLength; }
        void setLft(Node* p) { left = p; }
        void setRht(Node* p) { right = p; }
        
        void setMajorAnc(Node* p) { majorAncestor = p; }
        void setMinorAnc(Node* p) { minorAncestor = p; }

        void setGammaLft(double g) { gammaLft = g; }
        void setGammaRht(double g) { gammaRht = g; }
        double getGammaLft(void) { return gammaLft; }
        double getGammaRht(void) { return gammaRht; }
        void setGamma(double g) { gamma = g; }
        double getGamma(void) { return gamma; }

        void setBootSupport(double supp) { bootSupport = supp; }
        double getBootSupport(void) { return bootSupport; }

        void setIndex(int x) { index = x; }
        void setName(std::string s) { name = s; }
        void setName(int i) { name = std::to_string(i); }
        void setBranchLength(double x) { branchLength = x; }
        void printInfo(void);
        void setTime(double);
        double getTime();
        bool touched = false;
    
    protected:
        Node* left;
        Node* right;
        Node* majorAncestor;
        Node* minorAncestor;
        int index;
        std::string name;
        double branchLength;
        double time;
        double gammaLft = 0;
        double gammaRht = 0;
        double bootSupport = -1;
        double gamma = 0;
};

#endif
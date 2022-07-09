#include "../core/Network.hpp"
#include "../SimSuite.hpp"
#include "CoreTests.hpp"

#include <iostream>
#include <string>
#include <chrono>
#include <ctime>

std::vector<std::string> ALL_NEWICKS{
    "((1,((2,(3,(4)Y#H1)g)e,(((Y#H1,5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;",
    "((1:0.1,((2:0.2,(3:0.3,(4:0.4)Y#H1:0.5)g:0.6)e:0.7,(((Y#H1:0.8,5:0.9)h:1.0,6:1.1)f:1.2)X#H2:1.3)c:1.4)a:1.5,((X#H2:1.6,7:1.7)d:1.8,8:1.9)b:2.0)r;",
    "(15:11.0,(1:10.0,((14:8.0,(((7:2.8,(10:1.6,(9:0.4,8:0.4):1.2):1.2):0.8,(11:2.8,(13:0.4,12:0.4):2.4):0.8):3.4,#H1:0.4::0.3):1.0):1.2,(((2:0.4,3:0.4):5.2,((4:3.6,5:3.6):1.2,6:4.8):0.8):1.0)#H1:2.6::0.7):0.8):1.0);"
};
int TESTS_RUN = 0;
int TESTS_PASSED = 0;

int main(int narg, char **argv) {
    auto start = std::chrono::system_clock::now();
    TEST_ALL();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> time = end-start;

    std::cout << TESTS_PASSED << "/" << TESTS_RUN << " tests passed in " << time.count() << "s.\n";
}

void TEST_isomorphicNewick_TRUE(void) {
    int nRand = 10;
    for(std::string newick : ALL_NEWICKS) {
        Network net(newick, "newick");
        std::vector<std::string> randomNewicks = net.getRandomNewickRepresentations(nRand);
        bool allPassed = true;

        for(int i=0; i < nRand-1; i++) {
            for(int j=i+1; j<nRand; j++) {
                if(!ASSERT(isomorphicNewick(randomNewicks[i], randomNewicks[j]), "isomorphicNewick_TRUE", randomNewicks[i] + " " + randomNewicks[j]))
                    allPassed = false;
            }
        }

        TESTS_RUN += 1;
        if(allPassed)
            TESTS_PASSED += 1;
    }
}

void TEST_isomorphicNewick_FALSE(void) {
    // Obvious falses, comparing totally different networks
    int nRand = 10;
    for(unsigned int n=0; n < ALL_NEWICKS.size(); n++) {
        std::vector<std::string> randoms1 = Network(ALL_NEWICKS[n], "newick").getRandomNewickRepresentations(nRand);
        std::vector<std::string> randoms2 = Network(ALL_NEWICKS[n+1], "newick").getRandomNewickRepresentations(nRand);

        bool allPassed = true;
        for(int i=0; i<nRand; i++) {
            for(int j=0; j<nRand; j++) {
                if(!ASSERT(isomorphicNewick(randoms1[i], randoms2[j]), "isomorphicNewick_FALSE", randoms1[i] + " " + randoms2[j]))
                    allPassed = false;
            }
        }

        TESTS_RUN += 1;
        if(allPassed)
            TESTS_PASSED += 1;
    }
}

// Returns total number of tests passed
void TEST_ALL(void) {
    TEST_isomorphicNewick_TRUE();
}

bool ASSERT(bool assertion, std::string title, std::string info) {
    if(!assertion) {
        std::cerr << "ASSERT FAILED FOR " + title + ": " + info;
        exit(0);
    }
    return assertion;
}
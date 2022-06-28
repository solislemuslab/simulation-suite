#include <iostream>
#include <string>


int main(int narg, char **argv) {
    
}

void TEST_ALL(void) {
    
}


void ASSERT(bool assertion, std::string failMsg) {
    if(!assertion) {
        std::cerr << failMsg;
        exit(0);
    }
}
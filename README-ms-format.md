### Read a network in with Newick format

This is in the `main.cpp` file right now, so `main.cpp` should just be adjusted for now to test different Newick strings (do not commit your local changes):

```cpp
#include "SimSuite.hpp"
#include <iostream>

int main(int narg, char **argv) {
    // Convert extended Newick to MS
    std::string msString = newickToMS("(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;");
    std::cout << msString << std::endl;
}
```

### Running the code

```shell
git clone https://github.com/solislemuslab/simulation-suite.git
cd simulation-suite/src/

# For now, msString must be adjusted in main.cpp to test different Newick strings
make
./main
```

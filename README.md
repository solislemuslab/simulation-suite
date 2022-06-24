# simulation-suite

### TODO:

1. Fix the `Network::toms(void)` function so that it is not destructive, and leaves the original network intact
2. Add times to the output of `Network::toms(void)`
4. ~~Fix the Newick parser to work with the extended Newick format described in [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl/wiki/Introduction)~~
5. Have `Network::toms(void)` return a string containing the *entire* ms command line command, not just the `-es` and `-ej` parameters
6. Lots of testing!!!!

### Read a network in with Newick format

This is in the `main.cpp` file right now, so `main.cpp` should just be adjusted for now to test different Newick strings:

```cpp
#include "Network.hpp"

int main(int narg, char *argv[]) {
  // Read the network from Newick
  std::string myNewick = "(((A:13.2:45:,(B::27:)#H1:::0.9),(C,#H1:::0.1)f)g,D);";
  Network myNetwork(myNewick);
}
```

### Get the `ms` parameters that correspond to your network

```cpp
#include "Network.hpp"

int main(int narg, char *argv[]) {
  ...
  
  // Print out the ms parameters that go with the network
  myNetwork.toms();
}
```

### Running the code

```shell
cd src/
make
./main
```

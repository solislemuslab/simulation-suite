# simulation-suite

### TODO:

1. Fix the `Network::toms(void)` function so that it is not destructive, and leaves the original network intact
2. Add times to the output of `Network::toms(void)`
3. Fix the Newick parser to work with the extended Newick format described in [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl/wiki/Introduction)
4. Have `Network::toms(void)` return a string containing the *entire* ms command line command, not just the `-es` and `-ej` parameters
5. Lots of testing!!!!

### Read a network in with Newick format

This is in the `main.cpp` file right now, so `main.cpp` should just be adjusted for now to test different Newick strings:

```cpp
#include "Network.hpp"

int main(int narg, char *argv[]) {
  // Read the network from Newick
  std::string myNewick = "((1,((2,(3,(4)Y#H1)g)e,(((Y#H1,5)h,6)f)X#H2)c)a,((X#H2,7)d,8)b)r;";
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
make
./main
```

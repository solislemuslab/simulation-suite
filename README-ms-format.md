# simulation-suite

### TODO:

1. Fix the `Network::toms(void)` function so that it is not destructive, and leaves the original network intact
2. ~~Add times to the output of `Network::toms(void)`~~
3. ~~Add gamma to the output of `Network::toms(void)`~~
4. ~~Fix the Newick parser to work with the extended Newick format described in [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl/wiki/Introduction)~~
5. Have `Network::toms(void)` return a string containing the *entire* ms command line command, not just the `-es` and `-ej` parameters
6. ~~Write functionality to grow a network from an ms command (primarily used for testing)~~
7. Lots of testing!!!!

### Read a network in with Newick format

This is in the `main.cpp` file right now, so `main.cpp` should just be adjusted for now to test different Newick strings:

```cpp
#include "Network.hpp"

int main(int narg, char *argv[]) {
  // Read the network from Newick
  std::string myNewick = "(((A:1,(B:0.3)#H1:0.7::0.9)e:1,(C:0.5,#H1:0.2::0.1)f:1.5)g:0.3,D:2.3)h;";
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

Running this with the network given in the Newick string above yields the output:

```shell
-es 0.3 2 0.9
-ej 1 2 1
-ej 0.5 5 3
-ej 2 3 1
-ej 2.3 1 4
```

### Running the code

```shell
cd src/
make
./main
```

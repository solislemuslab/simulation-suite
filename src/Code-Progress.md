# Code Progress TODO Checklist

### core code:

1. Fix the `Network::toms(void)` function so that it is not destructive, and leaves the original network intact
2. ~~Add times to the output of `Network::toms(void)`~~
3. ~~Add gamma to the output of `Network::toms(void)`~~
4. ~~Fix the Newick parser to work with the extended Newick format described in [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl/wiki/Introduction)~~
5. Have `Network::toms(void)` return a string containing the *entire* ms command line command, not just the `-es` and `-ej` parameters
6. ~~Write functionality to grow a network from an ms command (primarily used for testing)~~
7. Lots of testing!!!!
8. Add to function to `Network.cpp` that returns the extended Newick string for the network.
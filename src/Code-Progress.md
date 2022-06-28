# Code Progress TODO Checklist

### core code:

1. ~~Fix the `Network::toms(void)` function so that it is not destructive, and leaves the original network intact~~
2. Add functionality to read Newick from a file and then convert that to ms (`Network::newickFileToMS(std::string location)`)
3. ~~Add times to the output of `Network::toms(void)`~~
4. ~~Add gamma to the output of `Network::toms(void)`~~
5. ~~Fix the Newick parser to work with the extended Newick format described in [PhyloNetworks](https://github.com/crsl4/PhyloNetworks.jl/wiki/Introduction)~~
6. ~~Have `Network::toms(void)` return a string containing the *entire* ms command line command, not just the `-es` and `-ej` parameters~~
7. ~~Write functionality to grow a network from an ms command (primarily used for testing)~~
8. Lots of testing!!!!
9. ~~Add function to `Network.cpp` that returns the extended Newick string for the network.~~
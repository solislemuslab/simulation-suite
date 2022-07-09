# Code Progress TODO Checklist

### core code:

1. Adjust code so that if taxon are already named 1, 2, 3, 4, ... then they are not overwritten in the process of converting from Newick --> ms
2. Add functionality to read Newick from input instead of taking a hard-coded Newick string from main
3. Add functionality to read Newick from a file and then convert that to ms (`Network::newickFileToMS(std::string location)`)
4. Finish setting up automated testing
5. Lots of testing!!!!
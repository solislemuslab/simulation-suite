# Code Progress TODO Checklist

### core code:

1. ~~Adjust code so that if taxon are already named 1, 2, 3, 4, ... then they are not overwritten in the process of converting from Newick --> ms~~
2. ~~Add a warning for when a hybrid node is included without a gamma specified, or with a gamma specified as 0 or 1.~~
2. ~~Fix Newick parser to be more ammenable to whitespace~~
3. ~~Add a warning to the Newick parser when branch lengths are not specified~~
4. Add functionality to read Newick from input instead of taking a hard-coded Newick string from main
5. Add functionality to read Newick from a file and then convert that to ms (`SimSuite::newickFileToMS(std::string location)`)
6. Finish setting up automated testing
7. Lots of testing!!!!
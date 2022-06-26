# Simulation suite
## Compilation of scripts to simulate data on phylogenetic networks

I want to
- Simulate phylogenetic networks under a birth-death-hybridization model using `SiPhyNetworks` R package
    - Use `scripts/simulate-networks.R`
    - Dependencies:
        - R (version 4.1.2 or higher)
        - [SiPhyNetwork](https://github.com/jjustison/SiPhyNetwork): 
        ```
        install.packages("devtools")
        devtools::install_github("jjustison/SiPhyNetwork")
        ```
    - Input parameters:
        - N: number of networks to simulate
        - lambda: speciation rate
        - mu: extinction rate
        - nu: hybridization rate
        - hybprobs: 3d vector with frequencies of each of the three types of hybridizations (see [here](https://github.com/jjustison/SiPhyNetwork/blob/master/vignettes/figures/hybridization_types.png))
        - ntips: number of tips in the network
        - seed: random seed for reproducibility
    - Output:
        - N output files each with a simulated network in parenthetical format. Note that branches on the network depend on the units of the rates
        - text file with parameters selected for reproducibility

- Simulate gene trees on a given species network (or tree) under the multispecies coalescent using [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html)
    - Use `scripts/simulate-gene-trees-ms.jl`
    - Dependencies:
        - `ms`:
            - Download `ms.tar.gz` [here](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)
            - Move the folder to wherever you want it
            - Untar `tar -xvf ms.tar.gz` to create the `msdir` folder
            - Inside `msdir`, compile 
            ```
            gcc -o ms ms.c streec.c rand2.c -lm
            ```
            Note that we use `rand2.c` so that the number generator used is `rand()`, not `drand48()` which is supposed to be inferior.
            - Copy the `ms` executable in `/usr/local/bin` or somewhere in your PATH, or add the path to it to PATH
        - Julia (version 1.8 or higher)
        - PhyloNetworks package, see [here](https://github.com/crsl4/PhyloNetworks.jl). In Julia:
        ```julia
        ] add PhyloNetworks
        ```
    - Input:
        - msnet: text file with the network in ms format, see [ms docs](http://home.uchicago.edu/~rhudson1/source/mksamples.html). Branch lengths are assumed to be in coalescent units. This file containts the -ej and -es events
        - net: text file with the network in parenthetical format.
        - numalleles: number of alleles per taxon. This script assumes all taxa have the same
               number of alleles, but see notes in the jl script).
        - numgt: number of gene trees to simulate
    - Output:
        - treefile: text file with simulated gene trees (one per line) in parenthetical format.
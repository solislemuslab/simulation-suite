# Simulation suite: compilation of scripts to simulate data on phylogenetic networks

I want to
- Simulate phylogenetic networks under a birth-death-hybridization model using `SiPhyNetworks` R package
    - Use `scripts/simulate-networks.R`
    - Dependencies:
        - R (version 4.1.2 or newer)
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
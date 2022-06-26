# R script to simulate N phylogenetic networks under a birth-death-hybridization
# model implemented in the R package of SiPhyNetwork
# Input parameters:
# - N: number of networks to simulate
# - lambda: speciation rate
# - mu: extinction rate
# - nu: hybridization rate
# - hybprobs: 3d vector with frequencies of each of the three types of
#             hybridizations (see https://github.com/jjustison/SiPhyNetwork/blob/master/vignettes/figures/hybridization_types.png)
# - ntips: number of tips in the network (there are also options to simulate networks of a given height)
# - seed: random seed for reproducibility
# Output:
# - N output files each with a simulated network in parenthetical format (note that branches on the network depend on the units of the rates)
# - txt file with parameter choices
# Claudia (June 2022)


library(SiPhyNetwork)

## Input parameters
N = 10
lambda = 0.9
mu = 0
nu = 0.02
hybprops = c(1,1,1)
ntips = 15
seed = 1234

## Where to save the network output files:
outfolder = "../data/"
setwd(outfolder)

## Name of output files:
rootnetworks = "network"
outfile = "parameters.txt"

set.seed(seed)

# From SiPhyNetworks docs: 
# The `sim.bdh.taxa.ssa` function uses the Simple Sampling Approach (SSA) to simulate to a specified number of taxa `n`, that is, the simulation will stop once the phylogeny has `n` tips.
# * `frac = 1` and `stochsampling = FALSE` assumes that all extant taxa are sampled in the phylogeny.
# * `mrca=TRUE` means that we start with a single lineage rather than the MRCA of two lineages 
# * `complete = TRUE` leaves extinct species on the phylogeny.
# * `hyb.rate.fxn = NULL` assumes that successful hybridization events are not a function of the genetic distance between taxa.
# * `trait.model = NULL` assumes that successful hybridization events do not depend on a trait value between taxa.
networks <- sim.bdh.taxa.ssa(n = ntips,
                             numbsim = N,
                             lambda = lambda,
                             mu = mu ,
                             nu = nu,
                             hybprops = hybprops,
                             hyb.inher.fxn = make.beta.draw(1, 1),
                             frac = 1,
                             mrca = FALSE,
                             complete = TRUE,
                             stochsampling = FALSE,
                             hyb.rate.fxn = NULL,
                             trait.model = NULL)

# GAB: Code for removing bad networks before writing
# get rid of null trees which go extinct=0:
networks <- networks[!sapply(X = networks, FUN = is.null)]
# and networks with no extinct tips are sampled=1:
#networks <- networks[sapply(X = networks, FUN = is.phylo)]
# fixit: we cannot use this bc is.phylo is in phylosim package which is not available for R 4.1.2


# Writing networks to files:
for(i in 1:length(networks)){
  SiPhyNetwork::write.net(net = networks[i], file = paste0(rootnetworks,"-",i,".net"))
}

# Writing parameters to reproduce to file:
txt = paste0("N = ", N, "\n lambda = ", lambda, "\n mu = ", mu, "\n nu = ", nu, "\n hybprops = ", hybprobs[1],",",hybprobs[2],"," ,hybprobs[3], "\n ntips = ", ntips, 
             "\n seed = ", seed)
writeLines(txt, outfile)

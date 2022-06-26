# Julia script that runs ms to simulate gene trees under
# a given species network (or tree).
#
# Input:
# - msnet: text file with the network in ms format, see ms docs: http://home.uchicago.edu/~rhudson1/source/mksamples.html
#          Branch lengths are assumed to be in coalescent units.
#          This file containts the -ej and -es events
# - net: text file with the network in parenthetical format.
#        We need this file simply for the number of taxa
# - numalleles: number of alleles per taxon (this script assumes all taxa have the same
#               number of alleles, but see notes below)
# - numgt: number of gene trees to simulate
#
# Output:
# - treefile: text file with simulated gene trees (one per line) in parenthetical format.
#             Note that ms includes the ms command in this file.
#
# NOTES:
# - The script does not check that both networks (parenthetical format and ms format) match!
# - The script also assumes that the ms executable is in your path
# - If different taxa should have different number of alleles, this needs to be manually 
#   modified in the alleles vector
# - We wanted to use the -seeds x1 x2 x3 option, but we always got an error. The seed is now saved in seedms.
# - The script needs to move into the folder with the input and output since ms saves seedms in the working directory.
# Claudia (June 2022)

using PhyloNetworks

# Input -----------------------
folder = "../data/"
msnet = "msnet1.txt"
net = "network-1.net"
numgt = 100
#seed1 = 123
#seed2 = 456
#seed3 = 789
numalleles = 1
treefile = "gene-trees-1.tre"
# -----------------------------

# Move into the folder
cd(folder)

# Reading the network on parenthetical format:
network = readTopology(net)
n = length(tipLabels(network))
alleles = numalleles * ones(Int,n)


# Writing the ms command (as Vector{String}):
mscommand = ["ms", "$n", "$numgt", "-T", "-I", "$n"]
for a in alleles
    push!(mscommand, "$a")
end

# Reading the network in ms format:
l = readlines(msnet)

# Adding it to the ms command:
for ll in split(l[1])
    push!(mscommand,ll)
end

# Running the ms command
out = read(Cmd(mscommand), String)

# Saving the gene trees
io = open(treefile,"w")
write(io,out)
close(io)


# for testing
#mscommand = ["ms", "6", "100", "-T", "-I", "6", "1", "1", "1", "1", "1", "1", "-ej", "0.1", "1", "2", "-ej", "0.5", "3", "4", "-es", "0.55", "4", "0.5", "-ej", "0.55", "7", "5", "-ej", "0.6", "2", "4", "-ej", "1.6", "4", "5", "-ej", "2.1", "5", "6"]
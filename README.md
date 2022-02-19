# MTBCPhyDyn
Phylodynamic simulations and inference of inter-species Mycobacterium spp transmission. 

determ_asymm.R
Deterministic approximations of tree asymmetry metrics for a two-host SI model, adopted from Frost & Volz (2013). Metrics modelled are: i) normalized number of cherries, ii) normalized Sackin's index, and iii) assortativity coefficient. 

stochasticsim.R
Simulating transmission in a two-host SI model using the adaptivetau R package. 

simtophylo.R
Generating transmission tree as a phylo object from simulated transmission history.

sampledPhylo.R
Given a specific sampling time, sample extant tips from whole transmission tree and shorten all terminal branch lengths to reflect time of sampling. This sampled and shortened tree is returned as a newick file and phylo object. 
Code for visualizing sampled and shortened tree is also available.
Given a transmission tree, calculate asymmetry metrics (Sackin's index and number of cherries).

seqGen.R
Given a transmission tree, simulate sequence evolution for its lineages according to a specified model of nucleotide substitution. 
Simulated sequences are aligned and returned as a fasta file.

trans_analysis.R
Analysing transmission trees from MASTER simulations, to i) generate a log file of all simulation events ii) identifying chains of transmission and iii) annotate transmission tree branches by their 'parent' infector.

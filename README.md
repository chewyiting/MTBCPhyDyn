# MTBCPhyDyn
Phylodynamic simulations and inference of inter-species Mycobacterium spp transmission. 

determ_asymm.R
Deterministic approximations of tree asymmetry metrics for a two-host SI model, adopted from Frost & Volz (2013). Metrics modelled are: i) normalized number of cherries, ii) normalized Sackin's index, and iii) assortativity coefficient. 

stochasticsim.R
Simulating transmission in a two-host SI model using the adaptivetau R package. 

simtophylo.R
Generating transmission tree as a phylo object from simulated transmission history.

seqgen_gen.sh
Generating XML files for SeqGen simulations of sequence evolution via BEAST from newick files.

trans_analysis.R
Analysing transmission trees from MASTER simulations, to i) generate a log file of all simulation events ii) identifying chains of transmission and iii) annotate transmission tree branches by their 'parent' infector.

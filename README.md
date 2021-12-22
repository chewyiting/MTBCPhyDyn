# MTBCPhyDyn
Phylodynamic simulations and inference of inter-species Mycobacterium spp transmission. 

trans_LHS.R
Exploring parameter space for the two-host SI model of bovine and zoonotic TB transmission outlined by Ellen Brooks-Pollock et al (2015).
To adapt this model in the context of endemic TB, equilibrium prevalence in both host species and turnover rates used are constrained to realistic values observed in high-burden TB countries experiencing both bovine and zoonotic TB.

trans_scenarios.sh 
Generating XML files for MASTER simulations via BEAST.

trans_analysis.R
Analysing transmission trees from MASTER simulations, to i) generate a log file of all simulation events ii) identifying chains of transmission and iii) annotate transmission tree branches by their 'parent' infector.

trans_phytrees.sh
Generating XML files for SeqGen simulations of sequence evolution via BEAST


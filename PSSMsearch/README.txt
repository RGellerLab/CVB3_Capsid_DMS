This folder contains the data and scripts for performing the PSSMsearch analysis based on the AA preferences of the 3CDpro cleavage sites
1. script_write_peptides.R: reads in data of AA preferences (p1_mean_prefs.csv) and writes peptides that correspond to the frequency of each amino acid at each position
2. peptides folder: contains the results of the R script_write_peptides
3. pssm_results.tdt: the output of PSSMsearch when querying with peptide sequences from both cleavage sites together
4. PSSM_from_PSSMsearch.csv: the PSSM calculated by PSSMsearch

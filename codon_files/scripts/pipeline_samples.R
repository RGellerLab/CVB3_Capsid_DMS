## this script will take the virvarseq output and generate filtered codon tables
## using the virvarseq_to_codon_table function
## in addition, it can filter single mutants in codons by running the 
## filter_single_codons

## requirements library(Biostrings)

## load functions
source("./functions_virvar_analysis_cod_table.R")


## get codon tables form virvarseq output
data.dir= "../virvarseq/"# input where the codon files are found
fs=list.files(pattern = "*.codon$",data.dir,include.dirs = F)

## run script on all files
## produces folder with codon_tables (after filtering by SOR if selected)
## and folders with the filtered tables and rejected mutations
for (f in fs){ 
  virvarseq_to_codon_table(file = paste0(data.dir,f) ,
                           name = unlist(strsplit (f,split = "d3"))[1], 
                           out.directory="./processed_tables/",
                           strand.test="SOR", SOR.cut=4,
                           ref="./p1_nuc.fasta",
                           f.length=851);
  rm (f)}


rm(fs)


##########################
###FILTER SINGLE MUTS

filter_single_codons(min.cov=0,
                     codon.identity.matrix="./codon_identity_matrix.csv", 
                     data.dir="./processed_tables/codon_tables/",
                     out.dir="./filter_single_mutants/")

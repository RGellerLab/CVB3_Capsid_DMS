#cleavage sites
library(tidyverse)
##### 
df=read.csv("./p1_mean_prefs.csv",stringsAsFactors = F)
df=df %>% select(-(site:wildtype))


#############
library(Biostrings)
wt=readAAStringSet("./p1_aa.fasta")
wt=unname(unlist(strsplit(as.character(wt),"")))


### set width of window around cut site
length_window=5
sites=length_window-1
vp23=df[(332-sites):(333+sites),]
vp31=df[(570-sites):(571+sites),]

########## write peptides
get_peptides=function(df,pept_num=1000, col.rm=F, col.num=1, wt.seq){
  ## data frame should be rows as position in peptide and columns as AA.
  if (col.rm==T){    df=df[,-col.num]}
  ## gets number of peptides to make for each one
  test=round(df*pept_num,digits = 0)
  res=NULL
  ## go by row and then by column and make peptides
  for (r in 1:nrow(test)){
    # r=4 # go by row
    for (c in 1:ncol(test)){ ## go by aa which is column
      # c=17
      aa=colnames(test)[c] # get new aa to put in
      temp.seq2=wt.seq # take sequence to mutate
      temp.seq2[r]=aa # replace the position with the new aa
      temp.seq2=paste0(temp.seq2,collapse = "") # collapse into string
      if (is.na(test[r,c])==0){
      res=c(res, rep(temp.seq2,test[r,c])) # add to list after making x many copies as required
    }}
  }
return(res)
}

## get peptides

vp23_pept=get_peptides(vp23, wt.seq=wt[(332-sites):(333+sites)])
vp31_pept=get_peptides(vp31, wt.seq=wt[(570-sites):(571+sites)])

write(vp23_pept,"./vp23_peptides.txt")
write(vp31_pept,"./vp31_peptides.txt")
##########

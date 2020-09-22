virvarseq_to_codon_table=function(file, name="l1", f.length=851, strand.test=c("fisher","SOR","none"),
                                  SOR.cut=4,
                                  out.directory= "", ref="../p1_nuc.fasta"){
  # takes vir var seq output file, filters it, removes results that fail strand bias, and writes both
  # cleaned file and a codon table of the results)
  
  ############### PARAMETERS
  # name="l1" ## sample name
  # out.directory= ## directory of where to write files
  # f.length=851 ## length of AA sequence 
  # file="../july19/pld2c0.7t10.codon" ## vir var seq ouput file
  # strand.test=c("fisher","SOR","none") ## which test wanted
  # SOR.cut=4 # where do you want to filter strand bias
  
  library(Biostrings)
  
  ### read wt sequence and get codons
  refseq=as.character(readDNAMultipleAlignment(ref))
  refseq=strsplit(refseq,"")[[1]]
  ref.codons=sapply(seq(from = 1,to = (length(refseq)-2), by = 3), function(x) paste0(refseq[x:(x+2)],collapse = ""))
  rm(ref,refseq)
  
  ## get genetic code
  gen_code=as.data.frame(cbind(codons=as.character(names(GENETIC_CODE)),AA=GENETIC_CODE),stringsAsFactors=F)
  gen_code=gen_code[order(gen_code$codons),]
  
  ##make results df
  res.df=data.frame(matrix(data=0, nrow=f.length,ncol=66))
  colnames(res.df)=c("site","wildtype",gen_code$codons)
  res.df$site=1:f.length
  rm(f.length)
  
  #### read vir var seq data
  vvs=read.delim(file,stringsAsFactors = F) 
  vvs$REF_CODON=toupper(vvs$REF_CODON)
  
  # remove amino acids that are X
  vvs=vvs[vvs$CODON%in%gen_code$codons,]
  
  ######## remove failure from strand bias, need to recalculate denom after removing X
  positions=unique(vvs$POSITION)
  
  print("fixing counts")
  
  #### redo denom
  for (i in positions){ 
    print(i)
    # i=831############
    ## fix up denominator
    vvs[vvs$POSITION==i,"FWD_DENOM"]=sum(vvs[vvs$POSITION==i,"FWD_CNT"])
    vvs[vvs$POSITION==i,"REV_DENOM"]=sum(vvs[vvs$POSITION==i,"REV_CNT"])
    vvs[vvs$POSITION==i,"DENOM"]=sum(vvs[vvs$POSITION==i,"CNT"])
    rm(i)
    
  }
  rm(positions)  
  
  
  ### function for SOR from GATK documentation, StrandOddsRatio
  SOR=function(wtF, wtR, mutF, mutR){
    # Add pseudocount
    wtF=wtF+1;
    wtR=wtR+1
    mutF=mutF+1
    mutR=mutR+1
    r= (wtF*mutR)/(wtR*mutF)
    R= r+1/r
    refRatio=(min(wtF,wtR)/ max(wtF,wtR))
    mutRatio=(min(mutF,mutR)/ max(mutF,mutR))
    log(R)+log(refRatio)- log(mutRatio)
  }
  
  
  #### do strand test for strand bias
  positions=unique(vvs$POSITION)
  vvs$std.bias=1
  vvs$id=1:nrow(vvs) #give an ID for each row 
  
  print("fixing strand bias")
  for (i in positions){
    print(i)
    vvs2=vvs[vvs$POSITION==i,] ## vvs of individual AA
    
    if(nrow(vvs2)>1){  ### if there are mutants to test get WT and mutant rows
      wt=vvs2[vvs2$CNT==max(vvs2$CNT),] ## WT is defined as most common
      muts=vvs2[vvs2$CNT!=max(vvs2$CNT),]  
      
      if (strand.test=="fisher"){
        ### iterate overall mutants to test vs wt
        for (j in 1:nrow(muts)){
          # j=1
          vvs$std.bias[vvs$id==muts$id[j]]=fisher.test(rbind(c(wt$FWD_CNT, wt$REV_CNT),  #### pvalue calculation for bias
                                                             c(muts$FWD_CNT[j], muts$REV_CNT[j])))$p.value
          rm(j)
        }
        rm(i,vvs2,wt,muts)
      }
      
      if (strand.test=="SOR"){
        ### iterate overall mutants to test vs wt
        for (j in 1:nrow(muts)){
          # j=1
          vvs$std.bias[vvs$id==muts$id[j]]=SOR(wt$FWD_CNT, wt$REV_CNT,muts$FWD_CNT[j], muts$REV_CNT[j])
          rm(j)
        }
        rm(i,vvs2,wt,muts)
      }
      
      if (strand.test=="none"){print("not fixing strand")
      }
      
    }
  }
  
  rm(positions)
  
  ## correct frequency
  vvs$FREQ=vvs$CNT/vvs$DENOM
  vvs=vvs[,-ncol(vvs)] ## remove ID column
  
  ## strand.test 
  
  if (strand.test==T){
    p.val_strand_bias=0.05/nrow(vvs)
    rejected=vvs[vvs$std.bias<=p.val_strand_bias,] # those that would get kicked out normally
    vvs=vvs[vvs$std.bias>p.val_strand_bias,] #strand bias if chosen
    rm(p.val_strand_bias)
  } 
  
  ## fdr correction
  if (strand.test=="SOR"){
    rejected=vvs[vvs$std.bias>=SOR.cut,] # those that would get kicked out normally
    vvs=vvs[vvs$std.bias<SOR.cut,] #strand bias if chosen
  } 
  
  
  
  # iterate over rows of data and put in for each position/codon put in count data
  for (r in 1:nrow(vvs)){res.df[vvs$POSITION[r],vvs$CODON[r]]=vvs$CNT[r]; rm(r)}
  
  ## add wildtype to codon file
  res.df$wildtype=ref.codons
  fix.mut=sum(ref.codons!=apply(res.df[3:ncol(res.df)], 1 , function (x) colnames(res.df[3:ncol(res.df)])[which(x== max(x))]))
  if (fix.mut>1){write(name,paste0(out.directory, "fixed_mutations.txt"),row.names = F,append = T)}
  
  #### write results     ##########
  ### create directories
  dir.create(  out.directory, recursive = T)
  dir.create(  paste0(out.directory,"/codon_tables/"), recursive = T)
  dir.create(  paste0(out.directory,"/rejected_codon_tables/"), recursive = T)
  dir.create(  paste0(out.directory,"/filtered_tables/"), recursive = T)
  ### write files
  write.csv(res.df,paste0(out.directory,"/codon_tables/",name, "_codon_table.csv"),row.names = F)
  write.csv(rejected,paste0(out.directory,"/rejected_codon_tables/",name, "_rej_codon_table.csv"),row.names = F)
  write.csv(vvs,paste0(out.directory,"/filtered_tables/",name, "_filtered.csv"),row.names = F)
  write(c(SOR.cut,strand.test),paste0(out.directory,"/codon_tables/", "strandtest.txt"))
}


filter_single_codons=function(min.cov=0,
                              codon.identity.matrix="./codon_identity_matrix.csv", 
                              data.dir="./codon_tables/",
                              out.dir="./filter_single_mutants/"){
  
  ## filter codon table for more than 1 mutation, uses output 
  ## from virvarseq_to_codon_table function
  
  # PARAMETERS##
  # min.cov=0 # minimum coverage for accepting mutation, if below put 0
  # codon.identity.matrix="./codon_identity_matrix.csv"
  # data.dir="./" # where the data for the codon files is located
  # out.dir="./filter_single_mutants/") # where you want to write the files
    
  ## load data on which codons are a single mutation away
  cod.tbl=read.csv(codon.identity.matrix,stringsAsFactors = F,row.names = 1)

  ## create output director
  dir.create(out.dir,recursive = T)
  ## iterate on all codon files in data directory
  fs=list.files(data.dir,pattern = "*.csv",full.names = T)
  for (f in fs){
    df=read.csv(f,stringsAsFactors = F)
    ## iterate on row, get columns that are single mutants and change to 0
    for( x in 1:nrow(df)){
      df[x, colnames(cod.tbl)[which(cod.tbl[df$wildtype[x],]==1)]]=0
      rm(x)
    }
    ### put at 0 all below threshold
    df2=df[,3:ncol(df)]
    df2[df2<min.cov]=0
    df[,3:ncol(df)]=df2
    rm(df2)
    
    # write results, name obtained by taking first element split by _ from
    # codon file name
    write.csv(x = df, file = paste0(out.dir,unlist(strsplit(basename(f),"_"))[1],
                                    "_2_3mut_codon_table.csv"),row.names = F)
    rm(f,df)
  }
}




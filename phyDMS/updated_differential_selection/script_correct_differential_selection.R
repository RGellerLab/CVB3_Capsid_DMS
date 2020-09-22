library(tidyverse)
## dnds
phy=read_delim("../program_output/out_ExpCM_p1_mean_prefs_diffprefsbysite.txt",delim = "\t",comment = "#")
phy=phy %>% arrange(site)
colnames(phy)=gsub("dpi_","",colnames(phy))


### read positions to remove single mutants in codons
aa.rm=read_csv("./p1_aa_positions_1away.csv")

for (i in 1:nrow(phy)){
  # i=1
  site.i=aa.rm$site[i]
  row.i=aa.rm[i,-c(1:2)]
  phy[phy$site==site.i,colnames(row.i)[row.i==1]]=NA
  rm(site.i,row.i,i)

}
rm(aa.rm)
### update sum of absoulte
phy$sum_absolute_differential_selection=apply(select(phy,-site,-half_sum_abs_dpi),1,function(x)  (sum(abs(x),na.rm = T)))

phy=phy %>%  select(-half_sum_abs_dpi)

write.csv(phy, ".abs_differential_selection.csv", row.names = F)

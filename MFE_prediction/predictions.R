library(tidyverse)
library(hexbin)
library(RColorBrewer)
library(randomForest)



## read data
data_fac=read_csv("./data_predictions.csv")%>% mutate_if(is.character, as.factor)

#### linear model


linmod=(lm(data = data_fac, p1_avg_log2effect~
             entB+
             wt*
             mutation+
             pent_ddg+
             RSApent,na.action = "na.omit"
))

cor.test(data_fac$p1_avg_log2effect,predict.lm(linmod,data_fac),method="pearson") ## 0.67
cor.test(data_fac$p1_avg_log2effect,predict.lm(linmod,data_fac),method="spearman") ## 0.67


###### graph###############################
## colors
hmcol<-rev(brewer.pal(9,"YlGnBu"))#RdBu and 11 nice
# predicted values
predicted.lm=predict.lm(linmod,data_fac)
## make data frame for plotting
df.lm=data.frame(MFE=data_fac$p1_avg_log2effect,predicted.MFE=predicted.lm)

ggplot(df.lm,aes(x=MFE,y=predicted.MFE))+
  geom_point(size=1)+
  scale_y_continuous(limits=c(-8,4),expand = c(0,0))+
  scale_x_continuous(limits=c(-8,4),expand = c(0,0))+
  xlab("MFE")+ylab("Predicted MFE")+
  geom_hex(bins = 80)+
  scale_fill_gradientn(limits=c(0,50), breaks=seq(0, 50, by=10), colours=hmcol)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10, family="Arial", color="black"),
        axis.title=element_text(size=10, color="black"))
rm(linmod,predicted.lm,df.lm, hmcol, data_fac)

#############################################################################
########### Random forest########################################

## split into train and test data sets
train=data_fac %>% filter(type=="train") %>% select(-position, -type)
test=data_fac %>% filter(type=="test") %>% select(-position, -type)
rm(data_fac)


# to generate the model uncheck and run the code below
#modelRF=randomForest(p1_avg_log2effect~.,data=train, importance = TRUE,mtry=10)#54.84

## load model 
load("./modelRF.RData")

### correlation between model and observed
cor.test(train$p1_avg_log2effect,modelRF$predicted,method="pearson")#0.75
cor.test(train$p1_avg_log2effect,modelRF$predicted,method="spearman")#0.74


## get variable importance
import=data.frame(importance(modelRF)) %>% arrange(desc(X.IncMSE))
import$parameters=factor(row.names(import), levels=row.names(import))


## plot var importance##############

## full data set
ggplot(import) + geom_col(aes(y = parameters,x=X.IncMSE))+
  xlab("% Increase MSE")+ylab("Parameters")+ggtitle("")+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.text=element_text(size=10),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=10, family="Arial", color="black"),
    axis.title=element_text(size=10, color="black"))
 


## top 10
import=import[order(import[,1],decreasing = T),]
graph.imp=import[1:10,] %>% arrange(X.IncMSE)
graph.imp$attribute=row.names(graph.imp)
graph.imp$attribute=factor(graph.imp$attribute,level=graph.imp$attribute)
library(RColorBrewer)
ggplot(graph.imp,aes(y=attribute,x=X.IncMSE,fill=brewer.pal(n = 10, name = "Set3")))+geom_bar(stat="identity")+
  xlab("% increase MSE")+ylab("")+
  scale_x_continuous(limits=c(0,80),expand = c(0,0))+
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.text=element_text(size=10),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=10, family="Arial", color="black"),
    axis.title=element_text(size=10, color="black"))
rm(graph.imp)
######################################


###### predictions:
pred.test=predict(modelRF,test)
cor.test(test$p1_avg_log2effect,pred.test,method="pearson")### 0.76
cor.test(test$p1_avg_log2effect,pred.test,method="spearman")### 0.75


# colors
hmcol<-rev(brewer.pal(9,"YlGnBu"))
# df for graphing
test.df=bind_cols(MFE=test$p1_avg_log2effect,predicted.MFE=pred.test)

ggplot(test.df,aes(x=MFE,y=predicted.MFE))+
  geom_point(size=1)+
  scale_y_continuous(limits=c(-7.5,2.5),expand = c(0,0))+
  scale_x_continuous(limits=c(-7.5,2.5),expand = c(0,0))+
  geom_hex(bins = 50)+
  xlab("MFE")+ylab("Predicted MFE")+ggtitle("Random forest model 54 predictors")+
  scale_fill_gradientn(limits=c(0,30), breaks=seq(0, 30, by=10), colours=hmcol)+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=10, family="Arial", color="black"),
    axis.title=element_text(size=10, color="black"))

#############################################

### select top 5 predictors
updated.fact=c("p1_avg_log2effect",rownames(import[1:5,]))
model5=randomForest(p1_avg_log2effect~.,data=train[,updated.fact])#50.4
pred.test=predict(model5,test)

cor.test(test$p1_avg_log2effect,pred.test,method="pearson")
cor.test(test$p1_avg_log2effect,pred.test,method="spearman")
save(model5,file = "scripts_github/model5.RData")

## graph
test.df=bind_cols(MFE=test$p1_avg_log2effect,predicted.MFE=pred.test)

hmcol<-rev(brewer.pal(9,"YlGnBu"))

ggplot(test.df,aes(x=MFE,y=predicted.MFE))+
  geom_point(size=1)+
  scale_y_continuous(limits=c(-7.5,2.5),expand = c(0,0))+
  scale_x_continuous(limits=c(-7.5,2.5),expand = c(0,0))+
  xlab("MFE")+ylab("Predicted MFE")+
  geom_hex(bins = 50)+
  xlab("MFE")+ylab("Predicted MFE")+ggtitle("Random forest model top 5 predictors")+
  scale_fill_gradientn(limits=c(0,30), breaks=seq(0, 30, by=10), colours=hmcol)+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text=element_text(size=10, family="Arial", color="black"),
    axis.title=element_text(size=10, color="black"))



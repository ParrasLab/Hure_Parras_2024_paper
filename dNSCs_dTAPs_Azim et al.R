# Auteur : Rihab AZMANI 
# script Azim data
setwd("/home/orlab/Desktop/Rihab/Dataset Azim et al") # pathway of dataset 

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

# install.packages('class', dependencies=TRUE)
# install.packages('codetools', dependencies=TRUE)
# install.packages('affy', dependencies=TRUE)
# install.packages('BiocGenerics', dependencies=TRUE)
# install.packages('rma', dependencies=TRUE)
library(class)
library(codetools)
library(affy)

# importer et visualiser les donnees avant normalisation 
Mydata <-ReadAffy() 
dim(Mydata)
boxplot(Mydata) 
hist(Mydata)

## normalisation intra-lames 

# normalize.methods(Mydata)
# MydataIntra=normalize(Mydata, method='quantiles')
# # write.csv(Mydata2, file="dSVZvs_lSVZ_quantile2.txt")
# hist(MydataIntra)
# boxplot(MydataIntra)
# normalisation inter-lames 

## normalisation avec le package rma 

library(gcrma)  
Mydata<-ReadAffy() 
eset2 <- rma(Mydata) 
# #(choisir rma ou gcrma )
write.exprs(eset2, file="dSVZvs_lSVZ_rma.xlsx", sep="\t") 
MyData=read.table("dSVZvs_lSVZ_rma.xlsx", header=TRUE) 
dim(MyData) # 45101 23


library(stringr)
tmp =NULL
tabtemp = NULL
for(i in colnames(MyData)){
  tabtemp = rbind(tabtemp,str_split_fixed(i,"_",5)[2:4])}

tabtemp[which(tabtemp[,2] =="Dorsal"),2]<-"d"
tabtemp[which(tabtemp[,2] =="Lateral"),2]<-"l"

tmpnames = NULL

for(u in 1:dim(tabtemp)[1]){
  tmpnames = c(tmpnames,paste0(tabtemp[u,1],"_",tabtemp[u,2],tabtemp[u,3]))
}

colnames(MyData) <- tmpnames

# ACP 
## make CPA with ExPosition package 
# install.packages('prettyGraphs', dependencies=TRUE)
# install.packages('ExPosition', dependencies=TRUE)

library('prettyGraphs')
library('ExPosition')
MydataQ=apply(MyData, 1, IQR)
##  barplot to visualize the quantiles from 0 to 1 by 10%
barplot(quantile(MydataQ, probs=seq(0,1,.01)))
# 
##  Filter the data for > or = to .4 (40%) quantiles
MyFiltData50=MyData[MydataQ >= quantile(MydataQ, .4), ]
write.csv(MyFiltData50, file="dSVZvs_lSVZ_quantile1.txt")  # 27061
MyEpPCA = epPCA(t(MyFiltData50), scale = TRUE, center = TRUE, DESIGN = NULL,
                make_design_nominal = TRUE, graphs = T, k = 0)
dim(MyFiltData50)
write.table(MyEpPCA$ExPosition.Data$fj, file='MyEpPCA_filt40_IDs_FsR.txt', sep="\t",
            quote=F, col.names=F, row.names=F)
DataFit=read.table("MyEpPCA_filt40_IDs_FsR.txt", sep = "\t", header=TRUE, row.names=F) # 
dim(DataFit)
Fs2=read.table('MyEpPCA_filt40_IDs_FsR.txt', sep = "\t", header=TRUE, row.names=F)
dim(Fs2)
colnames (Fs2) = c (paste0('F',1:22))  ### to add names to the columns
q=sapply(Fs2, function(x) quantile(x, probs=seq(from=0, to=1, by=.01), type = 1))
barplot (q)
sapply(as.data.frame(q), barplot)
sapply(1:22, function(i) barplot(q[,i]))
q[paste0(c(10,90), '%'),] #10%s extreme values
F2tens = subset(Fs2, F2 >=1.4 | F2 <=-2.3) #subset for the 10%s extreme values of F2
dim (F2tens) #5184   22

## make CPA with Ade4 package 

# install.packages(c('ade4', 'factoextra', 'ggplot2'),
# repos = c("http://rstudio.org/_packages",
# "http://cran.rstudio.com"))
library('ade4')
library('ggplot2')
library('factoextra')
MydataQ=apply(MyData, 1, IQR)
##  barplot to visualize the quantiles from 0 to 1 by 10%
barplot(quantile(MydataQ, probs=seq(0,1,.01)))
MyFiltData50=MyData[MydataQ >= quantile(MydataQ, .4), ]
write.csv(MyFiltData50, file="dSVZvs_lSVZ_quantile2.txt")  # 27061
MyadPCA = dudi.pca(t(MyFiltData50), scale = TRUE , scan = FALSE , nf=10)
fviz_eig(MyadPCA)
fviz_pca_ind(MyadPCA,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB"),
             repel = TRUE)

##ID to gene names
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
# BiocManager::install("annotate")
# BiocManager::install("mouse4302.db", dependencies=TRUE)

require(org.Mm.eg.db)
library("annotate")
library("mouse4302.db")

x=keys(mouse4302.db, keytype="PROBEID")
addNames = select(mouse4302.db,keys = x , c("SYMBOL","GENENAME","ENTREZID"))  # 47653 
dim(addNames)
removeDuplicat=split(addNames, addNames$PROBEID) 
n=sapply( removeDuplicat, nrow) 
removeDuplicat.uniq=do.call(rbind, removeDuplicat[n==1]) 
removeDuplicat.duplic=do.call(rbind, removeDuplicat[n!=1]) 
lapply(removeDuplicat.uniq, nrow) # see the number of rows in yy.uniq objects
lapply(removeDuplicat.duplic, nrow) 
removeDuplicat.all = (rbind(removeDuplicat.uniq, removeDuplicat.duplic)) 
dim(removeDuplicat.all)
## to get the "SYMBOL","GENENAME" of a list of probesets

ids <- rownames(F2tens)
IDS=intersect(ids, removeDuplicat.uniq$PROBEID)
F2tens=F2tens[IDS, ]
IDS.uniq=removeDuplicat.uniq[match(IDS, removeDuplicat.uniq$PROBEID), ]
all(IDS.uniq$PROBEID == IDS) 
F2tens.ids=cbind(F2tens, subset(IDS.uniq, select=-PROBEID))
dim(F2tens.ids)
head(F2tens.ids)
write.table(Fs2.ids, file='MyEpPCA_F2_10p_IDsR.txt', sep="\t", quote=F, col.names=T, row.names=F)


## ANALYSE Limmma

library(limma)
library(reshape)

dim (MyData) # 45101   23
summary(MyData)
head (MyData)
tail(MyData)
test=MyData 

### comparaison 

test=MyData
test.group=factor(with(melt(sapply(c('dS','dN','dT', 'lS', 'lN', 'lT'), function(x) grep(x, names(test), value=FALSE),
                                   simplify=FALSE)), L1[order(value)]), levels=c('dS','dN','dT', 'lS', 'lN', 'lT'), labels=paste0('Group', 1:6))
#design <- model.matrix(~test.group) ## TO BE CHECKED == in this case WT is the reference group
design <- model.matrix(~test.group-1) ## TO BE CHECKED AGAIN == no reference here
colnames(design)=levels(test.group)
contrast.matrix=makeContrasts(contrasts=c('(Group2+Group3)-(Group5+Group6)'), levels=design)
fit <- lmFit(MyData, design)
fit2c <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2c)
dim (fit2) # 45101    3
tt0=topTable(fit2, number=nrow(MyData), coef=NULL) 
i = match(rownames(MyData), rownames(tt0)) 
o = cbind(Id=rownames(MyData), MyData, tt0[i, ]) 
write.table(x=o, file='dNSCsTAPs_vs_lNSCsTAPs_rma_limma.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.')

new0=tt0
names(new0)=paste0('NSCsTAPs', names(tt0))
p=cbind(Id=rownames(MyData), MyData, new0[i, ])
write.table(x=p, file='dNSCsTAPs_vs_lNSCsTAPs_rma_limma.txt', row.names=F, col.names=T, quote=F, sep="\t", dec='.')

ids= rownames(p)
IDS=intersect(ids, removeDuplicat.uniq$PROBEID)
newp=p[IDS, ]
IDS.uniq=removeDuplicat.uniq[match(IDS, removeDuplicat.uniq$PROBEID), ]
#In order to check ids matches one can do:
all(IDS.uniq$PROBEID == IDS)    #should return TRUE
p_ids=cbind(newp, subset(IDS.uniq, select=-PROBEID))
dim (p_ids) 
head (p_ids)
write.table(p_ids, file='dNSCsTAPs_vs_lNSCsTAPs_rma_limma_IDs.txt', sep="\t",
            quote=F, col.names=T, row.names=F)

IDS=intersect(ids, removeDuplicat.uniq$PROBEID)
newp=p[IDS, ]
IDS.uniq=removeDuplicat.uniq[match(IDS, removeDuplicat.uniq$PROBEID), ]
#In order to check ids matches one can do:
all(IDS.uniq$PROBEID == IDS)    #should return TRUE
pp_ids=cbind(newp, subset(IDS.uniq, select=-PROBEID))
dim (pp_ids) 
head (pp_ids)
write.table(pp_ids, file='dNSCsTAPs_vs_lNSCsTAPs_rma_limma_IDs.txt', sep="\t",
            quote=F, col.names=T, row.names=F)

Data2=read.table("dNSCsTAPs_vs_lNSCsTAPs_rma_limma_IDs.txt", header=TRUE , sep="\t")

#install.packages(c("tidyverse" ,"dplyr") , dependencies = T)

library(tidyverse)
library(dplyr)
FC_Max= ifelse((Data2$NSCsTAPslogFC < 0) ,(-(0.5)^(Data2$NSCsTAPslogFC)) , (2^(Data2$NSCsTAPslogFC)))
FCRealMax=round(FC_Max,2)
Data2<- cbind(Data2 , FCRealMax)
up2 <- filter(Data2 , Data2$FCRealMax >= 1.2 & Data2$NSCsTAPsP.Value < 0.05)
x=unique(up2$SYMBOL)
dim(up2) # 3143
head(up2)


library(stats)
library(ggplot2)
library(RColorBrewer)
donnes <- (up2[c(1,4,5,6,8,7,9,12,13,14,15,16)])
coul = colorRampPalette(brewer.pal(10, "PiYG"))(25)
# coul=sort(c("red", "green") ,decreasing = TRUE)
heatmap(as.matrix(donnes[,-1]) , col = coul ,  Colv = NA, Rowv = NA  ,  main ="heatmap up")
write.table(up2, file='RdNSCsTAPs_vs_lNSCsTAPs_rma_limma_IDs_up.txt', sep="\t",
            quote=F, col.names=T, row.names=F)
down <- filter(Data2 , Data2$FCRealMax <= -1.2 & Data2$NSCsTAPsP.Value < 0.05)
dim(down) # 2858
write.table(down, file='RdNSCsTAPs_vs_lNSCsTAPs_rma_limma_IDs_down.txt', sep="\t",
            quote=F, col.names=T, row.names=F)
library(stats)
library(ggplot2)
library(RColorBrewer)
donnes <- (down[c(1,4,5,6,8,7,9,12,13,14,15,16)])
# coul=sort(c("red", "green") ,decreasing = TRUE)
coul = colorRampPalette(brewer.pal(10, "PiYG"))(25)
heatmap(as.matrix(donnes[,-1]) , col = coul ,  Colv = NA, Rowv = NA  ,  main ="heatmap down")

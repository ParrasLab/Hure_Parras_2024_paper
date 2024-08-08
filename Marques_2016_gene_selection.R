
setwd("...")

library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(data.table)
library(stringr) 

# Load the dataset

data <-read.table("GSE75330_Marques_et_al_mol_counts2.txt", 
                  sep="", header = TRUE, row.names = 1, quote="")
dim(data) #23556  5069


# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes

mOLglia <- CreateSeuratObject(counts = data, min.cells = 10, min.features = 200, project = "Marques")
head(mOLglia)

# store mitochondrial percentage in object meta data    
mOLglia <- PercentageFeatureSet(mOLglia, pattern = "^MT-", col.name = "percent.mt")

#Apply sctransform normalization
mOLglia <- SCTransform(mOLglia, vars.to.regress = "percent.mt", verbose = FALSE)

#Run PCA
mOLglia <- RunPCA(mOLglia, verbose = FALSE)

#DETERMINE THE DIMENSIONALITY OF THE DATASET

ElbowPlot(mOLglia) #here, one can observe an 'elbow' around PC16-17, suggesting that the majority of true signal is captured in the first 10 PCs.
ggsave("mOLglia_ElbowPlot_PCs.pdf", width = 10, height = 6.00)



#Find Neighbours
mOLglia <- FindNeighbors(mOLglia, dims = 1:20, verbose = FALSE)


mOLglia <- FindClusters(mOLglia, resolution = 0.9, verbose = FALSE)
#it produces 18 clusters

#Run Non-linear dimensional reduction (tSNE)
mOLglia <- RunTSNE(object = mOLglia, dims.use = 1:20, do.fast = TRUE, perplexity = 80)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = mOLglia)

mOLglia.markers <- FindAllMarkers(object = mOLglia, only.pos = TRUE)
head(mOLglia.markers)

dim(mOLglia.markers) # 21852     7
length(unique(mOLglia.markers$gene))#6502
table(mOLglia.markers$cluster)
#  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 536  801 1025  707 1179  611 1469 1116 2083 1319 2655 1810  942 1079  923  665 2570  362 
head(mOLglia.markers)

top150 <- mOLglia.markers %>% group_by(cluster) %>% top_n(150, avg_log2FC)

dim(top150) #2701    7
length(unique(top150$gene)) #1615
selected_genes <- unique(top150$gene)


##################################################################################################
#Selecting genes from Marques 2016 supplemental table ‘aaf6463 Table S1’, sheet ‘specific genes’ 
#pooling all genes except those of VLMC column
##################################################################################################
 
dat1 <- c("Plp1", "Mbp", "Cnp", "Mal", "Enpp2", "Ugt8a", "Mog", "Cldn11", "Mobp", "Tspan2", "Mag", "Trf", "Ermn", "Car2", "Cryab", "Cmtm5", "Dbndd2", "Pllp", "Opalin", "Tmem88b", "Gsn", "Gpr37", "Sept4", "Slc12a2", "Qdpr", "Sirt2", "Tmeff2", "Olig1", "Plekhb1", "Grb14", "Tubb4a", "Gatm", "Cers2", "Psat1", "Csrp1", "Aspa", "Phgdh", "Fa2h", "Cyp51", "Elovl7", "2810468N07Rik", "Lgi3", "Pdlim2", "Slc44a1", "Gjc3", "Phldb1", "Cd9", "S100b", "Ddr1", "Omg")
dat2 <- c("Ptprz1", "Pdgfra", "Serpine2", "Cspg5", "Vcan", "Cspg4", "Atp1a2", "3110035E14Rik", "Fabp7", "Ednrb", "Sdc3", "Lhfpl3", "Bcan", "Ccnd1", "Ntm", "Sox11", "Zfp36l1", "Neu4", "Pcdh15", "Gpr17", "Kcnip3", "Gpr37l1", "Matn4", "Tmem176b", "Slc1a2", "Apoe", "Rgcc", "Tmem100", "Sox6", "Spon1", "Id2", "Gria3", "Traf4", "Rlbp1", "Ncald", "Ncan", "Slc1a3", "Tnr", "Cdo1", "C1ql1", "Midn", "S100a1", "Tpm1", "Tril", "Nnat", "Alcam", "Grin3a", "Fxyd6", "Gm2a", "Mtss1l")
dat3 <- c("Plp1", "Trf", "Mbp", "Mal", "Cnp", "Cldn11", "Mog", "Mobp", "Mag", "Ugt8a", "Apod", "Tspan2", "Ermn", "Cryab", "Car2", "Enpp2", "Gsn", "42617", "Dbndd2", "Opalin", "Cmtm5", "Tmem88b", "Gpr37", "Qdpr", "Cers2", "Pdlim2", "Slc12a2", "Grb14", "Aspa", "Csrp1", "Tmeff2", "Ptgds", "Gatm", "Sirt2", "Cyp51", "Elovl7", "Ndrg1", "Fa2h", "Epb4.1l3", "Tubb4a", "Efnb3", "Lpar1", "Lgi3", "Kcna1", "Psat1", "Plekhb1", "Arpc1b", "Tmem151a", "Pcyt2", "Phgdh")
dat4 <- c("9630013A20Rik", "Itpr2", "Gpr17", "Bcan", "Fyn", "Slc1a1", "Frmd4a", "Tmem2", "Mpzl1", "Sirt2", "Sema5a", "Col9a3", "Rap2a", "Idh1", "Tcf7l2", "Kank1", "Tns3", "Gpr37l1", "Enpp6", "1810041L15Rik", "Rras2", "Scrg1", "Luzp2", "Cnksr3", "Pik3r3", "Tmem163", "Rnf122", "Trio", "Mfsd2a", "Rims2", "Peli1", "Prom1", "Dock9", "Klhl5", "Rasgef1b", "Bcas1", "Tnr", "Parvb", "Rictor", "Bfsp2", "Csad", "Epb4.1l2", "Sema6a", "Gnb4", "Vcan", "Neu4", "Atp6v0a2", "Pcdh7", "Capn5", "Etv6")
dat5 <- c("Mal", "Apod", "Trf", "Car2", "Ermn", "Mog", "Gsn", "Aspa", "Opalin", "Ndrg1", "Pdlim2", "Cryab", "Sepp1", "Cntn2", "Sept4", "Mobp", "Glul", "Grb14", "Phgdh", "Omg", "Pla2g16", "Qdpr", "Litaf", "Plp1", "Fnbp1", "Tppp3", "Epb4.1l3", "Gpr37", "Evi2a-evi2b", "Kcna1", "Gng11", "Edil3", "Anln", "Ptgds", "Cmtm5", "Nrbp2", "Sez6l2", "Desi1", "Endod1", "Hapln2", "Ttll7", "Fez1", "Enpp2", "Tmem88b", "Nmral1", "Csrp1", "Fa2h", "Plekhh1", "Gjb1", "Cers2")
dat6 <- c("Cd9", "Neu4", "3110035E14Rik", "Bmp4", "Gpr17", "Vcan", "Bcas1", "AI414108", "S100a1", "Cyp2j6", "S100a13", "Sox6", "Epb4.1l2", "Ptprz1", "Sgk1", "Tubb2b", "Fyn",  "S100b", "Phyhipl", "Ppfibp1", "Fabp7", "Chn2", "Sepw1", "Sept8", "Arsb", "Omg", "Ptpre", "Timp4", "Apod", "Sulf2", "Ncald", "Zfp365", "Tnr", "Sh3bp4", "Gp1bb", "Sox8", "Mycl", "Trio", "S100a6", "Edil3", "Traf4", "Car2", "Cask", "Camsap2", "Scamp2", "Opcml", "Cdv3", "Lims2", "Serpine2")
dat7 <- c("9630013A20Rik", "Arpc1b", "Tmem2", "Sema6a", "Mag", "Mobp", "Mbp", "Rras2", "Prom1", "Eml1", "Rims2", "Idh1", "Ttyh1", "Tmem141", "Adamts4", "Tspan2", "Dhcr24", "Kndc1", "Plekha1", "Rap2a", "Sox2ot", "Plat", "Slc12a2", "Cldn11", "Rnf122", "Cnksr3", "Ugt8a", "Gjc2", "Fam107b", "Tmem163", "Myo6", "Mar8", "Plp1", "Elovl7", "Cyp51", "Lpar1", "Ddr1", "Ddc", "Gng12", "Acat2", "Tmem125", "Gpr37", "Capn5", "Sema4d", "Rell1", "Lap3", "Bace1", "Sema6d", "Idi1", "Ephb1")
dat8 <- c("Ctps", "Tmem141", "Opalin", "Ttyh1", "Plekha1", "Ddr1", "Cd9", "Sirt2", "Arpc1b", "Fam214a", "Mcam", "Nudt4", "Clic4", "Tspan2", "Serinc5", "Prom1", "Sh3gl3", "Plat", "Kndc1", "Olig1", "Lpar1", "Grb14", "Sema6a", "Rap2a", "Mag", "Slk", "Slc9a3r2", "9630013A20Rik", "Gpr17", "Tspan15", "Igsf8", "Tssc4", "Lgi3", "Arap2", "Klf13", "St8sia3", "Rims2", "Prickle1", "Ick", "Clmn", "Daam1", "Phyhipl", "Onecut2", "Bmp1", "Adamts4", "Ablim2", "Mobp", "Thbs3", "Zdhhc9", "Capn5")
dat9 <- c("Apod", "Sepp1", "S100b", "Trf", "Anln", "Hapln2", "Glul", "Ndrg1", "Sez6l2", "Klk6", "Gng11", "Sept8", "Sec11c", "Car2", "Spock3", "Serpinb1a", "Pmp22", "Zdhhc20", "Cpd", "Anxa5", "Slco3a1", "Enpp6", "Gatm", "Cd59a", "B3galt5", "Mgst3", "Tubb4a", "Etv1", "Ttyh2", "Psat1", "Elovl5", "Aspa", "Plin3", "Usp31", "Litaf", "Gjb1", "Cntn2", "Kctd3", "S100a1", "Prr5l", "Mid1ip1", "Pla2g16", "Dbndd2", "S100a16", "Dpy19l1", "S100a13", "S100a6", "Edil3", "Gstm7", "Kctd13")
dat10 <- c("Chn2", "Mpzl1", "Frmd4a", "Gjc3", "S100b", "Bcas1", "Pik3r3", "Cd9", "Epb4.1l2", "Enpp6", "Ptpre", "Sgk1", "Pdcd4", "Dct", "Tmem108", "Fam13c", "Fyn", "Ppfibp1", "Trio", "Grlf1", "Neu4", "Ust", "Tnr", "Vcan", "Phyhipl", "Sox4", "Sox2", "Lrrc42", "Zfp365", "Snx18", "Ptgds", "Rgs16", "Sez6l", "Ttyh2", "Slc1a1", "Panx1", "Rnd2", "Tcf7l2", "Tril", "Ctnnal1", "Pnpla3", "Rasgef1b", "Opcml", "Lcorl", "Ece1", "Mycl", "Brca1", "Dock4", "Tmem132b", "Adam11")
dat11 <- c("Mobp", "Ddr1", "Tspan2", "Cyp51", "Plat", "Tmem141", "Tmem2", "Mag", "Gnb4", "Rims2", "Cmtm5", "Trf", "Sema6a", "Rras2", "Enpp2", "Slc12a2", "Dhcr24", "Marcksl1", "Fam107b", "Ctps", "Gjc2", "Eml1", "Cers2", "Grb14", "Elovl6", "Gpr37", "Gamt", "Mmd2", "Plp1", "Rap2a", "Gsn", "Ugt8a", "Clic4", "Tmem125", "Idi1", "Acat2", "Fdps", "Sox2ot", "Kndc1", "Ttyh1", "Qdpr", "9630013A20Rik", "Paqr8", "Lap3", "Lss", "Arap2", "Rap1a", "Tmeff2", "Il23a", "Elovl1")
dat12 <- c("9630013A20Rik", "Tmem141", "Ctps", "Sema6a", "Mag", "Prom1", "Sirt2", "Gng12", "Plat", "Clic4", "Ddr1", "Nudt4", "Mcam", "Rims2", "Arap2", "Birc2", "Marcksl1", "Serinc5", "Tspan15", "Mmd2", "Mbp", "Olig1", "Ugt8a", "Slc12a2", "Mobp", "Elovl7", "Gnb4", "Tssc4", "Eml1", "Rap2a", "Kndc1", "Rap1a", "Gamt", "Onecut2", "Sh3gl3", "Tppp3", "Luzp2", "Adamts4", "Fam214a", "Klf13", "Dhcr24", "Fscn1", "Cldn11", "Gjc2", "Mgll", "Sox2ot", "Rab33a", "Tmem125", "Col9a3", "Golga7")
dat13 <- c("Mal", "Ptgds", "Evi2a-evi2b", "Wfdc18", "Efhd1", "Gprasp2", "Car2", "Resp18", "Pla2g16", "Zcchc12", "Nmral1", "Apod", "Scg2", "Ndrg1", "Plxnb1", "Tro", "Aspa", "Sst", "6330403K07Rik", "Pnck", "Tmem130", "Prkar1b", "Spock3", "Gap43", "Unc80", "Rgs17", "Ncald", "Syt4", "Rab3c", "Syn2", "Pgm2l1", "Hap1", "Gng2", "Gabra1", "Gda", "Rit2", "Tmem59l", "Cd200", "Pcp4", "Nsg2", "Vsnl1", "Myt1l", "Pcsk2", "Erc2", "A030009H04Rik", "Nsg1", "Ndrg2", "Atp1a3", "Phgdh", "A730017C20Rik")
dat14 <- c("S100b", "Sepp1", "Klk6", "Apod", "Ugt8a", "Cyp51", "Anxa5", "Enpp6", "Gsn", "Anln", "Cd59a", "Mgst3", "Trf", "Dhcr24", "Tppp3", "Idi1", "Elovl7", "Mobp", "Pmp22", "Sccpdh", "Spock3", "Hapln2", "Cpd", "Plp1", "Acat2", "Gltp", "Dbndd2", "Sqle", "Elovl5", "Mmd2", "Elovl1", "Mob3b", "Chst2", "Fez1", "Mbp", "Mog", "Scrg1", "Plat", "Reep3", "Rap1a", "Fdps", "Sema5a", "Sept8", "Gng11", "Mgll", "Plin3", "Sc5d", "Efr3b", "Far1", "Idh1")
dat15 <- c("Ptgds", "Opalin", "Qdpr", "Il33", "Apoe", "Cd63", "Jph4", "C030029H02Rik", "Clmn", "Ccp110", "Tmod1", "Bin1", "Car2", "Csrp1", "Cd9", "Tubb4a", "Neat1", "Mcam", "Grm3", "Erbb2ip", "Spock1", "Faim2", "Kcna1", "Dusp26", "Pls1", "Prickle1", "Grb14", "Fam13c", "Glul", "Ndrg2", "Nfia", "Kcna6", "Lgi3", "Pex5l", "Wfdc18", "Sema6d", "Cep97", "Gpt", "Olfml1", "Daam1", "Gpr37", "Abca8a", "Aebp1", "Synpr", "Npsr1", "Cerk", "Fbxo7", "Fbxl5", "Plekha1", "Sbf1")
dat16 <- c("Opalin", "Ptgds", "Fosb", "Dusp1", "Qdpr", "Car2", "Btg2", "Dnajb1", "Ccp110", "Grb14", "Klf4", "Tspan2", "Cd9", "Zfp36", "Tmod1", "Tnfaip6", "Hspa1a", "Cd63", "Kcna1", "Arrdc3", "Grm3", "Csrp1", "Mcam", "Sept4", "Serpinb1a", "Tob1", "Arc", "Jund", "Phldb1", "Clmn", "Junb", "Ddit4", "Cntn2", "Hspa1b", "Trim59", "Erbb2ip", "Erf", "Tubb4a", "Plekha1", "Ddit3", "Bin1", "Desi1", "Ier2", "Egr2", "Ttyh1", "Dusp26", "Daam1", "Lgi3", "Pim3", "Ninj2", "Ā")
dat17 <- c("S100b", "Klk6", "Anxa5", "Mgst3", "Sepp1", "Pmp22", "Cd59a", "Chst2", "Scrg1", "Apod", "Mmd2", "Mob3b", "Dhcr24", "Gsn", "Cpd", "Enpp6", "Tppp3", "Acat2", "Cyp51", "Ugt8a", "Sccpdh", "Spock3", "Sema5a", "Elovl7", "Trf", "Elovl5", "Idi1", "Cdkn1c", "Hapln2", "Rab37", "Aacs", "Plat", "Eml1", "Sqle", "Anln", "Gpd1", "Hn1l", "Dbndd2", "Lss", "Me1", "Tmem254b_loc3", "Far1", "Snca", "Tmem254b_loc1", "B3galt5", "Apbb2", "Tmem254b_loc2", "Hopx", "Gjb1", "Kcnj10")
dat18 <- c("Cyp51", "Dhcr24", "Pdlim2", "Elovl1", "Tmem125", "Ugt8a", "Gamt", "Elovl7", "Srd5a1", "Sqle", "Mgll", "Hcn2", "Mobp", "Plekha1", "Pcyt2", "Idi1", "Tsix", "Car14", "Gjb1", "Mmd2", "Elovl6", "Gsn", "Arpc1b", "0610007P14Rik", "Scg2", "Rhoc", "Acat2", "Hsd17b7", "S100a6", "Fam13c", "Slco3a1", "Atp1a2", "Pdxdc1", "Nceh1", "Marcksl1", "Glul", "Pop4", "Commd6", "Slc29a3", "Gltp", "Eno2", "Eml1", "R3hdm1", "Aacs", "Fdps", "Atp1a1", "Plp1", "Sc5d", "Ndufaf2", "Mrpl43")
dat19 <- c("Il33", "Apoe", "Ptgds", "Car2", "Edil3", "Gpr37", "Grm3", "Sepw1", "Fgfr2", "Qdpr", "Neat1", "Erbb2ip", "C030029H02Rik", "Ndrg1", "Plekhb1", "Tubb4a", "Cryab", "Sema6d", "Dpy19l1", "Dnajb2", "Sept4", "Npc1", "Scd1", "Ttyh2", "Cdk19", "Olfml1", "Spock1", "Gnai1", "Cyp2j12", "Serpinb1a", "Trim59", "Carhsp1", "Sez6l2", "Gatm", "Kcna6", "Clic4", "Ndrg2", "Pls1", "Tmcc3", "Dock9", "Abca8a", "Slain1", "Ppp1cc", "Tkt", "Cd82", "Galnt6", "Efnb3", "Ddx39b", "Npsr1", "Tmem229a")
dat20 <- c("Fosb", "Dusp1", "Dnajb1", "Btg2", "Hspa1a", "Hspa1b", "Tob1", "Mal", "Ddit4", "Klf4", "Jund", "Zfp36", "Arc", "Ypel2", "Nudt4", "Phldb1", "Tspan2", "Reep3", "Lpar1", "Egr2", "Ddit3", "Endod1", "Pim3", "Gsn", "Cd9", "Phyhipl", "Mbp", "Ier2", "Ppp1r16b", "Insig1", "Pcyt2", "Rap1a", "Srcin1", "Clic4", "Mcam", "Dhcr24", "Olig1", "Junb", "Marcksl1", "Cyp51", "Lgi3", "Ugt8a", "Erf", "Idi1", "Plekha1", "2810468N07Rik", "Sema6a", "Plekhb1", "Erbb2ip", "Abhd17b")
dat21 <- c("Serpinb1a", "Neat1", "Sepp1", "Pmp22", "Glul", "Eml1", "Pex5l", "Gng11", "Sec11c", "Ttr", "Hapln2", "Car2", "Slc1a2", "Klk6", "Pls1", "Hcn2", "S100a1", "Fam13c", "Spock3", "Grm3", "Aspa", "Atp1a2", "Apoe", "Arsg", "Kctd4", "Slco3a1", "Pcolce2", "E330020D12Rik", "Etv1", "B3galt5", "Cystm1", "Dock9", "Mt3", "C630043F03Rik", "S100b", "Serpind1", "Apln", "Skap2", "Aig1", "9330117O12Rik", "Tmem109", "Urod", "Schip1", "Fam213b", "A330049N07Rik", "Tmem144", "Atg4c", "Lhfpl2", "Polr3e", "Cdk5rap2")
dat22 <- c("Anxa5", "Klk6", "Mgst3", "S100b", "Tppp3", "Cyp51", "Sepp1", "Plekhb1", "Apod", "Gatm", "Ugt8a", "Gltp", "Mbp", "Hn1l", "Sccpdh", "Sirt2", "Gsn", "Tubb4a", "Pcyt2", "Tmem254b_loc2", "Rap1gds1", "Pdlim2", "Frmd8", "Reep3", "Tmem254b_loc1", "Sqle", "Tmem254b_loc3", "Ndrg1", "Cdkn1c", "Fa2h", "Cmtm5", "Golga7", "Acat2", "Sec11c", "S100a16", "Cpd", "Fdps", "Plp1", "Josd2", "Rab37", "Plat", "Phgdh", "Chst2", "Cnp", "Vamp3", "Elovl1", "Trf", "Pllp", "Elovl7", "Pmp22")
dat23 <- c("Car2", "Cntn2", "Gad2", "Dqx1", "Nrxn3", "Mef2c", "Arpp21", "Ryr2", "Gad1", "Nrgn", "Camk2a", "Rgs4", "Rgs7bp", "Rasgrp1", "C030029H02Rik", "Hspa12a", "Pbx1", "Gabrg2", "Gabra1", "Chst1", "Vcl", "Bsn", "Slc32a1", "Kalrn", "2900060B14Rik", "Gabrb2", "Npcd", "Prrt2", "Actl6b", "Gm410", "Brinp1", "Snord68", "5031426D15Rik", "Cwc22", "D430041D05Rik", "Zbtb16", "Slc8a2", "Unc13a", "Camkk2", "Sdr42e1", "Zbtb38", "Pcdh15", "Rnf150", "Atp2b4", "Cers4", "Csmd1", "Ttn", "Ppm1h", "Sv2b", "Zfp71-rs1")

# Retrieve all vectors into a list
all_vectors <- mget(paste0("dat", 1:23))

# Merge and remove duplicates using Reduce and union
merged_genes <- Reduce(union, all_vectors)
length(merged_genes) #532


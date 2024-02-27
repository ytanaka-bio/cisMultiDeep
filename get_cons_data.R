#module load r/4.0.5
#get conserved genes
mart <- read.table("mart_export.txt.gz",sep="\t",header=T)
mart <- mart[which(mart[,2] != "" & mart[,4] != "" & mart[,6] != "" & mart[,8] != ""),]

mart <- mart[which(duplicated(mart[,2])==F),]
mart <- mart[which(duplicated(mart[,4])==F),]
mart <- mart[which(duplicated(mart[,6])==F),]
mart <- mart[which(duplicated(mart[,8])==F),]

#max score
max <- 300

#RNA
cons_rna_list <- mart[,c(6,2,4,8)]
colnames(cons_rna_list) <- c("Mouse","Human","Macaque","Marmoset")

#RNA
Human_rna <- read.table("gene_expression/Human_rna_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Macaque_rna <- read.table("gene_expression/Macaque_rna_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Marmoset_rna <- read.table("gene_expression/Marmoset_rna_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Mouse_rna <- read.table("gene_expression/Mouse_rna_dif.csv",sep=',',header=T,row.names=1,check.names=F)

#common celltype across species
#celltype_rna <- intersect(colnames(Human_rna),colnames(Macaque_rna))
#celltype_rna <- intersect(celltype_rna,colnames(Marmoset_rna))
#celltype_rna <- intersect(celltype_rna,colnames(Mouse_rna))
#celltype_rna <- sort(celltype_rna)
celltype_rna <- colnames(Mouse_rna)

#get list of orthologus gene
rownames(cons_rna_list) <- cons_rna_list[,2]
cons_rna_list <- cons_rna_list[intersect(rownames(cons_rna_list),rownames(Human_rna)),]
rownames(cons_rna_list) <- cons_rna_list[,3]
cons_rna_list <- cons_rna_list[intersect(rownames(cons_rna_list),rownames(Macaque_rna)),]
rownames(cons_rna_list) <- cons_rna_list[,4]
cons_rna_list <- cons_rna_list[intersect(rownames(cons_rna_list),rownames(Marmoset_rna)),]
rownames(cons_rna_list) <- cons_rna_list[,1]
cons_rna_list <- cons_rna_list[intersect(rownames(cons_rna_list),rownames(Mouse_rna)),]

#Focus on conserved peaks
Human_rna_cons <- Human_rna[cons_rna_list[,2],celltype_rna]
Macaque_rna_cons <- Macaque_rna[cons_rna_list[,3],celltype_rna]
Marmoset_rna['Sst Chodl'] <- 0     #Add Sst Chodl columns, since this cell type is absent in Marmoset data
Marmoset_rna_cons <- Marmoset_rna[cons_rna_list[,4],celltype_rna]
Mouse_rna_cons <- Mouse_rna[cons_rna_list[,1],celltype_rna]
Human_rna_cons[Human_rna_cons > max] <- max
Human_rna_cons[Human_rna_cons < -1*max] <- -1*max
Macaque_rna_cons[Macaque_rna_cons > max] <- max
Macaque_rna_cons[Macaque_rna_cons < -1*max] <- -1*max
Marmoset_rna_cons[Marmoset_rna_cons > max] <- max
Marmoset_rna_cons[Marmoset_rna_cons < -1*max] <- -1*max
Mouse_rna_cons[Mouse_rna_cons > max] <- max
Mouse_rna_cons[Mouse_rna_cons < -1*max] <- -1*max


#mCG
cons_mCG_list <- mart[,c(6,2,4,8)]
colnames(cons_mCG_list) <- c("Mouse","Human","Macaque","Marmoset")

#mCG
Human_mCG <- read.table("methylation/Human_mCG_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Macaque_mCG <- read.table("methylation/Macaque_mCG_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Marmoset_mCG <- read.table("methylation/Marmoset_mCG_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Mouse_mCG <- read.table("methylation/Mouse_mCG_dif.csv",sep=',',header=T,row.names=1,check.names=F)

#common celltype across species
celltype_mCG <- intersect(colnames(Human_mCG),colnames(Macaque_mCG))
celltype_mCG <- intersect(celltype_mCG,colnames(Marmoset_mCG))
celltype_mCG <- intersect(celltype_mCG,colnames(Mouse_mCG))
celltype_mCG <- sort(celltype_mCG)

#get list of orthologus gene
rownames(cons_mCG_list) <- cons_mCG_list[,2]
cons_mCG_list <- cons_mCG_list[intersect(rownames(cons_mCG_list),rownames(Human_mCG)),]
rownames(cons_mCG_list) <- cons_mCG_list[,3]
cons_mCG_list <- cons_mCG_list[intersect(rownames(cons_mCG_list),rownames(Macaque_mCG)),]
rownames(cons_mCG_list) <- cons_mCG_list[,4]
cons_mCG_list <- cons_mCG_list[intersect(rownames(cons_mCG_list),rownames(Marmoset_mCG)),]
rownames(cons_mCG_list) <- cons_mCG_list[,1]
cons_mCG_list <- cons_mCG_list[intersect(rownames(cons_mCG_list),rownames(Mouse_mCG)),]

#Focus on conserved peaks
Human_mCG_cons <- Human_mCG[cons_mCG_list[,2],celltype_mCG]
Macaque_mCG_cons <- Macaque_mCG[cons_mCG_list[,3],celltype_mCG]
Marmoset_mCG_cons <- Marmoset_mCG[cons_mCG_list[,4],celltype_mCG]
Mouse_mCG_cons <- Mouse_mCG[cons_mCG_list[,1],celltype_mCG]
Human_mCG_cons[Human_mCG_cons > max] <- max
Human_mCG_cons[Human_mCG_cons < -1*max] <- -1*max
Macaque_mCG_cons[Macaque_mCG_cons > max] <- max
Macaque_mCG_cons[Macaque_mCG_cons < -1*max] <- -1*max
Marmoset_mCG_cons[Marmoset_mCG_cons > max] <- max
Marmoset_mCG_cons[Marmoset_mCG_cons < -1*max] <- -1*max
Mouse_mCG_cons[Mouse_mCG_cons > max] <- max
Mouse_mCG_cons[Mouse_mCG_cons < -1*max] <- -1*max


#mCH
cons_mCH_list <- mart[,c(6,2,4,8)]
colnames(cons_mCH_list) <- c("Mouse","Human","Macaque","Marmoset")

#mCH
Human_mCH <- read.table("methylation/Human_mCH_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Macaque_mCH <- read.table("methylation/Macaque_mCH_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Marmoset_mCH <- read.table("methylation/Marmoset_mCH_dif.csv",sep=',',header=T,row.names=1,check.names=F)
Mouse_mCH <- read.table("methylation/Mouse_mCH_dif.csv",sep=',',header=T,row.names=1,check.names=F)

#common celltype across species
celltype_mCH <- intersect(colnames(Human_mCH),colnames(Macaque_mCH))
celltype_mCH <- intersect(celltype_mCH,colnames(Marmoset_mCH))
celltype_mCH <- intersect(celltype_mCH,colnames(Mouse_mCH))
celltype_mCH <- sort(celltype_mCH)

#get list of orthologus gene
rownames(cons_mCH_list) <- cons_mCH_list[,2]
cons_mCH_list <- cons_mCH_list[intersect(rownames(cons_mCH_list),rownames(Human_mCH)),]
rownames(cons_mCH_list) <- cons_mCH_list[,3]
cons_mCH_list <- cons_mCH_list[intersect(rownames(cons_mCH_list),rownames(Macaque_mCH)),]
rownames(cons_mCH_list) <- cons_mCH_list[,4]
cons_mCH_list <- cons_mCH_list[intersect(rownames(cons_mCH_list),rownames(Marmoset_mCH)),]
rownames(cons_mCH_list) <- cons_mCH_list[,1]
cons_mCH_list <- cons_mCH_list[intersect(rownames(cons_mCH_list),rownames(Mouse_mCH)),]

#Focus on conserved peaks
Human_mCH_cons <- Human_mCH[cons_mCH_list[,2],celltype_mCH]
Macaque_mCH_cons <- Macaque_mCH[cons_mCH_list[,3],celltype_mCH]
Marmoset_mCH_cons <- Marmoset_mCH[cons_mCH_list[,4],celltype_mCH]
Mouse_mCH_cons <- Mouse_mCH[cons_mCH_list[,1],celltype_mCH]
Human_mCH_cons[Human_mCH_cons > max] <- max
Human_mCH_cons[Human_mCH_cons < -1*max] <- -1*max
Macaque_mCH_cons[Macaque_mCH_cons > max] <- max
Macaque_mCH_cons[Macaque_mCH_cons < -1*max] <- -1*max
Marmoset_mCH_cons[Marmoset_mCH_cons > max] <- max
Marmoset_mCH_cons[Marmoset_mCH_cons < -1*max] <- -1*max
Mouse_mCH_cons[Mouse_mCH_cons > max] <- max
Mouse_mCH_cons[Mouse_mCH_cons < -1*max] <- -1*max

#Get conserved ATAC peaks across species
Mouse_Human <- read.table("Mouse_atac_Human_overlap.bed",sep="\t",header=F)
Mouse_Macaque <- read.table("Mouse_atac_Macaque_overlap.bed",sep="\t",header=F)
Mouse_Marmoset <- read.table("Mouse_atac_Marmoset_overlap.bed",sep="\t",header=F)

Mouse_Human <- Mouse_Human[which(duplicated(Mouse_Human[,8]) == F),]
Mouse_Human <- Mouse_Human[which(duplicated(Mouse_Human[,4]) == F),]
Mouse_Macaque <- Mouse_Macaque[which(duplicated(Mouse_Macaque[,4]) == F),]
Mouse_Macaque <- Mouse_Macaque[which(duplicated(Mouse_Macaque[,8]) == F),]
Mouse_Marmoset <- Mouse_Marmoset[which(duplicated(Mouse_Marmoset[,8]) == F),]
Mouse_Marmoset <- Mouse_Marmoset[which(duplicated(Mouse_Marmoset[,4]) == F),]

rownames(Mouse_Human) <- Mouse_Human[,4]
rownames(Mouse_Macaque) <- Mouse_Macaque[,4]
rownames(Mouse_Marmoset) <- Mouse_Marmoset[,4]

cons_atac <- intersect(rownames(Mouse_Human),rownames(Mouse_Macaque))
cons_atac <- intersect(cons_atac,rownames(Mouse_Marmoset))
cons_atac_list <- data.frame(Mouse=cons_atac,Human=Mouse_Human[cons_atac,8],Macaque=Mouse_Macaque[cons_atac,8],Marmoset=Mouse_Marmoset[cons_atac,8])
write.table(cons_atac_list,'cons_atac_list.csv',sep=",",quote=F,row.names=F)

#Read cell type-specificity data of ATAC
Human_atac <- read.csv("chromatin_access/Human_atac_dif.csv",header=T,row.names=1,check.names=F)
Macaque_atac <- read.csv("chromatin_access/Macaque_atac_dif.csv",header=T,row.names=1,check.names=F)
Marmoset_atac <- read.csv("chromatin_access/Marmoset_atac_dif.csv",header=T,row.names=1,check.names=F)
Mouse_atac <- read.csv("chromatin_access/Mouse_atac_dif.csv",header=T,row.names=1,check.names=F)

#common celltypes across species
#celltype_atac <- intersect(colnames(Human_atac),colnames(Mouse_atac))
#celltype_atac <- intersect(colnames(Human_atac),colnames(Macaque_atac))
#celltype_atac <- intersect(celltype_atac,colnames(Marmoset_atac))
#celltype_atac <- intersect(celltype_atac,colnames(Mouse_atac))
#celltype_atac <- sort(celltype_atac)
celltype_atac <- colnames(Mouse_atac)

#get common celltypes and genes
celltype <- intersect(celltype_rna,celltype_mCG)
celltype <- intersect(celltype,celltype_mCH)
celltype <- intersect(celltype,celltype_atac)

common_gene <- intersect(rownames(cons_rna_list),rownames(cons_mCG_list))
common_gene <- intersect(common_gene,rownames(cons_mCH_list))

cons_rna_list <- cons_rna_list[common_gene,]
cons_mCG_list <- cons_mCG_list[common_gene,]
cons_mCH_list <- cons_mCH_list[common_gene,]

Human_rna_cons <- Human_rna_cons[cons_rna_list[,2],celltype]
Macaque_rna_cons <- Macaque_rna_cons[cons_rna_list[,3],celltype]
Marmoset_rna_cons <- Marmoset_rna_cons[cons_rna_list[,4],celltype]
Mouse_rna_cons <- Mouse_rna_cons[cons_rna_list[,1],celltype]

Human_mCG_cons <- Human_mCG_cons[cons_mCG_list[,2],celltype]
Macaque_mCG_cons <- Macaque_mCG_cons[cons_mCG_list[,3],celltype]
Marmoset_mCG_cons <- Marmoset_mCG_cons[cons_mCG_list[,4],celltype]
Mouse_mCG_cons <- Mouse_mCG_cons[cons_mCG_list[,1],celltype]

Human_mCH_cons <- Human_mCH_cons[cons_mCH_list[,2],celltype]
Macaque_mCH_cons <- Macaque_mCH_cons[cons_mCH_list[,3],celltype]
Marmoset_mCH_cons <- Marmoset_mCH_cons[cons_mCH_list[,4],celltype]
Mouse_mCH_cons <- Mouse_mCH_cons[cons_mCH_list[,1],celltype]

Human_atac <- Human_atac[,celltype]
Macaque_atac <- Macaque_atac[,celltype]
Marmoset_atac['Sst Chodl'] <- 0
Marmoset_atac <- Marmoset_atac[,celltype]
Mouse_atac <- Mouse_atac[,celltype]

#RNA output
all <- (Mouse_rna_cons + Human_rna_cons + Macaque_rna_cons + Marmoset_rna_cons)/4

#write CSV file
write.table(all,"all_rna_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_rna_list,"cons_rna_list.csv",sep=",",quote=F,row.names=F)

#mCG output
all <- (Mouse_mCG_cons + Human_mCG_cons + Macaque_mCG_cons + Marmoset_mCG_cons)/4

#write CSV file
write.table(all,"all_mCG_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_mCG_list,"cons_mCG_list.csv",sep=",",quote=F,row.names=F)

#mCH output
all <- (Mouse_mCH_cons + Human_mCH_cons + Macaque_mCH_cons + Marmoset_mCH_cons)/4

#write CSV file
write.table(all,"all_mCH_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_mCH_list,"cons_mCH_list.csv",sep=",",quote=F,row.names=F)

#Make slimed ATAC dif profile
num_peak <- 10000
selected_peak <- c()
for(i in 1:length(celltype)){
      dif <- Mouse_atac[,celltype[i]]
      names(dif) <- rownames(Mouse_atac)
      selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=T))[1:num_peak])
      #selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=F))[1:num_peak])
}
selected_peak <- unique(selected_peak)
write.table(Mouse_atac[selected_peak,],"Mouse_atac_dif_selected.csv",sep=",",quote=F,col.names=NA)

num_peak <- 10000
selected_peak <- c()
for(i in 1:length(celltype)){
      dif <- Human_atac[,celltype[i]]
      names(dif) <- rownames(Human_atac)
      selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=T))[1:num_peak])
      #selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=F))[1:num_peak])
}
selected_peak <- unique(selected_peak)
write.table(Human_atac[selected_peak,],"Human_atac_dif_selected.csv",sep=",",quote=F,col.names=NA)

num_peak <- 10000
selected_peak <- c()
for(i in 1:length(celltype)){
      dif <- Macaque_atac[,celltype[i]]
      names(dif) <- rownames(Macaque_atac)
      selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=T))[1:num_peak])
      #selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=F))[1:num_peak])
}
selected_peak <- unique(selected_peak)
write.table(Macaque_atac[selected_peak,],"Macaque_atac_dif_selected.csv",sep=",",quote=F,col.names=NA)

num_peak <- 10000
selected_peak <- c()
for(i in 1:length(celltype)){
      dif <- Marmoset_atac[,celltype[i]]
      names(dif) <- rownames(Marmoset_atac)
      selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=T))[1:num_peak])
      #selected_peak <- c(selected_peak,names(sort(rank(dif),decreasing=F))[1:num_peak])
}
selected_peak <- unique(selected_peak)
write.table(Marmoset_atac[selected_peak,],"Marmoset_atac_dif_selected.csv",sep=",",quote=F,col.names=NA)

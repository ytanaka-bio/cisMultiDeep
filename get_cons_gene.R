#module load r/4.0.5
#get conserved genes
mart <- read.table("mart_export.txt.gz",sep="\t",header=T)
mart <- mart[which(mart[,2] != "" & mart[,4] != "" & mart[,6] != "" & mart[,8] != ""),]

mart <- mart[which(duplicated(mart[,2])==F),]
mart <- mart[which(duplicated(mart[,4])==F),]
mart <- mart[which(duplicated(mart[,6])==F),]
mart <- mart[which(duplicated(mart[,8])==F),]

#RNA
cons_rna_list <- mart[,c(6,2,4,8)]
colnames(cons_rna_list) <- c("Mouse","Human","Macaque","Marmoset")

#max score
max <- 300

#RNA
Human_rna <- read.table("gene_expression/Human_rna_dif.csv",sep=',',header=T,row.names=1)
Macaque_rna <- read.table("gene_expression/Macaque_rna_dif.csv",sep=',',header=T,row.names=1)
Marmoset_rna <- read.table("gene_expression/Marmoset_rna_dif.csv",sep=',',header=T,row.names=1)
Mouse_rna <- read.table("gene_expression/Mouse_rna_dif.csv",sep=',',header=T,row.names=1)

#common celltypes across species
celltype <- intersect(colnames(Human_rna),colnames(Macaque_rna))
celltype <- intersect(celltype,colnames(Marmoset_rna))
celltype <- intersect(celltype,colnames(Mouse_rna))
celltype <- sort(celltype)

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
Human_rna_cons <- Human_rna[cons_rna_list[,2],celltype]
Macaque_rna_cons <- Macaque_rna[cons_rna_list[,3],celltype]
Marmoset_rna_cons <- Marmoset_rna[cons_rna_list[,4],celltype]
Mouse_rna_cons <- Mouse_rna[cons_rna_list[,1],celltype]
Human_rna_cons[Human_rna_cons > max] <- max
Human_rna_cons[Human_rna_cons < -1*max] <- -1*max
Macaque_rna_cons[Macaque_rna_cons > max] <- max
Macaque_rna_cons[Macaque_rna_cons < -1*max] <- -1*max
Marmoset_rna_cons[Marmoset_rna_cons > max] <- max
Marmoset_rna_cons[Marmoset_rna_cons < -1*max] <- -1*max
Mouse_rna_cons[Mouse_rna_cons > max] <- max
Mouse_rna_cons[Mouse_rna_cons < -1*max] <- -1*max

all <- (Mouse_rna_cons + Human_rna_cons + Macaque_rna_cons + Marmoset_rna_cons)/4

#write CSV file
write.table(all,"all_rna_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_rna_list,"cons_rna_list.csv",sep=",",quote=F,row.names=F)

#mCG
cons_mCG_list <- mart[,c(6,2,4,8)]
colnames(cons_mCG_list) <- c("Mouse","Human","Macaque","Marmoset")

Human_mCG <- read.table("methylation/Human_mCG_dif.csv",sep=',',header=T,row.names=1)
Macaque_mCG <- read.table("methylation/Macaque_mCG_dif.csv",sep=',',header=T,row.names=1)
Marmoset_mCG <- read.table("methylation/Marmoset_mCG_dif.csv",sep=',',header=T,row.names=1)
Mouse_mCG <- read.table("methylation/Mouse_mCG_dif.csv",sep=',',header=T,row.names=1)

#common celltypes across species
celltype <- intersect(colnames(Human_mCG),colnames(Macaque_mCG))
celltype <- intersect(celltype,colnames(Marmoset_mCG))
celltype <- intersect(celltype,colnames(Mouse_mCG))
celltype <- sort(celltype)

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
Human_mCG_cons <- Human_mCG[cons_mCG_list[,2],celltype]
Macaque_mCG_cons <- Macaque_mCG[cons_mCG_list[,3],celltype]
Marmoset_mCG_cons <- Marmoset_mCG[cons_mCG_list[,4],celltype]
Mouse_mCG_cons <- Mouse_mCG[cons_mCG_list[,1],celltype]
Human_mCG_cons[Human_mCG_cons > max] <- max
Human_mCG_cons[Human_mCG_cons < -1*max] <- -1*max
Macaque_mCG_cons[Macaque_mCG_cons > max] <- max
Macaque_mCG_cons[Macaque_mCG_cons < -1*max] <- -1*max
Marmoset_mCG_cons[Marmoset_mCG_cons > max] <- max
Marmoset_mCG_cons[Marmoset_mCG_cons < -1*max] <- -1*max
Mouse_mCG_cons[Mouse_mCG_cons > max] <- max
Mouse_mCG_cons[Mouse_mCG_cons < -1*max] <- -1*max

all <- (Mouse_mCG_cons + Human_mCG_cons + Macaque_mCG_cons + Marmoset_mCG_cons)/4

#write CSV file
write.table(all,"all_mCG_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_mCG_list,"cons_mCG_list.csv",sep=",",quote=F,row.names=F)

#mCH
cons_mCH_list <- mart[,c(6,2,4,8)]
colnames(cons_mCH_list) <- c("Mouse","Human","Macaque","Marmoset")

#Read cell type-specificity data of mCH
Human_mCH <- read.table("methylation/Human_mCH_dif.csv",sep=',',header=T,row.names=1)
Macaque_mCH <- read.table("methylation/Macaque_mCH_dif.csv",sep=',',header=T,row.names=1)
Marmoset_mCH <- read.table("methylation/Marmoset_mCH_dif.csv",sep=',',header=T,row.names=1)
Mouse_mCH <- read.table("methylation/Mouse_mCH_dif.csv",sep=',',header=T,row.names=1)

#common celltypes across species
celltype <- intersect(colnames(Human_mCH),colnames(Macaque_mCH))
celltype <- intersect(celltype,colnames(Marmoset_mCH))
celltype <- intersect(celltype,colnames(Mouse_mCH))
celltype <- sort(celltype)

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
Human_mCH_cons <- Human_mCH[cons_mCH_list[,2],celltype]
Macaque_mCH_cons <- Macaque_mCH[cons_mCH_list[,3],celltype]
Marmoset_mCH_cons <- Marmoset_mCH[cons_mCH_list[,4],celltype]
Mouse_mCH_cons <- Mouse_mCH[cons_mCH_list[,1],celltype]
Human_mCH_cons[Human_mCH_cons > max] <- max
Human_mCH_cons[Human_mCH_cons < -1*max] <- -1*max
Macaque_mCH_cons[Macaque_mCH_cons > max] <- max
Macaque_mCH_cons[Macaque_mCH_cons < -1*max] <- -1*max
Marmoset_mCH_cons[Marmoset_mCH_cons > max] <- max
Marmoset_mCH_cons[Marmoset_mCH_cons < -1*max] <- -1*max
Mouse_mCH_cons[Mouse_mCH_cons > max] <- max
Mouse_mCH_cons[Mouse_mCH_cons < -1*max] <- -1*max

all <- (Mouse_mCH_cons + Human_mCH_cons + Macaque_mCH_cons + Marmoset_mCH_cons)/4

#write CSV file
write.table(all,"all_mCH_dif_cons.csv",sep=",",quote=F,col.names=NA)
write.table(cons_mCH_list,"cons_mCH_list.csv",sep=",",quote=F,row.names=F)


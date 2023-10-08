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

#max score
max <- 300

#Read cell type-specificity data of ATAC
Human_atac <- read.csv("chromatin_access/Human_atac_dif.csv",header=T,row.names=1)
Macaque_atac <- read.csv("chromatin_access/Macaque_atac_dif.csv",header=T,row.names=1)
Marmoset_atac <- read.csv("chromatin_access/Marmoset_atac_dif.csv",header=T,row.names=1)
Mouse_atac <- read.csv("chromatin_access/Mouse_atac_dif.csv",header=T,row.names=1)

#common celltypes across species
celltype <- intersect(colnames(Human_atac),colnames(Mouse_atac))
celltype <- intersect(colnames(Human_atac),colnames(Macaque_atac))
celltype <- intersect(celltype,colnames(Marmoset_atac))
celltype <- intersect(celltype,colnames(Mouse_atac))
celltype <- sort(celltype)

#Focus on conserved peaks
Human_atac_cons <- Human_atac[cons_atac_list[,2],celltype]
Macaque_atac_cons <- Macaque_atac[cons_atac_list[,3],celltype]
Marmoset_atac_cons <- Marmoset_atac[cons_atac_list[,4],celltype]
Mouse_atac_cons <- Mouse_atac[cons_atac_list[,1],celltype]
Human_atac_cons[Human_atac_cons > max] <- max
Human_atac_cons[Human_atac_cons < -1*max] <- -1*max
Macaque_atac_cons[Macaque_atac_cons > max] <- max
Macaque_atac_cons[Macaque_atac_cons < -1*max] <- -1*max
Marmoset_atac_cons[Marmoset_atac_cons > max] <- max
Marmoset_atac_cons[Marmoset_atac_cons < -1*max] <- -1*max
Mouse_atac_cons[Mouse_atac_cons > max] <- max
Mouse_atac_cons[Mouse_atac_cons < -1*max] <- -1*max

all <- (Mouse_atac_cons + Human_atac_cons + Macaque_atac_cons + Marmoset_atac_cons)/4

#write CSV file
write.table(all,"all_atac_dif_cons.csv",sep=",",quote=F,col.names=NA)
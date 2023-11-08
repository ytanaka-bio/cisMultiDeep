#ATAC
Human_atac_bed <- read.table("Human_atac.bed",sep="\t",header=F)
Macaque_atac_bed <- read.table("Macaque_atac.bed",sep="\t",header=F)
Marmoset_atac_bed <- read.table("Marmoset_atac.bed",sep="\t",header=F)
Mouse_atac_bed <- read.table("Mouse_atac.bed",sep="\t",header=F)

Human_peak <- read.table("Human_atac_10000_deep_DR_meanshap.csv",sep=",",header=T,row.names=1,check.names=F)
Human_peak <- rownames(Human_peak)
Macaque_peak <- read.table("Macaque_atac_10000_deep_DR_meanshap.csv",sep=",",header=T,row.names=1,check.names=F)
Macaque_peak <- rownames(Macaque_peak)
Marmoset_peak <- read.table("Marmoset_atac_10000_deep_DR_meanshap.csv",sep=",",header=T,row.names=1,check.names=F)
Marmoset_peak <- rownames(Marmoset_peak)
Mouse_peak <- read.table("Mouse_atac_10000_deep_DR_meanshap.csv",sep=",",header=T,row.names=1,check.names=F)
Mouse_peak <- rownames(Mouse_peak)

rownames(Human_atac_bed) <- Human_atac_bed[,4]
rownames(Macaque_atac_bed) <- Macaque_atac_bed[,4]
rownames(Marmoset_atac_bed) <- Marmoset_atac_bed[,4]
rownames(Mouse_atac_bed) <- Mouse_atac_bed[,4]


write.table(Mouse_atac_bed[intersect(rownames(Mouse_atac_bed),Mouse_peak),],"Mouse_atac_10000.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Human_atac_bed[intersect(rownames(Human_atac_bed),Human_peak),],"Human_atac_10000.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Macaque_atac_bed[intersect(rownames(Macaque_atac_bed),Macaque_peak),],"Macaque_atac_10000.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Marmoset_atac_bed[intersect(rownames(Marmoset_atac_bed),Marmoset_peak),],"Marmoset_atac_10000.bed",sep="\t",quote=F,row.names=F,col.names=F)

#RNA
Human_gene_bed <- read.table("Human_gene.bed",sep="\t",header=F)
Macaque_gene_bed <- read.table("Macaque_gene.bed",sep="\t",header=F)
Marmoset_gene_bed <- read.table("Marmoset_gene.bed",sep="\t",header=F)
Mouse_gene_bed <- read.table("Mouse_gene.bed",sep="\t",header=F)

cons_rna_list <- read.csv("cons_rna_list.csv",sep=",",header=T)
input <- read.csv("rna_600_input.csv",sep=",",header=T,check.names=F,nrow=2,row.names=1)

rownames(Human_gene_bed) <- Human_gene_bed[,4]
rownames(Macaque_gene_bed) <- Macaque_gene_bed[,4]
rownames(Marmoset_gene_bed) <- Marmoset_gene_bed[,4]
rownames(Mouse_gene_bed) <- Mouse_gene_bed[,4]
rownames(cons_rna_list) <- cons_rna_list[,1]

write.table(Mouse_gene_bed[intersect(rownames(Mouse_gene_bed),cons_rna_list[colnames(input),1]),],"Mouse_rna_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Human_gene_bed[intersect(rownames(Human_gene_bed),cons_rna_list[colnames(input),2]),],"Human_rna_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Macaque_gene_bed[intersect(rownames(Macaque_gene_bed),cons_rna_list[colnames(input),3]),],"Macaque_rna_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Marmoset_gene_bed[intersect(rownames(Marmoset_gene_bed),cons_rna_list[colnames(input),4]),],"Marmoset_rna_600.bed",sep="\t",quote=F,row.names=F,col.names=F)

#mCG
cons_mCG_list <- read.csv("cons_mCG_list.csv",sep=",",header=T)
input <- read.csv("mCG_600_input.csv",sep=",",header=T,check.names=F,nrow=2,row.names=1)
rownames(cons_mCG_list) <- cons_mCG_list[,1]

write.table(Mouse_gene_bed[intersect(rownames(Mouse_gene_bed),cons_mCG_list[colnames(input),1]),],"Mouse_mCG_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Human_gene_bed[intersect(rownames(Human_gene_bed),cons_mCG_list[colnames(input),2]),],"Human_mCG_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Macaque_gene_bed[intersect(rownames(Macaque_gene_bed),cons_mCG_list[colnames(input),3]),],"Macaque_mCG_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Marmoset_gene_bed[intersect(rownames(Marmoset_gene_bed),cons_mCG_list[colnames(input),4]),],"Marmoset_mCG_600.bed",sep="\t",quote=F,row.names=F,col.names=F)

#mCH
cons_mCH_list <- read.csv("cons_mCH_list.csv",sep=",",header=T)
input <- read.csv("mCH_600_input.csv",sep=",",header=T,check.names=F,nrow=2,row.names=1)
rownames(cons_mCH_list) <- cons_mCH_list[,1]

write.table(Mouse_gene_bed[intersect(rownames(Mouse_gene_bed),cons_mCH_list[colnames(input),1]),],"Mouse_mCH_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Human_gene_bed[intersect(rownames(Human_gene_bed),cons_mCH_list[colnames(input),2]),],"Human_mCH_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Macaque_gene_bed[intersect(rownames(Macaque_gene_bed),cons_mCH_list[colnames(input),3]),],"Macaque_mCH_600.bed",sep="\t",quote=F,row.names=F,col.names=F)
write.table(Marmoset_gene_bed[intersect(rownames(Marmoset_gene_bed),cons_mCH_list[colnames(input),4]),],"Marmoset_mCH_600.bed",sep="\t",quote=F,row.names=F,col.names=F)

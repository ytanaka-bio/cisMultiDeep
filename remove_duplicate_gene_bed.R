bed <- read.table("Human_gene.bed",sep="\t",header=F)
bed <- bed[which(duplicated(bed[,4])==F),]
write.table(bed,"Human_gene.bed",sep="\t",quote=F,row.names=F,col.names=F)

bed <- read.table("Macaque_gene.bed",sep="\t",header=F)
bed <- bed[which(duplicated(bed[,4])==F),]
write.table(bed,"Macaque_gene.bed",sep="\t",quote=F,row.names=F,col.names=F)

bed <- read.table("Marmoset_gene.bed",sep="\t",header=F)
bed <- bed[which(duplicated(bed[,4])==F),]
write.table(bed,"Marmoset_gene.bed",sep="\t",quote=F,row.names=F,col.names=F)

bed <- read.table("Mouse_gene.bed",sep="\t",header=F)
bed <- bed[which(duplicated(bed[,4])==F),]
write.table(bed,"Mouse_gene.bed",sep="\t",quote=F,row.names=F,col.names=F)
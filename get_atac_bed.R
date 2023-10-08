data <- read.csv("chromatin_access/Human_atac_dif.csv",header=T,row.names=1)
bed <- matrix(unlist(strsplit(rownames(data),"[:-]")),nrow(data),3,byrow=T)
bed <- cbind(bed,rownames(data))
write.table(bed,"Human_atac.bed",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.csv("chromatin_access/Macaque_atac_dif.csv",header=T,row.names=1)
bed <- matrix(unlist(strsplit(rownames(data),"[:-]")),nrow(data),3,byrow=T)
bed <- cbind(bed,rownames(data))
write.table(bed,"Macaque_atac.bed",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.csv("chromatin_access/Marmoset_atac_dif.csv",header=T,row.names=1)
bed <- matrix(unlist(strsplit(rownames(data),"[:-]")),nrow(data),3,byrow=T)
bed <- cbind(bed,rownames(data))
write.table(bed,"Marmoset_atac.bed",sep="\t",quote=F,row.names=F,col.names=F)

data <- read.csv("chromatin_access/Mouse_atac_dif.csv",header=T,row.names=1)
bed <- matrix(unlist(strsplit(rownames(data),"[:-]")),nrow(data),3,byrow=T)
bed <- cbind(bed,rownames(data))
write.table(bed,"Mouse_atac.bed",sep="\t",quote=F,row.names=F,col.names=F)

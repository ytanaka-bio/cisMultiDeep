sum_SHAP_by_HiC_v3 <- function(dif_file,atac_file,atac_cons_file,gene_files,inLoop_dir,outprefix,celltype,num_gene=600,sufix=c("rna","mCG","mCH"),orgs=c("Mouse","Human","Macaque","Marmoset"),...){
		loop_name <- function(x){
			  return(paste(x,collapse=" "))
		}

		#HiC prefix list
		list <- read.table("HiC_prefix.txt",sep="\t",header=T,row.names=1)

		dif <- read.table(dif_file,sep=",",header=T,row.names=1,check.names=F)
		#celltype <- colnames(dif)

		atac_cons <- read.table(atac_cons_file,sep=",",header=T)
		rownames(atac_cons) <- atac_cons[,1]

		atac_shap <- read.table(atac_file,sep=",",header=T,row.names=1,check.names=F)
		other_shap <- vector("list",length(gene_files))
		names(other_shap) <- gene_files
		for(i in 1:length(other_shap)){
		      other_shap[[i]] <- read.table(gene_files[i],sep=",",header=T,row.names=1,check.names=F)
		}

		#analyze in each cell type

		print(celltype)
		rank <- dif[,celltype]
		names(rank) <- rownames(dif)
		rank <- sort(rank(rank),decreasing=T)
		peak <- names(rank)[1:num_gene]

		#SHAP in each celltype
		score <- atac_shap[peak,celltype]
		names(score) <- peak
		other_shap_type <- other_shap
		for(j in 1:length(other_shap)){
		      other_shap_type[[j]] <- other_shap[[j]][,celltype]
		      names(other_shap_type[[j]]) <- rownames(other_shap[[j]])
		}

		cons_peak <- intersect(names(score),rownames(atac_cons))

		#Read HiC
		loop_atac <- vector("list",length(orgs))
		names(loop_atac) <- orgs
		HiC_atac <- loop_atac
		for(l in 1:length(orgs)){
		      HiC_file <- paste(inLoop_dir,"/",orgs[l],"_",list[celltype,1],".loop_atac.txt",sep="")
		      HiC_data <- read.table(HiC_file,sep="\t",header=F)
		      loop_atac[[l]] <- apply(HiC_data[,1:6],1,loop_name)
		      HiC_atac[[l]] <- data.frame(LOOP=loop_atac[[l]],ATAC=HiC_data[,11])
		      loop_atac[[l]] <- unique(loop_atac[[l]])
	        }

		for(j in 1:length(other_shap)){
		      print(sufix[j])
		      print(org[1])
		      HiC_file <- paste(inLoop_dir,"/",orgs[1],"_",list[celltype,1],".loop_",sufix[j],".txt",sep="")
		      HiC_data <- read.table(HiC_file,sep="\t",header=F)
		      loop <- apply(HiC_data[,1:6],1,loop_name)
		      for(k in 1:length(loop_atac[[1]])){
			     select <- which(loop == loop_atac[[1]][k])
			     select2 <- HiC_atac[[1]][which(HiC_atac[[1]][,1] == loop_atac[[1]][k]),2]
			     select2 <- intersect(select2,peak)
			     if(length(select) != 0 & length(select2) != 0){
				   #print(score[select2])
				   score[select2] = score[select2] + sum(HiC_data[select,7] * other_shap_type[[j]][HiC_data[select,11]])
			     }
		      }
		      for(l in 2:length(orgs)){
		      	    print(org[2])
		      	    names(cons_peak) <- atac_cons[cons_peak,l]
		      	   HiC_file <- paste(inLoop_dir,"/",orgs[l],"_",list[celltype,1],".loop_",sufix[j],".txt",sep="")
			   HiC_data <- read.table(HiC_file,sep="\t",header=F)
			   loop <- apply(HiC_data[,1:6],1,loop_name)
			   for(k in 1:length(loop_atac)){
			     	 select <- which(loop == loop_atac[[l]][k])
			     	 select2 <- HiC_atac[[l]][which(HiC_atac[[l]][,1] == loop_atac[[l]][k]),2]
			     	 select2 <- intersect(select2,names(cons_peak))
			     	 if(length(select) != 0 & length(select2) != 0){
				   		   #print(score[select2])
				   		   score[cons_peak[select2]] = score[cons_peak[select2]] + sum(HiC_data[select,7] * other_shap_type[[j]][HiC_data[select,11]])
			     	 }
		      	  }
		      }
		}
		#print(score)
		output <- paste(outprefix,"_",list[celltype,1],".txt",sep="")
		result <- data.frame(sort(score,decreasing=T))
		result <- cbind(result,rep(celltype,nrow(result)))
		write.table(result, output,sep=",",quote=F,col.names=F)
}
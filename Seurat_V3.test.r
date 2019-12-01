#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
file=args[1]
name=args[2]
#library(Seurat,lib.loc="/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/Clusteringtest/Seurat3")
library(EZsinglecell)
require(RColorBrewer)
require(dplyr)
library(raster)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#### DGE file preprocess
CheckContamination_V3<-function(object,genes=c("INS","GCG","PPY","SST"),tocheck="INS",yrange=c(0,0.1),name="contamination_violinplot"){
	as.data.frame(t(as.matrix(GetAssayData(object = object, slot = "counts")[genes,]))) %>% Tomerge_v2(object@meta.data,.) -> datatoplot
	datatoplot$RNA_snn_res.0.8<-factor(datatoplot$RNA_snn_res.0.8,levels=0:12)
	all.list<-list()
	count=0
	for(tocheck in genes){
		count=count+1
		p<-cbind(datatoplot,ratio=datatoplot[,tocheck]/as.numeric(as.character(datatoplot$nCount_RNA))) %>% ggplot(.)+aes(RNA_snn_res.0.8,log10(ratio),fill=RNA_snn_res.0.8)+geom_violin()+theme(legend.position="null")+labs(x="Clusters",y=paste("% of ",tocheck,"transcripts\n Log scale"),title=name)
		all.list[[count]]<-p
	}
	p<-as.data.frame(datatoplot[,c("INS","GCG")]) %>% ggplot(.)+aes(INS,GCG)+geom_point(size=0.5)
	grid.arrange(all.list[[1]],all.list[[2]],all.list[[3]],all.list[[4]],nrow=2)
	print(p)
}
dgepreprocess<-function(dge_data,Txcutoff=1000,norowname=T){
	if(norowname){
		rownames(dge_data)<-dge_data[,1]
		dge_data<-dge_data[,-1]
	}else{
		dge_data<-as.data.frame(dge_data)
	}
	dge_data<-rbind(dge_data,Sum=apply(dge_data,2,sum))
	new_dge<-dge_data[,which(dge_data["Sum",]>Txcutoff)]
	new_dge<-new_dge[1:(nrow(new_dge)-2),]
	new_dge<-cbind(new_dge,genesum=apply(new_dge,1,sum))
	new_dge<-subset(new_dge,genesum!=0)
	new_dge<-new_dge[,1:ncol(new_dge)-1]
	return(new_dge)
}
Tomerge<-function(A,B){
	mergeAB<-merge(as.matrix(A),as.matrix(B),by="row.names",all=TRUE)
	row.names(mergeAB)<-mergeAB[,1]
	mergeAB<-mergeAB[,-1]
	return(mergeAB)
}

Tomerge_v2<-function(A,B,leavex=T,leavey=F,withnumber=F){
        mergeAB<-merge(as.matrix(A),as.matrix(B),by="row.names",all.x=leavex,all.y=leavey)
        row.names(mergeAB)<-mergeAB[,1]
        mergeAB<-mergeAB[,-1]
        if(withnumber){
                for (i in 1: ncol(mergeAB)){
                        tmp<-as.numeric(as.character(mergeAB[,i]))
                        if(all(is.na(tmp))){
                                mergeAB[,i]<-mergeAB[,i]
                        }else{
                                mergeAB[,i]<-tmp
                        }
                }
        }
        return(mergeAB)
}
format2<-function(p,data,x,y,nolegend=T){
	data<-data[complete.cases(data),]
	p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_fixed( (max(data[,x])-min(data[,x]))/(max(data[,y])-min(data[,y])))
	if(nolegend){
		p<-p+theme(legend.position="none")
	}
	return(p)
}

#### Get top genes
Gettopgenes<-function(object,number){
	pbmc.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
	pbmc.markers %>% group_by(cluster) %>% top_n(number, avg_logFC) -> topgene
	return(topgene)
}
#### Get gene signature plots
Gettsnesignature<-function(object,signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"),toreturn=T,doUMAP=F,doPCA=F,dotSNE=F,dotsize=0.3,buttomgrey=T,islegend=T,highcolor="red",dotsize2=1.5,extratitle=""){
	geneNA<-signiture[which(!signiture %in% row.names(GetAssayData(object = object, slot = "scale.data")))]
        signiture.expre<-setdiff(signiture,geneNA)
	tSNE_signature<-as.data.frame(Tomerge_v2(object@reductions$tsne@cell.embeddings,t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])))
	UMAP_signature<-as.data.frame(Tomerge_v2(object@reductions$umap@cell.embeddings,t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])))
	for (gene in geneNA){
                tSNE_signature<-cbind(gene=0,tSNE_signature)
                names(tSNE_signature)[1]<-gene
        }
	toPlottSNEandPCA<-Tomerge(tSNE_signature,as.data.frame(object@reductions$pca@cell.embeddings)[,c(1,2,3,4)])
	toPlotUMAPandPCA<-Tomerge(UMAP_signature,as.data.frame(object@reductions$pca@cell.embeddings)[,c(1,2,3,4)])
	all.list<-list()
	if(dotSNE){
                for (i in 1:length(signiture)){
                        genename<-signiture[i]
                        if(genename %in% geneNA){
                                lowlim<-0
                                highlim<-0
                        }else{
                                lowlim<-min(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                                highlim<-max(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                        }
                        if(buttomgrey){
                                if(genename %in% geneNA){
                                	p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high="grey",limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)
                                }else{
					p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize2)
				}
                        }else{
                                p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string("tSNE_1","tSNE_2",color=genename),size=dotsize)
                        }
                        p<-format2(p,toPlottSNEandPCA,"tSNE_1","tSNE_2",nolegend=islegend)
                        all.list[[i]]<-p
                }
        }
	if(doUMAP){
                for (i in 1:length(signiture)){
                        genename<-signiture[i]
                        if(genename %in% geneNA){
                                lowlim<-0
                                highlim<-0
                        }else{
                                lowlim<-min(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                                highlim<-max(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                        }
                        if(buttomgrey){
                                if(genename %in% geneNA){
                                        p<-ggplot(subset(toPlotUMAPandPCA,get(genename)<=0))+geom_point(aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high="grey",limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"UMAP_space",extratitle))+geom_point(data=subset(toPlotUMAPandPCA,get(genename)>0),aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize)
                                }else{
                                        p<-ggplot(subset(toPlotUMAPandPCA,get(genename)<=0))+geom_point(aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"UMAP_space",extratitle))+geom_point(data=subset(toPlotUMAPandPCA,get(genename)>0),aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize2)
                                }
                        }else{
                                p<-ggplot(subset(toPlotUMAPandPCA,get(genename)<=0))+geom_point(aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"UMAP_space",extratitle))+geom_point(data=subset(toPlotUMAPandPCA,get(genename)>0),aes_string("UMAP_1","UMAP_2",color=genename),size=dotsize)
                        }
                        p<-format2(p,toPlotUMAPandPCA,"UMAP_1","UMAP_2",nolegend=islegend)
                        all.list[[i]]<-p
                }
        }

	if(doPCA){
		for (i in 1:length(signiture)){
			genename<-signiture[i]
                        if(genename %in% geneNA){
                                lowlim<-0
                                highlim<-0
                        }else{
                                lowlim<-min(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                                highlim<-max(t(GetAssayData(object = object, slot = "scale.data")[signiture.expre,])[,genename])
                        }
			if(buttomgrey){
                                if(genename %in% geneNA){
                                	p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high="grey",limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"PC_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
                                }else{
                                	p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
				}
                        }else{
                        	p<-ggplot(subset(toPlottSNEandPCA,get(genename)<=0))+geom_point(aes_string(PC1used,PC2used,color=genename),size=dotsize)+scale_color_gradient(low="grey",high=highcolor,limit=c(lowlim,highlim))+ggtitle(paste(signiture[i],"tSNE_space",extratitle))+geom_point(data=subset(toPlottSNEandPCA,get(genename)>0),aes_string(PC1used,PC2used,color=genename),size=dotsize)
                        }
                        p<-format2(p,toPlottSNEandPCA,PC1used,PC2used)
        #                print(p)
                        all.list<-c(all.list,list(p))
                }
        }
        return(all.list)
}

#### Plot ALL 
Plot_cluster<-function(object,name,dotsize=1,resolusion="RNA_snn_res.0.8",color="Paired"){
	titlename=name
	pdf(paste(name,".Seurat3.pdf",sep=''))
	PC_df<-merge(as.data.frame(object@reductions$pca@cell.embeddings),object@meta.data,by="row.names")
	getPalette = colorRampPalette(brewer.pal(9, "RdYlBu"))
	tSNE_df<-merge(object@reductions$tsne@cell.embeddings,object@meta.data,by="row.names")
	topgene<-Gettopgenes(object,10)
	UMAP_df<-merge(object@reductions$umap@cell.embeddings,object@meta.data,by="row.names")

	p1<-ggplot(tSNE_df)+aes(tSNE_1,tSNE_2,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename) # tSNE plot
	#DoHeatmap(object, genes.use = topgene$gene, order.by.ident = TRUE, col.use = rev(getPalette(10)),slim.col.label = TRUE, remove.key = TRUE, cexRow=0.5,cex.col=0.5)
	p2<-DoHeatmap(object, features = topgene$gene, group.by = "ident",group.bar = TRUE, disp.min = -2.5, disp.max = NULL,slot = "scale.data", assay = NULL, label = TRUE, size = 5.5,hjust = 0, angle = 45, raster = TRUE, draw.lines = TRUE,lines.width = NULL, group.bar.height = 0.02, combine = TRUE) # heatmap
        p3<-ggplot(PC_df)+aes(PC_1,PC_2,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename) # PCA  PC1 PC2 plot
        p4<-ggplot(PC_df)+aes(PC_3,PC_4,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename) # PCA  PC3 PC4 plot
	p5<-ggplot(UMAP_df)+aes(UMAP_1,UMAP_2,color=eval(parse(text=resolusion)))+geom_point(size=dotsize)+theme_bw()+coord_fixed()+scale_color_brewer(palette=color,guide_legend(title="Cell clustering"))+guides(color=guide_legend(override.aes=list(size=5)))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(titlename) # UMAP plot

	print(p1)
	print(p5)
	print(p2)
	print(p3)
	print(p4)
	sig.plot<-Gettsnesignature(object,dotSNE=T,doUMAP=F)
	cycles<-as.integer(length(sig.plot)/6)
	for (i in 0:(cycles-1)){
		grid.arrange(grobs=sig.plot[(6*i+1):(6*i+6)])
	}
	grid.arrange(grobs=sig.plot[(6*cycles+1):length(sig.plot)])
#	sig.plot<-Gettsnesignature(object,dotSNE=F,doUMAP=T)
 #       cycles<-as.integer(length(sig.plot)/6)
  #      for (i in 0:(cycles-1)){
   #             grid.arrange(grobs=sig.plot[(6*i+1):(6*i+6)])
    #    }
     #   grid.arrange(grobs=sig.plot[(6*cycles+1):length(sig.plot)])
	CheckContamination_V3(object,genes=c("INS","GCG","PPY","SST"))
	
#	for(i in 1:4){
#	DimHeatmap(object, dims = i, nfeatures = 500, cells = NULL,reduction = "pca", disp.min = -2.5, disp.max = NULL,balanced = TRUE, projected = FALSE, ncol = NULL, combine = TRUE,fast = TRUE, raster = TRUE, slot = "scale.data", assays = NULL)
#	}
	dev.off()
}
detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}
cutoff=500
data=read.table(file,stringsAsFactors=FALSE,header=T)
CLUSTER=docluster(dgepreprocess(data,cutoff,norowname=T),GetinformativeGene(dgepreprocess(data,cutoff,norowname=T),cutoff),"dropseq",reso=0.6)
Fullplot_v2(CLUSTER,paste(name,"Seurat2.pdf",sep='.'),topgene=NULL,resolusion="res.0.6",signiture=c("INS","GCG","SST","PPY","KRT19","COL1A2","REG1A","DNAJB1","GHRL"))
#detach("package:Seurat",unload=T)
detach_package("Seurat", TRUE)
library(Seurat,lib.loc="/mnt/rstor/genetics/JinLab/cxw486/Chip-seq/ATAC/scATAC/Clusteringtest/Seurat3")
DGE=dgepreprocess(data,500,T)
pbmc <- CreateSeuratObject(counts = DGE)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
pbmc <- RunUMAP(object = pbmc,dims = 1:5)
Plot_cluster(pbmc,name)


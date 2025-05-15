##SIGLNCOS
SIGLNCOS <- function(file = NULL, Tissue = NULL, Cell_annotation = NULL, GENCODE = ftq, nFeature_RNA_min = 200, nFeature_RNA_max = 2500, percent.mt_max = 10, P_threshold = 0.05){
  # ---- Confirm data integrity ----
  if (length(file) != 1) {
    if (length(file) == 0) {
      cli::cli_abort("{.arg file} cannot be empty.")
    }
    len <- length(file)
    file <- file[1]
    cli::cli_warn(c("{.arg file} must have length 1, not length {len}.", 
                    `!` = "Only the first, {.file {filename}}, will be used."))
  }
  
  if (is.null(Tissue)){
    cli::cli_abort("{.arg Tissue} cannot be empty.")
  }
  
  # ---- library packages ----
  lapply(c("Seurat","dplyr","psych","clusterProfiler","SingleR","celldex"), 
         library, character.only = T)
  
  # ---- load scRNA-seq and reference data ----
  scdata<-Read10X(file) ##barcodes.tsv.gz,features.tsv.gz,matrix.mts.gz
  
  # ---- Stage1:Constructing lncRNA Visualization landscape ----
  
  SO<-CreateSeuratObject(scdata,min.cells = 3, min.features = 200)##Create the Seurat object
  
  SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")##QC
  SO <- subset(SO, subset  = nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max & percent.mt < percent.mt_max)
  
  
  ## data pre-processing
  
  SO <- NormalizeData(object = SO,normalization.method = "LogNormalize", scale.factor = 10000)
  SO <- FindVariableFeatures(object = SO,selection.method = "vst", nfeatures = 2000)
  SO <- ScaleData(object = SO)
  SO <- RunPCA(object = SO, features = VariableFeatures(object = SO))
  SO <- FindNeighbors(object = SO,dims = 1:15)
  SO <- FindClusters(object = SO,resolution = 0.5)
  SO <- RunUMAP(object = SO,dims = 1:15)
  
  ## cell annotation
  
  if(!is.null(Cell_annotation)){
    print("Running with cell types")
    rownames(Cell_annotation) <- Cell_annotation[,1]
    SO$celltype <- Cell_annotation[colnames(SO),2]
  }else{
    print("Running with seurat clusters")
    SO$celltype <- SO$seurat_clusters
  }
  
  saveRDS(SO,"./Seurat_Object.rds")
  
  ##gene annotation with GENCODE
  ftq<-readRDS("./contaxt/GENCODE.rds")
  lncRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "lncRNA")%>%arrange("gene_name")
  lncRNA<-intersect(lncRNA$gene_name,rownames(SO))##lncRNA
  mRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "protein_coding")%>%arrange("gene_name")
  mRNA<-intersect(mRNA$gene_name,rownames(SO))##mRNA
  
  

  # ---- Stage2:Identifying cellular-signature lncRNA ----
  
  Idents(SO)<-SO$celltype##Calculated according to cell type
  SO.markers<-FindAllMarkers(SO,only.pos = T,)##Only high expression data were selected(avglog2FC>0.25,min.pct = 0.1)
  SO.markers$tumor<-tumor
  SO.markers<-subset(SO.markers,p_val_adj < P_threshold)##Genes with significant corrected p-values were selected
  SO.markers$tumor<-tumor
  lnc.markers<-SO.markers[which(SO.markers$gene %in% lncRNA),]##cellular_signature_lncRNA
  write.table(lnc.markers,"cellular_signature_lncRNA.txt",sep = "\t",quote = F,col.names = T,row.names = F)
  mRNA.markers<-SO.markers[which(SO.markers$gene %in% mRNA),]##cellular_signature_mRNA
  
  # ---- Stage3:Functional enrichment of clncRNAs ----
  
  ##mRNA correlations were calculated based on lncrnas（spearman）
  clncRNA<-unique(lnc.markers$gene)
  
  package_type=substr(packageVersion("Seurat"), 1, 1)
  
  if(package_type==5){
    exp.data<-GetAssayData(object = SO, slot = "data")##Extracting expression profiles
  }
  else{
    exp.data<-GetAssayData(object = SO, layer  = "data")##Extracting expression profiles
  }

  clnc.exp<-exp.data[which(rownames(exp.data) %in% clncRNA),]##Extract the clncRNA expression matrix
  mRNA.exp<-exp.data[which(rownames(exp.data) %in% mRNA),]##Extract the mRNA expression matrix
  options(warn=0)
  cor_data<-lapply(clncRNA,function(x){
    cell<-as.vector(subset(lnc.markers,gene == x)$cluster)
    aa<-lapply(cell,function(z){
      cellid<-colnames(subset(SO,cell_type == z))
      clnc.exp2<-clnc.exp[,cellid]
      mRNA.exp2<-mRNA.exp[,cellid]
      cmRNA<-subset(mRNA.markers,cluster == z)$gene
      bb<-lapply(mRNA, function(y){##mRNA for caculating all mRNA、cmRNA for mRNA with DE
        cor<-corr.test(clnc.exp2[x,],mRNA.exp2[y,],method = "spearman")
        return(c(x,y,z,cor$r,cor$p.adj))
      })
      return(do.call(rbind,bb))
    })
    do.call(rbind,aa)
  })
  cor_data<-as.data.frame(do.call(rbind,cor_data))
  colnames(cor_data)<-c("lncRNA","mRNA","cell_type","cor","p.adj")
  cor_data<-na.omit(cor_data)
  cor_data$tumor<-tumor
  write.table(cor_data,"lncRNA_mRNA_cor.txt",sep = "\t",quote = F,col.names = T,row.names = F)
  
  ##GSEA enrichment analysis was performed based on the sorted mRNA
  ##Read the functions used for enrichment
  GO<-read.gmt("./contaxt/c5.go.v2022.1.Hs.symbols.gmt")
  imm<-read.gmt("./contaxt/c7.all.v2022.1.Hs.symbols.gmt")
  kegg<-read.gmt("./contaxt/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
  ##function annotation
  cor_data$cor<-as.numeric(cor_data$cor)
  cor_data$p.adj<-as.numeric(cor_data$p.adj)##Convert to a numeric type
  cell<-as.vector(unique(cor_data$cell_type))
  lnc_function<-c()
  for(i in cell){
    subcor<-subset(cor_data,cell_type == i)
    clncRNA2<-unique(subcor$lncRNA)
    celldata<-lapply(clncRNA2,function(x){
      clncRNA_cor<-subset(subcor,lncRNA == x)$cor
      clncRNA_p<-subset(cor_data,lncRNA == x)$p.adj
      ##Replace the 0 value in padj with the minimum value
      minp<-min(clncRNA_p[-which(clncRNA_p==0)])
      clncRNA_p<-ifelse(clncRNA_p == 0,minp,clncRNA_p)
      Rank<- -log(as.numeric(clncRNA_p))*sign(as.numeric(clncRNA_cor))
      names(Rank)<-subset(subcor,lncRNA == x)$mRNA
      Rank<-sort(Rank,decreasing = T)##sort
      ##GSEA
      gseaGO<-GSEA(Rank ,TERM2GENE =  GO)
      gseaimmune<-GSEA(Rank,TERM2GENE = imm)
      gseaKEGG<-GSEA(Rank,TERM2GENE = kegg)
      ##Select significant results
      GOresult<-gseaGO@result[which(gseaGO@result$p.adjust<0.05),]
      immgsearesult<-gseaimmune@result[which(gseaimmune@result$p.adjust<0.05),]
      keggsearesult<-gseaKEGG@result[which(gseaKEGG@result$p.adjust<0.05),]
      ##Add lncRNA information
      if(nrow(GOresult)>0){
        GOresult$lncRNA<-x
      }
      if(nrow(immgsearesult)>0){
        immgsearesult$lncRNA<-x
      }
      if(nrow(keggsearesult)>0){
        keggsearesult$lncRNA<-x
      }
      ##Add cell information
      if(nrow(GOresult)>0){
        GOresult$cell_type<-i
      }
      if(nrow(immgsearesult)>0){
        immgsearesult$cell_type<-i
      }
      if(nrow(keggsearesult)>0){
        keggsearesult$cell_type<-i
      }
      ##Add classification labels
      if(nrow(GOresult)>0){
        GOresult$func_type<-"GO"
      }
      if(nrow(immgsearesult)>0){
        immgsearesult$func_type<-"IMMUNE"
      }
      if(nrow(keggsearesult)>0){
        keggsearesult$func_type<-"KEGG"
      }
      return(rbind(GOresult,immgsearesult,keggsearesult))
    })
    celldata<-do.call(rbind,celldata)
    lnc_function<-rbind(lnc_function,celldata)
  }
  lnc_function$tumor<-tumor
  write.table(lnc_function,"lncRNA_function_regulation.txt",sep = "\t",quote = F,col.names = T,row.names = F)
  
  
  # ---- Stage4:Identitying co-lncRNA pairs based on functional co-regulated ----
  ##分别根据三种功能计算
  func<-list(
    "GO" = subset(lnc_function,func_type == "GO"),
    "IMMUNE" = subset(lnc_function,func_type == "IMMUNE"),
    "KEGG" = subset(lnc_function,func_type == "KEGG")
  )
  ##计算每一对lncRNA的杰卡尔系数
  jacmatrix<-lapply(func,function(i){
    cell<-unique(i$cell_type)
    celljac<-c()
    for(j in cell){
      subdata<-subset(i,cell_type == j)
      lncRNA<-unique(subdata$lncRNA)
      if(length(lncRNA)<2){next}
      aa<-data.frame(t(apply(combn(lncRNA,2),2,function(x){
        func1<-subset(subdata,lncRNA == x[1])$ID
        func2<-subset(subdata,lncRNA == x[2])$ID
        jac<-length(intersect(func1,func2))/length(union(func1,func2))
        return(c(x,jac,j))
      })))
      celljac<-rbind(celljac,aa)
    }
    colnames(celljac)<-c("lncRNA1","lncRNA2","jaccord","cell_type")
    return(celljac)
  })
  
  
  ##对于每一种功能，根据lncRNA调控情况计算其权重
  W<-list()
  for(f in 1:3){
    fun<-func[[f]]##依次对三种功能计算权重
    lncnum<-length(unique(fun$lncRNA))##总的lncRNA数量
    fun2<-split(fun$lncRNA,fun$ID)##根据功能id分离lncRNA
    fun_lncnum<-unlist(lapply(fun2,function(x){length(unique(x))}))##调控每个功能的lncRNA数量
    w<- -log10(fun_lncnum/lncnum)##每个功能的权重
    W[[f]]<-w
  }
  ##Calculate the co-regulatory coefficients among lncrnas
  jacmatrix2<-list()
  for(f in 1:3){##Calculate the co-regulation coefficient based on the three functions
    fun<-func[[f]]##The regulation of this type of function by lncRNA
    w<-W[[f]]##The weight of this type of function
    jac<-jacmatrix[[f]]##The jaccard coefficients between lncrnas calculated based on this type of function
    jac2<-jac[which(jac$jaccord>=0.5),]##The lncRNA filtering relationship is correct according to jaccard
    if(nrow(jac2)<1){next}
    coindex<-c()##The co-regulatory coefficient of each pair of lncrnas
    for(i in 1:nrow(jac2)){##For each lncRNA relationship pair
      lncRNA1<-jac2$lncRNA1[i]
      lncRNA2<-jac2$lncRNA2[i]
      celltype<-jac2$cell_type[i]
      fun_in_lncRNA1<-fun[which(fun$lncRNA == lncRNA1&fun$cell_type == celltype),]
      NES1<-as.numeric(fun_in_lncRNA1$NES)
      names(NES1)<-fun_in_lncRNA1$ID
      fun_in_lncRNA2<-fun[which(fun$lncRNA == lncRNA2&fun$cell_type == celltype),]
      NES2<-as.numeric(fun_in_lncRNA2$NES)
      names(NES2)<-fun_in_lncRNA2$ID
      interfun<-intersect(fun_in_lncRNA1$ID,fun_in_lncRNA2$ID)##The overlapping functions between two lncrnas
      ##Calculate the co-regulation coefficient
      a<-NES1[interfun]
      b<-NES2[interfun]
      s<-ifelse(sign(a)==sign(b),1,-1)
      x<-s*((abs(a)+abs(b))/2)
      X<-sum(x*w[interfun])/sum(w[interfun])
      
      coindex<-c(coindex,X)
    }
    jac2$coindex<-coindex
    jac2$type<-names(func)[f]
    jac2$tumor<-tumor
    jacmatrix2[[f]]<-jac2
  }
  co_lncRNA<-do.call(rbind,jacmatrix2)
  write.table(co_lncRNA,"co-lncRNA.txt",sep = "\t",quote = F,col.names = T,row.names = F)
}








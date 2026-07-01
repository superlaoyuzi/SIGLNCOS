

# ===== From: 整体识别流程.R =====

##导包
library(Seurat)##V4
library(dplyr)
library(psych)##计算相关性
library(clusterProfiler)##功能富集
##全局变量
tumor<-"HNSC"##癌症类型
setwd("I:/pan-cancer lncRNA课题/数据库扩展/test")##路径
##读取注释后单细胞数据
load("./GSE188737_hnscc_seu.RData/GSE188737_hnscc_seu.RData")
SO<-UpdateSeuratObject(hnscc_seu)

##差异表达识别clncRNA
ftq<-readRDS("./ftq.rds")
lncRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "lncRNA")%>%arrange("gene_name")
lncRNA<-intersect(lncRNA$gene_name,rownames(SO))##lncRNA
mRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "protein_coding")%>%arrange("gene_name")
mRNA<-intersect(mRNA$gene_name,rownames(SO))##mRNA
Idents(SO)<-SO$cell_type##根据细胞类型计算
SO.markers<-FindAllMarkers(SO,only.pos = T,)##仅选择高表达数据(avglog2FC>0.25,min.pct = 0.1)
SO.markers$tumor<-tumor
SO.markers<-subset(SO.markers,p_val_adj<0.05)##选择校正后p值显著的gene
SO.markers$tumor<-tumor
write.table(SO.markers,"all_DEG.txt",sep = "\t",quote = F,col.names = T,row.names = F)
lnc.markers<-SO.markers[which(SO.markers$gene %in% lncRNA),]##差异高表达的lncRNA
write.table(lnc.markers,"DElnc.txt",sep = "\t",quote = F,col.names = T,row.names = F)
mRNA.markers<-SO.markers[which(SO.markers$gene %in% mRNA),]##差异高表达的mRNA
##lncRNA功能调控分析
##根据lncRNA计算mRNA相关性（spearman）
clncRNA<-unique(lnc.markers$gene)
exp.data<-GetAssayData(object = SO, layer  = "data")##提取表达谱(5.0版本换成slot = "data")
clnc.exp<-exp.data[which(rownames(exp.data) %in% clncRNA),]##clncRNA表达矩阵
mRNA.exp<-exp.data[which(rownames(exp.data) %in% mRNA),]##mRNA表达矩阵
options(warn=0)
cor_data<-lapply(clncRNA,function(x){
  cell<-as.vector(subset(lnc.markers,gene == x)$cluster)
  aa<-lapply(cell,function(z){
    cellid<-colnames(subset(SO,cell_type == z))
    clnc.exp2<-clnc.exp[,cellid]
    mRNA.exp2<-mRNA.exp[,cellid]
    cmRNA<-subset(mRNA.markers,cluster == z)$gene
    bb<-lapply(mRNA, function(y){##mRNA为全部mRNA、cmRNA为差异表达的mRNA
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
##根据排序的mRNA进行GSEA富集分析
##读取富集用到的功能
GO<-read.gmt("./function_old/c5.go.v2022.1.Hs.symbols.gmt")
imm<-read.gmt("./function_old/c7.all.v2022.1.Hs.symbols.gmt")
kegg<-read.gmt("./function_old/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
##功能注释
cor_data$cor<-as.numeric(cor_data$cor)
cor_data$p.adj<-as.numeric(cor_data$p.adj)##转换成数值类型
cell<-as.vector(unique(cor_data$cell_type))
lnc_function<-c()
for(i in cell){
  subcor<-subset(cor_data,cell_type == i)
  clncRNA2<-unique(subcor$lncRNA)
  celldata<-lapply(clncRNA2,function(x){
    clncRNA_cor<-subset(subcor,lncRNA == x)$cor
    clncRNA_p<-subset(cor_data,lncRNA == x)$p.adj
    ##padj中的0值替换为最小值
    minp<-min(clncRNA_p[-which(clncRNA_p==0)])
    clncRNA_p<-ifelse(clncRNA_p == 0,minp,clncRNA_p)
    Rank<- -log(as.numeric(clncRNA_p))*sign(as.numeric(clncRNA_cor))
    names(Rank)<-subset(subcor,lncRNA == x)$mRNA
    Rank<-sort(Rank,decreasing = T)##排序
    ##GSEA
    gseaGO<-GSEA(Rank ,TERM2GENE =  GO)
    gseaimmune<-GSEA(Rank,TERM2GENE = imm)
    gseaKEGG<-GSEA(Rank,TERM2GENE = kegg)
    ##选取显著结果
    GOresult<-gseaGO@result[which(gseaGO@result$p.adjust<0.05),]
    immgsearesult<-gseaimmune@result[which(gseaimmune@result$p.adjust<0.05),]
    keggsearesult<-gseaKEGG@result[which(gseaKEGG@result$p.adjust<0.05),]
    ##添加lncRNA信息
    if(nrow(GOresult)>0){
      GOresult$lncRNA<-x
    }
    if(nrow(immgsearesult)>0){
      immgsearesult$lncRNA<-x
    }
    if(nrow(keggsearesult)>0){
      keggsearesult$lncRNA<-x
    }
    ##添加cell信息
    if(nrow(GOresult)>0){
      GOresult$cell_type<-i
    }
    if(nrow(immgsearesult)>0){
      immgsearesult$cell_type<-i
    }
    if(nrow(keggsearesult)>0){
      keggsearesult$cell_type<-i
    }
    ##添加分类标签
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

##co-lncRNA识别
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
##计算lncRNA之间的共调控系数
jacmatrix2<-list()
for(f in 1:3){##根据三种功能计算共调控系数
  fun<-func[[f]]##lncRNA对该类功能的调控情况
  w<-W[[f]]##该类功能的权
  jac<-jacmatrix[[f]]##基于该类功能计算的lncRNA之间的jaccard系数
  jac2<-jac[which(jac$jaccord>=0.5),]##根据jaccard过滤lncRNA关系对
  if(nrow(jac2)<1){next}
  coindex<-c()##每对lncRNA的共调控系数
  for(i in 1:nrow(jac2)){##对于每个lncRNA关系对
    lncRNA1<-jac2$lncRNA1[i]
    lncRNA2<-jac2$lncRNA2[i]
    celltype<-jac2$cell_type[i]
    fun_in_lncRNA1<-fun[which(fun$lncRNA == lncRNA1&fun$cell_type == celltype),]
    NES1<-as.numeric(fun_in_lncRNA1$NES)
    names(NES1)<-fun_in_lncRNA1$ID
    fun_in_lncRNA2<-fun[which(fun$lncRNA == lncRNA2&fun$cell_type == celltype),]
    NES2<-as.numeric(fun_in_lncRNA2$NES)
    names(NES2)<-fun_in_lncRNA2$ID
    interfun<-intersect(fun_in_lncRNA1$ID,fun_in_lncRNA2$ID)##两个lncRNA之间重叠的功能
    ##计算共调控系数
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













# ===== From: 识别colnc.R =====

##导包
library(Seurat)##V4
library(dplyr)
library(psych)##计算相关性
library(clusterProfiler)##功能富集
##全局变量
setwd("I:/pan-cancer lncRNA课题/数据库扩展/test")##路径
##读取注释后单细胞数据

SO<-UpdateSeuratObject(hnscc_seu)

lnc_function<-read.table("I:/pan-cancer lncRNA课题/数据库扩展/all_lnc_function.txt",sep = "\t",header = T)
##co-lncRNA识别
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
##计算lncRNA之间的共调控系数
jacmatrix2<-list()
for(f in 1:3){##根据三种功能计算共调控系数
  fun<-func[[f]]##lncRNA对该类功能的调控情况
  w<-W[[f]]##该类功能的权
  jac<-jacmatrix[[f]]##基于该类功能计算的lncRNA之间的jaccard系数
  jac2<-jac[which(jac$jaccord>=0.5),]##根据jaccard过滤lncRNA关系对
  if(nrow(jac2)<1){next}
  coindex<-c()##每对lncRNA的共调控系数
  for(i in 1:nrow(jac2)){##对于每个lncRNA关系对
    lncRNA1<-jac2$lncRNA1[i]
    lncRNA2<-jac2$lncRNA2[i]
    celltype<-jac2$cell_type[i]
    fun_in_lncRNA1<-fun[which(fun$lncRNA == lncRNA1&fun$cell_type == celltype),]
    NES1<-as.numeric(fun_in_lncRNA1$NES)
    names(NES1)<-fun_in_lncRNA1$ID
    fun_in_lncRNA2<-fun[which(fun$lncRNA == lncRNA2&fun$cell_type == celltype),]
    NES2<-as.numeric(fun_in_lncRNA2$NES)
    names(NES2)<-fun_in_lncRNA2$ID
    interfun<-intersect(fun_in_lncRNA1$ID,fun_in_lncRNA2$ID)##两个lncRNA之间重叠的功能
    ##计算共调控系数
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




func<-list(
  "GO" = subset(lnc_function,func_type == "GO"),
  "IMMUNE" = subset(lnc_function,func_type == "IMMUNE"),
  "KEGG" = subset(lnc_function,func_type == "KEGG")
)
file<-names(func)
##对于每个功能计算jaccard系数
{
  jacmatrix<-list()
  for(f in 1:3){##分别基于GO、KEGG以及免疫计算jaccard系数
    fun<-func[[f]]##读取功能信息
    ##所有癌症
    tumor<-unique(fun$tumor)
    tumorline<-c()##对肿瘤构建矩阵
    ##对所有癌症
    ##计算lncRNA之间功能的overlap（杰卡尔系数）
    for(i in tumor){##对所有癌症
      tumordata<-fun[which(fun$tumor == i),]
      celltype<-unique(tumordata$cell_type)
      cellline<-c()##对细胞构建矩阵
      for(j in celltype){##对所有细胞类型
        celldata<-tumordata[which(tumordata$cell_type == j),]
        lncRNA<-unique(celldata$lncRNA)
        if(length(lncRNA)<2){##去除lncRNA少于两个的情况
          next
        }
        lnclinex<-c()##对lncRNA x 构建矩阵
        for(x in 1:(length(lncRNA)-1)){##对于lncRNA x 的功能信息
          funcx<-unique(celldata[which(celldata$lncRNA==lncRNA[x]),"ID"])
          lncliney<-c()##对lncRNA y构建矩阵
          for(y in (x+1):length(lncRNA)){##对于lncRNA y 的功能信息
            funcy<-unique(celldata[which(celldata$lncRNA==lncRNA[y]),"ID"])
            inter<-length(intersect(funcx,funcy))
            uni<-length(union(funcx,funcy))
            jaccord<-inter/uni##计算overlap（jaccord）
            lncliney<-rbind(lncliney,c(lncRNA[x],lncRNA[y],jaccord,j,i))
          }
          lnclinex<-rbind(lnclinex,lncliney)
        }
        cellline<-rbind(cellline,lnclinex)
      }
      tumorline<-rbind(tumorline,cellline)
    }
    colnames(tumorline)<-c("lncRNA1","lncRNA2","jaccord","cell_type","tumor")
    rownames(tumorline)<-1:nrow(tumorline)
    jacmatrix[[f]]<-as.data.frame(tumorline)
    
  }
  names(jacmatrix)<-file
}




##基于加权平均的思想识别co-lncRNA
{
  ##对于每一种功能，根据lncRNA调控情况计算其权重
  W<-list()
  for(f in file){
    fun<-func[[f]]##依次对三种功能计算权重
    lncnum<-length(unique(fun$lncRNA))##总的lncRNA数量
    fun2<-split(fun$lncRNA,fun$ID)##根据功能id分离lncRNA
    fun_lncnum<-unlist(lapply(fun2,function(x){length(unique(x))}))##调控每个功能的lncRNA数量
    w<- -log10(fun_lncnum/lncnum)##每个功能的权重
    W[[f]]<-w
  }
  ##计算lncRNA之间的共调控系数
  jacmatrix2<-list()
  for(f in file){##根据三种功能计算共调控系数
    fun<-func[[f]]##lncRNA对该类功能的调控情况
    w<-W[[f]]##该类功能的权
    jac<-jacmatrix[[f]]##基于该类功能计算的lncRNA之间的jaccard系数
    jac2<-jac##根据jaccard过滤lncRNA关系对
    coindex<-c()##每对lncRNA的共调控系数
    for(i in 1:nrow(jac2)){##对于每个lncRNA关系对
      lncRNA1<-jac2$lncRNA1[i]
      lncRNA2<-jac2$lncRNA2[i]
      celltype<-jac2$cell_type[i]
      tumor<-jac2$tumor[i]
      fun_in_lncRNA1<-fun[which(fun$lncRNA == lncRNA1&fun$cell_type == celltype&fun$tumor == tumor),]
      NES1<-as.numeric(fun_in_lncRNA1$NES)
      names(NES1)<-fun_in_lncRNA1$ID
      fun_in_lncRNA2<-fun[which(fun$lncRNA == lncRNA2&fun$cell_type == celltype&fun$tumor == tumor),]
      NES2<-as.numeric(fun_in_lncRNA2$NES)
      names(NES2)<-fun_in_lncRNA2$ID
      interfun<-intersect(fun_in_lncRNA1$ID,fun_in_lncRNA2$ID)##两个lncRNA之间重叠的功能
      ##计算共调控系数
      a<-NES1[interfun]
      b<-NES2[interfun]
      s<-ifelse(sign(a)==sign(b),1,-1)
      x<-s*((abs(a)+abs(b))/2)
      X<-sum(x*w[interfun])/sum(w[interfun])
      
      coindex<-c(coindex,X)
    }
    jac2$coindex<-coindex
    jac2$func_type<-f
    jacmatrix2[[f]]<-jac2
  }
}

result<-do.call(rbind,jacmatrix2)
write.csv(result,"I:/pan-cancer lncRNA课题/数据库扩展/colnc_no_filter.csv",quote = F,row.names = F)





# ===== From: 识别colnc脚本.R =====

##脚本
library(Seurat)
library(dplyr)
library(psych)##计算相关性
library(clusterProfiler)##功能富集
setwd("I:/pan-cancer lncRNA课题/数据库扩展/test")##上传的路径
data<-data.frame(data.table::fread("./GSM3770932_SyS.tumors10x_counts.csv.gz",header = T),row.names = 1)##读取表达谱
cell_anno<-data.frame(data.table::fread("./GSM3770932_SyS.tumors10x_cell.annotations.csv"),row.names = 1)##读取注释
##转换为Seurat对象
SO<-CreateSeuratObject(data)
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")
SO <- subset(SO, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
SO <- NormalizeData(SO)
SO <- FindVariableFeatures(SO, selection.method = "vst", nfeatures = 2000)
SO <- ScaleData(SO, features = rownames(SO))
SO <- RunPCA(SO, features = VariableFeatures(object = SO),reduction.name = "pca")
SO <- FindNeighbors(SO, reduction = "pca", dims = 1:30)
SO <- FindClusters(SO, resolution = 1)

##注释细胞类型
rownames(cell_anno)<-cell_anno[,1]
SO$cell_type<-cell_anno[colnames(SO),2]

##分析
##差异表达识别clncRNA
dir.create("result")
ftq<-readRDS("./context/ftq.rds")
lncRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "lncRNA")%>%arrange("gene_name")
lncRNA<-intersect(lncRNA$gene_name,rownames(SO))##lncRNA
mRNA<-ftq%>%dplyr::select(c("type","gene_type","gene_name"))%>%filter(type == "gene" & gene_type == "protein_coding")%>%arrange("gene_name")
mRNA<-intersect(mRNA$gene_name,rownames(SO))##mRNA
Idents(SO)<-SO$cell_type##根据细胞类型计算
SO.markers<-FindAllMarkers(SO,only.pos = T,)##仅选择高表达数据(avglog2FC>0.25,min.pct = 0.1)
SO.markers<-subset(SO.markers,p_val_adj<0.05)##选择校正后p值显著的gene
#write.table(SO.markers,"./result/all_DEG.txt",sep = "\t",quote = F,col.names = T,row.names = F)
lnc.markers<-SO.markers[which(SO.markers$gene %in% lncRNA),]##差异高表达的lncRNA
#write.table(lnc.markers,"./result/DElnc.txt",sep = "\t",quote = F,col.names = T,row.names = F)
mRNA.markers<-SO.markers[which(SO.markers$gene %in% mRNA),]##差异高表达的mRNA
##lncRNA功能调控分析
##根据lncRNA计算mRNA相关性（spearman）
clncRNA<-unique(lnc.markers$gene)
exp.data<-GetAssayData(object = SO, layer  = "data")##提取表达谱(5.0版本换成slot = "data")
clnc.exp<-exp.data[which(rownames(exp.data) %in% clncRNA),]##clncRNA表达矩阵
mRNA.exp<-exp.data[which(rownames(exp.data) %in% mRNA),]##mRNA表达矩阵
options(warn=0)
cor_data<-lapply(clncRNA,function(x){
  cell<-as.vector(subset(lnc.markers,gene == x)$cluster)
  aa<-lapply(cell,function(z){
    cellid<-colnames(subset(SO,cell_type == z))
    clnc.exp2<-clnc.exp[,cellid]
    mRNA.exp2<-mRNA.exp[,cellid]
    cmRNA<-subset(mRNA.markers,cluster == z)$gene
    bb<-lapply(mRNA, function(y){##mRNA为全部mRNA、cmRNA为差异表达的mRNA
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
#write.table(cor_data,"./result/lncRNA_mRNA_cor.txt",sep = "\t",quote = F,col.names = T,row.names = F)
##根据排序的mRNA进行GSEA富集分析
##读取富集用到的功能
GO<-read.gmt("./context/c5.go.v2022.1.Hs.symbols.gmt")
imm<-read.gmt("./context/c7.all.v2022.1.Hs.symbols.gmt")
kegg<-read.gmt("./context/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
##功能注释
cor_data$cor<-as.numeric(cor_data$cor)
cor_data$p.adj<-as.numeric(cor_data$p.adj)##转换成数值类型
cell<-as.vector(unique(cor_data$cell_type))
lnc_function<-c()
for(i in cell){
  subcor<-subset(cor_data,cell_type == i)
  clncRNA2<-unique(subcor$lncRNA)
  celldata<-lapply(clncRNA2,function(x){
    clncRNA_cor<-subset(subcor,lncRNA == x)$cor
    clncRNA_p<-subset(cor_data,lncRNA == x)$p.adj
    ##padj中的0值替换为最小值
    minp<-min(clncRNA_p[-which(clncRNA_p==0)])
    clncRNA_p<-ifelse(clncRNA_p == 0,minp,clncRNA_p)
    Rank<- -log(as.numeric(clncRNA_p))*sign(as.numeric(clncRNA_cor))
    names(Rank)<-subset(subcor,lncRNA == x)$mRNA
    Rank<-sort(Rank,decreasing = T)##排序
    ##GSEA
    gseaGO<-GSEA(Rank ,TERM2GENE =  GO)
    gseaimmune<-GSEA(Rank,TERM2GENE = imm)
    gseaKEGG<-GSEA(Rank,TERM2GENE = kegg)
    ##选取显著结果
    GOresult<-gseaGO@result[which(gseaGO@result$p.adjust<0.05),]
    immgsearesult<-gseaimmune@result[which(gseaimmune@result$p.adjust<0.05),]
    keggsearesult<-gseaKEGG@result[which(gseaKEGG@result$p.adjust<0.05),]
    ##添加lncRNA信息
    if(nrow(GOresult)>0){
      GOresult$lncRNA<-x
    }
    if(nrow(immgsearesult)>0){
      immgsearesult$lncRNA<-x
    }
    if(nrow(keggsearesult)>0){
      keggsearesult$lncRNA<-x
    }
    ##添加cell信息
    if(nrow(GOresult)>0){
      GOresult$cell_type<-i
    }
    if(nrow(immgsearesult)>0){
      immgsearesult$cell_type<-i
    }
    if(nrow(keggsearesult)>0){
      keggsearesult$cell_type<-i
    }
    ##添加分类标签
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
#write.table(lnc_function,"lncRNA_function_regulation.txt",sep = "\t",quote = F,col.names = T,row.names = F)

##co-lncRNA识别
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
##计算lncRNA之间的共调控系数
jacmatrix2<-list()
for(f in 1:3){##根据三种功能计算共调控系数
  fun<-func[[f]]##lncRNA对该类功能的调控情况
  w<-W[[f]]##该类功能的权
  jac<-jacmatrix[[f]]##基于该类功能计算的lncRNA之间的jaccard系数
  jac2<-jac[which(jac$jaccord>=0.5),]##根据jaccard过滤lncRNA关系对
  if(nrow(jac2)<1){next}
  coindex<-c()##每对lncRNA的共调控系数
  for(i in 1:nrow(jac2)){##对于每个lncRNA关系对
    lncRNA1<-jac2$lncRNA1[i]
    lncRNA2<-jac2$lncRNA2[i]
    celltype<-jac2$cell_type[i]
    fun_in_lncRNA1<-fun[which(fun$lncRNA == lncRNA1&fun$cell_type == celltype),]
    NES1<-as.numeric(fun_in_lncRNA1$NES)
    names(NES1)<-fun_in_lncRNA1$ID
    fun_in_lncRNA2<-fun[which(fun$lncRNA == lncRNA2&fun$cell_type == celltype),]
    NES2<-as.numeric(fun_in_lncRNA2$NES)
    names(NES2)<-fun_in_lncRNA2$ID
    interfun<-intersect(fun_in_lncRNA1$ID,fun_in_lncRNA2$ID)##两个lncRNA之间重叠的功能
    ##计算共调控系数
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

aa<-apply(co_lncRNA,1,function(x){
  cell<-unlist(c(x[4]))
  lncRNA1<-unlist(c(x[1]))
  lncRNA2<-unlist(c(x[2]))
  functype<-unlist(c(x[7]))
  fun1<-unique(subset(newlncfun,lncRNA == lncRNA1 & cell_type == cell & func_type == functype)[,3])
  fun2<-unique(subset(newlncfun,lncRNA == lncRNA2 & cell_type == cell & func_type == functype)[,3])
  func<-intersect(fun1,fun2)
  return(paste(func,collapse = ";"))
})
co_lncRNA$func<-aa
write.table(co_lncRNA,"./result/co-lncRNA.txt",sep = "\t",quote = F,col.names = T,row.names = F)

# ===== From: 功能富集.R =====

##library packages
{
  library(Seurat)
  library(clusterProfiler)
  library(dplyr)
  library(enrichplot)
}

##gene注释(Gencode)
{
  ftq<-readRDS("F:/课题/context/ftq.rds")
  lncRNA<-unique(ftq$gene_name[which(ftq$gene_type == "lncRNA")])
  mRNA<-unique(ftq$gene_name[which(ftq$gene_type == "protein_coding")])
  
}

##泛癌路径和名称
{
  pan_cancer_file<-c("ATC/ATC1","BC/BC1","CRC/CRC1_Tumor","ESCC/ESCC1","GBM/GBM2","GC/GC1","HCC/HCC2","HNSCC/HNSCC1","LUAD/LUAD2","MCC","NBL","OV/OV1","STAD","TNBC/TNBC4","UCEC/EC1")
  name<-c("ATC","BC","CRC","ESCC","GBM","GC","HCC","HNSCC","LUAD","MCC","NBL","OV","STAD","TNBC","UCEC")
  
}

##注释基因集
{
  GO<-read.gmt("F:/课题/context/c5.go.v2022.1.Hs.symbols.gmt")
  ##hallmark注释
  ##hmt<-read.gmt("E:/work/课题/gsea/h.all.v7.5.1.entrez.gmt")
  #免疫相关注释
  immmt<-read.gmt("F:/课题/context/c7.all.v2022.1.Hs.symbols.gmt")
  ##KEGG注释
  keggmt<-read.gmt("F:/课题/context/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
  
}
##差异表达
{
  setwd("F:/20221016试验记录")
  alldegene<-c()
  for(a in 1:length(pan_cancer_file)){
    SO<-readRDS(paste("F:/20220314试验记录/",pan_cancer_file[a],"/SO.rds",sep = ""))
    setwd(paste("F:/20221016试验记录/",name[a],sep = ""))
    celltype<-unique(SO@meta.data$celltype)
    allcelltype<-SO@meta.data$celltype
    ##差异表达分析
    tumordegene<-c()
    for(i in 1:length(celltype)){
      degene<-FindMarkers(SO,group.by = "celltype",ident.1 = celltype[i],only.pos = T)
      degene$cluster<-rep(celltype[i],nrow(degene))
      degene$gene<-rownames(degene)
      tumordegene<-rbind(tumordegene,degene)
    }
    tumordegene<-cbind(tumordegene,tumor = name[a])
    alldegene<-rbind(alldegene,tumordegene)
  }
  write.csv(alldegene,"alldegene.txt",sep = "\t",quote = F,row.names = F)
  alldegene<-read.csv2("F:/20221016试验记录/alldegene.txt",header = T,sep = ",")
}
##spearman相关性分析
  alldegene<-subset(alldegene,as.numeric(alldegene$p_val_adj)<0.05)
  tumorresult<-c()
  setwd("F:/20221016试验记录")
  for(i in 1:15){##每一种癌症
    ##读取该癌症的表达谱
    SO<-readRDS(paste("F:/20220314试验记录/",pan_cancer_file[i],"/SO.rds",sep = ""))
    expdata<-SO@assays$RNA@data
    ##差异表达矩阵
    tumordata<-subset(alldegene,alldegene$tumor == name[i])
    celltype<-unique(tumordata$cluster)
    cellresult<-c()
    for(j in 1:length(celltype)){##每种细胞类型
      celldata<-subset(tumordata,tumordata$cluster == celltype[j])
      ##RNA注释
      cellRNA<-unique(celldata$gene)
      celllncRNA<-intersect(cellRNA,lncRNA)
      cellmRNA<-intersect(cellRNA,mRNA)
      if(length(celllncRNA) == 0){##确保存在lncRNA
        next
      }
      ##该细胞类型表达谱
      cellexpdata<-as.matrix(expdata[,which(SO$celltype == celltype[j])])
      lncexpdata<-cellexpdata[which(rownames(cellexpdata)%in%celllncRNA),]
      mexpdata<-cellexpdata[which(rownames(cellexpdata)%in%cellmRNA),]
      if(length(celllncRNA) == 1){
        lncexpdata<-t(as.matrix(lncexpdata))
        rownames(lncexpdata)<-celllncRNA
      }
      ##计算矩阵之间相关性(spearman)
      result<-as.data.frame(t(data.frame(apply(lncexpdata, 1, function(x){
        data.frame(
          apply(mexpdata, 1, function(y){
            co<-cor.test(as.numeric(x),as.numeric(y),method = "spearman")
            return(c(co$estimate,co$p.value))
          })
        )
      }))
      ))
      m<-rep(rownames(mexpdata),nrow(lncexpdata))
      lnc<-c()
      for (n in rownames(lncexpdata)) {
        lnc<-c(lnc,rep(n,nrow(mexpdata)))
      }
      fdr<-p.adjust(result[,2],method = "fdr")
      rownames(result)<-1:nrow(result)
      result$lncRNA<-lnc
      result$mRNA<-m
      result$fdr<-fdr
      result<-cbind(result,celltype = celltype[j])
      cellresult<-rbind(cellresult,result)
    }
    cellresult<-cbind(cellresult,tumor = name[i])
    tumorresult<-rbind(tumorresult,cellresult)
  }
  write.csv2(tumorresult,"lncRNA_mRNA_cor.txt",quote = F,row.names = F,sep = "\t")
  
  tumorresult<-read.csv2("F:/20221016试验记录/lncRNA_mRNA_cor.txt",header = T,sep = "\t")
  ##GSEA
  {
    allGO<-c()
    allKEGG<-c()
    allImmune<-c()
    for(i in 1:15){##对每个癌症
      tumordata<-subset(tumorresult,tumorresult$tumor == name[i])
      celltype<-unique(tumordata$celltype)
      cellGOresult<-cellKEGGresult<-cellimmresult<-c()
      
      for(j in 1:length(celltype)){##对每个细胞类型
        celldata<-subset(tumordata,tumordata$celltype == celltype[j])
        ##lncRNA
        celllnc<-unique(celldata$lncRNA)
        xGOresult<-ximmuneresult<-xKEGGresult<-c()
        for(x in 1:length(celllnc)){##对每一个lncRNA
          xdata<-subset(celldata[,c(1,4,5)],celldata$lncRNA == celllnc[x])
          xdata$score<- -log(as.numeric(xdata$fdr))*sign(as.numeric(xdata$rho))
          xdata_rank<-xdata[order(xdata$score,decreasing = T),c(2,4)]
          colnames(xdata_rank)<-c("SYMBOL","RANK")
          genelist<-xdata_rank$RANK
          names(genelist)<-xdata_rank$SYMBOL
          genelist<-sort(genelist,decreasing = T)
          gseaGO<-GSEA(genelist ,TERM2GENE =  GO)
          gseaimmune<-GSEA(genelist,TERM2GENE = immmt)
          gseaKEGG<-GSEA(genelist,TERM2GENE = keggmt)
          GOresult<-gseaGO@result[which(gseaGO@result$p.adjust<0.05),]
          if(nrow(GOresult) > 0){
            GOresult$lncRNA<-rep(celllnc[x],nrow(GOresult))
            xGOresult<-rbind(xGOresult,GOresult)
          }
          
          immgsearesult<-gseaimmune@result[which(gseaimmune@result$p.adjust<0.05),]
          if(nrow(immgsearesult) > 0){
            immgsearesult$lncRNA<-rep(celllnc[x],nrow(immgsearesult))
            ximmuneresult<-rbind(ximmuneresult,immgsearesult)
          }
          keggsearesult<-gseaKEGG@result[which(gseaKEGG@result$p.adjust<0.05),]
          if(nrow(keggsearesult) > 0){
            keggsearesult$lncRNA<-rep(celllnc[x],nrow(keggsearesult))
            xKEGGresult<-rbind(xKEGGresult,keggsearesult)
          }
          
        }
        xGOresult<-cbind(xGOresult,celltype = celltype[j])
        xKEGGresult<-cbind(xKEGGresult,celltype = celltype[j])
        ximmuneresult<-cbind(ximmuneresult,celltype = celltype[j])
        cellGOresult<-rbind(cellGOresult,xGOresult)
        cellKEGGresult<-rbind(cellKEGGresult,xKEGGresult)
        cellimmresult<-rbind(cellimmresult,ximmuneresult)
      }
      cellGOresult<-cbind(cellGOresult,tumor = name[i])
      cellKEGGresult<-cbind(cellKEGGresult,tumor = name[i])
      cellimmresult<-cbind(cellimmresult,tumor = name[i])
      allGO<-rbind(allGO,cellGOresult)
      allKEGG<-rbind(allKEGG,cellKEGGresult)
      allImmune<-rbind(allImmune,cellimmresult)
    }
    write.table(allGO,"GO_gsea.txt",quote = F,row.names = F,sep = "\t")
    write.table(allImmune,"immune_gsea.txt",quote = F,row.names = F,sep = "\t")
    write.table(allKEGG,"kegg_gsea.txt",quote = F,row.names = F,sep = "\t")
  }
  
  
     
    
    


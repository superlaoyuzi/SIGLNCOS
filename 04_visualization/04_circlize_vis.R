library(circlize)
library(dplyr)
setwd("I:/pan-cancer lncRNA课题/数据库扩展")

name<-c("ATC","BC","CRC","ESCC","GBM","GC","HCC","HNSCC","LUAD","MCC","NBL","OV","TNBC","UCEC",
        "STAD","BLCA","cSCC","ICC","LSCC","meningioma","OCCC","OS","OSCC","Pan-NET","PDAC","PRAD","PTC","RB","RMS","SyS","UVM")
colors<-c("#3b6291","#943c39","#779043","#624c7c","#388498","#bf7334","#3f6899","#9c403d","#7d9847","#675083","#3b8ba1","#c97937","#7FFFD4","#FFA500",
          "#FDF6DC","#8D5D47","#BD7F57","#4A6D8C","#2E4559","#CC3A40","#D4AFA5","#77967B","#730217","#D99441","#686973","#A6896F","#CFD9C5","#A37792","#F2BF80","#152937","#011640")

names(colors)<-name
delncRNA<-read.table("alldelncRNA2.txt",sep = "\t",header = T)[,-10]
delncRNA<-delncRNA%>%group_by(tumor)%>%top_n(30,avg_log2FC)
ftq<-readRDS("I:/pan-cancer lncRNA课题/context/ftq.rds")
chrominfo<-ftq%>%dplyr::select(c("seqnames","start","end","type","gene_name"))%>%filter(type == "gene")
chrominfo<-chrominfo[!duplicated(chrominfo$gene_name),]
##构建染色体
human_chromInfo = read.chromInfo(species = "hg38")$df
tumor<-unique(delncRNA$tumor)
tumorinfo<-c("tumor",0,1200000000)
chromInfo = rbind(human_chromInfo, tumorinfo)
chromosome.index = c(paste0("chr", c(1:22, "X","Y")),"tumor")
##添加基因名
lncRNA<-unique(delncRNA$gene)
lncRNA_index<-subset(chrominfo,gene_name %in% lncRNA)
lncRNA_index$seqnames<-as.vector(lncRNA_index$seqnames)
tumorprop<-table(delncRNA$tumor)
tumorprop2<-(tumorprop/sum(tumorprop))*1200000000
tumorlength<-tumorprop2[tumor[1]]
gap<-1200000000*0.001
for(i in 1:length(tumor)){
  if(i == 1){
    lncRNA_index<-rbind(lncRNA_index,c("tumor",gap,tumorlength,"tumor",tumor[i]))
  }else{
    lncRNA_index<-rbind(lncRNA_index,c("tumor",tumorlength2+gap,tumorlength,"tumor",tumor[i]))
  }
  tumorlength2<-tumorlength
  tumorlength<-tumorlength+tumorprop2[tumor[i+1]]
  
}
lncRNA_index$start<-as.numeric(lncRNA_index$start)
lncRNA_index$end<-as.numeric(lncRNA_index$end)
##添加连线
link<-lncRNA_index
rownames(link)<-lncRNA_index$gene_name
link1<-link[delncRNA$gene,]
link2<-link[delncRNA$tumor,]


##绘图
pdf(paste(tumor,"circos.pdf",sep = "_"),width = 50,height = 50)
circos.clear()
circos.initializeWithIdeogram(chromInfo, chromosome.index = chromosome.index,plotType = c("axis", "labels"))
circos.genomicLabels(lncRNA_index, labels.column = 5, side = "inside",cex = 0.8,connection_height = mm_h(11))
circos.genomicLink(link1, link2, col = colors[link2$gene_name], 
                   border = NA)
dev.off()

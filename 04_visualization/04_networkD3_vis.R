library(networkD3)
colnc<-read.table("I:/pan-cancer lncRNA课题/数据库扩展/all_colncRNA.txt",sep = "\t",header = T)
tumor<-unique(colnc$tumor)
networkData <- aggregate(colnc$jaccord,list(colnc$lncRNA1,colnc$lncRNA2),sum)
colnames(networkData)<-c("lncRNA1","lncRNA2","weight")
nodeData <- data.frame(name=unique(c(networkData$lncRNA1, networkData$lncRNA2)),group = "lncRNA")
linkData <- data.frame(source = (match(networkData$lncRNA1, nodeData$name)-1),
                       target = (match(networkData$lncRNA2, nodeData$name)-1),
                       value = networkData$weight)
nodeData$size<-as.character(node[nodeData$name,3])
nodeData$size<-nodeData$size/sum(nodeData$size)
forceNetwork(
  Links = linkData,
  Nodes = nodeData,
  Source = "source",
  Target = "target",
  Value = "value",
  Nodesize = 'size',
  NodeID = "name",Group = "group",
  opacity = 0.8, zoom = T, bounded = T,
  clickAction = MyClickScript,
  charge=-80,linkColour = "lightgrey",fontSize = 16,
  opacityNoHover = 0)
MyClickScript <- 
  '      d3.select(this).select("circle").transition()
.duration(750)
.attr("r", 30)
.attr("color", "red")'
network <- htmlwidgets::onRender(network, '
  function(el, x) {
    var svg = d3.select(el).select("svg");
    var width = +svg.attr("width");
    var height = +svg.attr("height");

    // 使用力导向布局并限制迭代次数
    d3.forceSimulation(x.nodes)
      .force("link", d3.forceLink(x.links).id(function(d) { return d.index; }))
      .force("charge", d3.forceManyBody().strength(-30))
      .force("center", d3.forceCenter(width / 2, height / 2))
      .alphaDecay(0.05) // 调整衰减率来限制迭代次数
      .on("tick", function() {
        d3.select(el).selectAll(".link")
          .attr("x1", function(d) { return d.source.x; })
          .attr("y1", function(d) { return d.source.y; })
          .attr("x2", function(d) { return d.target.x; })
          .attr("y2", function(d) { return d.target.y; });

        d3.select(el).selectAll(".node")
          .attr("cx", function(d) { return d.x; })
          .attr("cy", function(d) { return d.y; });
      })
      .on("end", function() {
        // 在布局结束后固定节点位置以减少动画计算
        x.nodes.forEach(function(d) {
          d.fx = d.x;
          d.fy = d.y;
        });
      });
  }
')


network
saveNetwork(network, "network_with_click_action.html")



##lncRNA_function
color <- 'd3.scaleOrdinal() .domain(["Tumor cell","Epithelial cell","Endothelial cell","Mesenchymal cell","Fibroblast","T cell","CD4+T cell","CD8+T cell","NK","NKT","B cell","Plasma cell","Macrophage","Monocyte","Melanocyte","MDM","Megakaryocyte","DC","CMP","Unassign","Oligodendrocyte","Mast cell","Cancer stem cells","Keratinocytes","Pericyte","Neutrophils","Meninge derived cell","Keratinocytes","Pericyte","Alpha Cell","Beta Cell","Ductal Cell","Pancreatic Stellate Cell","Proliferating Beta cell","Stellate","Basal Intermediate cell","Luminal", "Neuron cell","Neutrophil","Smooth Muscle cell","Amacrine cell","Cone cell", "Melanocyte cell","Microglia cell","Retinal ganglion cell" ]) 
.range(["#8FC1B5","#A6E2D1","#47A1BA","#9DC4D3","#D9E7EE","#C8A4C4","#D8A7B1","#E8C4C6","#FFE4E1"
,"#FFB6C1","#FF9A82","#F98080","#CC9E7C","#D1AD8E","#799956","#E3C6A6","#F2D6C0","#F5DEBA"
,"#E7E896","#B7CA79","#B0CC99","#D2DBBF","#E0DCC3","#B7CA55","#9CC6E9","#E0B696","#B6C2D9"
,"#B7CA55","#9CC6E9","#FDF5E6","#B2DFEE","#E6E6FA","#66CDAA","#BC8F8F","#A2CD5A","#FFEC8B"
,"#EEB4B4","#DEB887","#F5DEB3","#D8BFD8","#CDAF95","#CD6090","#E0EEE0","#CDB7B5","#B4CDCD"])'
color2<-'d3.scaleOrdinal() .domain(["lncRNA","function"]) 
.range(["#057F76","#D44C3C"])'
lncfun<-read.table("I:/pan-cancer lncRNA课题/数据库扩展/all_lnc_function2.txt",sep = "\t",header = T)
lncfun$tumor[which(lncfun$tumor == "meningioma")]<-"MG"
lncfun<-lncfun%>%group_by(tumor)%>%top_n(800,-log10(p.adjust))
tumor<-unique(lncfun$tumor)
for(i in tumor){
  cancer<-lncfun[lncfun$tumor == i,]
  nodeData<-as.data.frame(rbind(cbind(name = unique(cancer$lncRNA),group = "lncRNA"),
                            cbind(name = unique(cancer$ID),group = "function")))
  linkData <- data.frame(source = (match(cancer$lncRNA, nodeData$name)-1),
                         target = (match(cancer$ID, nodeData$name)-1),
                         value = cancer$enrichmentScore,group2 = cancer$cell_type)
  degree<-table(c(cancer$lncRNA,cancer$ID))
  nodeData$size<-as.character(degree[nodeData$name])
  network<-forceNetwork(
    Links = linkData,
    Nodes = nodeData,
    Source = "source",
    Target = "target",
    Value = "value",
    Nodesize = 'size',
    NodeID = "name",Group = "group",
    opacity = 1, zoom = T, bounded = T,
    #clickAction = MyClickScript,
    legend = T,
    charge=-20,linkColour = "lightgrey",fontSize = 16,colourScale = color2,
    opacityNoHover = 0,)
  network$x$links$colour<-cellcol[linkData$group2]
  saveNetwork(network = network,paste("I:/pan-cancer lncRNA课题/数据库扩展/新加图/lnc_fun_html/",i,"_pathway.html",sep = ""))
}

MyClickScript <- 'alert("You clicked " + d.name + " is in " + d.group );'
colnc<-read.table("I:/pan-cancer lncRNA课题/数据库扩展/all_colncRNA.txt",sep = "\t",header = T)
tumor<-unique(colnc$tumor)
for(i in tumor){
  cancer<-colnc[which(colnc$tumor == i),]
  anno<-unique(as.data.frame(cbind(c(cancer$lncRNA1,cancer$lncRNA2),c(cancer$cell_type,cancer$cell_type))))
  anno<-anno[!duplicated(anno$V1),]
  rownames(anno)<-anno$V1
  cancer2<-aggregate(cancer$coindex,list(cancer$lncRNA1,cancer$lncRNA2,cancer$cell_type),sum)
  nodeData<-data.frame(name = unique(c(cancer$lncRNA1,cancer$lncRNA2)),group = "lncRNA")
  nodeData$cell<-anno[nodeData$name,2]
  linkData <- data.frame(source = (match(cancer2$Group.1, nodeData$name)-1),
                         target = (match(cancer2$Group.2, nodeData$name)-1),
                         value = cancer2$x,group2 = cancer2$Group.3)
  degree<-table(c(cancer2$Group.1,cancer2$Group.2))
  nodeData$size<-as.character(degree[nodeData$name])
  nodeData$name<-paste(nodeData$cell,nodeData$name,sep = " : ")
  network<-forceNetwork(
    Links = linkData,
    Nodes = nodeData,
    Source = "source",
    Target = "target",
    Value = "value",
    Nodesize = 'size',
    NodeID = "name",Group = "cell",
    opacity = 1, zoom = T, bounded = T,
    clickAction = MyClickScript,
    legend = F,
    charge=-20,linkColour = "lightgrey",fontSize = 16,colourScale = color,
    opacityNoHover = 0,)
  network$x$links$colour<-cellcol[linkData$group2]
  saveNetwork(network = network,paste("I:/pan-cancer lncRNA课题/数据库扩展/新加图/colnchtml/",i,"_colncnet.html",sep = ""))
  
}
network


lncfun<-read.table("I:/pan-cancer lncRNA课题/数据库扩展/all_lnc_function2.txt",sep = "\t",header = T)
lncfun$tumor[which(lncfun$tumor == "meningioma")]<-"MG"
lnc<-unique(lncfun$lncRNA)
for(i in lnc){
  lncdata<-lncfun[lncfun$lncRNA == i,]
  lncdata<-lncdata[,c(3,4,5,7,11,12)]
  lncdata$p<- -log(lncdata$p.adjust)
  lncdata$scale_NES<-scale(lncdata$enrichmentScore,center = T,scale = F)[,1]
  
  
  fig <- plot_ly(
    type = 'scatter',
    x = lncdata$p,
    y = lncdata$enrichmentScore,
    opacity = 0.5,
    text = paste("Function: ", lncdata$ID,
                 
                 "<br>Tumor: ", lncdata$tumor,
                 
                 "<br>Cell Type: ", lncdata$cell_type,
                 
                 "<br>EnrichmentScore: ", lncdata$enrichmentScore,
                 "<br>P.adj: ", lncdata$p.adjust),
    hoverinfo = 'text',
    mode = 'markers',
    transforms = list(
      list(
        type = 'groupby',
        groups = lncdata$func_type,
        styles = list(
          list(target = "GO", value = list(marker =list(color = '#A52A2A'))),
          list(target = "KEGG", value = list(marker =list(color = '#228B22'))),
          list(target = "IMMUNE", value = list(marker =list(color = '#1E90FF')))
        ))))%>%layout(title = i)
  saveNetwork(network = fig,paste("I:/pan-cancer lncRNA课题/数据库扩展/新加图/funchtml/",i,"_function_regulate.html",sep = ""))
  
}


for(i in lnc){
  lncdata<-lncfun[lncfun$lncRNA == i,]
  lncdata<-lncdata[,c(3,4,5,7,11,12)]
  lncdata$`-log10p`<- -log(lncdata$p.adjust)
  lncdata$scale_NES<-scale(lncdata$enrichmentScore,center = T,scale = F)[,1]
  fig <- plot_ly(lncdata, x = ~`-log10p`, y = ~enrichmentScore, text = ~ID, type = 'scatter', mode = 'markers',color = ~cell_type,colors = cellcol,
                 size = ~`-log10p`,sizes = c(5, 25),marker = list(sizemode = 'diameter', opacity = 0.5))
  
  fig <- fig %>% layout(title = i)
  saveNetwork(network = fig,paste("I:/pan-cancer lncRNA课题/数据库扩展/新加图/funchtml_buble/",i,"_function_regulate.html",sep = ""))
  
}

lnc<-unique(c(colnc$lncRNA1,colnc$lncRNA2))
for(i in lnc){
  lncdata<-colnc[which(colnc$lncRNA1 == i | colnc$lncRNA2 == i),]
  lncdata$lncRNA1<-paste(lncdata$cell_type,lncdata$lncRNA1,sep = " : ")
  lncdata$lncRNA2<-paste(lncdata$cell_type,lncdata$lncRNA2,sep = " : ")
  anno<-unique(as.data.frame(cbind(c(lncdata$lncRNA1,lncdata$lncRNA2),c(lncdata$cell_type,lncdata$cell_type))))
  anno<-anno[!duplicated(anno$V1),]
  rownames(anno)<-anno$V1
  lncdata2<-aggregate(lncdata$coindex,list(lncdata$lncRNA1,lncdata$lncRNA2,lncdata$cell_type),sum)
  nodeData<-data.frame(name = unique(c(lncdata$lncRNA1,lncdata$lncRNA2)),group = "lncRNA")
  nodeData$cell<-anno[nodeData$name,2]
  linkData <- data.frame(source = (match(lncdata2$Group.1, nodeData$name)-1),
                         target = (match(lncdata2$Group.2, nodeData$name)-1),
                         value = lncdata2$x,group2 = lncdata2$Group.3)
  degree<-table(c(lncdata2$Group.1,lncdata2$Group.2))
  nodeData$size<-as.character(degree[nodeData$name])
  network<-forceNetwork(
    Links = linkData,
    Nodes = nodeData,
    Source = "source",
    Target = "target",
    Value = "value",
    Nodesize = 'size',
    NodeID = "name",Group = "cell",
    opacity = 1, zoom = T, bounded = T,
    #clickAction = MyClickScript,
    legend = F,
    charge=-10,linkColour = "lightgrey",fontSize = 16,colourScale = color,
    opacityNoHover = 0,)
  network$x$links$colour<-cellcol[linkData$group2]
  saveNetwork(network = network,paste("I:/pan-cancer lncRNA课题/数据库扩展/新加图/lnc_colnc/",i,"_colnc.html",sep = ""))
}

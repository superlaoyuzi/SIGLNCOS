##tset
setwd("./test")
file<-"./testdata"

list.files(file)##Paths containing data from single cells
##[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"

cellannation<-read.csv("./Cell_annotation.csv",header = F)##cell annotation table

head(cellannation)
# V1           V2
# 1 AAACATACAACCAC-1   CD14+ Mono
# 2 AAACATTGAGCTAC-1            B
# 3 AAACATTGATCAGC-1   CD14+ Mono
# 4 AAACCGTGCTTCCG-1 Memory CD4 T
# 5 AAACCGTGTATGCG-1           NK
# 6 AAACGCACTGGTAC-1   CD14+ Mono

source("./SIGLNCOS.R")##load SIGLNCOS

SIGLNCOS(file = file,Tissue = "PBMC",Cell_annotation = cellannation)##operation SIGLNCOS


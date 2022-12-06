# https://www.jianshu.com/p/82787ddada38
# https://www.csdn.net/tags/OtTaYgxsMjAwNTUtYmxvZwO0O0OO0O0O.html
# https://blog.csdn.net/weixin_45822007/article/details/116870038

getwd()
rm(list = ls()) 

# 
up_prot<-read.csv("ratio-4.CSV")
down_prot<-read.csv("ratio-0.25.CSV")


# kegg
#install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
dupl<-duplicated(up_prot$id)
up_prot<-up_prot[!dupl,]
dupl<-duplicated(down_prot$id)
down_prot<-down_prot[!dupl,]

#write.table(res_kegg, file="res_kegg.xls",
#            quote=F, sep="\t", row.names=T, col.names=T)

up_prot.ENTREZID <- bitr(up_prot, 
                            fromType = c("SYMBOL"),
                            toType = c("ENTREZID"),
                            OrgDb = org.Mm.eg.db)
down_prot.ENTREZID <- bitr(down_prot, 
                              fromType = c("SYMBOL"),
                              toType = c("ENTREZID"),
                              OrgDb = org.Mm.eg.db)

up_kk <- enrichKEGG(gene = as.numeric(up_prot.ENTREZID$ENTREZID),
                     organism     = 'mouse',
                     #universe     = gene_all,
                     pvalueCutoff = 0.8,
                     qvalueCutoff =0.8)
down_kk <- enrichKEGG(gene = as.numeric(down_prot.ENTREZID$ENTREZID),
                    organism     = 'mouse',
                    #universe     = gene_all,
                    pvalueCutoff = 0.8,
                    qvalueCutoff =0.8)
barplot(up_kk,label_format=50,title="KEGG-Up")
barplot(down_kk,label_format=50,title="KEGG-Down")

#dotplot(down_kk,showCategory = 10, title="The KEGG enrichment analysis of all DEGs")+
#  scale_size(range=c(2, 12))+
#  scale_x_discrete(labels=function(down_kk) str_wrap(down_kk,width = 25))


# https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enrichGO
# https://blog.csdn.net/Eric_blog/article/details/119951736
# go Molecular_Function (MF), Biological_Process (BP), and Cellular Component (CC).
#barplot(go_up,showCategory=10,title="GO_Molecular_Function")


go_up<-enrichGO(gene = as.numeric(up_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "ALL",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_up,showCategory=20,font.size=12,label_format=70,title="GO: All")


go_up<-enrichGO(gene = as.numeric(up_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "MF",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_up,showCategory=20,font.size=12,label_format=70,title="GO: Molecular Function")


go_up<-enrichGO(gene = as.numeric(up_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "BP",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_up,showCategory=20,font.size=12,label_format=70,title="GO: Biological Process")

go_up<-enrichGO(gene = as.numeric(up_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "CC",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_up,showCategory=20,font.size=12,label_format=70,title="GO: Cellular Component")






# down
go_down<-enrichGO(gene = as.numeric(down_prot.ENTREZID$ENTREZID),
                OrgDb      = org.Mm.eg.db,
                #keyType    = 'ENSEMBL',
                ont        = "ALL",
                #pAdjustMethod = "BH",
                pvalueCutoff = 0.8,
                qvalueCutoff = 0.8)
barplot(go_down,showCategory=10,label_format=70,title="GO: All")


go_down<-enrichGO(gene = as.numeric(down_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "MF",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_down,showCategory=10,label_format=70,title="GO: Molecular Function")


go_down<-enrichGO(gene = as.numeric(down_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "BP",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_down,showCategory=10,label_format=70,title="GO: Biological Process")

go_down<-enrichGO(gene = as.numeric(down_prot.ENTREZID$ENTREZID),
                  OrgDb      = org.Mm.eg.db,
                  #keyType    = 'ENSEMBL',
                  ont        = "CC",
                  #pAdjustMethod = "BH",
                  pvalueCutoff = 0.8,
                  qvalueCutoff = 0.8)
barplot(go_down,showCategory=10,label_format=70,title="GO: Cellular Component")


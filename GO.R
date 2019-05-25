#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")

library(clusterProfiler)
library(org.Mm.eg.db)

RU<- as.data.frame(fc_significant)
names(RU)[names(RU) == colnames(RU)[1]]<- "UNIPROT"
columns(org.Mm.eg.db)

gene.RU <- bitr(RU[,1], fromType = "UNIPROT", 
                toType = c("ENSEMBL","ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
library(dplyr)
gene<- gene.RU %>% distinct(gene.RU$ENTREZID)

names(gene)[names(gene) ==colnames(gene)[1]]<- "ENTREZID"
library(dplyr)
combineRU<- left_join(RU, gene.RU, by = "UNIPROT")
combineRU.omit<- na.omit(combineRU)

combiRU<- combineRU.omit[,c(7,2)]
ctRU<- t(gene)
colnames(ctRU) <- as.character(unlist(ctRU[1,]))
ctRU = ctRU[-1, ]
combineRU.3<- list(ctRU)
names(ctRU)
GO.RU<- names(ctRU)
###GO Enrichment 
#for enrichment based on gene symbol
#ego4 <- enrichGO(gene.df1$SYMBOL,'org.Mm.eg.db', keytype = "SYMBOL",
#ont = "BP",pvalueCutoff=0.01)
ego <- enrichGO(GO.RU,'org.Mm.eg.db',
                 ont = "CC",pvalueCutoff=0.01)

head(summary(ego5))
ego5$BgRatio

#Visualizations

dotplot(ego5,x = "GeneRatio", color = "p.adjust", 
        font.size = 7, showCategory=NULL)

enrichMap(ego5, n=50, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

cnetplot(ego5, showcategory = 5, categorySize = "pvalue", font.size = 10)

#Enrichment geneset analysis
genr <- gseGO(geneList=GO, ont="CC", OrgDb=org.Mm.eg.db, verbose=F)
head(summary(genr))

library(clusterProfiler)
library(org.Mm.eg.db)

fil<- data.frame(KEGG_enrichment)
names(fil)[names(fil) == colnames(fil)[1]]<- "UNIPROT"
names(fil)[names(fil) == colnames(fil)[2]]<- "logFC"
names(fil)[names(fil) == colnames(fil)[3]]<- "Groups"

######compare cluster########

fil1<- data.frame(KEGG_enrichment)
names(fil1)[names(fil1) == colnames(fil1)[1]]<- "UNIPROT"
df<- list(fil1[,c(1,2)])
df
df2<- as.factor(df$UNIPROT)

gene.df1<- bitr(fil1[,1], fromType = "UNIPROT",
                toType = c("ENSEMBL","ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
library(dplyr)
combine<- left_join(fil1, gene.df1, by = "UNIPROT")
combine.omit<- na.omit(combine)

mydf<- data.frame(combine.omit[, c(1:3, 5)])

#mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$logFC) > 0,]
mydf$Regulation <- "upregulated"
mydf$Regulation[mydf$logFC < 0] <- "downregulated"
mydf
mydf1<- mydf[,3:5]
formula_res<- compareCluster(ENTREZID~Groups+Regulation,organism='mmu', 
                     data=mydf1, fun="enrichKEGG")

func<- as.data.frame(formula_res)
func1<- func[,c(4,5,12)]
write.csv(func, file = 'KEGGEnriched.csv')

#source("https://bioconductor.org/biocLite.R")
#biocLite("IRanges")


class(formula_res)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

dotplot(formula_res, font.size = 7, by='count',includeAll = T,showCategory = 100,
  x=~Regulation)+facet_grid(~Groups)+geom_point(size=1)

# input MS quantification files for Downstream Processing

imp<- as.data.frame(example)
medium<- as.data.frame(imp)
colnames(medium[,3:11]) <- c("U1","U2", "U3","R1","R2", "R3","H1","H2", "H3")
##OR
names(medium)[names(medium) == colnames(medium)[3]] <- "U1"
#......
m1<- medium[,1:11]

# Log transform
m1.log<-log(m1[,3:11],2) 
m2.log<- cbind(m1[,1:2], m1.log)
is.na(m1.log)

final.2 <- as.data.frame(m2.log[rowSums(is.na(m2.log)) <= 6, ]) #min 6 NA's
colSums(is.na(final.2))
imputed<- cbind(m4[,1:2], m2.complete)

# Missmap
library(Amelia)
rows1<- m2.log[,1]
missmap(m2.log[,2:10], legend=T, col = c("wheat", "darkred"), x.cex = 0.5, 
        y.cex = 0.3, y.lables = rows1 )

#Imputation
library(mice)# for more robustness
md.pattern(m2.log)
md.pairs(m2.log)# freq of missing vals b/w pairs 
#m=1 number of multiple imputations
#maxit=2 number of iterations. 10-20 is sufficient.
CleanData.imp <- mice(m2.log, m=5, maxit=2, printFlag=TRUE) 
summary(CleanData.imp)

# Normalization

quanData<- m2.complete[,c(1,2,4,5,7,8)]
b <- data.frame(name=c("U1","U2",
                       "R1" ,"R2",
                       "H1","H2"))

sampleData<- data.frame(Sample.names=c("U1","U2",
                                       "R1" ,"R2",
                                       "H1","H2"),
                        Condition = c("U","U","R","R","H","H"), 
                        Bio.rep = c(1,2,3,4,5,6))

quanData<-as.matrix(quanData)

#library(tibble)
#col<-column_to_rownames(sampleData, var = "Sample.names")
#labels=c("U","U","R" ,"R","H","H")
normalized<- normalizeD(quanData, labels, "Quantile Centering", "within conditions", 
                            quantile = 0.5)
normalized <- normalizeD(quanData, labels, "Median Centering", method = "overall")
library(preprocessCore)
norm1<- data.frame(normalize.quantiles(quanData,copy=TRUE))
norm2<- norm1[,c(1,2,4,5,7,8)]

colnames(norm2) <- c("U1","U2","R1","R2","H1","H2")
normalized<- cbind(m4[,1:2], norm2)
write.csv(imputed, file = "imputed.csv")

#HEATMAP
dat.clust <- dist(normalized[,3:8])
dat2.clust <- hclust(dat.clust, method = "complete")
dat.dend <- as.dendrogram(dat2.clust)

library(dendextend)
# Color the branches based on the clusters:
dat.dend <- color_branches(dat.dend) 
# hang the dendrogram:
dat.dend <- hang.dendrogram(dat.dend,hang_height=0.1)
# reduce size of the labels:
dat.dend <- set(dat.dend, "labels_cex", 0.5)
some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), 
                                                      l = c(30, 90), power = c(1/5, 1.5)))

#using package d3heatmap
d3heatmap::d3heatmap(n3,
                     dendrogram = "both",#labRow = rownames(normalized1[,1]),
                     #Rowv = dat.dend,
                     colors = colorRamp(c("blue","white")),
                     # scale = "row",
                     width = 1000, height = 600,
                     show_grid = F, theme ="light")

###
library(gplots)
clusters = dendextend::cutree(dat2.clust, k = 5)
heatmap.2(n3, 
          distfun = dist, Rowv = T, Colv=T,
          hclustfun = hclust, 
          sepwidth = c(50,50),
          labRow = F,
          sepcolor="white", 
          #labCol = colnames(normalized1[,3:11]),
          dendrogram = "both", trace= 'none',rowsep=clusters, 
          col=colorRampPalette(c("blue","white","red"))(50))

##Correlation plot
library(reshape2)
n2<-melt(normalized)
library(GGally)
ggscatmat(normalized, corMethod = "pearson")

### for condition R/U
condition1<- 'U'
condition2<- 'R'
Diff<- diffAnaLimma(normalized[,3:8], sampleData,labels, condition1, condition2)
Diff.RU<- data.frame(m4[,7], normalized1[,2:7],Diff)
names(Diff.RU)[names(Diff.RU) == colnames(Diff.RU)[1]]<- "Accession"
Diff.RU <- cbind(Diff.RU, adj.pval$adjusted.p) 
names(Diff.RU)[names(Diff.RU) == colnames(Diff.RU)[10]]<- "Adj.pval"
write.csv(Diff.RU, file = "Differential-RU.csv")

### for condition H/U
condition1<- 'U'
condition2<- 'H'
Diff1<- diffAnaLimma(normalized1[,2:7], sampleData,labels, condition1, condition2)
Diff.HU<- data.frame(m4[,7], normalized1[,2:7],Diff1)
names(Diff.HU)[names(Diff.HU) == colnames(Diff.HU)[1]]<- "Accession"
Diff.HU <- cbind(Diff.HU, adj.pval$adjusted.p) 
names(Diff.HU)[names(Diff.HU) == colnames(Diff.HU)[10]]<- "Adj.pval"
write.csv(Diff.HU, file = "Differential-HU.csv")

#### for condition R/H
condition1<- 'H'
condition2<- 'R'
Diff2<- diffAnaLimma(normalized1[,2:7], sampleData,labels, condition1, condition2)
Diff.RH<- data.frame(m4[,7], normalized1[,2:7],Diff2)
names(Diff.RH)[names(Diff.RH) == colnames(Diff.RH)[1]]<- "Accession"
Diff.RH <- cbind(Diff.RH, adj.pval$adjusted.p) 
names(Diff.RH)[names(Diff.RH) == colnames(Diff.RH)[10]]<- "Adj.pval"
write.csv(Diff.RH, file = "Differential-RH.csv")

##caliberation plot
library(cp4p)
adjusted<- adjust.p(pVal, pi0.method = 1, alpha = 0.05, nbins=10, pz = 0.05)
adj.pval<- adjusted[["adjp"]]
calibration.plot(pVal, pi0.method = "ALL", nbins = 20, pz = 0.05)

#volcano plot
logFC<- Diff.HU[,9]
pVal<- Diff.HU[,8]
pVal
diffAnaVolcanoplot(logFC, pVal, threshold_pVal = 2,
                   threshold_logFC = 2, conditions = 1)

#ggpot
library(ggplot2)
library(gridExtra)

#suppressPackageStartupMessages(library("plotly"))

diff_df<- data.frame(Diff.HU[,1], logFC,pVal, adj.pval[,2])
names(diff_df)[names(diff_df) == colnames(diff_df)[1]]<- "Accession"
names(diff_df)[names(diff_df) == colnames(diff_df)[4]]<- "adj.pval"
#a<- diff_df[,1]
diff_df$Accession<- as.character(diff_df$Accession)

diff_df["Category"] <- "Non-significant"

# highlight-
# FDR < 0.05 (significance level) & Fold Change > 1.5

# change the grouping for the entries with significance but not a large enough Fold change
diff_df[which(diff_df['pVal'] < 0.01 & 
                abs(diff_df['logFC']) < 2 ),"Category"] <- "ttest-significant"

# change the grouping for the entries a large enough Fold change but not a low enough p value
diff_df[which(diff_df['pVal'] > 0.01 & 
                abs(diff_df['logFC']) > 2 ),"Category"] <- "FC-significant"

# change the grouping for the entries with both significance and large enough fold change
diff_df[which(diff_df['pVal'] < 0.01 & 
                abs(diff_df['logFC']) > 2 ),"Category"] <- "All-significant"

# Find and label the top peaks..
top_peaks <- diff_df[with(diff_df, order(logFC, pVal)),][1:10,]
top_peaks <- rbind(top_peaks, diff_df[with(diff_df, 
                                           order(-logFC, pVal)),][1:10,])
top_peaks["Category"]<- "top-hits"

# Add gene labels for all of the top genes we found
# Create an empty list, and fill it with entries for each row in the dataframe
# each list entry is another list with named items that will be used by Plot.ly
a <- list()
for (i in seq_len(nrow(top_peaks))) {
  m <- top_peaks[i, ]
  a[[i]] <- list(
    a = m[["logFC"]],
    b = -log10(m[["pVal"]]),
    text = m[["Accession"]],
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 0.5,
    ax = 20,
    ay = -40
  )
}

#with ggplot
diff_df1<- diff_df[,c(1:3,5)]
write.csv(diff_df, file = "diff_dfHR.csv")

library(ggalt)
gg<- ggplot(data = diff_df, aes(x = logFC, y = -log10(pVal)))+
  geom_point(aes(colour= Category))+labs(title= "Volcano plot-HU") 
gg 

#######################################################

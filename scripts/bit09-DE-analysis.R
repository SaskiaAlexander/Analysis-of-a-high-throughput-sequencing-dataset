################################################################################
### LOAD REQUIRED PACKAGES/LIBRARIES
################################################################################
## Bioconductor previous way of installing
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("edgeR")
# biocLite("Glimma")
## When you get Error: With R version >=3.5 --> https://www.bioconductor.org/install/
## use installation instructions below
install.packages("BiocManager")
library(BiocManager)
## Install specific bioconductor packages:
BiocManager::install(c("limma","edgeR"))
# 2022-01 Glimma dependencies: 
# sudo dnf install libcurl-devel openssl-devel libxml2-devel
BiocManager::install("Glimma")

# Loading required packages
library("limma")
library("edgeR")
library("Glimma")
library("RColorBrewer") # to use brewer.pal
library("Rtsne")        # t-SNE

# Bioconductor packages for gene annotations
# This is OPTIONAL! --> lots of dependencies and probably updating required
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# biocLite("Homo.sapiens")
# biocLite("Mus.musculus")
# biomaRt primarily works off Ensembl gene IDs, whereas Mus.musculus and 
# Homo.sapiens packages information from various sources and allows users 
# to choose between many different gene IDs as the key.



################################################################################
### The typical workflow for analyzing RNA-seq data described below
### is mainly based on the paper by Law et al, F1000 Research 2016
################################################################################
### DATA PACKAGING: READING IN COUNT DATA & CREATING DGEList OBJECT
# Text files can be read into R separately and combined into matrix of counts
# --> using readDGE function
# If counts from all samples are stored in single file 
# --> read data into R and convert into DGEList-object with DGEList function
################################################################################
# To start we need counts summarised at gene-level:
# read table from all_counts.txt
# setwd("/media/sf_VMshare/BIT09-HTA/BIT09-DE")
countData <- read.table(file = "all_counts.txt",
                        sep = "\t", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
# Check first 5 rows
countData[1:5,1:3]

# If first column contains genes (not counts) --> use as rownames
rownames(countData) <- countData$X
# Check first 5 rows again
countData[1:5,]
# Remove filename column
countData <- countData[,-1]
# Check first 5 rows again
countData[1:5,]

# Create DGEList object (requires edgeR) which will contains count data
DGEobject <- DGEList(counts = countData)
# Check dimensions 
dim(DGEobject)
# [1] 24421  6 --> 24421 rows=genes and 6 columns=samples

# Alternatively if you would read all .txt sample files
# instead of using the bit09-merge-htseqcounts.R script
#pathFiles <- "/media/sf_VMshare/BIT09-HTA/BIT09-DE/htseqcount/"
#files <- dir(path = pathFiles, pattern = "*sorted.txt$")
#x <- readDGE(files, path = pathFiles)
# You will see a warning: Meta tags detected: __no_feature, __ambiguous, ...
#class(x)
# [1] "DGEList"
# attr(,"package")
# [1] "edgeR"
#dim(x)
# [1] 24425  6



################################################################################
### DATA PACKAGING: READING & ORGANISING SAMPLE INFORMATION (METADATA)
################################################################################
# We also need metadata describing the experiment and samples:
# e.g. read table SraRunTable.csv or your own metaData file
metaData <- read.table(file = "SraRunTable.csv",
                       sep = ",",
                       #quote = '"',
                       header = TRUE, 
                       stringsAsFactors = TRUE)
# Show colnames
colnames(metaData)
# Show most important columns e.g. "Run" and "genotype" columns
metaData[,c("Run","genotype")]
# metaData <- metaData[order(metaData$Run),]

# Add SraRunTable info to DGEobject$samples
DGEobject$samples <- cbind(DGEobject$samples, metaData[,c("Run","genotype")])
DGEobject$samples
#extracol <- c(rep("WT",3), rep("KO",3))
#DGEobject$samples <- cbind(DGEobject$samples, extracol)
extracol <- c(rep("WT",3), rep("KO",3))
extracol2 <- c("S1","S2","S3","S4","S5","S6")
DGEobject$samples <- cbind(DGEobject$samples, extracol)
DGEobject$samples <- cbind(DGEobject$samples, extracol2)

# Additional you can add a genes data frame in the DGEList object
# to store gene annotation information e.g. gene symbols, gene names,
# chromosome names and locations, Entrez gene IDs, Refseq gene IDs and 
# Ensembl gene IDs to name just a few.



################################################################################
### CHECKING LIBRARY SIZES AND SOME GENES (COUNTS)
################################################################################
# Check the library sizes of the samples
DGEobject$samples$lib.size

# BARPLOT LIB SIZES (all samples)
# savePath <- "/home/user/results"
# pdf(paste(savePath,"QC-barplot-libsizes.pdf",sep=""))
par(mar = c(8,6,4,4)) # more margin: bottom, left, top, right
bp <- barplot(DGEobject$samples$lib.size*1e-6,
              axisnames = FALSE,
              main = "Barplot library sizes",
              ylab = "Library size (millions)")
axis(1, labels = rownames(DGEobject$samples), 
     at = bp, las = 2, cex.axis = 0.8)
# dev.off()

# To check counts of genes
DGEobject$counts[1:5,]
# --> shows counts for first five rows (genes)
DGEobject$counts["Actb",]
# --> shows Brca1 counts for each sample 
DGEobject$counts["Tnfrsf1b",]
# --> shows Tnfrsf1b counts 

# To look for genes by gene symbol with grep
rownames(DGEobject$counts) 
# --> all the gene symbols
grep("Brca", rownames(DGEobject$counts), 
     perl = TRUE, value = TRUE)
# --> returns all gene symbols containing Brca --> Brca1, Brca2
# Now use this to see the counts for those genes
DGEobject$counts[grep("Brca", rownames(DGEobject$counts), 
                      perl = TRUE, value = TRUE),]
# For Tnf* genes
DGEobject$counts[grep("Tnf", rownames(DGEobject$counts), 
                      perl = TRUE, value = TRUE),]



################################################################################
### DATA PRE-PROCESSING: Transformations from the raw scale
################################################################################
# Common practice to transform raw counts onto a scale that accounts 
# for such library size differences. Popular transformations include 
# counts per million (CPM), log2-counts per million (log-CPM), 
# reads per kilobase of transcript per million (RPKM), and fragments per 
# kilobase of transcript per million (FPKM)
cpm <- cpm(DGEobject)
log.cpm <- cpm(DGEobject, log = TRUE)
# Genes that are not expressed
no.samples <- nrow(DGEobject$samples)
table(rowSums(DGEobject$counts==0)==no.samples)
# FALSE  TRUE 
# 19558  4863
# --> 4863 genes have no expression (counts=0) in any of the samples
DGEobject$counts[rowSums(DGEobject$counts==0)==no.samples,]
DGEobject$counts["Zswim2",]



################################################################################
### DATA PRE-PROCESSING: Reduce subset of genes
################################################################################
# Genes that are not expressed at a biologically meaningful level 
# in any condition should be discarded to reduce the subset of genes 
# to those that are of interest, and to reduce the number of tests carried out 
# downstream when looking at differential expression.
# Using a nominal CPM value of 1 (which is equivalent to a log-CPM value of 0 
# a gene is deemed to be expressed in a given sample if its transformed count 
# is above this threshold, and unexpressed otherwise. 
# Genes must be expressed in at least one group (or roughly any three samples 
# across the entire experiment) to be kept for downstream analysis.
keep.exprs <- rowSums(cpm>1)>=3
DGEfiltered <- DGEobject[keep.exprs,, keep.lib.sizes = FALSE]
dim(DGEfiltered)
# [1] 14471  6

# Now transform subset to get filtered+transformed cpmF (and log.cpmF)
cpmF <- cpm(DGEfiltered)
log.cpmF <- cpm(DGEfiltered, log = TRUE)
# Using this criterion, the number of genes is reduced
# DGEList-object "DGEobject" contains raw pre-filtered data
# DGEList-object "DGEfiltered" contains post-filtered data

# Density plots of log-CPM values for raw pre-filtered data
# and post-filtered data (yF1000R)
nsamples = ncol(DGEfiltered)
# brewer.pal Paired max 12 colors # library(RColorBrewer)
sample.col <- brewer.pal(nsamples, "Paired")
# alternative use own list of colors to randomly choose from

#pdf("/pathToResults/logCPM-densityplots.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
# Density plot part A: raw data pre-filtered
plot(density(log.cpm[,1]), 
     col = sample.col[1], lwd = 2, las = 1, ylim = c(0,0.25),
     main = "A. Raw data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
    den <- density(log.cpm[,i])
    lines(den$x, den$y, col = sample.col[i], lwd = 2)
}
legend("topright", rownames(DGEobject$samples), 
       text.col = sample.col, cex = 0.7, bty = "n")
# Density plot part B: post-filtered
plot(density(log.cpmF[,1]), 
     col = sample.col[1], lwd = 2, las = 2, ylim = c(0,0.25),
     main = "B. Filtered data", xlab = "Log-cpm")
abline(v = 0, lty = 3)
for (i in 2:nsamples){
  den <- density(log.cpmF[,i])
  lines(den$x, den$y, col = sample.col[i], lwd = 2)
}
legend("topright", rownames(DGEfiltered$samples), 
       text.col = sample.col, cex = 0.7, bty="n")
#dev.off()



################################################################################
### DATA PRE-PROCESSING: Normalising gene expression distributions
################################################################################
## Normalisation generally required to ensure that expression distributions
## of each sample are similar across the entire experiment.
## Any plot showing the per sample expression distributions, such as a density or boxplot,
## is useful in determining whether any samples are dissimilar to others. 
## Normalisation by the method of trimmed mean of M-values (TMM)
## is performed using the calcNormFactors function in edgeR.
# cpmF = filtered subset # log.cpmF = log cpm filtered subset
# Copy unnormalized, filtered data to show in boxplot
DGEf.notnorm <- DGEfiltered 
DGEf.norm <- calcNormFactors(DGEfiltered, method = "TMM")
DGEf.norm$samples$norm.factors
# [1] 0.9979337 0.9948348 1.0076473 0.9993417 1.0057461 0.9945725
# --> scaling factors relatively close to 1 --> mild TMM normalisation

## Boxplots of log-CPM values showing expression distributions 
## for unnormalised data (A) and normalised data (B) for each sample
par(mfrow = c(1,2))
par(mar = c(8,2,4,2))
# Boxplot unnormalised data
DGEf.notnorm$samples$norm.factors
# [1] 1 1 1 1 1 1
log.cpmF.notnorm <- cpm(DGEf.notnorm, log = TRUE)
boxplot(log.cpmF.notnorm, las = 2, col = sample.col, cex = 0.9, 
        main = "A. Unnormalised data", ylab = "Log-cpm")
# Boxplot normalised data
DGEf.norm$samples$norm.factors
# [1] 0.9979337 0.9948348 1.0076473 0.9993417 1.0057461 0.9945725
log.cpmF.norm <- cpm(DGEf.norm, log = TRUE)
boxplot(log.cpmF.norm, las = 2, col = sample.col, cex = 0.9, 
        main = "B. Normalised data", ylab = "Log-cpm")
# --> only when scaling factors are not close to 1 you will see a big difference



################################################################################
### DATA PRE-PROCESSING: Unsupervised clustering of samples
################################################################################
## Most important exploratory plots to examine for gene expression analyses 
## is the multidimensional scaling (MDS) plot.
## The plot shows similarities and dissimilarities between samples 
## in an unsupervised manner so that one can have an idea of the extent to which
## differential expression can be detected before carrying out formal tests.
## Ideally, samples would cluster well within the primary condition of interest.
## MDS plot of log-CPM values using limma plotMDS function
## Dimensions 1 and 2 are examined using color grouping defined by sample group
par(mfrow = c(1,1), mar = c(8,6,4,2))
col.group <- DGEfiltered$samples$genotype
# OR: col.group <- metaData$genotype
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
# Warning if you only have two groups
col.group
# [1] #377EB8 #377EB8 #377EB8 #E41A1C #E41A1C #E41A1C
# Levels: #E41A1C #377EB8 #4DAF4A
col.group <- as.character(col.group)
col.group
# [1] "#377EB8" "#377EB8" "#377EB8" "#E41A1C" "#E41A1C" "#E41A1C"
plotMDS(log.cpmF.norm, 
        # pch = 21,
        # labels = DGEfiltered$samples$genotype
        col = col.group,
        cex = 0.7, 
        xlim = c(-4.0,4.0), 
        ylim = c(-2.0,2.0),
        dim = c(1,2),
        main="MDS plot")
# OPTIONAL: interactive MDS plot (MDS-Plot.html) using glMDSPlot from Glimma
glMDSPlot(log.cpmF.norm, 
          labels = DGEfiltered$samples$genotype,
          # select groups (columns) from metaData
          groups = metaData[,c("genotype","Run")], 
          launch = TRUE)



################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Design matrix and contrasts
################################################################################
# Comparing Wildtype (WT) vs TNFR2-knockout (KO)
colnames(DGEf.norm)
DGEfiltered$samples$genotype
group <- as.factor(c("WT", "WT", "WT", "KO", "KO", "KO"))
# OR: group <- DGEfiltered$samples$genotype

## DESIGN MATRIX
design <- model.matrix(~0+group)
colnames(design)
# [1] "groupKO" "groupWT"
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)
# [1] "KO" "WT"
design

# sample_group <- c(rep("WT",3), rep("T1",3), rep("T2",3))
# sample_group
# design <- model.matrix(~0+sample_group)
# colnames(design) <- gsub("sample_group", "", colnames(design))
# design

## CONTRASTS for pairwise comparisons between groups 
## are set up in limma using the makeContrasts function.
contr.matrix <- makeContrasts(
  KOvsWT = KO-WT,
  levels = colnames(design))
contr.matrix

## If you have more conditions e.g. three types of samples 
## the model will be more complex ()
contrast3 <- matrix(c(-1,1,0,-1,0,1), ncol = 2)
dimnames(contrast3) <- list(c("cont","trt1","trt2"),
                            c("trt1-cont","trt2-cont"))
contrast3
# trt1-cont trt2-cont
# cont        -1        -1
# trt1         1         0
# trt2         0         1



################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Removing heteroscedasticity
################################################################################
## Removing heteroscedasticity from count data
## In limma, linear modelling is carried out 
## on the log-CPM values ... by the voom function.
par(mfrow = c(1,2), mar = c(8,6,4,2))
vDGEf <- voom(DGEf.norm, design, plot = TRUE)
vDGEf
## Fitting linear models for comparisons of interest
## Linear modelling in limma is carried out using 
## the lmFit and contrasts.fit functions 
vfitDGEf <- lmFit(vDGEf, design)
vfitDGEf <- contrasts.fit(vfitDGEf, 
                          contrasts = contr.matrix)
efitDGEf <- eBayes(vfitDGEf)
plotSA(efitDGEf, main = "Final model")



################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: Examining the DE genes
################################################################################
### Examining the number of DE genes
## Significance is defined using an adjusted p-value cutoff that is set at 5% by default.
summary(decideTests(efitDGEf))

## Some studies require more than an adjusted p-value cut-off. 
## When testing requires genes to have a log-FC that is significantly greater than 1 
## (equivalent to a 2-fold difference between cell types on the original scale).
tfitWTvsKO <- treat(vfitDGEf, lfc = 1)
dtWTvsKO <- decideTests(tfitWTvsKO)
summary(dtWTvsKO)


## topTable for results using eBayes
WT.vs.KO <- topTable(efitDGEf, 
                     coef = 1, 
                     n = Inf)
#T2.vs.WT <- topTable(efitDGEf, 
#                     coef = 2, 
#                     n = Inf)
head(WT.vs.KO, n = 10)
topgenes <- head(WT.vs.KO, n = 10)
topgenes
# Top genes up and down sorted by logFC
topgenes.ordered <- topgenes[order(topgenes$logFC),] 
head(topgenes.ordered, n = 5)
tail(topgenes.ordered, n = 5)
top10genes.ordered <- head(topgenes.ordered, n = 5)
top10genes.ordered <- rbind(tail(topgenes.ordered, n = 5),
                            top10genes.ordered)
top10genes.ordered <- top10genes.ordered[order(top10genes.ordered$logFC),] 
top10genes.ordered

## topTreat for results using treat 
# tfitWTvsKO <- treat(vfitDGEf, lfc = 1)
topTreat.WTvsKO <- topTreat(tfitWTvsKO, coef = 1, n = Inf)
head(topTreat.WTvsKO, n = 10)
tail(topTreat.WTvsKO, n = 10)
topTreatgenes <- head(topTreat.WTvsKO, n = 4)
topTreatgenes
# Top genes
topTreatgenes.ordered <- topTreatgenes[order(topTreatgenes$logFC),] 
topTreatgenes.ordered
# write.csv(top10genes.ordered, file = "top10genes-ordered.csv")

################################################################################
### DIFFERENTIAL EXPRESSION ANALYSIS: graphical representations of DE genes
################################################################################
par(mfrow = c(1,1), mar = c(8,6,4,2))
# Volcano plot 
volcanoplot(efitDGEf, coef = 1, highlight = 4, 
            names = rownames(efitDGEf$coefficients))
## MD plot
plotMD(tfitWTvsKO, column = 1, status = dtWTvsKO[,1], 
       main = "KOvsWT", xlim = c(-2,15))

## Create your own color palette to use in heatmap
color.palette <- colorRampPalette(c("red", "yellow", "green"))(n = 200)
color.palette.2  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(n = 200)
color.palette.3 <- colorRampPalette(c("#0d50b2", 
                                      "white", 
                                      "#c5081a"))(n = 200)

## Heatmap topTable genes
# Get names of genes to show
genes10topTable <- rownames(top10genes.ordered) # top
# Packages required to use heatmap and str_count
# install.packages("gplots")
library("gplots")
# install.packages("stringr")
library("stringr")
# Create heatmap
heatmap.2(vDGEf$E[genes10topTable,],
          ## DENDOGRAM CONTROL
          Rowv = FALSE, # do not reorder rows!
          # distfun = dist,
          # hclustfun = hclust,
          dendrogram = "column",
          # symm = FALSE,
          ## DATA SCALING
          scale = "row",
          ## image plot
          # revC = identical(Colv, "Rowv"),
          # breaks,
          col = color.palette.3, 
          ## BLOCK SEPARATION
          # colsep,
          # rowsep,
          # sepcolor="white",
          # sepwidth=c(0.05,0.05),
          ## CELL LABELING
          # cellnote,
          # notecex=1.0,
          # notecol="cyan",
          # na.color=par("bg"),
          ## LEVEL TRACE
          trace = "none", 
          ## ROW/COLUMN LABELING
          margins = c(12,8),
          # ColSideColors,
          # RowSideColors,
          cexRow = 0.9,
          cexCol = 0.9,
          labRow = genes10topTable, 
          labCol = str_c(colnames(DGEf.norm), group, sep = "-"),
          ## COLOR KEY + DENSITY INFO
          key = TRUE,
          density.info = "none",
          ## OTHER
          lmat = rbind(4:3,2:1), 
          lhei = c(1,4), 
          lwid = c(1,4)
)
# The position of each element in the heatmap.2 plot can be controlled using 
# the lmat, lhei and lwid parameters. These are passed by heatmap.2 to the layout command as:
# layout(mat = lmat, widths = lwid, heights = lhei)
# lmat is a 2x2 matrix 
# 1 Heatmap - 2 Row dendrogram - 3 Column dendrogram - 4 Key
# --> default: rbind(4:3,2:1)
#      [,1] [,2]
# [1,]    4    3
# [2,]    2    1
# --> key underneath the heatmap
# lmat = rbind(c(0,3),c(2,1),c(0,4))lmat = rbind(c(0,3),c(2,1),c(0,4))
# [,1] [,2]
# [1,]    0    3
# [2,]    2    1
# [3,]    0    4

## Heatmap topTreat genes
genestopTreat <- rownames(topTreatgenes.ordered)
heatmap.2(vDGEf$E[genestopTreat,], 
          scale = "row",
          labRow = genestopTreat, 
          labCol = str_c(colnames(DGEf.norm), 
                         group, sep = "-"),
          col = color.palette.3, 
          trace = "none", density.info = "none",
          Rowv = FALSE, # do not reorder rows!
          margin = c(15,8), lmat = rbind(4:3,2:1), 
          lhei = c(2,6), lwid = c(1.5,4), 
          dendrogram = "column")



################################################################################
### OPTIONAL: check some genes
################################################################################
## Check raw counts for some genes of interest
DGEobject$counts["Tnfrsf1b",]
DGEobject$counts["Fbxo44",]
DGEobject$counts["Slc37a2",]
DGEobject$counts["Gbp2b",]

par(mfrow = c(1,1), mar = c(6,9,4,2))
expTable <- vDGEf$E
# Tnfrsf1b
geneSymbol <- "Tnfrsf1b"
DGEobject$counts[geneSymbol,]
mp <- barplot(expTable[geneSymbol,], 
              xlim = c(0,6), horiz = TRUE, 
              col = "skyblue", 
              axes = F, axisnames = F, 
              main = geneSymbol)
axis(1)
axis(2, labels = colnames(expTable), 
     at = mp, las = 1, cex.axis = 0.8)
title(main = geneSymbol)
# Fbxo44
geneSymbol <- "Fbxo44"
DGEobject$counts[geneSymbol,]
mp <- barplot(expTable[geneSymbol,], xlim = c(0,4), horiz = TRUE, 
              col = "skyblue", axes = F, axisnames = F, main = geneSymbol)
axis(1)
axis(2, labels = colnames(expTable), at = mp, las=1, cex.axis=0.8)
title(main = geneSymbol)
# Slc37a2
geneSymbol <- "Slc37a2"
DGEobject$counts[geneSymbol,]
mp <- barplot(expTable[geneSymbol,], xlim = c(0,10), horiz = TRUE, 
              col = "skyblue", axes = F, axisnames = F, main = geneSymbol)
axis(1)
axis(2, labels = colnames(expTable), at = mp, las=1, cex.axis=0.8)
title(main = geneSymbol)
# Gbp2b
geneSymbol <- "Gbp2b"
DGEobject$counts[geneSymbol,]
mp <- barplot(expTable[geneSymbol,], xlim = c(-2,6), horiz = TRUE, 
              col = "skyblue", axes = F, axisnames = F, main = geneSymbol)
axis(1)
axis(2, labels = colnames(expTable), at = mp, las=1, cex.axis=0.8)
title(main = geneSymbol)
################################################################################
# Barplot raw counts versus log-CPM counts
################################################################################
par(mfrow = c(1,2), mar = c(6,9,4,2))
# DGEobject$counts[geneSymbol,]
geneSymbol <- "Tnfrsf1b"
bp1 <- barplot(DGEobject$counts[geneSymbol,], 
               xlim = c(0,2000), horiz = TRUE, 
               col = "skyblue", 
               axes = F, axisnames = F)
axis(1)
axis(2, labels = colnames(DGEobject$counts), 
     at = bp1, las = 1, cex.axis = 0.8)
title(main = paste0("Raw counts ",geneSymbol))
# expTable[geneSymbol,]
bp2 <- barplot(vDGEf$E[geneSymbol,], 
               xlim = c(0,6), horiz = TRUE, 
               col = "skyblue", 
               axes = F, axisnames = F)
axis(1)
axis(2, labels = colnames(vDGEf$E), 
     at = bp2, las = 1, cex.axis = 0.8)
title(main = paste0("log-CPM counts ",geneSymbol))



################################################################################
### OPTIONAL: Kallisto quantification
################################################################################
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
# biocLite("tximportData")
# biocLite("rhdf5")
library("tximport")
library("tximportData")
library("rhdf5")
samples <- read.table(file = "/home/pacoh/Dropbox/howest/server/exampledata/SraRunTable.csv",
                       sep = ",", 
                       header = TRUE, 
                       stringsAsFactors = TRUE)
dir <- "/home/pacoh/Dropbox/howest/server"

## Importing Kallisto with TSV files
files <- file.path(dir, "kallisto", paste0(samples$Run,"_1"), "abundance.tsv")
names(files) <- paste0(samples$Run,"_1")
txi.kallisto.tsv <- tximport(files, type = "kallisto", 
                             txOut = TRUE)
head(txi.kallisto.tsv$counts, n = 10)

## METHOD 1: "original counts and offset" 
## --> for edgeR you need to assign a matrix to y$offset
## Example creating DGEList for use with edgeR
# library(edgeR)
cts <- txi.kallisto.tsv$counts
normMat <- txi.kallisto.tsv$length
normMat
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
# Assign matrix to y$offset (DGEobject$offset)
y$offset <- t(t(log(normMat)) + o)
dim(y)
# [1] 39627 (=transcripts) 6 (=samples)
# y is now ready for estimate dispersion functions --> see edgeR User's Guide
# Given a table of counts or a DGEList object, the qCML common dispersion 
# and tagwise dispersions can be estimated using the estimateDisp() function

## METHOD 2: "bias corrected counts without an offset"
## --> uses the gene-level count matrix txi$counts
## Because limma-voom does not use the offset matrix stored in y$offset
## it is recommended to use the scaled counts generated from abundances
## either "scaledTPM" or "lengthScaledTPM" --> tximport parameter: 
## countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM")
txi <- tximport(files, type = "kallisto", 
                txOut = TRUE,
                countsFromAbundance = "scaledTPM")
# library(limma)
y <- DGEList(txi$counts)
y <- calcNormFactors(y)
# design <- model.matrix()
# v <- voom(y, design)
# v is then ready for lmFit()



################################################################################
### OPTIONAL: dimensionality reduction techniques (instead of MDS)
################################################################################
## PCA or t-SNE
set.seed(1)
samples <- colnames(log.cpmF.norm)
genes <- rownames(log.cpmF.norm)
pch.groups <- as.numeric(metaData$genotype)
pch.groups # pch num: 2 2 2 1 1 1
metaData$genotype # WT WT WT KO KO KO
col.group # colors: "#FB9A99" "#FB9A99" "#FB9A99" "#B2DF8A" "#B2DF8A" "#B2DF8A"
# Run tSNE
set.seed(1)
tsne.out <- Rtsne(t(log.cpmF.norm), 
                  dims = 2, 
                  initial_dims = 30, 
                  perplexity = 5)
# --> error: perplexity is too large
# --> when only few samples (and filtered genes)
# Plot the results
par(mfrow=c(1,1))
plot(tsne.out$Y,
     col = col.group,
     pch = pch.groups,
     main = "tSNE")
legend("topleft",
       legend = metaData$genotype,
       col = col.group,
       pch = pch.groups,
       cex = 0.8,
       bty = 'n')
################################################################################
################################################################################
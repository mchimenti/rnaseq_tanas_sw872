## Analysis of Munir Tanas RNA-Seq data in human SW872 cells
## Date: 03.13.19
## Author: Michael Chimenti
## Organism: hg38 / human 
## Aligners: hisat2 / salmon
## Design: Condition vs Empty Vector
## Reps: 3

##########
## Imports
##########

source("https://bioconductor.org/biocLite.R")
biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/tanas_lab/project_rnaseq_SW872_feb2019")

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 > PCAExplorer
#######################################


samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample
#samples$batch <- as.factor(samples$batch)

files <- file.path(getwd(), samples$sample, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read.csv(file.path(getwd(), "tx2gene.csv"), header = FALSE, as.is = c(1:2)) 

tx2gene$V1 <- tx2gene$V1 %>% 
  strsplit(split = '.', fixed = TRUE) %>%
  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ cond)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

##---------------launch PCA Explorer on dds object 
anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rldTxi <- rlog(ddsTxi, blind=FALSE)

rld_tab <- as.data.frame(assay(rldTxi))
rld_tab$gene <- anno[row.names(rld_tab), "gene_name"]
write.csv(x=rld_tab, file = "rlog_transform_counts_genes.csv")


pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

pcaExplorer::pcaplot(ddsTxi, intgroup = 'cond', ellipse = FALSE, text_labels = FALSE)
pcaExplorer::pcaplot(rldTxi, intgroup = "cond", ellipse = FALSE, text_labels = FALSE)

plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05, ylim = c(-10,10))

##############
## DE analysis
##############

## CAMTA1 vs EV_sub
res_camta1 <- results(ddsTxi, contrast = c("cond","CAMTA1","EV_sub"))
res_camta1 <- na.omit(res_camta1)  #drop NA rows
res_camta1_sig <- res_camta1[res_camta1$padj < 0.05 & res_camta1$baseMean > 5.0,]
res_camta1_ord <- res_camta1_sig[order(res_camta1_sig$padj),]
res_camta1_ord$ext_gene <- anno[row.names(res_camta1_ord), "gene_name"]

png("camta1_DE_volcano_zoom.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_camta1_ord, main = "Volcano Plot: DE genes across camta1 vs ev_sub", lfcthresh=2, sigthresh=0.05, textcx=.6, xlim=c(-4, 4), ylim = c(3,50))
dev.off()

plotCounts(ddsTxi,"ENSG00000128739", intgroup = "cond", normalized = TRUE)

p <- degPlot(ddsTxi, xs = "cond", res = res_camta1_ord, n = 6, group = "cond", 
        groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p

## TAZ4SA vs EV
res_taz <- results(ddsTxi, contrast = c("cond","TAZ4SA","EV"))
res_taz <- na.omit(res_taz)  #drop NA rows
res_taz_sig <- res_taz[res_taz$padj < 0.05 & res_taz$baseMean > 5.0,]
res_taz_ord <- res_taz_sig[order(res_taz_sig$padj),]
res_taz_ord$ext_gene <- anno[row.names(res_taz_ord), "gene_name"]

png("taz_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_taz_ord, main = "Volcano Plot: DE genes across taz vs ev", lfcthresh=2, sigthresh=0.05, textcx=.4, xlim=c(-8, 8), ylim = c(3,100))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_taz_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p

## TAZ/CAMTA1 vs EV
res_tc <- results(ddsTxi, contrast = c("cond","TC","EV"))
res_tc <- na.omit(res_tc)  #drop NA rows
res_tc_sig <- res_tc[res_tc$padj < 0.05 & res_tc$baseMean > 5.0,]
res_tc_ord <- res_tc_sig[order(res_tc_sig$padj),]
res_tc_ord$ext_gene <- anno[row.names(res_tc_ord), "gene_name"]

png("TC_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_tc_ord, main = "Volcano Plot: DE genes across tc vs ev", lfcthresh=2, sigthresh=0.05, textcx=.6, xlim=c(-8, 8), ylim = c(3,120))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_tc_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p


## TFE3 vs EV
res_tfe3 <- results(ddsTxi, contrast = c("cond","TFE3","EV"))
res_tfe3 <- na.omit(res_tfe3)  #drop NA rows
res_tfe3_sig <- res_tfe3[res_tfe3$padj < 0.05 & res_tfe3$baseMean > 5.0,]
res_tfe3_ord <- res_tfe3_sig[order(res_tfe3_sig$padj),]
res_tfe3_ord$ext_gene <- anno[row.names(res_tfe3_ord), "gene_name"]

png("tfe3_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_tfe3_ord, main = "Volcano Plot: DE genes across tfe3 vs ev", lfcthresh=2, sigthresh=0.05, textcx=0.6, xlim=c(-10, 10), ylim = c(3,140))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_tfe3_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p

## YAP5SA vs EV
res_yap5sa <- results(ddsTxi, contrast = c("cond","YAP5SA","EV"))
res_yap5sa <- na.omit(res_yap5sa)  #drop NA rows
res_yap5sa_sig <- res_yap5sa[res_yap5sa$padj < 0.05 & res_yap5sa$baseMean > 5.0,]
res_yap5sa_ord <- res_yap5sa_sig[order(res_yap5sa_sig$padj),]
res_yap5sa_ord$ext_gene <- anno[row.names(res_yap5sa_ord), "gene_name"]

png("yap5sa_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_yap5sa_ord, main = "Volcano Plot: DE genes across yap5sa vs ev", lfcthresh=2, sigthresh=0.05, textcx=.5, xlim=c(-8, 8), ylim = c(3,120))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_yap5sa_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p

## YAP/TFE vs EV
res_yt <- results(ddsTxi, contrast = c("cond","YT","EV"))
res_yt <- na.omit(res_yt)  #drop NA rows
res_yt_sig <- res_yt[res_yt$padj < 0.05 & res_yt$baseMean > 5.0,]
res_yt_ord <- res_yt_sig[order(res_yt_sig$padj),]
res_yt_ord$ext_gene <- anno[row.names(res_yt_ord), "gene_name"]

png("YT_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_yt_ord, main = "Volcano Plot: DE genes across YT vs ev", lfcthresh=2, sigthresh=0.05, textcx=.5, xlim=c(-8, 8), ylim = c(3,120))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_yt_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p

## TAZ/CAMTA1 vs CAMTA1
res_tc_camta1 <- results(ddsTxi, contrast = c("cond","tc_camta1","EV"))
res_tc_camta1 <- na.omit(res_tc_camta1)  #drop NA rows
res_tc_camta1_sig <- res_tc_camta1[res_tc_camta1$padj < 0.05 & res_tc_camta1$baseMean > 5.0,]
res_tc_camta1_ord <- res_tc_camta1_sig[order(res_tc_camta1_sig$padj),]
res_tc_camta1_ord$ext_gene <- anno[row.names(res_tc_camta1_ord), "gene_name"]

png("tc_camta1_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_tc_camta1_ord, main = "Volcano Plot: DE genes across tc_camta1 vs ev", lfcthresh=2, sigthresh=0.05, textcx=.5, xlim=c(-8, 8), ylim = c(3,120))
dev.off()

p <- degPlot(ddsTxi, xs = "cond", res = res_tc_camta1_ord, n = 6, group = "cond", 
             groupLab = "Condition")

p <- p + theme_dark()
p <- p + theme(axis.text.x = element_text(size=8, angle=45, face="bold"))
p


## TAZ/CAMTA1 vs TAZ





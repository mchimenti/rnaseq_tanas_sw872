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
plotMA(object = ddsTxi, alpha = 0.05)

##############
## DE analysis
##############

## Genotype
res_geno <- results(ddsTxi2, contrast = c("geno","KO","WT"))
res_geno <- na.omit(res_geno)  #drop NA rows
res_geno_sig <- res_geno[res_geno$padj < 0.05 & res_geno$baseMean > 5.0,]
res_geno_ord <- res_geno_sig[order(res_geno_sig$padj),]
res_geno_ord$ext_gene <- anno[row.names(res_geno_ord), "gene_name"]

png("geno_DE_volcano.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_geno_ord, main = "Volcano Plot: DE genes across genotypes, KO vs. WT (Male responses)", lfcthresh=1.0, sigthresh=0.1, textcx=.6, xlim=c(-12, 12), ylim = c(3,25))
dev.off()

degPlot(dds = ddsTxi2, res = res_sex_ord, n = 9, xs = "sex")


# DEseq

#Commenting this out since I already installed this software
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install("ggplot2")
#BiocManager::install("ggrepel")
#install.packages("UpSetR")
require(devtools) || install.packages("devtools",quietly = T)
require(MoMAColors) || devtools::install_github("BlakeRMills/MoMAColors")


#Loading the libraries (required for this analysis)
suppressMessages(library("DESeq2",quietly = T))
suppressMessages(library("apeglm",quietly = T))
suppressMessages(library("dplyr",quietly = T))
#library("MoMAColors")
suppressMessages(library("ggplot2",quietly = T))
suppressMessages(library("ggrepel",quietly = T))
suppressMessages(library("RColorBrewer",quietly = T))
suppressMessages(library("rlang",quietly = T))
suppressMessages(library("purrr",quietly = T))
suppressMessages(library("gplots", quietly = T))
suppressMessages(library("UpSetR", quietly = T))

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  stop("Six arguments are required", call.=FALSE)
} else if (length(args)==6) {
  geneCountFile = args[1]
  colDataFile = args[2]
  alpha = args[3]
  resultsPath = args[4]
  treatVsCtrl = args[5]
  ctrl = args[6]
}

#geneCountFile = 'combined_gene_count.csv'
#colDataFile = 'colData.csv'
#alpha = 0.05
#treatVsCtrl = 'treatVsControl.csv'
#resultsPath = 'results'
#ctrl = "CL"


# Printing filepaths
geneCountFile
colDataFile
alpha = as.numeric(alpha)
alpha

#Loading counts into a matric
countData <- as.matrix(read.csv(geneCountFile, row.names=1, stringsAsFactors = TRUE))

#Loading phenotype data 
coldata <- read.csv(colDataFile,row.names=1,stringsAsFactors = TRUE)

#Loading treatment vs contrl data - must be spelled the same as the DE
ctrlData <- read.csv(treatVsCtrl,row.names=1)

#Formatting fixes (matching columns to rows, in same order) - 
#first two are checking if they are all there and they match
if (all(rownames(coldata) %in% colnames(countData))) {
  print('All count columns have a corresponding row in colData')
} else {
  print('Your coldata file is incomplete. All count columns must have 
        a corresponding row in colData')
}
if (all(rownames(coldata) == colnames(countData))) {
  print('Rows and Columns are in matching order')
} else {
  coldata <- coldata[dput(colnames(countData)),]
  print('Rearranged Row and Columns so the order matches.')
}


# Changes any NA data in count matrix to 0
countData[is.na(countData)] <- 0

#loading data into Deseq2
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = coldata,
  design = ~ DE)
dds$DE <- as.factor(dds$DE)
dds$Name <- as.factor(dds$DE)

class(alpha)

dds <- DESeq(dds)
dds$DE <- relevel(dds$DE, ctrl)
dds <- DESeq(dds)
resultsNames(dds)

# empty dataframe to save l2FC (log2FoldChange) and q-value (padj)
# for comparing comparisons
# Q-value is FDR (False Discovery Rate) when running many comparisons 
# DESEq2 calculates padj (can be used as q-value) using the BH method (Benjamini-Hochberg) by default and this seems to be
# a reasonable method.
fc_q_table <-data.frame(gene_id = rownames(countData))
fc_q_table$l2FC <- 0
fc_q_table$qvalue <- 0
qvalue_cols <- vector()
fc_cols <- vector()

for (row in as.integer(rownames(ctrlData))) {
  #run comparisons
  treat = ctrlData$treatment[row]
  control = ctrlData$control[row]
  print(paste("Comparison of ", treat, " vs ", control, sep="" ))
  res <- results(dds, alpha=alpha, contrast=c("DE", treat, control))
  res <- res[order(res$padj),]
  name = paste(treat,"vs",control,sep="_")
  resultsTxtFile = paste(resultsPath, "/", name, "_results.txt", sep="")
  sink(resultsTxtFile, append = FALSE)
  print(paste("Comparison of ", treat, " vs ", control, sep="" ))
  summary(res)
  sink()
  
  # getting l2FC and padj values to be able to do comparisons 
  # outside this loop.
  resTable <- data.frame(res)
  resTable$gene_id <- row.names(resTable) 
  resTable <- resTable %>%
    select(gene_id,log2FoldChange,padj)
  
  names(resTable) <- c("gene_id", "l2FC", "qvalue")
  
  fc_q_table <- full_join(fc_q_table, resTable, by="gene_id", suffix = c("",paste(".",name,sep="")))
  qvalue_cols <- append(qvalue_cols, paste("qvalue.",name,sep=""))
  fc_cols <- append(fc_cols, paste("l2FC.",name,sep=""))
  
  # export to csv
  file = paste(resultsPath, "/", "DEseq2",name,"_res",alpha,".csv",sep="")
  write.csv(res, file = file)
  
  #Volcano Plot
  DE <- read.csv(paste(resultsPath, "/", "DEseq2",name,"_res",alpha,".csv",sep=""), row.names = 1)
  #plotDE <- plotDE[complete.cases(plotDE),]
  plotDE <- DE %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
  plotDE$Name <- rownames(plotDE)
  plotDE <- plotDE[complete.cases(plotDE),]
  cols <- densCols(plotDE$log2FoldChange, plotDE$pvalue)
  summary(plotDE)
  cols[plotDE$gene_type=='Up-regulated']<-"#FF9900"
  cols[plotDE$gene_type=='Down-regulated']<-"#9933CC"
  cols[plotDE$gene_type=='Unchanged']<-"#999999"
  sizes <- c("Up-regulated" = 3, "Down-regulated" = 3, "Unchanged" = 2)
  alphas <- c("Up-regulated" = 1, "Down-regulated" = 1, "Unchanged" = 0.5)
  sig <- plotDE[plotDE$pvalue<=0.05&plotDE$gene_type!='Unchanged',]
  if (nrow(sig) == 0) {
    write(paste("No unchanged genes. Can't do Volcano plot for ", name, sep=""), file = resultsTxtFile, append=TRUE)
  } else {
    vptitle <- paste("Volcano plot for ", name, sep="")
    ggplot(plotDE, aes(x=log2FoldChange, y=-log10(pvalue), size = gene_type, alpha = gene_type))+
      geom_point(col = cols)+
      #scale_y_continuous(expand = c(0,0), limits = c(0, 8.5))+
      geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
      geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
      geom_label_repel(data = sig, aes(label = Name), force = 2, nudge_y = 1) +
      scale_size_manual(values = sizes) + # Modify point size
      scale_alpha_manual(values = alphas) + # Modify point transparency
      ggtitle(vptitle) +
      xlab("Effect size: log2(fold-change)") +
      ylab("-log10(p-value)")+
      theme_bw()+
      theme(legend.position = "none")
    volcanoFile = paste(resultsPath, "/",name,"_volcano_plot.png",sep="")
    ggsave(filename = volcanoFile)
  }
  print(paste("Analysis of ", name, " is complete!", sep=""))
}

for (i in 1:length(resultsNames(dds))) {
  if (i != 1) {
    # prepare data for maplot
    resLFC <- lfcShrink(dds, coef =i, type="apeglm")
    resLFC <- resLFC[order(resLFC$padj),]

    #export MAplot
    pdfName = paste(resultsPath, "/", resultsNames(dds)[i],"_maplot.pdf",sep="")
    pdf(pdfName, width = 10, height = 10)
    plot.window(plotMA(resLFC, alpha=alpha, ylim=c(-2,2), xlim=c(1e-01,1e+07), main = resultsNames(dds)[i]), ylim=c(-2,2), xlim=c(1e-01,1e+07), log="")
    dev.off() 
  } else {}
}

print("Before Norm Counts")
# export normalized counts
dds2 <- estimateSizeFactors(dds)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized=TRUE)
normalizedFile = paste(resultsPath, "/", "normalized_counts.txt", sep="" )
write.table(normalized_counts, file=normalizedFile, sep="\t", quote=F, col.names=NA)

print("Before PCA")
# PCA plots
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("DE"), returnData=TRUE)
percentVarCondition <- round(100 * attr(pcaData, "percentVar"))
PCA1 <- paste("PCA1 (", percentVarCondition[1],"% variance)", sep="")
PCA2 <- paste("PCA2 (", percentVarCondition[2],"% variance)", sep="")
PCATitle <- paste("PCA", sep="")
PCA <- ggplot(pcaData, aes(x=PC1, y=PC2)) +
  theme_bw()+
  ggtitle(PCATitle)+
  theme(legend.position = "right") +
  xlab(PCA1) + 
  ylab(PCA2) +
  guides(size = "none")+
  geom_point(aes(colour = group, shape = group, size = 3), show.legend = T) +  #different groups have different color
  geom_text(aes(label=name),hjust="inward", vjust="inward", size = 2, position=position_jitter(width=0,height=4))+  #adding labels to the points
  scale_color_moma_d("OKeeffe", direction = -1)
pcaFile <- paste(resultsPath, "/PCAplot.png",sep="")
ggsave(filename = pcaFile, plot = PCA, width=6, height=4)

print("Before Heatmap1")
#plotting all DE genes
#create a subtable of interest
fc_q_table1 <- subset(fc_q_table, select = -c(l2FC,qvalue))
rownames(fc_q_table1) <- fc_q_table1$gene_id

for (num in 1:length(qvalue_cols)) {
  col1 = qvalue_cols[num]
  fc_q_table1[,col1] <- case_when(fc_q_table1[,col1] <= alpha ~ "*", TRUE ~ "")
  fc_q_table1$combo <- paste(fc_q_table1$combo, fc_q_table1[,col1])
}

fc_q_table1$combo <- gsub(' ','',fc_q_table1$combo)

fc_q_table1 <- fc_q_table1 %>% filter(grepl('\\*+', combo))

fc_q_table1 <- subset(fc_q_table1, select = -c(combo))

#subset table to FC and qvalues
FC <- subset(fc_q_table1, select = fc_cols)
rownames(FC) <- FC$gene_id
qvalue <- subset(fc_q_table1, select = qvalue_cols)
rownames(qvalue) <- qvalue$gene_id

FC <- FC[complete.cases(FC),]

colors <- seq(-8,8,0.2)

colors_flash <- moma.colors("Flash", n=80, type="continuous")

pdfName = paste(resultsPath, "/heatmap_DE.pdf",sep="")
pdf(pdfName, width = 10, height = 10)
heatmap.2(as.matrix(FC), key=T, keysize=1, key.xlab= "log2(FC)", 
          trace="none", symm=F, symkey=F, cex.main=0.3, main="DE genes", 
          breaks=colors, scale="none", col=colors_flash, density.info="none", 
          cellnote = qvalue, notecol="black", 
          dendrogram = "both", notecex = 0.8, Rowv=T, cexRow=0.2,  cexCol=0.8, Colv=T)


print("Before Heatmap2")
#To plot Expression results as a Heatmap:
#load normalized counts that we exported earlier in this script
expression <- read.table(normalizedFile, header=TRUE, row.names=1, sep="\t")

expression <- expression %>% filter(row.names(expression) %in% row.names(fc_q_table1))

#export significant normalized genes
normalizedFileSig = paste(resultsPath, "/", "normalized_counts_significant_",alpha, ".txt", sep="" )
write.table(normalized_counts, file=normalizedFileSig, sep="\t", quote=F, col.names=NA)

# setting colors
colors = seq(0,10,0.1)
colors_kippenberger <- moma.colors("Kippenberger", n=100, type="continuous")

# creating and exporting heatmap to a pdf
pdfName = paste(resultsPath, "/heatmap_expression.pdf",sep="")
pdf(pdfName, width = 10, height = 10)
heatmap.2(as.matrix(expression), key=T, keysize= 1, key.xlab= "Log2(FPKM)", 
                                 trace="none", symm=F, symkey=F, cex.main=0.3, main = "Gene expression", 
                                 breaks=colors, scale="none", col=colors_kippenberger , density.info="none",
                                 dendrogram = "both", notecex = 0.3, Rowv=T, cexRow=0.4,  cexCol=0.8, Colv=T,  margin=c(10, 10))



print("Before UpsetPlot")
# Make an UpsetPlot

#transforming data to work for upsetplot
upsetPlot <- as.data.frame(qvalue)
colnames(upsetPlot) <- gsub("qvalue.","",colnames(upsetPlot))
for (num in 1:length(qvalue_cols)) {
  upsetPlot[,num] <- case_when(upsetPlot[,num] == "*" ~ 1, TRUE ~ 0)
}

#plotting UpSet plot
pdfName = paste(resultsPath, "/upset_plot.pdf",sep="")
pdf(pdfName, width = 12, height = 10, onefile=FALSE)
upset(upsetPlot, sets = colnames(upsetPlot), 
      order.by = "freq",
      keep.order = T,
      empty.intersections = T,
      number.angles = 0, point.size = 4, line.size = 1, 
      mainbar.y.label = "DE genes", sets.x.label = "DE genes",
      text.scale = c(1.5, 1.5, 1.2, 1.5, 1.5, 1.8), mb.ratio = c(0.7, 0.3))
dev.off() 



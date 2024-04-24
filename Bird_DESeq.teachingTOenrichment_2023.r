####  DESeq2 Script for Differential Gene Expression Analysis in 
      # Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported


## Use the Session menu to set working directory To Source File Directory

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
  ### How you inport this will depend on what your final output was from the mapper/counter that you used.
  ## this works with output from PrepDE.py from Ballgown folder.

countdata <- as.matrix(read.csv("Birdgene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)

#### IF necessary,depending on what program made the count matrix. Remove the unwanted row (the * and zero row)
#countdata<- countdata[-1,c(-1)]
# OR   Remove the unwanted column (length)
#countdata<- countdata[,c(-1)]
#dim(countdata)
#head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata

coldata <-(read.table("phenobird_Data.txt", header=TRUE, row.names=1))
# Assuming your data frame is named 'df' and you want to rename the first two columns
colnames(coldata)[0] <- "Treatment"
colnames(coldata)[1] <- "type"

dim(coldata)
head(coldata)

library(dplyr)

# Assuming your data frame is named 'df' and you want to rename the first two columns
coldata <- coldata %>%
  rename(Treatment = name_of_first_column, # replace with the actual column name
         Type = name_of_second_column)     # replace with the actual column name

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))
colnames(coldata)

colnames(coldata)[colnames(coldata) == "type"] <- "Treatment"

## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata, design = ~Treatment)
coldata <- as.data.frame(coldata)

#look at it
dds

#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
## You can play around with this number to see how it affects your results!
dds <- dds[ rowSums(counts(dds)) > 20, ]
# look.  How many genes were filtered out?
dds

## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds$condition <- factor(dds$Treatment, levels=c("control","heat-treatment"))

######     1.4 Differential expression analysis
### Question 2. Look at the manual - what is happening at this point?
dds <- DESeq(dds)
res <- results(dds)
res

###  Question 3. What does each column mean?
# We can order our results table by the smallest adjusted p value:
  resOrdered <- res[order(res$padj),]
  resOrdered
# We can summarize some basic tallies using the summary function the default is p<0.1.
  summary(res)
#How many adjusted p-values were less than 0.1?
  sum(res$padj < 0.01, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
  res05 <- results(dds, alpha=0.01)
  summary(res05)
  #sum(res05$padj < 0.05, na.rm=TRUE)
###    1.5.1 MA-plot
  ##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
  ## Points will be colored red if the adjusted p value is less than 0.1. 
  ## Points which fall out of the window are plotted as open triangles pointing either up or down
  plotMA(res, alpha = 0.1, colNonSig = "grey", colSig = "red",colLine = "grey4",main="DESeq2", ylim=c(-8,8))
  
  #Add a legend to the plot
  # Note: Adjust the legend parameters according to what you specifically want to denote in your plot.
  legend("topright",  # Position of the legend
         legend=c("Significantly down", "Significantly up", "line"),  # Legend text
         col=c("grey", "red", "grey4"),  # Colors corresponding to your legend text
         pch=20, cex = 0.75 )  # Symbol type
  #After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
  # One can then recover the gene identiers by saving the resulting indices:
  idx <- identify(res$baseMean, res$log2FoldChange)
    # after selecting a gene. You need to press escape to move on
  rownames(res)[idx]
    
##  1.5.2 Plot counts - sanity check!
  
  # You can select the gene to plot by rowname or by numeric index.
  plotCounts(dds, gene="TMEM86A|TMEM86A", intgroup="Treatment")
  # You can plot the gene with th lowest adjuated P-value
  plotCounts(dds, gene=which.min(res$padj), intgroup="Treatment")
  dds
  
  ##  Write your results to a file 
  write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  
  
  ## 2.1.2 Extracting transformed values
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  head(assay(rld), 3)
  
  ### Heatmap of the count matrix
  #library("genefilter")
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  
  library("pheatmap")
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  colnames(colData(vsd))
  anno <- as.data.frame(colData(vsd)[, c("Treatment", "type")])
  df <- as.data.frame(colData(dds)[,c("Treatment","type")])
    pheatmap(mat, annotation_col = anno)
  
  
 ## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
#  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
#  nt <- normTransform(dds) # defaults to log2(x+1)
#  log2.norm.counts <- assay(nt)[select,]
#   df <- as.data.frame(colData(dds)[,c("treatment","type")])
#  pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#          cluster_cols=TRUE, annotation_col=df)
  
# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df) 

  #2.2.2 Heatmap of the sample-to-sample distances
  sampleDists <- dist(t(assay(rld)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$treatment)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           main = "Heatmap of Sample-to-Sample Distances")
  
 # 2.2.3 Principal component plot of the samples
  plotPCA(rld, intgroup=c("Treatment"))
  
############ Preparing Data for GSEA and Cytoscape.  #############
  
### Merge 'gene names' with DGE results by Gene Model
  
## Import Annotation file with results from Blast to databases
Anno <- read.csv("Birdgene_count_matrix.csv", stringsAsFactors = FALSE, na.strings=c(""))
summary(Anno)
dim(Anno)
  
## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
DGEresults <- read.csv("Birdgene_count_matrix.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)

## Rename first column so it matches "gene_id" in annotation file
names(DGEresults)[1]<- "gene_id" 

#Merge anno with DGE results
DGE_Anno <- merge(Anno,DGEresults,by="gene_id",all.y = FALSE)
dim(DGE_Anno)
summary(DGE_Anno)
colnames(DGE_Anno)
############################# Make ranked list for GSEA ####################
   ## example aa <- within(resOrdered, z <- x + y - 2)
DGE_Anno_Rank <-  within(DGE_Anno, rank <- sign(DGE_Anno$ERR2984230.x) * -log10(DGE_Anno$ERR2984230.y))
DGE_Anno_Rank 
write.table(as.data.frame(DGE_Anno_Rank), file="DGE_Anno_Rank.txt", quote = FALSE, row.names=FALSE, sep = "\t")  

 #subset the results so only Gene Name and rank
DGErank = subset(DGE_Anno_Rank, select = c(gene_id, rank))
DGErank

#sebset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)

write.table(as.data.frame(DGErank_withName), file="DGErank_withName.rnk.txt", quote = FALSE, row.names=FALSE, sep = "\t")  
# If 'data' is a vector of the text lines you provided
# Convert to character vector if it's a factor or a similar non-character data type
# Assuming DGErank_withName is a data frame and the column you want to split is called 'gene_info'


DGErank_withName$gene_id  <- as.character(DGErank_withName$gene_id)
split_dataDGE <- strsplit(DGErank_withName$gene_id, "\\|")

# Assuming the result of the split is two elements, gene name and another identifier
new_data <- do.call(rbind, lapply(split_dataDGE, function(x) data.frame(gene_id=x[1], Identifier=x[2])))
new_data
# Add the numeric values to the new data frame
new_data$rank <- DGErank_withName$rank
#delete the identifier column
new_data$Identifier <- NULL

write.table(as.data.frame(new_data), file="Birddata_withName.rnk.txt", quote = FALSE, row.names=FALSE, sep = "\t")  

### IF NEEDED....remove the "gene-" from row names
  ## https://stackoverflow.com/questions/39897155/replace-and-remove-part-of-string-in-rownames/39897315  "URS000075AF9C-snoRNA_GTATGTGTGGACAGCACTGAGACTGAGTCT"    to   "snoRNA"
  ## We can use gsub to match one of more characters that are not a - ([^-]+) from the start (^) of the string followed by 
  ## a - or (|) one or more characters that are not an underscore ([^_]+) until the end of the string ($)  and replace it with blanks ("").
  ## gsub("^[^-]+-|_[^_]+$", "", rownames(df))

gene <-gsub("^[^-]+-", "", rownames(DGErank))
DGErankIDs  <-cbind(gene,DGErank)
head(DGErankIDs)
#summary(DGErankIDs)  

write.csv(as.data.frame(DGErankIDs), file="DGErankIDs.csv", row.names=FALSE)  

####  We also need the normalized expression DATA
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)
# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
gene <-gsub("^[^-]+-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(gene,NormTransExp)
head(NormTransExpIDs)

write.csv(as.data.frame(NormTransExpIDs), file="NormTransExpressionData.csv", row.names=FALSE)  


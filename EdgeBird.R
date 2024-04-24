#### Install the EdgeR package if you have not already
## try http:// if https:// URLs are not supported

setwd("/Users/temitopefolorunso/Desktop/Functional Genomics/Bird_Project")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

## Load the edgeR library 
library(limma)
library(edgeR)
library(dplyr)
library(stringr)
## Use the Session menu to set working directory To Source File Directory

# Reading in the feature count file as "counts.df"
countdata1 <- read.csv("2_4WC_3HBirdgene_count_matrix.csv", row.names="gene_id")
# Clean the gene IDs (currently your row names)
clean_gene_ids <- str_replace(rownames(countdata1), "gene-", "")
clean_gene_ids <- str_remove(clean_gene_ids, "\\|.+") 

# Update the row names with the cleaned versions
rownames(countdata1) <- clean_gene_ids


# Printing the start of the counts.df object in R...
dim(countdata1)
head(countdata1)
head(clean_gene_ids)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table(file.choose(), header=TRUE,row.names = 1))
dim(coldata)
head(coldata)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata1))
countdata1 <- countdata1[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata1))

# Subsetting gene counts according to experimental condition
counts_control.df  <- countdata1[,c("ERR2984226", "ERR2984227","ERR2984246")]
counts_treated.df <- countdata1[,c("HT_3hr_1", "HT_3hr_2", "HT_3hr_4")]

# Printing the structure of the gene counts set and subsets
str(countdata1)
str(counts_control.df)
str(counts_treated.df)

#Relative standard deviation (RSD), in edgeR, is calculated as the standard
#deviation of a gene's expression values, across samples, divided by the gene's 
#mean expression. We use RSD filtering to identify and remove genes with highly 
#variable expression across samples.  This filtering is important because genes 
#with high RSD may be unduly influenced by technical noise or sample-specific outliers, 
#making it harder to identify true biological differences. 
#By filtering out these genes, we improve the reliability and statistical power of our differential expression analysis.

# Defining function "RSD.test()"
RSD.test <- function(dataframe){
  # This function tests whether the relative standard deviation (RSD) is less
  # than or equal to one for each row in a data frame.
  # It adds the result to a new variable in the data frame called "RSD.test".
  # For a given row, if data.frame$RSD.test is TRUE, that row has an RSD less
  # than or equal to one, i.e. RSD <= 1.
  # If data.frame$RSD.test is FALSE, that row has an RSD outside of this range.
  RSD_tests = dataframe[,1]
  for (row_index in 1:nrow(dataframe)){
    row = as.numeric(dataframe[row_index,])
    RSD = sd(row) / mean(row)
    RSD_tests[row_index] = RSD <= 1 || is.na(RSD)
  }
  dataframe$RSD.test <- as.factor(RSD_tests)
  levels(dataframe$RSD.test) <- c(FALSE, TRUE)
  return(dataframe)
}

# Applying RSD.test() to gene count subsets
counts_control.df  <- RSD.test(counts_control.df)
counts_treated.df <- RSD.test(counts_treated.df)

# Printing the structure of the gene counts subsets
str(counts_control.df)
str(counts_treated.df)

# Creating list of genes which failed RSD test
RSD_failed_genes <- rownames(counts_control.df[
  which(counts_control.df$RSD.test == FALSE),])
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_treated.df[
  which(counts_treated.df$RSD.test == FALSE),]))
RSD_failed_genes <- unique(RSD_failed_genes)
length(RSD_failed_genes)

# Filtering gene counts
filtered_counts.df <- countdata1[
  which(!rownames(countdata1) %in% RSD_failed_genes),]

# Printing the structure of the filtered gene counts
str(filtered_counts.df)

# Checking that gene counts were correctly filtered
nrow(countdata1) - length(RSD_failed_genes) == nrow(filtered_counts.df)

# Creating a DGEList object using the filtered gene counts
counts.DGEList <- DGEList(counts = filtered_counts.df,
                          genes = rownames(filtered_counts.df))

#If you dont want to use RSD test, then creating a DGEList object using the original unfiltered gene counts
#counts.DGEList <- DGEList(counts = countdata, genes = rownames(countdata))

# Printing the design table
print(coldata)


# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df) == rownames(coldata))

# for unfiltered(with noRSD filtering) data use this
#summary(colnames(countdata) == rownames(coldata))

# Add grouping information to DGEList object
counts.DGEList$samples$group <- as.factor(coldata$type)

# Printing counts.DGEList
counts.DGEList
dim(counts.DGEList)
# save into a csv file 
# Assuming counts.DGEList is already a data frame
write.csv(counts.DGEList, "2_4_3Hcounts_DGEList.csv", row.names = FALSE)

# Creating an object to filter genes with low expression
counts.keep <- filterByExpr(counts.DGEList)
summary(counts.keep)

# Filtering lowly expressed genes
filtered.DGEList <- counts.DGEList[counts.keep, , keep.lib.sizes = FALSE]
dim(filtered.DGEList)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in counts.keep
length(counts.keep[counts.keep == TRUE]) == dim(filtered.DGEList)[1]

# Printing the normalisation factors for the libraries
filtered.DGEList$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
norm.DGEList <- calcNormFactors(filtered.DGEList)
norm.DGEList$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition <- coldata$type
disp.DGEList <- estimateDisp(norm.DGEList, design = model.matrix(~condition))
disp.DGEList

# Assuming counts.DGEList is already a data frame
write.csv(disp.DGEList, "2_4_3Hdisp.DGEList.csv", row.names = FALSE)

# Exact tests for differences between experimental conditions
std_treated.DGEExact <- exactTest(disp.DGEList, pair = c("control",
                                                         "heat-treatment")) 

# Extracting most differentially expressed genes from exact tests
std_treated.topTags <- topTags(std_treated.DGEExact,n = nrow(std_treated.DGEExact$table))

# Printing the most differentially expressed genes
std_treated.topTags
# Assuming counts.DGEList is already a data frame
write.csv(std_treated.topTags, "2_4_3Hstd_treated.topTags.csv", row.names = FALSE)

# extract significant differentially expressed genes, sort, & write to csv
resOrdered1 <-std_treated.topTags$table[order(std_treated.topTags$table$FDR),]
resSig <- resOrdered1[resOrdered1$FDR<0.05,]
resSig
print(resOrdered1)

out <- resOrdered1 %>% dplyr::select(logFC, logCPM, PValue, FDR) %>% dplyr::rename(meanExpr = logCPM, pval = PValue, adj.pval = FDR)
write.csv(as.data.frame(out),file = paste0("/Users/temitopefolorunso/Desktop/Functional Genomics/Bird_Project/edgeR3.csv")) #write results to a new csv


# Create the MA-plot 
edgeR::plotMD.DGEExact(std_treated.DGEExact)


#Calculating the total number of differentially expressed genes at FDR< 0:05 
de1 <- decideTestsDGE(std_treated.DGEExact , adjust.method="BH", p.value=0.1)
summary(de1)

#..............heatmap construction.........................................
library("RColorBrewer")
library("pheatmap")

logcpm <- cpm(norm.DGEList, log=TRUE) 
## Making heatmap of count matrix for significantly expressed genes
#Prepare Expression Matrix by subsetting the log-CPM matrix to include only significant genes
topVarGenes <- order(apply(logcpm, 1, var), decreasing=TRUE)[1:20] 
topVarGenes
expression_matrix <- logcpm[topVarGenes, ]
expression_matrix
#expression_matrix <- logcpm[rownames(resSig), ]

# Choose Color Palette
color_palette <- brewer.pal(9, "RdYlBu")  # Adjust the palette name if desired

treatment_groups <-  coldata$type 
annotation_df <- data.frame(Treatment = treatment_groups) 
rownames(annotation_df) <- colnames(expression_matrix)

# Generate Heatmap
pheatmap(logcpm,
         color = color_palette,
         scale = "row",          # Scale expression values by row
         cluster_rows = TRUE,    
         cluster_cols = TRUE,   
         show_rownames = TRUE, 
         show_colnames = TRUE,
         main = "Heatmap of all Genes",
         annotation_col = annotation_df) # Add a title
## If you want to make heatmap for all the genes after filtering steps, use logcpm instead of expression matrix
#pheatmap(expression_matrix,
#         color = color_palette,
#       cluster_rows = TRUE,    
#         cluster_cols = TRUE,  
#         show_rownames = TRUE, 
#        show_colnames = TRUE,
#         main = "Heatmap of Significant Genes",
#         annotation_col = annotation_df) # Add a title

#Heatmap for sample to sample distances

expression_matrix2 <- t(logcpm[rownames(resSig), ])
# Calculate sample distances
sampleDists <- dist(expression_matrix2, method = "euclidean") 
sampleDistMatrix <- as.matrix(sampleDists)

# Color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Heatmap of Sample-to-Sample Distances") 

df <- data.frame(coldata1$Treatment, logcpm[topVarGenes,])

# Perform PCA using limma's plotMDS function
plotMDS(expression_matrix, main = "PCA Plot", col = as.numeric(coldata1$Treatment)) 

#GSEA rank file making

DGE_Anno_Rank <-  within(std_treated.topTags$table, rank <- sign(logFC) * -log10(PValue))
DGE_Anno_Rank 

#subset the results so only Gene Name and rank
DGErank = subset(DGE_Anno_Rank, select = c(genes,rank) )
DGErank

#sebset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)



### IF NEEDED....remove the "gene-" from row names
## https://stackoverflow.com/questions/39897155/replace-and-remove-part-of-string-in-rownames/39897315  "URS000075AF9C-snoRNA_GTATGTGTGGACAGCACTGAGACTGAGTCT"    to   "snoRNA"
## We can use gsub to match one of more characters that are not a - ([^-]+) from the start (^) of the string followed by 
## a - or (|) one or more characters that are not an underscore ([^_]+) until the end of the string ($)  and replace it with blanks ("").
## gsub("^[^-]+-|_[^_]+$", "", rownames(df))

#gene <-gsub("^[^-]+-", "", rownames(DGErank))
#DGErankIDs  <-cbind(gene,DGErank)
#head(DGErankIDs)
#summary(DGErankIDs)  

#write.csv(as.data.frame(DGErankIDs), file="DGErankIDs.csv", row.names=FALSE)
#write.table(as.data.frame(DErank_withName),file="DEGrank_withName.txt",quote=FALSE, row.names=FALSE,sep="\t")
write.table(as.data.frame(DGErank_withName),file="BirdDGErankedge_withName1.rnk",quote=FALSE, row.names=FALSE,sep="\t")

# Create the data frame
NormTransExp <- data.frame(DGErank$genes, logcpm) 

# Export as CSV
write.csv(NormTransExp, file="BirdNormTransExpressionData.txt", row.names=FALSE)


#2. Reading second dataset 3hr control and 2-4wk Htreatment the feature count file as "counts.df"
countdata2 <- read.csv("3C_2_4W_Birdgene_count_matrix.csv", row.names="gene_id")
# Clean the gene IDs (currently your row names)
clean_gene_ids2 <- str_replace(rownames(countdata2), "gene-", "")
clean_gene_ids2 <- str_remove(clean_gene_ids2, "\\|.+") 

# Update the row names with the cleaned versions
rownames(countdata2) <- clean_gene_ids2


# Printing the start of the counts.df object in R...
dim(countdata2)
head(countdata2)
head(clean_gene_ids2)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata2 <-(read.table(file.choose(), header=TRUE,row.names = 1))
dim(coldata2)
head(coldata2)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata2) %in% colnames(countdata2))
countdata2 <- countdata2[, rownames(coldata2)]
all(rownames(coldata2) == colnames(countdata2))

# Subsetting gene counts according to experimental condition
counts_control.df1  <- countdata2[,c("ERR2984230", "ERR2984231", "ERR2984232")]
counts_treated.df1 <- countdata2[,c("ERR2984239", "ERR2984238","ERR2984249")]

# Printing the structure of the gene counts set and subsets
str(countdata2)
str(counts_control.df1)
str(counts_treated.df1)

#Relative standard deviation (RSD), in edgeR, is calculated as the standard
#deviation of a gene's expression values, across samples, divided by the gene's 
#mean expression. We use RSD filtering to identify and remove genes with highly 
#variable expression across samples.  This filtering is important because genes 
#with high RSD may be unduly influenced by technical noise or sample-specific outliers, 
#making it harder to identify true biological differences. 
#By filtering out these genes, we improve the reliability and statistical power of our differential expression analysis.

# Defining function "RSD.test()"
RSD.test <- function(dataframe){
  # This function tests whether the relative standard deviation (RSD) is less
  # than or equal to one for each row in a data frame.
  # It adds the result to a new variable in the data frame called "RSD.test".
  # For a given row, if data.frame$RSD.test is TRUE, that row has an RSD less
  # than or equal to one, i.e. RSD <= 1.
  # If data.frame$RSD.test is FALSE, that row has an RSD outside of this range.
  RSD_tests = dataframe[,1]
  for (row_index in 1:nrow(dataframe)){
    row = as.numeric(dataframe[row_index,])
    RSD = sd(row) / mean(row)
    RSD_tests[row_index] = RSD <= 1 || is.na(RSD)
  }
  dataframe$RSD.test <- as.factor(RSD_tests)
  levels(dataframe$RSD.test) <- c(FALSE, TRUE)
  return(dataframe)
}

# Applying RSD.test() to gene count subsets
counts_control.df1  <- RSD.test(counts_control.df1)
counts_treated.df1 <- RSD.test(counts_treated.df1)

# Printing the structure of the gene counts subsets
str(counts_control.df1)
str(counts_treated.df1)

# Creating list of genes which failed RSD test
RSD_failed_genes <- rownames(counts_control.df1[
  which(counts_control.df1$RSD.test == FALSE),])
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_treated.df1[
  which(counts_treated.df1$RSD.test == FALSE),]))
RSD_failed_genes <- unique(RSD_failed_genes)
length(RSD_failed_genes)

# Filtering gene counts
filtered_counts.df1 <- countdata2[
  which(!rownames(countdata2) %in% RSD_failed_genes),]

# Printing the structure of the filtered gene counts
str(filtered_counts.df1)

# Checking that gene counts were correctly filtered
nrow(countdata2) - length(RSD_failed_genes) == nrow(filtered_counts.df1)

# Creating a DGEList object using the filtered gene counts
counts.DGEList2 <- DGEList(counts = filtered_counts.df1,
                          genes = rownames(filtered_counts.df1))

#If you dont want to use RSD test, then creating a DGEList object using the original unfiltered gene counts
#counts.DGEList <- DGEList(counts = countdata, genes = rownames(countdata))

# Printing the design table
print(coldata2)


# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df1) == rownames(coldata2))

# for unfiltered(with noRSD filtering) data use this
#summary(colnames(countdata) == rownames(coldata))

# Add grouping information to DGEList object
counts.DGEList2$samples$group <- as.factor(coldata2$type)

# Printing counts.DGEList
counts.DGEList2
dim(counts.DGEList2)
# save into a csv file 
# Assuming counts.DGEList is already a data frame
write.csv(counts.DGEList2, "3C_2_4Hcounts_DGEList.csv", row.names = FALSE)

# Creating an object to filter genes with low expression
counts.keep <- filterByExpr(counts.DGEList2)
summary(counts.keep)

# Filtering lowly expressed genes
filtered.DGEList1 <- counts.DGEList2[counts.keep, , keep.lib.sizes = FALSE]
dim(filtered.DGEList1)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in counts.keep
length(counts.keep[counts.keep == TRUE]) == dim(filtered.DGEList1)[1]

# Printing the normalisation factors for the libraries
filtered.DGEList1$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
norm.DGEList1 <- calcNormFactors(filtered.DGEList1)
norm.DGEList1$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition <- coldata$type
disp.DGEList1 <- estimateDisp(norm.DGEList1, design = model.matrix(~condition))
disp.DGEList1

# Assuming counts.DGEList is already a data frame
write.csv(disp.DGEList1, "3C_2_4Hdisp.DGEList.csv", row.names = FALSE)

# Exact tests for differences between experimental conditions
std_treated.DGEExact1 <- exactTest(disp.DGEList1, pair = c("control",
                                                         "heat-treatment")) 

# Extracting most differentially expressed genes from exact tests
std_treated.topTags1 <- topTags(std_treated.DGEExact1,n = nrow(std_treated.DGEExact1$table))

# Printing the most differentially expressed genes
std_treated.topTags1
# Assuming counts.DGEList is already a data frame
write.csv(std_treated.topTags1, "3C_2_4Hstd_treated.topTags.csv", row.names = FALSE)

# extract significant differentially expressed genes, sort, & write to csv
resOrdered1 <-std_treated.topTags$table[order(std_treated.topTags$table$FDR),]
resSig <- resOrdered1[resOrdered1$FDR<0.05,]
resSig
print(resOrdered1)

out <- resOrdered1 %>% dplyr::select(logFC, logCPM, PValue, FDR) %>% dplyr::rename(meanExpr = logCPM, pval = PValue, adj.pval = FDR)
write.csv(as.data.frame(out),file = paste0("/Users/temitopefolorunso/Desktop/Functional Genomics/Bird_Project/edgeR3.csv")) #write results to a new csv

#3. Reading second dataset 3hr control and 3Control-3Htreatment the feature count file as "counts.df"
countdata3 <- read.csv("3C_3HBirdgene_count_matrix.csv", row.names="gene_id")
# Clean the gene IDs (currently your row names)
clean_gene_ids3 <- str_replace(rownames(countdata3), "gene-", "")
clean_gene_ids3 <- str_remove(clean_gene_ids3, "\\|.+") 

# Update the row names with the cleaned versions
rownames(countdata3) <- clean_gene_ids3


# Printing the start of the counts.df object in R...
dim(countdata3)
head(countdata3)
head(clean_gene_ids3)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata3 <-(read.table(file.choose(), header=TRUE,row.names = 1))
dim(coldata3)
head(coldata3)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata3) %in% colnames(countdata3))
countdata3 <- countdata3[, rownames(coldata3)]
all(rownames(coldata3) == colnames(countdata3))

# Subsetting gene counts according to experimental condition
counts_control.df3  <- countdata3[,c("ERR2984230", "ERR2984231", "ERR2984232","ERR2984233")]
counts_treated.df3 <- countdata3[,c("ERR2984242", "ERR2984243","ERR2984244","ERR2984245")]

# Printing the structure of the gene counts set and subsets
str(countdata3)
str(counts_control.df3)
str(counts_treated.df3)

#Relative standard deviation (RSD), in edgeR, is calculated as the standard
#deviation of a gene's expression values, across samples, divided by the gene's 
#mean expression. We use RSD filtering to identify and remove genes with highly 
#variable expression across samples.  This filtering is important because genes 
#with high RSD may be unduly influenced by technical noise or sample-specific outliers, 
#making it harder to identify true biological differences. 
#By filtering out these genes, we improve the reliability and statistical power of our differential expression analysis.

# Defining function "RSD.test()"
RSD.test <- function(dataframe){
  # This function tests whether the relative standard deviation (RSD) is less
  # than or equal to one for each row in a data frame.
  # It adds the result to a new variable in the data frame called "RSD.test".
  # For a given row, if data.frame$RSD.test is TRUE, that row has an RSD less
  # than or equal to one, i.e. RSD <= 1.
  # If data.frame$RSD.test is FALSE, that row has an RSD outside of this range.
  RSD_tests = dataframe[,1]
  for (row_index in 1:nrow(dataframe)){
    row = as.numeric(dataframe[row_index,])
    RSD = sd(row) / mean(row)
    RSD_tests[row_index] = RSD <= 1 || is.na(RSD)
  }
  dataframe$RSD.test <- as.factor(RSD_tests)
  levels(dataframe$RSD.test) <- c(FALSE, TRUE)
  return(dataframe)
}

# Applying RSD.test() to gene count subsets
counts_control.df3  <- RSD.test(counts_control.df3)
counts_treated.df3 <- RSD.test(counts_treated.df3)

# Printing the structure of the gene counts subsets
str(counts_control.df3)
str(counts_treated.df3)

# Creating list of genes which failed RSD test
RSD_failed_genes <- rownames(counts_control.df1[
  which(counts_control.df3$RSD.test == FALSE),])
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_treated.df3[
  which(counts_treated.df3$RSD.test == FALSE),]))
RSD_failed_genes <- unique(RSD_failed_genes)
length(RSD_failed_genes)

# Filtering gene counts
filtered_counts.df3 <- countdata3[
  which(!rownames(countdata3) %in% RSD_failed_genes),]

# Printing the structure of the filtered gene counts
str(filtered_counts.df3)

# Checking that gene counts were correctly filtered
nrow(countdata3) - length(RSD_failed_genes) == nrow(filtered_counts.df3)

# Creating a DGEList object using the filtered gene counts
counts.DGEList3 <- DGEList(counts = filtered_counts.df3,
                           genes = rownames(filtered_counts.df3))

#If you dont want to use RSD test, then creating a DGEList object using the original unfiltered gene counts
#counts.DGEList <- DGEList(counts = countdata, genes = rownames(countdata))

# Printing the design table
print(coldata3)


# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df3) == rownames(coldata3))

# for unfiltered(with noRSD filtering) data use this
#summary(colnames(countdata) == rownames(coldata))

# Add grouping information to DGEList object
counts.DGEList3$samples$group <- as.factor(coldata3$type)

# Printing counts.DGEList
counts.DGEList3
dim(counts.DGEList3)
# save into a csv file 
# Assuming counts.DGEList is already a data frame
write.csv(counts.DGEList3, "3C_3Hcounts_DGEList.csv", row.names = FALSE)

# Creating an object to filter genes with low expression
counts.keep1 <- filterByExpr(counts.DGEList3)
summary(counts.keep1)

# Filtering lowly expressed genes
filtered.DGEList3 <- counts.DGEList3[counts.keep1, , keep.lib.sizes = FALSE]
dim(filtered.DGEList3)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in counts.keep
length(counts.keep1[counts.keep == TRUE]) == dim(filtered.DGEList3)[1]

# Printing the normalisation factors for the libraries
filtered.DGEList3$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
norm.DGEList3 <- calcNormFactors(filtered.DGEList3)
norm.DGEList3$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition <- coldata3$type
disp.DGEList3 <- estimateDisp(norm.DGEList3, design = model.matrix(~condition))
disp.DGEList3

# Assuming counts.DGEList is already a data frame
write.csv(disp.DGEList3, "3C_3Hdisp.DGEList.csv", row.names = FALSE)

# Exact tests for differences between experimental conditions
std_treated.DGEExact3 <- exactTest(disp.DGEList3, pair = c("control",
                                                           "heat-treatment")) 

# Extracting most differentially expressed genes from exact tests
std_treated.topTags3 <- topTags(std_treated.DGEExact3,n = nrow(std_treated.DGEExact3$table))

# Printing the most differentially expressed genes
std_treated.topTags3
# Assuming counts.DGEList is already a data frame
write.csv(std_treated.topTags3, "3C_3CHstd_treated.topTags.csv", row.names = FALSE)

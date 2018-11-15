# ---
# title: "Shakeel's PCA script"
# author: Shakeel Jessa
# date: 18/Jun/2018
# output: pdf_document
# ---

# Import packages
library('ggplot2')
library(dplyr)
library(gplots)
library(getopt)

spec = matrix(c(
  'gene', 'g', 1, "character",
  'annotation', 'a', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
cat(opt$gene)

# Collect file names from user
#input_counts <- readline(prompt = "What is the name of your gene counts file?: ")
#input_annot <- readline(prompt = "What is the name of your annotation file?: ")
input_counts <- opt$gene
input_annot <- opt$annotation
# Read a txt file and convert gene names (first column) into indeces
my_data <- read.delim(input_counts, header = TRUE, row.names = 1, check.names = FALSE)

# Read a text file and create indeces
annot <- read.delim(input_annot, header = TRUE, row.names = 1, check.names = FALSE)

# Transpose data table
gene_df <- t(my_data)

# Convert transposed table into a data frame
gene_df <- data.frame(gene_df)

# Remove all non-expressed and 0-variance genes from data frame
express_df <- gene_df[,sapply(gene_df,function(x){!all(x==0)})]
express_final_df <- express_df[vapply(express_df, function(x) length(unique(x)) > 1, logical(1L))]

# Scale data frame
scale_df <- scale(express_final_df)

# Perform PCA analysis on scaled data and store
pca <- prcomp(scale_df)

# Plot PC1 vs PC2 and convert pca data into a data frame
pca_df <- data.frame(pca$x)

# Count number of columns in pca_df and store it
pre_ncol <- ncol(pca_df)

# Add annotations to pca dataframe
pca_df <- data.frame(pca$x, annot)

# Count number of columns in pca_df after adding annot and store
post_ncol <- ncol(pca_df)

# Generate range of values to iterate through for plot labeling
label_ncol <- pre_ncol + 1

# Generate Importance of Components and store it
components <- summary(pca)

# Pull values of the two highest Proportion of Variance's
Prop_ofPC1 <- components$importance[2, 1] * 100
Prop_ofPC2 <- components$importance[2, 2] * 100

# Convert Variance values into strings
Prop_ofPC1 <- toString(Prop_ofPC1)
Prop_ofPC2 <- toString(Prop_ofPC2)

# Create x and y axis labels with variance %
x_label <- paste("PC1,", Prop_ofPC1, sep=' Var %: ')
y_label <- paste("PC2,", Prop_ofPC2, sep=' Var %: ')

# Capture the rotation matrix in a data frame
rotation_data <- data.frame(pca$rotation, variable=row.names(pca$rotation))

# Sort PC1/PC2 columns from G->L and store in new data frames
sorted_PC1 <- rotation_data[order(rotation_data$PC1, decreasing=TRUE),, drop=FALSE]
sorted_PC2 <- rotation_data[order(rotation_data$PC2, decreasing=TRUE),, drop=FALSE]

# Store gene names and values for top 10 genes
PC1_names <- rownames(sorted_PC1[1:10,])
PC1_values <- sorted_PC1[1:10,1]
PC1_10list <- paste(PC1_names, PC1_values, sep=': ')
PC2_names <- rownames(sorted_PC2[1:10,])
PC2_values <- sorted_PC2[1:10,2]
PC2_10list <- paste(PC2_names, PC2_values, sep=': ')

#Output PCA plot(s) to pdf
pdf("pca_plots.pdf")
for (n in c(label_ncol:post_ncol)) {print(ggplot(pca_df, aes(x=PC1, y=PC2, color= pca_df[, n]))
                                          + geom_point(size = 2) + labs(x =x_label, y = y_label, color=colnames(pca_df[n])))}
textplot(PC1_10list, cex = 1.65)
title("10 most highly contributing genes to PC1")
textplot(PC2_10list, cex= 1.65)
title("10 most highly contributing genes to PC2")
dev.off()

#write ("The number of genes removed during filtering is: ", stderr() )
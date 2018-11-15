#install.packages("Rtsne") # Install Rtsne package from CRAN
library("Rtsne")

input_counts = "genecounts.variancestabilized.csv"
input_annot = "meta_data_with_paths_all_patients_final.csv"

# let's get the file

# Read a txt file and convert gene names (first column) into indeces
my_data <- read.csv(input_counts, header = TRUE, row.names = 1, check.names = FALSE)

# Read a text file and create indeces
annot <- read.csv(input_annot, header = TRUE, row.names = 1, check.names = FALSE )

# Transpose data table
gene_df <- t(my_data)

# Convert transposed table into a data frame
gene_df <- data.frame(gene_df)

# remove stuff we don't want
gene_df <- gene_df[rownames(annot),]
gene_df <- gene_df[,sapply(gene_df,function(x){!all(x==0)})]
gene_df <- gene_df[vapply(express_df, function(x) length(unique(x)) > 1, logical(1L))]

gene_df <- gene_df[complete.cases(gene_df), ] # remove rows with NaN, I think there is 1 like that
head(gene_df)

png("t-snePlot perplexity 10.png")
labels <- rownames(gene_df)
colors = rainbow(length(unique(labels)))
names(colors) = unique(labels)
tsne_model = Rtsne(gene_df[2,], perplexity=10, verbose=TRUE, check_duplicates=FALSE) # be careful with downstream

plot.new()
plot(tsne_model$Y,  t='n',main='Tsne')
text(tsne_model$Y, labels=labels, col=colors[labels])

dev.off()
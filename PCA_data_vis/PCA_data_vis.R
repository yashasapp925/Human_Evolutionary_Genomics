# Install packages
install.packages("umap")
install.packages ("ggplot2")
install.packages ("magrittr")
install.packages ("dplyr")

# Load packages
library(umap)
library(ggplot2)
library(magrittr)
library(dplyr)

# Load data
Global_data<-read.table("Global_50kSNPs_64pops.txt",sep="\t",header=TRUE)

# View first 10 rows of data
Global_data[1:10,]

# Determine number of SNPs
nrow(Global_data)

# Extract allele frequency data
freqs=Global_data [,2:65]

# Transpose the matrix of allele frequencies
freqs_t=t(freqs)

# Use prcomp function to run PCA
my_pca=prcomp(freqs_t)

# Find percent variance explained by each principal component
summary(my_pca)

# Plot PC1 and PC2
plot(my_pca$x[,1],my_pca$x[,2])

# Obtain the list of the population names
PopNames=names(Global_data [,2:65])

# Plot PC1 and PC2 using population names
plot(my_pca$x[,1],my_pca$x[,2],cex=0,xlab="PC1",ylab="PC2")
text(my_pca$x[,1],my_pca$x[,2],PopNames,cex=0.5)

# Save plot
pdf("PCA_data_vis/PCA.pdf")
plot(my_pca$x[,1],my_pca$x[,2],cex=0,xlab="PC1",ylab="PC2")
text(my_pca$x[,1],my_pca$x[,2],PopNames,cex=0.5)
dev.off()

# Plot a different pair of PCs
pdf("PCA_data_vis/PCA_diffpair.pdf")
plot(my_pca$x[,3],my_pca$x[,4],cex=0,xlab="PC1",ylab="PC2")
text(my_pca$x[,3],my_pca$x[,4],PopNames,cex=0.5)
dev.off()

# Repeat above analysis using only 100 SNPs
LessSNPs= Global_data[sample(nrow(Global_data),100),]
freqs_less=LessSNPs[,2:65]
freqs_less_t=t(freqs_less)
my_pca_less=prcomp(freqs_less_t)
pdf("PCA_data_vis/PCA_100SNPs.pdf")
plot(my_pca_less$x[,1],my_pca_less$x[,2],cex=0,xlab="PC1",ylab="PC2")
text(my_pca_less$x[,1],my_pca_less$x[,2],PopNames,cex=0.5)
dev.off()

# Compute UMAP embedding
umap_result <- umap(freqs_t, n_neighbors=10)

# Read population codes file
popcodes <- read.table("1M1KG_population_codes.txt", sep="\t", header=TRUE)

# Convert UMAP result to data frame and rename columns
umap_df <- as.data.frame(umap_result$layout) %>% rename(UMAP1 = V1, UMAP2 = V2)

# Add Population column as characters
umap_df$Population <- as.character(popcodes$Population)

# Join with popcodes using the Population column
umap_df <- inner_join(umap_df, popcodes, by = "Population")

# Plot UMAP plot with ggplot2
umap_plot <- umap_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Superpop, shape = Subsistence_pattern)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2", subtitle = "UMAP plot")
ggsave("PCA_data_vis/UMAP_plot_points.png", plot = umap_plot, width = 10, height = 8, units = "in")

umap_plot <- umap_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Superpop, shape = Subsistence_pattern, label=Population)) +
  geom_text() +
  labs(x = "UMAP1", y = "UMAP2", subtitle = "UMAP plot")
ggsave("PCA_data_vis/UMAP_plot_text.png", plot = umap_plot, width = 10, height = 8, units = "in")



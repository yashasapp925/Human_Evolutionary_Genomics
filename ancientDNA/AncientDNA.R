# Install and load packages
install.packages("mapview")
install.packages ("ggplot2")
install.packages ("magrittr")
install.packages ("dplyr")

library(mapview)
library(ggplot2)
library(magrittr)
library(dplyr)

# Load data
data<-read.table("aDNA_MetaData_EffectAlleleCalls.txt",sep="\t",header=TRUE)
data[1:10,1:20]

# Find the number of rows and columns
nrow(data)
ncol(data)

# Use the mapview package to generate an interactive map that shows the geographic location of each ancient sample
mapview(data, xcol="Longitude", ycol="Latitude",crs=4269, grid=FALSE)


# Examine the age distribution of the ancient samples
pdf("ancientDNA/age_dist.pdf")
hist(data[,'SampleDate_YearsBP'], breaks=60, xlim=c(0,60000), main="Sample age distribution", xlab="Years before present", ylab="Counts")
dev.off()

# Examine how missingness varies for different amounts of sequence coverage
pdf("ancientDNA/missingness.pdf")
plot(data[,'Coverage'], data[,'Missingness'], main="Genotyping quality", xlab="Sequence coverage", ylab="Missingness")
dev.off()

# Find the rs# of a random SNP and generate a loess smoothed plot of allele frequency vs time, 
# focusing only on samples from the last 10,000 years
rs_num <- colnames(data)[1000]
ggplot(data=data, aes(x=SampleDate_YearsBP, y= data[,rs_num])) + 
                                                    geom_point() + geom_smooth(method = "loess", se = FALSE) + xlim(0,10000) + ylim(0,1)

# Split data into two parts
sample_info=data[,1:9]
calls= data[,10:3617]

# Impute missing data by replacing missing genotype calls with the mean allele frequency in the dataset
for(i in c(1:3608)){
  column_i <- calls[,i]
  mean_i <- (mean(column_i, na.rm=TRUE))
  NAs_i <- which(is.na(column_i))
  column_i[NAs_i] <- mean_i
  calls[,i] <- column_i
}

# Use prcomp() function to run PCA
aDNA_pca=prcomp(calls)
summary(aDNA_pca)

# Convert PCs 1 and 2 to a data frame
pca_df <- as.data.frame(aDNA_pca$x[,1:2]) 

# Add SampleID column as characters
pca_df$SampleID <- as.character(sample_info[,'SampleID'])

# Join with sample_info using the SampleID column
pca_df <- inner_join(pca_df, sample_info, by = "SampleID")
pca_df[1:10,]

# Plot first two PCs using ggplot2 with different colors indicating different longitude coordinates
# Each point corresponds to a different ancient sample
pca_plot <- pca_df %>%
  ggplot(aes(x = PC1, y = PC2, color = Longitude)) + geom_point() + labs(x = "PC1", y = "PC2", subtitle = "PCA plot")
ggsave("ancientDNA/PCA_plot_Longitude_Points.png", plot = pca_plot, width = 10, height = 8, units = "in")

pca_plot <- pca_df %>%
  ggplot(aes(x = PC1, y = PC2, color = Longitude, label=SampleID)) +
  geom_text() + labs(x = "PC1", y = "PC2", title = "PCA plot")
ggsave("ancientDNA/PCA_plot_Longitude_SampleID.png", plot = pca_plot, width = 10, height = 8, units = "in")

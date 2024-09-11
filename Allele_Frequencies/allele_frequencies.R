# Load data
data=read.table("chr22_phase3_DAF.txt",header=TRUE,sep="\t")

# Select a random 10000 SNPs out of the dataset
subset=data[sample(nrow(data),10000),]

# View the first five rows of this subset of the data
subset[1:5,]

# Generate a histogram of allele frequencies in Europe
pdf("Allele_Frequencies/hist_EUR.pdf")
hist(subset[,'EUR_DAF'],breaks=100)
dev.off()

# Zoom in on low frequency alleles
pdf("Allele_Frequencies/hist_EUR_subset.pdf")
hist(subset[,'EUR_DAF'],breaks=1000,xlim=c(0,0.05))
dev.off()

# Generate a scatterplot of African allele frequencies vs. East Asian allele frequencies
pdf("Allele_Frequencies/scatter_AFR_EAS.pdf")
plot(subset[,'AFR_DAF'],subset[,'EAS_DAF'],cex=0.25)
dev.off()

# Convert derived allele frequencies (DAF) to minor allele frequencies (MAF)
SAS_MAF=pmin(data[,'SAS_DAF'],1-data[,'SAS_DAF'])

# Generate a folded frequency spectrum
pdf("Allele_Frequencies/fold_freq_spec.pdf")
hist(SAS_MAF,breaks=100)
dev.off()

# Generate a folded spectrum for only polymorphic sites
pdf("Allele_Frequencies/fold_freq_spec_polymorphic.pdf")
hist((SAS_MAF[SAS_MAF>0]),breaks=100)
dev.off()

# Plot the DAF of all the variants on a chromosome 22
pdf("Allele_Frequencies/DAF_all_var_chr22.pdf")
plot(data[,'POS_HG19'],data[,'AMR_DAF'], cex=0.25)
dev.off()

# Plot the MAF of all the variants on a chromosome 22
pdf("Allele_Frequencies/MAF_all_var_chr22.pdf")
plot(data[,'POS_HG19'],SAS_MAF, cex=0.25)
dev.off()

# Plot the absolute value of the difference in allele frequency between two populations across all of chromosome 22
pdf("Allele_Frequencies/absval_AFdif.pdf")
plot(data[,'POS_HG19'],abs(data[,'EAS_DAF']-data[,'EUR_DAF']), cex=0.25)
dev.off()

# Plot the absolute value of the difference in allele frequency between two populations across a subset of chromosome 22
pdf("Allele_Frequencies/absval_AFdif_subset.pdf")
plot(data[,'POS_HG19'],abs(data[,'EAS_DAF']-data[,'EUR_DAF']), xlim=c(35000000, 36000000),cex=0.25)
dev.off()
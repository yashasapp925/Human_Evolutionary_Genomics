# Load data
data=read.table("AlleleFreqs_1KG_chrom16.txt",header=TRUE)
data[1:5,]

# Select three populations and obtain the allele frequencies for each
p1=data[,'VarFreq_CEU']
p2=data[,'VarFreq_YRI']
p3=data[,'VarFreq_JPT']

# Calculate the mean and variance of allele frequencies
mean12=(p1+p2)/2
mean13=(p1+p3)/2
mean23=(p2+p3)/2
var12=0.5*(p1-mean12)^2+0.5*(p2-mean12)^2
var13=0.5*(p1-mean13)^2+0.5*(p3-mean13)^2
var23=0.5*(p2-mean23)^2+0.5*(p3-mean23)^2

# Calculate FST statistics for each SNP between each pair of populations
Fst12=var12/(mean12*(1-mean12))
Fst13=var13/(mean13*(1-mean13))
Fst23=var23/(mean23*(1-mean23))

# Force any undefined values to equal zero
Fst12[is.nan(Fst12)] <- 0
Fst13[is.nan(Fst13)] <- 0
Fst23[is.nan(Fst23)] <- 0

# Rescale genetic distances between populations 
T12=-log(1-Fst12)
T13=-log(1-Fst13)
T23=-log(1-Fst23)
PBS1=(T12+T13-T23)/2

# Ensure no negative PBS scores
PBS1[PBS1<0] <--0

# Merge PBS score calculations with rs# information and genomic position
rs=data[,'ID']
chr= data[,'CHROM']
pos= data[,'POS']
PBS_unsorted=cbind(as.data.frame(rs),chr,pos,PBS1)

# Sort scores
PBS_sorted=PBS_unsorted[order(PBS_unsorted[,'PBS1'],decreasing=TRUE),]
PBS_sorted[1:10,]

# Plot PBS scores vs. chromosome position
options(scipen=999)
plot(PBS_sorted[,'pos'],PBS_sorted[,'PBS1'],xlab="Genomic position",ylab="PBS",cex=0.25)

pdf("PBS_stats/PBSvsChrPos.pdf")
plot(PBS_sorted[,'pos'],PBS_sorted[,'PBS1'],xlab="Genomic position",ylab="PBS",cex=0.25,xlim=c(20000000,22000000))
dev.off()



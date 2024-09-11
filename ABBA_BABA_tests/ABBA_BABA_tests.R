# Install packages
install.packages("bootstrap")

# Load packages
library(bootstrap)

# Function to calculate Patterson’s D across an entire chromosome
D_chromosome_wide <- function(data, pop_1, pop_2, pop_3, pop_4) {
  cdata = data[, c("POS", pop_1, pop_2, pop_3, pop_4)]
  cdata = na.omit(cdata)
  abba = sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 1 & cdata[[5]] == 0) + 
    sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 0 & cdata[[5]] == 1)
  baba = sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 1 & cdata[[5]] == 0) +
    sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 0 & cdata[[5]] == 1)
  D = (abba - baba) / (abba + baba)
  return(D)
}

# Function to calculate Patterson’s D for sliding windows
D_sliding_window <- function(data, pop_1, pop_2, pop_3, pop_4, window_size) {
  n_windows = max(data$POS) %/% window_size + 1
  abba_windows = rep(0, n_windows)
  baba_windows = rep(0, n_windows)
  data = data[, c("POS", pop_1, pop_2, pop_3, pop_4)]
  window_start = 1
  for (i in 1:n_windows) { 
    cdata = data[data$POS >= window_start & data$POS < window_start + window_size, 
                 c("POS", pop_1, pop_2, pop_3, pop_4)]
    cdata = na.omit(cdata) # drop na
    abba_windows[i] = sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 1 & cdata[[5]] == 0) + 
      sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 0 & cdata[[5]] == 1)
    baba_windows[i] = sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 1 & cdata[[5]] == 0) +
      sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 0 & cdata[[5]] == 1)
    window_start = window_start + window_size
  }
  print(sum(abba_windows))
  print(sum(baba_windows))
  return(list(abba_windows, baba_windows,
              ((abba_windows + 1) - (baba_windows + 1)) / ((abba_windows + 1) + (baba_windows + 1))))
}

# Load data
data = read.table('chr7_1KG_altai_vindija_chimp_final.tab', sep='\t', header=T)

# Define 4 populations
pop_1 = 'YRI'
pop_2 = 'CEU'
pop_3 = 'altai' 
pop_4 = 'chimp' # Outgroup

# Calculate counts of ABBA SNPs, counts of BABA SNPs, Patterson’s D statistic, and a jack-knifed p-value
cdata = data[, c("POS", pop_1, pop_2, pop_3, pop_4)]
cdata = na.omit(cdata)
abba = sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 1 & cdata[[5]] == 0) + 
  sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 0 & cdata[[5]] == 1)
baba = sum(cdata[[2]] == 1 & cdata[[3]] == 0 & cdata[[4]] == 1 & cdata[[5]] == 0) +
  sum(cdata[[2]] == 0 & cdata[[3]] == 1 & cdata[[4]] == 0 & cdata[[5]] == 1)
D = D_chromosome_wide(data, pop_1, pop_2, pop_3, pop_4)
window_size = 1e7
windows = D_sliding_window(data, pop_1, pop_2, pop_3, pop_4, window_size)
D_windows = windows[[3]]
jackknife_results <- jackknife(D_windows, mean) # Apply jackknife (i.e., leave one window out) to obtain estimate of standard error
mean_jackknife = mean(jackknife_results$jack.values)
se_jackknife = jackknife_results$jack.se
Z_jackknife = mean_jackknife / se_jackknife
pval_jackknife = 2 * pnorm(-abs(Z_jackknife))

# View key stats
abba
baba
D
pval_jackknife

# View how Patterson’s D varies along chromosome 7
pdf("ABBA_BABA_tests/PatD_variation_chr7.pdf")
plot((seq(from=window_size, to=max(data$POS) + window_size, by=window_size) -
        window_size %/% 2) / 1e6, D_windows,
     type='b', xlab='Chromosome position in Mb',
     ylab=paste("D(", pop_1, ",", pop_2,",",  pop_3, ",", pop_4, ")", sep=''), ylim=c(-1, 1))
abline(h = 0, lty=2)
arrows(x0=(seq(from=window_size, to=max(data$POS) + window_size, by=window_size) - window_size %/% 2) / 1e6, 
       y0=D_windows - jackknife_results$jack.se,
       x1=(seq(from=window_size, to=max(data$POS) + window_size, by=window_size) - window_size %/% 2) / 1e6,
       y1=D_windows + jackknife_results$jack.se, code=3, angle=90, length=0.05, col="black", lwd=1)
abline(h = mean(jackknife_results$jack.values), lty=3)
legend(1, 1, legend=c("D=0", "E[D]"),
       col=c("black", "black"), lty=2:3, cex=0.8)
dev.off()

# Show counts of possible configuration of alleles for four populations
data_merged=paste(cdata[,2], cdata[,3], cdata[,4], cdata[,5],sep="",colapse="")
table(data_merged)
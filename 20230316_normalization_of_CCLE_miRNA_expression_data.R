# This script is to try various normalization of CCLE miRNA expression data
# 2023/03/16 made

# activate package to draw meanSdplot
library(vsn)
library(ggplot2)

# create new directory
setwd('C:/Rdata')
dir.create("20230316_normalization_of_CCLE_miRNA_expression_data")

# import CCLE miRNA expression data
# this data is located at "//fsz-p21.naist.jp/okamura-lab/Files_related_to_M1_Projects/Hirota/CCLE_Data"
setwd("C:/Rdata/CCLE_data")
CCLE.miRNA.data <-read.table("CCLE_miRNA_20181103.gct.txt",sep="\t",header = T,stringsAsFactors = F)
d <-duplicated(CCLE.miRNA.data[,2])
CCLE.miRNA.data <-CCLE.miRNA.data[!d,]

# change directory
setwd("C:/Rdata/20230316_normalization_of_CCLE_miRNA_expression_data")

# draw meansd plot of log2 normalization
log2.df <-apply(CCLE.miRNA.data[c(-1,-2)], 2, log2)
pdf("meansd_plot_of_log2_normalization.pdf")
l <-meanSdPlot(log2.df)
l$gg+scale_y_continuous(limits = c(0, 6))
dev.off()

# activate package to do median normalization
library(limma)

# do median normalization
median.df <-as.data.frame(normalizeMedianValues(CCLE.miRNA.data[,c(-1,-2)]))
median.df[,955] <-CCLE.miRNA.data[,2]
median.df <-median.df[,c(955,1:954)]
colnames(median.df)[1] <-"miRNA"
write.table(median.df,"CCLE_miRNA_median_normalization.txt",sep="\t",row.names = F,quote = F)

# draw meansd plot of median normalization
log2.median <-apply(median.df[,c(2:955)],2,log2)
pdf("meansd_plot_of_median_normalization.pdf")
m <-meanSdPlot(log2.median)
m$gg+scale_y_continuous(limits = c(0, 6))
dev.off()

# activate package to do quantile normalization
library(preprocessCore)

# do quantile normalization 
quantile.df <-as.data.frame(normalize.quantiles(as.matrix(CCLE.miRNA.data[,c(-1,-2)])))
quantile.df[,955] <-CCLE.miRNA.data[,2]
quantile.df <-quantile.df[,c(955,1:954)]
colnames(quantile.df) <-colnames(median.df)
write.table(quantile.df,"CCLE_miRNA_quantile_normalization.txt",sep="\t",row.names = F,quote = F)

# draw meansd plot of quantile normalization
log2.quant <-apply(quantile.df[,c(2:955)], 2, log2)
pdf("meansd_plot_of_quantile_normalization.pdf")
q <-meanSdPlot(log2.quant)
q$gg+scale_y_continuous(limits = c(0, 6))
dev.off()

# activate package to do Variance stabilization normalization
library(vsn)

# do variance normalization ()
variance.df <-as.data.frame(justvsn(as.matrix(CCLE.miRNA.data[,c(-1,-2)])))
variance.df[,955] <-CCLE.miRNA.data[,2]
variance.df <-variance.df[,c(955,1:954)]
colnames(variance.df) <-colnames(median.df)
write.table(variance.df,"CCLE_miRNA_variance_normalization_.txt",sep="\t",row.names = F,quote = F)

# draw meansd plot of variance stabilizing normalization
pdf("meansd_plot_of_variance_stabilizing_normalization.pdf")
v <-meanSdPlot(as.matrix(variance.df[,c(2:955)]))
v$gg+scale_y_continuous(limits = c(0, 6))  
dev.off()

# do mean normalization
mean.df <-CCLE.miRNA.data[,c(-1,-2)]
mean <-as.numeric(apply(mean.df,1,mean))
for (i in 1:ncol(mean.df)) {
 mean.df[,i] <-mean.df[,i]/mean 
}
mean.df[,955] <-CCLE.miRNA.data[,2]
mean.df <-mean.df[,c(955,1:954)]
colnames(mean.df)[1] <-"miRNA"
write.table(mean.df,"CCLE_miRNA_mean_normalization_.txt",sep="\t",row.names = F,quote = F)

# draw meansd plot of variance stabilizing normalization
pdf("meansd_plot_of_mean_normalization.pdf")
m <-meanSdPlot(as.matrix(mean.df[,c(2:955)]))
dev.off()

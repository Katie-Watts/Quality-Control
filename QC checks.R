library(GWASinspector)
library(qqman)
library(ggplot2)
library(QCGWAS)

## import the QC-configuration file 
job <- setup.inspector("GWASinspector/config.ini")

## check the created instance
## input result files that will be inspected are also displayed
job

## run the algorithm 
job <- run.inspector(job)

## check the results
## comprehensive report and result file are already saved in the output folder
result.inspector(job)

# Load back in for final checks

data_GWAS <-load_GWAS("results.txt.gz")

# Check P-values match other statistics

check_P(data_GWAS,
        plot_correlation = TRUE, plot_if_threshold = FALSE,
        save_name = "STATS_after")

# To calculate a correlation between predicted and actual p-values and plot the correlation:

calc_kurtosis(data_GWAS$EFFECT)
calc_kurtosis(data_GWAS$EFF_ALL_FREQ)

calc_skewness(data_GWAS$EFFECT)
calc_skewness(data_GWAS$EFF_ALL_FREQ)


QC_plots(data_GWAS,
         plot_QQ = TRUE, plot_Man = TRUE,
         filter_NA = TRUE,
         plot_cutoff_p = 1, plot_names = FALSE,
         QQ_colors = c("red", "blue", "orange", "green3", "yellow"),
         plot_QQ_bands = TRUE,
         save_name = "dataset", save_dir = getwd(),
         use_log = FALSE,
         check_impstatus = FALSE, ignore_impstatus = TRUE,
         NA_strings = c(NA, "NA", ".", "-"))


Cairo::Cairo(
  80, #length
  30, #width
  file = paste("Manhattan", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
par(mar=c(5,6,4,1)+.1)
manhattan(data_GWAS,
          bp = "POSITION", 
          snp = 'MARKER',
          p = "PVALUE", 
          ylim = c(0, 10), 
          ylab = expression(-log[10]*"("*italic("P")*")"), 
          cex.lab=1.7,
          cex.axis=1.7,
          col = c("navyblue", "dodgerblue2"))
dev.off()

Cairo::Cairo(
  50, #length
  50, #width
  file = paste("QQ", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
par(mar=c(5,6,4,1)+.1)
qq(data_GWAS$PVALUE)
dev.off()


Cairo::Cairo(
  40, #length
  20, #width
  file = paste("Imputation", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
par(mar=c(5,6,4,1)+.1)
plot <- hist(data_GWAS$IMP_QUALITY, col = "grey", xlab="Imputation quality score", main = "Distribution of imputation quality scores") 
dev.off()


Cairo::Cairo(
  40, #length
  20, #width
  file = paste("Effect size", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
par(mar=c(5,6,4,1)+.1)
plot <- hist(data_GWAS$EFFECT, col = "grey", xlab="Effect size", main = "Distribution of Effect sizes")
dev.off()


Cairo::Cairo(
  40, #length
  20, #width
  file = paste("Effect size", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
par(mar=c(5,6,4,1)+.1)
plot <- hist(data_GWAS$PVALUE, col = "grey", xlab="P-value", main = "Distribution of P-values")
dev.off()


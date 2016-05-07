# Code to plot coverage from the output of processed dcov files.
# Input must have 3 columns and header Gene\tnumMissing\tfracMissing

# Example Useage:
# Rscript --vanilla plot_coverage.R coverage/file.coverage_metrics.txt coverage_min bq_min mq_min

# Outputs a pdf with 2 plots

library(ggplot2)
require(gridExtra)


args = commandArgs(trailingOnly=TRUE)
file = args[1]
coverage = args[2]
baseQ = args[3]
mapQ = args[4]

data<- read.table(file, header=T) 
  plot.num = ggplot(data, aes(Gene, numMissing)) + 
          geom_point(size=5, colour=c("darkblue")) + theme_bw(20) +
          theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("Bases Not Covered\n") + xlab("")  +
          ggtitle(paste("Coverage Threshold: ", coverage, "\nMin Base Quality: ", baseQ, "\nMin Mapping Quality: ", mapQ))
  


  data$percMissed = (100 * data$fracMissing)
  plot.percent = ggplot(data, aes(Gene, percMissed)) +geom_point(size=5, colour=c("darkblue")) + theme_bw(20) + 
          theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ylab("% of Bases Not Covered\n")  +
          ggtitle("")

pdf(paste0("coverageOfKeyGenes.",coverage,"cov.",baseQ,"baseQ.", mapQ, "mapQ.pdf"), height=14, round(length(unique(data$Gene))/4))
grid.arrange(plot.num, plot.percent, ncol=1)
dev.off()

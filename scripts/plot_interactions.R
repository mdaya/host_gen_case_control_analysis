source("plot_functions.R")
#TODO: check that the packages that are loaded at the top of plot_functions.R is installed

tb.data <- read.delim("../data/input/TB_assoc.txt")
genenames <- read.delim("../data/input/gene_info.txt")

pdf("../data/output/genotype_combinations.pdf", width=5.08, height=7, onefile = TRUE)
grid.arrange(
  GetGenoCombinationPlot("rs4833095", "TLR2_GT", tb.data),
  GetGenoCombinationPlot("rs4833095", "rs3804100", tb.data) )
dev.off()

pdf("../data/output/allele_combinations.pdf", width=5.08, height=9, onefile = TRUE)
grid.arrange(
  GetAlleleCombinationPlot("rs4833095", "TLR2_GT", tb.data, T, lev=c("G-T", "G-A", "A-T", "A-A")),
  GetAlleleCombinationPlot("rs4833095", "rs3804100", tb.data, F, lev=c("G-T", "G-C", "A-T", "A-C")) )
dev.off()

pdf("../data/output/logits.pdf", width=5.08, height=9, onefile = TRUE)
grid.arrange(
  GetLogitPlot("rs4833095", "TLR2_GT", tb.data),
  GetLogitPlot("rs4833095", "rs3804100", tb.data) )
dev.off()

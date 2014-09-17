infile.name <- "../data/input/TB_assoc.txt"
outfile.name <- "../data/output/haplo_summary.txt"

source("haplo_stats.R")

#TODO: Use haploview to visualize potential haplotype blocks in your gene.
#If a haplotype is present, specify the names of the SNPs in the haplotype,
#according to their physical order.
snp.names <- c("TLR4_rs1927914", "TLR4_rs10759932", "TLR4_rs2770148")
RunHaplo(infile.name, outfile.name, snp.names)

#TODO: At the top of the data/output/haplo_summary.txt, the haplotype
#frequencies are printed. Paste this into the HaploAnalysis 
#spreadsheet. Ensure that the "haplo.base" genotype is the first
#row in your table. At the bottom of the file, the unadjusted and adjusted
#P-values of the model are printed. Update your spreadsheet to report 
#these values. If the adjusted P value is significant, (<0.05), then 
#you need to look at the Coefficients part of the output and work out
#ORs and CIs for the haplotypes that are significant.

#TODO: update the effect size according to the coeff column
effect.size <-   0.381787
#TODO: update the standard error according to the se column
se <- 0.180260

#TODO: run the below to get the confidence interval for the OR
or <- exp(effect.size)
ci.lower <- exp(effect.size-1.96*se)
ci.upper <- exp(effect.size+1.96*se)
print(or)
print(ci.lower)
print(ci.upper)

#TODO: update the HaploAnalysis spreadsheet with the OR and confidence interval
#See the following page for advice on how to interpret the odds ratio:
#https://tethys.mb.sun.ac.za/share/page/site/hostgen/wiki-page?title=How_to_run_and_interpret_a_logistic_regression_model

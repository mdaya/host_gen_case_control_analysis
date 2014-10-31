# Initialize variables
#######################
tb.assoc <- read.delim("../data/input/TB_assoc.txt")
genenames <- read.delim("../data/input/gene_info.txt")
x.chr.genes <- c()
#x.chr.genes <- c("TLR8")   #TODO: If any of the genes are on the X chromosme, and you want the analysis to be stratified by gender, use this variable to list the genes on the X-chromosome
outfile.name <- "../data/output/interaction_summary.txt"
snp.start.col <- 10 #Column number of first SNP
p.vals <- c()

# Function to fit and report model
##################################
FitInteractionModel <- function(model.name, test.data, nr.tests, adjust.gender=T) {
  run.name <- model.name
  if (dim(test.data)[1] > 0) {
    nr.cases <- sum(test.data$Class == 1)
    nr.controls <- sum(test.data$Class == 0)
    snp1 <- test.data[,11]
    snp2 <- test.data[,12]
    if (adjust.gender) {
      model.out <- try(glm(Class ~ Age + Gender + SanAncestry  + AfricanAncestry  + 
                             EuropeanAncestry+ SouthAsianAncestry + snp1*snp2, 
                           family="binomial", data=test.data))
    } else {
      model.out <- try(glm(Class ~ Age  + SanAncestry  + AfricanAncestry  + 
                             EuropeanAncestry+ SouthAsianAncestry + snp1*snp2,
                           family="binomial", data=test.data))
    }
    if(inherits(model.out,"try-error")) {
      cat("ERROR! Failed to run model",
          model.name, "\n", file=paste("../data/output/interaction_errors_", run.name, ".txt", sep=""), append=T)
    } else {
      anova.out <- anova(model.out, test="Chisq")
      p.val <- anova.out[dim(anova.out)[1],dim(anova.out)[2]]
      if (!is.na(p.val) & (p.val <= 0.05)) {
        cat(model.name, nr.cases, nr.controls, p.val, 
            "\n", sep="\t", file=outfile.name, append=T)
      }
    }
  } else {
    cat("ERROR! No data for model",
        model.name, "\n", file=paste("../data/output/interaction_errors_", run.name, ".txt", sep=""), append=T)
  }
  return (p.val)
}

# Run the models
################
snp.names <- names(tb.assoc)
for (i in snp.start.col:(dim(tb.assoc)[2]-1)) {
  snp1 <- tb.assoc[,i]
  name.snp1 <- snp.names[i]
  gene1.name <-
    genenames[genenames$SnpName==name.snp1, 1]
  for (j in (i+1):dim(tb.assoc)[2]) {
    snp2 <- tb.assoc[,j]
    name.snp2 <- snp.names[j]
    gene2.name <-
      genenames[genenames$SnpName==name.snp2, 1]
    if (gene1.name != gene2.name) {
      test.data <- tb.assoc[,c(1:10, which(names(tb.assoc) == name.snp1), which(names(tb.assoc) == name.snp2))]
      test.data <- test.data[complete.cases(test.data),]
      if (!(gene1.name %in% x.chr.genes) & !(gene2.name %in% x.chr.genes)) {
        model.name <- paste(gene1.name, "_", name.snp1, " - ", 
                            gene2.name, "_", name.snp2, sep="")
        p.val <- FitInteractionModel(model.name, test.data)
        if (!is.na(p.val)) {
          p.vals <- c(p.vals, as.numeric(p.val))
        }
      } else {
        model.name <- paste(gene1.name, "_", name.snp1, " - ", 
                            gene2.name, "_", name.snp2, " (Males)", sep="")
        test.data.m <- test.data[test.data$Gender == 1,]
        p.val <- FitInteractionModel(model.name, test.data.m, F)
        if (!is.na(p.val)) {
          p.vals <- c(p.vals, as.numeric(p.val))
        }  
        model.name <- paste(gene1.name, "_", name.snp1, " - ", 
                            gene2.name, "_", name.snp2, " (Females)", sep="")
        test.data.f <- test.data[test.data$Gender == 2,]
        p.val <- FitInteractionModel(model.name, test.data.f, F)
        if (!is.na(p.val)) {
          p.vals <- c(p.vals, as.numeric(p.val))
        }
      }    
    } 
  }
}

print(paste("Nr tests:", length(p.vals)))
#source("http://bioconductor.org/biocLite.R") #TODO: Uncomment to install qvalue package
#biocLite("qvalue") #TODO: Uncomment to install qvalue package
library(qvalue)
q.vals <- sort(qvalue(p.vals, pi0.method="bootstrap")$qvalues)
results <- read.delim("../data/output/interaction_summary.txt", head=F)
names(results) <- c("Model", "NrCases", "NrControls", "P", "Q")
results <- results[order(results$P),]
results$Q <- q.vals[1:dim(results)[1]]
write.table(results, "../data/output/interaction_summary.txt", row.names=F, quote=F, sep='\t')

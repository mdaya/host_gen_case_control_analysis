#Note: for now this script does not adjust for ancestry, but that may
#change in future

################################################################################
# Setup parameters, load packages and read the input file
################################################################################
input.file.name <- "../data/input/TB_assoc.txt"
outfile.name <- paste("../data/output/snp_summary.txt")
outfile.monosnp.name <- "../data/output/monomorphic_snps.txt" #If there are monomorphich SNPs they will be saved in this file
excl.snp.names <- c() #To exclude polymorphims that are not SNPs, e.g. indels. If you just have SNPs in your data set, keep this empty
begin.col <- 11 # where the SNP starts (col nr). This might change if you add/remove covariates

library(genetics)
sample.geno.data <- read.delim(input.file.name)

################################################################################
# Functions to get allele and genotype count from the genetics genotype()
# object
################################################################################
GetGenotypeCount <- function(freq, geno.name) {
  if (!is.null(rownames(freq))) {
    geno.count <- 
      freq[
        rownames(freq) %in% geno.name, 1]
    geno.count <- ifelse(length(geno.count) == 0, 0, geno.count)
  }
  #If the rownames are null, there is only 1 type of allele
  if (is.null(rownames(freq))) {
    if (geno.name[1] == paste(wild.a, wild.a, sep="/")) {
      geno.count <- as.numeric(freq[1])
    }
    if (geno.name[1] != paste(wild.a, wild.a, sep="/")) {
      geno.count <- 0
    }
  }
  return (geno.count)
}

GetGenotypeFreq <- function(freq, geno.name) {
  if (!is.null(rownames(freq))) {
    geno.freq <- 
      freq[
        rownames(freq) %in% geno.name, 2]
    geno.freq <- ifelse(length(geno.freq) == 0, 0, geno.freq)
  }
  #If the rownames are null, there is only 1 type of allele
  if (is.null(rownames(freq))) {
    if (geno.name[1] == paste(wild.a, wild.a, sep="/")) {
      geno.freq <- as.numeric(freq[2])
    }
    if (geno.name[1] != paste(wild.a, wild.a, sep="/")) {
      geno.freq <- 0
    }
  }
  return (geno.freq)
}

GetAlleleCount <- function(freq, allele.name) {
  if (!is.null(rownames(freq))) {
    allele.count <- freq[rownames(freq) == allele.name, 1]
    allele.count <- ifelse(length(allele.count) == 0, 0, allele.count)
  }
  if (is.null(rownames(freq))) {
    if (allele.name == wild.a) {
      allele.count <- as.numeric(freq[1])			
    }
    if (allele.name == rare.a) {
      allele.count <- 0
    }
  }
  return (allele.count)
}

GetAlleleFreq <- function(freq, allele.name) {
  if (!is.null(rownames(freq))) {
    allele.freq <- freq[rownames(freq) == allele.name, 2]
    allele.freq <- ifelse(length(allele.freq) == 0, 0, allele.freq)
  }
  if (is.null(rownames(freq))) {
    if (allele.name == wild.a) {
      allele.freq <- as.numeric(freq[2])  			
    }
    if (allele.name == rare.a) {
      allele.freq <- 0
    }
  }
  return (allele.freq)
}


################################################################################
# Create summary statistics
################################################################################
mono.frame <- data.frame()
i <- begin.col
for (i in (begin.col:dim(sample.geno.data)[2])) {
#for (i in (begin.col:begin.col)) {
  #Get SNP name and SNP
  snp.name = colnames(sample.geno.data)[i]
  snp <- sample.geno.data[,i]
  is.monomorphic <- (length(levels(snp)) == 1)
  
  if (is.monomorphic & !(snp.name %in% excl.snp.names) ) {
    is.control <- sample.geno.data$Class==0 & !is.na(snp)
    is.case <- sample.geno.data$Class==1 & !is.na(snp)
    nr.controls <- length(snp[is.control])
    nr.cases <- length(snp[is.case])
    genotype.val <-  levels(snp)
    if (dim(mono.frame)[1] > 0) {
      mono.frame <- rbind(mono.frame, data.frame(
        snp.name=snp.name,
        controls.count=nr.controls,
        cases.count=nr.cases,
        genotype=genotype.val,
        stringsAsFactors = F) )
    }  
    if (dim(mono.frame)[1] == 0) {
      mono.frame <- data.frame(
        snp.name=snp.name,
        controls.count=nr.controls,
        cases.count=nr.cases,
        genotype=genotype.val,
        stringsAsFactors = F) 
    }
  }
  
  if (!(snp.name %in% excl.snp.names) & !is.monomorphic) {
    
    #Get Case and Control Count  
    is.control <- sample.geno.data$Class==0 & !is.na(snp)
    is.case <- sample.geno.data$Class==1 & !is.na(snp)
    nr.controls <- length(snp[is.control])
    nr.cases <- length(snp[is.case])
    
    #Get summary genotype objects, HWE and rare and wildtype alleles
    gene.snp.controls <- genotype(snp[is.control], sep="")
    gene.snp.cases <- genotype(snp[is.case], sep="")  
    gene.snp <- genotype(snp[!is.na(snp)], sep="")
    if (length(levels(gene.snp.controls)) > 1) {
      hwe.controls.p <- HWE.exact(gene.snp.controls)$p.value
    }
    if (length(levels(gene.snp.controls)) <= 1) {
      hwe.controls.p <- -999
    }
    if (length(levels(gene.snp.cases)) > 1) {
      hwe.cases.p <- HWE.exact(gene.snp.cases)$p.value
    }
    if (length(levels(gene.snp.cases)) <= 1) {
      hwe.cases.p <- -999
    }
    snp.info <- summary(gene.snp)
    snp.info.controls <- summary(gene.snp.controls)
    snp.info.cases <- summary(gene.snp.cases)
    rare.a <- snp.info$allele.names[snp.info$allele.freq[,2] == min(snp.info$allele.freq[,2])]
    wild.a <- snp.info$allele.names[snp.info$allele.freq[,2] == max(snp.info$allele.freq[,2])]
    het1.g <- paste(rare.a, wild.a, sep="")
    het2.g <- paste(wild.a, rare.a, sep="")
    homo.rare.g <- paste(rare.a, rare.a, sep="")
    homo.wild.g <- paste(wild.a, wild.a, sep="")  
    
    #Get Geno and Allele P vals
    trait <- sample.geno.data$Class
    age <- sample.geno.data$Age
    gender <- sample.geno.data$Gender
    san <- sample.geno.data$SanAncestry
    afr <- sample.geno.data$AfricanAncestry
    eur <- sample.geno.data$EuropeanAncestry
    sas <- sample.geno.data$SouthAsianAncestry
    trait <- trait[!(is.na(snp))]
    age <- age[!(is.na(snp))]
    gender <- gender[!(is.na(snp))]
    san <- san[!(is.na(snp))]
    afr <- afr[!(is.na(snp))]
    eur <- eur[!(is.na(snp))]
    sas <- sas[!(is.na(snp))]
    snp <- snp[!(is.na(snp))]
    
    snp.a <- rep(0, length(snp))
    snp.a[(snp == het1.g) | (snp == het2.g)] <- 1
    snp.a[snp == homo.rare.g] <- 2
    test <- anova(glm(trait ~ snp, family="binomial"), test="Chisq")
    p.unadj.geno <- test[2,5]
    test <- anova(glm(trait ~ age + gender + san + afr + eur + sas + snp, family="binomial"), test="Chisq")
    p.adj.geno <- test[8,5]
    test <- anova(glm(trait ~ snp.a, family="binomial"), test="Chisq")
    p.unadj.a <- test[2,5]
    test <- anova(glm(trait ~ age + gender + san + afr + eur + sas + snp.a, family="binomial"), test="Chisq")
    p.adj.a <- test[8,5]
    
    #Write the SNP names and case and control counts
    if (i == begin.col) {
      output.frame <- data.frame(
        snp.name=snp.name,
        controls.count="",
        controls.freq="",
        controls.hwe.p="",
        cases.count="",
        cases.freq="",
        cases.hwe.p="",
        unadj.p="",
        adj.p="",
        stringsAsFactors = F)
    }
    if (i > begin.col) {
      output.frame <- rbind(output.frame, data.frame(
        snp.name=snp.name,
        controls.count="",
        controls.freq="",
        controls.hwe.p="",
        cases.count="",
        cases.freq="",
        cases.hwe.p="",
        unadj.p="",
        adj.p="",
        stringsAsFactors = F) )
    }
    output.frame <- rbind(output.frame, data.frame(
      snp.name="",
      controls.count=nr.controls,
      controls.freq="",
      controls.hwe.p=round(hwe.controls.p,3),
      cases.count=nr.cases,
      cases.freq="",
      cases.hwe.p=round(hwe.cases.p,3),
      unadj.p="",
      adj.p="",
      stringsAsFactors = F) )
    
    #Get Geno Counts and write them
    geno.names <- rownames(snp.info$genotype.freq)
    for (k in 1:length(geno.names)) {
      geno.name.l <- c(
        geno.names[k], 
        paste(substr(geno.names[k],3,3),"/",substr(geno.names[k],1,1), sep=""))
      geno.name <- geno.name.l[1]
      geno.controls.count <- GetGenotypeCount(snp.info.controls$genotype.freq, geno.name.l)
      geno.controls.freq <- GetGenotypeFreq(snp.info.controls$genotype.freq, geno.name.l)
      geno.cases.count <- GetGenotypeCount(snp.info.cases$genotype.freq, geno.name.l)
      geno.cases.freq <- GetGenotypeFreq(snp.info.cases$genotype.freq, geno.name.l) 
      if (k == 1) {
        unadj.p=round(p.unadj.geno,3)
        adj.p=round(p.adj.geno,3)    
      }
      if (k > 1 ) {
        unadj.p=""
        adj.p="" 
      }      
      output.frame <- rbind(output.frame, data.frame(
        snp.name=geno.name,
        controls.count=geno.controls.count,
        controls.freq=round(geno.controls.freq,2),
        controls.hwe.p="",
        cases.count=geno.cases.count,
        cases.freq=round(geno.cases.freq,2),
        cases.hwe.p="",
        unadj.p=unadj.p, 
        adj.p=adj.p, 
        stringsAsFactors = F) )
    }
    
    allele.names <- rownames(snp.info$allele.freq)
    for (k in 1:length(allele.names)) {
      allele.name <- allele.names[k]
      allele.controls.count <- GetAlleleCount(snp.info.controls$allele.freq, allele.name)
      allele.controls.freq <- GetAlleleFreq(snp.info.controls$allele.freq, allele.name)
      allele.cases.count <- GetAlleleCount(snp.info.cases$allele.freq, allele.name)
      allele.cases.freq <- GetAlleleFreq(snp.info.cases$allele.freq, allele.name) 
      if (k == 1) {
        unadj.p=round(p.unadj.a,3)
        adj.p=round(p.adj.a,3)    
      }
      if (k > 1 ) {
        unadj.p=""
        adj.p="" 
      }      
      output.frame <- rbind(output.frame, data.frame(
        snp.name=allele.name,
        controls.count=allele.controls.count,
        controls.freq=round(allele.controls.freq,2),
        controls.hwe.p="",
        cases.count=allele.cases.count,
        cases.freq=round(allele.cases.freq,2),
        cases.hwe.p="",
        unadj.p=unadj.p, 
        adj.p=adj.p, 
        stringsAsFactors = F) )
    }
  
  }
}

################################################################################
# Write the SNP summary to an output file
################################################################################
#TODO: paste this an excel file
write.table(output.frame, outfile.name,sep="\t",row.names=F, col.names=F, quote=F)

if (dim(mono.frame)[1] > 0) {
  write.table(mono.frame, outfile.monosp.name, sep="\t",row.names=F, col.names=F, quote=F)
}
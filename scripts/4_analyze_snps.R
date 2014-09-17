input.file.name <- "../data/input/TB_assoc.txt"
library(genetics)
sample.geno.data <- read.delim(input.file.name)
#TODO: change the SNP name according to the SNP you want to analyze
snp <- sample.geno.data$TLR9_rs352139

#Set the variables to put in the model
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

#Get rare allele and code snp.a, snp.d, snp.r as 0, 1, or 2 according to the number of copies of the rare variant
gene.snp <- genotype(snp, sep="")
snp.info <- summary(gene.snp)
rare.a <- snp.info$allele.names[snp.info$allele.freq[,2] == min(snp.info$allele.freq[,2])]
wild.a <- snp.info$allele.names[snp.info$allele.freq[,2] == max(snp.info$allele.freq[,2])]
het1.g <- paste(rare.a, wild.a, sep="")
het2.g <- paste(wild.a, rare.a, sep="")
homo.rare.g <- paste(rare.a, rare.a, sep="")
snp.a <- rep(0, length(snp))
snp.a[(snp == het1.g) | (snp == het2.g)] <- 1
snp.a[snp == homo.rare.g] <- 2
snp.d <- rep(0, length(snp))
snp.d[(snp == het1.g) | (snp == het2.g)] <- 1
snp.d[snp == homo.rare.g] <- 1
snp.r <- rep(0, length(snp))
snp.r[(snp == het1.g) | (snp == het2.g)] <- 0
snp.r[snp == homo.rare.g] <- 1

#TODO: Check whether the genotype needs to be releveled.
#The wild type homozygote should be the first genotype level
#The relevel command should be commented out and updated with the relevant reference
#genotyp if required
table(snp)
levels(snp)
#snp <- relevel(snp, ref="GG")

#Fit the models
#We use logistic regression to model the relationship between TB case/control 
#status and the predictor variables. 
geno.test <- glm(trait ~ age + gender + san + afr + eur + sas + snp, family="binomial")
allele.test.a <- glm(trait ~ age + gender + san + afr + eur + sas + snp.a, family="binomial")
allele.test.d <- glm(trait ~ age + gender + san + afr + eur + sas + snp.d, family="binomial")
allele.test.r <- glm(trait ~ age + gender + san + afr + eur + sas + snp.r, family="binomial")

#TODO: run the below lines one by one and decide which allelic model 
#(a = additive, d=dominant, r=recessive) is best. That would be the
#model with the smallest P-value of the snp.x row (column=Pr(>Chi))
anova(allele.test.a, test="Chisq")
anova(allele.test.d, test="Chisq")
anova(allele.test.r, test="Chisq")

#TODO: update this to contain the best model (this could also be
#the geno.test model)
best.model <- allele.test.a

#TODO: run the below and get the effect size of the SNP, from the 
#Coefficients: Estimate column
summary(best.model)

#TODO: update the effect size according to the summary above
effect.size <-  0.3609568

#TODO: run the below to get the confidence interval for the OR
or <- exp(effect.size)
ci <- exp(confint(best.model))
print(or)
print(ci)

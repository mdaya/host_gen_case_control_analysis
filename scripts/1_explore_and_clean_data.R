snp.in.data <- read.delim("../data/raw/SNPs.txt")

#Check that the data has been read correctly 
head(snp.in.data)

#Check for incorrect values and fix them appropriately
summary(snp.in.data[,-1])

#Below commented out commands is examples of fixes that could
#be made to the data
# levels(snp.in.data$TLR4_rs4986791)
# snp.in.data$SampleID[snp.in.data$TLR4_rs4986791=='GG'] 
# levels(snp.in.data$TLR4_rs4986791)[3] <- 'TT'
# snp.in.data$TLR4_rs4986791[snp.in.data$SampleID==169] <- 'TT'
# summary(snp.in.data[,-1])
# 
# levels(snp.in.data$TLR4_rs11536889)
# snp.in.data$TLR4_rs11536889[snp.in.data$TLR4_rs11536889==" "] <- NA
# snp.in.data$TLR4_rs11536889[snp.in.data$TLR4_rs11536889=="  "] <- NA
# summary(snp.in.data[,-1])
# 
# levels(snp.in.data$TLR4_rs11536891) 
# snp.in.data$TLR4_rs11536891[snp.in.data$TLR4_rs11536891=="T"] <- 'TT'
# snp.in.data$TLR4_rs11536891[snp.in.data$TLR4_rs11536891=="GG"] <- NA
# summary(snp.in.data[,-1])

#Merge the data with sample_info.txt 
sample.info <- read.delim("../data/raw/sample_info.txt")
head(sample.info) #Check that the data has been read correctly
snp.out.data <- merge(sample.info, snp.in.data, all.y=T)
#Check that the number of rows of the input and output files match
dim(snp.in.data)
dim(snp.out.data)
#Check that the new data have no missing Case IDs.
#The sum should be zero. If not, investigate and remove the duplicate
#sample/or update the sample nr
sum(is.na(snp.out.data$CaseID)) 
snp.out.data[is.na(snp.out.data$CaseID),]

#Check for duplicate individuals and remove them if necessary.
#The sum should be zero. If not, investigate and remove the duplicate
#sample.
sum(duplicated(snp.out.data$CaseID))

#Remove samples with missing covariates
snp.out.data <- snp.out.data[complete.cases(snp.out.data[,1:10]),]

#Write the analysis input file
write.table(snp.out.data, "../data/input/TB_assoc.txt", row.names=F, sep="\t", quote=F)

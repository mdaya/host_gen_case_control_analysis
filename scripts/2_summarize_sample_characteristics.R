################################################################################
# Setup parameters and read the file
################################################################################
input.file.name <- "../data/input/TB_assoc.txt"
output.file.name <- "../data/output/sample_info_summary.txt"
output.image.gender.age <- "../data/output/sample_gender_age.pdf"
output.image.ancestry <- "../data/output/sample_ancestry.pdf"
cols <- c(3:10)

sample.info <- read.delim(input.file.name)[,cols]

################################################################################
# Calculate statistics
################################################################################
is.case <- sample.info$Class==1
is.control <- sample.info$Class==0
nr.cases <- sum(is.case, na.rm=T)
nr.controls <- sum(is.control, na.rm=T)

is.male <- !is.na(sample.info$Gender) & (sample.info$Gender == 1)
nr.male.cases <- sum(is.case & is.male)
nr.male.controls <- sum(is.control & is.male )
prop.male.cases <- nr.male.cases/nr.cases
prop.male.controls <- nr.male.controls/nr.cases

mean.age.cases <- mean(sample.info$Age[is.case], na.rm=T)
sd.age.cases <- sd(sample.info$Age[is.case], na.rm=T)
mean.age.controls <- mean(sample.info$Age[is.control], na.rm=T)
sd.age.controls <- sd(sample.info$Age[is.control], na.rm=T)

median.san.cases <- quantile(sample.info[is.case,4])[3]
q1.san.cases <- quantile(sample.info[is.case,4])[2]
q3.san.cases <- quantile(sample.info[is.case,4])[4]
median.san.controls <- quantile(sample.info[is.control,4])[3]
q1.san.controls <- quantile(sample.info[is.control,4])[2]
q3.san.controls <- quantile(sample.info[is.control,4])[4]

median.afr.cases <- quantile(sample.info[is.case,5])[3]
q1.afr.cases <- quantile(sample.info[is.case,5])[2]
q3.afr.cases <- quantile(sample.info[is.case,5])[4]
median.afr.controls <- quantile(sample.info[is.control,5])[3]
q1.afr.controls <- quantile(sample.info[is.control,5])[2]
q3.afr.controls <- quantile(sample.info[is.control,5])[4]

median.eur.cases <- quantile(sample.info[is.case,6])[3]
q1.eur.cases <- quantile(sample.info[is.case,6])[2]
q3.eur.cases <- quantile(sample.info[is.case,6])[4]
median.eur.controls <- quantile(sample.info[is.control,6])[3]
q1.eur.controls <- quantile(sample.info[is.control,6])[2]
q3.eur.controls <- quantile(sample.info[is.control,6])[4]

median.sas.cases <- quantile(sample.info[is.case,7])[3]
q1.sas.cases <- quantile(sample.info[is.case,7])[2]
q3.sas.cases <- quantile(sample.info[is.case,7])[4]
median.sas.controls <- quantile(sample.info[is.control,7])[3]
q1.sas.controls <- quantile(sample.info[is.control,7])[2]
q3.sas.controls <- quantile(sample.info[is.control,7])[4]

median.eas.cases <- quantile(sample.info[is.case,8])[3]
q1.eas.cases <- quantile(sample.info[is.case,8])[2]
q3.eas.cases <- quantile(sample.info[is.case,8])[4]
median.eas.controls <- quantile(sample.info[is.control,8])[3]
q1.eas.controls <- quantile(sample.info[is.control,8])[2]
q3.eas.controls <- quantile(sample.info[is.control,8])[4]

model <- glm(Class ~ Age + Gender + SanAncestry + AfricanAncestry + 
               EuropeanAncestry + SouthAsianAncestry, data=sample.info, 
             family="binomial")
age.p <- summary(model)$coefficients[2,4]
gender.p <- summary(model)$coefficients[3,4]
san.p <- summary(model)$coefficients[4,4]
afr.p <- summary(model)$coefficients[5,4]
eur.p <- summary(model)$coefficients[6,4]
sas.p <- summary(model)$coefficients[7,4]

################################################################################
# Functions to format output
################################################################################
FormatEst <- function(nr) {
  return (formatC(nr, 2, format="f"))
}

FormatP <- function(nr) {
  if (nr >= 0.0001 ) {
    return (formatC(nr, 4, format="f"))
  } else {
    return ("< 0.0001")
  }
}

################################################################################
# Write the output
################################################################################
output.frame <- data.frame(
    col1="",
    col2="TB cases",
    col3="Controls",
    col4=""
  )
output.frame <- rbind(output.frame, data.frame(
  col1="",
  col2=paste("(n=", nr.cases, ")", sep=""),
  col3=paste("(n=", nr.controls, ")", sep=""),
  col4="P-value"
))
output.frame <- rbind(output.frame, data.frame(
  col1="Age (mean +- SD)",
  col2=paste(FormatEst(mean.age.cases), " +- ", FormatEst(sd.age.cases), sep=""),
  col3=paste(FormatEst(mean.age.controls), " +- ", FormatEst(sd.age.controls), sep=""),
  col4=FormatP(age.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="Nr males (prop)",
  col2=paste(nr.male.cases, " (", FormatEst(prop.male.cases), ")", sep=""),
  col3=paste(nr.male.controls, " (", FormatEst(prop.male.controls), ")", sep=""),
  col4=FormatP(gender.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="African San [IQR]",
  col2=paste(FormatEst(median.san.cases), " [", FormatEst(q1.san.cases), "-", FormatEst(q3.san.cases), "]", sep=""),
  col3=paste(FormatEst(median.san.controls), " [", FormatEst(q1.san.controls), "-", FormatEst(q3.san.controls), "]", sep=""),
  col4=FormatP(san.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="African non-San [IQR]",
  col2=paste(FormatEst(median.afr.cases), " [", FormatEst(q1.afr.cases), "-", FormatEst(q3.afr.cases), "]", sep=""),
  col3=paste(FormatEst(median.afr.controls), " [", FormatEst(q1.afr.controls), "-", FormatEst(q3.afr.controls), "]", sep=""),
  col4=FormatP(afr.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="European [IQR]",
  col2=paste(FormatEst(median.eur.cases), " [", FormatEst(q1.eur.cases), "-", FormatEst(q3.eur.cases), "]", sep=""),
  col3=paste(FormatEst(median.eur.controls), " [", FormatEst(q1.eur.controls), "-", FormatEst(q3.eur.controls), "]", sep=""),
  col4=FormatP(eur.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="South Asian [IQR]",
  col2=paste(FormatEst(median.sas.cases), " [", FormatEst(q1.sas.cases), "-", FormatEst(q3.sas.cases), "]", sep=""),
  col3=paste(FormatEst(median.sas.controls), " [", FormatEst(q1.sas.controls), "-", FormatEst(q3.sas.controls), "]", sep=""),
  col4=FormatP(sas.p)
))
output.frame <- rbind(output.frame, data.frame(
  col1="East Asian [IQR]",
  col2=paste(FormatEst(median.eas.cases), " [", FormatEst(q1.eas.cases), "-", FormatEst(q3.eas.cases), "]", sep=""),
  col3=paste(FormatEst(median.eas.controls), " [", FormatEst(q1.eas.controls), "-", FormatEst(q3.eas.controls), "]", sep=""),
  col4=""
))
#TODO: paste the output file contents into an excel file
write.table(output.frame, output.file.name, quote=F, sep="\t", col.names = F, row.names=F)

################################################################################
# Create summary graphs
################################################################################
pdf(output.image.gender.age)
par(mfrow=c(2,2))

#Barplot of gender
class <- sample.info$Class
class[class==1] <- -1
gender.table <- table(sample.info$Gender,class) 
rownames(gender.table) <- c("Male", "Female")
colnames(gender.table) <- c("TB cases", "Controls")
barplot(gender.table, beside=T, col=c("grey30", "grey80"), ylab="Number")
legend(1, 300, rownames(gender.table), col=c("grey30", "grey80"), pch=c(15,15), cex=0.8)

#Boxplot of age
boxplot(sample.info$Age[is.case], sample.info$Age[is.control], 
  names=c("TB cases", "Controls"), ylab="Age", col=c("grey30", "grey80"))
dev.off()

#Boxplots of ancestry proportions
pdf(output.image.ancestry)
par(mfrow=c(1,1))
boxplot(sample.info[is.case,4:8], outline=T,ylab="Proportion Ancestry", boxwex=0.25, at=1:5-0.2, col="grey30", xaxt="n")
boxplot(sample.info[is.control,4:8], outline=T, add=T, boxwex=0.25, at=1:5+0.2, col="grey80", xaxt="n")
axis(1,1:5, labels=c("San","African","European","South Asian","East Asian"), cex.axis=0.5)
legend(3.5,0.8,c("TB cases","Controls"), fill =c("grey30","grey80"), cex=0.7)
dev.off()



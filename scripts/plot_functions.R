library(effects)
library(ggplot2)
library(gridExtra)
library(haplo.stats)
library(reshape2)

GetGenoCombinationPlot <- function(snp1, snp2, tb.data) {
  gene1 <- genenames[genenames$SnpName==snp1, 1]
  gene2 <- genenames[genenames$SnpName==snp2, 1]
  xlab <- paste(gene2, snp2)
  flab <- paste(gene1, snp1)
  prop.data <- tb.data[,which(names(tb.data) %in% c(snp1, snp2, "Class"))]
  prop.data <- prop.data[complete.cases(prop.data),c(2,3,1)]
  d <- as.data.frame(prop.table(table(prop.data),3))
  table.1 <- table(prop.data[,1])
  if ( (length(table.1) == 3) & (table.1[3] > table.1[1]) ) {
    lev <- levels(d[,1])
    d[,1]<- factor(d[,1], levels = c(lev[3], lev[2], lev[1]))
  }
  table.2 <- table(prop.data[,2])
  if ( (length(table.2) == 3) & (table.2[3] > table.2[1]) ) {
    lev <- levels(d[,2])
    d[,2]<- factor(d[,2], levels = c(lev[3], lev[2], lev[1]))
  } 
  if (sum(levels(d[,1]) %in% c("TA", "TC", "TG", "GA", "GC")) > 0) {
    levels(d[,1])[levels(d[,1]) == "TA"] <- "AT"
    levels(d[,1])[levels(d[,1]) == "TC"] <- "CT"
    levels(d[,1])[levels(d[,1]) == "TG"] <- "GT"
    levels(d[,1])[levels(d[,1]) == "GA"] <- "AG"
    levels(d[,1])[levels(d[,1]) == "GC"] <- "CG"
  }
  if (sum(levels(d[,2]) %in% c("TA", "TC", "TG", "GA", "GC")) > 0) {
    levels(d[,2])[levels(d[,2]) == "TA"] <- "AT"
    levels(d[,2])[levels(d[,2]) == "TC"] <- "CT"
    levels(d[,2])[levels(d[,2]) == "TG"] <- "GT"
    levels(d[,2])[levels(d[,2]) == "GA"] <- "AG"
    levels(d[,2])[levels(d[,2]) == "GC"] <- "CG"
  }
  class <- rep(NA, dim(d)[1])
  class[d$Class == 0] <- "Controls"
  class[d$Class == 1] <- "TB cases"
  geno.comb <- data.frame(geno1 = d[,which(names(d) == snp1)], geno2 = d[,which(names(d) == snp2)], prop = d$Freq, class)
  g <- ggplot(data=geno.comb, aes(x=geno2, y=prop, fill=geno1))
  b <- geom_bar(stat="identity", position=position_dodge())
  f <- facet_wrap(~class, nrow=1, ncol=4)
  l <- labs(title="", y="", x=xlab, fill=flab)
  c <- scale_fill_brewer(palette="Set1")
  t <- theme_bw() + theme(legend.position="top", 
                          text = element_text(size=8, face="plain"),
                          legend.title = element_text(size=8, face="plain"), 
                          legend.key.size = unit(2.5, "mm"))
  p <- g + b + f + l + c + t
  return (p)
}


GetAlleleCombinationPlot <- function(snp1, snp2, tb.data, draw.legend=F, lev=c()) {
  gene1 <- genenames[genenames$SnpName==snp1, 1]
  gene2 <- genenames[genenames$SnpName==snp2, 1]
  snp.names <- c(snp1, snp2)
  prop.data <- tb.data[,which(names(tb.data) %in% c(snp1, snp2, "Class"))]
  prop.data <- prop.data[complete.cases(prop.data),]
  single.alleles <- cbind(substr(prop.data[,snp.names[1]], 1, 1),   
                          substr(prop.data[,snp.names[1]], 2, 2) )  
  for (k in (2:length(snp.names)) ) {
    single.alleles <- cbind(single.alleles,
                            substr(prop.data[,snp.names[k]], 1, 1),   
                            substr(prop.data[,snp.names[k]], 2, 2) )  
  }
  geno <- setupGeno(single.alleles, miss.val=c(NA), locus.label=snp.names)
  Class <- prop.data$Class
  group.bin <- haplo.group(Class, geno, locus.label=snp.names)
  alleles <- attributes(geno)$unique.alleles
  allele.pairs <- c(
    paste(alleles[[1]][as.numeric(group.bin$group.df[1,1])],
          alleles[[2]][as.numeric(group.bin$group.df[1,2])], 
          sep="-"),
    paste(alleles[[1]][as.numeric(group.bin$group.df[2,1])],
          alleles[[2]][as.numeric(group.bin$group.df[2,2])], 
          sep="-"),
    paste(alleles[[1]][as.numeric(group.bin$group.df[3,1])],
          alleles[[2]][as.numeric(group.bin$group.df[3,2])], 
          sep="-"),
    paste(alleles[[1]][as.numeric(group.bin$group.df[4,1])],
          alleles[[2]][as.numeric(group.bin$group.df[4,2])], 
          sep="-"))
  allele.comb <- data.frame(
    alleles = rep(allele.pairs,2),
    prop = unlist(c(group.bin$group.df[4], group.bin$group.df[5])),
    class =  c(rep("Controls", 4), rep("TB cases", 4)) )
  if (length(lev) > 0) {
    allele.comb$alleles <- factor(allele.comb$alleles,lev)
  }
  title.text = ""
  xlab = paste(gene1, snp1, "-", gene2, snp2)
  g <- ggplot(data=allele.comb, aes(x=alleles, y=prop, fill=class))
  b <- geom_bar(stat="identity", position=position_dodge())
  l <- labs(y="", x=xlab, fill="", title=title.text)
  c <- scale_fill_brewer(palette="Set1")
  t <-  theme_bw() + theme(text = element_text(size=8, face="plain"))
  if (draw.legend) {
    t <- theme_bw() + theme(text = element_text(size=8, face="plain"),
                            legend.key.size = unit(2.5, "mm"),
                            legend.position = c(0.72, 0.7))
    p <- g + b + l + c + t 
  } else {
    p <- g + b + l + c + t + guides(fill=F)
  }
  return(p)
}

GetLogitLabel <- function(prop) {
  logit <- round(log(prop/(1-prop)), 1)
  return (formatC(logit, 1, format="f"))
}

GetLogitPlot <- function(snp1, snp2, tb.data) {
  gene1 <- genenames[genenames$SnpName==snp1, 1]
  gene2 <- genenames[genenames$SnpName==snp2, 1]
  test.data <- tb.data[,c(2:10,which(names(tb.data) %in% c(snp1, snp2)))]
  names(test.data)[10] <- "snp1"
  names(test.data)[11] <- "snp2"
  model <- glm(Class ~ Age + Gender + SanAncestry + AfricanAncestry + EuropeanAncestry + SouthAsianAncestry + snp1*snp2, family="binomial", data=test.data)
  e <- summary(effect("snp1:snp2", model))$effect
  d <- melt(e)
  table.1 <- table(test.data$snp1)
  if ( (length(table.1) == 3) & (table.1[3] > table.1[1]) ) {
    lev <- levels(d$snp1)
    d$snp1 <- factor(d$snp1, levels = c(lev[3], lev[2], lev[1]))
  }
  table.2 <- table(test.data$snp2)
  if ( (length(table.2) == 3) & (table.2[3] > table.2[1]) ) {
    lev <- levels(d$snp2)
    d$snp2 <- factor(d$snp2, levels = c(lev[3], lev[2], lev[1]))
  }
  if (sum(levels(d[,1]) %in% c("TA", "TC", "TG", "GA", "GC")) > 0) {
    levels(d[,1])[levels(d[,1]) == "TA"] <- "AT"
    levels(d[,1])[levels(d[,1]) == "TC"] <- "CT"
    levels(d[,1])[levels(d[,1]) == "TG"] <- "GT"
    levels(d[,1])[levels(d[,1]) == "GA"] <- "AG"
    levels(d[,1])[levels(d[,1]) == "GC"] <- "CG"
  }
  if (sum(levels(d[,2]) %in% c("TA", "TC", "TG", "GA", "GC")) > 0) {
    levels(d[,2])[levels(d[,2]) == "TA"] <- "AT"
    levels(d[,2])[levels(d[,2]) == "TC"] <- "CT"
    levels(d[,2])[levels(d[,2]) == "TG"] <- "GT"
    levels(d[,2])[levels(d[,2]) == "GA"] <- "AG"
    levels(d[,2])[levels(d[,2]) == "GC"] <- "CG"
  }
  g <- ggplot(data=d, aes(x=snp2, y=value, group=snp1, colour=snp1, linetype=snp1, shape=snp1)) + 
    geom_line() + geom_point()
  t <- theme_bw() + theme(legend.position="top",
                          text = element_text(size=8, face="plain"),
                          plot.title = element_text(size=8, face="plain"), 
                          legend.key = element_blank(),
                          legend.margin = unit(-2, "mm"),
                          legend.key.height = unit(2, "mm"),
                          legend.key.width = unit(10, "mm"))
  c <- scale_color_brewer(palette="Set1")
  l <- labs(title=paste(gene1, " ", snp1, sep=""), 
            y="logit", x=paste(gene2, snp2),
            colour="", linetype="", shape="") #, colour=paste(gene1, snp1.name))
  lower <- range(e, na.rm=T)[1]
  upper <- range(e, na.rm=T)[2] 
  if (upper > 0.975) {
    upper <- 0.975
  }
  incr <- (upper - lower)/4
  breaks <- seq(lower, upper, incr)
  s <- scale_y_continuous(breaks = breaks, 
                          labels =  GetLogitLabel(breaks)) 
  lt <- scale_linetype_manual(values=c("solid", "dashed", "dotted"))
  p <- g + t + c + l + s + lt
  return(p)
}
library(genetics)
library(haplo.stats)
seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

#-------------------------------------------------------------------------------
SetupGeno <- function(data, snp.names) {
	# Create dataframe with columns of single alleles. E.g. 2 columns with 
	# genotypes AG CT will become 4 columns with alleles A G C T
	single.alleles <- cbind(substr(data[,snp.names[1]], 1, 1),   
									substr(data[,snp.names[1]], 2, 2) )  
	for (i in (2:length(snp.names)) ) {
		single.alleles <- cbind(single.alleles,
										substr(data[,snp.names[i]], 1, 1),   
										substr(data[,snp.names[i]], 2, 2) )  
	}

	# Return geno data structure required by haplo.glm
	return (setupGeno(single.alleles, miss.val=c(NA), locus.label=snp.names))
}

#-------------------------------------------------------------------------------
RunHaplo <- function(infile.name, outfile.name, snp.names, x.chr=F) {
	# Given a dataset of SNPs, if SNPs are located in different genes, 
	# run some haplotype summary statistics, and fit a haplotype model
	# Args:
	#	in.file.name: The input file name. The input file should
	#		contain the following columns:
	#			Class: Set to either 0 (control) or 1 (case)
	#			Age
	#			Gender: Male=1, Female=2
	#			SnpName1, SnpName2, ... : SNPs should be 2 letter genotypes,
	#				e.g. AC 
	#	out.file.name: The output file name 
	#	snp.names: The names of the SNPs that should form a haplotype
	#	x.chr: Whether this gene is on the X chromosome
 
 	#Set up
	sample.geno.data <- read.delim(infile.name)
	nr.snps <- length(snp.names)
	
	#Fit model
	geno <- SetupGeno(sample.geno.data, snp.names)
	y <- sample.geno.data$Class
	age <- sample.geno.data$Age
	gender <- sample.geno.data$Gender
	san <- sample.geno.data$SanAncestry
	afr <- sample.geno.data$AfricanAncestry
	eur <- sample.geno.data$EuropeanAncestry
	sas <- sample.geno.data$SouthAsianAncestry
	model.data <- data.frame(geno, age=age, gender=gender, san=san, afr=afr, eur=eur, sas=sas, y=y)
	fit.bin <- haplo.glm(y ~ age + gender + san + afr + eur + sas + geno,
		family=binomial,
		data=model.data,
		locus.label=snp.names)
		
	#Get P values for the model (adjusted and unadjusted)
	score.unadj <- haplo.score(
		y, 
		geno,
		trait.type="binomial",
		haplo.effect="additive",
		locus.label=snp.names,
		miss.val=0,
		#sim.control=score.sim.control(min.sim=20000),
		simulate=T)
	unadj.p <- score.unadj$score.global.p
	score.adj <- haplo.score(
		y, 
		geno,
		trait.type="binomial",
		x.adj=data.frame(age=sample.geno.data$Age, gender=sample.geno.data$Gender,
                     san=sample.geno.data$SanAncestry, afr=sample.geno.data$AfricanAncestry,
                     eur=sample.geno.data$EuropeanAncestry, sas=sample.geno.data$SouthAsianAncestry),
		haplo.effect="additive",
		locus.label=snp.names,
		miss.val=0,
		#sim.control=score.sim.control(min.sim=20000),
		simulate=T)
	adj.p <- score.adj$score.global.p
		
	#Run frequency summaries
	if (x.chr == F) {
		group.bin <- haplo.group(y, geno, locus.label=snp.names)
		alleles <- attr(geno,"unique.alleles")
		a.lookup <- c()
		for (k in (1:length(alleles))) {
			a1 <-  substr(alleles[k], 4, 4)
			a2 <-  substr(alleles[k], 9, 9)
			a.lookup <- rbind(a.lookup, c(a1, a2))
		}
		for (row.nr in (1:dim(group.bin$group.df)[1]) ) {
			haplo.name <- ""
			for (i in (1:nr.snps)) {
				val <- a.lookup[i, as.numeric(group.bin$group.df[row.nr, i]) ]
				if (i < nr.snps) {
					haplo.name <- paste(haplo.name, val, "/", sep="")
				}
				else {
					haplo.name <- paste(haplo.name, val, sep="")
				}
			}
			controls.freq <- group.bin$group.df[row.nr, nr.snps+2]
			cases.freq <- group.bin$group.df[row.nr, nr.snps+3]
			if (row.nr == 1) {
				output.frame <- data.frame(
					haplo.name=haplo.name,
					controls.freq=round(controls.freq,2),
					cases.freq=round(cases.freq,2),
					stringsAsFactors = F)
			}
			else {
				output.frame <- rbind(output.frame, data.frame(
					haplo.name=haplo.name,
					controls.freq=round(controls.freq,2),
					cases.freq=round(cases.freq,2),
					stringsAsFactors = F) )
			}
		}
	}
	if (x.chr == T) {
		sample.geno.data.f <- sample.geno.data[sample.geno.data$Gender==2, ]
		sample.geno.data.m <- sample.geno.data[sample.geno.data$Gender==1, ]

		output.frame <- data.frame(
       		haplo.name="Female",
        	controls.freq="",
        	cases.freq="",
        	stringsAsFactors = F)
		geno <- SetupGeno(sample.geno.data.f, snp.names)
		y <- sample.geno.data.f$Class
		group.bin <- haplo.group(y, geno, locus.label=snp.names)
		alleles <- attr(geno,"unique.alleles")
		a.lookup <- c()
		for (k in (1:length(alleles))) {
			a1 <-  substr(alleles[k], 4, 4)
        	a2 <-  substr(alleles[k], 9, 9)
        	a.lookup <- rbind(a.lookup, c(a1, a2))
		}
		for (row.nr in (1:dim(group.bin$group.df)[1]) ) {
        	haplo.name <- ""
        	for (i in (1:nr.snps)) {
                val <- a.lookup[i, as.numeric(group.bin$group.df[row.nr, i]) ]
                if (i < nr.snps) {
                        haplo.name <- paste(haplo.name, val, "/", sep="")
                }
                else {
                        haplo.name <- paste(haplo.name, val, sep="")
                }
        	}
        	controls.freq <- group.bin$group.df[row.nr, nr.snps+2]
        	cases.freq <- group.bin$group.df[row.nr, nr.snps+3]
        	output.frame <- rbind(output.frame, data.frame(
                haplo.name=haplo.name,
                controls.freq=round(controls.freq,2),
                cases.freq=round(cases.freq,2),
                stringsAsFactors = F) )
		}

		output.frame <- rbind(output.frame, data.frame(
        	haplo.name="Male",
        	controls.freq="",
        	cases.freq="",
        	stringsAsFactors = F) )
		geno <- SetupGeno(sample.geno.data.m, snp.names)
		y <- sample.geno.data.m$Class
		group.bin <- haplo.group(y, geno, locus.label=snp.names)
		alleles <- attr(geno,"unique.alleles")
		a.lookup <- c()
		for (k in (1:length(alleles))) {
       		a1 <-  substr(alleles[k], 4, 4)
        	a2 <-  substr(alleles[k], 9, 9)
        	a.lookup <- rbind(a.lookup, c(a1, a2))
		}
		for (row.nr in (1:dim(group.bin$group.df)[1]) ) {
        	haplo.name <- ""
        	for (i in (1:nr.snps)) {
                val <- a.lookup[i, as.numeric(group.bin$group.df[row.nr, i]) ]
                if (i < nr.snps) {
                        haplo.name <- paste(haplo.name, val, "/", sep="")
                }
                else {
                        haplo.name <- paste(haplo.name, val, sep="")
                }
	        }
        	controls.freq <- group.bin$group.df[row.nr, nr.snps+2]
    	    cases.freq <- group.bin$group.df[row.nr, nr.snps+3]
        	output.frame <- rbind(output.frame, data.frame(
                haplo.name=haplo.name,
                controls.freq=round(controls.freq,2),
                cases.freq=round(cases.freq,2),
                stringsAsFactors = F) )
		}

	}
	
	write.table(output.frame, outfile.name,sep="\t",row.names=F, col.names=T, quote=F)

	sink(outfile.name, append=T)
	print(summary(fit.bin))
	print(paste("Overall haplotype unadjusted P value:", unadj.p))
	print(paste("Overall haplotype adjusted P value:", adj.p))
	sink()
}






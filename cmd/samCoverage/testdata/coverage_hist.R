library(ggplot2)
library(ggpubr)
library(fitdistrplus)
df <- read.table("/Users/nikitha/kellisUrop/software/gonomics/cmd/samCoverage/testdata/test.txt",header=TRUE,sep="\t")
ggline(df, x = "Coverage", y = "Pileups", numeric.x.axis = TRUE, color = "Group", plot_type = "l", xlab = "Coverage", ylab = "Pileups")
freqTable <- rep(df$Coverage, times = df$Pileups)
totalPileups <- sum(subset(df, Group == "Empirical")$Pileups)

########################  FIT NEGATIVE BINOMIAL DISTRIBUTION TO DATA  ############################

negbinom.params <- fitdistr(freqTable,"negative binomial", method = "SANN")$estimate
negbinom <- totalPileups * dnbinom(x = 0:60, size = negbinom.params[1], mu = negbinom.params[2], log = FALSE)
coverage <- seq_along(negbinom) - 1
rfit_df <- data.frame(Coverage = coverage, Pileups = negbinom)
rfit_df$Group <- "rFit"
merged_df <- rbind(df, rfit_df)
ggline(merged_df, x = "Coverage", y = "Pileups", numeric.x.axis = TRUE, color = "Group", plot_type = "l", xlab = "Coverage", ylab = "Pileups")

########################  FIT POISSON DISTRIBUTION TO DATA  #############################

poisson.param <- fitdistr(freqTable,"Poisson", method = "SANN")$estimate
poissondistr = totalPileups * dpois(x = 0:60, lambda = poisson.param)
coverage <- seq_along(poissondistr) - 1
p_rfit_df <- data.frame(Coverage = coverage, Pileups = poissondistr)
p_rfit_df$Group <- "rFitPois"
p_merged_df <- rbind(df, p_rfit_df)
ggline(p_merged_df, x = "Coverage", y = "Pileups", numeric.x.axis = TRUE, color = "Group", plot_type = "l", xlab = "Coverage", ylab = "Pileups")



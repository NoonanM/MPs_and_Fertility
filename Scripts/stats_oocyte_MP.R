setwd("~/Desktop/R/Oocyte + MP")

library(reshape2)
library(nlme)
library(multcomp)
library(viridis)
library(ggplot2)
library(gridExtra)
library(survival)
library(survminer)
library(lme4)
library(MASS)
library("dplyr")
library("ggpubr")

#data input
data <- read.csv("all_data_2.csv")
data$Treatment <- as.factor(data$Treatment)
data$Nuclear_stage <- factor(data$Nuclear_stage, ordered = T, levels = c("BZP", "DG", "MT"))
data$replicate <- as.factor(data$Replicate)


head(data)



##############################################
####     Basic Summary Statistics     ########
##############################################

# Means
means <- aggregate(Percentage ~ Nuclear_stage*Treatment, data = data, FUN = "mean")

# SDs

sds <- aggregate(Percentage ~ Nuclear_stage*Treatment, data = data, FUN = "sd")

means$SD <- sds$Percentage
names(means)[3] <- c("Mean")

write.csv(means, file = "Summary_Statistics_oocyte.csv")



#####################
## plot  Oocyte total##
#####################


png(filename="maturation_final.png",
    width = 5, height = 4, units = "in",
    res = 600)




boxplot(Percentage ~ Treatment*Nuclear_stage,
        data = data,
        cex.axis = 0.5,
        col = viridis(3, alpha = 0.7),
        lwd = 0.7,
        ylab = "",
        ylim = c(0,100),
        axes = FALSE,
        frame.plot = TRUE,
        xlab = "")

TIME <- rev(unique(data$Nuclear_stage))
groups <- 3
numbox <- 3
total <- groups * numbox
xpoints <- seq(median(1:numbox),total,numbox)

axis(2, pretty(data$Percentage), las = 1, cex.axis = 0.7)
axis(1, labels = TIME, at = xpoints,las = 1, cex.axis = 0.7)
abline(v = 3.5, lty = "dashed", col = "grey70")
abline(v = 6.5, lty = "dashed", col = "grey70")



stripchart(data$Percentage ~ data$Treatment*data$Nuclear_stage,               # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = "black",           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE,
           cex = 0.3)        # Add it over


title(ylab = expression("% of total oocytes"), line = 1.7)

legend("top", fill = viridis(3), legend = c( "B0.3", "B1.1", "CT"), horiz = T, cex = 0.7)

dev.off()



#######################################
#general linear model on Maturation####
######################################


#Test if the data are under/overdispersed
test1 <- glm(Percentage/100 ~ Treatment + Nuclear_stage, data = data, family=binomial)
test2 <- glm(Percentage/100 ~ Treatment + Nuclear_stage, data = data, family=quasibinomial)

#Return results
summary(test1)
summary(test2)


fit.quasi.binom <- MASS::glmmPQL(Percentage/100 ~ Treatment + Nuclear_stage,
                                 random = ~ 1|Replicate,
                                 family=quasibinomial,
                                 data = data)

#Return results
summary(fit.quasi.binom)


# Tukey post-hoc test
TUKEY <- glht(fit.quasi.binom, mcp(Treatment = "Tukey"), data = data)
summary(TUKEY)



# Differences at certain maturation stages

####MT####
fit.T1 <- MASS::glmmPQL(Percentage/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Nuclear_stage == "MT"),])

summary(fit.T1)

# Tukey post-hoc test
TUKEY1 <- glht(fit.T1,
              mcp(Treatment = "Tukey"),
              data = data[which(data$Nuclear_stage == "MT"),])
summary(TUKEY1)

###BZP###
fit.T2 <- MASS::glmmPQL(Percentage/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Nuclear_stage == "BZP"),])

summary(fit.T2)

# Tukey post-hoc test
TUKEY2 <- glht(fit.T2,
               mcp(Treatment = "Tukey"),
               data = data[which(data$Nuclear_stage == "BZP"),])
summary(TUKEY2)



###DG###
fit.T3 <- MASS::glmmPQL(Percentage/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Nuclear_stage == "DG"),])

summary(fit.T3)

# Tukey post-hoc test
TUKEY3 <- glht(fit.T3,
               mcp(Treatment = "Tukey"),
               data = data[which(data$Nuclear_stage == "DG"),])
summary(TUKEY3)


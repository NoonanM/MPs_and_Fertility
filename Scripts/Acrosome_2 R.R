setwd("~/Desktop/R/Sperm + MP")

library(reshape2)
library(nlme)
library(multcomp)
library(viridis)
library(ggplot2)
library(gridExtra)
library(survival)
library(survminer)
library(lme4)
library(MuMIn)
library(MASS)

#data input
data <- read.csv("Acrosome.csv") 
data$Treatment <- as.factor(data$Treatment)
data$Replicate <- as.factor(data$Replicate)



##############################################
####     Basic Summary Statistics     ########
##############################################

# Means
means <- aggregate(acrosome ~ Time*Treatment, data = data, FUN = "mean")

# SDs

sds <- aggregate(acrosome ~ Time*Treatment, data = data, FUN = "sd")

means$SD <- sds$acrosome
names(means)[3] <- c("Mean")

write.csv(means, file = "Summary_Statistics_acrosome.csv")



#####################
## plot  Acrosome total##
#####################


png(filename="acrosome.png",
    width = 6.86, height = 4, units = "in",
    res = 600)

par(mgp=c(1.5,0.5,0),
    mar = c(4, #bottom
            4, #left
            3, #top
            2),
    cex.main = 0.8,
    family = "sans" #right
)


boxplot(acrosome ~ Treatment*Time,
        data = data,
        cex.axis = 0.5,
        col = viridis(5, alpha = 0.7),
        lwd = 0.7,
        ylab = "",
        #ylab = expression(Area ~ (mu * m^2)), # "Area (um2)"
        ylim = c(0,100),
        axes = FALSE,
        frame.plot = TRUE,
        xlab = "Time (h)")

TIME <- unique(data$Time)
groups <- 5
numbox <- 5
total <- groups * numbox
xpoints <- seq(median(1:numbox),total,numbox)

axis(2, pretty(data$acrosome), las = 1, cex.axis = 0.7)
axis(1, labels = TIME, at = xpoints,las = 1, cex.axis = 0.7)
abline(v = 5.5, lty = "dashed", col = "grey70")
abline(v = 10.5, lty = "dashed", col = "grey70")
abline(v = 15.5, lty = "dashed", col = "grey70")
abline(v = 20.5, lty = "dashed", col = "grey70")
#abline(v = 25.5, lty = "dashed", col = "grey70")

stripchart(data$acrosome ~ data$Treatment*data$Time,               # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = "black",           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE,
           cex = 0.3)        # Add it over

#title("All treatments")
title(ylab = expression("Intact acrosome (%)"), line = 1.5)

legend("top", fill = viridis(5), legend = c("B0.05", "B0.1", "B0.3", "B1.1", "CT"), horiz = T, cex = 0.5)

dev.off()



############################
#general linear model on acrosome
############################


#Test if the data are under/overdispersed
test1 <- glm(acrosome/100 ~ Treatment + Time, data = data, family=binomial)
test2 <- glm(acrosome/100 ~ Treatment + Time, data = data, family=quasibinomial)

#Return results
summary(test1)
summary(test2)


fit.quasi.binom <- MASS::glmmPQL(acrosome/100 ~ Treatment + Time,
                                 random = ~ 1|Replicate,
                                 family=quasibinomial,
                                 data = data)

#Return results
summary(fit.quasi.binom)


# Tukey post-hoc test
TUKEY <- glht(fit.quasi.binom, mcp(Treatment = "Tukey"), data = data)
summary(TUKEY)



# Differences at certain timepoints
fit.T0.5 <- MASS::glmmPQL(acrosome/100 ~ Treatment,
                                 random = ~ 1|Replicate,
                                 family=quasibinomial,
                                 data = data[which(data$Time == 0.5),])

summary(fit.T0.5)

# Tukey post-hoc test
TUKEY0.5 <- glht(fit.T0.5,
              mcp(Treatment = "Tukey"),
              data = data[which(data$Time == 0.5),])
summary(TUKEY0.5)



# Differences at certain timepoints
fit.T1 <- MASS::glmmPQL(acrosome/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Time == 1.0),])

summary(fit.T1)

# Tukey post-hoc test
TUKEY1 <- glht(fit.T1,
              mcp(Treatment = "Tukey"),
              data = data[which(data$Time == 1.0),])
summary(TUKEY1)


# Differences at certain timepoints
fit.T1.5 <- MASS::glmmPQL(acrosome/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Time == 1.5),])

summary(fit.T1.5)

# Tukey post-hoc test
TUKEY1.5 <- glht(fit.T1.5,
               mcp(Treatment = "Tukey"),
               data = data[which(data$Time == 1.5),])
summary(TUKEY1.5)


# Differences at certain timepoints
fit.T2 <- MASS::glmmPQL(acrosome/100 ~ Treatment,
                          random = ~ 1|Replicate,
                          family=quasibinomial,
                          data = data[which(data$Time == 2.0),])

summary(fit.T2)

# Tukey post-hoc test
TUKEY2 <- glht(fit.T2,
                 mcp(Treatment = "Tukey"),
                 data = data[which(data$Time == 2.0),])
summary(TUKEY2)


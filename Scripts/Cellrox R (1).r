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
data <- read.csv("CR_3.0.csv") 
#data$treatment <- as.factor(data$Treatment, ordered = T, levels = c("CT+", "CT", "b0.05", "b0.1", "b0.3", "b1.1"))
data$Treatment <- as.factor(data$Treatment)
data$Replicate <- as.factor(data$Replicate)



##############################################
####     Basic Summary Statistics     ########
##############################################

# Means
means <- aggregate(CR ~ Time*Treatment, data = data, FUN = "mean")

# SDs

sds <- aggregate(CR ~ Time*Treatment, data = data, FUN = "sd")

means$SD <- sds$CR
names(means)[3] <- c("Mean")

write.csv(means, file = "Summary_Statistics_cellrox.csv")



#####################
## plot  Cellrox total considering time##
#####################


png(filename="cellrox2.png",
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


boxplot(CR ~ Treatment*Time,
        data = data,
        cex.axis = 0.5,
        col = c(viridis(4, alpha = 0.7), "grey75", "grey50"),
        lwd = 0.7,
        ylab = "",
        #ylab = expression(Area ~ (mu * m^2)), # "Area (um2)"
        ylim = c(0,75),
        axes = FALSE,
        frame.plot = TRUE,
        xlab = "Time (h)")

TIME <- unique(data$Time)
groups <- 5
numbox <- 6
total <- groups * numbox
xpoints <- seq(median(1:numbox),total,numbox)

axis(2, pretty(data$CR), las = 1, cex.axis = 0.7)
axis(1, labels = TIME, at = xpoints,las = 1, cex.axis = 0.7)
abline(v = 6.5, lty = "dashed", col = "grey70")
abline(v = 12.5, lty = "dashed", col = "grey70")
abline(v = 18.5, lty = "dashed", col = "grey70")
abline(v = 24.5, lty = "dashed", col = "grey70")

stripchart(data$CR ~ data$Treatment*data$Time,               # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = "black",           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE,
           cex = 0.2)        # Add it over

#title("All treatments")
title(ylab = expression("Cellrox + (%)"), line = 1.5)

legend("top", fill = viridis(6), legend = c("B0.05", "B0.1", "B0.3", "B1.1", "CT", "CT+"), horiz = T, cex = 0.5)


dev.off()




#####################
## plot  Cellrox total no time##
#####################



png(filename="cellrox2notime.png",
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


boxplot(CR ~ Treatment,
        data = data,
        cex.axis = 0.5,
        col = viridis(6, alpha = 0.7),
        lwd = 0.7,
        ylab = "",
        #ylab = expression(Area ~ (mu * m^2)), # "Area (um2)"
        ylim = c(0,100),
        axes = FALSE,
        frame.plot = TRUE,
        xlab = )

stripchart(data$CR ~ data$Treatment,               # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = "black",           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE,
           cex = 0.2)        # Add it over

#title("All treatments")
title(ylab = expression("Cellrox + (%)"), line = 1.5)

legend("top", fill = viridis(6), legend = c("B0.05", "B0.1", "B0.3", "B1.1", "CT", "CT+"), horiz = T, cex = 0.5)

dev.off()



############################
#general linear model on cellrox
############################


#Test if the data are under/overdispersed
test1 <- glm(CR/100 ~ Treatment + Time, data = data, family=binomial)
test2 <- glm(CR/100 ~ Treatment + Time, data = data, family=quasibinomial)

#Return results
summary(test1)
summary(test2)


fit.quasi.binom <- MASS::glmmPQL(CR/100 ~ Treatment + Time,
                                 random = ~ 1|Replicate,
                                 family=quasibinomial,
                                 data = data)

#Return results
summary(fit.quasi.binom)


# Tukey post-hoc test
TUKEY <- glht(fit.quasi.binom, mcp(Treatment = "Tukey"), data = data)
summary(TUKEY)



# Differences at certain timepoints
fit.T0.5 <- MASS::glmmPQL(CR/100 ~ Treatment,
                                 random = ~ 1|Replicate,
                                 family=quasibinomial,
                                 data = data[which(data$Time == 0.5),])

summary(fit.T0.5)

# Tukey post-hoc test
TUKEY <- glht(fit.T1,
              mcp(Treatment = "Tukey"),
              data = data[which(data$Time == 0.5),])
summary(TUKEY)


# Differences at certain timepoints
fit.T1 <- MASS::glmmPQL(CR/100 ~ Treatment,
                        random = ~ 1|Replicate,
                        family=quasibinomial,
                        data = data[which(data$Time == 1.0),])

summary(fit.T1)

# Tukey post-hoc test
TUKEY <- glht(fit.T1,
              mcp(Treatment = "Tukey"),
              data = data[which(data$Time == 0.5),])
summary(TUKEY)


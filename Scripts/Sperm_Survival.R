library(survival)
library(ggplot2)
library(survminer)
library(coxme)
library(cmprsk)
library(viridis)
library(gridExtra)

# Import the sperm motility data
data <- read.csv("Data/Motility.csv")

#Some data carpentry to get the data in the correct format for survival analysis
data$Dead <- data$Total - data$Motile
data$Proportion <- data$Motile/data$Total

#Convert the times from time points to minutes
data[which(data$Time == 0.5),"Time"] <- 30
data[which(data$Time == 1),"Time"] <- 60
data[which(data$Time == 1.5),"Time"] <- 90
data[which(data$Time == 2),"Time"] <- 120

#Rename the groups
data[which(data$Treatment == "CT"),"Treatment"] <- "Control"
data[which(data$Treatment == "B0.3"),"Treatment"] <- "MP_0.3"
data[which(data$Treatment == "B1.1"),"Treatment"] <- "MP_1.1"
data[which(data$Treatment == "B0.05"),"Treatment"] <- "MP_0.05"
data[which(data$Treatment == "B0.1"),"Treatment"] <- "MP_0.1"


DATA <- list()
for(i in 1:nrow(data)){
  
  #Create a new dataframe with an entry for each sperm that was counted
  data_i <- data.frame(Replicate = rep(data[i,"Replicate"], data[i,"Total"]),
                       bull = rep(data[i,"bull"], data[i,"Total"]),
                       Treatment = rep(data[i,"Treatment"], data[i,"Total"]),
                       Time = rep(data[i,"Time"], data[i,"Total"]),
                       Status = rep(0, data[i,"Total"]))
  
  # Set up which sperm should be marked as dead
  START <- data[i,"Motile"]+1
  END <- as.numeric(data[i,"Total"])
  
  data_i[START:END,"Status"] <- 1
  
  #Save the results in the list
  DATA[[i]] <- data_i
}

# convert from a list to a data frame
DATA <- do.call(rbind, DATA)


#----------------------------------------------------------------------
# Survival analysis
#----------------------------------------------------------------------

# Function for extracting information from the cox model
extract_coxme_table <- function (mod){
  beta <- beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

# Generate the Surv object
surv_data <- Surv(DATA$Time, DATA$Status, type = "right")

# Fit a mixed effects cox model to the data
test <- coxme(surv_data ~ Treatment + (1|Replicate) + (1|bull), data = DATA)
summary(test)


#Table of the results
results <- data.frame(HR = round(exp(coef(test)),2),
                      HR_min = round(exp(confint(test)),2)[,1],
                      HR_max = round(exp(confint(test)),2)[,2],
                      P = extract_coxme_table(test)[,4])

#write.csv(results, file = "Results/Motility_Results.csv")




# Fit a survival model (for plotting purposes only)
# Not used to obtain significance as it doesn't account for the effects of the
# Bulls and Replicates
f1 <- survfit(surv_data ~ Treatment, data = DATA)

png(filename="Figures/Survival_Fig.png", width = 2.25, height = 2, units = 'in', res=600)     

B <- 
  ggsurvplot(f1,
             conf.int = T,
             palette = c("black", viridis(4)),
             size = 0.3,
             xlab = "Time (min)", 
             ylab = "Probability of motility",
             legend.title = "",
             legend.labs = c("Control",
                             "0.05 \U03BCm",
                             "0.1 \U03BCm",
                             "0.3 \U03BCm",
                             "1.1 \U03BCm"),
             
             ggtheme = theme_bw(base_size=6, base_family = "sans"),
             font.family = "sans",
             censor = F
  )

B$plot <- B$plot + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = c(0.15,0.2),
        legend.title = element_text(vjust = 5, size = 4, family = "sans", face = "bold"),
        legend.text = element_text(size = 3, family = "sans", face = "bold"),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.y = unit(-0.1, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))

B$plot <- B$plot + 
  ggtitle("B")

B

dev.off()

data$Time2 <- as.factor(data$Time)

A <- 
ggplot(data = data, 
       aes(x = Time2,
               y = Proportion,
               fill = Treatment)) +
  ggtitle("A") +
  geom_point(aes(col = Treatment),
             position = position_jitterdodge(jitter.width = 0.05),
             size = 0.1) +
  geom_boxplot(
               alpha = 0.5,
               size = 0.2,
               outlier.size = 0) +
  scale_fill_manual(labels = c("Control",
                               "0.05 \U03BCm",
                               "0.1 \U03BCm",
                               "0.3 \U03BCm",
                               "1.1 \U03BCm"),
                    values = c("black", viridis(4))) +
  scale_colour_manual(labels = c("Control",
                               "0.05 \U03BCm",
                               "0.1 \U03BCm",
                               "0.3 \U03BCm",
                               "1.1 \U03BCm"),
                    values = c("black", viridis(4))) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=5, family = "sans", face = "bold"),
        axis.title.x = element_text(size=5, family = "sans", face = "bold"),
        axis.text.y = element_text(size=4, family = "sans"),
        axis.text.x  = element_text(size=4, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 6, family = "sans", face = "bold"),
        legend.position = c(0.1,0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 3, family = "sans", face = "bold"),
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.15, "cm"),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.y = unit(-0.1, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  
  ylab(expression(bold(Proportion~motile)))+
  xlab(expression(bold(Time~(min))))



FIG <-
  grid.arrange(A,
               B$plot,
               ncol=2,
               nrow=1)

#Save the figures
ggsave(FIG,
       width = 5.75, height = 2, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Figures/Motility_Fig.png")


# Quick test of the cumulative incidence rate
# Not used for reporting, just as an additional check
test2 <- cuminc(DATA$Time, DATA$Status, DATA$Treatment)
test2
plot(test2,
     xlab = "Time (min)",
     ylab = "Probability of death")

library(lme4)
library(ggResidpanel)
library(dplyr)
library(ggplot2)
library(emmeans)
library(arm)
library(tidyverse)
library(MASS)
library(bestNormalize)


setwd("C:/Users/mfakude/OneDrive - Iowa State University/Desktop/HFF_4")

#Read in HFF data

data1 <- read.csv("./HFF4.csv")
data1$PI = as.factor(data1$PI)
data1$Rep = as.factor(data1$Rep)
data1$HFF = as.numeric(data1$HFF)


hist(data1$HFF)
hist(data1$HFF, col = "Cyan", xlab = "Haploid Female Fertility", main = "Trait distribution")
Boxplot(data1$HFF, col = "Cyan", ylab = "Haploid Female Fertility", main = "Haploid Female Fertility")


# Check for assumptions: P-Value less than 0.05 violates assumption. Can be seen with histogram
shapiro.test(data1$HFF) 



# the code below will find a transformation that works for your data
bestNormalize::bestNormalize(data1$HFF) # read its output
tmp <- bestNormalize::bestNormalize(data1$HFF)
shapiro.test(tmp$chosen_transform$x.t) # best transformation
data1$transformed_data <- tmp$chosen_transform$x.t # adding transformed data
hist(data1$transformed_data) # better?
hist(data1$transformed_data, col = "Cyan", xlab = "Haploid Female Fertility", main = "Trait distribution")

# ANOVA model
str(data1)
table(data1$PI, data1$Rep)
model <- lm(transformed_data~Rep+PI, data = data1)
#model <- lmer(HFF ~ Genotype + (1 | Rep), data = data1) 

ggResidpanel::resid_panel(model)
anova(model) 

#Calculate Entry-mean heritability from ANOVA output
#PI_variance = Genetic variance
#Genetic variance= (PI MSquare - Residual MSquare)/2
heritability <- PI_variance / (PI_variance + (residual_variance / n))

###################################################################### 

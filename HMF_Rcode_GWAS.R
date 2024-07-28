library(lme4)
library(ggResidpanel)
library(dplyr)
library(ggplot2)
library(emmeans)
library(tidyverse)
library("car")
library(qqman)


setwd("C:/Users/mfakude/OneDrive - Iowa State University/Desktop/HMF_L")

#Read in HMF data
da <- read.csv("./HMF_LAST.csv")
#Ensure other variables are categorical(#Transform into Factor)
da$PI = as.factor(da$PI)
da$Rep = as.factor(da$Rep)
da$HMF = as.numeric(da$HMF)


str(da)
#Check for normality
hist(da$HMF) #See distribution
#visualize your data
hist(da$HMF, col = "Cyan", xlab = "Haploid male fertility", main = "Trait distribution")
Boxplot(da$HMF, col = "Cyan", xlab = "Haploid male fertility", main = "Haploid male fertility")
qqnorm(da$HMF, pch = 1, frame = FALSE)
qqline(da$HMF, col = "steelblue", lwd = 2)
qqPlot(da$HMF)

shapiro.test(da$HMF) #violates assumption 

#da = read_xlsx("HMF_raw.xlsx",sheet = 1)
da.agg = da %>% group_by(PI, Rep) %>% summarize_at(.vars = c("Fertile","Sterile"),.funs = sum) %>% 
  mutate(HMF = Fertile/(Fertile + Sterile))
da.agg$HMF.mod = da.agg$HMF
da.agg$HMF.mod[da.agg$HMF == 0] = min(da.agg$HMF[da.agg$HMF > 0])/2
#Emmean inverse logit from glm
glm_GWAS = glm(cbind(Fertile, Sterile) ~ Rep + PI, da = da.agg, family = binomial)
resid_panel(glm_GWAS)
glm_eblues_GWAS <- data.frame(emmeans(glm_GWAS, "PI"))
glm_eblues_GWAS$inverse_logit <- arm::invlogit(glm_eblues_GWAS$emmean)*100
view(glm_eblues_GWAS)
summary(glm_eblues_GWAS$inverse_logit)
view(glm_eblues_GWAS$inverse_logit)
write.csv(x = , glm_eblues_GWAS, file = "HMF_BLUES.csv",row.names = F,quote = F)#writes the csv file with BLUES

############### Run BS39 and BS39+SHGD HMF GWAS############################################################
Pheno.data = read.csv("./HMF_BLUES.csv",header = TRUE)
Geno.data = read.delim("./filtered_genotype_data1.txt", sep ='\t',header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#GWAS
myGAPIT=GAPIT(
  Y=Pheno.data[,c(1,2)], #fist column is ID. Currently using one trait. HMF
  G=Geno.data,
  PCA.total=2,
  model=c("FarmCPU", "MLMM"),
  Multiple_analysis=TRUE,SNP.MAF = 0.01,
  file.output=T)


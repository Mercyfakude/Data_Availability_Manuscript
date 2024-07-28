
#For HMF Entry-mean heritability
library(lme4)

setwd("C:/Users/mfakude/OneDrive - Iowa State University/Desktop/HMF_L")

#Read in HFF data
da <- read.csv("./HMF_LAST.csv")

# Install and load necessary package
install.packages("lme4")
library(lme4)

# Fit a binomial GLMM
glmm_model <- glmer(cbind(Fertile, Sterile) ~ (1 | PI) + (1 | Rep), 
                    family = binomial, data = da.agg)

# Print the summary of the model
summary(glmm_model)
resid_panel(glmm_model)

# Extract variance components
varcomp <- as.data.frame(VarCorr(glmm_model))

# Variance for PI
pi_variance <- varcomp[varcomp$grp == "PI", "vcov"]

# Variance for Rep
rep_variance <- varcomp[varcomp$grp == "Rep", "vcov"]

# Calculate residual variance
deviance_residuals <- residuals(glmm_model, type = "deviance")
residual_variance <- var(deviance_residuals)

# Print variance components
cat("Variance for PI:", pi_variance, "\n")
cat("Variance for Rep:", rep_variance, "\n")
cat("Residual Variance:", residual_variance, "\n")


# Variance components 
pi_variance <- 2.172924
residual_variance <- 4.572366
rep_variance <- 0.03924474

# Number of environments 
n <- 2

# Calculate entry-mean heritability
heritability <- pi_variance / (pi_variance + (residual_variance / n))

# Print the heritability
cat("Entry-Mean Heritability (H²):", heritability, "\n")


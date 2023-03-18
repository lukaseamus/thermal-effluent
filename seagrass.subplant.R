###########################################################################
## Project: Effect of thermal effluent on Posidonia oceanica communities ##
## Script purpose: Analysis of detailed seagrass data                    ##
## Author: Luka Seamus Wright                                            ##
###########################################################################

#### 1.    Data exploration ####
#### 1.1   Load data ####
seagrass <- read.csv("~/Desktop/University of Malta/Posidonia/Data/seagrass.subplant.csv")
seagrass <- seagrass[1:240,]

#### 1.2   Rename variables ####
site <- seagrass$Site
dista <- seagrass$Distance
stage <- seagrass$Stage

leaves <- seagrass$Leaves # number of leaves per shoot
area <- seagrass$Area # leaf area per shoot (cm2)
biomass <- seagrass$Biomass # leaf biomass per shoot (g)
SLA <- area/biomass # specific leaf area (cm2 g-1)

SLA <- SLA[-c(152, 201, 209)]
ssite <- site[-c(152, 201, 209)]
sdist <- dista[-c(152, 201, 209)]
sstage <- stage[-c(152, 201, 209)]

#### 1.3   Test for assumptions ####
#### 1.3.1 Normality ####
par(mfrow = c(2,1), mar = c(2,2,2,1)) # modify plot environment

hist(leaves)
boxplot(leaves, horizontal = T) # somewhat right-skewed

hist(area)
boxplot(area, horizontal = T) # right-skewed

hist(biomass)
boxplot(biomass, horizontal = T) # right-skewed

hist(SLA)
boxplot(SLA, horizontal = T) # right-skewed

#### 1.3.2 Homogeneity ####
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1) # return to default

boxplot(leaves ~ site * dista * stage) # slightly heterogenous

boxplot(area ~ site * dista * stage) # heterogenous

boxplot(biomass ~ site * dista * stage) # heterogenous

boxplot(SLA ~ ssite * sdist * sstage) # heterogenous

#### 2.    Data anaylsis ####
#### 2.1   Number of leaves ####
#### 2.1.1 Build model ####
m1 <- lm(leaves ~ site * dista * stage)

#### 2.1.2 Determine fixed components ####
m2 <- lm(leaves ~ site * stage) # remove interaction term and distance variable
anova(m1, m2) # models are not significantly different
# -> proceed with the simpler model: m2

#### 2.1.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m2, pch = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m2)) # residuals are quite normally distributed
qqnorm(resid(m2))
qqline(resid(m2)) # and only deviate slightly at the upper distribution edge

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m2 is chosen as the optimal model

#### 2.1.4 Interpret model ####
require(car)
# type 3 sums of squares for overall p-values
Anova(m2, type = 3)
# Response: leaves
#             Sum Sq  Df   F value    Pr(>F)    
# (Intercept) 429.03   1 1103.8136 < 2.2e-16 ***
# site          1.25   1    3.2161  0.074211 .  
# stage        67.95   2   87.4123 < 2.2e-16 ***
# site:stage    5.03   2    6.4750  0.001832 ** 
# Residuals    90.95 234 

# pairwise contrasts
require(emmeans)
emmeans(m2, pairwise ~ site * stage)$contrasts
# contrast                                 estimate    SE  df t.ratio p.value
# Kbira Adult - Zghira Adult                   0.25 0.139 234  1.793  0.4720 
# Kbira Adult - Kbira Intermediate             1.73 0.139 234 12.374  <.0001 
# Kbira Adult - Zghira Intermediate            1.73 0.139 234 12.374  <.0001 
# Kbira Adult - Kbira Juvenile                 1.43 0.139 234 10.222  <.0001 
# Kbira Adult - Zghira Juvenile                2.12 0.139 234 15.243  <.0001 
# Zghira Adult - Kbira Intermediate            1.48 0.139 234 10.581  <.0001 
# Zghira Adult - Zghira Intermediate           1.48 0.139 234 10.581  <.0001 
# Zghira Adult - Kbira Juvenile                1.18 0.139 234  8.429  <.0001 
# Zghira Adult - Zghira Juvenile               1.88 0.139 234 13.450  <.0001 
# Kbira Intermediate - Zghira Intermediate     0.00 0.139 234  0.000  1.0000 
# Kbira Intermediate - Kbira Juvenile         -0.30 0.139 234 -2.152  0.2646 
# Kbira Intermediate - Zghira Juvenile         0.40 0.139 234  2.869  0.0506 
# Zghira Intermediate - Kbira Juvenile        -0.30 0.139 234 -2.152  0.2646 
# Zghira Intermediate - Zghira Juvenile        0.40 0.139 234  2.869  0.0506 
# Kbira Juvenile - Zghira Juvenile             0.70 0.139 234  5.021  <.0001


#### 2.2   Leaf biomass ####
#### 2.2.1 Build model ####
m3 <- lm(biomass ~ site * dista * stage)

#### 2.2.2 Determine fixed components ####
m4 <- lm(biomass ~ site * stage) # remove interaction term and distance variable
anova(m3, m4) # m3 fits significantly better

#### 2.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, pch = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m3) ~ stage)
plot(resid(m3) ~ dista)
# residual variance is most contrasted between stages
# with nlme::gls() heterogeneity can be modelled for stages

# test assumption of normality
hist(resid(m3)) # residuals seem quite normally distributed
qqnorm(resid(m3))
qqline(resid(m3)) # but deviate at the distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving
# first, let's try improving normality with gamma,
# which may simultaneously fix homogeneity

#### 2.2.4 Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(biomass + 1e-3, "gamma") 
norm <- fitdist(biomass + 1e-3, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits better

#### 2.2.5 Build gamma model ####
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))

#### 2.2.6 Determine fixed components ####
m6 <- glm(biomass + 1e-3 ~ site * stage, family = Gamma(link = "log"))
anova(m5, m6, test = "Chisq") # m5 fits significantly better

#### 2.2.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m5, pch = site)
# homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m5))
qqnorm(resid(m5))
qqline(resid(m5)) # residuals are normal 
# except for one outlier

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# gamma fixed normality and homogeneity
# m5 is chosen as the optimal model

#### 2.2.8 Interpret model ####
# type 3 sums of squares for overall p-values
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.000  1     0.9891    
# dista               1.224  1     0.2686    
# stage             193.258  2     <2e-16 ***
# site:dista          2.283  1     0.1308    
# site:stage          0.605  2     0.7389    
# dista:stage         0.511  2     0.7745    
# site:dista:stage    2.109  2     0.3483  

# pairwise contrasts (only significant contrasts shown)
summary(m5)
# Kbira:Adult vs Kbira:Intermediate, p < 0.001 ***
# Kbira:Adult vs Kbira:Juvenile, p < 0.001 ***

coef(m5)
# y = exp(0.004060946*x - 0.907064503)

stage <- factor(stage, levels = c("Intermediate", "Juvenile", "Adult"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.923  1   0.336776    
# dista               0.509  1   0.475412    
# stage             193.258  2  < 2.2e-16 ***
# site:dista          7.178  1   0.007379 ** 
# site:stage          0.605  2   0.738917    
# dista:stage         0.511  2   0.774519    
# site:dista:stage    2.109  2   0.348308   

summary(m5)
# Kbira:Intermediate vs Kbira:Juvenile, p < 0.001 ***

coef(m5)
# y = exp(0.002626632*x - 3.167281493)

stage <- factor(stage, levels = c("Juvenile", "Intermediate", "Adult"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.874  1     0.3497    
# dista               0.012  1     0.9112    
# stage             193.258  2     <2e-16 ***
# site:dista          0.448  1     0.5034    
# site:stage          0.605  2     0.7389    
# dista:stage         0.511  2     0.7745    
# site:dista:stage    2.109  2     0.3483 

coef(m5)
# y = exp(0.0004045533*x - 4.8147319484)

site <- factor(site, levels = c("Zghira", "Kbira"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.874  1     0.3497    
# dista               0.700  1     0.4027    
# stage              76.364  2     <2e-16 ***
# site:dista          0.448  1     0.5034    
# site:stage          0.605  2     0.7389    
# dista:stage         3.703  2     0.1570    
# site:dista:stage    2.109  2     0.3483  

summary(m5)
# Zghira:Juvenile vs Zghira:Intermediate, p < 0.001 ***
# Zghira:Juvenile vs Zghira:Adult, p < 0.001 ***

coef(m5)
# y = exp(0.005254916*x - 5.192194106)

stage <- factor(stage, levels = c("Intermediate", "Juvenile", "Adult"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.923  1   0.336776    
# dista              12.134  1   0.000495 ***
# stage              76.364  2  < 2.2e-16 ***
# site:dista          7.178  1   0.007379 ** 
# site:stage          0.605  2   0.738917    
# dista:stage         3.703  2   0.156980    
# site:dista:stage    2.109  2   0.348308 

summary(m5)
# Zghira:Intermediate vs Zghira:Adult, p < 0.001 ***

coef(m5)
# y = exp(0.022362431*x - 3.559028296)

stage <- factor(stage, levels = c("Adult", "Intermediate", "Juvenile"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass + 0.001
#                  LR Chisq Df Pr(>Chisq)    
# site                0.000  1    0.98908    
# dista               5.776  1    0.01625 *  
# stage              76.364  2    < 2e-16 ***
# site:dista          2.283  1    0.13079    
# site:stage          0.605  2    0.73892    
# dista:stage         3.703  2    0.15698    
# site:dista:stage    2.109  2    0.34831 

coef(m5)
# y = exp(0.014909765*x - 0.912492835)

site <- factor(site, levels = c("Kbira", "Zghira"))
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))

#### 2.3   Leaf area ####
#### 2.3.1 Build model ####
m7 <- lm(area ~ site * dista * stage)

#### 2.3.2 Determine fixed components ####
m8 <- lm(area ~ site * stage) # remove distance variable
anova(m7, m8) # m7 fits better

#### 2.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m7, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m7) ~ stage)
plot(resid(m7) ~ dista)
# residual variance is most contrasted between stages
# with nlme::gls() heterogeneity can be modelled for stages

# test assumption of normality
hist(resid(m7)) # residuals seem quite normally distributed
qqnorm(resid(m7))
qqline(resid(m7)) # but deviate dramatically at the distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving
# first, let's try improving normality with gamma,
# which may simultaneously fix homogeneity

#### 2.3.4 Fit gamma distribution ####
gamma <- fitdist(area + 1, "gamma") 
norm <- fitdist(area + 1, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits better

#### 2.3.5 Build gamma model ####
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))

#### 2.3.6 Determine fixed components ####
m10 <- glm(area + 1 ~ site * stage, family = Gamma(link = "log"))
anova(m9, m10, test = "Chisq") # m9 fits better

#### 2.3.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m9, col = site)
# homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m9))
qqnorm(resid(m9))
qqline(resid(m9)) # residuals are normal 
# except for one outlier

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# gamma fixed normality and homogeneity
# m9 is chosen as the optimal model

#### 2.3.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                0.201  1    0.65375    
# dista               0.322  1    0.57012    
# stage             146.866  2    < 2e-16 ***
# site:dista          2.724  1    0.09886 .  
# site:stage          0.911  2    0.63409    
# dista:stage         1.055  2    0.59015    
# site:dista:stage    4.339  2    0.11425  

# pairwise contrasts (only significant contrasts shown)
summary(m9)
# Kbira:Adult vs Kbira:Intermediate, p < 0.001 ***
# Kbira:Adult vs Kbira:Juvenile, p < 0.001 ***

coef(m9)
# y = exp(0.001822496*x + 4.451824305)

stage <- factor(stage, levels = c("Intermediate", "Juvenile", "Adult"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                3.110  1  0.0778150 .  
# dista               0.697  1  0.4038283    
# stage             146.866  2  < 2.2e-16 ***
# site:dista         12.166  1  0.0004867 ***
# site:stage          0.911  2  0.6340925    
# dista:stage         1.055  2  0.5901493    
# site:dista:stage    4.339  2  0.1142480  

summary(m9)
# Kbira:Intermediate vs Kbira:Juvenile, p < 0.001 ***

coef(m9)
# y = exp(-0.002708839*x + 2.824796102)

stage <- factor(stage, levels = c("Juvenile", "Intermediate", "Adult"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                1.083  1     0.2981    
# dista               0.207  1     0.6490    
# stage             146.866  2     <2e-16 ***
# site:dista          0.390  1     0.5324    
# site:stage          0.911  2     0.6341    
# dista:stage         1.055  2     0.5901    
# site:dista:stage    4.339  2     0.1142 

coef(m9)
# y = exp(-0.001449520*x + 1.488421157)

site <- factor(site, levels = c("Zghira", "Kbira"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                1.083  1    0.29806    
# dista               0.210  1    0.64701    
# stage              56.208  2  6.233e-13 ***
# site:dista          0.390  1    0.53239    
# site:stage          0.911  2    0.63409    
# dista:stage         4.945  2    0.08439 .  
# site:dista:stage    4.339  2    0.11425  

summary(m9)
# Zghira:Juvenile vs Zghira:Intermediate, p = 0.01 *
# Zghira:Juvenile vs Zghira:Adult, p < 0.001 ***

coef(m9)
# y = exp(0.002514847*x + 1.120449690)

stage <- factor(stage, levels = c("Intermediate", "Juvenile", "Adult"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                3.110  1  0.0778150 .  
# dista              12.548  1  0.0003967 ***
# stage              56.208  2  6.233e-13 ***
# site:dista         12.166  1  0.0004867 ***
# site:stage          0.911  2  0.6340925    
# dista:stage         4.945  2  0.0843919 .  
# site:dista:stage    4.339  2  0.1142480  

summary(m9)
# Zghira:Intermediate vs Zghira:Adult, p < 0.001 ***

coef(m9)
# y = exp(0.019780860*x + 2.193761703)

stage <- factor(stage, levels = c("Adult", "Intermediate", "Juvenile"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
Anova(m9, type = 3)
# Response: area + 1
#                  LR Chisq Df Pr(>Chisq)    
# site                0.201  1    0.65375    
# dista               5.035  1    0.02484 *  
# stage              56.208  2  6.233e-13 ***
# site:dista          2.724  1    0.09886 .  
# site:stage          0.911  2    0.63409    
# dista:stage         4.945  2    0.08439 .  
# site:dista:stage    4.339  2    0.11425  

coef(m9)
# y = exp(0.012248067*x + 4.294896324)

site <- factor(site, levels = c("Kbira", "Zghira"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))



#### 2.4   Specific leaf area ####
#### 2.4.1 Build model ####
m11 <- lm(SLA ~ ssite * sdist * sstage)

#### 2.4.2 Determine fixed components ####
m12 <- lm(SLA ~ ssite * sstage) # remove distance variable
anova(m11, m12) # m12 fits better but retain extra variable for hypothesis test
m12b <- lm(SLA ~ ssite + sstage)
anova(m12, m12b) # m12b fits better but retain interaction for hypothesis test

#### 2.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m11, pch = ssite)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m11) ~ sstage)
plot(resid(m11) ~ sdist)
# residual variance is most contrasted between stages
# with nlme::gls() heterogeneity can be modelled for stages

# test assumption of normality
hist(resid(m11)) # residuals seem quite balanced
qqnorm(resid(m11))
qqline(resid(m11)) # but deviate dramatically at distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving
# first, let's try improving normality with gamma,
# which may simultaneously fix homogeneity

#### 2.4.4 Fit gamma distribution ####
gamma <- fitdist(SLA, "gamma") 
norm <- fitdist(SLA, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits better

#### 2.4.5 Build gamma model ####
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))

#### 2.4.6 Determine fixed components ####
m14 <- glm(SLA ~ ssite * sstage, family = Gamma(link = "log"))
anova(m13, m14, test = "Chisq") # m14 fits better but retain extra variable for hypothesis test
m14b <- glm(SLA ~ ssite + sstage, family = Gamma(link = "log"))
anova(m14, m14b, test = "Chisq") # m14b fits better but retain interaction for hypothesis test

#### 2.4.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m13, col = site)
# more homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m13))
qqnorm(resid(m13))
qqline(resid(m13)) # residuals are somewhat more normal 

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# gamma improved normality and homogeneity
# m13 is chosen as the optimal model

#### 2.4.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)   
# ssite                1.7797  1   0.182181   
# sdist                5.0818  1   0.024178 * 
# sstage              10.3761  2   0.005583 **
# ssite:sdist          0.4653  1   0.495140   
# ssite:sstage         0.8180  2   0.664299   
# sdist:sstage         1.9609  2   0.375136   
# ssite:sdist:sstage   0.5052  2   0.776776 

# pairwise contrasts (only significant contrasts shown)
summary(m13)
# Kbira:Adult vs Kbira:Juvenile, p = 0.002 **

coef(m13)
# y = exp(-0.005700669*x + 5.577933726)

sstage <- factor(sstage, levels = c("Intermediate", "Juvenile", "Adult"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)   
# ssite                0.1336  1   0.714775   
# sdist                3.2493  1   0.071453 . 
# sstage              10.3761  2   0.005583 **
# ssite:sdist          0.2180  1   0.640605   
# ssite:sstage         0.8180  2   0.664299   
# sdist:sstage         1.9609  2   0.375136   
# ssite:sdist:sstage   0.5052  2   0.776776  

summary(m13)
# no significant pairwise contrasts

coef(m13)
# y = exp(-0.004636059*x + 5.921731885)

sstage <- factor(sstage, levels = c("Juvenile", "Intermediate", "Adult"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)   
# ssite                0.0117  1   0.913747   
# sdist                0.1283  1   0.720163   
# sstage              10.3761  2   0.005583 **
# ssite:sdist          0.0815  1   0.775236   
# ssite:sstage         0.8180  2   0.664299   
# sdist:sstage         1.9609  2   0.375136   
# ssite:sdist:sstage   0.5052  2   0.776776  

coef(m13)
# y = exp(-0.0009103437*x + 6.1915341712)

ssite <- factor(ssite, levels = c("Zghira", "Kbira"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)  
# ssite                0.0117  1    0.91375  
# sdist                0.2828  1    0.59487  
# sstage               7.5890  2    0.02249 *
# ssite:sdist          0.0815  1    0.77524  
# ssite:sstage         0.8180  2    0.66430  
# sdist:sstage         0.0002  2    0.99991  
# ssite:sdist:sstage   0.5052  2    0.77678 

summary(m13)
# Zghira:Juvenile vs Zghira:Adult, p = 0.007 **

coef(m13)
# y = exp(-0.002381412*x + 6.160182)

sstage <- factor(sstage, levels = c("Intermediate", "Juvenile", "Adult"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)  
# ssite                0.1336  1    0.71478  
# sdist                0.2956  1    0.58668  
# sstage               7.5890  2    0.02249 *
# ssite:sdist          0.2180  1    0.64061  
# ssite:sstage         0.8180  2    0.66430  
# sdist:sstage         0.0002  2    0.99991  
# ssite:sdist:sstage   0.5052  2    0.77678 

summary(m13)
# no significant pairwise contrasts

coef(m13)
# y = exp(-0.002312997*x + 5.821348)

sstage <- factor(sstage, levels = c("Adult", "Intermediate", "Juvenile"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: SLA
#                    LR Chisq Df Pr(>Chisq)  
# ssite                1.7797  1    0.18218  
# sdist                0.2910  1    0.58955  
# sstage               7.5890  2    0.02249 *
# ssite:sdist          0.4653  1    0.49514  
# ssite:sstage         0.8180  2    0.66430  
# sdist:sstage         0.0002  2    0.99991  
# ssite:sdist:sstage   0.5052  2    0.77678   

coef(m13)
# y = exp(-0.002307157*x + 5.208123)

ssite <- factor(ssite, levels = c("Kbira", "Zghira"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))

# type II sums of squares for overall p-values
Anova(m14b, type = 2)
# Response: SLA
#        LR Chisq Df Pr(>Chisq)    
# ssite     5.226  1    0.02225 *  
# sstage  182.560  2    < 2e-16 ***

summary(m14b)
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.26128    0.05270  99.842  < 2e-16 ***
# ssiteZghira        -0.12122    0.05295  -2.289    0.023 *  
# sstageIntermediate  0.50182    0.06464   7.763  2.6e-13 ***
# sstageJuvenile      0.88783    0.06485  13.690  < 2e-16 ***
  
sstage <- factor(sstage, levels = c("Intermediate", "Juvenile", "Adult"))
m14b <- glm(SLA ~ ssite + sstage, family = Gamma(link = "log"))
summary(m14b)
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     5.76310    0.05278 109.195  < 2e-16 ***
# ssiteZghira    -0.12122    0.05295  -2.289    0.023 *  
# sstageJuvenile  0.38601    0.06505   5.934 1.06e-08 ***
# sstageAdult    -0.50182    0.06464  -7.763 2.60e-13 ***

sstage <- factor(sstage, levels = c("Adult", "Intermediate", "Juvenile"))
m14b <- glm(SLA ~ ssite + sstage, family = Gamma(link = "log"))

#### 3.    Data visualisation ####
#### 3.1   Calculate descriptive statistics ####
leaves.stat <- aggregate(leaves ~ site + stage, FUN = mean)
leaves.stat$se <- aggregate(leaves ~ site + stage, FUN = function(x) sd(x)/sqrt(length(x)))[,3]
colnames(leaves.stat)[3] <- "mean"

biomass.stat <- aggregate(biomass ~ site + dista + stage, FUN = mean)
biomass.stat$se <- aggregate(biomass ~ site + dista + stage, FUN = function(x) sd(x)/sqrt(length(x)))[,4]
colnames(biomass.stat)[4] <- "mean"
biomass.stat <- biomass.stat[rep(seq_len(nrow(biomass.stat)), each = 103), ]
rownames(biomass.stat) <- NULL
biomass.stat <- biomass.stat[-c(820:824), ]
rownames(biomass.stat) <- NULL
biomass.stat <- biomass.stat[-c(1639:1643), ]
rownames(biomass.stat) <- NULL
biomass.stat <- biomass.stat[-c(2458:2462), ]
rownames(biomass.stat) <- NULL

area.stat <- aggregate(area ~ site + dista + stage, FUN = mean)
area.stat$se <- aggregate(area ~ site + dista + stage, FUN = function(x) sd(x)/sqrt(length(x)))[,4]
colnames(area.stat)[4] <- "mean"
area.stat <- area.stat[rep(seq_len(nrow(area.stat)), each = 103), ]
rownames(area.stat) <- NULL
area.stat <- area.stat[-c(820:824), ]
rownames(area.stat) <- NULL
area.stat <- area.stat[-c(1639:1643), ]
rownames(area.stat) <- NULL
area.stat <- area.stat[-c(2458:2462), ]
rownames(area.stat) <- NULL

SLA.stat <- aggregate(SLA ~ ssite + sdist + sstage, FUN = mean)
SLA.stat$se <- aggregate(SLA ~ ssite + sdist + sstage, FUN = function(x) sd(x)/sqrt(length(x)))[,4]
colnames(SLA.stat)[4] <- "mean"
SLA.stat <- SLA.stat[rep(seq_len(nrow(SLA.stat)), each = 103), ]
rownames(SLA.stat) <- NULL
SLA.stat <- SLA.stat[-c(820:824), ]
rownames(SLA.stat) <- NULL
SLA.stat <- SLA.stat[-c(1639:1643), ]
rownames(SLA.stat) <- NULL
SLA.stat <- SLA.stat[-c(2458:2462), ]
rownames(SLA.stat) <- NULL

#### 3.2   Calculate model predictions ####
# re-run models
m5 <- glm(biomass + 1e-3 ~ site * dista * stage, family = Gamma(link = "log"))
m9 <- glm(area + 1 ~ site * dista * stage, family = Gamma(link = "log"))
m13 <- glm(SLA ~ ssite * sdist * sstage, family = Gamma(link = "log"))

# generate new data frames for model predictions to work on
new <- data.frame(site = c(rep("Kbira", 496), rep("Zghira", 323),
                           rep("Kbira", 496), rep("Zghira", 323),
                           rep("Kbira", 496), rep("Zghira", 323)),
                  stage = c(rep("Adult", 819), rep("Intermediate", 819),
                            rep("Juvenile", 819)),
                  sig = c(rep("NS", 496), rep("S", 323), rep("NS", 496),
                          rep("S", 323), rep("NS", 819)),
                  dista = c(seq(21.22, 70.78, by = 0.1),
                            seq(36.35, 68.63, by = 0.1),
                            seq(21.22, 70.78, by = 0.1),
                            seq(36.35, 68.63, by = 0.1),
                            seq(21.22, 70.78, by = 0.1),
                            seq(36.35, 68.63, by = 0.1)))

snew <- data.frame(ssite = c(rep("Kbira", 496), rep("Zghira", 323),
                           rep("Kbira", 496), rep("Zghira", 323),
                           rep("Kbira", 496), rep("Zghira", 323)),
                   sstage = c(rep("Adult", 819), rep("Intermediate", 819),
                              rep("Juvenile", 819)),
                   sig = c(rep("S", 496), rep("NS", 1961)),
                   sdist = c(seq(21.22, 70.78, by = 0.1),
                             seq(36.35, 68.63, by = 0.1),
                             seq(21.22, 70.78, by = 0.1),
                             seq(36.35, 68.63, by = 0.1),
                             seq(21.22, 70.78, by = 0.1),
                             seq(36.35, 68.63, by = 0.1)))

# calculate predicted values and intervals
inv <- family(m5)$linkinv # calculate inverse of link function
pred.b <- data.frame(predict(m5, type = "link", newdata = new,
                             se.fit = TRUE))
new$b.fit <- inv(pred.b$fit) # add line fit to dataset
new$b.hi <- inv(pred.b$fit + pred.b$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$b.lo <- inv(pred.b$fit - pred.b$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m9)$linkinv # calculate inverse of link function
pred.a <- data.frame(predict(m9, type = "link", newdata = new,
                             se.fit = TRUE)) 
new$a.fit <- inv(pred.a$fit) # add line fit to dataset
new$a.hi <- inv(pred.a$fit + pred.a$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$a.lo <- inv(pred.a$fit - pred.a$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m13)$linkinv # calculate inverse of link function
pred.sla <- data.frame(predict(m13, type = "link", newdata = snew,
                               se.fit = TRUE)) 
snew$sla.fit <- inv(pred.sla$fit) # add line fit to dataset
snew$sla.hi <- inv(pred.sla$fit + pred.sla$se.fit * qnorm(0.975)) # upper confidence interval boundary
snew$sla.lo <- inv(pred.sla$fit - pred.sla$se.fit * qnorm(0.975)) # lower confidence interval boundary

#### 3.3   Define theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 3.4   Plot ####
leaves.stat$site <- factor(leaves.stat$site, levels = c("Zghira", "Kbira"))

lp <- ggplot(leaves.stat, aes(site, mean)) +
        geom_col(aes(fill = site, alpha = stage), position = "dodge", width = 0.7) +
        # geom_point(data = seagrass, aes(Site, Leaves, colour = Stage), 
        #            position = position_dodge(width = 0.7)) +
        geom_errorbar(aes(ymin = mean - se*qnorm(0.975), 
                          ymax = mean + se*qnorm(0.975), 
                          colour = stage),
                      width = .1, lwd = .4, position = position_dodge(width = 0.7)) +
        geom_text(aes(x = site, y = mean + se*qnorm(0.975) + 0.1,
                      label = c("a","a","bc","bc","b","c")),
                      nudge_x = c(-0.233,-0.233,0,0,0.233,0.233)) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          guide = "none") +
        scale_colour_manual(values = rep("#000000", 3),
                            guide = "none") +
        scale_alpha_manual(values = c(0.3, 0.6, 1)) +
        ylab(expression("Number of leaves (shoot"^-1*")")) +
        scale_x_discrete(labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira")) +
        coord_cartesian(ylim = c(0, 5)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = c(.21, .92),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()) +
        mytheme
lp # dimensions: 4.5 x 4 in

new$site <- factor(new$site, levels = c("Zghira", "Kbira"))
colnames(seagrass)[4] <- "stage"

bp <- ggplot() +
        geom_line(data = new, aes(dista, b.fit - 1e-3, colour = site, lty = sig), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = b.lo - 1e-3, ymax = b.hi - 1e-3, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Biomass, colour = Site),
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = biomass.stat, aes(dista, mean, ymin = mean - se,
                                                 ymax = mean + se, fill = site),
                        size = 0.5, shape = 21, 
                        colour = rep(c(rep("#417892",206),rep("#ee850d",412),rep("#417892",201)),3)) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = "none") +
        ylab(expression("Leaf biomass (g shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 2.5), xlim = c(10, 80)) +
        scale_x_continuous(breaks = seq(10, 80, by = 10), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), expand = c(0, 0)) +
        facet_grid(~stage) +
        theme(legend.position = c(.13, .9),
              strip.background = element_rect(colour = "black", size = 1,
                                              fill = NA),
              strip.text = element_text(size = 15, hjust = 0),
              panel.spacing = unit(1, "cm")) +
        mytheme
bp # dimensions: 4.5 x 8 in

ap <- ggplot() +
        geom_line(data = new, aes(dista, a.fit - 1, colour = site, lty = sig), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = a.lo - 1, ymax = a.hi - 1, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Area, colour = Site),
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = area.stat, aes(dista, mean, ymin = mean - se,
                                              ymax = mean + se, fill = site),
                        size = 0.5, shape = 21,
                        colour = rep(c(rep("#417892",206),rep("#ee850d",412),rep("#417892",201)),3)) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = "none") +
        ylab(expression("Leaf area (cm"^2*" shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 300), xlim = c(10, 80)) +
        scale_x_continuous(breaks = seq(10, 80, by = 10), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 300, by = 50), expand = c(0, 0)) +
        facet_grid(~stage) +
        theme(legend.position = c(.13, .9),
              strip.background = element_rect(colour = "black", size = 1,
                                              fill = NA),
              strip.text = element_text(size = 15, hjust = 0),
              panel.spacing = unit(1, "cm")) +
        mytheme
ap # dimensions: 4.5 x 8 in

snew$ssite <- factor(snew$ssite, levels = c("Zghira", "Kbira"))
colnames(seagrass)[4] <- "sstage"

slap <- ggplot() +
        geom_line(data = snew, aes(sdist, sla.fit, colour = ssite, lty = sig), size = 0.5) +
        geom_ribbon(data = snew, aes(sdist, ymin = sla.lo, ymax = sla.hi, fill = ssite),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Area/Biomass, colour = Site),
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = SLA.stat, aes(sdist, mean, ymin = mean - se,
                                             ymax = mean + se, fill = ssite),
                        size = 0.5, shape = 21,
                        colour = rep(c(rep("#417892",206),rep("#ee850d",412),rep("#417892",201)),3)) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = "none") +
        ylab(expression("Specific leaf area (cm"^2*" g"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 1400), xlim = c(10, 80)) +
        scale_x_continuous(breaks = seq(10, 80, by = 10), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 1400, by = 200), expand = c(0, 0)) +
        facet_grid(~sstage) +
        theme(legend.position = c(.13, .9),
              strip.background = element_rect(colour = "black", size = 1,
                                              fill = NA),
              strip.text = element_text(size = 15, hjust = 0),
              panel.spacing = unit(1, "cm")) +
        mytheme
slap # dimensions: 4.5 x 8 in

#### 3.5   Combine plots ####
require(cowplot)
absp <- plot_grid(ap + theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank()), 
                  bp + theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank(),
                             legend.position = "none",
                             strip.text = element_blank(),
                             strip.background = element_blank()), 
                  slap + theme(legend.position = "none",
                               strip.text = element_blank(),
                               strip.background = element_blank()), align = "v", rel_heights = c(1, 0.91, 1.02),
                  ncol = 1, labels = "auto", label_size = 15, label_fontfamily = "Helvetica Neue")
absp # dimensions: 10 x 8 in

#### 4.    Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")

###########################################################################
## Project: Effect of thermal effluent on Posidonia oceanica communities ##
## Script purpose: Analysis of seagrass data                             ##
## Author: Luka Seamus Wright                                            ##
###########################################################################

#### 1.    Data exploration ####
#### 1.1   Load data ####
seagrass <- read.csv("~/Desktop/University of Malta/Posidonia/Data/seagrass.plant.csv")
epiphytes <- seagrass[81:160,]
seagrass <- seagrass[1:80,]

#### 1.2   Rename variables ####
site <- seagrass$Site
dista <- seagrass$Distance

leaves <- seagrass$Leaves # number of leaves (shoot-1)
area <- seagrass$Area # leaf area (cm2 shoot-1)
biomass <- seagrass$Biomass # leaf biomass (g shoot-1)
SLA <- area/biomass # specific leaf area (cm2 g-1)
shoots <- seagrass$Shoots # number of shoots (m-2)

esite <- epiphytes$Site
edist <- epiphytes$Distance
ebiomass <- epiphytes$Biomass # epiphyte biomass (shoot-1)

#### 1.3   Test for assumptions ####
#### 1.3.1 Normality ####
par(mfrow = c(2,1), mar = c(2,2,2,1)) # modify plot environment

hist(leaves)
boxplot(leaves, horizontal = T) # almost normal

hist(area)
boxplot(area, horizontal = T) # almost normal

hist(biomass)
boxplot(biomass, horizontal = T) # slightly right-skewed

hist(SLA)
boxplot(SLA, horizontal = T) # two outliers
SLA <- SLA[-c(9, 22)] # remove outliers
hist(SLA)
boxplot(SLA, horizontal = T) # normal

hist(shoots)
boxplot(shoots, horizontal = T) # normal

hist(ebiomass)
boxplot(ebiomass, horizontal = T) # right-skewed

#### 1.3.2 Homogeneity ####
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1) # return to default

plot(leaves ~ dista, pch = site) # heterogenous
plot(area ~ dista, pch = site) # homogenous
plot(biomass ~ dista, pch = site) # heterogenous

ssite <- site[-c(9, 22)]
sdist <- dista[-c(9, 22)]
plot(SLA ~ sdist, pch = site) # heterogenous
plot(shoots ~ dista, pch = site) # heterogenous

plot(ebiomass ~ edist, pch = esite) # heterogenous


#### 2.    Data anaylsis ####
#### 2.1   Number of leaves ####
#### 2.1.1 Build model ####
m1 <- lm(leaves ~ site * dista)

#### 2.1.2 Determine fixed components ####
m2 <- update(m1, .~. - site : dista) # remove interaction term
anova(m1, m2) # m2 fits better but retain interaction for hypothesis test

#### 2.1.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, pch = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # quite normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m1 is chosen as the optimal model

#### 2.1.4 Interpret model ####
require(car)
# type III sums of squares for overall p-values
Anova(m1, type = 3)
# Response: leaves
#             Sum Sq Df  F value    Pr(>F)    
# (Intercept) 405.71  1 297.2431 < 2.2e-16 ***
# site          9.98  1   7.3085  0.008464 ** 
# dista         3.48  1   2.5494  0.114490    
# site:dista    4.08  1   2.9929  0.087692 .  
# Residuals   103.73 76 

site <- factor(site, levels = c("Zghira", "Kbira"))
m1 <- lm(leaves ~ site * dista)
Anova(m1, type = 3)
# Response: leaves
#              Sum Sq Df F value    Pr(>F)    
# (Intercept)  59.402  1 43.5211 4.992e-09 ***
# site          9.975  1  7.3085  0.008464 ** 
# dista         1.537  1  1.1264  0.291910    
# site:dista    4.085  1  2.9929  0.087692 .  
# Residuals   103.733 76      

site <- factor(site, levels = c("Kbira", "Zghira"))
m1 <- lm(leaves ~ site * dista)

# type II sums of squares for overall p-values
Anova(m2, type = 2)
# Response: leaves
#            Sum Sq Df F value   Pr(>F)   
# site       16.076  1 11.4807 0.001111 **
# dista       0.932  1  0.6656 0.417090   
# Residuals 107.818 77 

#### 2.2   Leaf biomass ####
#### 2.2.1 Build model ####
m3 <- lm(biomass ~ site * dista)

#### 2.2.2 Determine fixed components ####
m4 <- update(m3, .~. - site : dista) # remove interaction term
anova(m3, m4) # m3 fits better

#### 2.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, pch = site)

par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m3) ~ dista)
boxplot(resid(m3) ~ site)
# quite homogenous but residual heterogeneity could
# be modelled by site (greater residual spread at Zghira)

# test assumption of normality
hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normal except for one outlier

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# try improving the model with glm

#### 2.2.4 Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(biomass, "gamma") 
norm <- fitdist(biomass, "norm")
par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
          legendtext = c("Gamma", "Normal"), 
          fitlty = 1)
cdfcomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, gamma fits better 

#### 2.2.5 Build model ####
m5 <- glm(biomass ~ site * dista, family = Gamma(link = "log"))

#### 2.2.6 Determine fixed components ####
m6 <- update(m5, .~. - site : dista) # remove interaction term
anova(m5, m6, test = "Chisq") # m5 fits better

#### 2.2.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m5, pch = site)
# more homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m5))
qqnorm(resid(m5))
qqline(resid(m5)) # more normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m5 is chosen as the optimal model

#### 2.2.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m5, type = 3)
# Response: biomass
#            LR Chisq Df Pr(>Chisq)  
# site         0.0220  1    0.88221  
# dista        1.9099  1    0.16698  
# site:dista   4.4388  1    0.03513 *

coef(m5)
# y = exp(0.003900894*x - 0.796639077)

site <- factor(site, levels = c("Zghira", "Kbira"))
m5 <- glm(biomass ~ site * dista, family = Gamma(link = "log"))
Anova(m5, type = 3)
# Response: biomass
#            LR Chisq Df Pr(>Chisq)   
# site         0.0220  1   0.882212   
# dista       10.5714  1   0.001149 **
# site:dista   4.4388  1   0.035130 * 

coef(m5)
# y = exp(0.01555936*x - 0.84194074)

site <- factor(site, levels = c("Kbira", "Zghira"))
m5 <- glm(biomass ~ site * dista, family = Gamma(link = "log"))

#### 2.3   Leaf area ####
#### 2.3.1 Build model ####
m7 <- lm(area ~ site * dista)

#### 2.3.2 Determine fixed components ####
m8 <- update(m7, .~. - site : dista) # remove interaction term
anova(m7, m8) # m7 fits better

#### 2.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m7, col = site)

par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m7) ~ dista)
boxplot(resid(m7) ~ site)
# quite homogenous

# test assumption of normality
hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7)) # normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m7 is chosen as the optimal model

#### 2.3.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m7, type = 3)
# Response: area
#             Sum Sq Df F value    Pr(>F)    
# (Intercept)  82155  1 28.4294 9.649e-07 ***
# site          4978  1  1.7228  0.193285    
# dista          342  1  0.1184  0.731769    
# site:dista   24824  1  8.5901  0.004462 ** 
# Residuals   219624 76 

coef(m7)
# y = 0.1262619*x + 103.6155388

site <- factor(site, levels = c("Zghira", "Kbira"))
m7 <- lm(area ~ site * dista)
Anova(m7, type = 3)
# Response: area
#             Sum Sq Df F value    Pr(>F)    
# (Intercept)   6446  1  2.2306 0.1394414    
# site          4978  1  1.7228 0.1932855    
# dista        37862  1 13.1018 0.0005293 ***
# site:dista   24824  1  8.5901 0.0044621 ** 
# Residuals   219624 76

coef(m7)
# y = 2.223535*x + 51.592480

# intersection of the two lines
eq <- rbind(c(51.5925, 2.2235), c(103.6155, 0.1263))
c(-solve(cbind(eq[,2],-1)) %*% eq[,1])
# x = 24.80593, y = 106.74849

site <- factor(site, levels = c("Kbira", "Zghira"))
m7 <- lm(area ~ site * dista)

#### 2.4   Number of shoots ####
#### 2.4.1 Build model ####
m9 <- lm(shoots ~ site * dista)

#### 2.4.2 Determine fixed components ####
m10 <- update(m9, .~. - site : dista) # remove interaction term
anova(m9, m10) # m10 fits better but retain interaction for hypothesis test

#### 2.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m9, col = site[c(1:3, 11:13, 21:23, 31:33, 41:43, 51:53, 61:63, 71:73)])

par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m9) ~ dista[c(1:3, 11:13, 21:23, 31:33, 41:43, 51:53, 61:63, 71:73)])
boxplot(resid(m9) ~ site[c(1:3, 11:13, 21:23, 31:33, 41:43, 51:53, 61:63, 71:73)])
# somewhat heterogenous but residual spread 
# does not vary predictably

# test assumption of normality
hist(resid(m9))
qqnorm(resid(m9))
qqline(resid(m9)) # normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# based on good normality, this model is chosen as optimal

#### 2.4.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m9, type = 3)
# Response: shoots
#             Sum Sq Df F value    Pr(>F)    
# (Intercept)  36387  1 31.2383 1.803e-05 ***
# site            51  1  0.0436    0.8367    
# dista         2535  1  2.1764    0.1557    
# site:dista    1615  1  1.3863    0.2529    
# Residuals    23296 20  

site <- factor(site, levels = c("Zghira", "Kbira"))
m9 <- lm(shoots ~ site * dista)
Anova(m9, type = 3)
# Response: shoots
#              Sum Sq Df F value   Pr(>F)   
# (Intercept)  9826.8  1  8.4364 0.008765 **
# site           50.8  1  0.0436 0.836661   
# dista         279.9  1  0.2403 0.629353   
# site:dista   1614.8  1  1.3863 0.252855   
# Residuals   23296.4 20    

site <- factor(site, levels = c("Kbira", "Zghira"))
m9 <- lm(shoots ~ site * dista)

# type II sums of squares for overall p-values
Anova(m10, type = 2)
# Response: shoots
#            Sum Sq Df F value    Pr(>F)    
# site      21661.5  1 18.2605 0.0003381 ***
# dista      1200.2  1  1.0117 0.3259385    
# Residuals 24911.2 21

#### 2.5   Epiphyte biomass ####
#### 2.5.1 Build model ####
m11 <- lm(ebiomass ~ esite * edist)

#### 2.5.2 Determine fixed components ####
m12 <- update(m11, .~. - esite : edist) # remove interaction term
anova(m11, m12) # m11 fits better

#### 2.5.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m11, pch = esite)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m11) ~ esite)
plot(resid(m11) ~ edist)
# residual variance varies predictably with distance and site

# test assumption of normality
hist(resid(m11))
qqnorm(resid(m11))
qqline(resid(m11)) # non-normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving

#### 2.5.4 Fit gamma distribution ####
gamma <- fitdist(ebiomass, "gamma") 
norm <- fitdist(ebiomass, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits much better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits much better

#### 2.5.5 Build gamma model ####
m13 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))

#### 2.5.6 Determine fixed components ####
m14 <- update(m13, .~. - esite : edist) # remove interaction term
anova(m13, m14, test = "Chisq") # m13 fits better

#### 2.5.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m13, pch = esite)
# homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m13))
qqnorm(resid(m13))
qqline(resid(m13)) # residuals are normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m13 is chosen as the optimal model

#### 2.5.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m13, type = 3)
# Response: ebiomass
#             LR Chisq Df Pr(>Chisq)  
# esite         0.0298  1    0.86291  
# edist         3.4591  1    0.06290 .
# esite:edist   4.0845  1    0.04328 *

coef(m13)
# y = exp(0.01028275*x - 3.04316807)

esite <- factor(esite, levels = c("Zghira", "Kbira"))
m13 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))
Anova(m13, type = 3)
# Response: ebiomass
#             LR Chisq Df Pr(>Chisq)    
# esite         0.0298  1  0.8629057    
# edist        11.4220  1  0.0007258 ***
# esite:edist   4.0845  1  0.0432775 * 

coef(m13)
# y = exp(0.03229668*x - 3.14693442)

esite <- factor(esite, levels = c("Kbira", "Zghira"))
m13 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))


#### 2.6   Specific leaf area ####
#### 2.6.1 Build model ####
m15 <- lm(SLA ~ ssite * sdist)

#### 2.6.2 Determine fixed components ####
m16 <- update(m15, .~. - ssite : sdist) # remove interaction term
anova(m15, m16) # m16 fits better but retain interaction for hypothesis test

#### 2.6.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m15, pch = ssite)

par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m15) ~ sdist)
boxplot(resid(m15) ~ ssite)
# somewhat heterogenous but residual spread 
# does not vary predictably

# test assumption of normality
hist(resid(m15))
qqnorm(resid(m15))
qqline(resid(m15)) # normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# this model is chosen as optimal

#### 2.6.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m15, type = 3)
# Response: SLA
#             Sum Sq Df   F value    Pr(>F)    
# (Intercept) 374787  1 1104.6753 < 2.2e-16 ***
# ssite         2648  1    7.8061  0.006626 ** 
# sdist         8028  1   23.6615 6.297e-06 ***
# ssite:sdist    542  1    1.5980  0.210153    
# Residuals    25106 74

coef(m15)
# y = -0.6298916*x + 228.0910761

ssite <- factor(ssite, levels = c("Zghira", "Kbira"))
m15 <- lm(SLA ~ ssite * sdist)
Anova(m15, type = 3)
# Response: SLA
#             Sum Sq Df  F value    Pr(>F)    
# (Intercept)  87298  1 257.3085 < 2.2e-16 ***
# ssite         2648  1   7.8061  0.006626 ** 
# sdist          772  1   2.2752  0.135713    
# ssite:sdist    542  1   1.5980  0.210153    
# Residuals    25106 74 

coef(m15)
# y = -0.3174917*x + 189.8647657

ssite <- factor(ssite, levels = c("Kbira", "Zghira"))
m15 <- lm(SLA ~ ssite * sdist)

# type II sums of squares for overall p-values
Anova(m16, type = 2)
# Response: SLA
#            Sum Sq Df F value    Pr(>F)    
# ssite      8956.7  1  26.191 2.310e-06 ***
# sdist      8257.5  1  24.146 5.106e-06 ***
# Residuals 25648.4 75   

#### 3.    Data visualisation ####
#### 3.1   Calculate descriptive statistics ####
require(psych)
leaves.stat <- describeBy(leaves, list(site, dista), mat = T, digits = 10)
leaves.stat$group2 <- as.numeric(leaves.stat$group2)
l.stat <- describeBy(leaves, site, mat = T, digits = 10)

biomass.stat <- describeBy(biomass, list(site, dista), mat = T, digits = 10)
biomass.stat$group2 <- as.numeric(biomass.stat$group2)
b.stat <- describeBy(biomass, site, mat = T, digits = 10)

area.stat <- describeBy(area, list(site, dista), mat = T, digits = 10)
area.stat$group2 <- as.numeric(area.stat$group2)
a.stat <- describeBy(area, site, mat = T, digits = 10)

shoots.stat <- describeBy(shoots, list(site, dista), mat = T, digits = 10)
shoots.stat$group2 <- as.numeric(shoots.stat$group2)
s.stat <- describeBy(shoots, site, mat = T, digits = 10)

ebiomass.stat <- describeBy(ebiomass, list(esite, edist), mat = T, digits = 10)
ebiomass.stat$group2 <- as.numeric(ebiomass.stat$group2)
e.stat <- describeBy(ebiomass, esite, mat = T, digits = 10)

SLA.stat <- describeBy(SLA, list(ssite, sdist), mat = T, digits = 10)
SLA.stat$group2 <- as.numeric(SLA.stat$group2)
sla.stat <- describeBy(SLA, ssite, mat = T, digits = 10)

#### 3.2   Calculate model predictions ####
# re-run models
site <- factor(site, levels = c("Zghira", "Kbira"))
esite <- factor(esite, levels = c("Zghira", "Kbira"))
ssite <- factor(ssite, levels = c("Zghira", "Kbira"))

m1 <- lm(leaves ~ site * dista)
m5 <- glm(biomass ~ site * dista, family = Gamma(link = "log"))
m7 <- lm(area ~ site * dista)
m9 <- lm(shoots ~ site * dista)
m13 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))
m15 <- lm(SLA ~ ssite * sdist)

# generate new data frames for model predictions to work on
new <- data.frame(site = c(rep("Kbira", 496), rep("Zghira", 323)),
                  dista = c(seq(21.22, 70.78, by = 0.1),
                            seq(36.35, 68.63, by = 0.1)))
new$site <- as.factor(new$site)
new$site <- factor(new$site, levels = c("Zghira", "Kbira"))

enew <- data.frame(esite = c(rep("Kbira", 496), rep("Zghira", 323)),
                   edist = c(seq(21.22, 70.78, by = 0.1),
                             seq(36.35, 68.63, by = 0.1)))
enew$esite <- as.factor(enew$esite)
enew$esite <- factor(enew$esite, levels = c("Zghira", "Kbira"))

slanew <- data.frame(ssite = c(rep("Kbira", 496), rep("Zghira", 323)),
                     sdist = c(seq(21.22, 70.78, by = 0.1),
                               seq(36.35, 68.63, by = 0.1)))
slanew$ssite <- as.factor(slanew$ssite)
slanew$ssite <- factor(slanew$ssite, levels = c("Zghira", "Kbira"))

# calculate predicted values and intervals
pred.l <- data.frame(predict(m1, interval = "confidence", newdata = new)) 
new$l.fit <- pred.l$fit # add line fit to dataset
new$l.hi <- pred.l$upr # upper confidence interval boundary
new$l.lo <- pred.l$lwr # lower confidence interval boundary

inv <- family(m5)$linkinv # calculate inverse of link function
pred.b <- data.frame(predict(m5, type = "link", newdata = new,
                             se.fit = TRUE))

new$b.fit <- inv(pred.b$fit) # add line fit to dataset
new$b.hi <- inv(pred.b$fit + pred.b$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$b.lo <- inv(pred.b$fit - pred.b$se.fit * qnorm(0.975)) # lower confidence interval boundary

pred.a <- data.frame(predict(m7, interval = "confidence", newdata = new)) 
new$a.fit <- pred.a$fit # add line fit to dataset
new$a.hi <- pred.a$upr # upper confidence interval boundary
new$a.lo <- pred.a$lwr # lower confidence interval boundary

pred.s <- data.frame(predict(m9, interval = "confidence", newdata = new))
new$s.fit <- pred.s$fit # add line fit to dataset
new$s.hi <- pred.s$upr # upper confidence interval boundary
new$s.lo <- pred.s$lwr # lower confidence interval boundary

inv <- family(m13)$linkinv # calculate inverse of link function
pred.e <- data.frame(predict(m13, type = "link", newdata = enew,
                             se.fit = TRUE)) 

enew$fit <- inv(pred.e$fit) # add line fit to dataset
enew$hi <- inv(pred.e$fit + pred.e$se.fit * qnorm(0.975)) # upper confidence interval boundary
enew$lo <- inv(pred.e$fit - pred.e$se.fit * qnorm(0.975)) # lower confidence interval boundary

pred.sla <- data.frame(predict(m15, interval = "confidence", newdata = slanew))
slanew$sla.fit <- pred.sla$fit # add line fit to dataset
slanew$sla.hi <- pred.sla$upr # upper confidence interval boundary
slanew$sla.lo <- pred.sla$lwr # lower confidence interval boundary

#### 3.3   Calculate intersections of exponential curves ####
bdf <- data.frame(pred = predict(m5, type = "response", 
                                 newdata = data.frame(site = c(rep("Kbira", 3001), rep("Zghira", 3001)),
                                                      dista = c(seq(0, 30, by = 0.01),
                                                                seq(0, 30, by = 0.01)))),
                 site = c(rep("Kbira", 3001), rep("Zghira", 3001)),
                 dista = c(seq(0, 30, by = 0.01),
                           seq(0, 30, by = 0.01)))

bdf$diff <- rep(bdf$pred[3002:6002] - bdf$pred[1:3001], 2)
bdf$dista[which.min(abs(bdf$diff - 0))] # x = 3.89

edf <- data.frame(pred = predict(m13, type = "response", 
                                 newdata = data.frame(esite = c(rep("Kbira", 3001), rep("Zghira", 3001)),
                                                      edist = c(seq(0, 30, by = 0.01),
                                                                seq(0, 30, by = 0.01)))),
                  esite = c(rep("Kbira", 3001), rep("Zghira", 3001)),
                  edist = c(seq(0, 30, by = 0.01),
                            seq(0, 30, by = 0.01)))

edf$diff <- edf$pred[3002:6002] - edf$pred[1:3001] 
edf$edist[which.min(abs(edf$diff - 0))] # x = 4.71


#### 3.4   Define theme ####
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

#### 3.5   Plot ####
lp <- ggplot() +
        geom_line(data = new, aes(dista, l.fit, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = l.lo, ymax = l.hi, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Leaves, colour = Site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = leaves.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                                fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = l.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Number of leaves (shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(2, 12), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = c(0.35, 0.95)) +
        mytheme
lp # dimensions: 4.5 x 3.5

bp <- ggplot() +
        geom_line(data = new, aes(dista, b.fit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = b.lo, ymax = b.hi, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Biomass, colour = Site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = biomass.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                                 fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = b.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Leaf biomass (g shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 3), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 3, by = 0.5)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
bp # dimensions: 4.5 x 3.5

ap <- ggplot() +
        geom_line(data = new, aes(dista, a.fit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = a.lo, ymax = a.hi, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Area, colour = Site),
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = area.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                                 fill = group1), size = 0.5, shape = 21,
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = a.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Leaf area (cm"^2*" shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 350), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 350, by = 50)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
ap # dimensions: 4.5 x 3.5

sp <- ggplot() +
        geom_line(data = new, aes(dista, s.fit, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = s.lo, ymax = s.hi, fill = site),
                    alpha = 0.5) +
        geom_point(data = seagrass, aes(Distance, Shoots, colour = Site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = shoots.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                                fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = s.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Number of shoots (m"^-2*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 250), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
sp # dimensions: 4.5 x 3.5

ep <- ggplot() +
        geom_line(data = enew, aes(edist, fit, colour = esite, lty = esite), size = 0.5) +
        geom_ribbon(data = enew, aes(edist, ymin = lo, ymax = hi, fill = esite), 
                    alpha = 0.5) +
        geom_point(data = epiphytes, aes(Distance, Biomass, colour = Site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = ebiomass.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se, 
                                                  fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = e.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Epiphyte biomass (g shoot"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 1.2), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 1.2, by = 0.2)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
ep # dimensions: 4.5 x 3.5

slap <- ggplot() +
        geom_line(data = slanew, aes(sdist, sla.fit, colour = ssite, lty = ssite), size = 0.5) +
        geom_ribbon(data = slanew, aes(sdist, ymin = sla.lo, ymax = sla.hi, fill = ssite),
                    alpha = 0.5) +
        geom_point(aes(sdist, SLA, colour = ssite),
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = SLA.stat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                             fill = group1), size = 0.5, shape = 21,
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = sla.stat, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = F) +
        ylab(expression("Specific leaf area (cm"^2*" g"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(120, 260), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(120, 260, by = 20)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
slap # dimensions: 4.5 x 3.5

#### 3.6   Combine plots ####
require(cowplot)
final.plot <- plot_grid(lp + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        ap + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        bp + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        slap, ep, sp, rel_heights = c(0.928, 1), nrow = 2, align = "v",
                        labels = "auto", label_size = 15, label_fontfamily = "Helvetica Neue")
final.plot # dimensions: 9 x 10.5 in

#### 4.    Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:psych)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")


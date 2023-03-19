##########################################################################################
## Project: Impacts of thermal effluent on Posidonia oceanica and associated macrofauna ##
## Script purpose: Analysis of physico-chemical data                                    ##
## Author: Luka Seamus Wright                                                           ##
##########################################################################################

#### 1.    Data exploration ####
#### 1.1   Load data ####
physchem <- read.csv("~/PATH/physchem.csv")

#### 1.2   Test for relationships ####
require(ggplot2)
ggplot(physchem, aes(Distance, Temperature, colour = Site)) +
  geom_point() + facet_grid(cols = vars(Month)) 
# trend = temperature is higher at Zghira and declines with distance,
# while it does not change with distance at Kbira

ggplot(physchem, aes(Distance, Temperature, colour = Site)) +
  geom_smooth()
# trend still clear across months but large confidence intervals


ggplot(physchem, aes(Distance, pH, colour = Site)) +
  geom_point() + facet_grid(cols = vars(Month))
# trend = pH is generally lower at Zghira
# this becomes clearer when boxplotted
ggplot(physchem, aes(Site, pH)) +
  geom_boxplot()

ggplot(physchem, aes(Distance, Salinity, colour = Site)) +
  geom_point() + facet_grid(cols = vars(Month)) # no trend
# remove outliers
ggplot(physchem[-c(882, 1000),], aes(Distance, Salinity, colour = Site)) +
  geom_point() + facet_grid(cols = vars(Month)) # still no trend

# calculate site, month and overall means
require(psych)
with(physchem[-c(882, 1000),], describeBy(Salinity, list(Site, Month), mat = T))
with(physchem[-c(882, 1000),], describeBy(Salinity, Site), mat = T)

# conclusion: there are relationships for the variables 
# temperature and pH but these should be analysed separately
# for each month because of large intra-annual differences

#### 1.3   Split dataframe by month ####
July <- physchem[1:500,]
September <- physchem[501:1000,]
February <- physchem[1001:1500,]
June <- physchem[1501:2010,] 
# in June variables were measured in the meadow rather than at the sea surface

#### 1.4   Rename variables ####
Jsite <- factor(July$Site)
Jdist <- July$Distance
Jtemp <- July$Temperature

Ssite <- factor(September$Site)
Sdist <- September$Distance
Stemp <- September$Temperature
SpH <- September$pH

Fsite <- factor(February$Site)
Fdist <- February$Distance
Ftemp <- February$Temperature
FpH <- February$pH

Msite <- factor(June$Site) # M for meadow
Mdist <- June$Distance
Mtemp <- June$Temperature
MpH <- June$pH

#### 2.    Data analysis ####
#### 2.1   July Temperature ####
#### 2.1.1 Build models ####
m1 <- lm(Jtemp ~ Jsite * Jdist)

#### 2.1.2 Determine fixed components ####
m2 <- update(m1, .~. - Jsite : Jdist) # remove interaction term
m3 <- update(m1, .~. - Jdist - Jsite : Jdist) # remove interaction term and second variable
anova(m1, m2, m3) # interaction term is not significant
# but retain for hypothesis test

#### 2.1.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, pch = Jsite)
# not very homogenous

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m1) ~ Jsite) # more residual variance in Zghira
plot(resid(m1) ~ Jdist) # more resudual variance for lower distances
# -> residual variance varies with distance and site, so can potentially 
# be modelled with nlme package

# test assumption of normality
hist(resid(m1)) # residuals seem quite normally distributed
qqnorm(resid(m1))
qqline(resid(m1)) # but deviate from normality at the upper distribution edge

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# perhaps normality could be improved by fitting a gamma distribution

#### 2.1.4 Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(Jtemp, "gamma") 
norm <- fitdist(Jtemp, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# the gamma distribution only marginally fits better

#### 2.1.5 Model heterogeneity ####
require(nlme)
m4 <- gls(Jtemp ~ Jsite * Jdist, method = "REML")
m5 <- gls(Jtemp ~ Jsite * Jdist, weights = varIdent(form = ~ 1 | Jsite), 
           method = "REML")
m6 <- gls(Jtemp ~ Jsite * Jdist, weights = varExp(form = ~ Jdist | Jsite), 
           method = "REML")
m7 <- gls(Jtemp ~ Jsite * Jdist, 
           weights = varComb(varIdent(form = ~ 1 | Jsite),
                             varExp(form = ~ Jdist | Jsite)), method = "REML")
anova(m4, m5, m6, m7) # m7 has the lowers AIC and BIC scores
# -> proceed with m7

#### 2.1.6 Determine fixed components ####
m7 <- gls(Jtemp ~ Jsite * Jdist, 
           weights = varComb(varIdent(form = ~ 1 | Jsite),
                             varExp(form = ~ Jdist | Jsite)), method = "ML")

m8 <- update(m7, .~. - Jsite : Jdist) # remove interaction term
m9 <- update(m7, .~. - Jdist - Jsite : Jdist) # remove interaction term and second variable
anova(m8, m7, m9) # interaction is not significant
# but retain for hypothesis test

m7 <- gls(Jtemp ~ Jsite * Jdist, 
          weights = varComb(varIdent(form = ~ 1 | Jsite),
                            varExp(form = ~ Jdist | Jsite)), method = "REML")

#### 2.1.7 Test model fit ####
plot(m7, pch = Jsite)
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7)) # normality is the same as m2

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m8 is chosen as the optimal model

#### 2.1.8 Interpret model ####
# type II sums of squares for overall p-values
require(car)
Anova(m7, type = 3)
# Response: Jtemp
#             Df      Chisq Pr(>Chisq)    
# (Intercept)  1 3.5373e+05  < 2.2e-16 ***
# Jsite        1 4.3228e+01  4.872e-11 ***
# Jdist        1 1.8446e+02  < 2.2e-16 ***
# Jsite:Jdist  1 2.3865e+00     0.1224 

# pairwise contrasts not required


#### 2.2   September Temperature ####
#### 2.2.1 Build model ####
m10 <- lm(Stemp ~ Ssite * Sdist)

#### 2.2.2 Determine fixed components ####
m11 <- update(m10, .~. - Ssite : Sdist) # remove interaction term
m12 <- update(m10, .~. - Sdist - Ssite : Sdist) # remove interaction term and second variable
anova(m10, m11, m12) # m10 is still the best model

#### 2.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m10, pch = Ssite)
# clearly heterogenous

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m10) ~ Ssite) # more residual variance in Zghira
plot(resid(m10) ~ Sdist) # more resudual variance for lower distances
# -> residual variance varies with distance and site, so can potentially 
# be modelled with nlme package

# test assumption of normality
hist(resid(m10)) # residuals seem quite normally distributed
qqnorm(resid(m10))
qqline(resid(m10)) # but deviate from normality at the upper and lower distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# perhaps normality could be improved by fitting a gamma distribution

#### 2.2.4 Fit gamma distribution ####
gamma <- fitdist(Stemp, "gamma") 
norm <- fitdist(Stemp, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# the gamma distribution only marginally fits better

#### 2.2.5 Model heterogeneity ####
m13 <- gls(Stemp ~ Ssite * Sdist, method = "REML")
m14 <- gls(Stemp ~ Ssite * Sdist, weights = varIdent(form = ~ 1 | Ssite), 
           method = "REML")
m15 <- gls(Stemp ~ Ssite * Sdist, weights = varExp(form = ~ Sdist | Ssite), 
           method = "REML")
m16 <- gls(Stemp ~ Ssite * Sdist, 
           weights = varComb(varIdent(form = ~ 1 | Ssite),
                             varExp(form = ~ Sdist | Ssite)), method = "REML")
anova(m13, m14, m15, m16) # m16 has the lowers AIC and BIC scores
# -> proceed with m16

#### 2.2.6 Determine fixed components ####
m16 <- gls(Stemp ~ Ssite * Sdist, 
           weights = varComb(varIdent(form = ~ 1 | Ssite),
                             varExp(form = ~ Sdist | Ssite)), method = "ML")

m17 <- update(m16, .~. - Ssite : Sdist) # remove interaction term
m18 <- update(m16, .~. - Sdist - Ssite : Sdist) # remove interaction term and second variable
anova(m17, m16, m18) # m16 is the best model
# -> proceed with m16
m16 <- gls(Stemp ~ Ssite * Sdist, 
           weights = varComb(varIdent(form = ~ 1 | Ssite),
                             varExp(form = ~ Sdist | Ssite)), method = "REML")

#### 2.2.7 Test model fit ####
plot(m16, pch = Ssite)
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m16))
qqnorm(resid(m16))
qqline(resid(m16)) # normality is the same as m10

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m16 is chosen as the optimal model

#### 2.2.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m16, type = 3)
# Response: Stemp
#             Df      Chisq Pr(>Chisq)    
# (Intercept)  1 1.4186e+06  < 2.2e-16 ***
# Ssite        1 1.5338e+02  < 2.2e-16 ***
# Sdist        1 1.3518e+00      0.245    
# Ssite:Sdist  1 1.5931e+01  6.571e-05 *** 

# pairwise contrasts
summary(m16)
# slope for Kbira: p = 0.26
# Kbira vs. Zghira intercept: p < 0.001 ***
# Kbira vs. Zghira slope: p < 0.001 ***

Ssite <- factor(Ssite, levels = c("Zghira", "Kbira"))
m16 <- gls(Stemp ~ Ssite * Sdist, 
           weights = varComb(varIdent(form = ~ 1 | Ssite),
                             varExp(form = ~ Sdist | Ssite)), method = "REML")
summary(m16)
# slope for Zghira: p < 0.001 ***

Ssite <- factor(Ssite, levels = c("Kbira", "Zghira"))


#### 2.3   February Temperature ####
#### 2.3.1 Build model ####
m19 <- lm(Ftemp ~ Fsite * Fdist)

#### 2.3.2 Determine fixed components ####
m20 <- update(m19, .~. - Fsite : Fdist) # remove interaction term
m21 <- update(m19, .~. - Fdist - Fsite : Fdist) # remove interaction term and second variable
anova(m19, m20, m21) # m19 is still the best model

#### 2.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m19, pch = Fsite)
# heterogenous

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m19) ~ Fsite) # more residual variance in Zghira
plot(resid(m19) ~ Fdist) # resudual variance somewhat lower for larger distances
# -> residual variance varies with distance and site, so can potentially 
# be modelled with nlme package

# test assumption of normality
hist(resid(m19)) # residuals seem quite normally distributed
qqnorm(resid(m19))
qqline(resid(m19)) # but deviate strongly from normality at the upper and lower distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# perhaps normality could be improved by fitting a gamma distribution

#### 2.3.4 Fit gamma distribution ####
gamma <- fitdist(Ftemp, "gamma") 
norm <- fitdist(Ftemp, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# the gamma distribution fits better,
# but not enough to make a marked difference in normality

#### 2.3.5 Model heterogeneity ####
m22 <- gls(Ftemp ~ Fsite * Fdist, method = "REML")
m23 <- gls(Ftemp ~ Fsite * Fdist, weights = varIdent(form = ~ 1 | Fsite), 
           method = "REML")
m24 <- gls(Ftemp ~ Fsite * Fdist, weights = varExp(form = ~ Fdist | Fsite), 
           method = "REML")
m25 <- gls(Ftemp ~ Fsite * Fdist, 
           weights = varComb(varIdent(form = ~ 1 | Fsite),
                             varExp(form = ~ Fdist | Fsite)), method = "REML")
anova(m22, m23, m24, m25) # m25 has the lowers AIC and BIC scores
# -> proceed with m25

#### 2.3.6 Determine fixed components ####
m25 <- gls(Ftemp ~ Fsite * Fdist, 
           weights = varComb(varIdent(form = ~ 1 | Fsite),
                             varExp(form = ~ Fdist | Fsite)), method = "ML")

m26 <- update(m25, .~. - Fsite : Fdist) # remove interaction term
m27 <- update(m25, .~. - Fdist - Fsite : Fdist) # remove interaction term and second variable
anova(m25, m26, m27) # m25 is the best model
# -> proceed with m25
m25 <- gls(Ftemp ~ Fsite * Fdist, 
           weights = varComb(varIdent(form = ~ 1 | Fsite),
                             varExp(form = ~ Fdist | Fsite)), method = "REML")

#### 2.3.7 Test model fit ####
plot(m25, pch = Fsite)
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m25))
qqnorm(resid(m25))
qqline(resid(m25)) # normality is the same as m19

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# this model is chosen as optimal because it
# improves homogeneity substantially

#### 2.3.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m25, type = 3)
# Response: Ftemp
#             Df      Chisq Pr(>Chisq)    
# (Intercept)  1 3.9238e+06  < 2.2e-16 ***
# Fsite        1 6.4750e+02  < 2.2e-16 ***
# Fdist        1 2.8253e+01  1.065e-07 ***
# Fsite:Fdist  1 1.3241e+01  0.0002739 ***

# pairwise contrasts
summary(m25)
# slope for Kbira: p < 0.001 ***
# Kbira vs. Zghira intercept: p < 0.001 ***
# Kbira vs. Zghira slope: p < 0.001 ***

Fsite <- factor(Fsite, levels = c("Zghira", "Kbira"))
m25 <- gls(Ftemp ~ Fsite * Fdist, 
           weights = varComb(varIdent(form = ~ 1 | Fsite),
                             varExp(form = ~ Fdist | Fsite)), method = "REML")
summary(m25)
# slope for Zghira: p = 0.003 **

Fsite <- factor(Fsite, levels = c("Kbira", "Zghira"))



#### 2.4   September pH ####
#### 2.4.1 Build model ####
m28 <- lm(SpH ~ Ssite * Sdist)

#### 2.4.2 Determine fixed components ####
m29 <- update(m28, .~. - Ssite : Sdist) # remove interaction term
anova(m28, m29) # interaction term is not significant
# but retain for hypothesis test

#### 2.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))

plot(m28, col = Ssite)
# fairly homogenous

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m28) ~ Ssite) # variance is slightly larger at Kbira
plot(resid(m28) ~ Sdist) # and increases somewhat with distance
# however, overall homogeneity is given

# test assumption of normality
hist(resid(m28)) # residuals are not normally distributed
qqnorm(resid(m28))
qqline(resid(m28)) # and deviate from normality at the upper and lower distribution edges
# however, this cannot be fixed by modelling with a different distribution

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m28 is chosen as the optimal model

#### 2.4.4 Interpret model ####
# type II sums of squares for overall p-values
Anova(m28, type = 3)
# Response: SpH
#             Sum Sq  Df    F value    Pr(>F)    
# (Intercept) 4166.3   1 34533.2483 < 2.2e-16 ***
# Ssite          2.5   1    20.9211 6.053e-06 ***
# Sdist          0.0   1     0.0572    0.8111    
# Ssite:Sdist    0.1   1     0.9707    0.3250    
# Residuals     59.8 496  

Ssite <- factor(Ssite, levels = c("Zghira", "Kbira"))
m28 <- lm(SpH ~ Ssite * Sdist)
Anova(m28, type = 3)
# Response: SpH
#             Sum Sq  Df    F value    Pr(>F)    
# (Intercept) 3881.3   1 32170.9704 < 2.2e-16 ***
# Ssite          2.5   1    20.9211 6.053e-06 ***
# Sdist          0.3   1     2.6650    0.1032    
# Ssite:Sdist    0.1   1     0.9707    0.3250    
# Residuals     59.8 496   

Ssite <- factor(Ssite, levels = c("Kbira", "Zghira"))
m28 <- lm(SpH ~ Ssite * Sdist)

#### 2.5   February pH ####
#### 2.5.1 Build model ####
m30 <- lm(FpH ~ Fsite * Fdist)

#### 2.5.2 Determine fixed components ####
m31 <- update(m30, .~. - Fsite : Fdist) # remove interaction term
anova(m30, m31) # interaction is not significant
# but retain for hypothesis test

#### 2.5.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m30, col = Fsite)
# heterogenous

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m30) ~ Fsite) # more residual variance in Kbira
plot(resid(m30) ~ Fdist) # residual variance does not change with distance
# -> residual variance varies between sites, so can potentially 
# be modelled with nlme package

# test assumption of normality
hist(resid(m30)) # residuals seem balanced
qqnorm(resid(m30))
qqline(resid(m30)) # but deviate strongly from normality at the upper and lower distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# perhaps normality could be improved by fitting a gamma distribution

#### 2.5.4 Fit gamma distribution ####
gamma <- fitdist(FpH, "gamma") 
norm <- fitdist(FpH, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# the gamma distribution fits better,
# but not enough to make a marked difference in normality

#### 2.5.5 Model heterogeneity ####
m32 <- gls(FpH ~ Fsite * Fdist, weights = varIdent(form = ~ 1 | Fsite), 
           method = "ML")

#### 2.3.6 Determine fixed components ####
m33 <- update(m32, .~. - Fsite : Fdist) # remove interaction term
anova(m32, m33) # interaction term is not significant
# but retain for hypothesis test
m32 <- gls(FpH ~ Fsite * Fdist, weights = varIdent(form = ~ 1 | Fsite), 
           method = "REML")

#### 2.5.7 Test model fit ####
plot(m32, col = Fsite)
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m32))
qqnorm(resid(m32))
qqline(resid(m32)) # normality is the same as m30

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# this model is chosen as optimal because it
# improves homogeneity substantially

#### 2.5.8 Interpret model ####
# type II sums of squares for overall p-values
Anova(m32, type = 3)
# Response: FpH
#             Df      Chisq Pr(>Chisq)    
# (Intercept)  1 67280.9947     <2e-16 ***
# Fsite        1   854.8485     <2e-16 ***
# Fdist        1     0.2108     0.6461    
# Fsite:Fdist  1     0.9137     0.3391 

Fsite <- factor(Fsite, levels = c("Zghira", "Kbira"))
m32 <- gls(FpH ~ Fsite * Fdist, weights = varIdent(form = ~ 1 | Fsite), 
           method = "REML")
Anova(m32, type = 3)
# Response: FpH
#             Df      Chisq Pr(>Chisq)    
# (Intercept)  1 1.6807e+07     <2e-16 ***
# Fsite        1 8.5485e+02     <2e-16 ***
# Fdist        1 7.8794e+01     <2e-16 ***
# Fsite:Fdist  1 9.1370e-01     0.3391 



#### 2.6   June Temperature ####
#### 2.6.1 Build model ####
m34 <- lm(Mtemp ~ Msite * Mdist)

#### 2.6.2 Determine fixed components ####
m35 <- update(m34, .~. - Msite : Mdist) # remove interaction term
m36 <- update(m34, .~. - Mdist - Msite : Mdist) # remove interaction term and second variable
anova(m34, m35, m36) # m34 is the best model

#### 2.6.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m34, pch = Msite)
# nonlinear trend, try modelling with gamma
# or nonlinear exponential decay

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m34) ~ Msite) # more residual variance in Zghira
plot(resid(m34) ~ Mdist) # residual variance somewhat lower for larger distances
# -> residual variance varies with distance and site, so can potentially 
# be modelled with nlme package

# test assumption of normality
hist(resid(m34)) # residuals are right-skewed
qqnorm(resid(m34))
qqline(resid(m34))

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# perhaps normality and homogeneity could be improved by fitting a gamma distribution

#### 2.6.4 Fit gamma distribution ####
gamma <- fitdist(Mtemp, "gamma") 
norm <- fitdist(Mtemp, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# the gamma distribution fits better

#### 2.6.5 Determine fixed components ####
m37 <- glm(Mtemp ~ Msite * Mdist, 
           family = Gamma(link = "log"))

m38 <- update(m37, .~. - Msite : Mdist) # remove interaction term
m39 <- update(m37, .~. - Mdist - Msite : Mdist) # remove interaction term and second variable
anova(m37, m38, m39, test = "Chisq") # m37 is the best model

#### 2.6.6 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))

plot(m37, pch = Msite)
# homogeneity is somewhat improved but gamma isn't enough to describe 
# the curvature of the data

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m37))
qqnorm(resid(m37))
qqline(resid(m37)) # normality is the same as m34

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

# fit nonlinear exponential decay model

#### 2.6.7 Self-starting function ####
m40 <- nls(Temperature ~ SSasymp(Distance, A, I, logk), data = June)
summary(m40)
# starting parameters to be used are 
# A = 23.51206, C = 1.62355 (I-A), k = -0.0458656 (exp(logk))

#### 2.6.8 Fit model ####
June$Site <- factor(June$Site)
m41 <- nls(Temperature ~ C[Site] * exp(k[Site] * Distance) + A[Site],
           start = list(C = c(1.62355, 1.62355),
                        k = c(-0.0458656, -0.0458656),
                        A = c(23.51206, 23.51206)),
           data = June)

# convert to gnls for confidence interval purposes (note none of the paramaters changed)
m41 <- gnls(Temperature ~ C * exp(k * Distance) + A,
            start = list(C = c(1.62355, 1.62355),
                         k = c(-0.0458656, -0.0458656),
                         A = c(23.51206, 23.51206)),
            params = list(C ~ Site, k ~ Site, A ~ Site),
            data = June) 

#### 2.6.9 Test model fit ####
plot(resid(m41) ~ Mdist) # residual variance slightly decreases with distance
abline(0,0)
# fairly homogenous

m42 <- gnls(Temperature ~ C * exp(k * Distance) + A,
            start = list(C = c(1.62355, 1.62355),
                         k = c(-0.0458656, -0.0458656),
                         A = c(23.51206, 23.51206)),
            params = list(C ~ Site, k ~ Site, A ~ Site),
            weights = varExp(),
            data = June) 

plot(resid(m42, type = "normalized") ~ Mdist) # weighting overcorrects homogeneity
abline(0,0)
# stick with m41

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m41))
qqnorm(resid(m41))
qqline(resid(m41)) # normality is ok 
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 2.6.12 Interpret model ####
summary(m41)
#                   Value  Std.Error   t-value p-value
# C.(Intercept)  0.427599 0.17651141   2.42250  0.0158
# C.SiteZghira   2.401851 0.25991361   9.24096  0.0000
# k.(Intercept) -0.034447 0.04013022  -0.85838  0.3911
# k.SiteZghira  -0.013269 0.04079640  -0.32524  0.7451
# A.(Intercept) 23.168530 0.13944463 166.14860  0.0000
# A.SiteZghira   0.679681 0.16838848   4.03639  0.0001

June$Site <- factor(June$Site, levels = c("Zghira", "Kbira"))
m41 <- gnls(Temperature ~ C * exp(k * Distance) + A,
            start = list(C = c(1.62355, 1.62355),
                         k = c(-0.0458656, -0.0458656),
                         A = c(23.51206, 23.51206)),
            params = list(C ~ Site, k ~ Site, A ~ Site),
            data = June) 
summary(m41)
#                   Value  Std.Error   t-value p-value
# C.(Intercept)  2.829449 0.19078471  14.83059  0.0000
# C.SiteKbira   -2.401851 0.25991362  -9.24096  0.0000
# k.(Intercept) -0.047715 0.00734248  -6.49853  0.0000
# k.SiteKbira    0.013268 0.04079643   0.32524  0.7451
# A.(Intercept) 23.848211 0.09439215 252.65036  0.0000
# A.SiteKbira   -0.679681 0.16838813  -4.03639  0.0001

June$Site <- factor(June$Site, levels = c("Kbira", "Zghira"))
m41 <- gnls(Temperature ~ C * exp(k * Distance) + A,
            start = list(C = c(1.62355, 1.62355),
                         k = c(-0.0458656, -0.0458656),
                         A = c(23.51206, 23.51206)),
            params = list(C ~ Site, k ~ Site, A ~ Site),
            data = June) 

#### 2.7   June pH ####
#### 2.7.1 Build model ####
m43 <- lm(MpH ~ Msite * Mdist)

#### 2.7.2 Determine fixed components ####
m44 <- update(m43, .~. - Msite : Mdist) # remove interaction term
anova(m43, m44) # m43 is the better model

#### 2.7.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))

plot(m43, col = Msite)
# fairly homogenous but some noticeable nonlinearity

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m43) ~ Msite)
plot(resid(m43) ~ Mdist) # nonlinearity
# overall homogenous

# test assumption of normality
hist(resid(m43)) # residuals are normally distributed
qqnorm(resid(m43))
qqline(resid(m43))

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m43 is good but try polynomial regression for nonlinearity

#### 2.7.4 Build polynomial model ####
m45 <- lm(MpH ~ poly(Mdist, 3, raw = T) * Msite)
  
#### 2.7.5 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))

plot(m45, col = Msite)
# homogeneity is improved

# plot of residuals against site and distance variables
par(mfrow = c(1,2), mar = c(2,2,2,1))

boxplot(resid(m45) ~ Msite)
plot(resid(m45) ~ Mdist)
# overall homogenous

# test assumption of normality
hist(resid(m45)) # residuals are fairly normally distributed
qqnorm(resid(m45))
qqline(resid(m45))

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m45 is chosen as the optimal model


#### 2.7.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m45, type = 3)
# Response: MpH
#                                Sum Sq  Df    F value    Pr(>F)    
# (Intercept)                   1214.16   1 1.8169e+06 < 2.2e-16 ***
# poly(Mdist, 3, raw = T)          0.09   3 4.4969e+01 < 2.2e-16 ***
# Msite                            0.03   1 5.1532e+01 2.551e-12 ***
# poly(Mdist, 3, raw = T):Msite    0.07   3 3.2952e+01 < 2.2e-16 ***
# Residuals                        0.34 502    

Msite <- factor(Msite, levels = c("Zghira", "Kbira"))
m45 <- lm(MpH ~ poly(Mdist, 3, raw = T) * Msite)
Anova(m45, type = 3)
# Response: MpH
#                                Sum Sq  Df    F value    Pr(>F)    
# (Intercept)                   1195.94   1 1.7896e+06 < 2.2e-16 ***
# poly(Mdist, 3, raw = T)          0.02   3 1.1002e+01 5.236e-07 ***
# Msite                            0.03   1 5.1532e+01 2.551e-12 ***
# poly(Mdist, 3, raw = T):Msite    0.07   3 3.2952e+01 < 2.2e-16 ***
# Residuals                        0.34 502 

Msite <- factor(Msite, levels = c("Kbira", "Zghira"))
m45 <- lm(MpH ~ poly(Mdist, 3, raw = T) * Msite)


#### 3.    Data visualisation ####
#### 3.1   Calculate descriptive statistics ####
Stemp.stat <- describeBy(Stemp, list(Ssite, Sdist), mat = T, digits = 10)
Ftemp.stat <- describeBy(Ftemp, list(Fsite, Fdist), mat = T, digits = 10)
Mtemp.stat <- describeBy(Mtemp, list(Msite, Mdist), mat = T, digits = 10)

SpH.stat <- describeBy(SpH, list(Ssite, Sdist), mat = T, digits = 10)
FpH.stat <- describeBy(FpH, list(Fsite, Fdist), mat = T, digits = 10)
MpH.stat <- describeBy(MpH, list(Msite, Mdist), mat = T, digits = 10)

physchem.stat <- rbind(Stemp.stat, Mtemp.stat, Ftemp.stat, SpH.stat, MpH.stat, FpH.stat)
physchem.stat$group2 <- as.numeric(physchem.stat$group2)
physchem.stat$Month <- c(rep("September",100), rep("June",102), rep("February",100),
                         rep("September",100), rep("June",102), rep("February",100))

#### 3.2   Calculate model predictions ####
# re-run all optimal models
m16 <- gls(Stemp ~ Ssite * Sdist, 
           weights = varComb(varIdent(form = ~ 1 | Ssite),
                             varExp(form = ~ Sdist | Ssite)), method = "REML")

m25 <- gls(Ftemp ~ Fsite * Fdist, 
           weights = varComb(varIdent(form = ~ 1 | Fsite),
                             varExp(form = ~ Fdist | Fsite)), method = "REML")

m28 <- lm(SpH ~ Ssite * Sdist)

m32 <- gls(FpH ~ Fsite * Fdist, weights = varIdent(form = ~ 1 | Fsite), 
           method = "REML")

m41 <- gnls(Temperature ~ C * exp(k * Distance) + A,
            start = list(C = c(1.62355, 1.62355),
                         k = c(-0.0458656, -0.0458656),
                         A = c(23.51206, 23.51206)),
            params = list(C ~ Site, k ~ Site, A ~ Site),
            data = June)

m45 <- lm(MpH ~ poly(Mdist, 3, raw = T) * Msite)

# calculate predicted values and intervals
September$fit.temp <- predict(m16) # predicted values from m16
modmat <-  model.matrix(formula(m16)[-2]) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m16) %*% t(modmat))
September$lower.temp <- with(September, fit.temp - qnorm(0.975)*sqrt(int))
September$upper.temp <- with(September, fit.temp + qnorm(0.975)*sqrt(int))

February$fit.temp <- predict(m25) # predicted values from m25
modmat <-  model.matrix(formula(m25)[-2]) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m25) %*% t(modmat))
February$lower.temp <- with(February, fit.temp - qnorm(0.975)*sqrt(int))
February$upper.temp <- with(February, fit.temp + qnorm(0.975)*sqrt(int))

pred <- data.frame(predict(m28, interval = "confidence")) # predicted values and interval from m28
September$fit.pH <- pred$fit
September$lower.pH <- pred$lwr
September$upper.pH <- pred$upr

February$fit.pH <- predict(m32) # predicted values from m33
modmat <-  model.matrix(formula(m32)[-2]) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m32) %*% t(modmat))
February$lower.pH <- with(February, fit.pH - qnorm(0.975)*sqrt(int))
February$upper.pH <- with(February, fit.pH + qnorm(0.975)*sqrt(int))


# new <- data.frame(Distance = rep(seq(0, 100, by = 0.5), 2),
#                   Site = c(rep("Kbira", 201), rep("Zghira", 201))) # generated explanatory variable
# new$Site <- factor(new$Site) 
June$fit.temp <- predict(m41, newdata = June) # predicted values from m41

bootfun <- function(newdata) {
  start <- coef(m41)
  boot <- June[sample(nrow(June), size = nrow(June), replace = TRUE),]
  bootfit <- try(update(m41,
                        start = start,
                        data = boot),
                 silent = TRUE)
  if (inherits(bootfit, "try-error")) return(rep(NA, nrow(newdata)))
  predict(bootfit, newdata)
}

bmat <- replicate(1000, bootfun(June))
June$lower.temp <- apply(bmat, 1, quantile, 0.025, na.rm = TRUE)
June$upper.temp <- apply(bmat, 1, quantile, 0.975, na.rm = TRUE) # bootstrap 95% confidence interval

pred <- data.frame(predict(m45, interval = "confidence")) # predicted values and interval from m28
June$fit.pH <- pred$fit
June$lower.pH <- pred$lwr
June$upper.pH <- pred$upr

physchem <- rbind(September, June, February)
rownames(physchem) <- NULL
physchem$sig.temp <- c(rep("NS", 250), rep("S", 505), rep("NS", 255), rep("S", 500))
physchem$sig.pH <- c(rep("NS", 500), rep("S", 510), rep("NS", 250), rep("S", 250))

#### 3.3   Customise theme ####
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .4, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 3.4   Plot ####
physchem$Site <- factor(physchem$Site, levels = c("Zghira", "Kbira"))
physchem$Month <- factor(physchem$Month)
levels(physchem$Month) <- c("February (sea surface)","June (seagrass meadow)","September (sea surface)")
names(physchem.stat)[names(physchem.stat) == "group1"] <- "Site"
physchem.stat$Site <- factor(physchem.stat$Site, levels = c("Zghira", "Kbira"))
physchem.stat$Month <- factor(physchem.stat$Month)
levels(physchem.stat$Month) <- c("February (sea surface)","June (seagrass meadow)","September (sea surface)")

tp <- ggplot(data = physchem) +
        # geom_vline(xintercept = c(21.22, 70.78, 36.35, 68.63),
        #            colour = rep(c(rep("#6ea4be",2), rep("#f5a54a",2)),3)) +
        geom_pointrange(data = physchem.stat, aes(group2, mean, colour = Site,
                                                  ymin = mean-se, ymax = mean+se),
                        alpha = 0.5, size = 0.3, shape = 16) +
        geom_line(data = physchem, aes(Distance, fit.temp, colour = Site,
                      lty = sig.temp), size = 0.5) +
        geom_ribbon(data = physchem, aes(Distance, ymin = lower.temp, ymax = upper.temp,
                        fill = Site), alpha = 0.45) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = "none") +
        ylab("Temperature (°C)") +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(14.727, 29.273)) +
        scale_y_continuous(breaks = seq(14, 30, by = 4)) +
        scale_x_continuous(expand = c(0,0)) +
        facet_grid(~Month) +
        theme(legend.position = c(.1, .92),
              strip.background = element_rect(colour = "black", size = 1,
                                              fill = NA),
              strip.text = element_text(size = 15, hjust = 0),
              panel.spacing = unit(1, "cm")) +
        mytheme

tp # dimensions: 4.5 x 8 in

pHp <- ggplot(data = physchem) +
        # geom_vline(xintercept = c(21.22, 70.78, 36.35, 68.63),
        #            colour = rep(c(rep("#6ea4be",2), rep("#f5a54a",2)),3)) +
        geom_pointrange(data = physchem.stat, aes(group2, mean, colour = Site,
                                            ymin = mean-se, ymax = mean+se),
                        alpha = 0.5, size = 0.3, shape = 16) +
        geom_line(data = physchem, aes(Distance, fit.pH, colour = Site, lty = sig.pH),
                  size = 0.5) +
        geom_ribbon(data = physchem, aes(Distance, ymin = lower.pH, ymax = upper.pH,
                        fill = Site), alpha = 0.45) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = "none") +
        ylab("pH") +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(7.2545, 8.3455)) +
        scale_y_continuous(breaks = seq(7.2, 8.4, by = 0.4)) +
        scale_x_continuous(expand = c(0,0)) +
        facet_grid(~Month) +
        theme(legend.position = "none",
              strip.background = element_rect(colour = "black", size = 1,
                                              fill = NA),
              strip.text = element_text(size = 15, hjust = 0),
              panel.spacing = unit(1, "cm")) +
        mytheme

pHp # dimensions: 4.5 x 6 in


#### 3.5   Combined plot ####
require(cowplot)
tpH <- plot_grid(tp + theme(axis.text.x = element_blank(),
                            axis.title.x = element_blank()), 
                 pHp + theme(strip.text = element_blank(),
                             strip.background = element_blank()), 
                 rel_heights = c(1, 1.02), ncol = 1, 
                 labels = "auto", align = "v",
                 label_size = 15, label_fontfamily = "Helvetica Neue")
tpH # dimensions: 9 x 10.5 in

#### 4.    Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:nlme)
detach(package:psych)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")

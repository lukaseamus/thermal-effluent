##########################################################################################
## Project: Impacts of thermal effluent on Posidonia oceanica and associated macrofauna ##
## Script purpose: Analysis of seagrass habitat data                                    ##
## Author: Luka Seamus Wright                                                           ##
##########################################################################################

#### 1.    Data exploration ####
#### 1.1   Load data ####
seagrass <- read.csv("~/PATH/seagrass.habitat.csv")
epiphytes <- seagrass[81:160,]
seagrass <- seagrass[1:80,]

#### 1.2   Rename variables ####
site <- seagrass$Site
dista <- seagrass$Distance

leaves <- seagrass$Leaves # number of leaves (m-2)
area <- seagrass$Area # leaf area (m2 m-2)
biomass <- seagrass$Biomass # leaf biomass (g m-2)

esite <- epiphytes$Site
edist <- epiphytes$Distance
ebiomass <- epiphytes$Biomass # epiphyte biomass (m-2)

#### 1.3   Test for assumptions ####
#### 1.3.1 Normality ####
par(mfrow = c(2,1), mar = c(2,2,2,1)) # modify plot environment

hist(leaves)
boxplot(leaves, horizontal = T) # almost normal

hist(area)
boxplot(area, horizontal = T) # almost normal

hist(biomass)
boxplot(biomass, horizontal = T) # normal except for outlier

hist(ebiomass)
boxplot(ebiomass, horizontal = T) # right-skewed

#### 1.3.2 Homogeneity ####
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1) # return to default

plot(leaves ~ dista, pch = site) # heterogenous
plot(area ~ dista, pch = site) # homogenous
plot(biomass ~ dista, pch = site) # homogenous

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
qqline(resid(m1)) # normality could be better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 2.1.4 Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(leaves, "gamma") 
norm <- fitdist(leaves, "norm")

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

#### 2.1.5 Build model ####
m3 <- glm(leaves ~ site * dista, family = Gamma(link = "log"))

#### 2.1.6 Determine fixed components ####
m4 <- update(m3, .~. - site : dista) # remove interaction term
anova(m3, m4, test = "Chisq") # m4 fits better but retain interaction for hypothesis test

#### 2.1.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, pch = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normality is perfect

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model

#### 2.1.8 Interpret model ####
require(car)
# type III sums of squares for overall p-values
Anova(m3, type = 3)
# Response: leaves
#            LR Chisq Df Pr(>Chisq)  
# site         6.4400  1    0.01116 *
# dista        1.5417  1    0.21436  
# site:dista   1.1748  1    0.27843

site <- factor(site, levels = c("Zghira", "Kbira"))
m3 <- glm(leaves ~ site * dista, family = Gamma(link = "log"))
Anova(m3, type = 3)
# Response: leaves
#            LR Chisq Df Pr(>Chisq)  
# site         6.4400  1    0.01116 *
# dista        0.2776  1    0.59830  
# site:dista   1.1748  1    0.27843 

site <- factor(site, levels = c("Kbira", "Zghira"))
m3 <- glm(leaves ~ site * dista, family = Gamma(link = "log"))

# type II sums of squares for overall p-values
Anova(m4, type = 2)
# Response: leaves
#       LR Chisq Df Pr(>Chisq)    
# site   140.245  1     <2e-16 ***
# dista    0.637  1     0.4247 

#### 2.2   Leaf biomass ####
#### 2.2.1 Build model ####
m5 <- lm(biomass ~ site * dista)

#### 2.2.2 Determine fixed components ####
m6 <- update(m5, .~. - site : dista) # remove interaction term
anova(m5, m6) # m6 fits better but retain interaction for hypothesis test

#### 2.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m5, col = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m5))
qqnorm(resid(m5))
qqline(resid(m5)) # normal except for one outlier

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m5 is chosen as the optimal model

#### 2.2.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m5, type = 3)
# Response: biomass
#             Sum Sq Df F value    Pr(>F)    
# (Intercept)  21723  1 15.7677 0.0001613 ***
# site           413  1  0.2998 0.5856203    
# dista         9787  1  7.1041 0.0093895 ** 
# site:dista    1066  1  0.7738 0.3818065    
# Residuals   104706 76 

site <- factor(site, levels = c("Zghira", "Kbira"))
m5 <- lm(biomass ~ site * dista)
Anova(m5, type = 3)
# Response: biomass
#             Sum Sq Df F value Pr(>F)  
# (Intercept)   3552  1  2.5780 0.1125
# site           413  1  0.2998 0.5856  
# dista         9436  1  6.8491 0.0107 *
# site:dista    1066  1  0.7738 0.3818  
# Residuals   104706 76

site <- factor(site, levels = c("Kbira", "Zghira"))
m5 <- lm(biomass ~ site * dista)

# type II sums of squares for overall p-values
Anova(m6, type = 2)
# Response: biomass
#           Sum Sq Df F value    Pr(>F)    
# site        1214  1  0.8841 0.3500341    
# dista      18157  1 13.2182 0.0004986 ***
# Residuals 105772 77 

#### 2.3   Leaf area ####
#### 2.3.1 Build model ####
m7 <- lm(area ~ site * dista)

#### 2.3.2 Determine fixed components ####
m8 <- update(m7, .~. - site : dista) # remove interaction term
anova(m7, m8) # m8 fits better but retain interaction for hypothesis test

#### 2.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m7, col = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7)) # normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m7 is chosen as the optimal model

#### 2.3.4 Interpret model ####
# type III sums of squares for overall p-values
Anova(m7, type = 3)
# Response: area
#              Sum Sq Df F value    Pr(>F)    
# (Intercept) 12.3526  1 31.2520 3.394e-07 ***
# site         0.3161  1  0.7997   0.37401    
# dista        1.7899  1  4.5283   0.03658 *  
# site:dista   0.1743  1  0.4411   0.50860    
# Residuals   30.0397 76 

site <- factor(site, levels = c("Zghira", "Kbira"))
m7 <- lm(area ~ site * dista)
Anova(m7, type = 3)
# Response: area
#              Sum Sq Df F value  Pr(>F)  
# (Intercept)  1.7745  1  4.4895 0.03737 *
# site         0.3161  1  0.7997 0.37401  
# dista        1.6530  1  4.1819 0.04432 *
# site:dista   0.1743  1  0.4411 0.50860  
# Residuals   30.0397 76 

site <- factor(site, levels = c("Kbira", "Zghira"))
m7 <- lm(area ~ site * dista)

# type II sums of squares for overall p-values
Anova(m8, type = 2)
# Response: area
#            Sum Sq Df F value   Pr(>F)   
# site       0.2866  1  0.7304 0.395389   
# dista      3.2685  1  8.3296 0.005058 **
# Residuals 30.2140 77  

#### 2.4   Epiphyte biomass ####
#### 2.4.1 Build model ####
m9 <- lm(ebiomass ~ esite * edist)

#### 2.4.2 Determine fixed components ####
m10 <- update(m9, .~. - esite : edist) # remove interaction term
anova(m10, m9) # m9 fits better

#### 2.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m9, pch = esite)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m9) ~ esite)
plot(resid(m9) ~ edist)
# residual variance varies predictably with distance and site

# test assumption of normality
hist(resid(m9))
qqnorm(resid(m9))
qqline(resid(m9)) # non-normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving

#### 2.4.4 Fit gamma distribution ####
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

#### 2.4.5 Build gamma model ####
m11 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))

#### 2.4.6 Determine fixed components ####
m12 <- update(m11, .~. - esite : edist) # remove interaction term
anova(m11, m12, test = "Chisq") # m12 fits better but retain interaction for hypothesis test

#### 2.4.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m11, pch = esite)
# homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m11))
qqnorm(resid(m11))
qqline(resid(m11)) # residuals are normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m11 is chosen as the optimal model

#### 2.4.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m11, type = 3)
# Response: ebiomass
#             LR Chisq Df Pr(>Chisq)   
# esite         0.0576  1   0.810367   
# edist         7.5009  1   0.006167 **
# esite:edist   1.6972  1   0.192661

esite <- factor(esite, levels = c("Zghira", "Kbira"))
m11 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))
Anova(m11, type = 3)
# Response: ebiomass
#             LR Chisq Df Pr(>Chisq)   
# esite         0.0576  1   0.810367   
# edist         9.5394  1   0.002011 **
# esite:edist   1.6972  1   0.192661

esite <- factor(esite, levels = c("Kbira", "Zghira"))
m11 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))

# type II sums of squares for overall p-values
Anova(m12, type = 2)
# Response: ebiomass
#       LR Chisq Df Pr(>Chisq)    
# esite   9.6972  1  0.0018455 ** 
# edist  14.5806  1  0.0001343 ***

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

ebiomass.stat <- describeBy(ebiomass, list(esite, edist), mat = T, digits = 10)
ebiomass.stat$group2 <- as.numeric(ebiomass.stat$group2)
e.stat <- describeBy(ebiomass, esite, mat = T, digits = 10)

#### 3.2   Calculate model predictions ####
# re-run models
site <- factor(site, levels = c("Zghira", "Kbira"))
esite <- factor(esite, levels = c("Zghira", "Kbira"))

m3 <- glm(leaves ~ site * dista, family = Gamma(link = "log"))
m5 <- lm(biomass ~ site * dista)
m7 <- lm(area ~ site * dista)
m11 <- glm(ebiomass ~ esite * edist, family = Gamma(link = "log"))

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

# calculate predicted values and intervals
inv <- family(m3)$linkinv # calculate inverse of link function
pred.l <- data.frame(predict(m3, type = "link", newdata = new,
                             se.fit = TRUE))
new$l.fit <- inv(pred.l$fit) # add line fit to dataset
new$l.hi <- inv(pred.l$fit + pred.l$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$l.lo <- inv(pred.l$fit - pred.l$se.fit * qnorm(0.975)) # lower confidence interval boundary

pred.b <- data.frame(predict(m5, interval = "confidence", newdata = new)) 
new$b.fit <- pred.b$fit # add line fit to dataset
new$b.hi <- pred.b$upr # upper confidence interval boundary
new$b.lo <- pred.b$lwr # lower confidence interval boundary

pred.a <- data.frame(predict(m7, interval = "confidence", newdata = new)) 
new$a.fit <- pred.a$fit # add line fit to dataset
new$a.hi <- pred.a$upr # upper confidence interval boundary
new$a.lo <- pred.a$lwr # lower confidence interval boundary

inv <- family(m11)$linkinv # calculate inverse of link function
pred.e <- data.frame(predict(m11, type = "link", newdata = enew,
                             se.fit = TRUE)) 
enew$fit <- inv(pred.e$fit) # add line fit to dataset
enew$hi <- inv(pred.e$fit + pred.e$se.fit * qnorm(0.975)) # upper confidence interval boundary
enew$lo <- inv(pred.e$fit - pred.e$se.fit * qnorm(0.975)) # lower confidence interval boundary

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
        ylab(expression("Number of leaves (m"^-2*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(expand = c(0,0), breaks = seq(200, 1800, by = 400)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        coord_cartesian(ylim = c(200, 1800), xlim = c(10, 80)) +
        theme(legend.position = c(0.375, 0.95)) +
        mytheme
lp # dimensions: 4.5 x 3.5

bp <- ggplot() +
        geom_line(data = new, aes(dista, b.fit, colour = site), size = 0.5) +
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
        ylab(expression("Leaf biomass (g m"^-2*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        coord_cartesian(ylim = c(0, 250), xlim = c(10, 80), clip = "off") +
        theme(legend.position = "none") +
        mytheme
bp # dimensions: 4.5 x 3.5

ap <- ggplot() +
        geom_line(data = new, aes(dista, a.fit, colour = site), size = 0.5) +
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
        ylab(expression("Leaf area (m"^2*" m"^-2*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        coord_cartesian(ylim = c(0, 4), xlim = c(10, 80)) +
        theme(legend.position = "none") +
        mytheme
ap # dimensions: 4.5 x 3.5

ep <- ggplot() +
        geom_line(data = enew, aes(edist, fit, colour = esite), size = 0.5) +
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
        ylab(expression("Epiphyte biomass (g m"^-2*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0, 120), xlim = c(10, 80)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 120, by = 30)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
ep # dimensions: 4.5 x 3.5

#### 3.5   Combine plots ####
require(cowplot)
final.plot <- plot_grid(lp + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        ap + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        bp, ep, rel_heights = c(0.928, 1), align = "v",
                        nrow = 2, labels = "auto", label_size = 15, label_fontfamily = "Helvetica Neue")
final.plot # dimensions: 9 x 7 in

#### 4.    Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:psych)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")

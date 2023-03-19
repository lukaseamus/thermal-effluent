##########################################################################################
## Project: Impacts of thermal effluent on Posidonia oceanica and associated macrofauna ##
## Script purpose: Analysis of macrofauna data                                          ##
## Author: Luka Seamus Wright                                                           ##
##########################################################################################

#### 1.    Load data ####
fauna <- read.csv("~/PATH/fauna.csv")

#### 2.    Data wrangling ####
spec <- data.frame(fauna[,4:70], row.names = 1) # species abundance data
env <- fauna[,1:3] # environmental data

#### 3.    Multivariate data analysis ####
require(vegan)
#### 3.1   Calculate Bray-Curtis dissimilarity ####
specdist <- vegdist(spec, method = "bray")

#### 3.2   Test for homogeneity of dispersion ####
Sbetad <- with(env, betadisper(specdist, Site)) # dispersion by site
permutest(Sbetad, permutations = 9999) # homogenous

Pbetad <- with(env, betadisper(specdist, Plot)) # dispersion by plot
permutest(Pbetad, permutations = 9999) # homogenous

#### 3.3   PERMANOVA ####
adonis(spec ~ Site, env, permutations = 9999) # model incorporating site
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Site       1    1.1526 1.15258  8.8616 0.38762  1e-04 ***
# Residuals 14    1.8209 0.13007         0.61238           
# Total     15    2.9735                 1.00000   

adonis(spec[1:8,] ~ Plot, env[1:8,], permutations = 9999) # model incorporating Kbira plots  
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Plot       3   0.48658 0.16219  1.2768 0.48918 0.2587
# Residuals  4   0.50811 0.12703         0.51082       
# Total      7   0.99469                 1.00000  

adonis(spec[9:16,] ~ Plot, env[9:16,], permutations = 9999) # model incorporating Zghira plots  
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Plot       3   0.48869 0.162895  1.9304 0.59147 0.0175 *
# Residuals  4   0.33754 0.084384         0.40853         
# Total      7   0.82622                  1.00000 


#### 3.4   SIMPER ####
# SIMPER is only interesting where PERMANOVA detected significant dissimilarity
simp <- with(env, simper(spec, Site, permutations = 9999)) # SIMPER for site
summary(simp, ordered = T)
# Contrast: Kbira_Zghira 
#                                    ava    avb   cumsum p
# Hippolyte.inermis                  9.750 25.125 0.1860 **
# Jujubinus.exasperatus              1.875 14.250 0.3329 ***
# Rissoa.variabilis                  1.000  8.625 0.4357 ***
# Mysida.sp.A                        0.750  8.250 0.5173 ***
# Hyale.sp                           5.125  0.500 0.5707
# Alvania.discors                    1.250  5.125 0.6198 *
# Leptochelia.sp                     0.125  5.000 0.6621 **
# Rissoa.oriscalpium                 0.750  2.750 0.6895 *
# Amphitoe.ramondi                   2.250  2.000 0.7132
# Leucothoe.spinicapra               0.000  2.250 0.7359 ***
# Asterina.gibbosa                   0.125  1.750 0.7541 *
# Rissoa.sp..A                       0.000  1.375 0.7712 ***
# Gibbula.umbilicaris                1.500  0.625 0.7875
# Syllidae.sp.B                      1.250  0.125 0.8022
# Alvania.lineata                    0.750  1.250 0.8149
# Gnathia.sp                         0.250  0.875 0.8260
# Tricolia.pullus                    0.125  0.875 0.8364 *
# Dexamine                           0.625  0.625 0.8458
# Nereis.rava                        0.250  0.625 0.8551
# Tricolia.speciosa                  0.125  0.750 0.8637
# Rissoa.sp..B..juvenile.oriscalpum. 0.000  0.625 0.8718 ***
# Calliostoma.conulus                0.125  0.750 0.8798
# Photis.sp.                         0.375  0.625 0.8876
# Columbella.rustica                 0.250  0.500 0.8943
# Opeatogenys.gracilis               0.250  0.500 0.9008
# Jujubinus.striatus                 0.000  0.500 0.9062 ***
# Aoridae.sp                         0.375  0.250 0.9112
# Cheirocratus.sundevalli            0.000  0.375 0.9160 ***
# Ampeliska                          0.375  0.000 0.9205
# Nematonereis.unicornis             0.250  0.250 0.9249
# Amphilocus                         0.250  0.125 0.9291
# Apherusa.bispinosa                 0.250  0.250 0.9333
# Apseudes.sp                        0.125  0.250 0.9373
# Idotea.sp                          0.125  0.250 0.9411
# Lysianassa                         0.000  0.375 0.9449 ***
# Syllidae.sp.A                      0.125  0.250 0.9486
# Cestopagurus.timidus               0.125  0.250 0.9520
# Bittium.latreillii                 0.000  0.250 0.9553 ***
# Apherusa.sp                        0.250  0.000 0.9581
# Amphipholis.squamata               0.000  0.250 0.9607 ***
# Cirolana                           0.250  0.000 0.9633
# Calcinus.ornatus                   0.000  0.250 0.9657 ***
# Idotea.linearis                    0.125  0.125 0.9679
# Ophellidae.sp                      0.125  0.125 0.9700
# Anapagurus.sp                      0.000  0.250 0.9721 ***
# Isopoda.sp                         0.000  0.125 0.9739 ***
# Palaemon.sp                        0.125  0.000 0.9755
# Iphimedia                          0.000  0.125 0.9772 ***
# Philinopsis.sp.                    0.000  0.125 0.9788 ***
# Nemertea.sp                        0.000  0.125 0.9803 ***
# Nudibranchia.sp.                   0.000  0.125 0.9818 ***
# Alvania.mamillata                  0.125  0.000 0.9831
# Muricopsis.cristata                0.125  0.000 0.9845
# Cumacea.sp                         0.000  0.125 0.9859 ***
# Syllidae.sp.C                      0.000  0.125 0.9872 ***
# Mitrolumna.olivoidea               0.000  0.125 0.9886 ***
# Polyophthalmus.sp                  0.000  0.125 0.9899 ***
# Rissoa.sp..C                       0.000  0.125 0.9912 ***
# Gibbula.ardens                     0.000  0.125 0.9925 ***
# Phyllodocidae                      0.125  0.000 0.9937
# Gibberula.miliaria                 0.125  0.000 0.9949
# Eurydice                           0.000  0.125 0.9959 ***
# Opistobranchia                     0.000  0.125 0.9969 ***
# Monophorus.sp.                     0.000  0.125 0.9979 ***
# Trivia.pulex                       0.000  0.125 0.9990 ***
# Aspidosiphon.muelleri              0.000  0.125 1.0000 ***

Zsimp <- with(env[9:16,], simper(spec[9:16,], Plot, permutations = 9999)) # SIMPER for Zghira plots
summary(Zsimp, ordered = T)
# because of the abundance of data, only species and between-plot comparisons
# with significant contributions are shown

# A vs B
#                                     average       sd   ratio  ava  avb cumsum      p 
# Mysida.sp.A                        0.061109 0.022534  2.7119  3.5 13.5 0.3108 0.0148 * 
# Alvania.discors                    0.042044 0.022986  1.8291  2.0  8.5 0.6091 0.0331 * 
# Leucothoe.spinicapra               0.026824 0.023266  1.1529  0.0  4.5 0.6664 0.0449 * 
# Gnathia.sp                         0.009512 0.004012  2.3708  1.5  0.0 0.7266 0.0224 * 
# Nematonereis.unicornis             0.006282 0.000520 12.0808  0.0  1.0 0.8339 0.0050 **
# Syllidae.sp.B                      0.003348 0.003869  0.8652  0.0  0.5 0.9053 0.0228 * 
# Nemertea.sp                        0.003348 0.003869  0.8652  0.0  0.5 0.9339 0.0228 * 
# Nudibranchia.sp.                   0.003348 0.003869  0.8652  0.0  0.5 0.9410 0.0228 * 
# Amphilocus                         0.002935 0.003391  0.8654  0.0  0.5 0.9875 0.0251 * 
# Gibbula.ardens                     0.002935 0.003391  0.8654  0.0  0.5 1.0000 0.0251 * 

# A vs C
#                                     average       sd  ratio  ava  avb cumsum      p  
# Rissoa.variabilis                  0.099610 0.053200 1.8724 17.5  2.5 0.4110 0.0049 **
# Rissoa.sp..A                       0.020597 0.012789 1.6104  1.0  4.0 0.6141 0.0438 * 
# Rissoa.oriscalpium                 0.020017 0.006542 3.0599  3.5  0.5 0.6518 0.0431 * 
# Rissoa.sp..B..juvenile.oriscalpum. 0.017940 0.020737 0.8651  0.0  2.5 0.7211 0.0236 * 
# Calliostoma.conulus                0.009635 0.002858 3.3710  0.0  1.5 0.8304 0.0344 *
# Iphimedia                          0.003588 0.004147 0.8651  0.0  0.5 0.9277 0.0236 * 
# Philinopsis.sp.                    0.003588 0.004147 0.8651  0.0  0.5 0.9344 0.0236 * 
# Polyophthalmus.sp                  0.003023 0.003494 0.8654  0.0  0.5 0.9943 0.0226 * 
# Rissoa.sp..C                       0.003023 0.003494 0.8654  0.0  0.5 1.0000 0.0226 * 

# A vs D
#                                     average        sd  ratio  ava  avb cumsum      p 
# Hippolyte.inermis                  0.145254 0.0558733 2.5997 11.0 36.5 0.2307 0.0154 * 
# Leptochelia.sp                     0.099305 0.1147248 0.8656  0.0 20.0 0.3885 0.0225 * 
# Asterina.gibbosa                   0.029272 0.0166324 1.7599  0.0  5.0 0.6851 0.0140 *
# Hyale.sp                           0.012471 0.0144114 0.8653  0.0  2.0 0.7573 0.0246 * 
# Lysianassa                         0.008718 0.0043437 2.0071  0.0  1.5 0.7989 0.0049 **
# Calcinus.ornatus                   0.005600 0.0007526 7.4417  0.0  1.0 0.8691 0.0049 **
# Anapagurus.sp                      0.004965 0.0057362 0.8656  0.0  1.0 0.8946 0.0225 * 
# Cumacea.sp                         0.003118 0.0036029 0.8653  0.0  0.5 0.9194 0.0246 * 
# Syllidae.sp.C                      0.003118 0.0036029 0.8653  0.0  0.5 0.9243 0.0246 * 
# Mitrolumna.olivoidea               0.003118 0.0036029 0.8653  0.0  0.5 0.9293 0.0246 * 
# Idotea.linearis                    0.002483 0.0028681 0.8656  0.0  0.5 0.9763 0.0225 * 
# Ophellidae.sp                      0.002483 0.0028681 0.8656  0.0  0.5 0.9803 0.0225 * 
# Eurydice                           0.002483 0.0028681 0.8656  0.0  0.5 0.9842 0.0225 * 
# Opistobranchia                     0.002483 0.0028681 0.8656  0.0  0.5 0.9882 0.0225 * 
# Monophorus.sp.                     0.002483 0.0028681 0.8656  0.0  0.5 0.9921 0.0225 * 
# Trivia.pulex                       0.002483 0.0028681 0.8656  0.0  0.5 0.9961 0.0225 * 
# Aspidosiphon.muelleri              0.002483 0.0028681 0.8656  0.0  0.5 1.0000 0.0225 * 
  

# B vs C
#                                     average        sd  ratio  ava  avb cumsum      p 
# Rissoa.sp..A                       0.021881 0.0080909 2.7044  0.0  4.0 0.4568 0.0341 * 
# Rissoa.oriscalpium                 0.021541 0.0050980 4.2254  4.5  0.5 0.5149 0.0098 **
# Nematonereis.unicornis             0.005376 0.0005596 9.6068  1.0  0.0 0.8280 0.0404 * 
  
#### 4.   Multivariate data visualisation ####
#### 4.1  Calculate nMDS data ####
nMDS <- metaMDS(spec, distance = "bray", autotransform = FALSE)
nMDS$stress # Stress = 0.06751876

#### 4.2  Extract points and calculate 95% confidence ellipses ####
df <- data.frame(nMDS$points, group = env$Site, plot = env$Plot)
df$group <- as.factor(df$group)
df$group <- factor(df$group, levels = c("Zghira", "Kbira"))
ordiplot(nMDS)
ell <- ordiellipse(nMDS, env$Site, display = "sites", kind = "se",
                   conf = 0.95, label = F)

#### 4.3  Make data compatible with ggplot2 ####
# the following function was sourced from 
# https://github.com/jarioksa/vegan/blob/master/R/veganCovEllipse.R
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  Q <- chol(cov, pivot = TRUE)
  o <- attr(Q, "pivot")
  t(center + scale * t(Circle %*% Q[,o]))
}


df.ell <- data.frame()
for(g in levels(df$group)) {
  df.ell <- rbind(df.ell, 
                 cbind(as.data.frame(with(df[df$group==g,],
                                          veganCovEllipse(ell[[g]]$cov,
                                                          ell[[g]]$center,
                                                          ell[[g]]$scale))),
                       group = g))
}


#### 4.4  Define theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill = NA, size = 1),
                 axis.line = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.margin = margin(t = -7, r = -7, b = -7, l = -7, 
                                        unit = "pt"),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 4.5  Plot ####
specp <- ggplot(data = df, aes(MDS1, MDS2)) +
                annotate("segment", x = 0.05, xend = 0.87, y = -0.7, yend = 0.35, 
                         arrow = arrow(type = "closed", angle = 15), size = 1) +
                geom_point(aes(colour = group, shape = plot, fill = group),
                             size = 3) +
                geom_line(aes(group = interaction(plot, group), colour = group), size = 1) +
                geom_polygon(data = df.ell, aes(NMDS1, NMDS2, colour = group,
                                                fill = group), size = 0.7, alpha = 0.3) +
                scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                                    labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira")) +
                scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                                  labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira")) +
                scale_shape_manual(values = c(22, 23, 24, 25),
                                   labels = c("Plot A", "Plot B", "Plot C", "Plot D"),
                                   guide = guide_legend(ncol = 2)) +
                annotate("text", label = "Stress = 0.07", x = -0.845, y = -0.42,
                         size = 4.2) +
                coord_cartesian(ylim = c(0.55, -0.8)) +
                scale_y_reverse() +
                mytheme +
                theme(legend.position = c(0.22, 0.86))
specp # dimensions: 4.5 x 4.5 in


#### 5.    Univariate data exploration ####
#### 5.1   Calculate diversity indices ####
fauna$A <- rowSums(spec) # Abundance
fauna$S <- specnumber(spec) # Species richness
fauna$D <- diversity(spec, index = "invsimpson", base = exp(1)) # Simpson's reciprocal index (1/D)
fauna$D2 <- diversity(spec, index = "simpson", base = exp(1)) # Simpson's diversity index (1-D)
fauna$H <- diversity(spec, index = "shannon", base = exp(1)) # Shannon-Wiener diversity index
fauna$E <- fauna$D/fauna$S # Simpson's evenness
fauna$J <- fauna$H/log(fauna$S) # Pielou's evenness

#### 5.2   Look for trends ####
with(fauna, plot(A ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(S ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(D ~ Distance, pch = Site))
with(fauna, plot(D2 ~ Distance, pch = Site))
with(fauna, plot(H ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(E ~ Distance, pch = Site))
with(fauna, plot(J ~ Distance, pch = Site)) # potentially interesting

with(fauna, plot(Crustacea ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Crustacea.S ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Mollusca ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Mollusca.S ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Echinodermata ~ Distance, pch = Site)) # this is probably just Asterina gibbosa
with(fauna, plot(Echinodermata.S ~ Distance, pch = Site))
with(fauna, plot(Polychaeta ~ Distance, pch = Site))
with(fauna, plot(Polychaeta.S ~ Distance, pch = Site))
with(fauna, plot(Sipunculida ~ Distance, pch = Site))
with(fauna, plot(Sipunculida.S ~ Distance, pch = Site))
with(fauna, plot(Vertebrata ~ Distance, pch = Site))
with(fauna, plot(Vertebrata.S ~ Distance, pch = Site))

with(fauna, plot(Hippolyte.inermis ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Jujubinus.exasperatus ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Rissoa.variabilis ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Mysida.sp.A ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Leptochelia.sp ~ Distance, pch = Site))
with(fauna, plot(Asterina.gibbosa ~ Distance, pch = Site)) # potentially interesting
with(fauna, plot(Rissoa.sp..A ~ Distance, pch = Site))
with(fauna, plot(Alvania.discors ~ Distance, pch = Site)) # potentially interesting

#### 5.3   Rename interesting variables ####
A <- fauna$A
S <- fauna$S
H <- fauna$H
J <- fauna$J
C <- fauna$Crustacea
CS <- fauna$Crustacea.S
M <- fauna$Mollusca
MS <- fauna$Mollusca.S
Hi <- fauna$Hippolyte.inermis
Je <- fauna$Jujubinus.exasperatus
Rv <- fauna$Rissoa.variabilis
My <- fauna$Mysida.sp.A
Ag <- fauna$Asterina.gibbosa
Ad <- fauna$Alvania.discors

site <- fauna$Site
dista <- fauna$Distance

#### 6.    Univariate index data analysis ####
#### 6.1   Adundance ####
#### 6.1.1 Build model ####
m1 <- lm(A ~ site * dista)

#### 6.1.2 Determine fixed components ####
m2 <- update(m1, .~. - site : dista) # remove interaction term
anova(m1, m2) # models are not different but retain interaction for hypothesis test
# continue with m1

#### 6.1.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, pch = site)
# slightly heterogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 6.1.4 Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(A, "gamma") 
norm <- fitdist(A, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits slightly better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits slightly better

#### 6.1.5 Build gamma model ####
m3 <- glm(A ~ site * dista, family = Gamma(link = "log"))

#### 6.1.6 Determine fixed components ####
m4 <- update(m3, .~. - site : dista) # remove interaction term
anova(m3, m4, test = "Chisq") # m4 fits best but retain interaction term for hypothesis est

#### 6.1.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, pch = site)
# homogeneity seems similar

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normality is slightly better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model

#### 6.1.8 Interpret model ####
# type III sums of squares for overall p-values
require(car)
Anova(m3, type = 3)
# Response: A
#            LR Chisq Df Pr(>Chisq)  
# site         5.1378  1    0.02341 *
# dista        1.0674  1    0.30153  
# site:dista   0.2190  1    0.63978   

site <- factor(site, levels = c("Zghira", "Kbira"))
m3 <- glm(A ~ site * dista, family = Gamma(link = "log"))
Anova(m3, type = 3)
# Response: A
#            LR Chisq Df Pr(>Chisq)  
# site         5.1378  1    0.02341 *
# dista        1.4113  1    0.23484  
# site:dista   0.2190  1    0.63978 

site <- factor(site, levels = c("Kbira", "Zghira"))
m3 <- glm(A ~ site * dista, family = Gamma(link = "log"))

# type II sums of squares for overall p-values
Anova(m4, type = 2)
# Response: A
#       LR Chisq Df Pr(>Chisq)    
# site    34.611  1  4.026e-09 ***
# dista    2.428  1     0.1192 

#### 6.2   Species richness ####
#### 6.2.1 Build model ####
m5 <- lm(S ~ site * dista)

#### 6.2.2 Determine fixed components ####
m6 <- update(m5, .~. - site : dista) # remove interaction term
anova(m5, m6) # m5 is the best model

#### 6.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m5, pch = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m5) ~ site)
plot(resid(m5) ~ dista, pch = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m5)) # residuals are slightly left-skewed
qqnorm(resid(m5))
qqline(resid(m5)) # but could pass as normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity could be improved with gls

#### 6.2.4 Build model ####
require(nlme)
m7 <- gls(S ~ site * dista, weights = varIdent(form = ~1|site),
          method = "ML")

#### 6.2.5 Determine fixed components ####
m8 <- update(m7, .~. - site : dista) # remove interaction term
anova(m7, m8) # m7 fits best
m7 <- gls(S ~ site * dista, weights = varIdent(form = ~1|site),
          method = "REML")

#### 6.2.6 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m7, pch = site)
# homogeneity is slightly improved

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7)) # normality is the same

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m7 is chosen as the optimal model

#### 6.2.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m7, type = 3)
# Response: S
# Df   Chisq Pr(>Chisq)    
# (Intercept)  1 30.7065  3.002e-08 ***
# site         1 23.3985  1.317e-06 ***
# dista        1  1.2105   0.271243    
# site:dista   1 10.2210   0.001388 ** 

summary(m7)
# y = -0.060481x + 16.129944

site <- factor(site, levels = c("Zghira", "Kbira"))
m7 <- gls(S ~ site * dista, weights = varIdent(form = ~1|site),
          method = "REML")
Anova(m7, type = 3)
# Response: S
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 127.857  < 2.2e-16 ***
# site         1  23.398  1.317e-06 ***
# dista        1  29.150  6.698e-08 ***
# site:dista   1  10.221   0.001388 ** 

summary(m7)
# y = -0.31810x + 37.46297

# intersection of the two lines
eq <- rbind(c(16.129944, -0.060481), c(37.46297, -0.31810))
c(-solve(cbind(eq[,2],-1)) %*% eq[,1])
# x = 82.80843, y = 11.12161

site <- factor(site, levels = c("Kbira", "Zghira"))

#### 6.3   Shannon-Wiener diversity ####
#### 6.3.1 Build model ####
m9 <- lm(H ~ site * dista)

#### 6.3.2 Determine fixed components ####
m10 <- update(m9, .~. - site : dista) # remove interaction term
anova(m9, m10) # models are not different
# continue with m10

#### 6.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m10, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m10) ~ site)
plot(resid(m10) ~ dista, col = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m10)) 
qqnorm(resid(m10))
qqline(resid(m10)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# could be improved with gls

#### 6.3.4 Build model ####
m11 <- gls(H ~ site * dista, weights = varIdent(form = ~1|site),
           method = "ML")

#### 6.3.5 Determine fixed components ####
m12 <- update(m11, .~. - site : dista) # remove interaction term
anova(m11, m12) # models are not different
# continue with m12
m12 <- gls(H ~ site + dista, weights = varIdent(form = ~1|site),
           method = "REML")

#### 6.3.6 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
plot(m12, pch = site)
# more homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m12)) 
qqnorm(resid(m12))
qqline(resid(m12)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m12 is chosen as optimal model

#### 6.3.7 Interpret model ####
# type II sums of squares for overall p-values
Anova(m12, type = 2)
# Response: H
#       Df  Chisq Pr(>Chisq)  
# site   1 2.0451    0.15270  
# dista  1 3.0499    0.08074 .

# pairwise contrasts not required

#### 6.4   Pielou's evenness ####
#### 6.4.1 Build model ####
m13 <- lm(J ~ site * dista)

#### 6.4.2 Determine fixed components ####
m14 <- update(m13, .~. - site : dista) # remove interaction term
anova(m13, m14) # models are not different
# continue with m14

#### 6.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m14, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m14) ~ site)
plot(resid(m14) ~ dista, col = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m14))
qqnorm(resid(m14))
qqline(resid(m14)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# model heterogeneity

#### 6.4.4 Build model ####
m15 <- gls(J ~ site * dista, weights = varIdent(form = ~1|site),
           method = "ML")

#### 6.4.5 Determine fixed components ####
m16 <- update(m15, .~. - site : dista) # remove interaction term
anova(m15, m16) # m15 fits best
m15 <- gls(J ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")

#### 6.4.6 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
plot(m15, col = site)
# homogeneity is improved

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m15))
qqnorm(resid(m15))
qqline(resid(m15)) # normality is similar

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m15 is chosen as the optimal model

#### 6.4.7 Interpret model ####
# type III sums of squares for overall p-values
Anova(m15, type = 3)
# Response: J
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 212.3555  < 2.2e-16 ***
# site         1   9.3022   0.002289 ** 
# dista        1   0.0193   0.889554    
# site:dista   1   3.7290   0.053475 .

summary(m15)
# y = -0.0001512x + 0.8404082

site <- factor(site, levels = c("Zghira", "Kbira"))
m15 <- gls(J ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
Anova(m15, type = 3)
# Response: J
#              Df    Chisq Pr(>Chisq)    
# (Intercept)  1 274.5214  < 2.2e-16 ***
# site         1   9.3022  0.0022888 ** 
# dista        1  11.8229  0.0005851 ***
# site:dista   1   3.7290  0.0534746 . 

summary(m15)
# y = 0.0023240x + 0.6297496

# intersection of the two lines
eq <- rbind(c(0.8404082, -0.0001512), c(0.6297496, 0.0023240))
c(-solve(cbind(eq[,2],-1)) %*% eq[,1])
# x = 85.1077085, y = 0.8275399

site <- factor(site, levels = c("Kbira", "Zghira"))

#### 6.5   Crustacea abundance ####
#### 6.5.1 Build model ####
m17 <- lm(C ~ site * dista)

#### 6.5.2 Determine fixed components ####
m18 <- update(m17, .~. - site : dista) # remove interaction term
anova(m17, m18) # models are different
# continue with m17

#### 6.5.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m17, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m17) ~ site)
plot(resid(m17) ~ dista, col = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m17))
qqnorm(resid(m17))
qqline(resid(m17)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# could be improved with gamma

#### 6.5.4 Fit gamma distribution ####
gamma <- fitdist(C, "gamma") 
norm <- fitdist(C, "norm")

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

#### 6.5.5 Build gamma model ####
m19 <- glm(C ~ site * dista, family = Gamma(link = "log"))

#### 6.5.6 Determine fixed components ####
m20 <- update(m19, .~. - site : dista) # remove interaction term
anova(m19, m20, test = "Chisq") # m20 fits best

#### 6.5.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m20, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m20))
qqnorm(resid(m20))
qqline(resid(m20)) # normality improved

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m20 is chosen as the optimal model

#### 6.5.8 Interpret model ####
# type II sums of squares for overall p-values
Anova(m20, type = 3)
# Response: C
#       LR Chisq Df Pr(>Chisq)   
# site   10.3976  1   0.001262 **
# dista   2.1344  1   0.144031 

# pairwise contrasts not required

#### 6.6   Crustacea richness ####
#### 6.6.1 Build model ####
m21 <- lm(CS ~ site * dista)

#### 6.6.2 Determine fixed components ####
m22 <- update(m21, .~. - site : dista) # remove interaction term
anova(m21, m22) # m21 fits best

#### 6.6.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m21, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m21) ~ site)
plot(resid(m21) ~ dista, col = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m21)) # residuals are slightly right-skewed
qqnorm(resid(m21))
qqline(resid(m21))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity can be improved with gls

#### 6.6.4 Build model ####
m23 <- gls(CS ~ site * dista, weights = varIdent(form = ~1|site),
           method = "ML")

#### 6.6.5 Determine fixed components ####
m24 <- update(m23, .~. - site : dista) # remove interaction term
anova(m23, m24) # m23 fits best
m23 <- gls(CS ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")

#### 6.6.6 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
plot(m23, col = site)
# homogeneity is improved
# m23 is chosen as the optimal model

#### 6.6.7 Interpret model ####
# type III sums of squares for overall p-values
Anova(m23, type = 3)
# Response: CS
# Df   Chisq Pr(>Chisq)    
# (Intercept)  1  8.1517  0.0043022 ** 
# site         1 11.7802  0.0005986 ***
# dista        1  0.2848  0.5935944    
# site:dista   1 10.0709  0.0015063 ** 

# pairwise contrasts
site <- factor(site, levels = c("Zghira", "Kbira"))
m23 <- gls(CS ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
Anova(m23, type = 3)
# Response: CS
# Df  Chisq Pr(>Chisq)    
# (Intercept)  1 41.463  1.201e-10 ***
# site         1 11.780  0.0005986 ***
# dista        1 13.209  0.0002786 ***
# site:dista   1 10.071  0.0015063 ** 

summary(m23)
# y = -0.177280x + 17.662671

site <- factor(site, levels = c("Kbira", "Zghira"))


#### 6.7   Mollusca abundance ####
#### 6.7.1 Build model ####
m25 <- lm(M ~ site * dista)

#### 6.7.2 Determine fixed components ####
m26 <- update(m25, .~. - site : dista) # remove interaction term
anova(m25, m26) # models are different
# continue with m25

#### 6.7.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m25, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m25) ~ site)
plot(resid(m25) ~ dista, col = site)
# residual variance differs between sites
# but not with distance
# heterogeneity could be modelled by site

# test assumption of normality
hist(resid(m25))
qqnorm(resid(m25))
qqline(resid(m25)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# could be improved with gls

#### 6.7.4 Build model ####
m27 <- gls(M ~ site * dista, weights = varIdent(form = ~1|site),
           method = "ML")

#### 6.7.5 Determine fixed components ####
m28 <- update(m27, .~. - site : dista) # remove interaction term
anova(m27, m28) # m29 fits best

m27 <- gls(M ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")

#### 6.7.6 Test model fit ####
plot(m27, col = site) # homogeneity is improved
# m29 is chosen as the optimal model

#### 6.7.7 Interpret model ####
# type III sums of squares for overall p-values
Anova(m27, type = 3)
# Response: M
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 16.8802  3.981e-05 ***
# site         1  0.1587    0.69031    
# dista        1  4.2358    0.03958 *  
# site:dista   1  3.8781    0.04892 * 

site <- factor(site, levels = c("Zghira", "Kbira"))
m27 <- gls(M ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
Anova(m27, type = 3)
# Response: M
#             Df  Chisq Pr(>Chisq)  
# (Intercept)  1 2.3024    0.12918  
# site         1 0.1587    0.69031  
# dista        1 2.1325    0.14421  
# site:dista   1 3.8781    0.04892 *

site <- factor(site, levels = c("Kbira", "Zghira"))

#### 6.8   Mollusca richness ####
#### 6.8.1 Build model ####
m29 <- lm(MS ~ site * dista)

#### 6.8.2 Determine fixed components ####
m30 <- update(m29, .~. - site : dista) # remove interaction term
anova(m29, m30) # m30 fits best

#### 6.8.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m30, col = site)
# quite homogenous

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m30))
qqnorm(resid(m30))
qqline(resid(m30)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m30 is chosen as the optimal model

#### 6.8.4 Interpret model ####
# type II sums of squares for overall p-values
Anova(m30, type = 2)
# Response: MS
#           Sum Sq Df F value    Pr(>F)    
# site      84.657  1  46.298 1.254e-05 ***
# dista     37.104  1  20.292 0.0005922 ***
# Residuals 23.771 13 

# pairwise contrasts are not required


#### 7.    Univariate species data analysis ####
#### 7.1   Hippolyte inermis abundance ####
#### 7.1.1 Build model ####
m31 <- lm(Hi ~ site * dista)

#### 7.1.2 Determine fixed components ####
m32 <- update(m31, .~. - site : dista) # remove interaction term
anova(m31, m32) # m31 fits best

#### 7.1.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m31, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m31) ~ site)
plot(resid(m31) ~ dista, col = site)
# residual variance varies predictably with site and
# distance

# test assumption of normality
hist(resid(m31)) 
qqnorm(resid(m31))
qqline(resid(m31)) # residuals are quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.1.4 Fit gamma distribution ####
gamma <- fitdist(Hi, "gamma") 
norm <- fitdist(Hi, "norm")

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

#### 7.1.5 Build gamma model ####
m33 <- glm(Hi ~ site * dista, family = Gamma(link = "log"))

#### 7.1.6 Determine fixed components ####
m34 <- update(m33, .~. - site : dista) # remove interaction term
anova(m33, m34, test = "Chisq") # m33 fits best

#### 7.1.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m33, col = site)
# homogeneity is slightly better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m33))
qqnorm(resid(m33))
qqline(resid(m33)) # normality is worse (left skew)

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# try modelling heterogeneity with gls

#### 7.1.8 Model heterogeneity ####
m35 <- gls(Hi ~ site * dista,
           weights = varIdent(form = ~1|site),
           method = "REML")
m36 <- gls(Hi ~ site * dista,
           weights = varExp(form = ~dista|site),
           method = "REML")
m37 <- gls(Hi ~ site * dista,
           weights = varComb(varIdent(form = ~1|site),
                             varExp(form = ~dista|site)),
           method = "REML")
anova(m35, m36, m37) # proceed with m36
m36 <- gls(Hi ~ site * dista,
           weights = varExp(form = ~dista|site),
           method = "ML")

#### 7.1.9 Determine fixed components ####
m38 <- update(m36, .~. - site : dista) # remove interaction term
anova(m36, m38) # m36 fits best
m36 <- gls(Hi ~ site * dista,
           weights = varExp(form = ~dista|site),
           method = "REML")

#### 7.1.10 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
plot(m36, pch = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m36))
qqnorm(resid(m36))
qqline(resid(m36)) # normality is the same as m31

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m36 is chosen as the optimal model

#### 7.1.11 Interpret model ####
# type III sums of squares for overall p-values
Anova(m36, type = 3)
# Response: Hi
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 14.6516  0.0001293 ***
# site         1 33.6440  6.618e-09 ***
# dista        1  0.0292  0.8642406    
# site:dista   1 11.9803  0.0005377 ***

summary(m36)
# y = 0.01027x + 9.30709

site <- factor(site, levels = c("Zghira", "Kbira"))
m36 <- gls(Hi ~ site * dista,
           weights = varExp(form = ~dista|site),
           method = "REML")
Anova(m36, type = 3)
# Response: Hi
#             Df  Chisq Pr(>Chisq)    
# (Intercept)  1 52.533  4.231e-13 ***
# site         1 33.644  6.618e-09 ***
# dista        1 13.151  0.0002874 ***
# site:dista   1 11.980  0.0005377 ***

summary(m36)
# y = -0.58994x + 57.22562

# intersection of the two lines
eq <- rbind(c(9.30709, 0.01027), c(57.22562, -0.58994))
c(-solve(cbind(eq[,2],-1)) %*% eq[,1])
# x = 79.83627, y = 10.12701

site <- factor(site, levels = c("Kbira", "Zghira"))

#### 7.2   Jujubinus exasperatus abundance ####
#### 7.2.1 Build model ####
m39 <- lm(Je ~ site * dista)

#### 7.2.2 Determine fixed components ####
m40 <- update(m39, .~. - site : dista) # remove interaction term
anova(m39, m40) # models are not different
# continue with m40

#### 7.2.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m40, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m40) ~ site)
plot(resid(m40) ~ dista, col = site)
# residual variance varies predictably with site

# test assumption of normality
hist(resid(m40))
qqnorm(resid(m40))
qqline(resid(m40)) # residuals deviate at distribution edges
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.2.4 Fit gamma distribution ####
gamma <- fitdist(Je + 1, "gamma") 
norm <- fitdist(Je + 1, "norm")

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

#### 7.2.5 Build gamma model ####
m41 <- glm(Je + 1 ~ site * dista, family = Gamma(link = "log"))

#### 7.2.6 Determine fixed components ####
m42 <- update(m41, .~. - site : dista) # remove interaction term
anova(m41, m42, test = "Chisq") # m42 fits best

#### 7.2.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m42, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m42))
qqnorm(resid(m42))
qqline(resid(m42)) # normality is better (but one outlier)

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.2.8 Interpret model ####
# type II sums of squares for overall p-values
Anova(m42, type = 2)
# Response: Je + 1
#       LR Chisq Df Pr(>Chisq)    
# site    52.506  1  4.289e-13 ***
# dista    1.117  1     0.2906  

# no pairwise contrasts required

#### 7.3   Rissoa variabilis abundance ####
#### 7.3.1 Build model ####
m43 <- lm(Rv ~ site * dista)

#### 7.3.2 Determine fixed components ####
m44 <- update(m43, .~. - site : dista) # remove interaction term
anova(m43, m44) # continue with m43

#### 7.3.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m43, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m43) ~ site)
plot(resid(m43) ~ dista, col = site)
# residual variance varies predictably with site

# test assumption of normality
hist(resid(m43))
qqnorm(resid(m43))
qqline(resid(m43)) # residuals deviate at distribution edges
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.3.4 Fit gamma distribution ####
gamma <- fitdist(Rv + 1, "gamma") 
norm <- fitdist(Rv + 1, "norm")

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

#### 7.3.5 Build gamma model ####
m45 <- glm(Rv + 1 ~ site * dista, family = Gamma(link = "log"))

#### 7.3.6 Determine fixed components ####
m46 <- update(m45, .~. - site : dista) # remove interaction term
anova(m45, m46, test = "Chisq") # m45 fits best

#### 7.3.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m45, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m45))
qqnorm(resid(m45))
qqline(resid(m45)) # normality is better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.3.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m45, type = 3)
# Response: Rv + 1
#            LR Chisq Df Pr(>Chisq)   
# site         0.8403  1   0.359300   
# dista        1.8592  1   0.172721   
# site:dista   7.9051  1   0.004929 **

coef(m45)
# y = exp(-0.01103323*x + 1.18636843)

site <- factor(site, levels = c("Zghira", "Kbira"))
m45 <- glm(Rv + 1 ~ site * dista, 
           family = Gamma(link = "log"))
Anova(m45, type = 3)
# Response: Rv + 1
#            LR Chisq Df Pr(>Chisq)   
# site         0.8403  1   0.359300   
# dista        6.0551  1   0.013866 * 
# site:dista   7.9051  1   0.004929 **

coef(m45)
# y = exp(0.03197992*x + 0.41657740)

site <- factor(site, levels = c("Kbira", "Zghira"))


#### 7.4   Opossum shrimp (Mysida) abundance ####
#### 7.4.1 Build model ####
m47 <- lm(My ~ site * dista)

#### 7.4.2 Determine fixed components ####
m48 <- update(m47, .~. - site : dista) # remove interaction term
anova(m47, m48) # continue with m48

#### 7.4.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m48, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m48) ~ site)
plot(resid(m48) ~ dista, col = site)
# residual variance varies predictably with site

# test assumption of normality
hist(resid(m48))
qqnorm(resid(m48))
qqline(resid(m48)) # residuals deviate at distribution edges
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.4.4 Fit gamma distribution ####
gamma <- fitdist(My + .6, "gamma") 
norm <- fitdist(My + .6, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better



par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.4.5 Build gamma model ####
m49 <- glm(My + .6 ~ site * dista, family = Gamma(link = "log"))

#### 7.4.6 Determine fixed components ####
m50 <- update(m49, .~. - site : dista) # remove interaction term
anova(m49, m50, test = "Chisq") # m50 fits best

#### 7.4.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m50, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m50))
qqnorm(resid(m50))
qqline(resid(m50)) # normality is better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.4.8 Interpret model ####
# type II sums of squares for overall p-values
Anova(m50, type = 2)
# Response: My + 0.6
#       LR Chisq Df Pr(>Chisq)    
# site    53.211  1  2.996e-13 ***
# dista    0.002  1     0.9664 

# pairwise contrasts are not required


#### 7.5   Asterina gibbosa abundance ####
#### 7.5.1 Build model ####
m51 <- lm(Ag ~ site * dista)

#### 7.5.2 Determine fixed components ####
m52 <- update(m51, .~. - site : dista) # remove interaction term
anova(m51, m52) # continue with m51

#### 7.5.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m51, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m51) ~ site)
plot(resid(m51) ~ dista, col = site)
# residual variance varies predictably with site 
# and distance

# test assumption of normality
hist(resid(m51))
qqnorm(resid(m51))
qqline(resid(m51)) # residuals deviate at distribution edges
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.5.4 Fit gamma distribution ####
gamma <- fitdist(Ag + 1, "gamma") 
norm <- fitdist(Ag + 1, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.5.5 Build gamma model ####
m53 <- glm(Ag + 1 ~ site * dista, family = Gamma(link = "log"))

#### 7.5.6 Determine fixed components ####
m54 <- update(m53, .~. - site : dista) # remove interaction term
anova(m53, m54, test = "Chisq") # m53 fits best

#### 7.5.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m53, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m53))
qqnorm(resid(m53))
qqline(resid(m53)) # normality is better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.5.8 Interpret model ####
# type III sums of squares for overall p-values
Anova(m53, type = 3)
# Response: Ag + 1
#            LR Chisq Df Pr(>Chisq)    
# site         54.427  1  1.613e-13 ***
# dista         0.887  1     0.3462    
# site:dista   35.016  1  3.269e-09 ***

coef(m53)
# y = exp(-0.004219811*x + 0.314099917)

site <- factor(site, levels = c("Zghira", "Kbira"))
m53 <- glm(Ag + 1 ~ site * dista, family = Gamma(link = "log"))
Anova(m53, type = 3)
# Response: Ag + 1
#            LR Chisq Df Pr(>Chisq)    
# site         54.427  1  1.613e-13 ***
# dista        55.062  1  1.168e-13 ***
# site:dista   35.016  1  3.269e-09 ***

coef(m53)
# y = exp(-0.05500731*x + 3.72308366)

#### 7.6   Alvania discors abundance ####
#### 7.6.1 Build model ####
m55 <- lm(Ad ~ site * dista)

#### 7.6.2 Determine fixed components ####
m56 <- update(m55, .~. - site : dista) # remove interaction term
anova(m55, m56) # continue with m56

#### 7.6.3 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m56, col = site)
# heterogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
boxplot(resid(m56) ~ site)
plot(resid(m56) ~ dista, col = site)
# residual variance varies predictably with site

# test assumption of normality
hist(resid(m56))
qqnorm(resid(m56))
qqline(resid(m56)) # residuals are quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# should be improved
# -> try gamma

#### 7.6.4 Fit gamma distribution ####
gamma <- fitdist(Ad + .6, "gamma") 
norm <- fitdist(Ad + .6, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.6.5 Build gamma model ####
m57 <- glm(Ad + .6 ~ site * dista, family = Gamma(link = "log"))

#### 7.6.6 Determine fixed components ####
m58 <- update(m57, .~. - site : dista) # remove interaction term
anova(m57, m58, test = "Chisq") # m58 fits best

#### 7.6.7 Test model fit ####
# test assumption of homogeneity
# plot of residuals against fitted values
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m58, col = site)
# homogeneity is better

# test assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m58))
qqnorm(resid(m58))
qqline(resid(m58)) # normality is better

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.6.8 Interpret model ####
# type II sums of squares for overall p-values
Anova(m58, type = 2)
# Response: Ad + 0.6
#       LR Chisq Df Pr(>Chisq)  
# site     6.612  1    0.01013 *
# dista    0.007  1    0.93342 

# pairwise contrasts are not required

#### 8.    Data visualisation ####
#### 8.1   Calculate descriptive statistics ####
require(psych)
#### 8.1.1 Indices ####
Astat <- with(fauna, describeBy(A, list(Site, Distance), mat = T, digits = 10))
Astat$group2 <- as.numeric(Astat$group2)
Astat2 <- with(fauna, describeBy(A, Site, mat = T, digits = 10))

Sstat <- with(fauna, describeBy(S, list(Site, Distance), mat = T, digits = 10))
Sstat$group2 <- as.numeric(Sstat$group2)
Sstat2 <- with(fauna, describeBy(S, Site, mat = T, digits = 10))

Hstat <- with(fauna, describeBy(H, list(Site, Distance), mat = T, digits = 10))
Hstat$group2 <- as.numeric(Hstat$group2)
Hstat2 <- with(fauna, describeBy(H, Site, mat = T, digits = 10))

Jstat <- with(fauna, describeBy(J, list(Site, Distance), mat = T, digits = 10))
Jstat$group2 <- as.numeric(Jstat$group2)
Jstat2 <- with(fauna, describeBy(J, Site, mat = T, digits = 10))

Cstat <- with(fauna, describeBy(Crustacea, list(Site, Distance), mat = T, digits = 10))
Cstat$group2 <- as.numeric(Cstat$group2)
Cstat2 <- with(fauna, describeBy(Crustacea, Site, mat = T, digits = 10))

CSstat <- with(fauna, describeBy(Crustacea.S, list(Site, Distance), mat = T, digits = 10))
CSstat$group2 <- as.numeric(CSstat$group2)
CSstat2 <- with(fauna, describeBy(Crustacea.S, Site, mat = T, digits = 10))

Mstat <- with(fauna, describeBy(Mollusca, list(Site, Distance), mat = T, digits = 10))
Mstat$group2 <- as.numeric(Mstat$group2)
Mstat2 <- with(fauna, describeBy(Mollusca, Site, mat = T, digits = 10))

MSstat <- with(fauna, describeBy(Mollusca.S, list(Site, Distance), mat = T, digits = 10))
MSstat$group2 <- as.numeric(MSstat$group2)
MSstat2 <- with(fauna, describeBy(Mollusca.S, Site, mat = T, digits = 10))

#### 8.1.2 Species ####
Histat <- with(fauna, describeBy(Hippolyte.inermis, list(Site, Distance), mat = T, digits = 10))
Histat$group2 <- as.numeric(Histat$group2)
Histat2 <- with(fauna, describeBy(Hippolyte.inermis, Site, mat = T, digits = 10))

Jestat <- with(fauna, describeBy(Jujubinus.exasperatus, list(Site, Distance), mat = T, digits = 10))
Jestat$group2 <- as.numeric(Jestat$group2)
Jestat2 <- with(fauna, describeBy(Jujubinus.exasperatus, Site, mat = T, digits = 10))

Rvstat <- with(fauna, describeBy(Rissoa.variabilis, list(Site, Distance), mat = T, digits = 10))
Rvstat$group2 <- as.numeric(Rvstat$group2)
Rvstat2 <- with(fauna, describeBy(Rissoa.variabilis, Site, mat = T, digits = 10))

Mystat <- with(fauna, describeBy(Mysida.sp.A, list(Site, Distance), mat = T, digits = 10))
Mystat$group2 <- as.numeric(Mystat$group2)
Mystat2 <- with(fauna, describeBy(Mysida.sp.A, Site, mat = T, digits = 10))

Agstat <- with(fauna, describeBy(Asterina.gibbosa, list(Site, Distance), mat = T, digits = 10))
Agstat$group2 <- as.numeric(Agstat$group2)
Agstat2 <- with(fauna, describeBy(Asterina.gibbosa, Site, mat = T, digits = 10))

Adstat <- with(fauna, describeBy(Alvania.discors, list(Site, Distance), mat = T, digits = 10))
Adstat$group2 <- as.numeric(Adstat$group2)
Adstat2 <- with(fauna, describeBy(Alvania.discors, Site, mat = T, digits = 10))

#### 8.2   Calculate model predictions ####
#### 8.2.1 Re-run optimal models ####
site <- factor(site, levels = c("Zghira", "Kbira"))

m3 <- glm(A ~ site * dista, family = Gamma(link = "log"))
m7 <- gls(S ~ site * dista, weights = varIdent(form = ~1|site),
          method = "REML")
m12 <- gls(H ~ site + dista, weights = varIdent(form = ~1|site),
           method = "REML")
m15 <- gls(J ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
m20 <- glm(C ~ site + dista, family = Gamma(link = "log"))
m23 <- gls(CS ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
m27 <- gls(M ~ site * dista, weights = varIdent(form = ~1|site),
           method = "REML")
m30 <- lm(MS ~ site + dista)
m36 <- gls(Hi ~ site * dista, weights = varExp(form = ~dista|site),
           method = "REML")
m42 <- glm(Je + 1 ~ site + dista, family = Gamma(link = "log"))
m45 <- glm(Rv + 1 ~ site * dista, family = Gamma(link = "log"))
m50 <- glm(My + .6 ~ site + dista, family = Gamma(link = "log"))
m53 <- glm(Ag + 1 ~ site * dista, family = Gamma(link = "log"))
m58 <- glm(Ad + .6 ~ site + dista, family = Gamma(link = "log"))

#### 8.2.2 Generate new data frame ####
new <- data.frame(site = c(rep("Kbira", 496), rep("Zghira", 323)),
                  dista = c(seq(21.22, 70.78, by = 0.1),
                            seq(36.35, 68.63, by = 0.1)))
new$site <- as.factor(new$site)
new$site <- factor(new$site, levels = c("Zghira", "Kbira"))

#### 8.2.3 Calculate fits and 95% confidence intervals ####
inv <- family(m3)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m3, type = "link", newdata = new,
                           se.fit = TRUE))
new$Afit <- inv(pred$fit) # line fit
new$Ahi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Alo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

new$Sfit <- predict(m7, newdata = new) # line fit
modmat <-  model.matrix(formula(m7)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m7) %*% t(modmat))
new$Slo <- with(new, Sfit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$Shi <- with(new, Sfit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

new$Hfit <- predict(m12, newdata = new) # line fit
modmat <-  model.matrix(formula(m12)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m12) %*% t(modmat))
new$Hlo <- with(new, Hfit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$Hhi <- with(new, Hfit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

new$Jfit <- predict(m15, newdata = new) # line fit
modmat <-  model.matrix(formula(m15)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m15) %*% t(modmat))
new$Jlo <- with(new, Jfit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$Jhi <- with(new, Jfit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

inv <- family(m20)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m20, type = "link", newdata = new,
                           se.fit = TRUE))
new$Cfit <- inv(pred$fit) # line fit
new$Chi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Clo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

new$CSfit <- predict(m23, newdata = new) # line fit
modmat <-  model.matrix(formula(m23)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m23) %*% t(modmat))
new$CSlo <- with(new, CSfit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$CShi <- with(new, CSfit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

new$Mfit <- predict(m27, newdata = new) # line fit
modmat <-  model.matrix(formula(m27)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m27) %*% t(modmat))
new$Mlo <- with(new, Mfit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$Mhi <- with(new, Mfit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

pred <- data.frame(predict(m30, interval = "confidence", newdata = new))
new$MSfit <- pred$fit # line fit
new$MShi <- pred$upr # upper confidence interval boundary
new$MSlo <- pred$lwr # lower confidence interval boundary

new$Hifit <- predict(m36, newdata = new) # line fit
modmat <-  model.matrix(formula(m36)[-2], new) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m36) %*% t(modmat))
new$Hilo <- with(new, Hifit - qnorm(0.975)*sqrt(int)) # lower confidence interval boundary
new$Hihi <- with(new, Hifit + qnorm(0.975)*sqrt(int)) # upper confidence interval boundary

inv <- family(m42)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m42, type = "link", newdata = new,
                           se.fit = TRUE))
new$Jefit <- inv(pred$fit) # line fit
new$Jehi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Jelo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m45)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m45, type = "link", newdata = new,
                           se.fit = TRUE))
new$Rvfit <- inv(pred$fit) # line fit
new$Rvhi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Rvlo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m50)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m50, type = "link", newdata = new,
                           se.fit = TRUE))
new$Myfit <- inv(pred$fit) # line fit
new$Myhi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Mylo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m53)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m53, type = "link", newdata = new,
                           se.fit = TRUE))
new$Agfit <- inv(pred$fit) # line fit
new$Aghi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Aglo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

inv <- family(m58)$linkinv # calculate inverse of link function
pred <- data.frame(predict(m58, type = "link", newdata = new,
                           se.fit = TRUE))
new$Adfit <- inv(pred$fit) # line fit
new$Adhi <- inv(pred$fit + pred$se.fit * qnorm(0.975)) # upper confidence interval boundary
new$Adlo <- inv(pred$fit - pred$se.fit * qnorm(0.975)) # lower confidence interval boundary

#### 8.3   Define theme ####
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

#### 8.4   Plot ####
Ap <- ggplot() +
        geom_line(data = new, aes(dista, Afit, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Alo, ymax = Ahi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, A, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Astat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Astat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Abundance (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(9, 191), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = c(0.35, 0.95)) +
        mytheme
Ap # dimensions: 4.5 x 3.5

Sp <- ggplot() +
        geom_line(data = new, aes(dista, Sfit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Slo, ymax = Shi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, S, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Sstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Sstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Species richness (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(1.36, 28.64), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Sp # dimensions: 4.5 x 3.5

Hp <- ggplot() +
        geom_line(data = new, aes(dista, Hfit, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Hlo, ymax = Hhi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, H, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Hstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Hstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Diversity ("*italic(H)*" net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(1.568, 2.932), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Hp # dimensions: 4.5 x 3.5

Jp <- ggplot() +
        geom_line(data = new, aes(dista, Jfit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Jlo, ymax = Jhi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, J, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Jstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Jstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Evenness ("*italic(J)*" net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0.618, 0.982), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Jp # dimensions: 4.5 x 3.5

Cp <- ggplot() +
        geom_line(data = new, aes(dista, Cfit, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Clo, ymax = Chi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, C, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Cstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Cstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Crustacea abundance (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(5.65, 119.35), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Cp # dimensions: 4.5 x 3.5

CSp <- ggplot() +
        geom_line(data = new, aes(dista, CSfit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = CSlo, ymax = CShi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, CS, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = CSstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = CSstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression("Crustacea species richness (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0.68, 14.32), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
CSp # dimensions: 4.5 x 3.5

Mp <- ggplot() +
        geom_line(data = new, aes(dista, Mfit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Mlo, ymax = Mhi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, M, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Mstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                          fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Mstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1), guide = F) +
        ylab(expression("Mollusca abundance (net"^-1*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(breaks = seq(0, 60, by = 20)) +
        coord_cartesian(ylim = c(2.7, 57.3), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Mp # dimensions: 4.5 x 3.5

MSp <- ggplot() +
        geom_line(data = new, aes(dista, MSfit, colour = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = MSlo, ymax = MShi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, MS, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = MSstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = MSstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Mollusca richness (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0.54, 11.46), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
MSp # dimensions: 4.5 x 3.5


Hip <- ggplot() +
        geom_line(data = new, aes(dista, Hifit, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Hilo, ymax = Hihi, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, Hi, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Histat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Histat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression(italic("Hippolyte inermis")*" (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(2.27, 47.73), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Hip # dimensions: 4.5 x 3.5

Jep <- ggplot() +
        geom_line(data = new, aes(dista, Jefit - 1, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Jelo - 1, ymax = Jehi - 1, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, Je, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Jestat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Jestat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression(italic("Jujubinus exasperatus")*" (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(1.36, 28.64), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Jep # dimensions: 4.5 x 3.5

Rvp <- ggplot() +
        geom_line(data = new, aes(dista, Rvfit - 1, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Rvlo - 1, ymax = Rvhi - 1, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, Rv, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Rvstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Rvstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression(italic("Rissoa variabilis")*" (net"^-1*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 25, by = 5)) +
        coord_cartesian(ylim = c(0, 25), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Rvp # dimensions: 4.5 x 3.5

Myp <- ggplot() +
        geom_line(data = new, aes(dista, Myfit - .6, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Mylo - .6, ymax = Myhi - .6, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, My, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Mystat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Mystat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression("Mysid shrimps (net"^-1*")")) +
        xlab("Distance from source (m)") +
        coord_cartesian(ylim = c(0.9, 19.1), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Myp # dimensions: 4.5 x 3.5

Agp <- ggplot() +
        geom_line(data = new, aes(dista, Agfit - 1, colour = site, lty = site), size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Aglo - 1, ymax = Aghi - 1, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, Ag, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Agstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Agstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5), guide = F) +
        ylab(expression(italic("Asterina gibbosa")*" (net"^-1*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(expand = c(0,0), breaks = seq(0, 8, by = 2)) +
        coord_cartesian(ylim = c(0, 8), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Agp # dimensions: 4.5 x 3.5

Adp <- ggplot() +
        geom_line(data = new, aes(dista, Adfit - .6, colour = site), lty = 5, size = 0.5) +
        geom_ribbon(data = new, aes(dista, ymin = Adlo - .6, ymax = Adhi - .6, fill = site),
                    alpha = 0.5) +
        geom_point(aes(dista, Ad, colour = site), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_pointrange(data = Adstat, aes(group2, mean, ymin = mean - se, ymax = mean + se,
                                           fill = group1), size = 0.5, shape = 21, 
                        colour = rep(c("#417892","#ee850d"),8)) +
        geom_rug(data = Adstat2, aes(27.94, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f5a54a", "#6ea4be"),
                            labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#f5a54a", "#6ea4be"),
                          labels = c("Il-Ħofra ż-Żgħira", "Il-Ħofra l-Kbira"),
                          guide = guide_legend()) +
        ylab(expression(italic("Alvania discors")*" (net"^-1*")")) +
        xlab("Distance from source (m)") +
        scale_y_continuous(breaks = seq(-2, 12, by = 2)) +
        coord_cartesian(ylim = c(-1.37, 11.37), xlim = c(10, 80)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(10, 80, by = 10)) +
        theme(legend.position = "none") +
        mytheme
Adp # dimensions: 4.5 x 3.5

#### 8.5   Combine plots ####
require(cowplot)
plot1 <- plot_grid(Ap, Sp, Hp, Jp,
                   nrow = 1, labels = "auto")
plot1 # dimensions: 4.5 x 14 in

plot2 <- plot_grid(Hip, Jep, Rvp, Agp,
                   nrow = 1, labels = "auto")
plot2 # dimensions: 4.5 x 14 in

final.plot <- plot_grid(Ap + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        Sp + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        Jp + theme(axis.text.x = element_blank(),
                                   axis.title.x = element_blank()), 
                        Hip, Agp, Rvp, rel_heights = c(0.928, 1),
                        nrow = 2, align = "v", labels = "auto",
                        label_size = 15, label_fontfamily = "Helvetica Neue")
final.plot # dimensions: 9 x 10.5 in

#### 9.    Clean up ####
detach(package:vegan)
detach(package:fitdistrplus)
detach(package:car)
detach(package:nlme)
detach(package:psych)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")

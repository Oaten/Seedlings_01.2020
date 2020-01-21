 # Script to explore and analyse tree-island data.
 # Started January 2019
 # working directory set here. File paths will be relative.
setwd("D://Ecosystems_01.2020")

 # define data
input_file <- "ti_analysis_ecosystems.csv"

 # single functions
source("vif function.R")  # for another VIF function (for mixed models)
source("overdispersion function.R")  # for overdispersion function (lme4). See Ben Bolker's GLMM FAQ
source("confint function.R")  # for generating confidence interval for glmer objects for plotting. See Ben Bolker's GLMM FAQ

 # load required packages
library(ggplot2)  # plotting the spatial data
library(lme4)  # for modeling with mixed effects
library(lmerTest)  # for extracting p-values, such as they are, from lme4
library(psych)  # for pairplots
library(MuMIn)
library(piecewiseSEM)  # for SEM, and for r2 for mixed models: sem.model.fits().
library(dplyr)  # for pipes
library(car)  # for VIF function (for simple linear model's only)
library(tidyr)  # for "spread" function (for ICC)
library(lattice)  # for lattice plots (for visualising possible interactions)
library(broom)  # for better viewing of model output (coefficients etc in table form)
library(ggpubr)
library(nlme)  # for gls() function for ACF
library(aods3) # has an overdispersion-checking function, recommended by Ben Bolker

#############################
## Tree island Exploration ##-----------------------------------------------------------
#############################
 
 # Script for initial exploration of tree-island data

ti_analysis <- read.csv(input_file, header = TRUE) 
head(ti_analysis)
str(ti_analysis)
ti_analysis$plot_id <- as.factor(ti_analysis$plot_id)
ti_analysis$year <- as.factor(ti_analysis$year)
levels(ti_analysis$zone)

 # look at pairs plots for all zone-scale predictors
ti_analysis %>% 
 dplyr::select(zone, swe:htc, zone_density:zone_height) %>%
 pairs.panels()

 # look at density response
hist(ti_analysis$density)
median(ti_analysis$density)
sort(ti_analysis$density)
mean(ti_analysis$density)
var(ti_analysis$density)
dotchart(ti_analysis$density, group = ti_analysis$zone)  # greater heterogeneity in interior

 # look at percent mortality response
hist(ti_analysis$percent_mort)
sort(ti_analysis$percent_mort)
dotchart(ti_analysis$percent_mort, group = ti_analysis$zone)

 # look at germination response
hist(ti_analysis$germ)  
median(ti_analysis$germ)
mean(ti_analysis$germ)
sort(ti_analysis$germ) # lots of zeros
 # look at zeroes
x1 <- 0:10
y1 <- dpois(x1, 1.45)
plot(x1, y1) # prob of seeing a zero is 0.23
 # of zeros in data
 ti_analysis %>% count(germ == 0) # 126 zeros, 247 total
247 * 0.23  # expect 56 zeros
dotchart(ti_analysis$germ, group = ti_analysis$zone)

 # Look for outliers in predictors. 
 # No obvious outliers
dotchart(ti_analysis$swe, main = "swe", group = ti_analysis$zone)
dotchart(ti_analysis$htc, main = "htc", group = ti_analysis$zone)  # much greater variability in windward and deflation zones
dotchart(ti_analysis$zone_recruitment, main = "recruitment", group = ti_analysis$zone)
dotchart(ti_analysis$recruit_decade, main = "recruitment", group = ti_analysis$zone)
dotchart(ti_analysis$zone_age, main = "age", group = ti_analysis$zone)
dotchart(ti_analysis$zone_basal, main = "basal area", group = ti_analysis$zone)
dotchart(ti_analysis$zone_density, main = "density", group = ti_analysis$zone)
dotchart(ti_analysis$zone_height, main = "height density", group = ti_analysis$zone)

 # Check how many plots/years within plots we have
ti_analysis %>% group_by(site, zone, plot_id) %>% add_tally() %>% filter(n==3) %>% 
summarise(count = n_distinct(plot_id)) 

#-------------------------------------------------------------------------------

 # Make table of dispersion parameters for variables
head(ti_analysis)
unique(ti_analysis$year)

summary_tbl <- ti_analysis %>% 
 dplyr::select(density, germ, percent.mort = percent_mort, htc, swe, recruit.decade = recruit_decade, 
 zone.age = zone_age, zone.basal = zone_basal, zone.height = zone_height, zone.density = zone_density) %>%
  summarise_all(funs(mean, min, max, sd), na.rm = T)

str(summary_tbl)
df.stats.tidy <- summary_tbl %>% gather(stat, val) %>%
 separate(stat, into = c("var", "stat"), sep = "_") %>%
 spread(stat, val) %>% dplyr::select(var, mean, min, max, sd) %>% arrange(c(1,2,4,5,3,6,7:10))
print(df.stats.tidy)

#-------------------------------------------------------------------------------

  # Check to see how many samples (treeling plots) we have per site and zone.
ti_analysis %>% 
 count(site, zone, plot_id) %>% 
 mutate(median = median(n)) # ti11 leeward only has one plot.

#------------------------------------------------------------------------------- 
 # Calculate ICC for site/zone/plot 

 # first run unconditional model
lm_icc <- lmer(density ~ 1 + (1 | site/zone/plot_id), data = ti_analysis, REML = F)  # model is singular
summary(lm_icc) # variance of zone random effect is very small relative to other groupings

 # then get variance of random effects
icc_var <- as.data.frame(VarCorr(lm_icc))  # get variance of random effects and put into dataframe
icc_var <- icc_var %>% select(-c(var1, var2, sdcor)) %>% spread(grp, vcov) %>% # make a single row of variances 
 rename(zone = "zone:site", plot = "plot_id:(zone:site)")  # and rename the columns

 # look at ICC for plots
(icc_var$site + icc_var$zone + icc_var$plot) /
 (icc_var$site + icc_var$zone + icc_var$plot + icc_var$Residual) 
 
 # look at ICC for zones
(icc_var$site + icc_var$zone) /
 (icc_var$site + icc_var$zone + icc_var$plot + icc_var$Residual)
 
 # Look at ICC for sites
(icc_var$site) /
 (icc_var$site + icc_var$zone + icc_var$plot + icc_var$Residual) 
 
 # check to see, of the total variation between zones, how much is due to sites 
(icc_var$site) /
 (icc_var$site + icc_var$zone)  

########################
## Full model fitting ##--------------------------------------------------------------
########################

 # filter out rows where any variables are missing
ti_analysis <- ti_analysis %>% filter(!is.na(swe) | !is.na(htc))  

 # Standardize variables #
ti_analysis_std <- stdize(ti_analysis ,omit.cols = c("site", "zone", "plot_id", "year", "density", "germ", "percent_mort", "tot_prev"))

 # Now check how many plots/years within plots we have
ti_analysis %>% group_by(site, zone, plot_id) %>% add_tally() %>% filter(n==3) %>% 
summarise(count = n_distinct(plot_id)) 
 # and how many zones (= 18)
ti_analysis %>% group_by(site, zone) %>% tally()
 # how many plots present in at least one year (= 75)
ti_analysis %>% group_by(site, zone, plot_id) %>% tally() 
 # how many plots*years (= 162)
ti_analysis %>% group_by(site, zone, plot_id, year) %>% tally() 

######################
## Density response ##------------------------------------------------------------------
######################

  # Fit full model as LMM first:
lm_dens1 <- lmer(density ~ relevel(zone, ref = 4) + htc + swe + recruit_decade + zone_age + zone_basal + zone_height + zone_density + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis, na.action = na.omit, REML=F)
summary(lm_dens1)
vif.mer(lm_dens1)  # could remove height
plot(lm_dens1)  # very bad heterogeneity

 # Try full model with poisson distribution
 # did not converge; need to rescale variables
glm_dens1 <- glmer(density ~ relevel(zone, ref = 4) + htc + swe + recruit_decade + zone_age + zone_basal + zone_height + zone_density + (1|site/zone/plot_id)+ (1|year),
 data = ti_analysis, family = poisson, na.action = na.omit)

 # Try full model with poisson distribution with standardized predictor variables
 # not converging
glm_dens1 <- glmer(density ~ relevel(zone, ref = 4) + z.htc + z.swe + z.recruit_decade + z.zone_age + z.zone_basal + z.zone_height + z.zone_density + (1|site/zone/plot_id) + (1|year),
data = ti_analysis_std, family = poisson, na.action = na.omit)
 # check VIF
vif.mer(glm_dens1)  # could remove height

 # Removed height; VIF values are all <3 (except for zone)
 glm_dens2 <- glmer(density ~ relevel(zone, ref = 4) + z.htc + z.swe + z.recruit_decade + z.zone_age + z.zone_basal + z.zone_density + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std, family = poisson, na.action = na.omit)
summary(glm_dens2)
vif.mer(glm_dens2)

 # Get help with convergence
 pars_dens <- getME(glm_dens2, c("theta","fixef"))  # get parameter values
glm_dens2.restart <- update(glm_dens2, start = pars_dens) 
summary(glm_dens2.restart)  #
vif.mer(glm_dens2.restart)

#-------------------------------------------------------------------------------

 # Model validation

 # check for overdispersion
gof(glm_dens2.restart) # no evidence of overdispersion
overdisp_fun(glm_dens2.restart) # no evidence of overdispersion

 # pearson residuals vs fitted
plot(glm_dens2.restart)

 # fitted vs observed
plot(fitted(glm_dens2.restart), ti_analysis_std$density, xlab = "Predicted density", ylab="Observed density")

 # get pearson residuals
Eglm_dens <- resid(glm_dens2.restart, type = "pearson")

 # Look at residuals versus independent variables
plot(Eglm_dens ~ ti_analysis_std$z.htc)  
plot(Eglm_dens ~ ti_analysis_std$z.swe)
plot(Eglm_dens ~ ti_analysis_std$z.recruit_decade) # no evidence of a pattern in residuals
plot(Eglm_dens ~ ti_analysis_std$z.zone_basal)  
plot(Eglm_dens ~ ti_analysis_std$z.zone_age)  
plot(Eglm_dens ~ ti_analysis_std$z.zone_height)  
plot(Eglm_dens ~ ti_analysis_std$z.zone_density)  # no evidence of a pattern in residuals
plot(Eglm_dens3 ~ ti_analysis_std$year) # no evidence of a pattern in residuals
plot(Eglm_dens3 ~ ti_analysis_std$site)
plot(Eglm_dens3 ~ ti_analysis_std$zone)

 # get F-test value for quadratic term of density
glm_dens_td <- glmer(density ~ z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_density + I(z.zone_density^2) + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std, family = poisson, na.action = na.omit)
pars_dens_td <- getME(glm_dens_td, c("theta","fixef"))  # get parameter values
glm_dens4_td <- update(glm_dens_td, start = pars_dens_td) 
anova(glm_dens4_td, glm_dens2.restart)  # so overall there is no difference; accept the simpler model (eg., zone is not significant)
anova(glm_dens2.restart, glm_dens4.restart) 

 # Check for temporal autocorrelation
 # first get residuals from final model
 # warning: this step will take a long time
#dd2 <- augment(glm_dens2.restart, ti_analysis_std)  
dd2 <- readRDS("residuals", refhook = NULL) 
 
# put in form  for ACF function
 g <- gls(.resid ~ 1, data = dd2)  

 # use ACF() function
acf_dens <- ACF(g, form = ~ year|plot_id)

 # Calculate CI. 
 # CI = 1.96/sqrt(sample size - T)  # where T = time lag.
ciline <- qnorm(1-.025/2)/sqrt(76)  # use 76 (number of plots) as sample size

 # make ggplot ACF figure
ggplot(acf_dens, aes(lag, ACF)) + 
geom_segment(aes(xend = lag, yend = 0)) + 
geom_hline(aes(yintercept = 0), linetype="dashed") +
geom_hline(aes(yintercept = ciline), linetype = 3, color = 'red') + 
geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'red') +
   theme_bw()

 # get F-test value for zone
glm_dens_zone <- glmer(density ~ z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_density + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std, family = poisson, na.action = na.omit)
pars_dens_zone <- getME(glm_dens_zone, c("theta","fixef"))  # get parameter values
glm_dens_zone.restart <- update(glm_dens_zone, start = pars_dens_zone) 
anova(glm_dens_zone.restart, glm_dens2.restart)  # so overall there is no difference; accept the simpler model (eg., zone is not significant)

#-------------------------------------------------------------------------------

 # Get r-squared for final density model
rsquared(glm_dens2.restart)  
 
 # So this is the final model for density:
ranef(glm_dens2.restart)[[1]] # these are the random intercepts for plot
ranef(glm_dens2.restart)[[2]] # these are the random intercepts for zone
ranef(glm_dens2.restart)[[3]] # these are the random intercepts for ste
summary(glm_dens2.restart)
tidy(glm_dens2.restart, conf.int = TRUE, conf.method = "Wald")

##########################
## germination response ##-------------------------------------------------------
##########################

 # Try poisson. Possibly need zero-inflated model.
glm_germ1 <- glmer(germ ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_height + z.zone_density + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std, family = poisson, na.action = na.omit)
summary(glm_germ1)  # not converging
vif.mer(glm_germ1)  # height and density have high VIF

 # try again with density removed (had highest VIF)
 # didn't converge
glm_germ2 <- glmer(germ ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_height + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std, family = poisson, na.action = na.omit)
summary(glm_germ2)
vif.mer(glm_germ2)  # just zone has high VIF
 # restart
pars_germ2 <- getME(glm_germ2, c("theta","fixef"))  # get parameter values
glm_germ2.restart <- update(glm_germ2, start = pars_germ2)  # restart process with these parameter values
 # restart 
pars_germ2.restart <- getME(glm_germ2.restart, c("theta","fixef"))  # get parameter values
glm_germ2.restart2 <- update(glm_germ2.restart, start = pars_germ2.restart)  # restart process with these parameter values
summary(glm_germ2.restart2)  # singular fit. Zone error variance is close to zero.
vif.mer(glm_germ2.restart2) 

#-------------------------------------------------------------------------------

 # model validation
plot(glm_germ2.restart2)
plot(fitted(glm_germ2.restart2), ti_analysis_std$germ, xlab = "Predicted germ", ylab="Observed germ")
overdisp_fun(glm_germ2.restart2)  # no evidence of overdispersion

 # Look at residuals vs independent variables
Elm_germ2 <- resid(glm_germ2.restart2, type = "pearson")
plot(Elm_germ2 ~ ti_analysis_std$z.htc)  
plot(Elm_germ2 ~ ti_analysis_std$z.swe)
plot(Elm_germ2 ~ ti_analysis_std$z.recruit_decade)  
plot(Elm_germ2 ~ ti_analysis_std$year)  
plot(Elm_germ2 ~ ti_analysis_std$site)
plot(Elm_germ2 ~ ti_analysis_std$zone)

#-------------------------------------------------------------------------------

 # check r-squared of model
rsquared(glm_germ2.restart2) 

 # final model:
summary(glm_germ2.restart2) 
tidy(glm_germ2.restart2, conf.int = TRUE, conf.method = "Wald")

##########################################
## mortality/percent mortality response ##---------------------------------------
##########################################

 # remove  NAs and plots where there were no seedlings in prvs year
ti_analysis_mort <- ti_analysis %>% filter(!is.na(percent_mort) & tot_prev != 0) 
ti_analysis_std_mort <- ti_analysis_std %>% filter(!is.na(percent_mort) & tot_prev != 0)

 # try binomial model.

glm_mort1 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_height + z.zone_density + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
summary(glm_mort1)
vif.mer(glm_mort1) # remove zone density

glm_mort2 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_age + z.zone_basal + z.zone_height + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
vif.mer(glm_mort2) # age still high

glm_mort3 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_basal + z.zone_height + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
vif.mer(glm_mort3) # age still high

glm_mort4 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_basal + z.zone_height + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
vif.mer(glm_mort4) # both height and SWE high; but SWE is of more relevance. Remove height

glm_mort5 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.recruit_decade + z.htc + z.swe + z.zone_basal + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
vif.mer(glm_mort5) # recruitment high now
 
 # try with recruitment removed
glm_mort6 <- glmer(percent_mort ~ relevel(zone, ref = 4) + z.htc + z.swe + z.zone_basal + (1|site/zone/plot_id) + (1|year),
 data = ti_analysis_std_mort, family = binomial, weights = tot_prev, na.action = na.omit)
vif.mer(glm_mort6) # now nothing higher than 3
summary(glm_mort6)

#-------------------------------------------------------------------------------

 # model validation
plot(glm_mort6)
plot(fitted(glm_mort6), ti_analysis_std_mort$percent_mort, xlab = "Predicted mort", ylab="Observed mort")

 # Look at residuals vs independent variables
Elm_mort <- resid(glm_mort6, type = "pearson")
plot(Elm_mort ~ ti_analysis_std_mort$z.htc)  
plot(Elm_mort ~ ti_analysis_std_mort$z.swe)
plot(Elm_mort ~ ti_analysis_std_mort$z.recruit_decade)  
plot(Elm_mort ~ as.factor(ti_analysis_std_mort$year))  
plot(Elm_mort ~ ti_analysis_std_mort$site)
plot(Elm_mort ~ ti_analysis_std_mort$zone)

#-------------------------------------------------------------------------------

 # check r-squared of model
rsquared(glm_mort5) 

 # final model:
summary(glm_mort5) 
tidy(glm_mort5, conf.int = TRUE, conf.method = "Wald")

#######################
## Coefficient plots ##-------------------------------------------------------------------
#######################

 # Seedling density coefficient plot 

summary(glm_dens2.restart)
confints_dens <- as.data.frame(confint(glm_dens2.restart, method = "Wald"))
glm_dens2.restartFrame <- data.frame(var = rownames(summary(glm_dens2.restart)$coefficients),
coef = fixef(glm_dens2.restart), lowCI = confints_dens[match(rownames(summary(glm_dens2.restart)$coefficients),
 rownames(confints_dens)), 1], highCI = confints_dens[match(rownames(summary(glm_dens2.restart)$coefficients), rownames(confints_dens)), 2])
glm_dens2.restartFrame$var <- reorder(glm_dens2.restartFrame$var, glm_dens2.restartFrame$coef) #reorders
glm_dens2.restartFrame$pos <- as.factor(ifelse(glm_dens2.restartFrame$coef > 0, 1, 0))

d_coef <-ggplot(glm_dens2.restartFrame[-1,]) + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
geom_linerange(aes(x = var, ymin = lowCI, ymax = highCI), size =  0.8, colour = gray(1/2)) +
geom_point(aes(x = var, y = coef)) + theme_classic() +
theme(plot.title = element_text(face = "bold", margin = margin(b = 10), hjust = 0.5)) +
theme(axis.title = element_text(face = "bold", size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.text = element_text(size = rel(1), colour = "black")) +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
theme(axis.title.x = element_text(margin = margin(t = 5)), axis.title.y = element_text(margin = margin(r = 1)))+
labs(title = "Seedling density", x= "Predictors", y = "standardized coefficient estimate") + 
geom_text(aes(x = var, y = coef, label = sprintf("%0.2f",round(coef, digits = 2))), nudge_x = 0.3) 
d_coef <- d_coef + scale_x_discrete(labels = c("relevel(zone, ref = 4)deflation" = "zone: deflation",
"relevel(zone, ref = 4)interior" = "zone: interior", "relevel(zone, ref = 4)leeward" = "zone: leeward", 
z.recruit_decade = "recruitment", z.htc = "htc", z.swe = "swe", z.zone_age = "age", z.zone_basal = "basal", 
z.zone_density = "density")) 
d_coef <- d_coef + coord_flip() #now flip coordinates

d_coef

#-------------------------------------------------------------------------------

 # Germination coefficient plot 

summary(glm_germ2.restart2)
confints_germ <- as.data.frame(confint(glm_germ2.restart2, method = "Wald"))
glm_germ2.restart2Frame <- data.frame(var = rownames(summary(glm_germ2.restart2)$coefficients),
coef = fixef(glm_germ2.restart2), lowCI = confints_germ[match(rownames(summary(glm_germ2.restart2)$coefficients),
 rownames(confints_germ)), 1], highCI = confints_germ[match(rownames(summary(glm_germ2.restart2)$coefficients), rownames(confints_germ)), 2])
glm_germ2.restart2Frame$var <- reorder(glm_germ2.restart2Frame$var, glm_germ2.restart2Frame$coef) #hey, this works to reorder stuff!
glm_germ2.restart2Frame$pos <- as.factor(ifelse(glm_germ2.restart2Frame$coef > 0, 1, 0))

g_coef <-ggplot(glm_germ2.restart2Frame[-1,]) + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
geom_linerange(aes(x = var, ymin = lowCI, ymax = highCI), size =  0.8, colour = gray(1/2)) +
geom_point(aes(x = var, y = coef)) + theme_classic() +
theme(plot.title = element_text(face = "bold", margin = margin(b = 10), hjust = 0.5)) +
theme(axis.title = element_text(face = "bold", size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.text = element_text(size = rel(1), colour = "black")) +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
theme(axis.title.x = element_text(margin = margin(t = 5)), axis.title.y = element_text(margin = margin(r = 1)))+
labs(title = "Seedling germination", x= "Predictors", y = "standardized coefficient estimate") + 
geom_text(aes(x = var, y = coef, label = sprintf("%0.2f",round(coef, digits = 2))), nudge_x = 0.3) 
g_coef <- g_coef + scale_x_discrete(labels = c("relevel(zone, ref = 4)deflation" = "zone: deflation", "relevel(zone, ref = 4)interior" = "zone: interior", "relevel(zone, ref = 4)leeward" = "zone: leeward", z.recruit_decade = "recruitment", z.htc = "htc", z.swe = "swe", z.zone_age = "age", z.zone_basal = "basal", z.zone_height = "height")) 
g_coef <- g_coef + coord_flip() #now flip coordinates
g_coef

#-------------------------------------------------------------------------------

 # Percent mortality coefficient plot

summary(glm_mort6)
confints_mort <- as.data.frame(confint(glm_mort6, method = "Wald"))
glm_mort6Frame <- data.frame(var = rownames(summary(glm_mort6)$coefficients),
coef = fixef(glm_mort6), lowCI = confints_mort[match(rownames(summary(glm_mort6)$coefficients),
 rownames(confints_mort)), 1], highCI = confints_mort[match(rownames(summary(glm_mort6)$coefficients), rownames(confints_mort)), 2])
glm_mort6Frame$var <- reorder(glm_mort6Frame$var, glm_mort6Frame$coef) #hey, this works to reorder stuff!
glm_mort6Frame$pos <- as.factor(ifelse(glm_mort6Frame$coef > 0, 1, 0))

m_coef <-ggplot(glm_mort6Frame[-1,]) + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
geom_linerange(aes(x = var, ymin = lowCI, ymax = highCI), size =  0.8, colour = gray(1/2)) +
geom_point(aes(x = var, y = coef)) + theme_classic() +
theme(plot.title = element_text(face = "bold", margin = margin(b = 10), hjust = 0.5)) +
theme(axis.title = element_text(face = "bold", size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.text = element_text(size = rel(1), colour = "black")) +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
theme(axis.title.x = element_text(margin = margin(t = 5)), axis.title.y = element_text(margin = margin(r = 1)))+
labs(title = "Seedling percent mortality", x= "Predictors", y = "standardized coefficient estimate") + 
geom_text(aes(x = var, y = coef, label = sprintf("%0.2f",round(coef, digits = 2))), nudge_x = 0.3) 
m_coef <- m_coef + scale_x_discrete(labels = c("relevel(zone, ref = 4)deflation" = "zone: deflation", "relevel(zone, ref = 4)interior" = "zone: interior", "relevel(zone, ref = 4)leeward" = "zone: leeward", z.htc = "htc", z.swe = "swe", z.zone_basal = "basal")) 
m_coef <- m_coef + coord_flip() #now flip coordinates
m_coef


#####################
## Predicted plots ##---------------------------------------------------------
#####################
 
 # Density vs recruitment metric

summary(glm_dens2.restart)
z.recruit_decade <- seq(min(ti_analysis_std$z.recruit_decade), max(ti_analysis_std$z.recruit_decade), 0.01)  # generate sequence of recruitment
dr <- data.frame(z.recruit_decade = rep(z.recruit_decade),
zone = rep(c("deflation", "interior", "leeward", "windward"), length(z.recruit_decade)*4),
z.htc = rep(mean(ti_analysis_std$z.htc), length(z.recruit_decade)),
z.swe = rep(mean(ti_analysis_std$z.swe), length(z.recruit_decade)),
z.zone_age = rep(mean(ti_analysis_std$z.zone_age), length(z.recruit_decade)),
z.zone_basal = rep(mean(ti_analysis_std$z.zone_basal), length(z.recruit_decade)),
z.zone_density = rep(mean(ti_analysis_std$z.zone_density), length(z.recruit_decade)))
head(dr)
dr$pred <- predict(glm_dens2.restart, newdata = dr, re.form = NA, type = "response")
pr <- easyPredCI(glm_dens2.restart, newdata = dr)
pr <- as.data.frame(pr)
dr <- cbind(dr, pr)
dr <- dr %>% filter(zone == "windward")

drp <- ggplot(dr, aes(z.recruit_decade, pred), size = 1) + geom_line(size = 1)
drp <- drp + geom_ribbon(data = dr, aes(ymin = lwr, ymax = upr), alpha = 0.1)
drp <- drp + geom_point(data = ti_analysis_std, aes(z.recruit_decade, density), shape = 1) 
drp <- drp + labs(x = "Recruitment", y = "Density") + theme_classic() +
theme(axis.title = element_text(face = "bold",size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 5))) +
theme(axis.text = element_text(size = rel(1), colour = "black"))+
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
drp
 # remove legend for multipanel plot
drp <- drp + theme(legend.position="none")

#-------------------------------------------------------------------------------

 # Density vs SWE

summary(glm_dens2.restart)
z.swe <- seq(min(ti_analysis_std$z.swe), max(ti_analysis_std$z.swe), by = 0.01)  # generate sequence of recruitment
ds <- data.frame(z.swe = rep(z.swe, 4),
zone = rep(c("deflation", "interior", "leeward", "windward"), each = length(z.swe)),
z.htc = rep(mean(ti_analysis_std$z.htc), length(z.swe)*4),
z.recruit_decade = rep(mean(ti_analysis_std$z.recruit_decade), length(z.swe)*4),
z.zone_age = rep(mean(ti_analysis_std$z.zone_age), length(z.swe)*4),
z.zone_basal = rep(mean(ti_analysis_std$z.zone_basal), length(z.swe)*4),
z.zone_density = rep(mean(ti_analysis_std$z.zone_density), length(z.swe)*4))
head(ds)
ds$pred <- predict(glm_dens2.restart, newdata = ds, re.form = NA, type = "response")
ps <- easyPredCI(glm_dens2.restart, newdata = ds)
ps <- as.data.frame(ps)
ds <- cbind(ds, ps)
ds <- ds %>% filter(zone == "windward")

dsp <- ggplot(ds, aes(z.swe, pred), size = 1) + geom_line(size = 1)
dsp <- dsp + geom_ribbon(data = ds, aes(ymin = lwr, ymax = upr), alpha = 0.1)
dsp <- dsp + geom_point(data = ti_analysis_std, aes(z.swe, density), shape = 1) 
dsp <- dsp + labs(x = "SWE", y = "Density") + theme_classic() +
theme(axis.title = element_text(face = "bold",size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 5))) +
theme(axis.text = element_text(size = rel(1), colour = "black"))+
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
dsp
 # remove y-axis label for multipanel plot
dsp <- dsp +theme(axis.title.y = element_blank())

#-------------------------------------------------------------------------------

 # Density vs basal area

summary(glm_dens2.restart)
z.zone_basal <- seq(min(ti_analysis_std$z.zone_basal), max(ti_analysis_std$z.zone_basal), 0.001)  # generate sequence of recruitment. Decreased increment size here bc larger increments didn't fit
db <- data.frame(z.zone_basal = rep(z.zone_basal, 4),
zone = rep(c("deflation", "interior", "leeward", "windward"), each = length(z.zone_basal)),
z.htc = rep(mean(ti_analysis_std$z.htc), length(z.zone_basal)*4),
z.recruit_decade = rep(mean(ti_analysis_std$z.recruit_decade), length(z.zone_basal)*4),
z.zone_age = rep(mean(ti_analysis_std$z.zone_age), length(z.zone_basal)*4),
z.swe = rep(mean(ti_analysis_std$z.swe), length(z.zone_basal)*4),
z.zone_density = rep(mean(ti_analysis_std$z.zone_density), length(z.zone_basal)*4))
head(db)
db$pred <- predict(glm_dens2.restart, newdata = db, re.form = NA, type = "response")
pb <- easyPredCI(glm_dens2.restart, newdata = db)
pb <- as.data.frame(pb)
db <- cbind(db, pb)
db <- db %>% filter(zone == "windward")
 
dbp <- ggplot(db, aes(z.zone_basal, pred))  + geom_line(size = 1)
dbp <- dbp + geom_ribbon(data = db, aes(ymin = lwr, ymax = upr), alpha = 0.1)
dbp <- dbp + geom_point(data = ti_analysis_std, aes(z.zone_basal, density), shape = 1) 
dbp <- dbp + labs(x = "Basal area", y = "Density") + theme_classic() +
theme(axis.title = element_text(face = "bold",size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 5))) +
theme(axis.text = element_text(size = rel(1), colour = "black"))+
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
dbp
 # remove y-axis label for multipanel plot
dbp <- dbp +theme(axis.title.y = element_blank())

 # adjust legend position
dbp <- dbp + theme(legend.position = c(0.80, 0.97), legend.justification = c(1, 1))
dbp

#-------------------------------------------------------------------------------

 # Density vs tree density

summary(glm_dens2.restart)
z.zone_density <- seq(min(ti_analysis_std$z.zone_density), max(ti_analysis_std$z.zone_density), 0.01)  # generate sequence of recruitment
dd <- data.frame(z.zone_density = rep(z.zone_density, 4),
zone = rep(c("deflation", "interior", "leeward", "windward"), each = length(z.zone_density)),
z.htc = rep(mean(ti_analysis_std$z.htc), length(z.zone_density)*4),
z.swe = rep(mean(ti_analysis_std$z.swe), length(z.zone_density)*4),
z.zone_age = rep(mean(ti_analysis_std$z.zone_age), length(z.zone_density)*4),
z.zone_basal = rep(mean(ti_analysis_std$z.zone_basal), length(z.zone_density)*4),
z.recruit_decade = rep(mean(ti_analysis_std$z.recruit_decade), length(z.zone_density)*4))
head(dd)
dd$pred <- predict(glm_dens2.restart, newdata = dd, re.form = NA, type = "response")
pd <- easyPredCI(glm_dens2.restart, newdata = dd)
pd <- as.data.frame(pd)
dd <- cbind(dd, pd)
dd <- dd %>% filter(zone == "windward")

ddp <- ggplot(dd, aes(z.zone_density, pred))  + geom_line(size = 1)
ddp <- ddp + geom_ribbon(data = dd, aes(ymin = lwr, ymax = upr), alpha = 0.1)
ddp <- ddp + geom_point(data = ti_analysis_std, aes(z.zone_density, density), shape = 1) 
ddp <- ddp + labs(x = "Tree Density", y = "Density") + theme_classic() +
theme(axis.title = element_text(face = "bold",size = rel(1)), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
theme(axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 5))) +
theme(axis.text = element_text(size = rel(1), colour = "black"))+
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
ddp

 # remove legend for multipanel plot
ddp <- ddp + theme(legend.position="none")

#-------------------------------------------------------------------------------

 # Combine plots

densityplots <- ggarrange(drp, dbp, ddp, dsp, ncol = 2, nrow = 2, labels=c("A","B","C","D"), common.legend = TRUE, legend = "right", align = "hv")
densityplots

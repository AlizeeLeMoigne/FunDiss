# This code was used to compute the analysis of the functional dissimilarity of genomes originated from different environments (microbiomes) in the manuscript "A global Eco-evolutionary pattern of bacterial functional dissimilarity"

## Load Packages
library(pbapply)
library(vegan)
library(ggplot2)
library(reshape2)
library(plyr)
library(car)
library(ggplot2)
library(RVAideMemoire)
library(betareg)
library(car)
library(emmeans)
library(multcomp)
library(effectsize)

# Load Files: Available in Figshare https://doi.org/10.6084/m9.figshare.c.6734349 in "Habitat specificity of functional dissimilarity profiles of bacterial genomes"
read.RDS("Habitat_genomes")
read.csv("all_env.csv")
# Load file df_GTDB.txt, in Figshare in "Species-Traits matrices"
read.table("df_GTDB.txt")

#1. Calculate the absolute Bray-Curtis (BC) dissimilarity of genomes from each environment ----
# How it works: Randomizing 100 genomes from each environment separately and look at the BC distribution
set.seed(07022023) 

# Select randomly 100 genomes 999 times for each environments
random_100_GU<-NULL
for (i in 1:999) {
  random_100_GU[[i]]<-sample(GU, 100, replace=T)
}

random_100_SS<-NULL
for (i in 1:999) {
  random_100_SS[[i]]<-sample(SS, 100, replace=T)
}

random_100_WA<-NULL
for (i in 1:999) {
  random_100_WA[[i]]<-sample(WA, 100, replace=T)
}

random_100_BE<-NULL
for (i in 1:999) {
  random_100_BE[[i]]<-sample(BE, 100, replace=T)
}


# you need the file df_GTDB where the traits of the selected genomes are stored
# this is the function to calculate BC for one category
f_null_GTDB<-function (genomes,category){
  sub_df<-df_GTDB3[row.names(df_GTDB) %in% category,genomes]
  FD_BC<-as.vector(vegdist(t(sub_df), method="bray"))
  return(FD_BC)
}

# this is a function that apply the previous function to a list of categories
f_null_list<-function (liste, envi){
  catx<-pblapply(envi, f_null_GTDB ,category=liste)
  cat_flat<-as.data.frame(do.call(rbind, catx))
  return(cat_flat)
}

## Generate the thing----
# apply both function
GU_res<-lapply(list_fun_cat, FUN=f_null_list,envi=random_100_GU)
SS_res<-lapply(list_fun_cat, FUN=f_null_list,envi=random_100_SS)
WA_res<-lapply(list_fun_cat, FUN=f_null_list,envi=random_100_WA)
BE_res<-lapply(list_fun_cat, FUN=f_null_list,envi=random_100_BE)


# for each category, I have 999 observations (the randomizations) of 4950 variables (the pairwise dissimilarities (100*99/2)
# compute 1 median for each randomization
f_med<-function(x){
  rand_med<-apply(x, 1, function (r) median(r, na.rm = TRUE) )
  return(rand_med)
}

GU_med<-lapply(GU_res, FUN=f_med)
SS_med<-lapply(SS_res, FUN=f_med)
WA_med<-lapply(WA_res, FUN=f_med)
BE_med<-lapply(BE_res, FUN=f_med)

# plot
GU_p<-as.data.frame(do.call(cbind, GU_med))
m_GU_p<-melt(GU_p)

SS_p<-as.data.frame(do.call(cbind, SS_med))
m_SS_p<-melt(SS_p)

WA_p<-as.data.frame(do.call(cbind, WA_med))
m_WA_p<-melt(WA_p)

BE_p<-as.data.frame(do.call(cbind, BE_med))
m_BE_p<-melt(BE_p)


all_env<-rbind.data.frame(m_GU_p,m_SS_p,m_WA_p,m_BE_p)
all_env$envi<-c(rep("GU",31968),rep("SS",31968),rep("WA",31968),rep("BE",31968))
write.csv(all_env, "all_env.csv")
# 127872 obs of 3 variables


## Plot Figure 2 of the manuscript  ----

### Panel a
# Densities per Habitat -customized categories
a<-subset(all_env, variable=="V1")
median_all_KO<-median(a$value)

panela<-paste0("V",1:11)
sub_panela<-subset(all_env, variable %in% panela)

ggplot(sub_panela, aes(x=value))+
  geom_density(aes(fill=envi), alpha=0.8)+
  geom_vline(aes(xintercept=median_all_KO),linetype=3)+
  facet_grid(rows=vars(variable), scales="free")+
  scale_fill_manual(values=c("#dcd300","#CC79A7","#661100","#88CCEE"))+
  theme_classic()

### Panel b
# Densities per Habitat - KEGG classical categories

panelb<-paste0("V",c(1,12:32))
sub_panelb<-subset(all_env, variable %in% panelb)

ggplot(sub_panelb, aes(x=value))+
  geom_density(aes(fill=envi), alpha=0.8)+
  geom_vline(aes(xintercept=median_all_KO),linetype=3)+
  facet_grid(rows=vars(variable), scales="free")+
  scale_fill_manual(values=c("#dcd300","#CC79A7","#661100","#88CCEE"))+
  theme_classic()


# Statistics - Betaregression ----

summary(all_env$value)

#1. Although I have a 0-1 dataset, I still need to transform the values to avoid 0 and 1
all_env$bornee<-p.beta(all_env$value)

# Models
# do a separate model : one on my customized categories and one on the Kegg defined categories

## Model 1----
Vx<-paste0("V",1:32)
sub_all_env_my<-subset(all_env, variable %in% Vx[1:11])

Mx_1<-betareg(bornee~variable*envi|variable+envi, data=sub_all_env_my)
sum_Mx_1<-summary(Mx_1)
sum_Mx_1

# Link=logit
# I used variable+envi as additional factors to improve the models in terms of modelling the variance
# I compared AIC with none of the 2 variable, with 1 variable and with both. Both is the best.

#Type of estimator: ML (maximum likelihood)
#Log-likelihood: 1.391e+05 on 58 Df
#Pseudo R-squared: 0.9871
#Number of iterations: 68 (BFGS) + 2 (Fisher scoring) 
# Pseudo R2 is the global explicative capacity.

aov_Mx_1<-Anova(Mx_1)

# validation of the model
# Independence of the residuals from the fitted values
plotresid(Mx_1)
# the y~x should be horizontal: no relationship - ok

# correct modelisation of the variance of the residuals 
plotresid(Mx_1)
# residuals dispersion must be more or less constant over the x axis - ok
# link function default is logit

## THE MODEL IS VALID

# Effect size calculated as phi=racine(chisquare/number of observation)
library(effectsize)

chisq_to_phi(3976508, n=43956)
chisq_to_phi(101085, n=43956)
chisq_to_phi(72101, n=43956)

## Model 2 ----
# Mx_3
sub_all_env_kegg<-subset(all_env, variable %in% Vx[c(1,12:32)])

Mx_3<-betareg(bornee~variable*envi|variable+envi, data=sub_all_env_kegg)
sum_Mx_3<-summary(Mx_3)
sum_Mx_3

plotresid(Mx_3)
Anova(Mx_3)
# I validated the model with the same assessements as for model Mx_1

chisq_to_phi(7260069, n=87912)
chisq_to_phi(46593, n=87912)
chisq_to_phi(135020, n=87912)

## Posthocs ----
#emmeans automatically adjusts for multiple comparisons

EEM_Mx_1_variable<-emmeans(Mx_1,~variable|envi)  # Are the categories (variable) different, for each environment
EEM_Mx_1_envi<-emmeans(Mx_1,~envi|variable) # Are the environments different, for each category

EEM_Mx_3_variable<-emmeans(Mx_3,~variable|envi)
EEM_Mx_3_envi<-emmeans(Mx_3,~envi|variable)

# Model 1 - effect of categories
a<-cld(EEM_Mx_1_variable, details=T)
a
df1<-a$emmeans
df2<-a$comparisons
write.csv(df1, "emmeans_Mx_1_variable.csv")
write.csv(df2, "comp_Mx_1_variable.csv")

# Model 1 - effect of Environments
a<-cld(EEM_Mx_1_envi, details=T)
a
df1<-a$emmeans
df2<-a$comparisons
write.csv(df1, "emmeans_Mx_1_envi.csv")
write.csv(df2, "comp_Mx_1_envi.csv")

# Model 2 - effect of categories
a<-cld(EEM_Mx_3_variable, details=T)
a
df1<-a$emmeans
df2<-a$comparisons
write.csv(df1, "emmeans_Mx_3_variable.csv")
write.csv(df2, "comp_Mx_3_variable.csv")

# Model 2 - effect of Environments
a<-cld(EEM_Mx_3_envi, details=T)
a
df1<-a$emmeans
df2<-a$comparisons
write.csv(df1, "emmeans_Mx_3_envi.csv")
write.csv(df2, "comp_Mx_3_envi.csv")

# Are the means different from 0.5?
summary(EEM_Mx_1_variable,infer=T, null=0.5)

## summary stats (mean/median) ----
stat_summary<-ddply(all_env, .(variable,envi),summarize,
                    a_mean=mean(value),
                    a_median=median(value))

write.csv(stat_summary, "stat_summary_meanmedian.csv")

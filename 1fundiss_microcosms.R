# This is the code used in the manuscript "A global eco-evolutionary pattern of bacterial functional diversity"
# This code calculates the functional dissimilarity (alpha, beta, gamma) of the bacterial communities generated in identical environmental conditions.

# Upload the files: Files available in Figshare https://doi.org/10.6084/m9.figshare.c.6734349 in the "Species-Traits matrices" folder

readRDS("List_fun_cat")
read.csv("MAGs_KO_input_alpha.csv")
read.csv("MAGs_ab_input_alpha.csv")
read.csv("pMAGs_KO_input_beta.csv")
read.csv("pMAGs_ab_input_beta.csv")
read.csv("MAGs_KO_input_gamma.csv")

## Upload the modified version of rao.diversity: Source code in modified_Rao.R (available on this Github repository)

## Load Packages
library(SYNCSA)
library(reshape2)
library(ggplot2)
library(gtools)
library(pbapply)
library(plyr)
library(dplyr)
library(rstatix)
library(RColorBrewer)
library(stats)
library(adiv)
library(vegan)
library(RVAideMemoire)

## 1. Alpha Functional Dissimilarity ----
### 1.1 Calculation and plot ----

## Calculate alpha diversity for each category

## this function compiles multiple diversity indices (Gini-Simpson, Rao, functional redundancy and normalized functinal redundancy) 
f1<-function(x){
  x2<-cbind.data.frame(names(x$Simpson),x$Simpson, x$FunRao, x$FunRedundancy)
  colnames(x2)<-c("microcosm","GiniS","Rao","funRed")
  x2$nFR<-x2$funRed/x2$GiniS
  return(x2)
}

A<-MAGs_KO_input_alpha ## allKO
A<-subset(MAGs_KO_input_alpha, row.names(MAGs_KO_input_alpha) %in% aa) ## CHANGE CATEGORY
div_MAG<-rao.diversity.mod(t(MAGs_ab_input_alpha), traits=t(A))
df_aa<-f1(div_MAG)  ## Repeat this chunk for each categories and change the name of the df

meta_alpha<-rbind.data.frame(df_allKO,df_aa,df_carbohydrate,df_covit,df_GI,
                             df_other, df_regulation,df_secmet,df_substrate,df_trans,df_signal)
meta_alpha$fun<-c(rep("allKO",14),rep("aa",14),rep("carbohydrate",14),
                  rep("covit",14),rep("GI",14),rep("other",14),
                  rep("regulation",14),rep("secmet",14),rep("substrate",14),
                  rep("trans",14),rep("signal",14))

## combine all observed data
write.csv(meta_alpha, "alpha_Obs.csv")
melt_obs_alpha<-melt(meta_alpha, id=c("microcosm","fun"))

#order the categories
O_fun_red<-rev(c("allKO","signal","regulation","trans","secmet","aa","covit","GI","substrate","carbohydrate","other"))

#calculate the median of allKO
sub_allKO_Rao<-subset(melt_obs_alpha, fun=="allKO" & variable =="Rao")
med_Rao<-ddply(sub_allKO_Rao, .(variable),summarize,
               med=median(value))

# Plot the alpha functional dissimilarity
meta_alpha %>%
  mutate(fun=factor(fun, levels=O_fun_red)) %>%
  ggplot(aes(x=fun, y=Rao))+
  geom_boxplot(outlier.shape = NA, lwd=1)+
  geom_jitter(width=0.1, color="#0072B2", alpha=0.8, size=2)+
  geom_hline(data=med_Rao,aes(yintercept=med),linetype=3)+
  coord_flip()+
  theme_classic()


### 1.2 Statistical analyses on alpha functional dissimilarity ----
# mean Rao per category
aggregate(meta_alpha$Rao, by=list(meta_alpha$fun),FUN=mean)

# do categories differ from median allKO? Are there differences between functional categories?

sub_alpha_red<-subset(melt_obs_alpha,variable=="Rao")
sub_alpha_red$logit_value<-logit(sub_alpha_red$value)

# Variance eqality (Yes)
bartlett.test(logit_value~fun, sub_alpha_red)  ## homoscedascticity ## 0.998 acceptable

# Normality of the residuals (No)
res.aov_alpha<-aov(logit_value ~ fun + Error(microcosm),data=sub_alpha_red)
summary(res.aov_alpha)
A<-res.aov_alpha$Within$residuals
plot(A)
qqnorm(A)
qqline(A)
shapiro.test(A) # not good 

## Non paramteric Friedman test
res.fried.alpha<-friedman_test(data=sub_alpha_red, logit_value ~ fun|microcosm)
res.fried.alpha

effect_size<-friedman_effsize(data=sub_alpha_red,logit_value ~ fun|microcosm)
effect_size
# large effect size, Kendall's W=0.900
# W goes from 0 to 1, >0.5 = large effect

# posthoc
# Comparisons AllKO and all the other categories.
pwc <- sub_alpha_red %>%
  pairwise_wilcox_test(
    logit_value ~ fun, paired = TRUE,
    ref.group="allKO",
    p.adjust.method = "BH"
  )
pwc


## 2. Beta Functional Dissimilarity ----
### 2.1 Calculation and Plot ----

# This function organises the output: transforms distance matrix in a vector
f_trans<-function (x) {
  x[lower.tri(x,diag=T)]<-NA
  m_x<-melt(x)
  m_x<- m_x[!is.na(m_x$value),]
  return(m_x)
}

index<-c(rep("beta_red",91),rep("fun_diss",91))

# I chose Bray-Curtis due to type of matrix (Sparse (a lot of zeros))
FD_BC<-vegdist(t(pMAGs_KO_input_beta), method="bray")
hist(FD_BC) #allKO

# Calculate for each category
KO_pMAG_f<-subset(pMAGs_KO_input_beta, row.names(pMAGs_KO_input_beta) %in% aa)  ## CHANGE CATEGORY
FD_BC<-vegdist(t(KO_pMAG_f), method="bray")
hist(FD_BC) # to verify that it changes according to the categories

results_beta<-betaUniqueness(t(pMAGs_ab_input_beta),FD_BC, Nind=10000) ## WAIT

beta_red_f<-results_beta$betaRedundancy
fun_diss_f<-results_beta$DKG
#tax_diss_f<-results_beta$DR # just compute once, it is the taxonomic dissimilarity as Bray-Curtis dissimilarity

m_beta_red_f<-f_trans(beta_red_f)
m_fun_diss_f<-f_trans(fun_diss_f)
#m_tax_diss_f<-f_trans(tax_diss_f)

# CHANGE DF NAME !!
df_allKO_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_aa_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_carbohydrate_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_covit_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_GI_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_other_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_regulation_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_secmet_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_substrate_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_trans_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)
df_signal_b<-cbind.data.frame(rbind.data.frame(m_beta_red_f,m_fun_diss_f),index)

meta_beta<-rbind.data.frame(df_allKO_b,df_aa_b,df_carbohydrate_b,df_covit_b,df_GI_b,
                            df_other_b, df_regulation_b, df_secmet_b, df_substrate_b, df_trans_b, df_signal_b)

meta_beta$fun<-c(rep("allKO",182),rep("aa",182),rep("carbohydrate",182),
                 rep("covit",182),rep("GI",182),rep("other",182),
                 rep("regulation",182),rep("secmet",182),rep("substrate",182),
                 rep("trans",182),rep("signal",182))
str(meta_beta)
write.csv(meta_beta,"beta_Obs.csv")

meta_beta_FD<-subset(meta_beta, index=="fun_diss")
#median all KO
med_beta<-ddply(subset(meta_beta_FD, fun=="allKO"), .(index),summarize,
                med=median(value))

# Plot
meta_beta_FD %>%
  mutate(fun=factor(fun, levels=O_fun_red)) %>%
  ggplot(aes(x=fun, y=value))+
  geom_boxplot(outlier.shape = NA, lwd=1)+
  geom_jitter(width=0.1, color="#009E73", alpha=0.8, size=2)+
  geom_hline(data=med,aes(yintercept=med),linetype=3)+
  coord_flip()+
  theme_classic()

### 2.2 Statistical Analysis of beta Functional Dissimilarity ----
aggregate(meta_beta_FD$value, by=list(meta_beta_FD$fun),FUN=mean)

# do categories differ from median allKO? Are there differences between functional categories?

ID_pair<-paste0("pair",1:91)
meta_beta_FD$pair<-rep(ID_pair,11)
meta_beta_FD$logit_value<-logit(meta_beta_FD$value)

# Variance equality (No)
bartlett.test(logit_value~fun, meta_beta_FD)  ## homoscedascticity ## 10-16!

## Non paramteric Friedman test
res.fried<-friedman_test(data=meta_beta_FD, logit_value ~ fun|pair)
res.fried

effect_size<-friedman_effsize(data=meta_beta_FD,logit_value ~ fun|pair)
effect_size
# large effect size, Kendall's W=0.830
# W goes from 0 to 1, >0.5 = large effect

# posthoc
# I had to do it individually as the code similar to alpha was not working
wilcox.test(meta_beta_FD[1:91,7], meta_beta_FD[911:1001,7], paired=T)
p_na<-c(7.52e-14,2.2e-16,2.2e-16,2.2e-16,0.2843,2.2e-16,8.121e-16,4.722e-13,2.2e-16,2.2e-16)
p.adjust(p_na, method="BH")

## 3. Gamma Functional Dissimilarity ----
### 3.1 Calculate and Plot ----

#repeat for each category
A<-MAGs_KO_input_gamma ## allKO
A<-subset(MAGs_KO_input_gamma, row.names(MAGs_KO_input_gamma) %in% aa) ## CHANGE CATEGORY

fun_diss_gamma<-vegdist(t(A),method="bray", binary=F)
fun_diss_gamma_tp<-as.vector(fun_diss_gamma)  
df_aa<-as.data.frame(fun_diss_gamma_tp) ## CHANGE NAME!!

meta_gamma_FD<-rbind.data.frame(df_allKO,df_aa,df_carbohydrate,df_covit,df_GI,
                                      df_other, df_regulation,df_secmet,df_substrate,df_trans,df_signal)

meta_gamma_FD$fun<-c(rep("allKO",7381),rep("aa",7381),rep("carbohydrate",7381),
                           rep("covit",7381),rep("GI",7381),rep("other",7381),
                           rep("regulation",7381),rep("secmet",7381),rep("substrate",7381),
                           rep("trans",7381),rep("signal",7381))

colnames(meta_gamma_FD)<-c("fun_diss_BC", "variable")
write.csv(meta_gamma_FD, "meta_gamma.csv")

gamma_stat<-ddply(meta_gamma_FD, .(variable), summarize,
                  mean_gamma=mean(fun_diss_BC),
                  median_gamma=median(fun_diss_BC, na.rm=T),
                  sd_gamma=sd(fun_diss_BC))

#plot
median_gamma<-gamma_stat[2,3]
meta_gamma_FD %>%
  mutate(variable=factor(variable, levels=O_fun_red)) %>%
  ggplot(aes(x=variable, y=fun_diss_BC))+
  geom_boxplot(outlier.shape = NA, fill="darkorchid", lwd=1)+
  geom_hline(aes(yintercept=median_gamma),linetype=3)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()


### 3.2 Statistical Analysis Gamma Functional Dissimilarity ----
aggregate(meta_gamma_FD$fun_diss_BC, by=list(meta_gamma_FD$variable),FUN=mean)

#add pairs of genomes ID
genome_pair<-NULL
for (i in 1:7381){
  genome_pair[[i]]<-paste0("G_pair",i)
}
genome_pair<-unlist(genome_pair)
meta_gamma_FD$genome_pair<-rep(genome_pair,11)

# issue with 0-1, logit transformed create infinity because I have values equal to 0 and 1. First transform with p.beta()
meta_gamma_FD$pbeta_value<-p.beta(meta_gamma_FD$fun_diss_BC)
meta_gamma_FD$logit_value<-logit(meta_gamma_FD$pbeta_value)

# Variance equalities (No)
bartlett.test(logit_value~variable, meta_gamma_FD)  ## homoscedascticity ## 10-16!

## Non parameteric Friedman test
res.fried<-friedman_test(data=meta_gamma_FD, logit_value ~ variable|genome_pair)
res.fried

effect_size<-friedman_effsize(data=meta_gamma_FD,logit_value ~ variable|genome_pair)
effect_size
# large effect size, Kendall's W=0.643
# W goes from 0 to 1, >0.5 = large effect

pwc_gamma<-pairwise_wilcox_test(data=meta_gamma_FD,
                                logit_value ~ variable,
                                ref.group="allKO",
                                paired = TRUE,
                                p.adjust.method = "BH" )
pwc_gamma

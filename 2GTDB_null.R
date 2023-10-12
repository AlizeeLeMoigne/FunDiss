# This code was used in the manuscript "A global Eco-evolutionary pattern of bacterial functional dissimilarity"

# Load Packages
library(ggplot2)
library(vegan)
library(pbapply)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(emmeans)
library(multcomp)
library(RVAideMemoire)


## Load the file: GTDB species-Trait matrice available in Figshare https://doi.org/10.6084/m9.figshare.c.6734349, in the folder "Species-Traits matrices"
## 28 714 GTDB genomes with 8376 genes
read.table(df_GTDB.txt)
dim(df_GTDB)

## 1. Compute the null model ----
genome_ID<-colnames(df_GTDB)

set.seed(777) 
# This function randomly selects n genomes from the df_GTDB matrix
# Repeat for each subset (n=10,25,50,122,200)
random<-NULL
for (i in 1:999) {
  random[[i]]<-sample(genome_ID, 50, replace=T) # example with 50
}

# This function calculates the Bray-Curtis dissimilarity for a specific set of sampled genomes for a specific category
f_null_GTDB<-function (genomes,category){
  sub_df<-df_GTDB[row.names(df_GTDB) %in% category,genomes]
  FD_BC<-as.vector(vegdist(t(sub_df), method="bray"))
  return(FD_BC)
}

# Calculates functional dissimilarity of 999 sets of n genomes
# Repeat for each category
GTDB_FD<-pblapply(random, f_null_GTDB ,category=allKO) ## CHANGE CATEGORY
GTDB_FD_flat<-as.data.frame(do.call(rbind, GTDB_FD))

# change name +++++
write.csv(GTDB_FD_flat, "GTDB_null_allKO.csv")

# Recombine all results in one list.
signal_GTDB<-read.table("GTDB_null_signal.csv", header=T, sep=",",row.names=1)
regulation_GTDB<-read.table("GTDB_null_regulation.csv", header=T, sep=",",row.names=1)
transporters_GTDB<-read.table("GTDB_null_trans.csv", header=T, sep=",",row.names=1)
secmet_GTDB<-read.table("GTDB_null_secmet.csv", header=T, sep=",",row.names=1)
aa_GTDB<-read.table("GTDB_null_aa.csv", header=T, sep=",",row.names=1)
covit_GTDB<-read.table("GTDB_null_covit.csv", header=T, sep=",",row.names=1)
GI_GTDB<-read.table("GTDB_null_GI.csv", header=T, sep=",",row.names=1)
substrate_GTDB<-read.table("GTDB_null_substrate.csv", header=T, sep=",",row.names=1)
carbohydrate_GTDB<-read.table("GTDB_null_carbohydrate.csv", header=T, sep=",",row.names=1)
other_GTDB<-read.table("GTDB_null_other.csv", header=T, sep=",",row.names=1)
allKO_GTDB<-read.table("GTDB_null_allKO.csv", header=T, sep=",",row.names=1)

categories_GTDB<-list(allKO_GTDB,aa_GTDB,carbohydrate_GTDB,covit_GTDB,GI_GTDB,
                      other_GTDB,regulation_GTDB,secmet_GTDB,substrate_GTDB,transporters_GTDB,
                      signal_GTDB)

# Calculates the median for each set of n genomes
median_GTDB<-lapply(categories_GTDB, function(x){apply(x,1,median,na.rm=T)})

median_GTDB_flat<-as.data.frame(do.call(cbind, median_GTDB))
colnames(median_GTDB_flat)<-c("allKO","aa","carbohydrate","covit","GI","other",
                              "regulation","secmet","substrate","trans","signal")
melt_median_GTDB_flat<-melt(median_GTDB_flat)

# Each point is the median of one randomisation (999 points). 
# Each median is calculated from the pairwise comparisons of genomes of one randomization (the number of pairwise comparisons depends on the subset of genomes : (n*(n-1))/2)

# calculate the deviation from allKO for each categories and for each subset of genomes (10 to 200)
diff_allKO_sub50<-summarise(median_GTDB_flat, 
                            allKO_fundiss=median(allKO),
                            aa_diff=allKO-aa,
                            cabrbo_diff=allKO-carbohydrate,
                            covit_diff=allKO-covit,
                            GI_diff=allKO-GI,
                            other_diff=allKO-other,
                            regulation_diff=allKO-regulation,
                            secmet_diff=allKO-secmet,
                            substrate_diff=allKO-substrate,
                            trans_diff=allKO-trans,
                            signal_diff=allKO-signal)

# plot the difference with allKO
m_10<-melt(diff_allKO_sub10[,-1])
m_25<-melt(diff_allKO_sub25[,-1])
m_50<-melt(diff_allKO_sub50[,-1])
m_122<-melt(diff_allKO_sub122[,-1])
m_200<-melt(diff_allKO_sub200[,-1])

df_diff<-rbind.data.frame(m_10,m_25,m_50,m_122,m_200)
df_diff$genome_subsample<-c(rep("010",9990),rep("025",9990),rep("050",9990),
                            rep("122",9990),rep("200",9990))

# Order the categories
O_fun_red2<-c("other_diff","carbo_diff","substrate_diff","GI_diff",
              "covit_diff","aa_diff","secmet_diff","trans_diff",
              "regulation_diff","signal_diff") 

df_diff %>%
  mutate(variable=factor(variable, levels=O_fun_red2)) %>%
  ggplot(aes(x=variable, y=-value, fill=genome_subsample))+
  geom_boxplot(outlier.shape=NA, position=position_dodge())+
  scale_fill_brewer(palette="Set2")+
  labs(title="Difference between allKO and categories for various number of subsampled genomes")+
  geom_hline(yintercept=0, linetype=3)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()

## 2. Statistical analysis ----
# the variable (delta functional dissimilarity between genomes) is confined but between -1,1, not 0-1, so I cannot use p.beta()
# I had to modify the value to lock it between ]0,1[

# "Zu" is a transformation of the absolute difference with 0 which transforms the limits from ]-1,1[ to ]0,1[.
# "theo" is the median to which we want to compare, corresponding to an absolute difference value of 0 between AllKo and the other categories.
df_diff$zu<-(df_diff$value+1)/2
theo<-0+1/2

summary(df_diff$zu)

# Betaregression model
M4<-betareg(zu ~ variable + genome_subsample, data=df_diff)
plotresid(M4)
summary(M4)
Anova(M4)

# Posthoc
# eemeans
EMM<-emmeans(M4,~variable)
res_stat_GTDB<-summary(EMM, infer=T, null=theo)
res_stat_GTDB

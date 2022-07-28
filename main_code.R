
######### Library #########
library(PRROC)
library(Biobase)
library(ggsignif)
library(ggpubr)
library(caret)
library(reshape2)
library(gridExtra)
library(ggplot2); theme_set(theme_bw())
library(UpSetR)
library(GGally)
library(rpart.plot)
library(Rtsne)
library(dplyr)
library(magrittr)
library(randomForest)
library(cluster)
library(pheatmap)
library(e1071)
library(lilikoi)
library(readxl)
library(dplyr)
library(tidyr)
library(gbm)
library(reticulate)
library(stringr)
library(princurve)
library(pathview)
library(affy)
library(BiocManager)
library(limma)
library(magrittr)
library(M3C)
library(stringi)
library(coxphw)
######### Data files ####
comp_data=read_excel("Metabolon_YAG.xlsx") # initial metabolite file provided
patient_full=comp_data[,c(1:7,1201:1203)] # separate patient information from metabolite data
outcome_full = patient_full[,c(1,6)]
######### Data Wrangling ########
## Metabolites and Names Data - Raw
data_full=comp_data[,-c(2:7,1201:1203)] # remove patient information from metabolite dataset
data_full=data_full %>%
  # make in long format to match name of metabolite to COMP_ID
  pivot_longer(cols=colnames(data_full)[2:1194],
               names_to="metabolites",
               values_to="value") %>%
  # remove _3 to match COMP_ID format
  separate(metabolites,
           c("COMP_ID","Visit"),
           sep = "_") %>%
  # remove M to match COMP_ID format
  separate(COMP_ID,
           c("M","COMP_ID"),
           sep = 1)
data_full=data_full[,-c(2,4)] # remove M and visit columns
# 8 metabolites have too little variation to be analysed and need to be excluded
data_full=data_full[!(data_full$COMP_ID=="37020" | data_full$COMP_ID=="37033" | data_full$COMP_ID=="37056" | data_full$COMP_ID=="38306" | data_full$COMP_ID=="38595" | data_full$COMP_ID=="40459" | data_full$COMP_ID=="41924" | data_full$COMP_ID=="42592"),]
# remove metabolites with too little variation to be analysed
key=read_excel("Biochemical-comp_id-key.xlsx") # key for metabolite ids
key=key[!(key$COMP_ID=="37020" | key$COMP_ID=="37033" | key$COMP_ID=="37056" |
            key$COMP_ID=="38306" | key$COMP_ID=="38595" | key$COMP_ID=="40459" |
            key$COMP_ID=="41924" | key$COMP_ID=="42592"),] # remove metabolites with too little variation to be analysed
key=key[1:829,]
sum(is.na(key$HMDB_ID)) # 349 missing values
key=key[,c(2,5)] # only keep comp_id and biochemical name
data_full=merge(data_full,key,  by = "COMP_ID") # merge key and data sets to match metabolites by COMP_ID
data_full=data_full[,-1] # remove COMP_ID column
# 356 metabolites have unknown structure (name starts with X -) and need to be excluded
data_full=data_full %>%
  pivot_wider(names_from="BIOCHEMICAL",values_from="value") %>% # pivot back to wide format
  select(-starts_with("X -")) # remove metabolites with unknown structure
## Patient Data
# pre-term category: 0 = term control from the comparator group (n=297), 1 = sPTD (n=98), 2=iPTD (n=4).
data_full=merge(data_full,outcome_full,  by = "id") # merge patient and data sets by patient id
data_full = data_full %>%
  mutate(PTD_cat1 = case_when(data_full$PTD_cat == "1" ~ "1", data_full$PTD_cat == "2" ~ "0", data_full$PTD_cat == "0" ~ "0" ))
data_full = data_full %>% select(PTD_cat1, everything())  # order PTD category first
data_full = data_full %>% select(PTD_cat, everything())  # order PTD category first
rownames(data_full)=data_full$id # set patient ids as rownames
data_full=data_full[,-3] # remove id column
colnames(data_full)[1]="Label" # rename PTD to label
colnames(data_full)[2]="Label1" # rename PTD to label
data_full$Label <- as.factor(data_full$Label) # make labels factor
data_full$Label1 <- as.factor(data_full$Label1) # make labels factor
## HMDB ids
# Transform the metabolite names to the HMDB ids using Lilikoi MetaTOpathway function
dataSet=list(cmpd=c(colnames(data_full)[3:831])) # create dataSet to fit lilikoi format
convertResults=lilikoi.MetaTOpathway('name') # lilikoi function
Metabolite_pathway_table = convertResults$table # create table of metabolite pathways
# Code all missing entries as NA ("NA","",NULL --> NA)
Metabolite_pathway_table[Metabolite_pathway_table=="NA"] <- NA
Metabolite_pathway_table[Metabolite_pathway_table==""] <- NA
Metabolite_pathway_table[is.null(Metabolite_pathway_table)==TRUE] <- NA
# Missing metabolite IDs after using Lilikoi
sum(is.na(Metabolite_pathway_table$HMDB)) # 516 missing HMDB IDs for metabolites
## We want to use the data we have on the metabolites to fill some of the
# missing values we obtained from using the lilikoi package functions.
key=read_excel("Biochemical-comp_id-key.xlsx") # key for metabolite ids
colnames(key)[2]="Query" # rename to query
# This code allows us to fill in the missing data with our own data for HMDB ID
table =  Metabolite_pathway_table %>%
  left_join(key, by='Query') %>% # join Metabolite_pathway_table with the key
  mutate(HMDB_f = ifelse(is.na(HMDB),HMDB_ID,HMDB))
# Missing metabolite IDs
sum(is.na(table$HMDB_f)) # 329 missing HMDB IDs for metabolites
# Conflicting metabolite IDs
length(which(table$HMDB != table$HMDB_ID)) # 24 conflicting HMDB IDs
# Replace metabolite with conflicting HMDB IDs with our own or NAs (?)
conf_HMDB=which(table$HMDB_f != table$HMDB_ID)
for (l in c(1,2,3,4,5,6,13,18,19,21)) {
  table$HMDB_f[conf_HMDB[l]]=table$HMDB[conf_HMDB[l]]
}

for (m in c(7,8,9,10,11,12,14,15,16,17,20,22,23,24)) {
  table$HMDB_f[conf_HMDB[m]]=table$HMDB_ID[conf_HMDB[m]]
}

#table$HMDB_f[conf_HMDB]=table$HMDB_ID[conf_HMDB]
table=table[,c(1,2,4,5,6,7,23)] # only keep needed columns
table=relocate(table, HMDB_f, .after = Match) # order HMDB column
colnames(table)=colnames(Metabolite_pathway_table) # rename table with Metabolite_pathway_table
data_mom=data_full[!(data_full$Label==2),]
data_mom=data_mom[,-1]
colnames(data_mom)[1]="Label"
data_mom$Label <- as.factor(data_mom$Label) # make labels factor
metData = data_mom[,colnames(data_mom)[colnames(data_mom)!='Label']]
metClass = data_mom$Label
######### Lilikoi preprocessing (0 nzv)  #########
# Create patient data frame
patient_full$BMI=patient_full$Weight/((patient_full$Height)*0.01)^2
patient=patient_full[!(patient_full$PTD_cat==2),] # remove the 4 iatrogenic preterm deliveries
patient = patient %>% select(PTD_cat, everything())  # order PTD category first
rownames(patient)=patient$id # set patient ids as rownames
patient=patient[,-c(2,3,4,5,6,7,10)] # remove non needed columns
colnames(patient)[1]="Label" # rename PTD to label
patient$Label <- as.factor(patient$Label) # make labels factor
data_mom=data_full[!(data_full$Label==2),]
data_mom=data_mom[,-1]
colnames(data_mom)[1]="Label"
# Log Transform
data_mom_log=log(data_mom[,2:830])
data_mom_log$Label=data_mom$Label
data_mom_log = data_mom_log %>% select(Label, everything()) 

# Perform standard normalization for metabolites using the lilikoi.preproc_norm
data_mom=lilikoi.preproc_norm(inputdata=data_mom_log, method="standard")
# Perform additional preprocessing to avoid error in anova: "sum of residuals is 0"
# Exclude Zero/Near Zero Variance Features
nzv <- nearZeroVar(data_mom, saveMetrics=T)
near.zero.variance=rownames(nzv[nzv$nzv==TRUE | nzv$zeroVar==TRUE,])
data_mom = data_mom[,!(colnames(data_mom)%in%near.zero.variance)]
## 17 non zero variance terms
# Re-add Label column
data_mom$Label=patient$Label
data_mom$Label <- as.factor(data_mom$Label) # make labels factor
data_mom = data_mom %>% select(Label, everything())  # order label  column first
data_mom2=data_mom[,-1]
######### Exploratory analysis #####
lilikoi.explr(data_mom2,patient[,2:4],pca=FALSE,tsne=FALSE)
tsne(t(data_mom2),labels=data_mom$Label)
pca(t(data_mom2),labels=data_mom$Label)
data_Tsne <- Rtsne(data_mom2, check_duplicates=FALSE, pca=TRUE, perplexity=50,
                   theta=0.5, dims=2)
data_Tsne$Y %>%
  as.data.frame() %>%
  rename(tSNE1=V1,  tSNE2=V2) %>%
  mutate( Type=as.character(data_mom$Label) )%>%
  ggplot() +
  scale_color_discrete(name = "Type of Birth", labels = c("Term", "sPTB"))+
  geom_point(mapping = aes(x=tSNE1, y=tSNE2, color=Type), alpha=0.5,size=1.8) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 

patient_full=patient_full %>%
  mutate(survtime=case_when(GAwk_del>=37 ~ 37-GAwk_3,
                            GAwk_del<37 ~ TTD))
######### Density plots #####

data_mom_stand=lilikoi.preproc_norm(inputdata=data_mom, method="standard")
#data_mom_quantile=lilikoi.preproc_norm(inputdata=data_mom, method="quantile")
#data_mom_median=lilikoi.preproc_norm(inputdata=data_mom, method="median")
#data_raw_stand=lilikoi.preproc_norm(inputdata=data_raw, method="standard")
#data_raw_quantile=lilikoi.preproc_norm(inputdata=data_raw, method="quantile")
#data_raw_median=lilikoi.preproc_norm(inputdata=data_raw, method="median")


df_long <- melt(data = as.matrix(data_mom_stand[,-813]), 
                id.vars = c("patient"),
                variable.name = "variable",
                value.name = "value")

df_long <- df_long %>%
  select(-Var2) %>%
  rename(patient = Var1, expression=value)

a=ggplot(df_long, aes(x=expression, color=as.factor(patient))) +
  geom_density(show.legend=FALSE) +
  ggtitle("Standard Normalisation")+
  xlab("Expression") + ylab("Density") + theme_bw() +
  theme(plot.title = element_text(size=25,hjust = 0.5), axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=22)) +
  xlim(-5,5)
a

######### Wilcox Tests ############
patient_full=patient_full %>%
  mutate(survtime=case_when(GAwk_del>=37 ~ 37-GAwk_3,
                            GAwk_del<37 ~ TTD))
sptb=subset(patient_full,patient_full$sPTD==1)
control=subset(patient_full,patient_full$sPTD==0)
control=control[!(control$PTD_cat==2),]
# GA
quantile(sptb$GAwk_del,c(1/4,1/2,3/4))
quantile(control$GAwk_del,c(1/4,1/2,3/4))
quantile(patient_full$GAwk_del,c(1/4,1/2,3/4))
t.test(sptb$GAwk_del,control$GAwk_del)
# Height
quantile(sptb$Height,c(1/4,1/2,3/4))
quantile(control$Height,c(1/4,1/2,3/4))
quantile(patient_full$Height,c(1/4,1/2,3/4))
wilcox.test(sptb$Height,control$Height) # we used wilcox test instead of t test because not normally distirbuted
# BMI
quantile(sptb$BMI,c(1/4,1/2,3/4))
quantile(control$BMI,c(1/4,1/2,3/4))
quantile(patient_full$BMI,c(1/4,1/2,3/4))
wilcox.test(sptb$BMI,control$BMI) # we used wilcox test instead of t test because not normally distirbuted
# Age
quantile(sptb$Age,c(1/4,1/2,3/4))
quantile(control$Age,c(1/4,1/2,3/4))
quantile(patient_full$Age,c(1/4,1/2,3/4))
wilcox.test(sptb$Age,control$Age) # we used wilcox test instead of t test because not normally distirbuted
# Follow up time
quantile(sptb$survtime,c(1/4,1/2,3/4))
quantile(control$survtime,c(1/4,1/2,3/4))
quantile(patient_full$survtime,c(1/4,1/2,3/4))
wilcox.test(sptb$survtime,control$survtime)
# Correlation ()
library(corrplot)
cordata=na.omit(patient_full[,c(5,8,9,11)])
colnames(cordata)=c("sPTB", "Age", "Height","BMI")
corr = cor(cordata)
testRes = cor.mtest(cordata, conf.level = 0.95)
corrplot(corr,type = 'lower', order = 'hclust', tl.col = 'black',diag=FALSE,
         p.mat = testRes$p, insig = 'p-value',sig.level = -1,
         cl.ratio = 0.2, tl.srt = 45, col = COL2('RdBu', 10))
# can change colors tp “RdBu”, “BrBG”, “PiYG”, “PRGn”, “PuOr”, “RdYlBu”
######### Account for Confounders ########
#clinical_df contains: case/control group, confounder1,counfounder2,...., counfounderX for each sample
library(limma)
design = model.matrix(~  Height ,  data = patient)
fit = lmFit(t(data_mom2), design)
fit = eBayes(fit)
residuals = residuals(fit, t(data_mom2))
adj_expr <-residuals+matrix(apply(t(data_mom2), 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
data_mom_adj = as.data.frame(t(adj_expr))
############ Differential metabolite tests ######
## Create ExpressionSet object for all before genes together
exp=data_mom[,-1]
phen=na.omit(patient_full[,c(4,5,8,9,10)])
rownames(phen)=rownames(exp)
exp=as.matrix(exp)
eset <- new("ExpressionSet", exprs=t(exp)) # set expression values (need transpose of count matrix to fit format)
phenoData(eset) <- new("AnnotatedDataFrame", data=phen) # set phenotype data


## All metabolites together univariate
design <- model.matrix(~sPTD,data=pData(eset)) # create design matrix 

## Fit the model
fit0 <- lmFit(eset, design)
fit0 <- eBayes(fit0)

## Results
results <- decideTests(fit0[, "sPTD"],adjust.method = "BH") # wow
summary(results) # 12 upregulated
tt <- topTable(fit0, number=812, coef="sPTD", adjust="BH") 
length(which(tt$adj.P.Val<=.05))
tt$met=rownames(tt)

## All metabolites together multivariate
design <- model.matrix(~sPTD+Height,data=pData(eset)) # create design matrix 

## Fit the model
fit1 <- lmFit(eset, design)
fit1 <- eBayes(fit1)

## Results
results <- decideTests(fit1[, "sPTD"],adjust.method = "BH") # wow
summary(results) # 7 upregulated
tt1 <- topTable(fit1, number=812, coef="sPTD", adjust="BH") 
length(which(tt1$adj.P.Val<=.05))
tt1$met=rownames(tt1)


## Merge 
ttfinal=merge(tt,tt1,by="met")
ttfinal=ttfinal[,-c(2:5,7:11,13)]
colnames(ttfinal)=c("Metabolite","Univariate Adj p-Value","Multivariate Adj p-Value")
ttfinal$`Univariate Adj p-Value`=round(ttfinal$`Univariate Adj p-Value`,3)
ttfinal$`Multivariate Adj p-Value`=round(ttfinal$`Multivariate Adj p-Value`,3)

write_xlsx(ttfinal,"/Users/yasmina/Desktop/Cambridge/Thesis/code/ttfinal2.xlsx")

# QQ Plot 
ggplot(tt1,aes(x=-log(seq(0,1,length.out=812)), y=-log(sort(adj.P.Val))))+
  geom_point()+
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=22)) +
  labs(x= "-log(Expected Adjusted P-Values)", y = "-log(Observed Adjusted P-Values)")+
  geom_abline() 

# Histogram
ggplot(tt1, aes(x=P.Value)) + 
  geom_histogram(bins=40)+
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=22)) +
  labs(x= "P-values", y = "Frequency")



## All metabolites together multivariate (after already adjusting)
exp=data_mom_adj
phen=na.omit(patient_full[,c(4,5,8,9,11)])
rownames(phen)=rownames(exp)
exp=as.matrix(exp)
eset <- new("ExpressionSet", exprs=t(exp)) # set expression values (need transpose of count matrix to fit format)
phenoData(eset) <- new("AnnotatedDataFrame", data=phen) # set phenotype data

design <- model.matrix(~sPTD+Height,data=pData(eset)) # create design matrix 

## Fit the model
fit2 <- lmFit(eset, design)
fit2 <- eBayes(fit2)

## Results
results <- decideTests(fit2[, "sPTD"],adjust.method = "BH") # wow
summary(results) # 7 upregulated
tt1 <- topTable(fit2, number=812, coef="sPTD", adjust="BH") 
length(which(tt1$adj.P.Val<=.05))
tt1$met=rownames(tt1)

######## Machine Learning - Individual Metabolites ########
# Machine learning
# cvnum	- number of cross-validation folds
# dlround	- epoch number for the deep learning method
# nrun - denotes the total number of runs of each method to get their averaged performance metrics
# Rpart	- TRUE if run Rpart method
# LDA	- TRUE if run LDA method
# SVM	- TRUE if run SVM method
# RF	- TRUE if run random forest method
# GBM	- TRUE if run GBM method
# PAM	- TRUE if run PAM method
# LOG	- TRUE if run LOG method
# DL	- TRUE if run deep learning method
data_mom_adj$Label=as.factor(data_mom$Label)
data_mom_adj$Label=recode_factor(data_mom_adj$Label, `0` = "control", `1` = "sptb")
Metadata <- data_mom_adj
dataSet=list(cmpd=c(colnames(data_mom_adj)[2:813]))

MLmatrix = Metadata
measurementLabels = Metadata$Label
significantPathways = 0
trainportion = 0.7
cvnum = 10
dlround=50
nrun=100
Rpart=TRUE
LDA=TRUE
SVM=TRUE
RF=TRUE
GBM=TRUE
PAM=FALSE
LOG=TRUE
DL=FALSE
# Machine Learning (Rpart, GBM, LOG, SVM, RF, LDA)
set.seed(42)
lilikoi.machine_learning(MLmatrix = Metadata, measurementLabels = Metadata$Label,
                         significantPathways = 0, trainportion = 0.7, cvnum = 10,
                         dlround=50, nrun=100, Rpart=TRUE, LDA=FALSE, SVM=FALSE,
                         RF=FALSE, GBM=FALSE, PAM=FALSE, LOG=FALSE, DL=FALSE)

## Machine Learning (DL)
# lilikoi.machine_learning(MLmatrix = Metadata[,-813], measurementLabels = Metadata$Label,
#                         significantPathways = 0, trainportion = 0.7, cvnum = 10,
#                         dlround=50,nrun=10, Rpart=FALSE,LDA=FALSE, SVM=FALSE,
#                         GBM=FALSE, RF=FALSE,LOG=FALSE, PAM=FALSE, DL= TRUE)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)
######## Machine Learning - PDS score ########
data_mom_adj$Label=as.factor(data_mom$Label)
data_mom_adj$Label=recode_factor(data_mom_adj$Label, `0` = "control", `1` = "sptb")
Metadata <- data_mom_adj
dataSet=list(cmpd=c(colnames(data_mom)[2:813]))
PDSmatrix=lilikoi.PDSfun(table)

# Select the most signficant pathway related to phenotype.
# selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= 0.5,method="gain")

PDSmatrix_df=as.data.frame(t(PDSmatrix))
PDSmatrix_df$id=rownames(PDSmatrix_df)
PDSmatrix_df=merge(outcome_full,PDSmatrix_df,  by = "id") # merge patient and data sets by patient id
PDSmatrix_df=PDSmatrix_df[!(PDSmatrix_df$PTD_cat==2),] # remove the 4 iatrogenic preterm deliveries
rownames(PDSmatrix_df)=PDSmatrix_df$id
PDSmatrix_df=PDSmatrix_df[,-1] # remove id column
colnames(PDSmatrix_df)[1]="Label" # rename PTD to label
PDSmatrix_df=as.data.frame(sapply(PDSmatrix_df,as.numeric))
PDSmatrix_df$Label <- as.factor(PDSmatrix_df$Label) # make labels factor
PDSmatrix_df$Label=recode_factor(PDSmatrix_df$Label, `0` = "control", `1` = "sptb")


### Debug
MLmatrix = PDSmatrix_df
measurementLabels = PDSmatrix_df$Label
significantPathways = 0
trainportion = 0.7
cvnum = 10
dlround=50
nrun=100
Rpart=TRUE
LDA=TRUE
SVM=TRUE
RF=TRUE
GBM=TRUE
PAM=FALSE
LOG=TRUE
DL=FALSE

## Machine Learning (Rpart, GBM, LOG, SVM, RF,LDA)
lilikoi.machine_learning(MLmatrix = PDSmatrix_df, measurementLabels = PDSmatrix_df$Label,
                         significantPathways = 0, trainportion = 0.7, cvnum = 10,
                         nrun=100, Rpart=TRUE,LDA=TRUE, SVM=TRUE, GBM=TRUE, RF=TRUE,
                         LOG=TRUE, PAM=FALSE, DL= FALSE)

########### Testing Proportional Hazard #############
data_full=data_full[,-2]
# Log Transform
data_full_log=log(data_full[,2:830])
data_full_log$Label=data_full$Label
data_full_log = data_full_log %>% select(Label, everything()) 

# Perform standard normalization for metabolites using the lilikoi.preproc_norm
data_full=lilikoi.preproc_norm(inputdata=data_full_log, method="standard")
data_full=data_full[,-830]
data_full$Label=comp_data$PTD_cat# readd label column
data_full = data_full %>% select(Label, everything())  # order label  column first

## Zero or Near Zero Variance
nzv <- nearZeroVar(data_full, saveMetrics=T)
near.zero.variance=rownames(nzv[nzv$nzv==TRUE | nzv$zeroVar==TRUE,])
data_full = data_full[,!(colnames(data_full)%in%near.zero.variance)]

outcome_full$PTD_cat=recode(outcome_full$PTD_cat, '2'= '0','0'='0','1'='1') # iatrogenic PTBs recoded to 0 to be censored at delivery
rownames(outcome_full)=outcome_full$id
outcome_full=outcome_full[,-1]
outcome_full=as.numeric(unlist(outcome_full))

patient_full=comp_data[,c(1:7,1201:1203)] # separate patient information from metabolite data
patient_full=patient_full %>%
  mutate(GAwk_del=case_when(GAwk_del>=37 ~ 37,
                            GAwk_del<37 ~ GAwk_del))
outcome_full2=patient_full[,-c(4,5,8:10)]
outcome_full2$PTD_cat=recode(outcome_full2$PTD_cat, '2'= '0','0'='0','1'='1') # iatrogenic PTBs recoded to 0 to be censored at delivery

library(limma)
design = model.matrix(~  Height,  data = patient_full)
fit = lmFit(t(data_full[,-1]), design)
fit = eBayes(fit)
residuals = residuals(fit, t(data_full[,-1]))
adj_expr <-residuals+matrix(apply(t(data_full[,-1]), 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
data_full_adj = as.data.frame(t(adj_expr))
data_full_adj$id=rownames(data_full_adj)

data_test=merge(data_full_adj,outcome_full2,by="id")
rownames(data_test)=data_test$id
data_test$PTD_cat=as.numeric(data_test$PTD_cat)
#data_test=data_test %>%
#  mutate(survtime=GAwk_del - GAwk_3)
data_test=data_test[,-c(1,817)]
# compute schoenfeld residuals
coxfit = coxph(Surv(GAwk_3,GAwk_del,PTD_cat) ~ ., data = data_test,
               method = "breslow", robust = T, id=rownames(data_test))
sresid = resid(coxfit, type = "schoenfeld")


# compute KM estimate
sfit = survfit(Surv(GAwk_3,GAwk_del,PTD_cat) ~ 1, data = data_test) 
sest = sfit$surv[sfit$n.event > 0]
ecnt = sfit$n.event[sfit$n.event > 0]
m=4512
n=325
km = rep(sest^(m/n), ecnt)

corrtests=as.data.frame(matrix(data=NA,nrow=812,ncol=4))
rownames(corrtests)=colnames(data_full_adj[,-813])
colnames(corrtests)=c("event_time","event_time_cor","km_estimates","minimum")  
for (i in 1:812) {
  # correlation test with event time 
  corrtests[i,1]=cor.test(sort(data_test$GAwk_del[data_test$PTD_cat==1]),as.data.frame(sresid)[,i],method = "pearson")$p.value
  # correlation test with rank order of event time cor.
  corrtests[i,2]=cor.test(rank(sort(data_test$GAwk_del[data_test$PTD_cat==1])),as.data.frame(sresid)[,i], method = "pearson")$p.value
  # correlation test with KM estimates 
  corrtests[i,3]=cor.test(km,as.data.frame(sresid)[,i],method = "pearson")$p.value
  corrtests[i,4]=min(corrtests[i,1],corrtests[i,2],corrtests[i,3])
}

length(which(corrtests$minimum <= 0.05)) # 41 metabolites violate the assumption of proportional hazards
rownames(corrtests[which(corrtests$minimum <= 0.05),])

########## RFE ##########
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # repeated cv
                      repeats = 5, # number of repeats
                      number = 10) # number of folds

# Features
data_test <- data_test %>%
  select(-Comparator, -GAwk_del, -GAwk_3, -survtime) %>%
  as.data.frame()

# Target variable
y <- as.factor(data_test$PTD_cat)

data_test=data_test[,-c(1,814)]
# Training: 80%; Test: 20%
set.seed(2021)
inTrain <- createDataPartition(y, p = .80, list = FALSE)[,1]

x_train <- data_test[ inTrain, ]
x_test  <- data_test[-inTrain, ]

y_train <- y[ inTrain]
y_test  <- y[-inTrain]

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

result_rfe1 <- rfe(x = x_train, 
                   y = y_train, 
                   sizes = c(1:100),
                   rfeControl = control)

# Print the results
result_rfe1

# Print the selected features
predictors(result_rfe1)

########### Prognosis Prediction Individual Metabolites ##############

patient_full=patient_full %>%
  mutate(survtime=case_when(GAwk_del>=37 ~ 37-GAwk_3,
                            GAwk_del<37 ~ TTD))

# Prognosis model other function!

########### Testing Proportional Hazard PDS #############
data_full_adj$Label=as.factor(outcome_full)
data_full_adj$Label=recode_factor(data_full_adj$Label, `0` = "control", `1` = "sptb")

Metadata <- data_full_adj
dataSet=list(cmpd=c(colnames(data_full_adj)[1:812]))
PDSmatrix=lilikoi.PDSfun(table)


PDSmatrix_df=as.data.frame(t(PDSmatrix))
PDSmatrix_df$id=rownames(PDSmatrix_df)

data_test=merge(PDSmatrix_df,outcome_full2,by="id")
rownames(data_test)=data_full_adj$id
data_test$PTD_cat=as.numeric(data_test$PTD_cat)
data_test=data_test %>%
  mutate(survtime=GAwk_del - GAwk_3)
data_test=as.data.frame(sapply(data_test,as.numeric))
data_test$Label <- as.factor(data_test$PTD_cat) # make labels factor
data_test$Label=recode_factor(data_test$Label, `0` = "control", `1` = "sptb")

# compute schoenfeld residuals
data_test=data_test[,-c(1,254,255,256)]
coxfit = coxph(Surv(GAwk_3,GAwk_del,PTD_cat) ~ ., data = data_test, robust = T, id=rownames(data_test))
sresid = resid(coxfit, type = "schoenfeld")


# compute KM estimate
sfit = survfit(Surv(GAwk_3,GAwk_del,PTD_cat) ~ 1, data = data_test) 
sest = sfit$surv[sfit$n.event > 0]
ecnt = sfit$n.event[sfit$n.event > 0]
m=4512
n=325
km = rep(sest^(m/n), ecnt)

corrtests=as.data.frame(matrix(data=NA,nrow=249,ncol=4))
rownames(corrtests)=colnames(PDSmatrix_df[,-250])
colnames(corrtests)=c("event_time","event_time_cor","km_estimates","minimum")  
for (i in 1:249) {
  # correlation test with event time 
  corrtests[i,1]=cor.test(sort(data_test$GAwk_del[data_test$PTD_cat==1]),as.data.frame(sresid)[,i],method = "pearson")$p.value
  # correlation test with rank order of event time cor.
  corrtests[i,2]=cor.test(rank(sort(data_test$GAwk_del[data_test$PTD_cat==1])),as.data.frame(sresid)[,i], method = "pearson")$p.value
  # correlation test with KM estimates 
  corrtests[i,3]=cor.test(km,as.data.frame(sresid)[,i],method = "pearson")$p.value
  corrtests[i,4]=min(corrtests[i,1],corrtests[i,2],corrtests[i,3])
}

length(which(corrtests$minimum <= 0.05)) # 26 pathways violate the assumption of proportional hazards
rownames(corrtests[which(corrtests$minimum<= 0.05),])

########### Prognosis Prediction PDS #############
# Prognosis model other function!



############# Downstream Pathway Analysis ###########
metamat <- t(t(Metadata[, -813]))
metamat <- log2(metamat)
sampleinfo <- Metadata$Label
names(sampleinfo) <- rownames(Metadata)
grouporder <- unique(Metadata$Label)

lilikoi.KEGGplot(metamat = metamat, sampleinfo = sampleinfo, grouporder = grouporder,
                 pathid = '01100', specie = 'hsa',
                 filesuffix = 'GSE01100',
                 Metabolite_pathway_table = Metabolite_pathway_table)


############# Correlation ####
to_keep=
  c("1-docosahexaenoyl-GPE (22:6)*","1-palmitoyl-2-palmitoleoyl-GPC (16:0/16:1)*","palmitoyl-arachidonoyl-glycerol (16:0/20:4) [2]*",
    "1-dihomo-linolenoyl-GPC (20:3n3 or 6)*","1-palmitoleoyl-GPE (16:1)*","1-palmitoyl-2-palmitoleoyl-GPE (16:0/16:1)*",
    "1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*","1-myristoyl-GPC (14:0)","1-oleoyl-GPC (18:1)",
    "1-palmitoleoyl-GPC (16:1)*", "1-palmitoleoyl-GPI (16:1)*", "1-palmitoleoylglycerol (16:1)*", "1-palmitoyl-2-palmitoleoyl-GPI (16:0/16:1)*",
    "1-stearoyl-2-arachidonoyl-GPE (18:0/20:4)", "1,2-dipalmitoyl-GPC (16:0/16:0)","5alpha-pregnan-3beta,20alpha-diol disulfate",
    "alanine", "andro steroid monosulfate (1)*", "inosine","leucylalanine",
    "palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*","phytanate","1-myristoyl-2-palmitoleoyl-GPC (14:0/16:1)*",
    "1-palmitoyl-2-arachidonoyl-GPE (16:0/20:4)*", "1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*","1-palmitoyl-GPC (16:0)", "1-stearoyl-GPE (18:0)","3-methoxytyrosine",
    "4-hydroxyhippurate","adenosine","ADSGEGDFXAEGGGVR*","allantoin","arachidate (20:0)",
    "arachidonate (20:4n6)","asparagine","guanine","isocitrate","isoleucine","linoleate (18:2n6)","lysine","malate","myo-inositol","N-palmitoyl-sphingosine (d18:1/16:0)",
    "palmitoyl-oleoyl-glycerol (16:0/18:1) [1]*", "suberate (octanedioate)","tyrosine","undecanedioate")

# Correlation ()
cordata=data_full_adj[,(colnames(data_full_adj)%in%to_keep)]
library(corrplot)
corr = cor(cordata)
testRes = cor.mtest(cordata, conf.level = 0.95)

# creating correlation matrix
corr_mat <- round(cor(cordata),2)

# reorder corr matrix
# using corr coefficient as distance metric
dist <- as.dist((1-corr_mat)/2)

# hierarchical clustering the dist matrix
hc <- hclust(dist)
corr_mat <-corr_mat[hc$order, hc$order]

# reduce the size of correlation matrix
melted_corr_mat <- melt(corr_mat)
#head(melted_corr_mat)

#plotting the correlation heatmap
library(ggplot2)
ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2,fill=value)) +
  geom_tile()+
  labs(fill="Correlation")+
  scale_fill_gradient2(midpoint=0, low="#B2182B", high="#2166AC",limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1,face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(face="bold"))



df=as.data.frame(matrix(data=NA,nrow=x,ncol=1))
colnames(df)="Binary"
for (i in 1:x) {
  if (is.null(intersect(list_n, column[i]))) {
    df$Binary[i]=0
  }
  else
    df$Binary[i]=1
}
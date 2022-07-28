#' Prognosis model - Adapted to case-cohort studies
#'
#' Performs 500 permutations to identify the significance level for the log-rank test of the 
#' cross-validated KM curves for case-cohort studies
#'
#' @param data_test metabolite or pds score data, with ID, comparator group, survival time and label
#' @param test empty vector to store 500 test statistics
#' @import reticulate survminer
#' @importFrom stats complete.cases median predict quantile
#' @importFrom utils write.table
#' @importFrom glmnet cv.glmnet
#' @importFrom survival Surv survfit survdiff coxph
#' @importFrom survey 
#' @importFrom surveyCV 
#' @importFrom jskm 


permutation.cv.km.ccs <- function(data_test,   k=10){


data_test=data_test %>%
  mutate(survtime=GAwk_del-GAwk_3)

k=10
set.seed(2)
CEsubset1 <- data_test[which(data_test$Comparator == 1), ]
CEsubset1$.foldID <- folds.svy(CEsubset1, nfolds = k, clusterID = "id")

set.seed(2)
CEsubset2 <- data_test[which(data_test$Comparator == 0), ]
CEsubset2$.foldID <- folds.svy(CEsubset2, nfolds = k, clusterID = "id")

CEsubset <- rbind(CEsubset1,CEsubset2)

ff=1
coxFeatures30 =character(0)
for (ff in 1:k) {
  
  met_train <- subset(CEsubset, .foldID != ff)
  met_test  <- subset(CEsubset, .foldID == ff)
  
  time_test=met_test$survtime
  event_test=met_test$PTD_cat
  
  dcch <- twophase(id=list(~id,~id), strata=list(NULL,~PTD_cat),
                   subset=~I(Comparator | PTD_cat), data=met_train)
  
  fit = svycoxph(Surv(survtime,PTD_cat)~`1-dihomo-linolenoyl-GPC (20:3n3 or 6)*` + 
                   `1-palmitoyl-2-palmitoleoyl-GPC (16:0/16:1)*` + 
                   `palmitoyl-arachidonoyl-glycerol (16:0/20:4) [2]*` +
                   `1-palmitoleoyl-GPC (16:1)*`+
                   `1-palmitoleoyl-GPE (16:1)*`+
                   `3-methoxytyrosine`+
                   `4-acetamidophenol`+
                   `alanine`+
                   `allantoin`+
                   `arachidate (20:0)`+
                   `arachidonate (20:4n6)`+
                   `beta-hydroxyisovalerate`+
                   `inosine`+
                   `isocitrate`+
                   `isoleucine`+
                   `linoleate (18:2n6)`+
                   `malate`+
                   `myo-inositol`+
                   `taurodeoxycholate`+
                   `threonine`+
                   `1-palmitoleoylglycerol (16:1)*`+
                   `1-docosahexaenoyl-GPE (22:6)*`+
                   `1-palmitoyl-2-palmitoleoyl-GPI (16:0/16:1)*`+
                   `1-palmitoleoyl-GPI (16:1)*`+
                   `1-arachidonoyl-GPE (20:4n6)*`+
                   `1,2-dipalmitoyl-GPC (16:0/16:0)`+
                   `13-HODE + 9-HODE`+
                   `1-(1-enyl-stearoyl)-2-docosahexaenoyl-GPE (P-18:0/22:6)*`+
                   `1-stearoyl-2-arachidonoyl-GPE (18:0/20:4)`+
                   `palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*`+
                   `1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*`+
                   `ADSGEGDFXAEGGGVR*`+
                   `1-dihomo-linoleoylglycerol (20:2)`+
                   `5alpha-pregnan-3beta,20alpha-diol disulfate`+
                   `N4-acetylcytidine`+
                   `undecanedioate`+
                   `1-palmitoyl-GPE (16:0)`+
                   `1-palmitoyl-2-palmitoleoyl-GPE (16:0/16:1)*`+
                   `1-palmitoyl-2-docosahexaenoyl-GPE (16:0/22:6)*`+
                   `1-stearoyl-GPE (18:0)`+ 
                   `1-palmitoyl-2-oleoyl-GPE (16:0/18:1)`+
                   `1-palmitoyl-2-oleoyl-GPC (16:0/18:1)`+
                   `palmitoyl-oleoyl-glycerol (16:0/18:1) [1]*`+
                   `1-palmitoyl-GPC (16:0)`+
                   `1-myristoyl-GPC (14:0)`+
                   `1-palmitoyl-2-linoleoyl-GPE (16:0/18:2)`+
                   `palmitoloelycholine`+
                   `1-palmitoyl-2-arachidonoyl-GPE (16:0/20:4)*`+
                   `1-palmitoyl-2-oleoyl-GPI (16:0/18:1)*`+
                   `1-myristoyl-2-linoleoyl-GPC (14:0/18:2)*`+
                   `sebacate (decanedioate)`+
                   `asparagine`+
                   `1-myristoyl-2-palmitoleoyl-GPC (14:0/16:1)*`+
                   `2-myristoyl-GPC (14:0)*`+
                   `1-palmitoyl-2-arachidonoyl-GPI (16:0/20:4)*`+
                   `palmitoyl-arachidonoyl-glycerol (16:0/20:4) [1]*`+
                   `1-palmitoyl-2-dihomo-linolenoyl-GPC (16:0/20:3n3 or 6)*`+
                   `1-linoleoyl-GPE (18:2)*`+
                   `1-palmitoylglycerol (16:0)`+
                   `1-stearoyl-2-oleoyl-GPE (18:0/18:1)`+
                   `N-palmitoyl-sphingosine (d18:1/16:0)`+
                   `1-myristoylglycerol (14:0)`+
                   `1-stearoyl-2-linoleoyl-GPE (18:0/18:2)*`+
                   `1-linolenoyl-GPC (18:3)*`+
                   `guanine`+
                   `estriol 3-sulfate`+
                   `N-acetyl-3-methylhistidine*`+
                   `phenylalanylphenylalanine`+
                   `3-hydroxy-3-methylglutarate`+
                   `4-hydroxyhippurate`+
                   `S-allylcysteine`+
                   `phosphatidylcholine (15:0/18:1, 17:0/16:1)*`+
                   `kynurenine`+
                   `4-androsten-3beta,17beta-diol disulfate (2)`+
                   `2-hydroxy-3-methylvalerate`+
                   `1-oleoyl-GPC (18:1)`+
                   `valylleucine`+
                   `I-urobilinogen`+
                   `1-linolenoylglycerol (18:3)`+
                   `lanthionine`+
                   `5-hydroxyhexanoate`+
                   `leucylalanine`+
                   `glycylvaline`+
                   `adenosine`+
                   `3-hydroxysebacate`+
                   `1-(1-enyl-stearoyl)-2-linoleoyl-GPC (P-18:0/18:2)*`+
                   `gamma-carboxyglutamate`+
                   `sphingosine 1-phosphate`+
                   `1-methylxanthine`+
                   `1-palmitoyl-2-linoleoyl-GPI (16:0/18:2)`+
                   `propyl 4-hydroxybenzoate sulfate`+
                   `bilirubin (E,E)*`+
                   `linoleoyl-arachidonoyl-glycerol (18:2/20:4) [2]*`+
                   `dihydroferulic acid`+
                   `1-nonadecanoyl-GPC (19:0)`+
                   `cystine`+
                   `ribose`+
                   `pimeloylcarnitine/3-methyladipoylcarnitine`+
                   `homoarginine`+
                   `behenoyl sphingomyelin (d18:1/22:0)*`+
                   `serylalanine`+
                   `1-linoleoyl-GPC (18:2)`+
                   `1-arachidoyl-2-docosahexaenoyl-GPC (20:0/22:6)*`+
                   `2-oleoyl-GPE (18:1)*`+
                   `andro steroid monosulfate (1)*`+
                   `phenylalanylleucine`+
                   `16a-hydroxy DHEA 3-sulfate`+
                   `kynurenate`+
                   `1-stearoyl-2-oleoyl-GPG (18:0/18:1)`+
                   `3-methylhistidine`+
                   `1-palmitoleoylglycerol (16:1)*`+
                   `2-palmitoleoyl-GPC (16:1)*`+
                   `N-methylproline`+
                   `N6,N6,N6-trimethyllysine`+
                   `phytanate`+
                   `2-palmitoylglycerol (16:0)`+
                   `4-ethylphenylsulfate`+
                   `glutamine`+
                   `isoleucylleucine/leucylisoleucine`+
                   `progesterone`+
                   `oleoyl-oleoyl-glycerol (18:1/18:1) [2]*`+
                   `5alpha-androstan-3beta,17beta-diol disulfate`+
                   `N-acetylalanine`+
                   `theobromine`+
                   `hyocholate`+
                   `2-acetamidophenol sulfate`+
                   `1-eicosapentaenoyl-GPE (20:5)*`+
                   `cysteine-glutathione disulfide`+
                   `1-oleoyl-2-docosahexaenoyl-GPE (18:1/22:6)*`+
                   `valylarginine`+
                   `1-linoleoyl-2-linolenoyl-GPC (18:2/18:3)*`+
                   `tartronate (hydroxymalonate)`+
                   `prolylproline`+
                   `pantothenate`+
                   `tetradecanedioate`+
                   `caproate (6:0)`+
                   `1-palmityl-2-stearoyl-GPC (O-16:0/18:0)*`+
                   `palmitoleate (16:1n7)`+
                   `1-(1-enyl-palmitoyl)-2-docosahexaenoyl-GPE (P-16:0/22:6)*`+
                   `glycochenodeoxycholate glucuronide (1)`+
                   `1-palmityl-GPE (O-16:0)*`,
                 design=dcch)
  
  score <- predict(fit, newdata = met_test[,-c(1,814,815,816,817,818)],type="risk")
  PIscore <- as.data.frame(score, drop = F)
  colnames(PIscore)  <- "score"
  PIscore$EVENT <- event_test
  PIscore$TIME <- time_test
  
  coxphFeat=as.data.frame(sort(fit$coefficients,decreasing=TRUE))
  coxphFeat$met=str_remove_all(rownames(coxphFeat),"`")
  colnames(coxphFeat)=c("coeff","met")
  coxphFeatures=head(coxphFeat[order(-coxphFeat$coeff),]$met,30)
  coxFeatures30 = c(coxFeatures30,coxphFeatures)
  
  ### for 4 risk groups
  PI = quantile(PIscore$score,probs=c(.25,.5,.75))
  for(i in 1:nrow(PIscore)){
    if (PIscore[i,"score"]>PI[[3]])
      PIscore[i,"Risk Quartile"]<-1
    
    else if (PIscore[i,"score"]>PI[[2]] && PIscore[i,"score"]<=PI[[3]] )
      PIscore[i,"Risk Quartile"]<-2
    
    else if (PIscore[i,"score"]>PI[[1]] && PIscore[i,"score"]<=PI[[2]] )
      PIscore[i,"Risk Quartile"]<-3
    
    
    else if(PIscore[i,"score"]<=PI[[1]])
      PIscore[i,"Risk Quartile"]<-4
  }
  
  ### for 2 risk groups
  #PI = quantile(PIscore$score,probs=.8)
  #for(i in 1:nrow(PIscore)){
  #  if (PIscore[i,"score"]>PI)
  #    PIscore[i,"Risk Quartile"]<-1
    
  #  else if (PIscore[i,"score"]<PI )
  #    PIscore[i,"Risk Quartile"]<-2
  #}
  
  if (ff==1) {
    PIfinal=PIscore}
  if  (ff!=1) {
    PIfinal=rbind(PIfinal,PIscore)}
  remove(PIscore)
}

PIscore_na <- PIfinal[complete.cases(PIfinal),]

# With svykm 
PIscore_na$id=rownames(PIscore_na)
fin=merge(CEsubset,PIscore_na,by="id")
dpbc <- twophase(id=list(~id,~id), strata=list(NULL,~PTD_cat),
                 subset=~I(Comparator | PTD_cat), data=fin,method="simple")

svykm=svykm(Surv(TIME, EVENT) ~ `Risk Quartile`, design=dpbc,se=TRUE)
newcox=svycoxph(Surv(TIME, EVENT) ~ `Risk Quartile`,design=dpbc)
AIC=2*sum(diag(solve(newcox$inv.info, newcox$var)))

svydifftest=svylogrank(Surv(TIME, EVENT) ~ `Risk Quartile`,design=dpbc)

survp=svyjskm(svykm,pval = TRUE,
        pval.size = 10, 
        pval.coord = c(2,.05),
        pval.testname = TRUE,
        legend = TRUE,
        ci=T,
        legendposition=c(.9,.15))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25)) 



coxFeatures30 =as.data.frame(sort(table(coxFeatures30),decreasing=TRUE)[1:30])[,1]
coxFeatures30 =as.character(coxFeatures30)

return(summary(newcox),svydifftest,survp,coxFeatures30)
}




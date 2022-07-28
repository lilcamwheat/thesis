####### Plotting ########

# AUC
output2=as.data.frame(PDS_ML_results[["output2"]])
colnames(output2)=c("Rpart","LDA","SVM","RF","GBM","LOG")
auc=pivot_longer(output2,cols=c(1:6), names_to = "model", values_to = "auc")
kruskal.test(auc ~ model, data = auc)
pairwise.wilcox.test(auc$auc, auc$model,
                     p.adjust.method = "BH")

# Sens
output4=as.data.frame(PDS_ML_results[["output4"]])
colnames(output4)=c("Rpart","LDA","SVM","RF","GBM","LOG")
sens=pivot_longer(output4,cols=c(1:6), names_to = "model", values_to = "sens")
kruskal.test(sens ~ model, data = sens)
pairwise.wilcox.test(sens$sens, sens$model,
                     p.adjust.method = "BH")

# Spec
output6=as.data.frame(PDS_ML_results[["output6"]])
colnames(output6)=c("Rpart","LDA","SVM","RF","GBM","LOG")
spec=pivot_longer(output6,cols=c(1:6), names_to = "model", values_to = "spec")
kruskal.test(spec ~ model, data = spec)
pairwise.wilcox.test(spec$spec, spec$model,
                     p.adjust.method = "BH")

# F1
output8=as.data.frame(PDS_ML_results[["output8"]])
colnames(output8)=c("Rpart","LDA","SVM","RF","GBM","LOG")
f1=pivot_longer(output8,cols=c(1:6), names_to = "model", values_to = "f1")
kruskal.test(f1 ~ model, data = f1)
pairwise.wilcox.test(f1$f1, f1$model,
                     p.adjust.method = "BH")

# BA
output10=as.data.frame(PDS_ML_results[["output10"]])
colnames(output10)=c("Rpart","LDA","SVM","RF","GBM","LOG")
BA=pivot_longer(output10,cols=c(1:6), names_to = "model", values_to = "BA")
kruskal.test(BA ~ model, data = BA)
pairwise.wilcox.test(BA$BA, BA$model,
                     p.adjust.method = "BH")

## ALL
all=cbind(auc,sens$sens,spec$spec,f1$f1,BA$BA)
colnames(all)=c("model","AUC","SENS","SPEC","F1","Balanced Accuracy")
all=pivot_longer(all,cols=c(2:6), names_to = "metric", values_to = "value")

all.summary <- all %>%
  group_by(model, metric) %>%
  summarise(
    sd = sd(value),
    value = mean(value))
all.summary

all.summary = all.summary %>%
  mutate(ymin=value-sd) %>%
  mutate(ymax=value+sd) %>%
  mutate(ymin=ifelse(ymin < 0, 0,ymin)) %>%
  mutate(ymax=ifelse(ymax > 1, 1,ymax))


all.summary$model <- factor(all.summary$model, levels = c('LOG', 'Rpart', 'RF','GBM','LDA','SVM'))

## Plot with p-values
ggplot(data = all.summary,
       aes(x = metric, y = value, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.2,
                position=position_dodge(.9)) +
  xlab("") + ylab("")  + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),legend.text=element_text(size=22)) +
  labs(fill = "") +
  scale_fill_manual(values = c("#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFF2AE"))+
  geom_signif(y_position = c(0.70, 0.75,0.8,0.85,0.9,
                             0.65, 0.70,0.75,0.8,0.85,
                             0.2,0.25,0.3,0.35,0.4,
                             0.5,0.55,0.6,0.65,0.7,
                             0.35,0.4,0.45,0.5,0.55), 
              xmin = c(0.62, 0.62, 0.62, 0.62,0.62,
                       1.62,1.62,1.62,1.62,1.62,
                       2.62,2.62,2.62,2.62,2.62,
                       3.62,3.62,3.62,3.62,3.62,
                       4.62,4.62,4.62,4.62,4.62), 
              xmax = c(0.78, 0.94,1.1,1.26,1.4,
                       1.78, 1.94,2.1,2.25,2.41,
                       2.78, 2.94,3.1,3.25,3.41,
                       3.78, 3.94,4.1,4.25,4.41,
                       4.78, 4.94,5.1,5.25,5.41),
              annotation = c("***", "ns","*","ns","ns",
                             "ns", "ns","ns","ns","ns",
                             "***", "***","***","*","***",
                             "***", "***","***","*","***",
                             "***", "***","***","ns","***"), tip_length = 0,textsize=5) 



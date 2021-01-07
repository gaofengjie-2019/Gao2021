library(cowplot)
library("tidyverse")
library("ggplot2")
library("ggsci")
library("ggpubr")

mydata<-data %>%
  gather(key="scale", value="score", PHQ15:IESR) %>%
  dplyr::select(SampleID, scale, score, everything())
head(mydata)
p<-ggboxplot(mydata, x = "Group", y= "score",
             color = "Group", palette = "jama",
             add = "jitter")
p+stat_compare_means()
p1<-ggboxplot(data, x = "group", y= "PHQ15",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="t.test")
p2<-ggboxplot(data, x = "group", y= "PHQ9",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="t.test")
p3<-ggboxplot(data, x = "group", y= "PSQI",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="t.test")
p4<-ggboxplot(data, x = "group", y= "GAD7",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="wilcox.test")
p5<-ggboxplot(data, x = "group", y= "SCL90",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="wilcox.test")
p6<-ggboxplot(data, x = "group", y= "IESR",
             color = "group", palette = "jama",
             add = "jitter")+
  stat_compare_means(method="wilcox.test")

p+stat_compare_means(method = "t.test", label = "p.signif",comparisons = my_comparisons)

my_comparisons<-list(c("V2","V3"),c("V1","V2"),c("V1","V3"))
ggboxplot(mydata, x="scale", y= "score",
          color="Group", add="jitter", palette = "jama")+
  stat_compare_means(method="wilcox.test", lable= "p.signif", comparisons = my_comparisons)

compare_means(score~Group, data=mydata, Group.by = "scale")
p<-ggboxplot(mydata,x = "Group", y= "score",
             color = "Group", palette = "jama",
             add = "jitter",
             facet.by = "scale", short.panel.labs = FALSE)
p+stat_compare_means(method="wilcox.test", lable= "p.signif", comparisons = my_comparisons)

p1<-ggboxplot(data, x="Group", y="IESR", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif") # Add pairwise comparisons p-value
p2<-ggboxplot(data, x="Group", y="PHQ9", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif") 
p3<-ggboxplot(data, x="Group", y="GAD7", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif")
p4<-ggboxplot(data, x="Group", y="PSQI", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif") 
p5<-ggboxplot(data, x="Group", y="PHQ15", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif") 
p6<-ggboxplot(data, x="Group", y="SCL90", color = "Group",palette = "jama", add="jitter",width=0.5)+
  stat_compare_means(comparisons=my_comparisons,step.increase= 0.15,label="p.signif") 
ggarrange(p6,p5,p4,p3,p2, p1, ncol=3,nrow=2,labels=c("a","b","c","d","e","f"))
ggsave("test4.pdf",width=118,height=178,units=c("mm"))

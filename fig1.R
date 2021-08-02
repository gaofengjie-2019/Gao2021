##Figure 1: Boxplot
library(reshape2)
library(ggpubr)
meta = read.csv("metadata.csv", row.names = 1)
meta1 = meta[meta$Group %in% c('Day0','Control'), ]
cbbPalette <- c("#87CEFA","#FFAEB9")
p <- ggplot(meta1, aes(x=Group, y=PHQ15, fill=Group))+#PHQ15 can be replaced by other factors
  stat_boxplot(geom='errorbar', width=0.15)+
  geom_boxplot(aes(fill=Group), outlier.size=0, outlier.alpha = 0)+
  scale_fill_manual(values = cbbPalette) +
  geom_jitter(colour = 'black',size=1, alpha=0.8)+
  #geom_jitter(colour='grey90',size=1.5)+
  stat_compare_means(method="wilcox.test", aes(label = paste0("p = ", ..p.format..)))+
  #scale_fill_brewer(palette='Set1')+
  theme_bw()+
  theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour='black'))
p
ggsave(paste0(prefile, 'PHQ15-boxplot.pdf'), p, width = 4, height = 5)#prefile is the the file that can save the picture

##Figure 1: Compute alpha diversity
library(vegan)
library(dplyr)
asv<-read.csv("asv_6174.csv",row.names = 1)
asv<-t(asv) %>% as.data.frame() #asv in col, asv is reads
richness<-estimateR(asv, index=c("chao", "Species"))
shannon<-diversity(asv, index="shannon")
simpson<-diversity(asv, index= "simpson")
data<-rbind(richness, shannon, simpson)
data<-t(data) %>% as.data.frame()
write.csv(data,"alpha.csv")

##Figure 1: Beta diversity
library(vegan)
library(dplyr)
library(ggplot2)
asv<-read.csv("asv_6174.csv",row.names = 1)
group<-read.csv("metadata.csv",row.names = 1)
asv<-asv/9655 #9655 is the number of reads that we sampled randomly from each sample
asv3<-t(asv) %>% as.data.frame()#sample in row
bray_dis <- vegdist(asv3, method = 'bray')
bray_dis <- as.matrix(bray_dis)
#write.csv(bray_dis,"bray_curtis.csv")
beta<-bray_dis
sub_design= subset(group, group$Group %in% c("Day0","Control"))
sub_design$group<-sub_design$Group
sub_beta=beta[rownames(sub_design),rownames(sub_design)]
# k is dimension, 3 is recommended; eig is eigenvalues
pcoa = cmdscale(sub_beta, k=3, eig=T)
# get coordinate string, format to dataframme
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
# rename group name
#levels(sub_design$group)=c("Baseline","Control")
mycolors<-c("#87CEFA","#FFAEB9")
points = cbind(points, sub_design$group)
colnames(points) = c("PC1", "PC2", "PC3","group") 
points$group=factor(points$group,levels = c("Day0","Control"))
p = ggplot(points, aes(x=PC1, y=PC2, color=group)) + geom_point(alpha=.7, size=2) +
  scale_color_manual(values=mycolors)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  stat_ellipse(level=0.95)+
  theme(plot.margin=unit(rep(0.5,4), 'lines'),
        panel.background=element_rect(fill='transparent', color='black'),
        panel.border=element_rect(fill='transparent', color='transparent'),
        panel.grid=element_blank())
p
ggsave("bray_curtis2.pdf", p, width = 5, height = 3)


cbbPalette <- c("#87CEFA","#FFAEB9")
p<-ggplot(points, aes(x=group, y=PC1)) + 
  geom_boxplot(aes(colour=group, fill= group),
               
               # custom boxes
               color=cbbPalette,
               fill=cbbPalette,
               alpha=0.5,
               
               # Notch?
               notch=TRUE,
               notchwidth = 0.2,
               
               outlier.alpha = 0
  )+
  #coord_flip()+
  expand_limits(y=c(-0.4, 0.4))+
  theme(#plot.margin=unit(rep(0.5,4), 'lines'),
    panel.background=element_rect(fill='transparent', color='black'),
    panel.border=element_rect(fill='transparent', color='transparent'),
    panel.grid=element_blank())+
  stat_compare_means(label = "p.format")


##Figure 4: Random forest regression
##day0 against deltascale(180-45)
genus <- read.csv("species.csv", row.names = 1)
meta <- read.csv("metadata.csv", row.names = 1)
genus <- t(genus) %>% as.data.frame()
meta <- meta[meta$CompleteSeries=="CS", ]
meta0 <- meta[meta$Group=="Day0", ]
meta45 <- meta[meta$Group=="Day45", ]
meta180 <- meta[meta$Group=="Day180", ]
genus0 <- genus[rownames(meta0), ]
meta45 <- meta45[,c(30:31)]
meta180 <- meta180[,c(30:31)]
deltascale <- meta180-meta45
prefile<-"D:/stress/analysis/Figure4/"
#regression 
i=2
train.x <- genus0
train.y <- deltascale[,i]
cv.fold <- 5
cv.step <- 0.9
cv.time <- 5
marker.num <- 0
prefix <- colnames(deltascale)[i]

# crossvalidation
pdf.dir <- paste0(prefile, prefix, "_randomForest.pdf")
pdf(pdf.dir, width = 20, height = 7)
par(mfrow = c(1, 2))
set.seed(0)
train.cv <- replicate(cv.time, rfcv1(train.x, train.y, cv.fold = cv.fold, step = cv.step), simplify = F)

#mse
error.cv <- sapply(train.cv, "[[", "error.cv")
error.cv.rm <- rowMeans(error.cv)
id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
error.cv[id, ]
marker.num=37
if (marker.num == 0) {
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
}
matplot(train.cv[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cv.time), main =  paste("select", marker.num, "Vars"), xlab = "Number of vars", 
        ylab = "CV mse", lty = 1)
lines(train.cv[[1]]$n.var, error.cv.rm, lwd = 2)
abline(v = marker.num, col = "pink", lwd = 2)

#nmse
nmse.cv <- sapply(train.cv, "[[", "nmse")
nmse.cv.rm <- rowMeans(nmse.cv)
id <- nmse.cv.rm < min(nmse.cv.rm) + sd(nmse.cv.rm)
nmse.cv[id, ]
if (marker.num == 0) {
  marker.num <- min(as.numeric(names(nmse.cv.rm)[id]))
}
matplot(train.cv[[1]]$n.var, nmse.cv, type = "l", log = "x", col = rep(1, cv.time), main = paste("select", marker.num, "Vars"), xlab = "Number of vars", 
        ylab = "CV nmse", lty = 1)
lines(train.cv[[1]]$n.var, nmse.cv.rm, lwd = 2)
abline(v = marker.num, col = "pink", lwd = 2)

dev.off()

# pick marker by corossvalidation
marker.t <- table(unlist(lapply(train.cv, function(x) {
  lapply(x$res, "[", 1:marker.num)
})))
marker.t <- sort(marker.t, d = T)
names(marker.t) <- colnames(train.x)[as.numeric(names(marker.t))]
marker.p <- names(marker.t)[1:marker.num]

# train model
set.seed(0)
train.rf <- randomForest(train.x[, marker.p], train.y, importance = T)
imp.dir <- paste0(prefile, prefix, "_marker.imp.csv")
write.csv(train.rf$importance, imp.dir)

train.p <- predict(train.rf)
R2 <- R_squa2(train.p, train.y)
write.table(R2,file = paste0(prefile, prefix,"R2.txt"))

##Figure 3: Boxplot
library(ggpubr)
library(ggsci)
library(reshape2)
prefile<-"F:/stress/analysis/Figure4/"
meta = read.csv("metadata.csv", row.names = 1)
species = read.csv("species_r.csv", row.names = 1)
species[species==0] =  1e-7
species = log10(species)
meta0 = meta[meta$Group %in% c("Day0", "Control"), ]
data0 = species[,rownames(meta0)]
data0 = t(data0) %>% as.data.frame()
data = data0[,c("Bacteroides eggerthii" , "[Eubacterium] hallii group_uncultured bacterium", "Intestinibacter bartlettii DSM 16795", "Stenotrophobacter_uncultured Acidobacteria bacterium")]#the differential species were obtained from Find_dif.R
data$Group = meta0$Group

data0 = melt(data, id.vars=c("Group"), measure.vars=colnames(data)[1:4], variable.name="bac", value.name="ra")
cbbPalette <- c("#87CEFA","#FFAEB9")
p=ggplot(data0, aes(x=bac, y=ra, fill=Group)) + 
  geom_boxplot(outlier.colour = "grey", outlier.size = 0.02)+
  scale_fill_manual(values = cbbPalette) +
  theme_classic()+
  theme(plot.margin=unit(rep(0.5,4), 'lines'),
        panel.background=element_rect(fill='transparent', color='black'),
        panel.border=element_rect(fill='transparent', color='transparent'),
        panel.grid=element_blank(), legend.position = "right")+
  stat_compare_means(method = "wilcox.test",label = "p.format",position = "identity")+
  ylab("Relative abundance (1og10)")+xlab("")+
  coord_flip()
p
ggsave(paste0(prefile, "inducer_boxplot2.pdf"), p, width = 160,height = 60, units = 'mm')

##Figure 3: linear mixed model
library(lmerTest)
library(lme4)
library(qvalue)
#loading data
meta <- read.csv("metadata.csv",row.names = 1)
meta1 <- meta[meta$Group3 %in% "V3F",]
meta2 <- meta[meta$Group3 %in% "V4F",]
pheno <- meta1[,c(11,12,17)] #extract data of age, sex, and BMI
meta1 <- meta1[,c(29,31)]
meta2 <- meta2[,c(29,31)]
meta1[meta1==0] <- 1
meta1 <- log10(meta1)#log10 transformation
meta1 <- scale(meta1)#z-score transformation
meta2[meta2==0] <- 1
meta2 <- log10(meta2)
meta2 <- scale(meta2)
deltameta <- meta2-meta1
deltameta <- as.data.frame(deltameta)#delta changes of IES-R score

genus0 <- read.csv("species.csv",row.names = 1)
genus0 <- t(genus0)
genus0 <- log10(genus0+1)
genus0 <- scale(genus0)
genus0 <- as.data.frame(genus0)
genus1 <- genus0[rownames(meta1),]
genus2 <- genus0[rownames(meta2),]
deltagenus <- genus2-genus1#delta changes of the relative abundance of species

data0 <- cbind(deltameta,deltagenus,pheno)
prefile <- "D:/stress/analysis/Figure4/"
write.csv(data0, paste0(prefile,"lm_data.csv"))

#regression
result=data.frame(matrix(NA,nrow=2*462,ncol=4))
colnames(result)=c("Y","X","beta","P")
k=1
for (i in 1:2){
  for (j in 3:464){
    association=summary(lm(data0[,i]~data0[,j]+Age+Sex+BMI,data0))
    result[k,1]=colnames(data0)[i]
    result[k,2]=colnames(data0)[j]
    result[k,3]=association[["coefficients"]][2]
    result[k,4]=association[["coefficients"]][17]
    k=k+1
  }
}
result<-na.omit(result)
write.csv(result,paste0(prefile,"result.csv"))
#correction
for (i in c("PSQI","IESR")){
  resultnodifphq15<-result[result$Y==i,]
  resultnodifphq15$FDR<-p.adjust(resultnodifphq15$P,method = "BH",n=length(resultnodifphq15$P))
  x1<-qvalue(resultnodifphq15$P, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL)
  resultnodifphq15$qvalue<-x1$qvalues
  resultnodifphq15$lfdr<-x1$lfdr
  write.csv(resultnodifphq15,paste0(prefile,"result",i,".csv",sep = ""))
}
##Figure 3: Barplot
library(ggplot2)
prefile = "D:/stress/analysis/Figure4/"
data0 = read.csv(paste0(prefile, "IESR.csv"))#read the result of mixed model of IES-R
data0 = data0[data0$frequent==1, ]#data0$frequent ==1 means the species present in more than 20% of samples
data0 = na.omit(data0)
data0 = data0[order(data0$abs.beta.,decreasing = T),]
data0$X = factor(data0$X, levels = rev(data0$X))
p = ggplot(data0, aes(X, abs.beta., fill=color))+
  geom_bar(stat = "identity", width = 0.8)+
  coord_flip()+
  scale_fill_manual(values = c("blue"="#87CEFA", "pink"= "#FFAEB9"))+
  geom_text(data=data0,aes(x=X, y= abs.beta., label =text, size=8, hjust = "outward"))+
  xlab("")+
  ylab("|beta|")+theme_bw()+theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6))
p
ggsave(paste0(prefile, "lm_barplot_top.pdf"), p, width = 160, height = 80, units = 'mm')


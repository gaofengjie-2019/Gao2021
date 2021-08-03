##Figure 3: Boxplot
library(ggpubr)
library(ggsci)
meta<-read.csv("metadata.csv",row.names = 1)
data<-meta[,c(1,3,7,10,26:31)]#extract scale and group for plotting
data$Group<-factor(data$Group,levels = c("Control","Day0","Day14","Day45","Day180"))
#plot the scale in turn
k=5
for (i in colnames(data)[5:10]){
  p1=ggboxplot(data, x = "Group", y = i, color="Group", palette = "npg", add="jitter", size=0.05) +
    geom_line(aes(group = id), color = 'gray', lwd = 0.2) + #绘制箱线图
    theme_classic()+
    theme(plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),legend.position = 'none')+
    stat_compare_means(method = "wilcox.test", method.args = list(alternative = "greater"),ref.group = "Control",label = "p.format")+
    ylab(i)+xlab("")
  ggsave(paste0(prefile, i, ".pdf"), p1, width = 80*2,height = 55*2, units = 'mm')#prefile is the file that can save the picture
  k=k+1
}

##Figure 3: Random forest regression
library(randomForest)
library(dplyr)
source("rfcv1.R")

R_squa2 <- function(y, y0)
{ 
  sse <- sum((y0-y)^2)
  sst <- sum((y0-mean(y0))^2)
  R_squa <- 1-(sse/sst)
  return(R_squa) 
}
#loading data
genus <- read.csv("species.csv", row.names = 1)
meta <- read.csv("metadata.csv", row.names = 1)
genus <- t(genus) %>% as.data.frame()
meta <- meta[meta$CompleteSeries=="CS", ]
genus0 <- genus[rownames(meta), ]
genus0 <- genus0[,!names(genus0) %in% c("Unclassified")]
genus0$tp <- meta$TP
meta0 <- meta[,c(27,28,29,31)]
prefile = "D:/stress/analysis/Figure3/"#prefile is the file that can save the picture
#regression using species against the psychological scale in turn
for (i in 1:4) {
train.x <- genus0
train.y <- meta0[,i]
cv.fold <- 5
cv.step <- 0.9
cv.time <- 5
marker.num <- 0
prefix <- colnames(meta0)[i]

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
if (marker.num == 0) {
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
}
matplot(train.cv[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cv.time), main = paste("select", marker.num, "Vars"), xlab = "Number of vars", 
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

##Figure 3: Heatmap
# get union of the optimal sets of four scales through random forest
a1<-read.csv( paste0(prefile, "GAD7_marker.imp.csv"), row.names = 1)
colnames(a1)[1] <- "GAD7"
a1$x<-rownames(a1)
a1 <- a1[,-2]
a2<-read.csv(paste0(prefile, "IESR_marker.imp.csv"), row.names = 1)
colnames(a2)[1] <- "IESR"
a2$y<-rownames(a2)
a2 <- a2[,-2]
merge1 <- merge(a1, a2, all = T, by.x = "x", by.y = "y")
for (i in c("PHQ9", "PSQI")){
  a1<-read.csv(paste0(prefile, i, "_marker.imp.csv"), row.names = 1)
  colnames(a1)[1] <- i
  a1$y<-rownames(a1)
  a1 <- a1[,-2]
  a1[,1]<- scale(a1[,1])
  merge1<-merge(merge1, a1, all = T, by.x = "x", by.y = "y")
}
merge1 <- merge1[-60,]#remove timepoint
write.csv(merge1,paste0(prefile, "union.csv") )

#get annotation
tax1 = read.csv(paste0(prefile, "species_to_genus.csv"), row.names = 1)
tax2<- tax1[merge1$x,]
write.csv(tax2, paste0(prefile, "union_anno.csv"))#If some species are uncultured or unclassified at the genus level, they have to be traced back to the family level, but my algorithm only traces back to the genus level, so we will get the Na value, which can be modified manually
anno_row <- read.csv(paste0(prefile,"union_anno.csv" ))
cluster <- read.csv(paste0(prefile, "mfuzz.csv"))
cluster1 <- cluster[cluster$id %in% anno_row$genus, ]
cluster1$genus <- cluster1$id
merge2<- merge(anno_row, cluster1, all.x = T, by.x = "genus", by.y = "genus")
rownames(merge2) = merge2$X.x
merge2 = merge2[,c(1,9)]
colnames(merge2)[2] = "cluster"
write.csv(merge2, paste0(prefile,"union_anno1.csv"))

##normalization to [0,1]
scale00<- read.csv(paste0(prefile, "union.csv"),row.names = 1)
a=min(scale00, na.rm = T)
b=max(scale00, na.rm = T)
scale01 = (scale00-a)/(b-a)
scale00 = scale01
anno_row<- read.csv(paste0(prefile, "union_anno1.csv"), row.names = 1)
anno_row<-anno_row[order(anno_row$cluster, decreasing = F),]
scale00 <- scale00[rownames(anno_row), ]

##heatmap
bk0<-seq(0, 1, 0.01)
pheatmap(scale00, 
         cluster_col = F,cluster_rows = F,
         color=c(colorRampPalette(colors = c("white","#FF6A6A"))(101)),
         legend_breaks = seq(0,1,0.5),
         cellwidth = 12,
         cellheight = 8,
         breaks = bk0,
         border_color="white",
         annotation_row=anno_row)
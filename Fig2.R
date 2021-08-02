library(tidyr)
library(RColorBrewer)
library(plyr)
library(ggplot2)
genus2<-read.csv("genus.csv",row.names = 1)
genus2<-t(genus2) %>% as.data.frame()
meta<-read.csv("metadata.csv",row.names = 1)
meta0<-meta[meta$CompleteSeries %in% c("CS","Control"),]
genus2<-genus2[rownames(meta0),]
genus2$group<-meta0$Group
#compute the mean of frontline healthcare workers at each time point
genus2$group[which(genus2$group=="Day0")] <-0#替换某列的值
genus2$group[which(genus2$group=="Day14")] <-14
genus2$group[which(genus2$group=="Day45")] <-45
genus2$group[which(genus2$group=="Day180")] <-180
genus2=as.data.frame(lapply(genus2,as.numeric))
mean<-aggregate(genus2[,1:178],list(genus2$group),mean)
test<-t(mean) %>% as.data.frame()
test<-test[,c(1,2,3,5,4)]
colnames(test)<-c("Control","Day0", "Day14","Day45","Day180")
test<-test[-1,]
write.csv(test,"mean.csv")

library(Mfuzz)
library(plyr)
test<-read.csv("mean.csv",row.names = 1)
mydata<-test[,-1]
mydata<-as.matrix(mydata)
df<-new("ExpressionSet",exprs=mydata)
df<-standardise(df)
m<-mestimate(df)#m was Fuzzy parameter
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))#wss (within-cluster sum of squares)
for (i in 2:15){
  wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")#Figure S5

cl<-mfuzz(df, c=6,m=m)#c is clusters
mfuzz.plot2(df, cl=cl, mfrow=c(3,2),ylab = "Relative abundance",new.window=F, time.labels=c("Day0","Day14", "Day45" ,"Day180")) #Figure 2
library(RColorBrewer)
output<-cbind(cl$cluster,cl[["membership"]])
output<-as.data.frame(output)
mydata1<-as.data.frame(mydata)
mydata1$id<-rownames(mydata1)
output$id<-rownames(output)
joindata <- join(mydata1, output, by = "id")
write.csv(joindata, file="mfuzz.csv")#each genus with the cluster they belong to
library(glmnet)
library(caret)
library(mlbench)
library(psych)
library(dplyr)
library(plyr)
library(data.table)

asv<-read.csv("3618ASV.csv", row.names = 1)
asv<-log2(asv+1)
asvz<-scale(asv,center=T,scale=T)
meta<-read.csv("metadata.csv",row.names = 1)
asv<-t(asvz) %>% as.data.frame()
asv$id<-rownames(asv)
meta$id<-rownames(meta)
joindata<-join(asv, meta, by = "id")
joindata<-joindata[joindata$Group!="Control", ]
x<-joindata[,1:3618]
x<-as.matrix(x)
y<-joindata$IESR
y<-log2(y+1)
y<-scale(y,center=T,scale=T)
set.seed(888)
fit <- glmnet(x, y, alpha=1,family = 'gaussian')
fit_cv <- cv.glmnet(x, y, alpha=1,family="gaussian",type.measure="deviance")
plot(fit_cv)
log(fit_cv$lambda.min)
get_coe <- function(the_fit,the_lamb){
  Coefficients <- coef(the_fit, s = the_lamb)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients <- Coefficients[Active.Index]
  re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients)
  re <- data.table('var_names'=rownames(Coefficients)[Active.Index],
                   'coef'=Active.Coefficients)
  re$expcoef <- exp(re$coef)
  return(re[order(expcoef)])
}
important<-get_coe(fit_cv,fit_cv$lambda.min)
get_plot<- function(the_fit,the_fit_cv,the_lamb,toplot = seq(1,50,2)){
  Coefficients <- coef(the_fit, s = the_lamb)
  Active.Index <- which(Coefficients != 0)
  coeall <- coef(the_fit, s = the_fit_cv$lambda[toplot])
  coe <- coeall[Active.Index[-1],]
  ylims=c(-max(abs(coe)),max(abs(coe)))
  sp <- spline(log(the_fit_cv$lambda[toplot]),coe[1,],n=100)
  plot(sp,type='l',col =1,lty=1, 
       ylim = ylims,ylab = 'Coefficient', xlab = 'log(lambda)') 
  abline(h=0) 
  for(i in c(2:nrow(coe))){
    lines(spline(log(the_fit_cv$lambda[toplot]),coe[i,],n=1000),
          col =i,lty=i)
  }
  legend("bottomright",legend=rownames(coe),col=c(1:nrow(coe)),
         lty=c(1:nrow(coe)),
         cex=0.5)
}
get_plot(fit,fit_cv,fit_cv$lambda.min)
write.csv(important,"importantasv.csv")

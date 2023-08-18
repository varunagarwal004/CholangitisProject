---
output:
  pdf_document: default
  html_document: default
---
# imports
library(tidyverse)
library(gpairs)
library('vioplot')
library(ggplot2)
library('leaps')
knitr::opts_chunk$set(echo = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE)

# get the data
df <- read.csv('cholangitis.csv')
df.dict <- read.csv('cholangitis_dictionary.csv')

# remove na values
df <- na.omit(df)
head(df)
#df.dict


# taking a look at the n_days distribution
hist(df$n_days,
     main='Distribution of n_days',
     breaks=15,
     freq=FALSE)
lines(density(df$n_days), col='red')

ggplot(df) +
  geom_histogram(aes(x=n_days, group=drug, fill=drug))

source("https://www.stat.berkeley.edu/~epurdom/RcodeForClasses/myvioplot.R")
par(mfrow=c(2,2))
vioplot2(df$n_days, df$sex, xlab='sex', ylab='n_days', main='n_days by sex')
vioplot2(df$n_days, df$ascites, xlab='ascites', ylab='n_days', main='n_days by ascites')
vioplot2(df$n_days, df$hepatomegaly, xlab='hepatomegaly', ylab='n_days', main='n_days by hepatomegaly')
vioplot2(df$n_days, df$spiders, xlab='spiders', ylab='n_days', main='n_days by spiders')


df$sex<-ifelse(df$sex=="M",1,0)
df$ascites<-ifelse(df$ascites=="Y",1,0)
df$hepatomegaly<-ifelse(df$hepatomegaly=="Y",1,0)
df$spiders<-ifelse(df$spiders=="Y",1,0)
head(df)

# not run to conserve space on pdf
#for (x in 11:19) {
#  print(colnames(df)[x])
#  print(sort(boxplot(x=df[[x]], breaks = 20, xlab=colnames(df)[x], plot=FALSE)$out), sep=' ')
#  hist(x=df[[x]], breaks = 20, xlab=colnames(df)[x])
#}

print('prothrombin')
print(sort(boxplot(x=df['prothrombin'], breaks = 20, xlab='prothrombin', plot=FALSE)$out), sep=' ')
hist(x=df[['prothrombin']], breaks = 20, xlab='prothrombin')

#remove outliers as shown by the plots
df <- df[df['bilirubin'] < 15,]
df <- df[df['cholesterol'] < 1600,]
#df <- df[df['albumin'] > 2.56 & df['albumin'] < 4.52,]
df <- df[df['copper'] < 550,]
df <- df[df['alk_phos'] < 13000,]
df <- df[df['sgot'] < 400,]
df <- df[df['tryglicerides'] < 500,]
#df <- df[df['platelets'] < 514,]
df <- df[df['prothrombin'] < 16,]
head(df)

panel.hist <- function(x, ...)
{
	#from help of pairs
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
	#from help of pairs
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
suppressWarnings(
  pairs(df[,c(2,5:9, 11:20)], lower.panel = panel.smooth, upper.panel = panel.cor,col=c("red","black")[df$drug],diag.panel = panel.hist)
)

suppressWarnings(corrgram(df[,-c(1,3,4, 10)]))

coplot(n_days ~ stage | albumin, data=df)

# not run to conserve pdf space
#for (i in 5:20) {
#  plt <- ggplot(df)+
#  geom_point(aes(y=df$n_days, x=df[[i]])) +
#  geom_smooth(aes(y=df$n_days, x=df[[i]])) +
#  labs(y='n_days', x=colnames(df)[i])
#  print(plt)
#}


plt <- ggplot(df)+
  geom_point(aes(y=df$n_days, x=log(alk_phos))) +
  geom_smooth(aes(y=df$n_days, x=log(alk_phos))) +
  labs(y='n_days', x='alk_phos')
print(plt)

# not run to conserve pdf space
#for (x in 11:19) {
#  par(mfrow=c(2,2))
#  hist(x=df[[x]], breaks = 20, xlab=colnames(df)[x])
#  hist(x=log(df[[x]]), breaks = 20, xlab=colnames(df)[x])
#  hist(x=df[[x]]^2, breaks = 20, xlab=colnames(df)[x])
#  hist(x=df[[x]]^(1/2), breaks = 20, xlab=colnames(df)[x])
#}

par(mfrow=c(2,2))
hist(x=df[['bilirubin']], breaks = 20, xlab='bilirubin')
hist(x=log(df[['bilirubin']]), breaks = 20, xlab='bilirubin')
hist(x=df[['bilirubin']]^2, breaks = 20, xlab='bilirubin')
hist(x=df[['bilirubin']]^(1/2), breaks = 20, xlab='bilirubin')

model1<-lm(n_days ~ . - id, data=df)
print(summary(model1))
plot(x=df$n_days, y=model1$fitted.values, ylab='fitted values',xlab='true values', main='regression fit')
abline(b = 1, a = 0, col='red')   

step(model1, trace=0, direction='both')

regs <- regsubsets(n_days ~ . - id, df)
summary(regs)

model2<-lm(n_days ~ log(bilirubin) + log(cholesterol) + albumin + log(copper) + log(alk_phos) + log(sgot) + log(tryglicerides) + platelets + log(prothrombin) + as.factor(status) + drug + age + ascites + hepatomegaly + spiders + edema + stage + sex, data=df)
summary(model2)

step(model2, trace=0, direction='both')

regs <- regsubsets(n_days ~ log(bilirubin) + log(cholesterol) + albumin + log(copper) + log(alk_phos) + log(sgot) + log(tryglicerides) + platelets + log(prothrombin) + as.factor(status) + drug + age + ascites + hepatomegaly + spiders + edema + stage + sex, df)
summary(regs)

modelFinal <- lm(n_days ~ log(bilirubin) + albumin + log(copper) +  log(alk_phos) + log(prothrombin) + as.factor(status) + as.factor(edema) + as.factor(stage), df)
print(summary(modelFinal))
par(mfrow=c(2,2))
plot(modelFinal)

df2 <- df[,c('n_days', "bilirubin", "albumin", "copper", "alk_phos", "prothrombin", "stage", "status", "edema")]
df2['bilirubin'] <- log(df2$bilirubin)
df2['copper'] <- log(df2$copper)
df2['alk_phos'] <- log(df2$alk_phos)
df2['prothrombin'] <- log(df2$prothrombin)
head(df2)

pairs(df2[,-c(8,9)])


set.seed(1249)
nTest <- 0.1 * nrow(df)
whTest <- sample(1:nrow(df), size = nTest)
dfTest <- df[whTest, ]
dfTrain <- df[-whTest, ]
df_subsets <- regsubsets(n_days ~ . - id, df)

modelTrain <- lm(n_days ~ log(bilirubin) + albumin + log(copper) +  log(alk_phos) + log(prothrombin) + as.factor(status) + as.factor(edema) + as.factor(stage), dfTrain)
preds <- predict(modelTrain, newdata = dfTest[,-2])
mean((dfTest$n_days - preds)^2)

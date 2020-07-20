library(readxl)
library(ggplot2)
library(nlme)
library(plyr)
library(dplyr)
library(stargazer)
library(metap)
library(metafor)
library(texreg)
library(lsmeans)
library(doBy)
library(reshape)
library(pastecs)


# Generate the table that returns the 2.5%, 50%, and 97.5%
# quantile(array, c(.025, .50, .975))

funSummaryBy <- function (x){
# c(LOW=quantile(x, c(0.025)), MEAN=mean(x), median=median(x), UP=quantile(x, c(0.975)))
  c(MEAN=mean(x), median=median(x), MIN=min(x,na.rm =TRUE), MAX=max(x,na.rm =TRUE))
}

#################################################################################################################
########################################### Increase precision of results  ######################################
#################################################################################################################


d <- as.data.frame(read_excel("Raw Data Testing.xlsx"))
d<-subset(subset(d, (TECHNIQUE=="BT" | TECHNIQUE=="EP")))
#d<-subset(subset(d, EFF<=1 & (TECHNIQUE=="BT" | TECHNIQUE=="EP")))

d$EXPERIMENT=as.factor(d$EXPERIMENT)
d$SUBJECT_ID=as.factor(d$SUBJECT_ID)
d$PROGRAM=as.factor(d$PROGRAM)
d$TECHNIQUE=as.factor(d$TECHNIQUE)
d$VERSION=as.factor(d$VERSION)

dUPM01 <-subset(d, EXPERIMENT=="UPM01")
dUPM02 <-subset(d, EXPERIMENT=="UPM02")
dUPM03 <-subset(d, EXPERIMENT=="UPM03")
dUPM04 <-subset(d, EXPERIMENT=="UPM04")
dUPM05 <-subset(d, EXPERIMENT=="UPM05")

# Run the t-tests
modelupm01 <- t.test(EFF ~ TECHNIQUE, data=dUPM01, paired = TRUE)
modelupm02 <- t.test(EFF ~ TECHNIQUE, data=dUPM02, paired = TRUE)
modelupm03 <- t.test(EFF ~ TECHNIQUE, data=dUPM03, paired = TRUE)
modelupm04 <- t.test(EFF ~ TECHNIQUE, data=dUPM04, paired = TRUE)
modelupm05 <- t.test(EFF ~ TECHNIQUE, data=dUPM05, paired = TRUE, na.action=na.omit)
modelupm01
modelupm02
modelupm03
modelupm04
modelupm05


d<- rbind(dUPM01,dUPM02,dUPM03,dUPM04,dUPM05)

experiences <- distinct(d, SUBJECT_ID, .keep_all = TRUE)

mdf<-do.call(data.frame, aggregate(EFF ~ EXPERIMENT, FUN = function(x) c(N = length(x)), data=experiences))
print(mdf)

mdf<-do.call(data.frame, aggregate(EFF ~ TECHNIQUE + EXPERIMENT, FUN = function(x) c(N = length(x), MEAN = mean(x), SD=sd(x) ,MDN=median(x)), data=d))
print(mdf)
stargazer(mdf[,c("EXPERIMENT", "TECHNIQUE", "EFF.N", "EFF.MEAN", "EFF.SD","EFF.MDN")], summary=FALSE, digits = 3, rownames = FALSE)


standardizedDf <- data.frame(EXPERIMENT=character(0), ESTIMATE=numeric(0), VARIANCE=numeric(0))
by(d,d$EXPERIMENT, function(dTemp){
  
  newd <-  dTemp %>% group_by(SUBJECT_ID) %>% filter(n()>1) #
  
  ITL <- newd[which(newd$TECHNIQUE=="BT"),]$EFF
  TDD <- newd[which(newd$TECHNIQUE=="EP"),]$EFF
  
  n1 <- length(ITL)
  n2 <- length(TDD)
  
  numerator=((n1-1)*(sd(ITL,na.rm=TRUE)^2))+((n2-1)*sd(TDD,na.rm=TRUE)^2)
  denominator=n1+n2-2
  sPooled=sqrt(numerator/denominator)
  cohenD <- (mean(TDD,na.rm=TRUE) - mean(ITL,na.rm=TRUE))/sPooled
  
  vd <-0
  n <- length(c(ITL))
  corITLTDD <-cor(TDD,ITL, use="pairwise.complete.obs")
  vd <- 2*((1-corITLTDD)/n)+ ((cohenD^2)/(2*n-2))
  df<- n-1
  J= 1- (3/(4*df-1))
  g=J*cohenD
  vg=(J^2)*vd
  
  rowDf <- data.frame(EXPERIMENT=dTemp$EXPERIMENT[1], yi=g, vi=vg)
  standardizedDf <<- rbind(standardizedDf, rowDf)
})

summary(res1<-rma(yi=yi,vi=vi,slab=standardizedDf$EXPERIMENT,method="REML", data=standardizedDf))

par(mfrow=c(1,1)) 
par(font=1)
forest(res1, xlab="Effect Size",mlab="")
par(font=2)
### add column headings to the plot
text(c(-0.3,0.3), 7, c("EP<BT", "EP>BT"))
text(-2.63, 7 , "Experiment",  pos=4)
text(2.4, 7, "Effect Size [95% CI]", pos=2)
text(-2.63, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
                                             .(formatC(res1$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res1$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res1$I2, digits=1, format="f")), "%)")))

confint(res1)

#################################################################################################################
########################################### Extract participant-level moderators ################################
#################################################################################################################

d <- as.data.frame(read_excel("Raw Data TDD.xlsx"))
d<- d[,c("SUBJECT_ID","EXPERIMENT","TREATMENT","PROD","SESSION","PROGRAMMING_LANGUAGE","SUBJECT_TYPE","EXPERIENCE_PROGRAMMING","EXPERIENCE_JAVA","EXPERIENCE_UNIT_TESTING","EXPERIENCE_JUNIT")]
d$TREATMENT[d$TREATMENT==1] <- "ITL"
d$TREATMENT[d$TREATMENT==2] <- "TDD"
d$TREATMENT<- as.factor(d$TREATMENT)
d$SUBJECT_ID<-as.factor(d$SUBJECT_ID)
d$EXPERIMENT=as.factor(d$EXPERIMENT)
d$SUBJECT_TYPE=as.factor(d$SUBJECT_TYPE)
d$PROGRAMMING_LANGUAGE=as.factor(d$PROGRAMMING_LANGUAGE)


# Select F-Secure and UPV
dFsecure_O <- subset(d, EXPERIMENT=="F-Secure O" & (SESSION==1 | SESSION==2))
dFsecure_H <- subset(d, EXPERIMENT=="F-Secure H" & (SESSION==1 | SESSION==2))
dFsecure_K <- subset(d, EXPERIMENT=="F-Secure K" & (SESSION==1 | SESSION==2))

# Run the t-tests
modelFsecure_O <- t.test(PROD ~ TREATMENT, data=dFsecure_O, paired = TRUE)
modelFsecure_H <- t.test(PROD ~ TREATMENT, data=dFsecure_H, paired = TRUE)
modelFsecure_K <- t.test(PROD ~ TREATMENT, data=dFsecure_K, paired = TRUE)
modelFsecure_O
modelFsecure_H
modelFsecure_K

d<- rbind(dFsecure_O,dFsecure_H, dFsecure_K)

d <- d[with(d, order(EXPERIMENT, TREATMENT, SUBJECT_ID)), ]


experiences <- distinct(d, SUBJECT_ID, .keep_all = TRUE)

by(experiences, experiences$EXPERIMENT, function(dTemp){
  nrow(dTemp)
})


experiences <- experiences[,c("EXPERIMENT","EXPERIENCE_PROGRAMMING","EXPERIENCE_JAVA","EXPERIENCE_UNIT_TESTING","EXPERIENCE_JUNIT")]
allExperiences <- melt(experiences,id.vars = "EXPERIMENT",variable_name = "EXPERIENCE")
summaryExperiences<-summaryBy(value ~ EXPERIENCE + EXPERIMENT, FUN=funSummaryBy, data=allExperiences)
stargazer(summaryExperiences, summary=FALSE, rownames = FALSE, digits = 2)


mdf<-do.call(data.frame, aggregate(PROD ~ TREATMENT + EXPERIMENT, FUN = function(x) c(N = length(x), MEAN = mean(x), SD=sd(x) ,MDN=median(x)), data=d))
stargazer(mdf[,c("EXPERIMENT","TREATMENT","PROD.N","PROD.MEAN","PROD.SD","PROD.MDN")], summary=FALSE, rownames = FALSE, digits = 2)

mdf <- summaryBy(. ~ EXPERIMENT, data=experiences, na.rm=TRUE)
longMdf<-melt(mdf, id.vars = c("EXPERIMENT"))
longMdf$variable <- as.character(longMdf$variable)
longMdf$variable[longMdf$variable=="EXPERIENCE_PROGRAMMING.mean"] <- "Programming"
longMdf$variable[longMdf$variable=="EXPERIENCE_JAVA.mean"] <- "Java"
longMdf$variable[longMdf$variable=="EXPERIENCE_UNIT_TESTING.mean"] <- "Unit Test"
longMdf$variable[longMdf$variable=="EXPERIENCE_JUNIT.mean"] <- "JUnit"

longMdf$variable<- as.factor(longMdf$variable)
longMdf$variable<-factor(longMdf$variable,levels(longMdf$variable)[c(3,1,4,2)])


#### Profile plot for experiences
ggplot(data=longMdf, aes(x=as.factor(variable), y=value, group = EXPERIMENT, colour = EXPERIMENT)) +
  geom_line() + geom_point( size=4, shape=21, fill="white")+ labs(color='Experiment') +
  labs(x="", y="Mean experience") + theme_bw() + scale_y_continuous(limits=c(1,4))




################################## AD
# Main Analysis
standardizedDf <- data.frame(EXPERIMENT=character(0), ESTIMATE=numeric(0), VARIANCE=numeric(0))
by(d,d$EXPERIMENT, function(dTemp){
  
  newd <-  dTemp %>% group_by(SUBJECT_ID) %>% filter(n()>1) #
  
  ITL <- newd[which(newd$TREATMENT=="ITL"),]$PROD
  TDD <- newd[which(newd$TREATMENT=="TDD"),]$PROD
  
  n1 <- length(ITL)
  n2 <- length(TDD)
  
  numerator=((n1-1)*(sd(ITL,na.rm=TRUE)^2))+((n2-1)*sd(TDD,na.rm=TRUE)^2)
  denominator=n1+n2-2
  sPooled=sqrt(numerator/denominator)
  cohenD <- (mean(TDD,na.rm=TRUE) - mean(ITL,na.rm=TRUE))/sPooled
  
  vd <-0
  n <- length(c(ITL))
  corITLTDD <-cor(TDD,ITL, use="pairwise.complete.obs")
  vd <- 2*((1-corITLTDD)/n)+ ((cohenD^2)/(2*n-2))
  df<- n-1
  J= 1- (3/(4*df-1))
  g=J*cohenD
  vg=(J^2)*vd
  
  rowDf <- data.frame(EXPERIMENT=dTemp$EXPERIMENT[1], yi=g, vi=vg)
  standardizedDf <<- rbind(standardizedDf, rowDf)
})

summary(res1<-rma(yi=yi,vi=vi,slab=standardizedDf$EXPERIMENT,method="REML", data=standardizedDf))

par(mfrow=c(1,1)) 
par(font=1)
forest(res1, xlab="Effect Size",mlab="")
par(font=2)
### add column headings to the plot
text(c(-0.7,0.7), 5, c("TDD<ITL", "TDD>ITL"))
text(-4.4, 5 , "Experiment",  pos=4)
text(6.9, 5, "Effect Size [95% CI]", pos=2)
text(-4.4, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
                                             .(formatC(res1$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res1$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res1$I2, digits=1, format="f")), "%)")))



confint(res1)


daux <- d

daux$TREATMENT=as.character(daux$TREATMENT)
daux$TREATMENT[daux$TREATMENT=="ITL"] <- 0
daux$TREATMENT[daux$TREATMENT=="TDD"] <- 1
daux$TREATMENT<-as.numeric(daux$TREATMENT)



# Participant-level variables
experienceProgrammingMean <- aggregate(EXPERIENCE_PROGRAMMING   ~ EXPERIMENT, FUN = mean, data=experiences)
colnames(experienceProgrammingMean)[2] <- "MEAN_EXPERIENCE_PROGRAMMING"
experienceJavaMean <- aggregate(EXPERIENCE_JAVA   ~ EXPERIMENT, FUN = mean, data=experiences)
colnames(experienceJavaMean)[2] <- "MEAN_EXPERIENCE_JAVA"
experienceUnitTestingMean <- aggregate(EXPERIENCE_UNIT_TESTING   ~ EXPERIMENT, FUN = mean, data=experiences)
colnames(experienceUnitTestingMean)[2] <- "MEAN_EXPERIENCE_UNIT_TESTING"
experienceJunitMean <- aggregate(EXPERIENCE_JUNIT   ~ EXPERIMENT, FUN = mean, data=experiences)
colnames(experienceJunitMean)[2] <- "MEAN_EXPERIENCE_JUNIT"

daux <- merge(daux,experienceProgrammingMean,by="EXPERIMENT")
daux <- merge(daux,experienceJavaMean,by="EXPERIMENT")
daux <- merge(daux,experienceUnitTestingMean,by="EXPERIMENT")
daux <- merge(daux,experienceJunitMean,by="EXPERIMENT")

daux$WITHIN_EXPERIENCE_PROGRAMMING <- daux$EXPERIENCE_PROGRAMMING-daux$MEAN_EXPERIENCE_PROGRAMMING
daux$WITHIN_EXPERIENCE_JAVA <- daux$EXPERIENCE_JAVA-daux$MEAN_EXPERIENCE_JAVA
daux$WITHIN_EXPERIENCE_UNIT_TESTING <- daux$EXPERIENCE_UNIT_TESTING-daux$MEAN_EXPERIENCE_UNIT_TESTING
daux$WITHIN_EXPERIENCE_JUNIT <- daux$EXPERIENCE_JUNIT-daux$MEAN_EXPERIENCE_JUNIT


summary(lmeResult1<- lme(PROD ~ TREATMENT + EXPERIENCE_PROGRAMMING + TREATMENT:MEAN_EXPERIENCE_PROGRAMMING + TREATMENT:WITHIN_EXPERIENCE_PROGRAMMING, random=list(EXPERIMENT=~TREATMENT, SUBJECT_ID=~1), data=daux, na.action = na.exclude, control=lmeControl(opt='optim')))
summary(lmeResult2<- lme(PROD ~ TREATMENT + EXPERIENCE_JAVA + TREATMENT:MEAN_EXPERIENCE_JAVA + TREATMENT:WITHIN_EXPERIENCE_JAVA, random=list(EXPERIMENT=~TREATMENT, SUBJECT_ID=~1), data=daux, na.action = na.exclude, control=lmeControl(opt='optim')))
summary(lmeResult3<- lme(PROD ~ TREATMENT + EXPERIENCE_UNIT_TESTING + TREATMENT:MEAN_EXPERIENCE_UNIT_TESTING + TREATMENT:WITHIN_EXPERIENCE_UNIT_TESTING, random=list(EXPERIMENT=~TREATMENT, SUBJECT_ID=~1), data=daux, na.action = na.exclude, control=lmeControl(opt='optim')))
summary(lmeResult4<- lme(PROD ~ TREATMENT + EXPERIENCE_JUNIT + TREATMENT:MEAN_EXPERIENCE_JUNIT + TREATMENT:WITHIN_EXPERIENCE_JUNIT, random=list(EXPERIMENT=~TREATMENT, SUBJECT_ID=~1), data=daux, na.action = na.exclude, control=lmeControl(opt='optim')))

intervals(lmeResult1)
intervals(lmeResult2)
intervals(lmeResult3)
intervals(lmeResult4)


a1<- summary(lmeResult1)
a2<- summary(lmeResult2)
a3<- summary(lmeResult3)
a4<- summary(lmeResult4)


int1<-a1$tTable[5,1]
int2<-a2$tTable[5,1]
int3<-a3$tTable[5,1]
int4<-a4$tTable[5,1]

int1<-c(1,2,3,4)*int1
int2<-c(1,2,3,4)*int2
int3<-c(1,2,3,4)*int3
int4<-c(1,2,3,4)*int4

int1df<-data.frame(EXPERIENCE=rep("Programming",4), X=seq(1,4), VALUE=int1 )
int2df<-data.frame(EXPERIENCE=rep("Java",4), X=seq(1,4), VALUE=int2 )
int3df<-data.frame(EXPERIENCE=rep("Unit Testing",4), X=seq(1,4), VALUE=int3 )
int4df<-data.frame(EXPERIENCE=rep("JUnit",4), X=seq(1,4), VALUE=int4 )

df <- rbind(int1df, int2df, int3df, int4df)


ggplot(data=df, aes(x=X, y=VALUE, group = EXPERIENCE, colour = EXPERIENCE)) +
  geom_line() +   geom_point( size=4, shape=21, fill="white")+ labs(color='Experience') +
  labs(x="Experience", y="QLTY") + scale_x_continuous(limits=c(1, 4)) + theme_bw() + scale_y_continuous(limits=c(0,100))


#################################################################################################################
########################################### Extract experiment-level moderators #################################
#################################################################################################################
# UPM, Ericsson and Paf

d <- as.data.frame(read_excel("Raw Data All Experiments TDD.xlsx"))
d<- d[,c("SUBJECT_ID","EXPERIMENT","TREATMENT","PROD","SESSION","PROGRAMMING_LANGUAGE","SUBJECT_TYPE","EXPERIENCE_PROGRAMMING","EXPERIENCE_JAVA","EXPERIENCE_UNIT_TESTING","EXPERIENCE_JUNIT")]
d$TREATMENT[d$TREATMENT==1] <- "ITL"
d$TREATMENT[d$TREATMENT==2] <- "TDD"
d$TREATMENT<- as.factor(d$TREATMENT)
d$SUBJECT_ID<-as.factor(d$SUBJECT_ID)
d$EXPERIMENT=as.factor(d$EXPERIMENT)
d$SUBJECT_TYPE=as.factor(d$SUBJECT_TYPE)
d$PROGRAMMING_LANGUAGE=as.factor(d$PROGRAMMING_LANGUAGE)


# Select data
dUPM <- subset(d, EXPERIMENT=="UPM" & (SESSION==1 | SESSION==2))
dEricsson <- subset(d, EXPERIMENT=="Ericsson" & (SESSION==1 | SESSION==2))
dPlayTech <- subset(d, EXPERIMENT=="PlayTech" & (SESSION==1 | SESSION==2))

d<- rbind(dUPM,dEricsson, dPlayTech)
d$EXPERIMENT <- droplevels(d$EXPERIMENT)

# t-tests
modelUPM <- t.test(PROD ~ TREATMENT, data=dUPM, paired = TRUE)
modelEricsson <- t.test(PROD ~ TREATMENT, data=dEricsson, paired = TRUE)
modelPlayTech <- t.test(PROD ~ TREATMENT, data=dPlayTech, paired = TRUE)
modelUPM
modelEricsson
modelPlayTech

################################## AD
# Main Analysis
standardizedDf <- data.frame(EXPERIMENT=character(0), ESTIMATE=numeric(0), VARIANCE=numeric(0))
by(d,d$EXPERIMENT, function(dTemp){
  
  newd <-  dTemp %>% group_by(SUBJECT_ID) %>% filter(n()>1) #
  
  ITL <- newd[which(newd$TREATMENT=="ITL"),]$PROD
  TDD <- newd[which(newd$TREATMENT=="TDD"),]$PROD
  
  n1 <- length(ITL)
  n2 <- length(TDD)
  
  numerator=((n1-1)*(sd(ITL,na.rm=TRUE)^2))+((n2-1)*sd(TDD,na.rm=TRUE)^2)
  denominator=n1+n2-2
  sPooled=sqrt(numerator/denominator)
  cohenD <- (mean(TDD,na.rm=TRUE) - mean(ITL,na.rm=TRUE))/sPooled
  
  vd <-0
  n <- length(c(ITL))
  corITLTDD <-cor(TDD,ITL, use="pairwise.complete.obs")
  vd <- 2*((1-corITLTDD)/n)+ ((cohenD^2)/(2*n-2))
  df<- n-1
  J= 1- (3/(4*df-1))
  g=J*cohenD
  vg=(J^2)*vd
  
  rowDf <- data.frame(EXPERIMENT=dTemp$EXPERIMENT[1], yi=g, vi=vg)
  standardizedDf <<- rbind(standardizedDf, rowDf)
})
standardizedDf$EXPERIMENT<-droplevels(standardizedDf)$EXPERIMENT
summary(res<-rma(yi=yi,vi=vi,slab=standardizedDf$EXPERIMENT,method="REML", data=standardizedDf))


par(mfrow=c(1,1)) 
par(font=1)
forest(res, xlab="Effect Size",mlab="")
par(font=2)
### add column headings to the plot
text(c(-0.3,0.3), 5, c("TDD<ITL", "TDD>ITL"))
text(-3.34, 5 , "Experiment",  pos=4)
text(2.7, 5, "Effect Size [95% CI]", pos=2)
text(-3.34, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
                                              .(formatC(res$QE, digits=2, format="f")),
                                              ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(res$I2, digits=1, format="f")), "%)")))


confint(res)


#Subgroup analysis technological environment

allEstimatesSE1 <- standardizedDf[which(standardizedDf$EXPERIMENT!="Ericsson"),]
allEstimatesSE2 <- standardizedDf[which(standardizedDf$EXPERIMENT=="Ericsson"),]
summary(res1<-rma(yi=allEstimatesSE1$yi,vi=allEstimatesSE1$vi,slab=allEstimatesSE1$EXPERIMENT,method="REML"))
summary(res2<-rma(yi=allEstimatesSE2$yi,vi=allEstimatesSE2$vi,slab=allEstimatesSE2$EXPERIMENT,method="REML"))


par(mar=c(4,4,1,2))
par(font=1)
forest(res, xlab="Effect Size",mlab="", ylim=c(-1, 10), rows=c(6:5,2), order=c(3,2,1))
par(font=2)
text(-2.85, 7, "Java and JUnit")
par(font=1)
addpoly(res1, row=4,  mlab="")
par(font=2)
text(-2.85, 3, "C++ and Boost")
par(font=1)
addpoly(res2, row=1,  mlab="")


par(font=2)
### add column headings to the plot
text(c(-0.5,0.5), 9, c("TDD<ITL", "TDD>ITL"))
text(-3.35, 9 , "Experiment",  pos=4)
text(2.7, 9, "Effect Size [95% CI]", pos=2)
text(-3.35, -1, pos=4, cex=0.75, bquote(paste("Joint result (Q = ",
                                              .(formatC(res$QE, digits=2, format="f")),
                                              ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(res$I2, digits=1, format="f")), "%)")))

par(font=1)
#text(-3.35, 4, pos=4, cex=0.75, bquote(paste("Subgroup joint result")))
#text(-3.35, 1, pos=4, cex=0.75, bquote(paste("Subgroup joint result")))

text(-3.35, 4, pos=4, cex=0.75, bquote(paste("RE Model subgroup (Q = ",
                                             .(formatC(res1$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res1$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res1$I2, digits=1, format="f")), "%)")))

text(-3.35, 1, pos=4, cex=0.75, bquote(paste("RE Model subgroup (Q = ",
                                             .(formatC(res2$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res2$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res2$I2, digits=1, format="f")), "%)")))


confint(res1)

#Subgroup analysis trainer experience

allEstimatesSE3 <- standardizedDf[which(standardizedDf$EXPERIMENT!="PlayTech"),]
allEstimatesSE4 <- standardizedDf[which(standardizedDf$EXPERIMENT=="PlayTech"),]
summary(res3<-rma(yi=allEstimatesSE3$yi,vi=allEstimatesSE3$vi,slab=allEstimatesSE3$EXPERIMENT,method="REML"))
summary(res4<-rma(yi=allEstimatesSE4$yi,vi=allEstimatesSE4$vi,slab=allEstimatesSE4$EXPERIMENT,method="REML"))


par(mar=c(4,4,1,2))
par(font=1)
forest(res, xlab="Effect Size",mlab="", ylim=c(-1, 10), rows=c(6:5,2), order=c(1,3,2))
par(font=2)
text(-2.85, 7, "High")
par(font=1)
addpoly(res3, row=4,  mlab="")
par(font=2)
text(-2.85, 3, "Low")
par(font=1)
addpoly(res4, row=1,  mlab="")


par(font=2)
### add column headings to the plot
text(c(-0.5,0.5), 9, c("TDD<ITL", "TDD>ITL"))
text(-3.35, 9 , "Experiment",  pos=4)
text(2.7, 9, "Effect Size [95% CI]", pos=2)
text(-3.35, -1, pos=4, cex=0.75, bquote(paste("Joint result (Q = ",
                                              .(formatC(res$QE, digits=2, format="f")),
                                              ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(res$I2, digits=1, format="f")), "%)")))

par(font=1)
#text(-3.35, 4, pos=4, cex=0.75, bquote(paste("Subgroup joint result")))
#text(-3.35, 1, pos=4, cex=0.75, bquote(paste("Subgroup joint result")))

text(-3.35, 4, pos=4, cex=0.75, bquote(paste("RE Model subgroup (Q = ",
                                             .(formatC(res3$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res3$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res3$I2, digits=1, format="f")), "%)")))

text(-3.35, 1, pos=4, cex=0.75, bquote(paste("RE Model subgroup (Q = ",
                                             .(formatC(res4$QE, digits=2, format="f")),
                                             ", p = ", .(formatC(res4$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                             .(formatC(res4$I2, digits=1, format="f")), "%)")))

confint(res3)


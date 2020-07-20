## @knitr all
library(texreg)
library(ggplot2)
library(metafor)
library(doBy)
library(orddom)
library(gridExtra)
library(cowplot)
library(stargazer)
library(effsize)
set.seed(123)

#################################################################################
############################ Functions Section ##################################
#################################################################################


listPlots <- list()

## 1. Function that plots a series of plots in the same row
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## 2. This function returns the simulated data for each experiment
# generateRandomData(typeExperiment="between", experimentID=1, numSubjectsPerArm=20, unbalance=0.5, corr=0.2, meanTDD=38, sdTDD=10, meanITL=30, sdITL=10)
generateRandomData <- function(typeExperiment, experimentID, numSubjectsPerArm, unbalance, corr, meanTDD, sdTDD, meanITL, sdITL){
  datITL <-rnorm(numSubjectsPerArm, mean=meanITL, sd=sdITL)
  datTDD <-rnorm(numSubjectsPerArm, mean=meanTDD, sd=sdTDD)
  # if balance > 0 #ITL>#TDD
  if (typeExperiment=="between"){
    if (unbalance>0){ # if unbalance>0 then subjectsITL>subjectsTDD
      numSubjectsITL= numSubjectsPerArm
      numSubjectsTDD= round(numSubjectsPerArm*(1-unbalance))
    }else if (unbalance<0){ # if balance < 0 subjectsTDD>subjectsITL
      numSubjectsTDD= numSubjectsPerArm
      numSubjectsITL= round(numSubjectsPerArm*(1+unbalance))
    }else{
      numSubjectsITL= numSubjectsPerArm
      numSubjectsTDD= numSubjectsPerArm
    }
  }
  else{
    numSubjectsITL= numSubjectsPerArm
    numSubjectsTDD= numSubjectsPerArm
  }
  
  ITL <- datITL[1:numSubjectsITL]
  TDD <- datTDD[1:numSubjectsTDD]
  
  
  PROD <- c(ITL,TDD)
  if (typeExperiment=="between"){
    SUBJECT_ID <-c(c(1:numSubjectsITL),c((numSubjectsITL+1):(numSubjectsITL+numSubjectsTDD)))
  }else if (typeExperiment=="within"){
    SUBJECT_ID <-c(c(1:numSubjectsITL),c(1:numSubjectsTDD))
  }
  TREATMENT <- c(rep("ITL",numSubjectsITL), rep("TDD",numSubjectsTDD))
  SUBJECTID <- as.factor(paste(experimentID,"_", SUBJECT_ID, sep=""))
  d <- data.frame(PROD=PROD, SUBJECT_ID=as.factor(SUBJECTID), TREATMENT=as.factor(TREATMENT), TREATMENT_NUM=TREATMENT, EXPERIMENT=as.factor(experimentID))
}


## This function returns (Cohen's d, var, lower, upper) following Borenstein et al.
cohenD <- function(data){
  
  ITL<-subset(data, TREATMENT=="ITL")$PROD
  TDD<-subset(data, TREATMENT=="TDD")$PROD
  
  val<-cohen.d(TDD, ITL)
  
  cohend<-val$estimate
  vd<-val$var
  low<-val$conf.int[1]
  up<-val$conf.int[2]
  
  return(c(cohend,vd,low,up))
}

# Generate the table that returns the 2.5%, 50%, and 97.5%
# quantile(array, c(.025, .50, .975))

funSummaryBy <- function (x){
  c(LOW=quantile(x, c(0.025)), MEAN=mean(x), median=median(x), UP=quantile(x, c(0.975)))
}

#######################################################################################################################################
############################################# 1. Simulate single experiment and run t-test ############################################
#######################################################################################################################################


cols <- c("LINE1"="red","LINE2"="blue")
df <- data.frame(x=seq(0, 100, 0.01))
ggplot(data=df, aes(x=x)) + 
  stat_function(aes(colour = "LINE1"), fun=dnorm, args=list(mean=50, sd=10)) +
  stat_function(aes(colour = "LINE2"), fun=dnorm, args=list(mean=58, sd=10)) +
  theme_bw(base_size = 14) +
  theme(legend.position ="top") +
  xlab("Performance") +
  ylab("Density") +
  scale_colour_manual("", values=c(LINE1="red",LINE2="blue"), labels = c("ITL","TDD")) + scale_x_continuous(breaks = c(0, 25, 50, 58, 75, 100)) + geom_vline(xintercept=50, linetype="dashed", color = "red") + geom_vline(xintercept = 58, linetype="dotted", color = "blue")

d <- generateRandomData("between", 1, 20, 0, 0, meanTDD=58, sdTDD=10, meanITL=50, sdITL=10)

pD1 <- ggplot(d, aes(x=TREATMENT, y=PROD))
pD2 <- pD1 + geom_violin(aes(fill = factor(TREATMENT)), width=0.9) + geom_boxplot(width=0.1, outlier.shape=NA, fill="white") + geom_jitter(width = 0.05) + stat_summary(fun.y=mean, geom="point", colour="black", shape=3, size=4, show_guide = FALSE)
pD3 <- pD2 # + facet_grid(. ~ EXPERIMENT) # + geom_text(data = means, aes(label = PROD, y = PROD + 4)) 
pD4 <- pD3  # + ggtitle("Productivity by treatment and experiment")
pD5 <- pD4 + theme_bw() + theme(legend.position="none")  +xlab("Treatment") +ylab("Performance")
print(pD5)


mdf<-do.call(data.frame, aggregate(PROD ~ TREATMENT, FUN = function(x) c(N = length(x), MEAN = mean(x), SD=sd(x) ,MDN=median(x)), data=d))
stargazer(mdf[,c("TREATMENT","PROD.N","PROD.MEAN","PROD.SD","PROD.MDN")], summary=FALSE, rownames = FALSE, digits = 2)

# Run the t-test
model1 <- lm(PROD ~ TREATMENT, data=d)
summary(model1)
confint(model1)

model2 <- t.test(PROD ~ TREATMENT, data=d, PAIRED=FALSE)
model2


cohenD(d)

######################################################################################################################
############################################## 2. Simulations with p-Values ##########################################
######################################################################################################################

print("P-values")

# effSizes <- seq(from = 1, to = 10, by =1)
effSizes <- c(2,5,8)
# numSubjectsPerArm <- seq(from = 10, to = 100, by = 10)
# Following Dyba
numSubjectsPerArm <- seq(from = 2, to = 74, by = 8)
maxSimulations=5000

resultsDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), P=numeric(0), SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0))
for(effSize in effSizes){
  for(numSubjects in numSubjectsPerArm){
    print(numSubjects)
      for(i in 1:maxSimulations){
        d <- generateRandomData("between", i, numSubjects, 0, 0, meanTDD=(50+effSize), sdTDD=10, meanITL=50, sdITL=10)
        # browser()
        res<-cohenD(d)
        estimate<-res[1]
        pValue<-2
        # if lower and upper bound have the same sign and result is positive then statistically significant
        if ((sign(res[3])==sign(res[4]))){
          pValue<-0  
        }else{
          pValue<-1
        }
        rowTemp <- data.frame(EXPERIMENT_ID=paste(effSize,"_",numSubjects,"_",i, sep=""), ESTIMATE=estimate, P=pValue, SAMPLE_SIZE=numSubjects*2, EFFECT_SIZE=effSize)
        resultsDf <- rbind(resultsDf, rowTemp)
      }
  }
  write.csv(resultsDf, paste("intermediatea",effSize,".csv", sep=""))
  resultsDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), P=numeric(0), SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0))
}

resultsDf2 <- read.csv("intermediatea2.csv")
resultsDf5 <- read.csv("intermediatea5.csv")
resultsDf8 <- read.csv("intermediatea8.csv")

resultsDf <- rbind(resultsDf2, resultsDf5, resultsDf8)

resultazos <- data.frame(EFFECT_SIZE=character(0), SAMPLE_SIZE=numeric(0),PSIGNIFICANT=numeric(0), PNEGATIVE=numeric(0), PPOSITIVE=numeric(0))

by(resultsDf, resultsDf$EFFECT_SIZE, function(subResultsDfFirst){
  by(subResultsDfFirst, subResultsDfFirst$SAMPLE_SIZE, function(subResultsDf){
    numSignificant <- length(which(subResultsDf$P==0))/length(subResultsDf$P)
    numNegative <- length(which(subResultsDf$ESTIMATE<0))/length(subResultsDf$ESTIMATE)
    numPositive <- length(which(subResultsDf$ESTIMATE>0))/length(subResultsDf$ESTIMATE)
    rowTemp <- data.frame(EFFECT_SIZE=subResultsDf$EFFECT_SIZE[1], SAMPLE_SIZE=subResultsDf$SAMPLE_SIZE[1],PSIGNIFICANT=numSignificant, PNEGATIVE=numNegative, PPOSITIVE=numPositive)
    resultazos <<- rbind(resultazos, rowTemp)
  })
})
resultazos$EFFECT_SIZE <- as.factor(resultazos$EFFECT_SIZE/10)

p<-ggplot(resultazos, aes(x=SAMPLE_SIZE, y=PSIGNIFICANT, group=EFFECT_SIZE))+
  geom_line(aes(color=EFFECT_SIZE))+ labs(color='True effect size') +
  geom_point(aes(color=EFFECT_SIZE)) + ylab("Percentage of significant results") + xlab("Sample size")  + scale_x_continuous(breaks=seq(from = 4, to = 148, by = 16)) +  theme_bw()

plot(p)


# Just with significant
dfTricked <- data.frame(EFFECT_SIZE=character(0), SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0), NUM_EXP=numeric(0), PIDENTICAL=numeric(0),PIDENTICALNEGATIVE=numeric(0), PIDENTICALPOSITIVE=numeric(0))
# Calculate probability that majority non-significant in N experiments
by(resultazos, 1:nrow(resultazos), function(rowTemp) {
  # print(rowTemp["PNONSIGNIFICANT"])
  probabilitySignificant <- rowTemp$PSIGNIFICANT
  probabilityNegative <- rowTemp$PNEGATIVE
  probabilityPositive <- rowTemp$PPOSITIVE
  for(numExperiments in 2:12){
    print(numExperiments)
    probabilityIdenticalPositive <- probabilityPositive^numExperiments
    probabilityIdenticalNegative <- probabilityNegative^numExperiments
    probabilityIdentical <- probabilitySignificant^numExperiments
    rowTempazo <- data.frame(EFFECT_SIZE=rowTemp$EFFECT_SIZE, SAMPLE_SIZE=rowTemp$SAMPLE_SIZE,NUM_EXP=numExperiments,PIDENTICAL=probabilityIdentical, PIDENTICALNEGATIVE=probabilityIdenticalNegative, PIDENTICALPOSITIVE=probabilityIdenticalPositive)
    dfTricked <<- rbind(dfTricked, rowTempazo)
  }
  print(paste("row"))
})

# For effect size 0.2
resultazosEffectsize2 <- subset(dfTricked, EFFECT_SIZE==0.2)
resultazosEffectsize2$SAMPLE_SIZE <- as.factor(resultazosEffectsize2$SAMPLE_SIZE)

p<-ggplot(resultazosEffectsize2, aes(x=NUM_EXP, y=PIDENTICAL, group=SAMPLE_SIZE))+
  geom_line(aes(color=SAMPLE_SIZE))+ labs(color='Sample size') +
  geom_point(aes(color=SAMPLE_SIZE)) + ylab("Probability of all significant results") +xlab("Number of replications") + guides(fill=guide_legend(title="Num Experiments"))  + scale_x_continuous(breaks=seq(2,12,1))  + theme_bw()

plot(p)

# For effect size 0.5
resultazosEffectsize5 <- subset(dfTricked, EFFECT_SIZE==0.5)
resultazosEffectsize5$SAMPLE_SIZE <- as.factor(resultazosEffectsize5$SAMPLE_SIZE)

p<-ggplot(resultazosEffectsize5, aes(x=NUM_EXP, y=PIDENTICAL, group=SAMPLE_SIZE))+
  geom_line(aes(color=SAMPLE_SIZE))+ labs(color='Sample size') +
  geom_point(aes(color=SAMPLE_SIZE)) + ylab("Probability of all significant results") +xlab("Number of replications") + guides(fill=guide_legend(title="Num Experiments"))  + scale_x_continuous(breaks=seq(2,12,1))  + theme_bw()

plot(p)

# For effect size 0.8
resultazosEffectsize8 <- subset(dfTricked, EFFECT_SIZE==0.8)
resultazosEffectsize8$SAMPLE_SIZE <- as.factor(resultazosEffectsize8$SAMPLE_SIZE)

p<-ggplot(resultazosEffectsize8, aes(x=NUM_EXP, y=PIDENTICAL, group=SAMPLE_SIZE))+
  geom_line(aes(color=SAMPLE_SIZE))+ labs(color='Sample size') +
  geom_point(aes(color=SAMPLE_SIZE)) + ylab("Probability of all significant results") +xlab("Number of replications") + guides(fill=guide_legend(title="Num Experiments"))  + scale_x_continuous(breaks=seq(2,12,1))  + theme_bw()

plot(p)

######################################################################################################################
############################################## 3. Simulations with effect sizes signs and magnitudes #################
######################################################################################################################

print("Effect size signs and magnitudes")

effSizes <- c(2,5,8)
# numSubjectsPerArm <- seq(from = 10, to = 100, by =10)
# Following Dyba
numSubjectsPerArm <- seq(from = 2, to = 74, by = 8)
maxSimulations=5000

resultsDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), P=numeric(0), SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0))
for(effSize in effSizes){
  for(numSubjects in numSubjectsPerArm){
    print(numSubjects)
    for(i in 1:maxSimulations){
      d <- generateRandomData("between", i, numSubjects, 0, 0, meanTDD=(50+effSize), sdTDD=10, meanITL=50, sdITL=10)
      res<-cohenD(d)
      estimate<-res[1]
      # if lower and upper bound have the same sign then statistically significant
      pValue<-2
      if ((sign(res[3])==sign(res[4]))){
        pValue<-0  
      }else{
        pValue<-1
      }
      rowTemp <- data.frame(EXPERIMENT_ID=paste(effSize,"_",numSubjects,"_",i, sep=""), ESTIMATE=estimate, P=pValue, SAMPLE_SIZE=numSubjects*2, EFFECT_SIZE=effSize)
      resultsDf <- rbind(resultsDf, rowTemp)
    }
  }
  write.csv(resultsDf, paste("intermediateb",effSize,".csv", sep=""))
  resultsDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), P=numeric(0), SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0))
}

resultsDf2 <- read.csv("intermediateb2.csv")
resultsDf5 <- read.csv("intermediateb5.csv")
resultsDf8 <- read.csv("intermediateb8.csv")

resultsDf <- rbind(resultsDf2,resultsDf5,resultsDf8)

obtainPlotForEffectSize <- function(resultsDf, effectSize){
  effSizeSample10 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==4),"ESTIMATE"]
  effSizeSample20 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==20),"ESTIMATE"]
  effSizeSample30 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==36),"ESTIMATE"]
  effSizeSample40 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==52),"ESTIMATE"]
  effSizeSample50 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==68),"ESTIMATE"]
  effSizeSample60 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==84),"ESTIMATE"]
  effSizeSample70 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==100),"ESTIMATE"]
  effSizeSample80 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==116),"ESTIMATE"]
  effSizeSample90 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==132),"ESTIMATE"]
  effSizeSample100 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectSize & resultsDf$SAMPLE_SIZE==148),"ESTIMATE"]
  dat <- data.frame(dens = c(effSizeSample10,effSizeSample20,effSizeSample30,effSizeSample40,effSizeSample50,effSizeSample60,effSizeSample70, effSizeSample80,effSizeSample90,effSizeSample100)
                  ,lines = rep(c("4","20","36","52","68","84","100","116","132","148"), each = length(effSizeSample10)))

  dat$lines <- ordered(dat$lines, levels = c("4","20","36","52","68","84","100","116","132","148"))
  if (effectSize==2){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.2), labels=c("0","0.2"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density") + guides(fill=guide_legend(title="Sample size")) + theme_bw() + geom_vline(xintercept = 0.2, linetype="solid", color = "black") 
  }else if (effectSize==5){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.5), labels=c("0","0.5"), limits= c(-3, 3.5)) + scale_y_continuous(name = "Density") + guides(fill=guide_legend(title="Sample size")) + theme_bw() + geom_vline(xintercept = 0.5, linetype="solid", color = "black")
  }else{
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.8), labels=c("0","0.8"), limits= c(-3, 3.8)) + scale_y_continuous(name = "Density") + guides(fill=guide_legend(title="Sample size")) + theme_bw() + geom_vline(xintercept = 0.8, linetype="solid", color = "black")
  }
  return(p)
}



summaryEstimates<-summaryBy(ESTIMATE ~ EFFECT_SIZE + SAMPLE_SIZE, FUN=funSummaryBy, data=resultsDf)
stargazer(summaryEstimates, summary=FALSE, rownames = FALSE, digits = 2)


p1<-obtainPlotForEffectSize(resultsDf,2)
p2<-obtainPlotForEffectSize(resultsDf,5)
p3<-obtainPlotForEffectSize(resultsDf,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))


######################################################################################################################
############################################## 4. Increase Precision of Estimates ####################################
######################################################################################################################

print("Increase precision of estimates")

effSizes <- c(2,5,8)
# numSubjectsPerArm <- seq(from = 10, to = 100, by =10)
# Following Dyba
numSubjectsPerArm <- seq(from = 2, to = 74, by = 8)
familySizes <- c(1, seq(from=2, to=12, by=2))
maxSimulations <- 1000

resultsDf <- data.frame(ESTIMATE=numeric(0),  SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0), FAMILY_SIZE=numeric(0), SE=numeric(0))
for(effSize in effSizes){
  for(numSubjects in numSubjectsPerArm){
    print(numSubjects)
    # for family's sample sizes
    for(familySize in familySizes){
      for(i in 1:maxSimulations){
        resultsFamilyDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), VI=numeric(0))
        #simulate as many experiments as required
        for(j in 1:familySize){
          d <- generateRandomData("between", j, numSubjects, 0, 0, meanTDD=(50+effSize), sdTDD=10, meanITL=50, sdITL=10)
          res<-cohenD(d)
          estimate<-res[1]
          # pass the variance!
          vi <- res[2]
          rowTempazo <- data.frame(EXPERIMENT_ID=j, ESTIMATE=estimate, VI=vi)
          resultsFamilyDf <- rbind(resultsFamilyDf, rowTempazo)
        }
        meta.FE = rma(yi=resultsFamilyDf$ESTIMATE,vi=resultsFamilyDf$VI, method="REML", data=resultsFamilyDf, control=list(maxiter=100000))
        estimateAD <- meta.FE$b[1]
        seAD <- meta.FE$se[1]
        rowTemp <- data.frame(EFFECT_SIZE=effSize, FAMILY_SIZE=familySize, SAMPLE_SIZE=numSubjects*2, ESTIMATE=estimateAD, SE=seAD )
        resultsDf <- rbind(resultsDf,rowTemp)
      }
    }
  }
  write.csv(resultsDf, paste("intermediatec",effSize,".csv", sep=""))
  resultsDf <- data.frame(ESTIMATE=numeric(0),  SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0), FAMILY_SIZE=numeric(0), SE=numeric(0))
}

resultsDf2 <- read.csv("intermediatec2.csv")
resultsDf5 <- read.csv("intermediatec5.csv")
resultsDf8 <- read.csv("intermediatec8.csv")

resultsDf <- rbind(resultsDf2,resultsDf5,resultsDf8)

# For certain sample sizes (4, 36, 100)
resultsDfA<- subset(resultsDf, SAMPLE_SIZE==4)
resultsDfB<- subset(resultsDf, SAMPLE_SIZE==36)
resultsDfC<- subset(resultsDf, SAMPLE_SIZE==100)

obtainPlotForEffectSize <- function(resultsDf, effectsize){
  effSizeSample1 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==1),"ESTIMATE"]
  effSizeSample2 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==2),"ESTIMATE"]
  effSizeSample4 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==4),"ESTIMATE"]
  effSizeSample6 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==6),"ESTIMATE"]
  effSizeSample8 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==8),"ESTIMATE"]
  effSizeSample10 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==10),"ESTIMATE"]
  effSizeSample12 <- resultsDf[which(resultsDf$EFFECT_SIZE==effectsize & resultsDf$FAMILY_SIZE==12),"ESTIMATE"]
  
  dat <- data.frame(dens = c(effSizeSample1,effSizeSample2,effSizeSample4,effSizeSample6,effSizeSample8, effSizeSample10, effSizeSample12)
                    ,lines = as.factor(rep(c("1","2","4","6","8","10","12"), each = length(effSizeSample1))))
  
  dat$lines <-factor(dat$lines, levels=c("1","2","4","6","8","10","12")) 
  
  
  if (effectsize==2){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.2), labels=c("0","0.2"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Group size")) + theme_bw() + geom_vline(xintercept = 0.2, linetype="solid", color = "black") 
  }else if (effectsize==5){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.5), labels=c("0","0.5"), limits = c(-3, 3.5)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Group size")) + theme_bw() + geom_vline(xintercept = 0.5, linetype="solid", color = "black") 
  }else{
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.8), labels=c("0","0.8"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Group size")) + theme_bw() + geom_vline(xintercept = 0.8, linetype="solid", color = "black") 
  }
  
  
  return(p)
}

summaryEstimates<-summaryBy(ESTIMATE ~ EFFECT_SIZE + FAMILY_SIZE, FUN=funSummaryBy, data=resultsDfA)
stargazer(summaryEstimates, summary=FALSE, rownames = FALSE, digits = 2)

p1<-obtainPlotForEffectSize(resultsDfA,2)
p2<-obtainPlotForEffectSize(resultsDfA,5)
p3<-obtainPlotForEffectSize(resultsDfA,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))




summaryEstimates<-summaryBy(ESTIMATE ~ EFFECT_SIZE + FAMILY_SIZE, FUN=funSummaryBy, data=resultsDfB)
stargazer(summaryEstimates, summary=FALSE, rownames = FALSE, digits = 2)

p1<-obtainPlotForEffectSize(resultsDfB,2)
p2<-obtainPlotForEffectSize(resultsDfB,5)
p3<-obtainPlotForEffectSize(resultsDfB,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))



summaryEstimates<-summaryBy(ESTIMATE ~ EFFECT_SIZE + FAMILY_SIZE, FUN=funSummaryBy, data=resultsDfC)
stargazer(summaryEstimates, summary=FALSE, rownames = FALSE, digits = 2)

p1<-obtainPlotForEffectSize(resultsDfC,2)
p2<-obtainPlotForEffectSize(resultsDfC,5)
p3<-obtainPlotForEffectSize(resultsDfC,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))


######################################################################################################################
############################################## 5. Big experiment vs. chunks ##########################################
######################################################################################################################

print("Big experiment vs. chunks")

effSizes <- c(2, 5, 8)
# For example, 72 subjects
numSubjectsPerArm <- c(72)
numChunks <- c(4, 8, 12)
maxSimulations <- 5000


resultsDf <- data.frame(SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0), NUM_CHUNKS=numeric(0), ESTIMATE_OVERALL=numeric(0), VI_OVERALL=numeric(0), ESTIMATE_CHUNKS=numeric(0), VI_CHUNKS=numeric(0))
for(effSize in effSizes){
  for(numSubjects in numSubjectsPerArm){
    for(i in 1:maxSimulations){
      # Generate the experiment with all the subjects
      d <- generateRandomData("between", i, numSubjects, 0, 0, meanTDD=(50+effSize), sdTDD=10, meanITL=50, sdITL=10)
      res<-cohenD(d)
      # Obtain Cohen's d and the variance
      estimate<-res[1]
      vi <- res[2]
      # Position for those subjects with TDD
      middleIndex <- numSubjects
      # Divide the experiment in 2 and 3 chunks
      for(numChunk in numChunks){
        # obtain the number of subjects per chunk
        numSubjectsPerChunk=numSubjects/numChunk
        # Index for the min and max in each iteration
        minIndex <- 1
        # create dataframe for meta-analysis
        indexChunk <- 1
        resultsFamilyDf <- data.frame(EXPERIMENT_ID=character(0), ESTIMATE=numeric(0), VI=numeric(0))
        while(indexChunk<=numChunk){
          chunkDf <- d[c((minIndex:(minIndex+numSubjectsPerChunk-1)), (minIndex + middleIndex):(minIndex+ middleIndex + numSubjectsPerChunk -1)), ]
          resazo <- cohenD(chunkDf)
          estimateChunk <- resazo[1]
          viChunk <- resazo[2]
          rowTempazo <- data.frame(EXPERIMENT_ID=indexChunk, ESTIMATE=estimateChunk, VI=viChunk)
          resultsFamilyDf <- rbind(resultsFamilyDf, rowTempazo)
          # Next group of developers
          minIndex <- minIndex + numSubjectsPerChunk
          indexChunk <- indexChunk + 1
          # print(chunkDf)
        }
        # print(resultsFamilyDf)
        meta.FE <- rma(yi=resultsFamilyDf$ESTIMATE,vi=resultsFamilyDf$VI, method="FE", data=resultsFamilyDf, control=list(maxiter=100000))
        # print(meta.FE)
        estimateAD <- meta.FE$b[1]
        viAD <- meta.FE$se[1]^2
        # print(estimate)
        # print(estimateAD)
        # print("*****")
        rowTemp <- data.frame(SAMPLE_SIZE=(numSubjects*2), EFFECT_SIZE=effSize, NUM_CHUNKS=numChunk, ESTIMATE_OVERALL=estimate, VI_OVERALL=vi, ESTIMATE_CHUNKS=estimateAD, VI_CHUNKS=viAD)
        resultsDf <- rbind(resultsDf,rowTemp)
      }
    }
  }
  write.csv(resultsDf, paste("intermediated",effSize,".csv", sep=""))
  resultsDf <- data.frame(SAMPLE_SIZE=numeric(0), EFFECT_SIZE=numeric(0), NUM_CHUNKS=numeric(0), ESTIMATE_OVERALL=numeric(0), VI_OVERALL=numeric(0), ESTIMATE_CHUNKS=numeric(0), VI_CHUNKS=numeric(0))
}


resultsDf2 <- read.csv("intermediated2.csv")
resultsDf5 <- read.csv("intermediated5.csv")
resultsDf8 <- read.csv("intermediated8.csv")

resultsDf <- rbind(resultsDf2,resultsDf5,resultsDf8)


obtainPlotChunks <- function(resultsDf, effectsize){
  allCurve <- subset(resultsDf,NUM_CHUNKS==4)$ESTIMATE_OVERALL
  chunkCurve2 <- subset(resultsDf,NUM_CHUNKS==4)$ESTIMATE_CHUNKS
  chunkCurve3 <- subset(resultsDf,NUM_CHUNKS==8)$ESTIMATE_CHUNKS
  chunkCurve4 <- subset(resultsDf,NUM_CHUNKS==12)$ESTIMATE_CHUNKS
  
  dat <- data.frame(dens = c(allCurve,chunkCurve2, chunkCurve3, chunkCurve4)
                    ,lines = as.factor(rep(c("Chunks == 1","Chunks == 4","Chunks == 8","Chunks == 12"), each = length(allCurve))))
  
  dat$lines <-factor(dat$lines, levels=c("Chunks == 1","Chunks == 4","Chunks == 8","Chunks == 12"))
  
  
  if (effectsize==2){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.2), labels=c("0","0.2"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.2, linetype="solid", color = "black") 
    p <- p + facet_wrap(~ lines, ncol=1) + scale_fill_discrete(labels = c("1", "4", "8", "12"))
  }else if (effectsize==5){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.5), labels=c("0","0.5"), limits = c(-3, 3.5)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.5, linetype="solid", color = "black") 
    p <- p + facet_wrap(~ lines, ncol=1) + scale_fill_discrete(labels = c("1", "4", "8", "12"))
  }else{
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.8), labels=c("0","0.8"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.8, linetype="solid", color = "black") 
    p <- p + facet_wrap(~ lines, ncol=1) + scale_fill_discrete(labels = c("1", "4", "8", "12"))
  }
  
  
  return(p)
}


# For certain effect sizes (2, 5, 8)
resultsDfA<- subset(resultsDf, EFFECT_SIZE==2)
resultsDfB<- subset(resultsDf, EFFECT_SIZE==5)
resultsDfC<- subset(resultsDf, EFFECT_SIZE==8)

p1<-obtainPlotChunks(resultsDfA,2)
p2<-obtainPlotChunks(resultsDfB,5)
p3<-obtainPlotChunks(resultsDfC,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))


obtainPlotChunks <- function(resultsDf, effectsize){
  allCurve <- subset(resultsDf,NUM_CHUNKS==4)$ESTIMATE_OVERALL
  chunkCurve2 <- subset(resultsDf,NUM_CHUNKS==4)$ESTIMATE_CHUNKS
  chunkCurve3 <- subset(resultsDf,NUM_CHUNKS==8)$ESTIMATE_CHUNKS
  chunkCurve4 <- subset(resultsDf,NUM_CHUNKS==12)$ESTIMATE_CHUNKS
  
  dat <- data.frame(dens = c(allCurve,chunkCurve2, chunkCurve3, chunkCurve4)
                    ,lines = as.factor(rep(c("1","4","8","12"), each = length(allCurve))))
  
  dat$lines <-factor(dat$lines, levels=c("1","4","8","12"))

  
  if (effectsize==2){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.2), labels=c("0","0.2"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.2, linetype="solid", color = "black") 
  }else if (effectsize==5){
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.5), labels=c("0","0.5"), limits = c(-3, 3.5)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.5, linetype="solid", color = "black") 
  }else{
    p<- ggplot(dat, aes(x = dens, fill = lines, order = as.numeric(lines))) + geom_density(alpha = 0.2) + scale_x_continuous("Effect size", breaks=c(0, 0.8), labels=c("0","0.8"), limits = c(-3, 3.2)) + scale_y_continuous(name = "Density", limits=c(0,4)) + guides(fill=guide_legend(title="Chunks")) + theme_bw() + geom_vline(xintercept = 0.8, linetype="solid", color = "black") 
  }
  
  
  return(p)
}


p1<-obtainPlotChunks(resultsDfA,2)
p2<-obtainPlotChunks(resultsDfB,5)
p3<-obtainPlotChunks(resultsDfC,8)

legend <- get_legend(p1)
p1<- p1 + theme(legend.position="none")  
p2<- p2 + theme(legend.position="none")  
p3<- p3 + theme(legend.position="none")  
p<-plot_grid(p1,p2,p3, ncol=3)
plot(plot_grid(p,legend, rel_widths = c(0.9,0.1)))



##Initializing the library
library(ggplot2)
library(MASS)
library(dplyr)
library(klaR)
library(caret)
source("http://bit.ly/theme_pub")
theme_set(theme_pub())

##Data loading
lythrum_data <- read.csv("/Users/sreevats/Downloads/Sima-paper_1-16s rRNA-sequences/ColauttiBarrett2013Data.csv", header = T)

head(lythrum_data)
str(lythrum_data)

##Data Inspection
ggplot(aes(x=Flwr07, data=lythrum_data) + geom_histogram(bins = 30)) 

## data looks more or less similar to normality

##Now separating the predictors and responses

pred_lythrum_data <- lythrum_data %>%
  dplyr::select(1:7)

resp_lythrum_data <- lythrum_data %>%
  dplyr::select(-c(1:7))

##Scaling of data

scaled <- resp_lythrum_data %>%
  mutate_all(scale)

## Eliminate the missing data 

scaled %>%
  select_if(function(x) any(is.na(x))) %>%
  names()

ScalComp <- scaled %>%
  mutate (Flwr07 = ifelse(is.na(Flwr07),0,Flwr07), FVeg07 = ifelse(is.na(FVeg07),0,FVeg07), InfMass07 = ifelse(is.na(InfMass07),0,InfMass07),
          Flwr08 = ifelse(is.na(Flwr08),0,Flwr08), FVeg08 = ifelse(is.na(FVeg08),0,FVeg08), InfMass08 = ifelse(is.na(InfMass08),0,InfMass08), HVeg08 = ifelse(is.na(HVeg08),0,HVeg08),
          Flwr09 = ifelse(is.na(Flwr09),0,Flwr09), FVeg09 = ifelse(is.na(FVeg09),0,FVeg09), InfMass09 = ifelse(is.na(InfMass09),0,InfMass09), HVeg09 = ifelse(is.na(HVeg09),0,HVeg09),
          Flwr10 = ifelse(is.na(Flwr10),0,Flwr10), FVeg10 = ifelse(is.na(FVeg10),0,FVeg10), InfMass10 = ifelse(is.na(InfMass10),0,InfMass10), HVeg10 = ifelse(is.na(HVeg10),0,HVeg10))

mean(ScalComp$Flwr07)
sd(ScalComp$Flwr07)

qplot(x=Flwr07,data = ScalComp)

##Dimension Reduction

dim(ScalComp)

##Running LDA model for genetic population

lda_gen_pop <- lda(x=ScalComp,grouping=pred_lythrum_data$Pop)
summary(lda_gen_pop)
lda_gen_pop$counts

##what does scaling show?

lda_gen_pop$scaling

##predict function

Pred_gen_pop <- predict(lda_gen_pop)
summary(Pred_gen_pop)

##confusion matrix

CatDat_gen_pop <- data.frame(Pop=as.factor(pred_lythrum_data$Pop), Predicted=Pred_gen_pop$class)
table(CatDat_gen_pop)

##Running LDA model for Regionr

lda_region <- lda(x=ScalComp,grouping=pred_lythrum_data$Region)
summary(lda_region)
lda_region$counts

##what does scaling show?

lda_region$scaling

##predict function

Pred_region <- predict(lda_region)
summary(Pred_region)

##confusion matrix

CatDat_region <- data.frame(Pop=as.factor(pred_lythrum_data$Region), Predicted=Pred_region$class)
table(CatDat_region)

##Running LDA model for site

lda_site <- lda(x=ScalComp,grouping=pred_lythrum_data$Site)
summary(lda_site)
lda_site$counts

##what does scaling show?

lda_site$scaling

##predict function

Pred_site <- predict(lda_site)
summary(Pred_site)

##confusion matrix

CatDat_site <- data.frame(Pop=as.factor(pred_lythrum_data$Site), Predicted=Pred_site$class)
table(CatDat_site)


###trying RDA for the data wrt population

##train control

CTL<-trainControl(method="repeatedcv",
                  number=4,
                  repeats=24,
                  classProbs=T,
                  verboseIter=F,
                  search="random")
set.seed(123)
randomRDA <- train(x=ScalComp, y=pred_lythrum_data$Pop, method = "rda", metric = "Accuracy", tuneLength = 24, trControl = CTL)

randomRDA

ggplot(randomRDA) + theme(legend.position="bottom")

RDAmod<-rda(x=ScalComp,grouping=pred_lythrum_data$Pop,
            regularization=c(gamma=0.5, lambda=0.5))
summary(RDAmod)

Pred_rda <- predict(RDAmod, ScalComp)
Pred_rda

Post<-as.data.frame(Pred_gen_pop$posterior)
head(Post)
Post



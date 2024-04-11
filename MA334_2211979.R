library(dplyr)
library(tidyr) 
library(moments)
library(reshape2)

par(mfrow=c(1, 1)) 

#Read the Data
Proj_data <-  read.csv("proportional_species_richness_V3.csv") # you use V2 or V3

Proj_data$period <- as.factor(Proj_data$period)  

Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)

names(Proj_data)

# 7Allocated Taxonomic Groups

#All Taxonomic Groups
all <- c(2:12)
#Selected Taxonomic Groups
eco_selected <- c(2,3,4,5,8,9,10)

eco_not_selected <- all[!(all%in%eco_selected)]
eco_names <- names(Proj_data[,2:12])
eco_selected_names <- names(Proj_data)[eco_selected]
eco_selected_names

# calculate the bio div measure over 7 taxonomic groups
mean_selected <- rowMeans(Proj_data[,eco_selected],na.rm=TRUE) # mean the 7 columns 
# Check the mean value
sum(is.na(mean_selected)) 

# add in the biodiversity measure which is the mean over 7 taxonomic groups
Proj_data_MA334 <- Proj_data%>%mutate(eco_status_7=mean_selected)
names(Proj_data_MA334)

# the data exploration phase (only some suggested approaches)
summary(Proj_data)
# you could split the data by period and compare these stats before and after 
summary_stats <- sapply(Proj_data[,eco_selected], function(x) c(mean=mean(x,na.rm=TRUE), sd=sd(x, na.rm=TRUE), skew=skewness(x,na.rm=TRUE)))

# doing a linear regression with only Northing as a predictor 
#
lin_mod <- lm(Proj_data_MA334$eco_status_7~Proj_data$Northing)
summary(lin_mod)
qqnorm(lin_mod$residuals)
qqline(lin_mod$residuals,col="red")

# following code splits between the two periods to find the BD7 change
# however it may be better to use period as a predictor 

# box plot comparisons for the two periods ignoring all other varaibles 
eco_status <- Proj_data_MA334%>%pull(eco_status_7)
eco_period <- Proj_data_MA334%>%pull(period)
plot(eco_status~eco_period,xlab="Ecological Periods",ylab = "Ecological Status", main="Period Vs Status")

names(Proj_data_MA334)
Proj_data_MA334_period <- Proj_data_MA334%>%select(Location,period,eco_status_7)
Proj_data_MA334_split <- Proj_data_MA334_period%>%pivot_wider(names_from =period,values_from=eco_status_7)
Proj_data_MA334_split <- Proj_data_MA334_split%>%mutate(BD7_change=Y00-Y70)
head(Proj_data_MA334_split)
hist(Proj_data_MA334_split$BD7_change,xlab = "BD7",main = "Histogram of BD7 Changes")  # the distribution of the BD7 change 

BD7_change <- Proj_data_MA334_split%>%pull(BD7_change)
t.test(BD7_change,mu=0)  # t test with H0: mu=0

# comparing the two distributions of bio div based on 7 and 11 taxonomic groups 

par(mfrow=c(1, 1))  # divide graph area in 1 columns

qqplot(Proj_data_MA334$eco_status_7,Proj_data_MA334$ecologicalStatus,main="QQ Plot",xlab = "Ecological Status7",ylab = "Ecological Status")
abline(0,1,col="red")

# both cdfs together  and do a kolmogorov test H0: distributions are the same
BD7_cdf <- ecdf(Proj_data_MA334$eco_status_7)
BD11_cdf <- ecdf(Proj_data_MA334$ecologicalStatus)

plot(BD11_cdf,col="red",main="Kolmogorov Test for Ecological Status")
lines(BD7_cdf,col="green")

ks.test(Proj_data_MA334$eco_status_7,Proj_data_MA334$ecologicalStatus)

# Simple linear regression part of the specified assignment

# regressions of eco_status_7 against ecologicalstatus based on all 11

plot(Proj_data_MA334$eco_status_7~Proj_data_MA334$ecologicalStatus,main="Scatter Plot with Regrssion line of Simple Linear Regrssion",xlab="Ecological Status",ylab ="Ecological Status of 7 Groups")
abline(0,1,col="red")

lin_mod <- lm(Proj_data_MA334$eco_status_7~Proj_data_MA334$ecologicalStatus)
summary(lin_mod)

abline(lin_mod,col="green")

qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")

# do the same for each period report and differences 

Proj_data_MA334_Y70 <- Proj_data_MA334%>%filter(period=="Y70")

lin_mody70 <- lm(Proj_data_MA334_Y70$eco_status_7~Proj_data_MA334_Y70$ecologicalStatus)
lin_mody70$coefficients
plot(lin_mody70)
# for later period 

Proj_data_MA334_Y00 <- Proj_data_MA334%>%filter(period=="Y00")
lin_mody00 <- lm(Proj_data_MA334_Y00$eco_status_7~Proj_data_MA334_Y00$ecologicalStatus)
lin_mody00$coefficients
plot(lin_mody00)

# linear regression of BD4 on BD7 
mean_selected <- rowMeans(Proj_data[,eco_not_selected ],na.rm=TRUE) # mean the rem 4 columns 
sum(is.na(mean_selected)) # check that there are no NAs in mean_selected
# add in the biodiversity measure which is the mean over 7 taxonomic groups
Proj_data_MA334 <- Proj_data_MA334%>%mutate(eco_status_4=mean_selected)
names(Proj_data_MA334)

# now multiple linear regression BD4 against the selected 7 

# Create Training and Test data 
trainingRowIndex <- sample(1:nrow(Proj_data_MA334), 0.8*nrow(Proj_data_MA334))  # row indices for 80% training data
trainingData <- Proj_data_MA334[trainingRowIndex, ]  # model training data
testData  <- Proj_data_MA334[-trainingRowIndex, ]%>%na.omit # for test data remove NAs 

# Build the model on training data
lmMod_train <- lm(eco_status_4~.,
                  data=trainingData[c(eco_selected_names,"eco_status_4")],
                  na.action=na.omit,y=TRUE)
summary(lmMod_train)
cor(lmMod_train$fitted.values,lmMod_train$y) # cor training data 


# Feature engineering

lmModdel_reduced_train2 <- lm(eco_status_4~.,
                              data=trainingData[c(eco_selected_names[-7],"eco_status_4","period","Northing")],
                              na.action=na.omit,y=TRUE)

summary (lmModdel_reduced_train2)
red_mod_err<-cor(lmModdel_reduced_train2$fitted.values,lmModdel_reduced_train2$y)

# compare the effect of each significant coefficient to that of period
mult_lin_mod$coefficients
as.numeric(mult_lin_mod$coefficients[3])*mean(Proj_data_MA334$Easting)
as.numeric(mult_lin_mod$coefficients[4])*mean(Proj_data_MA334$Northing)

#AIC
trainmodels<- AIC(lmMod_train)
reducemodel<-AIC(lmModdel_reduced_train2)
AICVlaues <- c(trainmodels,reducemodel)

# The following PCA method is an extension to the set book 
# PCA for visualizing the multi-dimensional spread of biodiversity values #######################

table(Proj_data_MA334_Y70$period); table(Proj_data_MA334_Y00$period) # check that these separate periods 
table(Proj_data_MA334_Y00$Location==Proj_data_MA334_Y70$Location) # check that Locations correspond between the two periods

eco_difference <- Proj_data_MA334_Y00[,eco_selected ]-Proj_data_MA334_Y70[,eco_selected ] # general differences between the two periods 
head(eco_difference)

# see ?prcomp the default here is the mean correct but not to scale 
pr.out=prcomp(na.omit(eco_difference)) # Principal Components 
pr.out$center  # gives the mean corrections the "centers"
pr.out$scale  # not scaled
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], Proj_data_MA334_Y00$dominantLandClass, cex=0.5, pos=4, col="red") # location labels

# label by location 
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], Proj_data_MA334_Y00$Location, cex=0.4, pos=4, col="red") # location labels


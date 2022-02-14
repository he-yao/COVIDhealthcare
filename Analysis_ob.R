# title: "Shandong Data Analysis - OB & GYN only"
# author: "Yao He"
# date: "09/17/2021"

# Load packages
library(knitr)
library(plyr)
library(readr)
library(dplyr) # mutate(), rename()
library(ggplot2)
library(lattice)
library(nlme)
library(lme4)
library(lspline)
library(stargazer)
library(UWbe536)
library(glarma)

library(sandwich)
library(MASS)
library(lmtest)
library(cowplot)
library(tseries) #adf.test Dickey-Fuller test; p<0.05, reject null that time series is not stationary

# Create the golden plot template
goldentheme <- theme(
  panel.background = element_rect(fill = "white"),
  aspect.ratio = ((1 + sqrt(5))/2)^(-1),
  axis.ticks.length = unit(0.5, "char"),
  axis.line.x.top = element_line(size = 0.2),
  axis.line.x.bottom = element_line(size = 0.2),
  axis.ticks.x = element_line(size = 0.2),
  axis.text.x = element_text(color = "black", size = 12),
  axis.title.x = element_text(size = 12,
                              margin = margin(t = 7.5, r = 0, b = 0, l = 0)),
  axis.ticks.y = element_blank(),
  axis.text.y = element_text(color = "black", size = 12,
                             margin = margin(t = 0, r = -4, b = 0, l = 0)),
  axis.title.y = element_text(size = 12,
                              margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
  legend.key = element_rect(fill = NA, color = NA),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "gray45", size = 0.2),
  strip.background = element_blank(),
  strip.text.x = element_text(size=12),
  strip.text.y = element_blank(), 
  strip.placement = "outside", 
  panel.spacing.x = unit(1.25, "lines"),
  panel.spacing.y = unit(1, "lines") 
)

class(goldentheme)

goldentheme_facet <- theme(
  panel.background = element_rect(fill = "white"),
  axis.ticks.length = unit(0.5, "char"),
  axis.line.x.top = element_line(size = 0.2),
  axis.line.x.bottom = element_line(size = 0.2),
  axis.ticks.x = element_line(size = 0.2),
  axis.text.x = element_text(color = "black", size = 12),
  axis.title.x = element_text(size = 12,
                              margin = margin(t = 7.5, r = 0, b = 0, l = 0)),
  axis.ticks.y = element_blank(),
  axis.text.y = element_text(color = "black", size = 12,
                             margin = margin(t = 0, r = -4, b = 0, l = 0)),
  axis.title.y = element_text(size = 12,
                              margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
  legend.key = element_rect(fill = NA, color = NA),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "gray45", size = 0.2),
  strip.background = element_blank(),
  strip.text.x = element_text(size=12),
  strip.text.y = element_blank(), 
  strip.placement = "outside", 
  panel.spacing.x = unit(1.25, "lines"),
  panel.spacing.y = unit(1, "lines") 
)

class(goldentheme_facet)

# create a function to replace any space in x axis labels with a new line
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}



# OBSTETRICS

## Outpatient

### Data & Descriptive

# read in data
datOS <- read_csv("Data/monthly_outpatient_surgery_inpatient.csv")

# subset for outpatient OBSTETRICS
datObOUT <- subset(datOS, grepl("产", datOS$keshi) & datOS$division=="outpatient")

# combine counts if a monthly count was divided among more than 1 subdivisions of the obstetrics department
datObOUT1 <- datObOUT %>% group_by(year, month) %>% summarise_at(vars(count), sum)

# create Time var: first ever month in dataset = 0
datObOUT1$time <- (datObOUT1$year-2015)*12+(datObOUT1$month-1)

# create Time1 var: 0 centered at Jan 2020 (1 time point before the pandemic actually happened in Feb 2020)
datObOUT1$time1 <- datObOUT1$time - 60

# create COVID1 var to indicate just Feb 2020 to make the drop in that month pop.
datObOUT1$COVID1 <- with(datObOUT1, ifelse(year == 2020 & month == 2, 1, 0))

# create COVID2 var to indicate the time points starting from Mar 2020. 
datObOUT1$COVID2 <- with(datObOUT1, ifelse(time1 >= 2, 1, 0))

# create Post var; 1 should start at April 2020
datObOUT1$post <- with(datObOUT1, ifelse(time1 >= 3, (year-2020)*12+(month-3), 0))

# initial chart of # cases by month; just raw data (2017-2020)
ggplot(data = datObOUT1, aes(x = time, y = count)) +
  geom_point() +
  geom_line() +
  ggtitle("Obstetric outpatient visits (Jan 2017-May 2021)") +
  labs(x = "Time", y = "Number of obstetric outpatient visits") +
  scale_x_continuous(limits=c(24,76), breaks=c(24,36,48,60,72,76),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","May 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=66.5, y=25000, label="COVID-19",size=5) +
  expand_limits(y = 0) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObOUT.17
pp.ObOUT.17

### Analysis

# adjust for monthly seasonality, January as the control，no offset for newborn population becuase 2021 data are not accurate
mod.ObOUT <- glm.nb(count ~ time1 + COVID1 + COVID2 + post 
                    + as.factor(month), 
                    data = subset(datObOUT1, year >= 2017))

lincom(mod.ObOUT)

# Extract the fitted values from this fitted model
datObOUT1$fit <- NA
datObOUT1$fit[datObOUT1$year>=2017] <- fitted(mod.ObOUT)

# Use AR(1)
coef.ObOUT <- coeftest(mod.ObOUT, vcov=NeweyWest(mod.ObOUT,lag=1,prewhite = F))

# put coef, CI, and p-value for time1, COVID, post together
m.ObOUT <- matrix(, nrow = 4, ncol = 4)
m.ObOUT[,1]<-exp(coef.ObOUT[c(2:5),1])
m.ObOUT[,2]<-exp(coef.ObOUT[c(2:5),1]-1.96*coef.ObOUT[c(2:5),2])
m.ObOUT[,3]<-exp(coef.ObOUT[c(2:5),1]+1.96*coef.ObOUT[c(2:5),2])
m.ObOUT[,4]<-coef.ObOUT[c(2:5),4]
rownames(m.ObOUT) <- c("time1", "COVID1","COVID2", "post")
colnames(m.ObOUT) <- c("exp(beta)", "95%low", "95%up", "p-value")
round(m.ObOUT, 5)

# Counterfactual or expected Y values based on the original model
X0.ObOUT <- X.ObOUT <- model.matrix(mod.ObOUT)
X0.ObOUT[,3:5] <- 0
datObOUT1$exp[datObOUT1$year>=2017] <- as.vector(exp(X0.ObOUT%*%coef.ObOUT)) # save expected Y values to the monthly dataset

####### 10,000 simulations
Sigma.ObOUT <- NeweyWest(mod.ObOUT,lag=1,prewhite = F)
set.seed(12)
beta.ObOUT <- mvrnorm(n = 10000, coef.ObOUT[,1], Sigma.ObOUT)

D.ObOUT <- Y_sum.ObOUT <- Y0_sum.ObOUT <- D_P.ObOUT <- RR.ObOUT <- NA

for (i in 1:10000){
  Y_sum.ObOUT[i] <- sum(as.vector(exp(X.ObOUT%*%beta.ObOUT[i,]))[38:53])
  Y0_sum.ObOUT[i] <- sum(as.vector(exp(X0.ObOUT%*%beta.ObOUT[i,]))[38:53])
  D.ObOUT[i] <- Y0_sum.ObOUT[i] - Y_sum.ObOUT[i]
  D_P.ObOUT[i] <- (Y0_sum.ObOUT[i] - Y_sum.ObOUT[i]) / Y0_sum.ObOUT[i]
  RR.ObOUT[i] <- Y_sum.ObOUT[i] / Y0_sum.ObOUT[i]
}

# Fitted post
mean(Y_sum.ObOUT)
quantile(Y_sum.ObOUT, c(0.025,0.975)) #95% CI

# Counterfactual post
mean(Y0_sum.ObOUT)
quantile(Y0_sum.ObOUT, c(0.025,0.975)) 

# !!! Lost visits (observed - counterfactual/forecasted)
mean(D.ObOUT) 
quantile(D.ObOUT, c(0.025,0.975))
p.D.ObOUT <- length(D.ObOUT[D.ObOUT < 0]) / length(D.ObOUT) # p-value; if p=0, write it as p<0.0001
p.D.ObOUT

mean(D_P.ObOUT)  ## percentage of decrease
quantile(D_P.ObOUT, c(0.025,0.975))
p.D_P.ObOUT <- length(D_P.ObOUT[D_P.ObOUT < 0]) / length(D_P.ObOUT) # p-value
p.D_P.ObOUT

mean(RR.ObOUT)  ## RR
quantile(RR.ObOUT, c(0.025,0.975))

### Chart for the model analysis
# CI around the line, just for the fitted line.
datObOUT1$fit_up[datObOUT1$year>=2017] <- exp(predict(mod.ObOUT) + 1.96*predict(mod.ObOUT, se.fit=T)$se.fit)
datObOUT1$fit_low[datObOUT1$year>=2017] <- exp(predict(mod.ObOUT) - 1.96*predict(mod.ObOUT,se.fit=T)$se.fit)

# chart of observed versus expected
ggplot(data = datObOUT1, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp, colour = "Expected")) +
  geom_line(aes(y = fit, colour = "Fitted")) + 
  geom_point(aes(y = count, colour = "Fitted")) +
  ggtitle("Obstetric outpatient visits: Observed versus expected (Jan 2017 - May 2021)") +
  labs(x = "Time", y = "Obstetric outpatient visits", color = "") +
  scale_x_continuous(limits=c(24,82), breaks=c(24,36,48,60,72,82),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Nov 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=67, y=11000, label="COVID-19",size=5) +
  expand_limits(y = c(8000:26000)) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObOUT.ove
pp.ObOUT.ove



## Surgery

### Data & Descriptive

# subset for surgery OBSTETRICS
datObSUR <- subset(datOS, grepl("产", datOS$keshi) & datOS$division=="surgery")

# combine counts if a monthly count was divided among more than 1 subdivisions of the obstetrics department
datObSUR1 <- datObSUR %>% group_by(year, month) %>% summarise_at(vars(count), sum)

# create Time var: first ever month in dataset = 0
datObSUR1$time <- (datObSUR1$year-2015)*12+(datObSUR1$month-1)

# create Time1 var: 0 centered at Jan 2020 (1 time point before the pandemic actually happened in Feb 2020)
datObSUR1$time1 <- datObSUR1$time - 60

# create COVID1 var to indicate just Feb 2020 to make the drop in that month pop.
datObSUR1$COVID1 <- with(datObSUR1, ifelse(year == 2020 & month == 2, 1, 0))

# create COVID2 var to indicate the time points starting from Mar 2020. 
datObSUR1$COVID2 <- with(datObSUR1, ifelse(time1 >= 2, 1, 0))

# create Post var; 1 should start at April 2020
datObSUR1$post <- with(datObSUR1, ifelse(time1 >= 3, (year-2020)*12+(month-3), 0))

# initial chart of # cases by month; just raw data (2017-2020)
ggplot(data = datObSUR1, aes(x = time, y = count)) +
  geom_point() +
  geom_line() +
  ggtitle("Obstetric surgeries (Jan 2017-May 2021)") +
  labs(x = "Time", y = "Number of obstetric surgeries") +
  scale_x_continuous(limits=c(24,76), breaks=c(24,36,48,60,72,76),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","May 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=66.5, y=1200, label="COVID-19",size=5) +
  expand_limits(y = 0) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObSUR.17
pp.ObSUR.17

### Analysis

# adjust for monthly seasonality, January as the control，no offset for newborn population becuase 2021 data are not accurate
mod.ObSUR <- glm.nb(count ~ time1 + COVID1 + COVID2 + post 
                    + as.factor(month), 
                    data = subset(datObSUR1, year >= 2017))

lincom(mod.ObSUR)

# Extract the fitted values from this fitted model
datObSUR1$fit <- NA
datObSUR1$fit[datObSUR1$year>=2017] <- fitted(mod.ObSUR)

# Use AR(1)
coef.ObSUR <- coeftest(mod.ObSUR, vcov=NeweyWest(mod.ObSUR,lag=1,prewhite = F))

# put coef, CI, and p-value for time1, COVID, post together
m.ObSUR <- matrix(, nrow = 4, ncol = 4)
m.ObSUR[,1]<-exp(coef.ObSUR[c(2:5),1])
m.ObSUR[,2]<-exp(coef.ObSUR[c(2:5),1]-1.96*coef.ObSUR[c(2:5),2])
m.ObSUR[,3]<-exp(coef.ObSUR[c(2:5),1]+1.96*coef.ObSUR[c(2:5),2])
m.ObSUR[,4]<-coef.ObSUR[c(2:5),4]
rownames(m.ObSUR) <- c("time1", "COVID1","COVID2", "post")
colnames(m.ObSUR) <- c("exp(beta)", "95%low", "95%up", "p-value")
round(m.ObSUR, 5)

# Counterfactual or expected Y values based on the original model
X0.ObSUR <- X.ObSUR <- model.matrix(mod.ObSUR)
X0.ObSUR[,3:5] <- 0
datObSUR1$exp[datObSUR1$year>=2017] <- as.vector(exp(X0.ObSUR%*%coef.ObSUR)) # save expected Y values to the monthly dataset

####### 10,000 simulations
Sigma.ObSUR <- NeweyWest(mod.ObSUR,lag=1,prewhite=F)
set.seed(12)
beta.ObSUR <- mvrnorm(n = 10000, coef.ObSUR[,1], Sigma.ObSUR)

D.ObSUR <- Y_sum.ObSUR <- Y0_sum.ObSUR <- D_P.ObSUR <- RR.ObSUR <- NA

for (i in 1:10000){
  Y_sum.ObSUR[i] <- sum(as.vector(exp(X.ObSUR%*%beta.ObSUR[i,]))[38:53])
  Y0_sum.ObSUR[i] <- sum(as.vector(exp(X0.ObSUR%*%beta.ObSUR[i,]))[38:53])
  D.ObSUR[i] <- Y0_sum.ObSUR[i] - Y_sum.ObSUR[i]
  D_P.ObSUR[i] <- (Y0_sum.ObSUR[i] - Y_sum.ObSUR[i]) / Y0_sum.ObSUR[i]
  RR.ObSUR[i] <- Y_sum.ObSUR[i] / Y0_sum.ObSUR[i]
}

# Fitted post
mean(Y_sum.ObSUR)
quantile(Y_sum.ObSUR, c(0.025,0.975)) #95% CI

# Counterfactual post
mean(Y0_sum.ObSUR)
quantile(Y0_sum.ObSUR, c(0.025,0.975)) 

# !!! Lost visits (observed - counterfactual/forecasted)
mean(D.ObSUR) 
quantile(D.ObSUR, c(0.025,0.975))
p.D.ObSUR <- length(D.ObSUR[D.ObSUR < 0]) / length(D.ObSUR) # p-value; if p=0, write it as p<0.0001
p.D.ObSUR

mean(D_P.ObSUR)  ## percentage of decrease
quantile(D_P.ObSUR, c(0.025,0.975))
p.D_P.ObSUR <- length(D_P.ObSUR[D_P.ObSUR < 0]) / length(D_P.ObSUR) # p-value
p.D_P.ObSUR

mean(RR.ObSUR)  ## RR
quantile(RR.ObSUR, c(0.025,0.975))

### Chart for the model analysis
# CI around the line, just for the fitted line.
datObSUR1$fit_up[datObSUR1$year>=2017] <- exp(predict(mod.ObSUR) + 1.96*predict(mod.ObSUR, se.fit=T)$se.fit)
datObSUR1$fit_low[datObSUR1$year>=2017] <- exp(predict(mod.ObSUR) - 1.96*predict(mod.ObSUR,se.fit=T)$se.fit)

# chart of observed versus expected
ggplot(data = datObSUR1, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp, colour = "Expected")) +
  geom_line(aes(y = fit, colour = "Fitted")) + 
  geom_point(aes(y = count, colour = "Fitted")) +
  ggtitle("Obstetric surgeries: Observed versus expected (Jan 2017 - May 2021)") +
  labs(x = "Time", y = "Obstetric surgeries", color = "") +
  scale_x_continuous(limits=c(24,82), breaks=c(24,36,48,60,72,82),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Nov 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=67, y=1050, label="COVID-19",size=5) +
  expand_limits(y = c(500:1200)) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObSUR.ove
pp.ObSUR.ove


## Inpatient

### Data & Descriptive

# read in data
datOB <- read_csv("Data/allob.csv")

# create Time var: first ever month in dataset = 0
datOB$time <- (datOB$year-2015)*12+(datOB$month-1)

# create Time1 var: 0 centered at Jan 2020 (1 time point before the pandemic actually happened in Feb 2020)
datOB$time1 <- datOB$time - 60

# create COVID1 var to indicate just Feb 2020 to make the drop in that month pop.
datOB$COVID1 <- with(datOB, ifelse(year == 2020 & month == 2, 1, 0))

# create COVID2 var to indicate the time points starting from Mar 2020. 
datOB$COVID2 <- with(datOB, ifelse(time1 >= 2, 1, 0))

# create Post var; 1 should start at April 2020
datOB$post <- with(datOB, ifelse(time1 >= 3, (year-2020)*12+(month-3), 0))

# change wrong age data
datOB$age <- with(datOB, ifelse(age==1000, 0, age))

# exclude newborns (age==0) from the data
datOBp <- subset(datOB, age>0)


# Summarize monthly counts
datOBp %>% 
  group_by(time, year, month, time1, COVID1, COVID2, post) %>% 
  count(time) -> datObIN

# june and july 2021 are very low because of opening new sub-hospital
# impute with average of may and aug 2021
datObIN$n <- with(datObIN, ifelse(year==2021 & month==6 | year==2021 & month==7,
                                  round((datObIN$n[datObIN$year==2021 & datObIN$month==5] +
                                           datObIN$n[datObIN$year==2021 & datObIN$month==8]) / 2, 0), 
                                  n))

# initial chart of # cases by month; just raw data (2017-2021)
ggplot(data = datObIN, aes(x = time, y = n)) +
  geom_point() +
  geom_line() +
  ggtitle("Obstetric inpatient admissions (Jan 2017-May 2021)") +
  labs(x = "Time", y = "Number of obstetric admissions") +
  scale_x_continuous(limits=c(24,76), breaks=c(24,36,48,60,72,76),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","May 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=66.5, y=2200, label="COVID-19",size=5) +
  expand_limits(y = 0) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObIN.17
pp.ObIN.17


### Analysis

# adjust for monthly seasonality, January as the control，no offset for newborn population becuase 2021 data are not accurate
mod.ObIN <- glm.nb(n ~ time1 + COVID1 + COVID2 + post 
                   + as.factor(month), 
                   data = subset(datObIN, year >= 2017))

lincom(mod.ObIN)

# Extract the fitted values from this fitted model
datObIN$fit <- NA
datObIN$fit[datObIN$year>=2017] <- fitted(mod.ObIN)

# Use AR(1)
coef.ObIN <- coeftest(mod.ObIN, vcov=NeweyWest(mod.ObIN,lag=1,prewhite = F))

# put coef, CI, and p-value for time1, COVID, post together
m.ObIN <- matrix(, nrow = 4, ncol = 4)
m.ObIN[,1]<-exp(coef.ObIN[c(2:5),1])
m.ObIN[,2]<-exp(coef.ObIN[c(2:5),1]-1.96*coef.ObIN[c(2:5),2])
m.ObIN[,3]<-exp(coef.ObIN[c(2:5),1]+1.96*coef.ObIN[c(2:5),2])
m.ObIN[,4]<-coef.ObIN[c(2:5),4]
rownames(m.ObIN) <- c("time1", "COVID1","COVID2", "post")
colnames(m.ObIN) <- c("exp(beta)", "95%low", "95%up", "p-value")
round(m.ObIN, 5)

# Counterfactual or expected Y values based on the original model
X0.ObIN <- X.ObIN <- model.matrix(mod.ObIN)
X0.ObIN[,3:5] <- 0
datObIN$exp[datObIN$year>=2017] <- as.vector(exp(X0.ObIN%*%coef.ObIN)) # save expected Y values to the monthly dataset

####### 10,000 simulations
Sigma.ObIN <- NeweyWest(mod.ObIN,lag=1,prewhite=F)
set.seed(12)
beta.ObIN <- mvrnorm(n = 10000, coef.ObIN[,1], Sigma.ObIN)

D.ObIN <- Y_sum.ObIN <- Y0_sum.ObIN <- D_P.ObIN <- RR.ObIN <- NA

for (i in 1:10000){
  Y_sum.ObIN[i] <- sum(as.vector(exp(X.ObIN%*%beta.ObIN[i,]))[38:53])
  Y0_sum.ObIN[i] <- sum(as.vector(exp(X0.ObIN%*%beta.ObIN[i,]))[38:53])
  D.ObIN[i] <- Y0_sum.ObIN[i] - Y_sum.ObIN[i]
  D_P.ObIN[i] <- (Y0_sum.ObIN[i] - Y_sum.ObIN[i]) / Y0_sum.ObIN[i]
  RR.ObIN[i] <- Y_sum.ObIN[i] / Y0_sum.ObIN[i]
}

# Fitted post
mean(Y_sum.ObIN)
quantile(Y_sum.ObIN, c(0.025,0.975)) #95% CI

# Counterfactual post
mean(Y0_sum.ObIN)
quantile(Y0_sum.ObIN, c(0.025,0.975)) 

# !!! Lost visits (observed - counterfactual/forecasted)
mean(D.ObIN) 
quantile(D.ObIN, c(0.025,0.975))
p.D.ObIN <- length(D.ObIN[D.ObIN < 0]) / length(D.ObIN) # p-value; if p=0, write it as p<0.0001
p.D.ObIN

mean(D_P.ObIN)  ## percentage of decrease
quantile(D_P.ObIN, c(0.025,0.975))
p.D_P.ObIN <- length(D_P.ObIN[D_P.ObIN < 0]) / length(D_P.ObIN) # p-value
p.D_P.ObIN

mean(RR.ObIN)  ## RR
quantile(RR.ObIN, c(0.025,0.975))

### Chart for the model analysis
# CI around the line, just for the fitted line.
datObIN$fit_up[datObIN$year>=2017] <- exp(predict(mod.ObIN) + 1.96*predict(mod.ObIN, se.fit=T)$se.fit)
datObIN$fit_low[datObIN$year>=2017] <- exp(predict(mod.ObIN) - 1.96*predict(mod.ObIN, se.fit=T)$se.fit)

# chart of observed versus expected
ggplot(data = datObIN, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp, colour = "Expected")) +
  geom_line(aes(y = fit, colour = "Fitted")) + 
  geom_point(aes(y = n, colour = "Fitted")) +
  ggtitle("Obstetric inpatient admissions: Observed vs. expected (Jan 2017 - Dec 2021)") +
  labs(x = "Time", y = "Obstetric admissions", color = "") +
  scale_x_continuous(limits=c(24,83), breaks=c(24,36,48,60,72,83),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Dec 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=67, y=1900, label="COVID-19",size=5) +
  expand_limits(y = c(1000:2000)) +
  goldentheme +
  theme(axis.text.x = element_text(size = 9.5)) -> pp.ObIN.ove
pp.ObIN.ove


#really lag 1?

library(tseries)

## outpatient
adf.test(datObOUT1$count) # larger than 0.05, dp not reject the null that the time series is not stationary. Since it is not stationary, it means we might not need to use autoregressive models.

# What's the lag?
acf(datObOUT1$count) #sort of decaying slowly, suitable for AR
pacf(datObOUT1$count) #seems lag = 2 is the last significant spike

## surgery
adf.test(datObSUR1$count) # larger than 0.05, dp not reject the null that the time series is not stationary. Since it is not stationary, it means we might not need to use autoregressive models.

# What's the lag?
acf(datObSUR1$count) #sort of decaying slowly, suitable for AR
pacf(datObSUR1$count) #seems lag = 2 is the last significant spike

## inpatient data
adf.test(datObIN$n) # larger than 0.05, dp not reject the null that the time series is not stationary. Since it is not stationary, it means we might not need to use autoregressive models.

# What's the lag?
acf(datObIN$n) #sort of decaying slowly, suitable for AR
pacf(datObIN$n) #seems lag = 1 is the last significant spike

# so lag really is 1



### Plot all 3 charts together

# outpatient
ggplot(data = datObOUT1, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp), color = "#0072B2", linetype = "dashed") +
  geom_line(aes(y = fit), color = "#0072B2",) + 
  geom_point(aes(y = count), color = "#0072B2", size = 1) +
  ggtitle("A. Outpatient visits") +
  labs(x = "", y = "", color = "") +
  scale_x_continuous(limits=c(24,83), breaks=c(24,36,48,60,72,83),
                     labels=addline_format(c("","","","","",""))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  expand_limits(y = c(8000:26000)) +
  scale_y_continuous(trans="log10") +
  goldentheme_facet +
  theme(axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold")) -> pp.ObOUT.ove1
pp.ObOUT.ove1

# surgeries
ggplot(data = datObSUR1, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp), color = "#000000", linetype = "dashed") +
  geom_line(aes(y = fit), color = "#000000") + 
  geom_point(aes(y = count), color = "#000000", size = 1) +
  ggtitle("B. Surgeries") +
  labs(x = "", y = "", color = "") +
  scale_x_continuous(limits=c(24,83), breaks=c(24,36,48,60,72,83),
                     labels=addline_format(c("","","","","",""))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  expand_limits(y = c(500:1200)) +
  scale_y_continuous(trans="log10") +
  goldentheme_facet +
  theme(axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold")) -> pp.ObSUR.ove1
pp.ObSUR.ove1

#inpatient
ggplot(data = datObIN, aes(x = time)) +
  geom_ribbon(aes(ymin = fit_low, ymax = fit_up), alpha = 0.15) +
  geom_line(aes(y = exp), color = "#D55E00", linetype = "dashed") +
  geom_line(aes(y = fit), color = "#D55E00") + 
  geom_point(aes(y = n), color = "#D55E00", size = 1) +
  ggtitle("C. Inpatient admissions") +
  labs(x = "Time", y = "", color = "") +
  scale_x_continuous(limits=c(24,83), breaks=c(24,36,48,60,72,83),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Dec 2021"))) +
  geom_vline(xintercept = 61, linetype = "dotted") +
  annotate(geom="text", x=67.5, y=1900, label="COVID-19",size=4) +
  expand_limits(y = c(1000:2000)) +
  scale_y_continuous(trans="log10") +
  goldentheme_facet +
  theme(axis.text.x = element_text(size = 10),
        axis.title.y=element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold")) -> pp.ObIN.ove1
pp.ObIN.ove1


threez <- plot_grid(pp.ObOUT.ove1, pp.ObSUR.ove1, pp.ObIN.ove1, nrow = 3,
                    rel_heights = c(1, 1, 1.15))

ggsave("~/Ob3in1.pdf", width = 5, height = 8)

threez

# all 3 lines in one chart
ggplot() +
  geom_vline(xintercept = 61, linetype = "dotted") +
  geom_ribbon(data = datObOUT1, aes(x = time, ymin = fit_low, ymax = fit_up), fill = "#0072B2", alpha = 0.15) +
  geom_ribbon(data = datObSUR1, aes(x = time, ymin = fit_low, ymax = fit_up), fill = "#000000", alpha = 0.15) +
  geom_ribbon(data = datObIN, aes(x = time, ymin = fit_low, ymax = fit_up), fill = "#D55E00", alpha = 0.15) +
  geom_line(data = datObOUT1, aes(x = time, y = exp), color = "#0072B2", linetype = "dashed") +
  geom_line(data = datObSUR1, aes(x = time, y = exp), color = "#000000", linetype = "dashed") +
  geom_line(data = datObIN, aes(x = time, y = exp), color = "#D55E00", linetype = "dashed") +
  geom_line(data = datObOUT1, aes(x = time, y = fit), color = "#0072B2",) + 
  geom_line(data = datObSUR1, aes(x = time, y = fit), color = "#000000") + 
  geom_line(data = datObIN, aes(x = time, y = fit), color = "#D55E00") + 
  geom_point(data = datObOUT1, aes(x = time, y = count), color = "#0072B2", size = 1) +
  geom_point(data = datObSUR1, aes(x = time, y = count), color = "#000000", size = 1) +
  geom_point(data = datObIN, aes(x = time, y = n), color = "#D55E00", size = 1) +
  ggtitle("") +
  labs(x = "Time", y = "", color = "") +
  scale_x_continuous(limits=c(24,83), breaks=c(24,36,48,60,72,83),
                     labels=addline_format(c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Dec 2021"))) +
  annotate(geom="text", x=67, y=3500, label="COVID-19",size=5) +
  annotate(geom="text", x=29, y=15000, label="Outpatient visits",size=4, color = "#0072B2") +
  annotate(geom="text", x=30.5, y=2400, label="Inpatient admissions",size=4, color = "#D55E00") +
  annotate(geom="text", x=27, y=800, label="Surgeries",size=4, color = "#000000") +
  expand_limits(y = c(100:30000)) +
  scale_y_continuous(trans="log10") +
  goldentheme +
  theme(axis.title.y=element_blank(),
        legend.position = "none",
        plot.title = element_blank()) -> pp.Ob.all.ove

ggsave("~/Ob3in1_v2.pdf", width = 8, height = 5)

pp.Ob.all.ove


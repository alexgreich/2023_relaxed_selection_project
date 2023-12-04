#####################
#alex sorts through all of the bullshit
#that is her masters thesis
#and organized results
#in preparation for publication and all that

#12/04/23 (yes, this has been ongoing for far too long)
###################

#NOTES to revisit
##p2 GSI has a large outlier, revisit that one (will impact results or not?)

#libraries
library(ggplot2)
library(cowplot)
library(lme4)
library(nlme)
library(tidyverse)
library(dplyr)

library(AICcmodavg)
library(MASS)
library(lattice)
library(GGally)
library(mgcv)
library(lmtest)
library(nlme)
library(visreg)

library(lmerTest)

#######################################
#####################################
#GSI
##a t-test.
#####################################
###########################################

#pink 2020
#load datas
p1 <- read.csv("Data/FemalePinkGSIotodataadd1_without_gaps copy.csv")

#prep GSI
p1$GSI <- (p1$Total.Egg.Mass.g./p1$Fish.Weight..g.) *100
p1.clean <- p1 %>% filter(Weird == "n", Otolith.results != "unknown") #getting rid of weird and unknown


p1.wild.GSI.ttest <- p1.clean %>% filter(Otolith.results == "no mark") #wild fish (they dont have an oto mark)
p1.hatch.GSI.ttest <- p1.clean %>% filter(Otolith.results == "PORT ARMSTRONG") #hatchery fish (they have PA's oto mark)

plot(p1.clean[c(2, 3, 4, 15, 11)]) #splorin

#testin
aov.p1.GSI <- aov(GSI~Otolith.results, p1.clean)
sum.p1.GSI <- summary.lm(aov(GSI~Otolith.results, p1.clean))
#(p1.GSI.t.test <- t.test(p1.hatch.GSI.ttest$GSI, p1.wild.GSI.ttest$GSI, alternative="greater")) 
(p1.GSI.t.test <- t.test(p1.hatch.GSI.ttest$GSI, p1.wild.GSI.ttest$GSI, alternative="greater", var.equal=T)) ##maybe I should use THIS one, if I assume the variances are equal

#equal variances?
plot(resid(aov.p1.GSI))

ggplot() + aes(y= resid(aov.p1.GSI), x=p1.clean$Otolith.results) + geom_point()
#identify(resid(aov.p1.GSI)) #point 35 is an outlier
#View(p1.clean)

plot(aov.p1.GSI) #says outlier is not sig

#sd 05/25/22
sd_p1_w <- sd(p1.wild.GSI.ttest$GSI)
sd_p1_h <- sd(p1.hatch.GSI.ttest$GSI) #excellent. Now do the GSI graphs

#length test 07/25/22 (dont use)
mod_p1_length <- lm(GSI~Otolith.results + Length..mm., p1.clean)
summary(mod_p1_length) #not sig diff lengths.


#################################################
#pink 2021 -length might be sig here
p2.GSI <- read.csv("Data/Female.p2.Rdata.3.csv")
names(p2.GSI)

p2.GSI <- p2.GSI %>%
  mutate(GSI.2 = (GSI.measure.2.g./Fish.weight.g.)*100, 
         GSI.1 = (GSI.measure.1.g./Fish.weight.g.)*100)  #we're going to use GSI.1 here though. It's the GSI before removing Ov fluid

p2.GSI.clean <- p2.GSI %>% filter(Weird =="n", Oto.reading != "No Oto", Oto.reading != "Overground")
p2.GSI.clean.relevant <- p2.GSI.clean[c(1, 2, 19,14,4)]

GSI.for.ttest.wild <- na.omit(p2.GSI.clean.relevant) %>% filter(Oto.reading == "No Mark")
GSI.for.ttest.hatch <- na.omit(p2.GSI.clean.relevant) %>% filter(Oto.reading == "PORT ARMSTRONG")

(p2.GSI.t.test <- t.test(GSI.for.ttest.hatch$GSI.1,GSI.for.ttest.wild$GSI.1, alternative="greater", var.equal=T))

#test resids
aov.p2.GSI <- aov(GSI.1 ~ Oto.reading, p2.GSI.clean.relevant)
plot(aov.p2.GSI) #51 looks like it may be an outlier, which is HF17
#View(p2.GSI.clean.relevant)
#View(p2.GSI.clean)
#NEED TO CHECK DATA ENTRY FOR HF17
#cooks says it is ok, but what about the unequal var?

aov.p2.GSI_testlength <- lm(GSI.1 ~ Oto.reading + Length.mm., p2.GSI.clean.relevant)
summary(aov.p2.GSI_testlength)

plot(p2.GSI.clean.relevant)


library(onewaytests)
?bf.test 
bf.test(GSI.1 ~ Oto.reading, p2.GSI.clean.relevant) #says variance is equal

##05/05/22
#test for date significance
#and maybe location significance
write.csv(p2.GSI.clean, "p2.GSI.clean.csv") #made this a csv. Now goitn in to edit the date, new date column with early, middle , late in the run
p2_GSI_clean_date_alt <- read.csv("Data/p2.GSI.clean_altdate.csv")

names(p2_GSI_clean_date_alt)
aov.p2.GSI_global <- lm(GSI.1 ~ Oto.reading + Date_alt, p2_GSI_clean_date_alt)
summary(aov.p2.GSI_global) #date not sig for GSI.

#sd 05/25/22
sd_w_p2 <- sd(GSI.for.ttest.wild$GSI.1)
sd_h_p2 <-sd(GSI.for.ttest.hatch$GSI.1)


#05/25/22
#LUBRIDATE
library(lubridate)
names(p2_GSI_clean_date_alt)
Date_date <- mdy(p2_GSI_clean_date_alt$Date)
class(Date_date)
Julian_GSI <- julian(Date_date, origin = as.Date("2021-01-01"))
p2_GSI_clean_date_alt$Julian <-  Julian_GSI

aov.p2.GSI_global <- lm(GSI.1 ~ Oto.reading + Julian, p2_GSI_clean_date_alt)
summary(aov.p2.GSI_global) #not sig at ALLL. Date is not signficant.

#REVISIT - this one has a werid outlier.
#10/21/22 Outlier check <- REVISIT
names(p2.GSI.clean.relevant)
max(p2.GSI.clean.relevant$GSI.1)  #26
range(p2.GSI.clean.relevant$GSI.1)
big_GSI <- p2.GSI.clean.relevant %>% filter(GSI.1>26)
## 1185 is fish weight
##fish id is HF17
##let's see if she has any weird notes
test_no_outlier <- p2.GSI.clean.relevant %>% filter(GSI.1<26)

GSI.for.ttest.wild_no_o <- na.omit(test_no_outlier) %>% filter(Oto.reading == "No Mark")
GSI.for.ttest.hatch_no_o <- na.omit(test_no_outlier) %>% filter(Oto.reading == "PORT ARMSTRONG")

(test_no_o_t <- t.test(GSI.for.ttest.hatch_no_o$GSI.1,GSI.for.ttest.wild_no_o$GSI.1, alternative="greater", var.equal=T))
#no different test results if we remove outlier.




#################################################################################################################33
#coho (full coho female GSI dataset)

c.GSI.dat <-read.csv("Data/MASTERFemaleCohoQCwitheggs_copy.csv")
c.GSI.dat$GSI <- (c.GSI.dat$Total.Egg.Mass.g./c.GSI.dat$Fish.Weight..g.)*100


#get GSI
c.GSI.clean <- na.omit(c.GSI.dat) %>% filter (Weird == "n")
c.GSI.wild.t.test <- c.GSI.clean %>% filter(Wild.or.Hatch=="wild")
c.GSI.hatch.t.test <- c.GSI.clean %>% filter(Wild.or.Hatch=="hatchery")

#write.csv(c.GSI.clean, "c.GSI.clean")

(coho.GSI.t.test <- t.test(c.GSI.hatch.t.test$GSI,c.GSI.wild.t.test$GSI, alternative="greater", var.equal=T))
#t = 1.4608, df = 51.395, p-value = 0.07507
#t = 1.4608, df = 52, p-value = 0.07504


hatchery_mean_value_GSI_fullcohodataset <- coho.GSI.t.test$estimate[1]
wild_mean_value_GSI_fullcohodataset <- coho.GSI.t.test$estimate[2]
coho_GSI_table_supplement <- data.frame(hatchery_mean_value_GSI_fullcohodataset, wild_mean_value_GSI_fullcohodataset)

sd_w_c <- sd(c.GSI.wild.t.test$GSI)
sd_h_c <- sd(c.GSI.hatch.t.test$GSI)


#var.equal?
c1.aov.GSI <- aov(GSI ~ Wild.or.Hatch, c.GSI.dat)
plot(c1.aov.GSI) #looks good


#####################################################
#results summary
p1.GSI.t.test #pink 2020
p2.GSI.t.test #pink 2021
coho.GSI.t.test #coho

GSI_results <- data.frame(name = c("pink 2020", "pink 2021", "coho"), 
                          t = c(p1.GSI.t.test$statistic, p2.GSI.t.test$statistic, coho.GSI.t.test$statistic),
                          df = c(p1.GSI.t.test$parameter, p2.GSI.t.test$parameter, coho.GSI.t.test$parameter),
                          p = c(p1.GSI.t.test$p.value, p2.GSI.t.test$p.value, coho.GSI.t.test$p.value),
                          hatchery_mean =c( p1.GSI.t.test$estimate[1], p2.GSI.t.test$estimate[1], hatchery_mean_value_GSI_fullcohodataset),
                          hatch_sd = c( sd_p1_h, sd_h_p2 ,sd_h_c),
                          wild_mean = c(p1.GSI.t.test$estimate[2], p2.GSI.t.test$estimate[2], wild_mean_value_GSI_fullcohodataset),
                          wild_sd = c( sd_p1_w, sd_w_p2 ,sd_w_c)
                          )



#length results summary

#####################################
#egg diameter
######################################

#pink 2020

#pink 2021

#coho

#results summary
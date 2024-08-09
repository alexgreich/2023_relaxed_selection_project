#####################
#organizing GSI/egg results for publication
#prupose: this turned into female code and results

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
library(patchwork)

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

#sample sizes
ss_p1_h <- length(p1.hatch.GSI.ttest$GSI)
ss_p1_w <- length(p1.wild.GSI.ttest$GSI)

#length test 07/25/22 (dont use)
mod_p1_length <- lm(GSI~Otolith.results + Length..mm., p1.clean)
summary(mod_p1_length) #not sig diff lengths.

#length t test for even yr female pinks
p1_f_L_ttest <- t.test(p1.hatch.GSI.ttest$Length..mm., p1.wild.GSI.ttest$Length..mm., var.equal=T)

mean_length_p1_h <- mean(p1.hatch.GSI.ttest$Length..mm.)
mean_length_p1_w <- mean(p1.wild.GSI.ttest$Length..mm.)
sd_length_p1_h <- sd(p1.hatch.GSI.ttest$Length..mm.)
sd_length_p1_w <- sd(p1.wild.GSI.ttest$Length..mm.)

#GSI but length is a covariate (added 1/22/24)
p1_GSI_linearmod <- lm(GSI~factor(Otolith.results) + Length..mm. , p1.clean)

p1_GSI_lin_sum <- summary(p1_GSI_linearmod) #will need to get the 1-tailed test for this later


#gsi hypothesis: hatch > wild
#so needs to be 1-p val to get the right tail, divided by 2 because one sided
p1_p_GSI_linmod <- 1-(p1_GSI_lin_sum$coefficients[2,4]/2) #having a fuck, how do I make this a one-sided test moment. It makes sense that the p-value doesnt change if the order of the factors change. Only the side (pos/neg) of t changes. How do I make it one-sided
p1_t_GSI_linmod <- p1_GSI_lin_sum$coefficients[2,3]
p1_df_GSI_linmod <-p1_GSI_lin_sum$df[2]

#AIC- added 5/20/24 - adj on 6/6/24 for the LRT values
##my models are: 
###
p1.clean$Otolith.results <- factor(p1.clean$Otolith.results)

p1_GSI_linearmod
p1_GSI_mod2 <- lm(GSI~Otolith.results, p1.clean)
p1_GSI_int <- lm(GSI~Otolith.results + Length..mm. +  Otolith.results:Length..mm., p1.clean)
p1_GSI_null <- lm(GSI~1, p1.clean) #NEW 6/6/24

AIC(p1_GSI_linearmod, p1_GSI_mod2, p1_GSI_int,p1_GSI_null ) #null mod wins. Nothing explains anything
summary(p1_GSI_mod2)
summary(p1_GSI_null)
library(lmtest)
?lrtest
lrtest(p1_GSI_null, p1_GSI_mod2) #small (null) mod wins, is what this is telling me
lrtest(p1_GSI_mod2, p1_GSI_linearmod) #new mod... is not better?
lrtest(p1_GSI_linearmod, p1_GSI_int ) #fuck, what does this mean?

anova(p1_GSI_linearmod, p1_GSI_int)



#0.84147/0.78953
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

#sample sizes
ss_p2_h <- length(GSI.for.ttest.hatch$GSI.1)
ss_p2_w <- length(GSI.for.ttest.wild$GSI.1)

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

#length t test:
p2_f_L_ttest <- t.test(GSI.for.ttest.hatch$Length.mm., GSI.for.ttest.wild$Length.mm., var.equal = T)

mean_length_p2_h <- mean(GSI.for.ttest.hatch$Length.mm.)
mean_length_p2_w <- mean(GSI.for.ttest.wild$Length.mm.)
sd_length_p2_h <- sd(GSI.for.ttest.hatch$Length.mm.)
sd_length_p2_w <- sd(GSI.for.ttest.wild$Length.mm.)




#GSI but length is a covariate (added 1/22/24)
p2_GSI_linearmod <- lm(GSI.1~factor(Oto.reading) + Length.mm. , p2.GSI.clean.relevant)

p2_GSI_lin_sum <- summary(p2_GSI_linearmod) #will need to get the 1-tailed test for this later


#gsi hypothesis: hatch > wild
#so needs to be 1-p val to get the right tail, divided by 2 because one sided
p2_p_GSI_linmod <- 1-(p2_GSI_lin_sum$coefficients[2,4]/2) 
p2_t_GSI_linmod <- p2_GSI_lin_sum$coefficients[2,3]
p2_df_GSI_linmod <-p2_GSI_lin_sum$df[2]

##added: 5/20/24: AIC tables
p2_GSI_linearmod_jdate <- lm(GSI.1~factor(Oto.reading) + Length.mm.+ Julian , p2_GSI_clean_date_alt)
p2_GSI_linearmod <- lm(GSI.1~factor(Oto.reading) + Length.mm., p2_GSI_clean_date_alt)
p2_GSI_linearmod_2 <- lm(GSI.1~factor(Oto.reading), p2_GSI_clean_date_alt)
p2_int_1 <- lm(GSI.1~factor(Oto.reading) * Length.mm.* Julian , p2_GSI_clean_date_alt)
p2_int_2 <- lm(GSI.1~factor(Oto.reading) + Length.mm.+ Julian + factor(Oto.reading):Length.mm.+ Length.mm.:Julian+ Julian:factor(Oto.reading) , p2_GSI_clean_date_alt)
p2_int_3 <- lm(GSI.1~factor(Oto.reading) + Length.mm.+ Julian + Length.mm.:Julian+ Julian:factor(Oto.reading), p2_GSI_clean_date_alt)
p2_int_3.5 <- lm(GSI.1~factor(Oto.reading) + Length.mm.+ Julian + Julian:factor(Oto.reading), p2_GSI_clean_date_alt)
p2_int_4 <- lm(GSI.1~factor(Oto.reading) + Julian + Julian:factor(Oto.reading), p2_GSI_clean_date_alt)
p2_int_5 <- lm(GSI.1~factor(Oto.reading) + Julian, p2_GSI_clean_date_alt)
p2_int_6 <- lm(GSI.1~factor(Oto.reading), p2_GSI_clean_date_alt)


summary(p2_int_2)
AIC(p2_GSI_linearmod, p2_GSI_linearmod_jdate, p2_GSI_linearmod_2, p2_int_1 ) #no length, no date wins!

AIC(p2_GSI_linearmod, p2_GSI_linearmod_jdate, p2_GSI_linearmod_2, p2_int_1,
    p2_int_2, p2_int_3, p2_int_3.5, p2_int_4, p2_int_5, p2_int_6
    ) 

#p2_GSI_linearmod
#p2_GSI_linearmod_jdate <- lm(GSI.1~factor(Oto.reading) + Length.mm.+ Julian , p2.GSI.clean.relevant)
#p2_GSI_linearmod_2 <- lm(GSI.1~factor(Oto.reading), p2.GSI.clean.relevant)
#AIC(p2_GSI_linearmod, p2_GSI_linearmod_jdate, p2_GSI_linearmod_2)


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

#test if length is sig for GSI
coho_mod <- lm(GSI ~ Length..mm. + Wild.or.Hatch, c.GSI.clean)
summary(coho_mod) #hmmm. length is barely non sig. Larger fish should have a larger GSI, but smaller fish have a larger GSI, for both wild and hatch
coho_mod_test <- lm(GSI ~ Wild.or.Hatch, c.GSI.clean)
AIC(coho_mod, coho_mod_test)

hatchery_mean_value_GSI_fullcohodataset <- coho.GSI.t.test$estimate[1]
wild_mean_value_GSI_fullcohodataset <- coho.GSI.t.test$estimate[2]
coho_GSI_table_supplement <- data.frame(hatchery_mean_value_GSI_fullcohodataset, wild_mean_value_GSI_fullcohodataset)

sd_w_c <- sd(c.GSI.wild.t.test$GSI)
sd_h_c <- sd(c.GSI.hatch.t.test$GSI)


#var.equal?
c1.aov.GSI <- aov(GSI ~ Wild.or.Hatch, c.GSI.dat)
plot(c1.aov.GSI) #looks good

#length t test is down by coho egg




#GSI but length is a covariate (added 1/22/24)
c_GSI_linearmod <- lm(GSI~factor(Wild.or.Hatch, levels=c("wild", "hatchery") ) + Length..mm. , c.GSI.clean)

c_GSI_lin_sum <- summary(c_GSI_linearmod) #will need to get the 1-tailed test for this later


#gsi hypothesis: hatch > wild
#so needs to be 1-p val to get the right tail, divided by 2 because one sided
c_p_GSI_linmod <- (c_GSI_lin_sum$coefficients[2,4]/2) 
c_t_GSI_linmod <- c_GSI_lin_sum$coefficients[2,3]
c_df_GSI_linmod <-c_GSI_lin_sum$df[2]

#AIC for tables- added 5/20/24
##model: c_GSI_linearmod
c_GSI_linearmod_2 <- lm(GSI~factor(Wild.or.Hatch, levels=c("wild", "hatchery") ) , c.GSI.clean)
c_GSI_linearmod_int <- lm(GSI~factor(Wild.or.Hatch, levels=c("wild", "hatchery")) + Length..mm. + factor(Wild.or.Hatch, levels=c("wild", "hatchery")):Length..mm. , c.GSI.clean)
AIC(c_GSI_linearmod ,c_GSI_linearmod_2, c_GSI_linearmod_int) #mod 4 wins but it is super borderline
##so need to get 


#####################################################
#results summary
p1.GSI.t.test #pink 2020
p2.GSI.t.test #pink 2021
coho.GSI.t.test #coho

GSI_results <- data.frame(#name = c("pink 2020", "pink 2021", "coho"), #code not working here, don't know why. R is being squirrely
                          t = c(p1.GSI.t.test$statistic, p2.GSI.t.test$statistic, coho.GSI.t.test$statistic),
                          df = c(p1.GSI.t.test$parameter, p2.GSI.t.test$parameter, coho.GSI.t.test$parameter),
                          p = c(p1.GSI.t.test$p.value, p2.GSI.t.test$p.value, coho.GSI.t.test$p.value),
                          p_lin_mod = c(p1_p_GSI_linmod,p2_p_GSI_linmod,c_p_GSI_linmod), #from when length is included in the models
                          hatchery_mean =c( p1.GSI.t.test$estimate[1], p2.GSI.t.test$estimate[1], hatchery_mean_value_GSI_fullcohodataset),
                          hatch_sd = c( sd_p1_h, sd_h_p2 ,sd_h_c),
                          wild_mean = c(p1.GSI.t.test$estimate[2], p2.GSI.t.test$estimate[2], wild_mean_value_GSI_fullcohodataset),
                          wild_sd = c( sd_p1_w, sd_w_p2 ,sd_w_c)
                          )

rownames(GSI_results) = c("pink 2020", "pink 2021", "coho")

p_lin_mod <- c(p1_p_GSI_linmod,p2_p_GSI_linmod,c_p_GSI_linmod)
t_lin_mod <- c(p1_t_GSI_linmod,p2_t_GSI_linmod,c_t_GSI_linmod)
df_lin_mod <- c(p1_df_GSI_linmod, p2_df_GSI_linmod, c_df_GSI_linmod )

datafram_lin_mod <- data.frame(p_lin_mod, t_lin_mod, df_lin_mod)
rownames(datafram_lin_mod) =c("pink 2020","pink 2021","coho")
colnames(datafram_lin_mod)=c("p", "t", "df")


#write.csv(x=datafram_lin_mod, file="Results/If we do GSI with length in the linear model.csv")


#5/20/24
##comittee had me re-do things with AIC. resutls from p1, p2 mods without length and the results from coho with length. Sigh.
##just took GSI_results p1 and p2 values and datafram_lin_mod c values

##
#08/08/24 - GSI part for the supplemental table from hell- p1 summary GSI
##so resukts from p1 and p2 are without length and coho is with length? Ok then

#p1
#GSI
#summary(p1_GSI_linearmod)
summary(p1_GSI_int) #full model p1 GSI
summary(p1_GSI_mod2) #selected model p1 GSI
#summary(p1_GSI_null)
#one way test t test:
p1.GSI.t.test 


#p2
p2_glob <- lm(GSI.1~factor(Oto.reading) * Length.mm. , p2.GSI.clean.relevant)
p2_sel <- lm(GSI.1~factor(Oto.reading) , p2.GSI.clean.relevant)

summary(p2_glob) #global p2 mod
summary(p2_sel) #selected p2 mod

p2.GSI.t.test #the one-way test for the p2 mod, (the hypothesis test)

#c
#summary(c_GSI_linearmod_2)
summary(c_GSI_linearmod_int) #global model
summary(c_GSI_linearmod) #selected model

c_p_GSI_linmod

##

###########################################################3
############################################################
#egg diameter
###################################################################
###################################################################
#pink 2020
p1.all.eggs <- read.csv("Data/P1.egg.diameters.csv")
names(p1.all.eggs)
head(p1.all.eggs)

#the other pink even data
p1.other.data<- read.csv("Data/FemalePinkGSIotodataadd1_without_gaps copy.csv")
p1.other.data$Fish.ID <- p1.other.data$ID

#extract the oto result...or Join by the ID!!!!
p1.df <- full_join(p1.all.eggs, p1.other.data, by="Fish.ID")
p1.df.no.unknown <-p1.df %>% filter(Otolith.results != "unknown") #hmm no weird column might have to adjust, make sure WEIRD fish are out.
p1.df.clean <- p1.df.no.unknown %>% filter(Weird =="n")


#conceptual model
fit.lm.p1 <- lm(Diameter..mm. ~ Length..mm. + Otolith.results + Length..mm.:Otolith.results, data=na.omit(p1.df.clean))

# use REML to find the best random structure
fit.p1.B <- lmer(Diameter..mm. ~ Length..mm. + Otolith.results +
                   Length..mm.:Otolith.results +
                   (1|ID), data=na.omit(p1.df.clean), REML=T)

AIC(fit.p1.B, fit.lm.p1) #random intercept wins, by a lot

#use ML to find best fitted structure
fit.p1.B1 <- lmer(Diameter..mm. ~ Length..mm. + Otolith.results +
                    Length..mm.:Otolith.results +
                    (1|ID), data=na.omit(p1.df.clean), REML=F)
fit.p1.C <- lmer(Diameter..mm. ~ Length..mm. + Otolith.results +(1|ID), data=na.omit(p1.df.clean), REML=F)
fit.p1.C1 <- lmer(Diameter..mm. ~ Otolith.results +(1|ID), data=na.omit(p1.df.clean), REML=F)
fit.p1.C2 <- lmer(Diameter..mm. ~ 1 + (1|ID), data=na.omit(p1.df.clean), REML=F)
AIC(fit.p1.B1,fit.p1.C, fit.p1.C1, fit.p1.C2) #simplest one wins

#best model for p1: fit.p1.D
fit.p1.D <- lmer(Diameter..mm. ~ 1 + (1|ID), data=na.omit(p1.df.clean), REML=T)
fit.p1.C1.REML <- lmer(Diameter..mm. ~ Otolith.results +(1|ID), data=na.omit(p1.df.clean), REML=T)
fit.p1.lme <-  lme(fixed = Diameter..mm. ~ 1, random = ~1|ID, data=p1.df.clean, method="REML")
summary(fit.p1.D)
summary(fit.p1.C1.REML) 

#t-test sig
##yeilds 0.162
(0.162)/2
#t=1.424 on 44 df, p= 0.081


#get those parameter estimates...
fit.p1.lme.param.noint <-  lme(fixed = Diameter..mm. ~ Otolith.results -1 , random = ~1|ID, data=p1.df.clean, method="REML")
fixef(fit.p1.lme.param.noint)
summary(fit.p1.lme.param.noint)

fit.p1.lme.param.yesint <-  lme(fixed = Diameter..mm. ~ Otolith.results , random = ~1|ID, data=p1.df.clean, method="REML")
summary(fit.p1.lme.param.yesint)


#wait but which way is that t-test?
0.162/2
1-0.162/2


#a graph
p1.df.clean$fitted <- fitted(fit.p1.D) 
(p1.goodplot <- ggplot(p1.df.clean) + aes() + geom_jitter(aes(y=Diameter..mm., x=Otolith.results)) + 
    geom_boxplot(aes(y=fitted, x=Otolith.results), alpha=0.7))

summary(fit.p1.D)

#05/23/24
##milo says add AIC tables
#glob_p1 <- 


#####################################################################################
#pink 2021
#load in them datas
p2.GSI <- read.csv("Data/Female.p2.Rdata.3.csv")
names(p2.GSI)

library(dplyr)
GSI.data <- p2.GSI %>%
  mutate(GSI.2 = (GSI.measure.2.g./Fish.weight.g.)*100, 
         GSI.1 = (GSI.measure.1.g./Fish.weight.g.)*100)

## y=GSI, x= length, color and shape= otolith resuls
GSI.data.2 <- GSI.data %>%
  mutate(Oto.reading.2 = ifelse(Oto.reading!="Overground", Oto.reading, "none"),)
GSI.data.3 <- GSI.data.2 %>%
  mutate(Oto.reading.3 = ifelse(Oto.reading.2!="No Oto", Oto.reading.2, "none"),)
GSI.data.4 <- GSI.data.3 %>%
  mutate(Oto.reading.4 = ifelse(is.na(Oto.reading.3), "none", Oto.reading.3),)

GSI.nogreen <- GSI.data.4 %>% filter(Weird == "n")

#ok , that was just GSI, don't actually need the above...
#EGG TIME!


#pink odd (p2)###########################
p2.all.eggs <- read.csv("Data/Results.p2.ImageJ.eggs.csv")
library(tidyverse)
df.hold <- data.frame(unique(p2.all.eggs$Label))
#I manually added in the photo ID's. They are under "Label" in both of the dataframes
length(na.omit(GSI.data.4$Label)) #and then QC'ed.because they didnt line up. Now they line up
length(unique(p2.all.eggs$Label))

order(unique(p2.all.eggs$Label))

##MERGE TIME
p2.df <-  full_join(p2.all.eggs, GSI.data.4, by="Label")
head(p2.df) 
names(p2.df)
p2.df.2 <- p2.df %>%   #get the diameter from the area
  mutate(Diameter = 2*sqrt(Area/pi))
head(p2.df.2)
names(p2.df.2)

#let's remove the unknowns
p2.df.no.unknown <- p2.df.2 %>% filter(Oto.reading.4 != "none")
#and remove weird (the green or spawned fish)
p2.df.clean <- p2.df.no.unknown %>% filter(Weird == "n")
#View(p2.df.clean)

##IN THE MIDDLE OF THIS NOW. TO FULLY INCLUDE DATE IN CONSIDERATION? I THINK YES....
##addition 05/24/22: test julian date
#View(p2.df.clean)
class(p2.df.clean$Date)
library(lubridate)
?lubridate #let's make the dates Julian dates, a numeric
Date_adj <- mdy(p2.df.clean$Date)
class(Date_adj) #seemed to work. Tranformed to date
Julian <- julian(Date_adj, origin = as.Date("2021-01-01"))
class(Julian) #neat! I did it!
p2.df.clean$Julian <- Julian

##ZURR steps
#1 select a full model that includes "important" fixed effects
fit.p2.A <- lm(Diameter ~ Length.mm. + Oto.reading.4 + Julian +  Length.mm.:Oto.reading.4, data=p2.df.clean)
summary(fit.p2.A)


#2 select/search for optimal random effect structure using REML
library(lme4)
fit.p2.B <- lmer(Diameter ~ Oto.reading.4 + Length.mm. + Julian + Oto.reading.4:Length.mm. + (1|ID), data=p2.df.clean, REML=T)
summary(fit.p2.B)

AIC(fit.p2.A, fit.p2.B) #B wins (with random effects) 

#3 Using optimal random-effect structure. select fixed components for best model using ML #COMPARE USING ML #QCING HERE CURRENTLY
fit.p2.B.ml <- lmer(Diameter ~ Oto.reading.4 + Length.mm. + Julian + Oto.reading.4:Length.mm. +(1|ID), data=p2.df.clean, REML=F)
fit.p2.Julian <- lmer(Diameter ~ Oto.reading.4 + Length.mm.+ Julian +(1|ID), data=p2.df.clean, REML=F)
fit.p2.Julian2 <- lmer(Diameter ~ Oto.reading.4 + Julian +(1|ID), data=p2.df.clean, REML=F)
fit.p2.B3.ml <- lmer(Diameter ~ Oto.reading.4 + Length.mm. +(1|ID), data=p2.df.clean, REML=F)
fit.p2.alpha <- lmer(Diameter ~ Oto.reading.4 + Length.mm. + Date + (1|ID), data=p2.df.clean, REML=F)
fit.p2.B4.ml <-lmer(Diameter ~ Oto.reading.4 + (1|ID), data=p2.df.clean, REML=F) #this one
fit.p2.B5.ml <-lmer(Diameter ~ 1 + (1|ID), data=p2.df.clean, REML=F) 

AIC(fit.p2.B.ml, fit.p2.B3.ml, fit.p2.B4.ml, fit.p2.alpha, fit.p2.B5.ml, fit.p2.Julian, fit.p2.Julian2) #fit.p2.B4.ml wins
#update 05/24/22. Fuck, Julian2 date is my best mod. Wait... nno it's not. B4 is better, actually, because it's simpler and within 1 AIC point. So no Julian date needed in this analysis. #NO JULIAN DATE NEEDED. JULIAN DATE DISPROVED. Seemed to have a small effect, but not below the signficance threshhold AND AIC didnt identify it as importatn. WOOO.
BIC(fit.p2.B.ml, fit.p2.B3.ml, fit.p2.B4.ml, fit.p2.alpha, fit.p2.B5.ml, fit.p2.Julian)


##TANGENT - DISPROVED ABOVE!
#summary(fit.p2.alpha)
#anova(fit.p2.alpha, fit.p2.B3.ml) #fuck, the date model is sig. vua the loglik but NOT via AIC... #for simplicity, I'll stick with NO DATE for now...
##TANGENT

#use the REML verstion of the best model to get param estimates
p2.relevant <- na.omit(p2.df.clean[c(14, 34, 35)])

fit.p2.B4.reml <-lmer(Diameter ~ Oto.reading.4 + (1|ID), data=p2.relevant, REML=T) #this one
summary(fit.p2.B4.reml) #attach lmerTest for one thing, detatch lmerTest for other #with lmerTest, gives t-value but no df or p-value. Wihtout lmer test, yes p-val
#what was the thing about REML not relevant for comparison... REML for parameter estimates, ML for model selection?
compare <- lmer(Diameter ~ (1|ID), data=p2.relevant, REML=T) #this one
anova(fit.p2.B4.reml, compare )
#Anova(fit.p2.B4.reml, compare )

fit.p2.lme <- lme(fixed = Diameter ~ Oto.reading.4, random = ~1|ID, data= p2.relevant, method="REML")
fit.p2.lme.noint <- lme(fixed = Diameter ~ Oto.reading.4-1, random = ~1|ID, data= p2.relevant, method="REML")

summary(fit.p2.B4.reml)
fixef(fit.p2.lme.noint)
summary(fit.p2.lme)
summary(fit.p2.lme.noint)

#p-val: 0.000252
1-(0.000252/2)
#(1-0.000252)/2

#


############################################################################################
#coho
#prepping data. And get rid of them NA's!
#setwd("/Users/alexandrareich/Desktop/THESIS TIME!/Compilation of GITHUB project code/Relaxed_selection_proj_reproducible_science/DATA") 
###what the hell is this (above)? some relic code fro when working directories were a grand mystery? - 12/04/23
c.all.eggs <- read.csv("Data/coho.egg.diameters.csv")
c.all.other.stuff <- read.csv("Data/MASTERFemaleCohoQCwitheggs_copy.csv")

names(c.all.eggs)
names(c.all.other.stuff)


c.df <- full_join(c.all.eggs, c.all.other.stuff, by="Fish.ID")
head(c.df)
#View(c.df)
names(c.df)

coho.clean <- c.df %>% filter(Weird == "n")
#write.csv(coho.clean, file="coho_clean_for_641_project")

#aside: lenght t-test
c_for_length <- c.all.other.stuff %>% filter(Weird=="n")
length(c_for_length$Length..mm.)
c_for_length_h <- c_for_length %>% filter(Wild.or.Hatch=="hatchery") %>%
  dplyr::select(Length..mm.)
c_for_length_w <- c_for_length %>% filter(Wild.or.Hatch=="wild") %>%
  dplyr::select(Length..mm.)

c_f_L_ttest <- t.test(c_for_length_h, c_for_length_w, var.equal=T)

mean_length_c_h <- mean(c_for_length_h$Length..mm.)
mean_length_c_w <- mean(c_for_length_w$Length..mm.)
sd_length_c_h <- sd(c_for_length_h$Length..mm.)
sd_length_c_w <- sd(c_for_length_w$Length..mm.)

#conceptual model (mixed model with interaction effects too)
fit.c.A.INT <- lm(Diameter..mm. ~ Length..mm. + Wild.or.Hatch + Length..mm.:Wild.or.Hatch, data=coho.clean)

#reml to select mixed effects
fit.c.GLOBAL <- lmer(Diameter..mm. ~ Length..mm. + Wild.or.Hatch+Length..mm.:Wild.or.Hatch +(1|ID), data=coho.clean, REML=T)
#fit.c.B.INT <- lmer(Diameter..mm. ~ Length..mm. + Wild.or.Hatch + Length..mm.:Wild.or.Hatch +(1|ID), data=coho.clean, REML=T)
AIC(fit.c.GLOBAL, fit.c.A.INT) #mixed effect model wins!

#ml to select fixed effects
fit.c.GLOBAL.ml <- lmer(Diameter..mm. ~ Length..mm. + Wild.or.Hatch+Length..mm.:Wild.or.Hatch +(1|ID), data=coho.clean, REML=F)

fit.c.C <- lmer(Diameter..mm. ~ Length..mm. + Wild.or.Hatch +(1|ID), data=coho.clean, REML = F)

fit.c.C1 <- lmer(Diameter..mm. ~ Wild.or.Hatch +(1|ID), data=coho.clean, REML = F)
fit.c.C2 <- lmer(Diameter..mm. ~ 1 + (1|ID), data=coho.clean, REML = F)
fit.c.C3 <- lmer(Diameter..mm. ~ Length..mm. +(1|ID), data=coho.clean, REML = F)
AIC(fit.c.C, fit.c.C1, fit.c.C2, fit.c.C3, fit.c.GLOBAL.ml) #C1 is winner.




#model that was selected
##(and fit to model without intercept for the standard errors)
fit.c.C1.REML <- lmer(Diameter..mm. ~ Wild.or.Hatch +(1|ID), data=coho.clean, REML = T)
summary(fit.c.C1.REML)

fit.lme <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch, random = ~1|ID, data=coho.clean, method="REML")
summary(fit.lme)

fit.lme.no.int <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch-1, random = ~1|ID, data=coho.clean, method="REML")
summary(fit.lme.no.int)


#parameter estimates
ranef(fit.lme.no.int)
fixef(fit.lme.no.int)
summary(fit.lme)
summary(fit.lme.no.int)

#graph that
(coho.egg.graph <- ggplot(coho.clean) + aes(x=Length..mm., y=Diameter..mm., color=Wild.or.Hatch) +
    geom_point(alpha=0.5) + scale_color_manual(values=c("blue", "orange"), name=element_blank(), labels=c("Wild origin", "Hatchery origin"), breaks=c("wild", "hatchery"))+
    labs(y="Egg diameter (mm)", x="MEHP length (mm)")+
    theme_cowplot() +
    theme(legend.position = c(0.03, 0.94))+
    coord_cartesian(ylim=c(5,8.5)) +
    scale_y_continuous(breaks=c(5,6,7,8), expand=c(0,0)))

(cf.fix <- fixef(fit.lme.no.int))
# and the corresponding random effects:
(cf.rand <- ranef(fit.lme.no.int))
Intervals <- intervals(fit.lme.no.int)
names(Intervals)
Intervals$fixed[1,1]
Intervals$fixed[1,3]


coho.egg.graph + 
  geom_ribbon(aes(ymin=Intervals$fixed[1,1], ymax=Intervals$fixed[1,3]), color="grey", alpha=0.3) + geom_ribbon(aes(ymin= Intervals$fixed[2,1], ymax=(Intervals$fixed[2,3])), color="grey", alpha=0.3)+
  geom_hline(yintercept = cf.fix[1], color="orange", size=1.5) + geom_hline(yintercept = cf.fix[2], color="blue", size=1.5)+
  coord_cartesian(ylim=c(5,8.3)) +
  scale_y_continuous(breaks=c(5,6,7,8), expand=c(0,0)) 


#QCing fixef values.
##what am I deleting when I na.omit(coho.clean) ???
names(coho.clean)
coho.clean.relevant <- coho.clean[c(1, 4, 3, 5, 13)]
names(coho.clean.relevant)


fit.lme.no.int.testna <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch-1, random = ~1|ID, data=coho.clean.relevant, method="REML")
fixef(fit.lme.no.int.testna)

fit.lme.no.int.testna2 <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch-1, random = ~1|ID, data=na.omit(coho.clean.relevant), method="REML")
fixef(fit.lme.no.int.testna2)

fit.lme.no.int.testna3 <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch-1, random = ~1|ID, data=na.omit(coho.clean), method="REML")
fixef(fit.lme.no.int.testna3)
#based on these tests, the original way of omiting all NA cut out things I didnt want to get rid of!
unique(coho.clean$Weird)



##AAAAND t-test that coho datas.
t.test.wild <- coho.clean %>% filter(Wild.or.Hatch=="wild")
t.test.hatch <- coho.clean %>% filter(Wild.or.Hatch=="hatchery")
#(t.test.coho.eggs <-t.test(t.test.hatch$Diameter..mm.,t.test.wild$Diameter..mm., alternative="less", var.equal = T))

summary(fit.lme) #this does not match
summary(fit.lme)
summary(fit.c.C1.REML)

#p-val for fitted model
fit.null.ML <- lme(fixed = Diameter..mm. ~ 1, random = ~1|ID, data=coho.clean, method="ML")
fit.lme.ML <- lme(fixed = Diameter..mm. ~ Wild.or.Hatch, random = ~1|ID, data=coho.clean, method="ML")
summary(fit.lme.ML)
anova(fit.lme, test=T, type="marginal")
anova(fit.lme.ML, fit.null.ML)


####################
#side function
#confient is: mean +- 1.98*se
sd_from_confint <- function(upper, lower){
  u <- mean(c(upper, lower))
  sd_calced <- (upper-u)/1.98 #ned to do this formula on paper. Ok I did. This looks right
}

###############################################################################################
# Eggs results summary
GSI_results #use for example

#the summary results
##p1
summary(fit.p1.C1.REML)
sum_reml_p1 <- summary(fit.p1.C1.REML)
fixef(fit.p1.lme.param.noint)
(sum_lme_noint_p1<-summary(fit.p1.lme.param.noint))
p_egg_p1 <- sum_reml_p1$coefficients[2,5]/2 #divided by 2 becsue the hypothesis is one-sided

##p2
(sum_reml_p2<-summary(fit.p2.B4.reml))
fixef(fit.p2.B4.reml)
fixef(fit.p2.lme.noint)
summary(fit.p2.lme)
(sum_lme_noint_p2<-summary(fit.p2.lme.noint))
names(summary(fit.p2.B4.reml))
confint(fit.p2.B4.reml)#seemed to have worked??
nlme::intervals(fit.p2.lme.noint)
p_egg_p2 <- 1-(sum_reml_p2$coefficients[2,5]/2) #1- because opposite of hypothesized. /2 becasue hypothesis is

##c
ranef(fit.lme.no.int)
fixef(fit.lme.no.int)
(sum_lme_noint_c<-summary(fit.lme.no.int))
(sum_reml_c <- summary(fit.c.C1.REML))
summary(fit.lme)
anova(fit.lme, test=T, type="marginal")
anova(fit.c.C1.REML, test=T, type="marginal")
p_egg_c <- sum_reml_c$coefficients[2,5]/2


#making a results table for export
Egg_results <- data.frame(
  #name=c("pink 2020", "pink 2021", "coho"),
  t=c(sum_reml_p1$coefficients[2,4], sum_reml_p2$coefficients[2,4], sum_reml_c$coefficients[2,4]),
  p=c( p_egg_p1, p_egg_p2, p_egg_c),
  df = c(sum_reml_p1$coefficients[2,3], sum_reml_p2$coefficients[2,3],sum_reml_c$coefficients[2,3]),
  hatchery_mean =c(fixef(fit.p1.lme.param.noint)[2],fixef(fit.p2.lme.noint)[2],fixef(fit.lme.no.int)[1]),
  hatch_sd=c (coef(sum_lme_noint_p1)[2,2], coef(sum_lme_noint_p2)[2,2],coef(sum_lme_noint_c)[1,2]),
  wild_mean=c(fixef(fit.p1.lme.param.noint)[1],fixef(fit.p2.lme.noint)[1],fixef(fit.lme.no.int)[2]),
  wild_sd=c(coef(sum_lme_noint_p1)[1,2], coef(sum_lme_noint_p2)[1,2], coef(sum_lme_noint_c)[2,2]),
  se_fixed=c(sum_reml_p1$sigma,sum_reml_p2$sigma, sum_reml_c$sigma), #within a fish
  se_random=c(0.233099,0.2411312,0.3157946)#, #among fish. Copy-pasted from the model summaries because it was taking too long to figure out how to input direct.
  #sample_size=c()
)

rownames(Egg_results) <- c("pink 2020", "pink 2021", "coho")


#results summary so far:
GSI_results
Egg_results

#spit out my final results in a csv
##write.csv(GSI_results, "Results/GSI_results.csv")
##write.csv(Egg_results, "Results/Egg_results.csv")

#08/09/24 supp table from hell
##summaries for eggs:

 #p1
fit.p1.B #glob model
sum_reml_p1  # selected model
p_egg_p1 # one-sided t test

 #p2
fit_p2_glob <- lmer(Diameter ~ Oto.reading.4 + Length.mm. + Oto.reading.4:Length.mm. + (1|ID), data=p2.df.clean, REML=T)#global mod
summary(fit_p2_glob)
sum_reml_p2 #selected mod
p_egg_p2     #one-sided t-test from selected mod

 #c
summary(fit.c.GLOBAL) #global mod
summary(fit.c.C1.REML) #selected mod
p_egg_c              # one-sided t-test from selected mod

#######################################################################
#########################################################################
#final graphs
############################################################################
#########################################################################


##########################
#how we get to plot_egg_final
#########################
##how we get these: p1.eggs/p2.eggs/coho.eggs
#p1
pink.egg.graph <- ggplot(p1.df.clean) + aes(x=Length..mm., y=Diameter..mm., color=Otolith.results) +
  geom_point(alpha=0.5) + scale_color_manual(values=c("blue", "orange"))+
  theme_cowplot() + 
  guides(color="none")
#theme(legend.position = c(0.03, 0.94))+
#coord_cartesian(ylim=c(5,8.5)) +
#scale_y_continuous(breaks=c(5,6,7,8), expand=c(0,0)))

sum.D <- summary(fit.p1.D)
summary(fit.p1.C1.REML)

fixef.p1 <-sum.D$coefficients

#p1.clean.wild <- p1.clean %>% filter(Otolith.results == "no mark")
#p1.clean.hatch <- p1.clean %>% filter(Otolith.results == "PORT ARMSTRONG")

upper.p1 <- fixef.p1[1] + (1.96*fixef.p1[2])
lower.p1 <- fixef.p1[1] - (1.96*fixef.p1[2])

p1.eggs <- pink.egg.graph + #geom_segment(y=fixef.p1[1], yend=fixef.p1[1], x= min(p1.clean$Length..mm.), xend=max(p1.clean$Length..mm.), color= "black", size=1.5)+
  #now I just need them error bars
  #geom_segment(y=upper.p1, yend=upper.p1, x= min(p1.clean$Length..mm.), xend=max(p1.clean$Length..mm.), color="black", linetype="dashed") + #upper
  #geom_segment(y=lower.p1, yend=lower.p1, x= min(p1.clean$Length..mm.), xend=max(p1.clean$Length..mm.), color="black", linetype= "dashed") + #lower
  labs(y="Egg diameter (mm)" , x="Length (mm)") +
  coord_cartesian(xlim=c(395,480)) +scale_x_continuous (breaks=c(400,420,440,460,480), expand=c(0,0)) +
  theme(plot.margin = margin(t=7,r=12,l=5,b=5))

#p2
names(p2.df.clean)

p2.egg.graph <- ggplot(p2.df.clean) + aes(x=Length.mm., y=Diameter, color=Oto.reading) +
  geom_point(alpha=0.5) + scale_color_manual(values=c("blue", "orange"))+
  theme_cowplot() + 
  guides(color="none")

wild.p2 <- fixef(fit.p2.lme.noint)[1]
hatch.p2 <- fixef(fit.p2.lme.noint)[2]

wild.subset.p2 <- p2.df.clean %>% filter(Oto.reading=="No Mark")
hatch.subset.p2 <- p2.df.clean %>% filter(Oto.reading=="PORT ARMSTRONG")

Intervals.p2 <- intervals(fit.p2.lme.noint)
hatchlow<- Intervals.p2$fixed[2,1]
hatchhigh <- Intervals.p2$fixed[2,3]
wildlow <- Intervals.p2$fixed[1,1]
wildhigh <- Intervals.p2$fixed[1,3]

p2.eggs <-p2.egg.graph + geom_segment(y=wild.p2, yend=wild.p2, x= min(wild.subset.p2$Length.mm.), xend=max(hatch.subset.p2$Length.mm.), color= "blue", size=1.5) +  #wild+
  geom_segment(y=hatch.p2, yend=hatch.p2, x= min(hatch.subset.p2$Length.mm.), xend=max(hatch.subset.p2$Length.mm.), color= "orange", size=1.5)+ #hatch
  #now just need to add error bars
  geom_segment(y= wildhigh,  yend=(wildhigh), x=  min(wild.subset.p2$Length.mm.), xend= max(wild.subset.p2$Length.mm.), color="black", linetype="dashed")+
  geom_segment(y= wildlow,  yend=(wildlow), x=  min(wild.subset.p2$Length.mm.), xend= max(wild.subset.p2$Length.mm.), color="black", linetype="dashed") +
  geom_segment(y= hatchlow,  yend=(hatchlow), x=  min(hatch.subset.p2$Length.mm.), xend= max(hatch.subset.p2$Length.mm.), color="black", linetype="dashed") +
  geom_segment(y= hatchhigh,  yend=(hatchhigh), x=  min(hatch.subset.p2$Length.mm.), xend= max(hatch.subset.p2$Length.mm.), color="black", linetype="dashed") + 
  labs(y= "Egg diameter (mm)",x="Length (mm)") +
  coord_cartesian(xlim=c(342,475), ylim=c(5,8)) +scale_x_continuous (breaks=c(350, 375,400,425,450, 475), expand=c(0,0)) +
  scale_y_continuous(breaks=c(5,6,7,8), expand = c(0,0))+
  theme(plot.margin = margin(t=7,r=12,l=5,b=5))

#coho
(coho.egg.graph <- ggplot(coho.clean) + aes(x=Length..mm., y=Diameter..mm., color=Wild.or.Hatch) +
    geom_point(alpha=0.5) + scale_color_manual(values=c("blue", "orange"), name=element_blank(), labels=c("Wild origin", "Hatchery origin"), breaks=c("wild", "hatchery"))+
    labs(y="Egg diameter (mm)", x="MEHP length (mm)")+
    theme_cowplot() +
    theme(legend.position = c(0.03, 0.94))+
    coord_cartesian(ylim=c(5,8.5)) +
    scale_y_continuous(breaks=c(5,6,7,8), expand=c(0,0)))

hatch.temp <- coho.clean %>% filter(Wild.or.Hatch=="hatchery")
wild.temp <- coho.clean %>% filter(Wild.or.Hatch=="wild")

hmin.length <- min(hatch.temp$Length..mm.)
hmax.length <- max(hatch.temp$Length..mm.)
wmin.length <- min(wild.temp$Length..mm.)
wmax.length <- max(wild.temp$Length..mm.)

coho.eggs <- coho.egg.graph + 
  geom_segment(y=Intervals$fixed[1,1], yend=Intervals$fixed[1,1], color="black", x=hmin.length, xend=hmax.length, linetype="dashed") +
  geom_segment(y=Intervals$fixed[1,3], yend=Intervals$fixed[1,3], color="black", x=hmin.length, xend=hmax.length, linetype="dashed") +
  #hatch
  geom_segment(y= Intervals$fixed[2,1],  yend=(Intervals$fixed[2,1]), x= wmin.length, xend=wmax.length, color="black", linetype="dashed")+ 
  geom_segment(y= Intervals$fixed[2,3],  yend=(Intervals$fixed[2,3]), x= wmin.length, xend=wmax.length, color="black", linetype="dashed")+ 
  #wild
  geom_segment(y = cf.fix[1], yend = cf.fix[1], color="orange", size=1.5, x= hmin.length, xend = hmax.length ) + geom_segment(y = cf.fix[2], yend= cf.fix[2], color="blue", size=1.5, x= wmin.length, xend = wmax.length )+
  coord_cartesian(ylim=c(5,8.3), xlim=c(473,650)) +
  scale_y_continuous(breaks=c(5,6,7,8), expand=c(0,0))+
  scale_x_continuous (breaks=c(500,550,600,650), expand=c(0,0))+
  guides(color="none") +
  labs(x="Length (mm)") +
  theme(plot.margin = margin(t=7,r=12,l=5,b=5))

##
xlab <-coho.eggs$labels$x
ylab <- coho.eggs$labels$y

coho.eggs_nolabs <- coho.eggs + labs(x=element_blank(), y=element_blank())
p1.eggs_nolabs <- p1.eggs + labs(x=element_blank(), y=element_blank())
p2.eggs_nolabs <- p2.eggs + labs(x=element_blank(), y=element_blank())

#nice. Now combine them?
eggs_combined_draft <- plot_grid(p1.eggs_nolabs,p2.eggs_nolabs,coho.eggs_nolabs, scale=0.95, nrow=3) #ok, that looks good. Good job

#now add axes labels
##set it up first
temp_egg <- plot_grid (NULL, eggs_combined_draft, ncol = 2, rel_widths =c(1,6))
egg_base <- plot_grid(eggs_combined_draft, NULL, ncol = 1, rel_heights = c(9,1)) #might need to make this smaller a lil bit
#and tehn actually add them labels
plot_egg_final <- ggdraw(egg_base) + draw_label(xlab, x = 0.55, y = 0.11, size = 15) + 
  draw_label ((ylab), angle= 90, x = 0.02, y = 0.6, size = 15) +
  draw_label ("Even-year pink", x = 0.85, y = 0.98, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Odd-year pink", x = 0.85, y = 0.67, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Coho", fontfamily = "Arial", fontface="bold", size = 13, x = 0.85, y = 0.37) 

plot_egg_final



#############################
#get plot GSI final
##########################

###GSI. Copied from Russia graphs and now edited. I'll remove the labels. And... how to panel these? Three across? All stacked on top of each other?
#coho  #c.GSI.clean - get those datsets
#coho gsi
ggGSI_coho <- ggplot(c.GSI.clean) + aes(x=Length..mm., y=GSI, color=Wild.or.Hatch) + geom_point(size=3) +
  #geom_smooth(method="lm") + #ADD THIS BACK IN if we decide to put length in the model
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(x="Length (mm)") +
  coord_cartesian(ylim=c(13,30), xlim=c(473,650)) +scale_x_continuous (breaks=c(500,550,600,650), expand=c(0,0))+
  scale_y_continuous(breaks=c(15, 20, 25, 30), expand=c(0,0)) #+
#theme(plot.margin = margin(t=7,r=12,l=5,b=5))

#pink2020  #p1.clean
range(p1.clean$Length..mm.) # 397 476
range(p1.clean$GSI) #14.43570 24.45902

ggGSI_p20 <-ggplot(p1.clean) + geom_point(aes(x=Length..mm., y=GSI, color=Otolith.results), size=3)+
  #geom_smooth(aes(x=Length..mm., y=GSI), method = "lm", color="black") +
  scale_color_manual(values=c("orange", "blue"), name=element_blank(), labels = c("Hatchery", "Wild"), limits=c("PORT ARMSTRONG", "no mark")) +
  theme_cowplot() +
  #guides(color="none") +
  labs(x="Length (mm)") +
  coord_cartesian(ylim=c(13.9,25.0), xlim=c(395,480))+ scale_x_continuous(breaks=c(400,420,440,460,480), expand=c(0,0))+
  scale_y_continuous(breaks=c(15, 20, 25), expand=c(0,0)) + #for legend.poition, can ise x=0.85 istead of 0.05
  theme(legend.position = c(0.02, 0.90), plot.margin= margin(t=7,r=12,l=5,b=5))  +
  theme(
    legend.title=element_blank(),
    legend.box.background = element_rect(),
    legend.box.margin = margin(-5, 3, 2, 2) #yes, finally I made a box.
  ) #this is looking ok. Note the negative number. This is where I'm adjusting the LEGEND POSITION. needs to adjust a tad down


#pink 2021
max(p2.GSI.clean$GSI.1) #26.7
min(p2.GSI.clean$GSI.1) #11.8
max(p2.GSI.clean$Length.mm.) #462
min(p2.GSI.clean$Length.mm.) #344
ggGSI_pODD<-ggplot(p2.GSI.clean) + geom_point(aes(x=Length.mm., y=GSI.1, color=Oto.reading), size=3)+
  #geom_smooth(aes(x=Length.mm., y=GSI.1), method = "lm", color="black") +
  scale_color_manual(values=c("blue", "orange")) + #bc gatch is listed first here
  theme_cowplot() + 
  guides(color="none") + #ooh. How to move a legend into the frame
  labs(y="GSI", x="Length (mm)")+
  coord_cartesian(ylim=c(11.0,28.0), xlim=c(340,475))+  #might try max of 465 instead
  scale_x_continuous(breaks=c(350,375,400,425,450,475), expand=c(0,0)) +  
  scale_y_continuous(breaks = c(12, 20, 28), expand=c(0,0))

##
ggGSI_p20_2 <- ggGSI_p20 + labs(x=element_blank(), y=element_blank())
ggGSI_p21_2 <- ggGSI_pODD + labs(x=element_blank(), y=element_blank())
ggGSI_c2 <- ggGSI_coho + labs(x=element_blank(), y=element_blank())

ggGSI_grid <- plot_grid(ggGSI_p20_2, ggGSI_p21_2, ggGSI_c2, scale=0.95, nrow=3) 

temp_ggGSI <- plot_grid(NULL, ggGSI_grid, ncol = 2, rel_widths =c(1,6))
ggGSI_base <- plot_grid(ggGSI_grid, NULL, ncol = 1, rel_heights = c(9,1)) 

xlab_GSI <- ggGSI_coho$labels$x
ylab_GSI <- ggGSI_coho$labels$y

plot_GSI_final <- ggdraw(ggGSI_base) + draw_label(xlab_GSI, x = 0.55, y = 0.11, size = 15) + 
  draw_label ((ylab_GSI), angle= 90, x = 0.02, y = 0.6, size = 15) +
  draw_label ("Even-year pink", x = 0.85, y = 0.98, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Odd-year pink", x = 0.85, y = 0.67, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Coho", fontfamily = "Arial", fontface="bold", size = 13, x = 0.85, y = 0.37) 

plot_GSI_final

#############################################
#plot_GSI_final + plot_egg_final
#combine plot_egg_final and plot_GSI_final
#############################################
library(patchwork)
design2 <- "
111111#22222222
111111#22222222
111111#22222222
"


plot_GSI_final + plot_egg_final + plot_layout(design=design2)
#8 by 10 works decent with design 2 spacing

dev.new (width = 10, height = 8, unit = "in", noRStudioGD = T); last_plot() #perfect
ggsave ("Plots/GSI_Eggs_combined1.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()



###################################################################
#####################################################################
#female length test and graph

##################################################################
#################################################################



#############################3
#lovers snag vs. seine size test, for Sam
###########################
#females caught from lovers
#weird fish included
p2_catch_data_f <- read.csv("Data/Snag/p2_female_snag.csv")

p2_female_lovers_snag <- p2_catch_data_f %>%
  filter(Location == "Lovers",
         Catch.method == "snag")
p2_female_lovers_seine <- p2_catch_data_f %>%
  filter(Location == "Lovers",
         Catch.method == "seine")

t.test(p2_female_lovers_snag$Length.mm., p2_female_lovers_seine$Length.mm.) #not sig different

seine_v_snag_f <- rbind(p2_female_lovers_snag, p2_female_lovers_seine) #just the lovers cove seine and snag data

snag_f <-ggplot(seine_v_snag_f) + aes(x=Catch.method, y=Length.mm.) + geom_boxplot() + geom_jitter() + theme_cowplot()+
  ggtitle("Lovers Cove - females")+
  labs(x="Catch method", y="Length")

#most of them were snagged. that's what this plot tells us



#males caught from lovers
p2_catch_data_m <- read.csv("Data/Snag/p2_male_snag.csv")

p2_male_lovers_snag <- p2_catch_data_m %>%
  filter(Location == "Lovers",
         Catch.method == "snag")

p2_male_lovers_seine <- p2_catch_data_m %>%
  filter(Location == "Lovers",
         Catch.method == "seine")

t.test(p2_male_lovers_snag$Length.mm., p2_male_lovers_seine$Length.mm.) #not sig diff

seine_v_snag_m <- rbind(p2_male_lovers_snag, p2_male_lovers_seine)

snag_m <- ggplot(seine_v_snag_m) + aes(x=Catch.method, y=Length.mm.) + geom_boxplot() + geom_jitter() + theme_cowplot()+
  ggtitle("Lovers Cove - males")+
  labs(x="Catch method", y="Length")

#plot of catch type comparison for lovers cove
snag_f + snag_m


#############################################################################
##############################################################################
#the lovers vs sashin test that charlie and sam wanted
##(re-visit notes please)
################################################################################
#################################################################################
#sashin lovers comparison that justifies pooling

#uknown fish are excluded from analysis, but included in these graphs



##p1
#GSI
names(p1.clean) #GSI
names(p1.df.clean) #eggs
p1.clean$Location
a <- ggplot(p1) + aes(y=GSI, x=Location, color=Otolith.results) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() + scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2020, separated by oto mark")

b <- ggplot(p1) + aes(y=GSI, x=Location) + #includes unknown
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter()+
  theme_cowplot()+
  ggtitle("Pink 2020")

#Eggs #does not include unknown
c <-ggplot(p1.df.clean) + aes(y=Diameter..mm., x=Location.x, color=Otolith.results) +
  geom_boxplot() + scale_color_manual(values=c("blue", "orange"))+ #+geom_violin()
  theme_cowplot()+
  ggtitle("Pink 2020, separated by oto mark")
d <- ggplot(p1.df.clean) + aes(y=Diameter..mm., x=Location.x) +
  geom_boxplot()+
  theme_cowplot()+
  ggtitle("Pink 2020")
#Length - includes uknown
L1 <- ggplot(p1) + aes(y=Length..mm., x=Location, color=Otolith.results) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() + scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2020, length separated by oto mark")

L11 <- ggplot(p1.clean) + aes(y=Length..mm., x=Location, color=Otolith.results) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() + scale_color_manual(values=c("blue", "orange"))+
  theme_cowplot()+
  ggtitle("Pink 2020, length separated by oto mark")


L2 <- ggplot(p1) + aes(y=Length..mm., x=Location) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() +
  theme_cowplot()+
  ggtitle("Pink 2020")

L22 <- ggplot(p1.clean) + aes(y=Length..mm., x=Location) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() +
  theme_cowplot()+
  ggtitle("Pink 2020")

##p2
#quick wrangle
p2.GSI.with.unknown <- p2.GSI %>%
  mutate(Oto.reading = ifelse(Oto.reading %in% c("Overground", "No Oto"), "unknown", Oto.reading),
         GSI.1 = GSI.measure.1.g./Fish.weight.g.)

#GSI-good - change so it includes uknown, werid fish
names(p2.GSI.with.unknown)
e <- ggplot(p2.GSI.with.unknown) + aes(y=GSI.1, x=Location, color=Oto.reading) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter() + scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2021, separated by oto mark") #does not work now becasu I 
f <- ggplot(p2.GSI.with.unknown) + aes(y=GSI.1, x=Location) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter()+
  theme_cowplot()+
  ggtitle("Pink 2021")

#eggs #meh, I dont really want to do this graph
names(p2.df.clean) 
g <- ggplot(p2.df.clean) +aes(y=Diameter, x=Location, color=Oto.reading) +
  geom_boxplot() + scale_color_manual(values=c("blue", "orange"))+ #+geom_violin()
  theme_cowplot()+
  ggtitle("Pink 2021, separated by oto mark")
h <- ggplot(p2.df.clean) +aes(y=Diameter, x=Location) +
  geom_boxplot()+ #+geom_violin()
  theme_cowplot()+
  ggtitle("Pink 2021")

#length
L3 <- ggplot(p2.GSI.with.unknown) + aes(x=Location, y=Length.mm., color=Oto.reading) +
  geom_boxplot() + geom_jitter() + scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2021, length separated by oto mark")

L33 <- ggplot(p2.GSI.clean) + aes(x=Location, y=Length.mm., color=Oto.reading) +
  geom_boxplot() + geom_jitter() + scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2021, length separated by oto mark")

L4 <- ggplot(p2.GSI.with.unknown) + aes(x=Location, y=Length.mm.) +
  geom_boxplot() + geom_jitter() +# scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2021")

L44 <- ggplot(p2.GSI.clean) + aes(x=Location, y=Length.mm.) +
  geom_boxplot() + geom_jitter() +# scale_color_manual(values=c("blue", "orange", "springgreen"))+
  theme_cowplot()+
  ggtitle("Pink 2021")

#length graphs for sharing
##NOTE: unknown and werid fish are incldued in these graphs, but not in analysis
(L1+L3)/(L2+L4) #with unknowns and weirds
(L11+L33)/(L22+L44) #without unknowns and werids. I think use this one


##01/23/24 Sam said use these graphs in supplement. So, clean it up. But first clear mem and activate
##p1.clean and p2.GSI.clean. Cleaning up plots L11+L33
L111 <- ggplot(p1.clean) + aes(y=Length..mm., x=Location, color=Otolith.results) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Even-year pink")+
  labs(y="Length (mm)", x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(390, 480), breaks=c(400, 440, 480), expand=c(0,0))
#?position_dodge()

L333 <- ggplot(p2.GSI.clean) + aes(x=Location, y=Length.mm., color=Oto.reading) +
  geom_boxplot() + geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Odd-year pink")+
  labs(y=element_blank(), x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(340, 480), breaks=c(360, 400, 440, 480), expand=c(0,0))


L444 <- L111+L333 + plot_layout(guides="collect")
L444

dev.new (width = 10, height = 5, unit = "in", noRStudioGD = T); last_plot()
ggsave(filename = "Plots/Justification for combining Lovers and Sashin.jpg", width=10, height=5)
dev.off()

##quick sashin vs lovers t.test, as charlie requested
p1_s <- p1.clean %>% filter(Location=="Sashin Creek", Otolith.results=="no mark")
p1_l <- p1.clean %>% filter(Location=="Lovers Cove", Otolith.results=="no mark")

p1_ttest <-t.test(p1_s$Length..mm.,p1_l$Length..mm.) #p val = 0.839


p2_s <- p2.GSI.clean %>% filter(Location=="Sashin", Oto.reading == "No Mark")
p2_l <-p2.GSI.clean %>% filter(Location=="Lovers", Oto.reading == "No Mark")

p2_ttest<-t.test(p2_s$Length.mm.,p2_l$Length.mm.)  #p val = 0.08907

names(p1_ttest)

sashin_lovers_length_comparison <- data.frame(p1=c(p1_ttest$statistic, p1_ttest$parameter, p1_ttest$p.value),
                                               p2=c(p2_ttest$statistic, p2_ttest$parameter, p2_ttest$p.value))

rownames(sashin_lovers_length_comparison)= c("t", "df", "pvalue")

write.csv(x=sashin_lovers_length_comparison, file="Results/sashin_lovers_length_comparison.csv")






#combine combine
(a+e)/(b+f) #GSI sashin/lovers/armstrong comparison plot
(c+g)/(d+h) #egg sashin/lovers/armstrong comparison plot



#####
#create csv for graphing fig3 R script
write.csv(c.GSI.clean, "Data/c.GSI.clean.csv")




################################################################
#length comparison
##################################################################
##t test for 3 female groups
p1_f_L_ttest #not sig
p2_f_L_ttest #not sig
c_f_L_ttest #sig


stat_cfl <- c_f_L_ttest$statistic
p_cfl <- c_f_L_ttest$p.value
df_cfl <- c_f_L_ttest$parameter
mean_length_c_h
mean_length_c_w
sd_length_c_h
sd_length_c_w


stat_p1fl <- p1_f_L_ttest$statistic
p_p1fl <- p1_f_L_ttest$p.value
df_p1fl <- p1_f_L_ttest$parameter
mean_length_p1_h 
mean_length_p1_w 
sd_length_p1_h 
sd_length_p1_w 

stat_p2fl <- p2_f_L_ttest$statistic
p_p2fl <- p2_f_L_ttest$p.value
df_p2fl <- p2_f_L_ttest$parameter
mean_length_p2_h 
mean_length_p2_w 
sd_length_p2_h 
sd_length_p2_w 

female_length_results <- data.frame(
  t= c(stat_p1fl, stat_p2fl, stat_cfl),
  df = c(df_p1fl, df_p2fl, df_cfl),
  p = c(p_p1fl, p_p2fl, p_cfl),
  hatch_mean = c(mean_length_p1_h, mean_length_p2_h, mean_length_c_h),
  hatch_sd = c(sd_length_p1_h, sd_length_p2_h, sd_length_c_h),
  wild_mean = c(mean_length_p1_w, mean_length_p2_w, mean_length_c_w),
  wild_sd = c(sd_length_p1_w, sd_length_p2_w, sd_length_c_w)
)

rownames(female_length_results) <- c("pink 2020", "pink 2021", "coho")

write.csv(female_length_results, "Results/Female length results.csv")

##t test for 3 male groups
###let's do this over in the male dataframe?

################################################################################
#FEMALE
##GSI - two sample, one tailed t-test, sometimes with same # of samples, sometimes with diff
library(pwr)
#library(effsize)
library(effectsize)
library(simr)
?pwr.t.test #one and two samples for equal sample sizes
?pwr.t2n.test #two samples of differnt sizes
################################################################################
#GSI- power analysis for a t-test
#Eggs - power analysis for a linear mixed effects model? is this more complex than a linear model?


#P1############################

#GSI
##data frame: p1.clean
#(p1.GSI.t.test <- t.test(p1.hatch.GSI.ttest$GSI, p1.wild.GSI.ttest$GSI, alternative="greater", var.equal=T))
#first sample will be hatch
length(p1.hatch.GSI.ttest$GSI) #36
#second sample will be wild
length(p1.wild.GSI.ttest) #15

#d is effect size. how do I choose effect size again? Ecalculate cohen's D?
##is this the regression coefficient?
#what is the effect size of a t-test?
#is the effect size different for a one-sided test? Is the power different for a one-sided test?
p1.GSI.t.test
effectsize(p1.GSI.t.test) #cohens D is -0.41

pwr_p1_GSI <- pwr.t2n.test(n1=36, n2=10, d= -0.41, sig.level=0.05, power = NULL, 
                           alternative="greater" #alternative is greater because hypothesis is that hatch fish can invest more in their GSI
                           ) #power is 0.26. Damn , low power. Even less for one-sided


#eggs
##model: fit.p1.C1.REM

pwr_p1_egg <- simr::powerSim(fit.p1.C1.REML,  #likelihood ratio test for power
                             test = simr::fcompare(Diameter..mm.~1) ,
                             nsim = 1000
) 

pwr_sum_p1_egg <- summary(pwr_p1_egg)

#experimenting with a power curve
test_p1 <- simr::powerCurve(fit.p1.C1.REML) #trying out a power curve
#plot(test_p1) #yeah that did not work


#P2###################################

#GSI
##data frame: p2.GSI.clean.relevant
#(p2.GSI.t.test <- t.test(GSI.for.ttest.hatch$GSI.1,GSI.for.ttest.wild$GSI.1, alternative="greater", var.equal=T))
p2.GSI.t.test
effectsize(p2.GSI.t.test)

length(GSI.for.ttest.hatch$GSI.1) #50 fish
length(GSI.for.ttest.wild$GSI.1) #28 fish 

pwr_p2_GSI <- pwr.t2n.test(n1=50, n2=28, d= -0.01, sig.level=0.05, power = NULL, alternative="greater") # power = 0.05020049. That's really low.

#eggs

pwr_p2_egg <- simr::powerSim(fit.p2.B4.reml,  #likelihood ratio test for power
                            test = simr::fcompare(Diameter~1) ,
                            nsim = 1000
) 

pwr_sum_p2_egg <- summary(pwr_p2_egg)


#C###################################

#GSI
##data frame: c.GSI.clean
#(coho.GSI.t.test <- t.test(c.GSI.hatch.t.test$GSI,c.GSI.wild.t.test$GSI, alternative="greater", var.equal=T))

length(c.GSI.hatch.t.test$GSI) #27
length(c.GSI.wild.t.test$GSI) #27 (for this sample). looks like I can use the regular pwr.t.test since sample size is equal
effectsize(coho.GSI.t.test) #0.40

pwr_c_GSI <- pwr.t.test(n=27, d= 0.40, sig.level=0.05, power = NULL, alternative="greater"
                        ) #power = 0.3027827. such low power. increases a bit with one-sided but still low

#eggs #data frame: coho.clean.relevant
##fit.c.C1.REML
## fit.lme.no.int
?simr::powerCurve #powersim over different sample sizes
?simr::powerSim #estimat power by simulation

#simr::powerSim(fit.c.C1.REML, 
 #              test = simr::fixed(xname= factor(coho.clean.relevant$Wild.or.Hatch), method="t") , #select fixed effect and test type. Seemed to work this time. runs inconsistently
               #when I include test, the power is 0. I do not think the power is 0.
               #test=fcompare(Diameter..mm.~1) #alternative attempt for test
               #seed =  ,
               #nsim = , #default is 1000. Takes a bit but not crazy long
               #  )

simr::powerSim(fit.c.C1.REML, 
              # test = simr::fixed(xname= Wild.or.Hatch, method="t") , 
)

pwr_c_egg <- simr::powerSim(fit.c.C1.REML,  #this also runs
                test = simr::fcompare(Diameter..mm.~1) ,
               nsim = 1000
)  #so I trust this one over the unknown test, but the t-test should work and idk why it does not.

pwr_sum_c_egg <- summary(pwr_c_egg)

## make csv########################################
power_females <- data.frame(GSI_power = c(pwr_p1_GSI$power, pwr_p2_GSI$power, pwr_c_GSI$power), Egg_power=c( pwr_sum_p1_egg$mean,pwr_sum_p2_egg$mean,pwr_sum_c_egg$mean)) #need to add in egg power still

rownames(power_females) = c("Pink even", "Pink odd", "Coho")

write.csv(power_females, "Results/Female post-hoc power analysis.csv")

############################################################################
#power analysis for linear mixed effect models: figuring out
#https://www.r-bloggers.com/2009/09/power-analysis-for-mixed-effect-models-in-r/
#"Data Analysis Using Regression and Multilevel/Hierarchical Models" by Andrew Gelman and Jennifer Hill
#For the current analysis we needed to know three things: effect size, sample size, and estimates of population variance
##effect size will be... the parameter estimate for wild/hatch??
##we know sample size
##we can find out population variance
###something about hyperparameters. Can we use pooled variance from nlme? What are hyperparameterss? Oh, fake parameters. I have real ones
#https://humburg.github.io/Power-Analysis/simr_power_analysis.html <- really good resource.

?lme4
?nlme
?nlme::pooledSD
?lmerTest
?lmtest
??simr #supposed to work with lme4
#?Hmisc::simRegOrd

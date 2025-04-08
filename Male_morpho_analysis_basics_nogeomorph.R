
#Phenotypic divergence paper
##just the males linear morph, removed geomorph
## Alex Reich
##4/7/25



#df.coho.simple <- df.coho %>% select(ID, snoutL, length, Wild.or.hatch, depth)

#write.csv(df.coho.simple, "Data/Coho_morpho_data.csv")

#libraries
library(tidyverse)


###coho analysis: snout
df.coho <- read.csv("Data/Coho_morpho_data.csv")

coho_all_snout_05 <- lm(snoutL ~ length + factor(Wild.or.hatch), df.coho)
sum_c_5<-summary(coho_all_snout_05)

#get values for table, mean, sd, ect. ADD THIS TO RESUTLS TABLE
coho_h <- df.coho %>% filter(Wild.or.hatch=="H")
coho_w <- df.coho %>% filter(Wild.or.hatch=="W")

mean_s_h <- 10* mean(coho_h$snoutL) #that measurement is in cm I think, times 10 is mm
mean_s_w <- 10*  mean(coho_w$snoutL)

sd_s_h <- 10*sd(coho_h$snoutL)
sd_s_w <- 10*sd(coho_w$snoutL) #qc'ed, checks out.

df_coho_snout_results <- data.frame(t=sum_c_5$coefficients[3,3], df= sum_c_5$df[2], p_one_sided=1-sum_c_5$coefficients[3,4]/2,
                                    length_t = sum_c_5$coefficients[2,3], length_p_two_Sided = sum_c_5$coefficients[2,4], 
                                    hatchery_mean= mean_s_h, hatchery_sd = sd_s_h,
                                    wild_mean=mean_s_w, wild_sd=sd_s_w)


###############################
#coho: depth

###########
W.depth.all.coho <- df.coho %>% filter (Wild.or.hatch=="W")
H.depth.all.coho <- df.coho %>% filter (Wild.or.hatch=="H")
#t.depth.c <- t.test(x=H.depth.all.coho$depth, y=W.depth.all.coho$depth, alternative= "less", var.equal = T)
#t.depth.c

#get sd and mean values for my tables
sd_d_hc <- sd(H.depth.all.coho$depth) #colo all depth sd
sd_d_wc <- sd(W.depth.all.coho$depth) #coho all depth sd
names(W.depth.all.coho)

mean_d_hc <- mean(H.depth.all.coho$depth)
mean_d_wc <- mean(W.depth.all.coho$depth)

aov.depth.c <- lm(depth ~ length + factor(Wild.or.hatch), df.coho)
sum_c_depth <- summary(aov.depth.c)
#plot(aov.depth.c) #looks fine


#prepare the table
df_coho_depth_results <- data.frame(t=sum_c_depth$coefficients[3,3], df= sum_c_depth$df[2], p_one_sided=1-sum_c_depth$coefficients[3,4]/2,
                                    length_t = sum_c_depth$coefficients[2,3], length_p_two_Sided = sum_c_depth$coefficients[2,4], 
                                    hatchery_mean= mean_d_hc, hatchery_sd = sd_d_hc,
                                    wild_mean=mean_d_wc, wild_sd=sd_d_wc)

#some density plots
ggplot(df.coho) + aes(x=snoutL) + geom_density()
################################

#05/23/24, I need AIC tables
#snout
coho_glob_snout <- lm(snoutL ~ length * factor(Wild.or.hatch), df.coho)
coho_base_s <- lm(snoutL ~factor(Wild.or.hatch), df.coho)
coho_base_d <- lm(depth ~factor(Wild.or.hatch), df.coho)
aov.glob.c <- lm(depth ~ length * factor(Wild.or.hatch), df.coho)

AIC(coho_glob_snout, coho_all_snout_05, coho_base_s)
AIC(aov.glob.c, aov.depth.c, coho_base_d)

#depth

#08/07/24- the BIG supplementary table from hell
summary(coho_glob_snout)
#summary(coho_base_s)
summary(coho_all_snout_05)
df_coho_snout_results

summary(aov.glob.c)
#summary(coho_base_d)
summary(aov.depth.c)
df_coho_depth_results



#PINK 1 SNOUTS AND DEPTHS

#created from former morpho code
#write.csv(df.pink2, "df.pink2.csv") #do I need this?

#df.pink2.simple <- df.pink2 %>% select(ID, snout, length, Wild.or.hatch, depth)
#write.csv(df.pink2.simple, "Data/pink_morpho_data.csv")

df.pink2 <- read.csv("Data/pink_morpho_data.csv")


#first separate wild and hatchery
wild_D.p.int <-  df.pink2 %>% filter(Wild.or.hatch=="W")
hatch_D.p.int <-  df.pink2 %>% filter(Wild.or.hatch=="H")

wild_D.p <- wild_D.p.int$depth
hatch_D.p<- hatch_D.p.int$depth


#first separate wild and hatchery
wild_S.p <-  df.pink2 %>% filter(Wild.or.hatch=="W")
hatch_S.p <-  df.pink2 %>% filter(Wild.or.hatch=="H")


#pink 2020 sd and mean
sd_s_hp1 <- sd(hatch_S.p$snout)*10
sd_s_wp1 <- sd(wild_S.p$snout)*10
mean_s_hp1 <- mean(hatch_S.p$snout)*10
mean_s_wp1 <- mean(wild_S.p$snout)*10

names(wild_S.p)
sd_d_hp1 <- sd(hatch_S.p$depth)
sd_d_wp1 <-  sd(wild_S.p$depth)
mean_d_hp1 <- mean(hatch_S.p$depth)
mean_d_wp1 <-  mean(wild_S.p$depth)

#find sample size:
length(hatch_S.p$depth)
length(wild_S.p$depth)



#exploratory
qqnorm(df.pink2$snout)
qqline(df.pink2$snout)
qqnorm(log(df.pink2$snout))
qqline(log(df.pink2$snout))
ggplot(df.pink2) + aes(x=log(snout)) + geom_density()
ggplot(df.pink2) + aes(x=log(depth)) + geom_density()
ggplot(df.pink2) + aes(x=depth) + geom_density()

ggplot(df.pink2) + aes(x=depth, y=length) + geom_point()+ geom_smooth()#added
ggplot(df.pink2) + aes(x=depth, y=log(length)) + geom_point()+ geom_smooth()
ggplot(df.pink2) + aes(x=log(depth), y=log(length)) + geom_point()+ geom_smooth()

#actual analysis - no more logging
p1_snout_mod <- lm(snout~ length + factor(Wild.or.hatch), data=df.pink2)
sum_p_snout <- summary(p1_snout_mod)
#plot(p1_snout_mod)

p1_depth_mod <- lm(depth~ length + factor(Wild.or.hatch), data=df.pink2)
sum_p_depth <- summary(p1_depth_mod)
#plot(p1_depth_mod)

#interaction effects?
mod.p.d.l.int <- aov(depth ~ length + Wild.or.hatch + length*Wild.or.hatch, df.pink2)
summary.lm(mod.p.d.l.int) #not sig

mod.p.s.l.int <- aov(snout ~ Wild.or.hatch+length + length*Wild.or.hatch, df.pink2)
summary.lm(mod.p.s.l.int) #no int effects

#added 5/23/24
base_s <- lm(snout~ factor(Wild.or.hatch), data=df.pink2)
null_s <- lm(snout~ 1, data=df.pink2)
base_d <- lm(depth~ factor(Wild.or.hatch), data=df.pink2)
null_d <- lm(depth~ 1, data=df.pink2)

AIC(mod.p.d.l.int, p1_depth_mod, base_d ) #depth
AIC(mod.p.s.l.int,p1_snout_mod, base_s) #snout
summary(null_s) #for 6/6/24 table edits snout and depth p1
summary(base_s)
#get the one-sided p value from the t-value
(1-0.429)/2
0.429/2
1-0.429/2 #this one
pt(q=-0.797, lower.tail=F, df=58)
summary(p1_snout_mod)
pt(q=1.237, lower.tail=T, df=58) #this is teh answer
summary.lm(mod.p.s.l.int) #0.890
pt(q=-1.282, lower.tail=F, df=58)

summary(null_d) #for 6/6/24 table edits depth p1
summary(base_d)
summary(p1_depth_mod) #selected mod for supp table from hell 08/2024
summary.lm(mod.p.d.l.int) #global mod for supp table from hell
#broom::tidy(mod.p.d.l.int)


summary(null_s) #for 6/6/24 table edits snout p1
summary(base_s)
summary(p1_snout_mod) #selected mod for supp table from hell 08/2024
summary.lm(mod.p.s.l.int) #global mod for supp table from hell

#results for table
#full_mods_male <- tibble(summary.lm(mod.p.s.l.int), summary.lm(mod.p.d.l.int) )
#selected_mods_male <- tibble() 


#results
df_pink_snout_results <- data.frame(t=sum_p_snout$coefficients[3,3], df= sum_p_snout$df[2], p_one_sided=1-sum_p_snout$coefficients[3,4]/2,
                                    length_t = sum_p_snout$coefficients[2,3], length_p_two_Sided = sum_p_snout$coefficients[2,4], 
                                    hatchery_mean= mean_s_hp1, hatchery_sd = sd_s_hp1,
                                    wild_mean=mean_s_wp1, wild_sd=sd_s_wp1)
  
df_pink_depth_results <- data.frame(t=sum_p_depth$coefficients[3,3], df= sum_p_depth$df[2], p_one_sided=1-sum_p_depth$coefficients[3,4]/2,
                                    length_t = sum_p_depth$coefficients[2,3], length_p_two_Sided = sum_p_depth$coefficients[2,4], 
                                    hatchery_mean= mean_d_hp1, hatchery_sd = sd_d_hp1,
                                    wild_mean=mean_d_wp1, wild_sd=sd_d_wp1)

###########################################

#line 1495- formal length tests, hatch vs wild, both for coho and pink 2020
#################################
#05/19/22
#coho
names(df.coho) #my dataframe, all
#df.coho$Wild.or.hatch
df.coho.W <- df.coho %>% filter(Wild.or.hatch=="W")
df.coho.H <- df.coho %>% filter(Wild.or.hatch=="H")
t.test(df.coho.H$length, df.coho.W$length, alternative="two.sided", var.equal=T)
sd(df.coho.H$length)
sd(df.coho.W$length)
#t.test(df.coho.H$length, df.coho.W$length, alternative="two.sided", var.equal=F)


#plot coho length
##likely do a violin wiht sina overlay
##or can also do a boxplot
library(patchwork)
ggplot(df.coho) + aes(x=Wild.or.hatch, y=length) + geom_boxplot()+ geom_jitter()


#even-year pinks length
#df.pink2
names(df.pink2)
df.evenpink.W <- df.pink2 %>% filter(Wild.or.hatch=="W")
df.evenpink.H <- df.pink2 %>% filter(Wild.or.hatch=="H")
t.test(df.evenpink.H$length, df.evenpink.W$length, alternative="two.sided", var.equal=T)
sd(df.evenpink.H$length)
sd(df.evenpink.W$length)

###################################

#odd-year pink males results
###

#pink 2021 linear morph data
p2.males <- read.csv("Data/Male.p2.Rdata.2.csv")


#exploratory graphs
#use graphs you already made, dummy!
#,TK can't find them
names(p2.males)
hist(p2.males$Body.depth.mm.)
hist(p2.males$Snout.length.mm.)
hist(p2.males$Length.mm.)
ggplot(p2.males) + aes(x=Length.mm., y=Body.depth.mm.) + geom_point() 
ggplot(p2.males) + aes(x=Length.mm., y=Snout.length.mm.) + geom_point()

#additions: 05/03/22 - more exploratory graphs
library(ggpubr)
ggqqplot(p2.males$Snout.length.mm.)
ggqqplot(log(p2.males$Snout.length.mm.))

ggqqplot(p2.males$Body.depth.mm.)
ggqqplot(log(p2.males$Body.depth.mm.)) #log that!1

ggqqplot(p2.males$Length.mm.)

ggplot(p2.males) + aes(x=Snout.length.mm.) + geom_density()
ggplot(p2.males) + aes(x=log(Snout.length.mm.)) + geom_density()

ggplot(p2.males) + aes(x=Body.depth.mm.) + geom_density()
ggplot(p2.males) + aes(x=log(Body.depth.mm.)) + geom_density() # I see no reason to log either of these


#some graphs I'm skipping over for now, will revisit later if they contribute to the final figures - #line 572 of Morpho_yr2_pinks.Rmd

#identyifying the strays, not sure if needed
p2.males1 <- p2.males %>% filter(Otolith.reading!= "Overground")
strays <- p2.males1 %>% filter(Location != "Armstrong")
strays <- strays %>% filter(Otolith.reading == "PORT ARMSTRONG")


#hypothesis tests
mod.p2021.depth <- lm(Body.depth.mm. ~ Length.mm. + factor(Otolith.reading),data=p2.males1)
sum_depth_p2 <- summary(mod.p2021.depth)
mod.p2021.snout <- lm(Snout.length.mm. ~ Length.mm. + factor(Otolith.reading),data=p2.males1)
sum_snout_p2 <- summary(mod.p2021.snout)


ggplot(p2.males1) + aes(x=Date, y=Body.depth.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))
ggplot(p2.males1) + aes(x=Date, y=Snout.length.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))
#maybe change above so date is in order...




##############################################

#means and sd's, for data table

w<-p2.males1 %>% filter(Otolith.reading=="No Mark") 
mean_d_wp2 <- mean(w$Body.depth.mm.) #mean of w male depth for p2

h<-p2.males1 %>% filter(Otolith.reading=="PORT ARMSTRONG") 
mean_d_hp2 <- mean(h$Body.depth.mm.) #mean of h male depth for p2

sd_d_wp2 <- sd(w$Body.depth.mm.) #sd's, for my table
sd_d_hp2 <- sd(h$Body.depth.mm.)

mean_s_wp2 <- mean(w$Snout.length.mm.) #mean of w snout for p2
length(w$Snout.length.mm.) #sample size

mean_s_hp2 <- mean(h$Snout.length.mm.)  #mean of hatch snout for p2
length(h$Snout.length.mm.) #sample size

sd_s_wp2 <- sd(w$Snout.length.mm.)#standard deviations, for my table, for snout lenght
sd_s_hp2 <- sd(h$Snout.length.mm.)



#results- QC please
df_pink2_odd_snout_results <- data.frame(t=-1*sum_snout_p2$coefficients[3,3], df= sum_snout_p2$df[2], p_one_sided=sum_snout_p2$coefficients[3,4]/2, #not 1-pval/2 because hatch is larger than wild. And it is divided by 2 to make it a one-sided t-test
                                         length_t = sum_snout_p2$coefficients[2,3], length_p_two_Sided = sum_snout_p2$coefficients[2,4], 
                                         hatchery_mean= mean_s_hp2, hatchery_sd = sd_s_hp2,
                                         wild_mean=mean_s_wp2, wild_sd=sd_s_wp2)

df_pink2_odd_depth_results <- data.frame(t=-1*sum_depth_p2$coefficients[3,3], df= sum_depth_p2$df[2], p_one_sided=sum_depth_p2$coefficients[3,4]/2, #not 1-pval/2 because hatch is larger than wild. And it is divided by 2 to make it a one-sided t-test
                                    length_t = sum_depth_p2$coefficients[2,3], length_p_two_Sided = sum_depth_p2$coefficients[2,4], 
                                    hatchery_mean= mean_d_hp2, hatchery_sd = sd_d_hp2,
                                    wild_mean=mean_d_wp2, wild_sd=sd_d_wp2)

#RESULTS RESULTS

SNOUT_RESULTS <- rbind(df_pink_snout_results, df_pink2_odd_snout_results, df_coho_snout_results)
SNOUT_RESULTS <-  data.frame(SNOUT_RESULTS, row.names=c("Pink 2020", "Pink 2021", "Coho"))

DEPTH_RESULTS <- rbind(df_pink_depth_results, df_pink2_odd_depth_results, df_coho_depth_results )
DEPTH_RESULTS <-  data.frame(DEPTH_RESULTS, row.names=c("Pink 2020", "Pink 2021", "Coho"))

#write.csv(SNOUT_RESULTS, "Results/morpho_snout_results.csv")
#write.csv(DEPTH_RESULTS, "Results/morpho_depth_results.csv")

#############################################################################################################
##############################################################################################################
#investigate date - line 275

p2.males_dateadj <- read.csv("Data/p2.males2_dateadj.csv") #line 721

#GRAPHS WITH DATE
ggplot(p2.males_dateadj) + aes(x=Date_adj, y=Body.depth.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))
ggplot(p2.males_dateadj) + aes(x=Date_adj, y=Snout.length.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))

#MODELS WITH DATE
#Snout_gl#GRAPHS WITH DATE
ggplot(p2.males_dateadj) + aes(x=Date_adj, y=Body.depth.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))
ggplot(p2.males_dateadj) + aes(x=Date_adj, y=Snout.length.mm.) + geom_boxplot() + geom_jitter(aes(color=Otolith.reading, shape=Location))

#MODELS WITH DATE
Snout_global_pinkodd <- lm(Snout.length.mm. ~ Length.mm. + Otolith.reading +  Location + Date_adj, data=p2.males_dateadj)
summary(Snout_global_pinkodd)


Depth_global_pinkodd <- lm(Body.depth.mm. ~ Length.mm. + Otolith.reading +  Location + Date_adj, data=p2.males_dateadj)
summary(Depth_global_pinkodd) #date not sig to this model, but lcoaiton is

#lubridate
library(lubridate)
Julian_pinkodd <- julian(mdy(p2.males_dateadj$Date), origin=as.Date("2021-01-01") )
Julian_pinkodd
p2.males_dateadj$Julian <- Julian_pinkodd

#nice. now test for Julian date:
Snout_global_pinkodd_julian <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
summary(Snout_global_pinkodd_julian) #ok well, that's important. DAte IS signficant here, barely

Depth_global_pinkodd_julian <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
summary(Depth_global_pinkodd_julian) 

Depth_global_pinkodd_julian_check <- lm(log(Body.depth.mm.) ~ Length.mm. + Julian + Otolith.reading, data=p2.males_dateadj)
summary(Depth_global_pinkodd_julian_check)

Depth_global_pinkodd_julian_global <- lm(log(Body.depth.mm.) ~ Length.mm. * Julian * Otolith.reading, data=p2.males_dateadj)
summary(Depth_global_pinkodd_julian_global)

names(p2.males_dateadj)
#View(p2.males_dateadj)

##################################3
#05/23/24 
##I need AIC tables 
Snout_global_pinkodd_julian <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
Snout_global_pinkodd_julian_glob <- lm(log(Snout.length.mm.) ~ Length.mm. * Otolith.reading * Julian, data=p2.males_dateadj)
S_2 <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Length.mm.:Julian + Otolith.reading:Julian, data=p2.males_dateadj)
S_3 <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Otolith.reading:Julian, data=p2.males_dateadj)
S_4 <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian + Otolith.reading:Julian, data=p2.males_dateadj)
S_5 <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
S_6 <- lm(log(Snout.length.mm.) ~ Length.mm. * Otolith.reading, data=p2.males_dateadj)
S_7 <- lm(log(Snout.length.mm.) ~ Otolith.reading, data=p2.males_dateadj)

B_glob <- lm(log(Body.depth.mm.) ~ Length.mm. * Otolith.reading * Julian, data=p2.males_dateadj)
B_2 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Length.mm.:Julian + Otolith.reading:Julian, data=p2.males_dateadj)
B_3 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Otolith.reading:Julian, data=p2.males_dateadj)
B_4 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Otolith.reading:Julian, data=p2.males_dateadj)
B_5 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
B_6 <- lm(log(Body.depth.mm.) ~ Length.mm. * Otolith.reading, data=p2.males_dateadj)
B_7 <- lm(log(Body.depth.mm.) ~ Otolith.reading, data=p2.males_dateadj)


AIC(B_glob, B_2, B_3, B_4, B_5, B_6, B_7)

#08/07/24
#summary for my tables supp from hell
summary(B_6) #global mod #right direction for signs and things??
sum_depth_p2 #TK is the sign in the right direction??
DEPTH_RESULTS

summary(S_6)
sum_snout_p2
SNOUT_RESULTS
################################################


#test date, early vs. late by location? Sashin, lovers, armstrong? 12/11/23
##Graph thus
#EXPLORATORY DATE PLOTS - these show not much diff from one to the next
ggplot(p2.males_dateadj) + aes(x=Julian, y=Body.depth.mm.) +geom_point(aes(color=Location, shape=Otolith.reading)) +
  geom_smooth()
ggplot(p2.males_dateadj) + aes(x=Julian, y=Snout.length.mm.) +geom_point(aes(color=Location, shape=Otolith.reading)) +
  geom_smooth()


##################################################################################################################
####################################################################################################################


#snaky fish tests -turns out they don't make a difference

#might be needed later: write.csv(p2.males1, "p2.males1.csv")


################################################################################3
################################################################################
#########   Male fish graphs: morpho and linear   #########################################
############    (and length too??)             ####################################
#################################################################################
#################################################################################


##########
#pink 2020 - spending my dec 2023 time unlogging things
##########

df.pink2$snoutmm <- df.pink2$snout*10

#range(log(df.pink2$snout)) #1.18, 1.91
#range(log(df.pink2$depth)) #4.6, 5.2
range(df.pink2$snoutmm) #32.6, 67.4
range(df.pink2$depth) #101, 175
range(df.pink2$length) #371, 515

#snout
ggSnout_p1 <-ggplot(df.pink2) +geom_point(size=2, aes(y=snoutmm, x=length, color=Wild.or.hatch) )+
  geom_smooth(method="lm", aes(y=snoutmm, x=length, color=Wild.or.hatch), color="black")+ #keep or get rid of this line?
  scale_color_manual(breaks =c("W","H"), values=c("blue", "orange"), name=element_blank(), labels=c("Wild origin", "Hatchery origin")) + theme_cowplot() + labs(x="Length (mm)", y= "Snout (mm)")+
  guides(color="none") + 
  coord_cartesian(xlim=c(367, 520), ylim=c(30, 70))+ #changed ylim to account for non-logging
  scale_x_continuous(breaks=c(400, 440, 480, 520), expand=c(0,0)) + 
  scale_y_continuous(breaks= c(30, 50, 70), expand=c(0,0)) #changed breaks to account for not logging . MAy want to change
#looks good - AR 12/11/23

#depth
ggDepth_p1 <- ggplot(df.pink2) + geom_point(size=2, aes(y=depth, x=length, color=Wild.or.hatch) ) +
  geom_smooth(aes(y=depth, x=length), method="lm", color="black")+ 
  scale_color_manual(breaks =c("W","H"), values=c("blue", "orange"), name=element_blank(), labels=c("Wild origin", "Hatchery origin")) + theme_cowplot() + labs(x="Length (mm)", y= "Depth (mm)")+
  guides(color="none")+
  coord_cartesian(xlim = c(367, 520), ylim=c(95, 180))+
  scale_x_continuous(breaks=c(400, 440, 480, 520), expand=c(0,0)) +
  scale_y_continuous(breaks=c(100, 140, 180), expand=c(0,0)) #ok?
#looks good-AR 12/11/23



##############
#pink 2021
###############
range(p2.males1$Length.mm.) #345, 492
#range(log(p2.males1$Body.depth.mm.)) #4.2, 5.0
#range(log(p2.males1$Snout.length.mm.)) #3.6, 4.3
range(p2.males1$Body.depth.mm.) #67, 143
range(p2.males1$Snout.length.mm.) #37 76

#p2 snout 
ggsnout_pinkodd <- ggplot(p2.males1) + aes(x=Length.mm., y=Snout.length.mm., color=Otolith.reading) + 
  scale_color_manual(values=c("blue","orange")) + theme_cowplot()+
  guides(color= "none") +
  geom_point(size=2) +
  geom_smooth(method = "lm") + 
  labs(x= "Length (mm)", y= "Snout (mm)") +
  coord_cartesian(xlim = c(340, 500), ylim=c(35, 80)) +
  scale_x_continuous(breaks=c(350, 400, 450, 500), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(40, 60, 80), expand=c(0,0))

#p2 depth
ggdepth_pinkodd<- ggplot(p2.males1) + aes(x=Length.mm., y=Body.depth.mm., color=Otolith.reading) +   scale_color_manual(values=c("blue","orange")) + theme_cowplot()+
  guides(color= "none")+
  geom_point(size=2) +
  geom_smooth(method = "lm") + 
  labs(x= "Length (mm)", y= "Depth (mm)") + 
  coord_cartesian(xlim = c(340, 500), ylim=c(65, 150)) +
  scale_x_continuous(breaks=c(350, 400, 450, 500), expand=c(0,0)) + 
  scale_y_continuous(breaks=c(75, 100, 125, 150), expand=c(0,0))

###################3
#coho 
#################
df.coho$snoutmm <- (df.coho$snoutL)*10

range(df.coho$length) #416, 631 #updated length
range(df.coho$snoutmm) #79, 194 #range stays same
range(df.coho$depth) #104, 192 #updated depth
mycolors.coho <- c("orange", "blue")

#snout
ggsnout_coho <- ggplot(data=df.coho) + aes(x=length, y=snoutmm, color=Wild.or.hatch) + geom_point(size=2) +
  geom_smooth(method="lm")+  
  scale_color_manual(values=mycolors.coho) + 
  theme_cowplot() + labs(y="Snout (mm)", x= "Length (mm)") + #changed snout to element blank
  guides(color="none") + 
  coord_cartesian(xlim=c(410, 650), ylim=c(75, 200)) + 
  scale_x_continuous(expand=c(0,0), breaks = c(450, 550, 650)) +
  scale_y_continuous(expand=c(0,0), breaks =c (100, 150, 200))


ggdepth_coho <-ggplot(df.coho)+aes(y=depth, x=length, color=Wild.or.hatch) + 
  geom_point(size=2) +
  geom_smooth(method="lm", key_glyph= "blank") + scale_color_manual(values=mycolors.coho, labels=c("Hatchery", "Wild")) + 
  theme_cowplot()+labs(y="Depth (mm)", x= "Length (mm)") +
  #guides(color="none") + 
  coord_cartesian(xlim=c(410, 650), ylim=c(100, 200)) + 
  scale_x_continuous(expand=c(0,0), breaks = c(450, 550, 650)) +
  scale_y_continuous(expand=c(0,0), breaks =c (100, 150, 200)) +
  theme(legend.title=element_blank(), legend.position = c(0.55, 0.12))

ggdepth_coho <- ggdepth_coho + theme(
  legend.box.background = element_rect(),
  legend.box.margin = margin(0, 3, 2, 2) #yes, finally I made a box.
)

###########################################
#combine the linear morpho graphs
##01/29/24 - I'm adding labels here.
###############################################
male_base <- plot_grid(ggSnout_p1, ggsnout_pinkodd, ggsnout_coho, ggDepth_p1, ggdepth_pinkodd, ggdepth_coho, nrow=2, ncol=3, scale=0.95) #megan suggests 2 x 3

male_base2 <- plot_grid(NULL, male_base, ncol = 1, rel_heights = c(0.6,9.4))


##01/29/24 insert below
###lol I guess this is code from before I discovered patchwork
#library(patchwork)
##xlab_GSI <- ggSnout_p1$labels$x

plot_male_lables <- ggdraw(male_base2) +
  draw_label ("Even-year pink", x = 0.18, y = 0.95, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Odd-year pink", x = 0.52, y = 0.95, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Coho", fontfamily = "Arial", fontface="bold", size = 13, x = 0.86, y = 0.95) 
plot_male_lables 


dev.new (width = 10, height =6.56, unit = "in", noRStudioGD = T); last_plot()
#dev.off()
ggsave ("Plots/Male_linear_labs_jan.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()
#nice!

#end 01/29/24

#argh. not quite right
##try patchwork. Jan comment: I prefer the graph above, plot_male_labels
(ggSnout_p1 + ggsnout_pinkodd + ggsnout_coho)/(ggDepth_p1 +  ggdepth_pinkodd + ggdepth_coho)
dev.new (width = 10, height =6.56, unit = "in", noRStudioGD = T); last_plot()
dev.off()
#ggsave ("Plots/Male_linear.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off() #one of these
##hastagged out because the output isnt consistent.

###########################################################
##########################################################
#male morpho graphs
#############################################################
############################################################

#pink even
#even pink (p1)
ggMorpho_pinkeven <- ggplot(data=df.pink2)+aes(x=Comp3, y=Comp1, color=Wild.or.hatch)+
  geom_point(size=2)+scale_color_manual(values=wildhatch) + stat_ellipse() + 
  theme_bw() + labs(x= "RW 3", y= "RW 1") +
  guides(color="none") +
  labs(x=NULL, y=NULL)
#RW 3 is hump. RW 1 is hump and snout. RW 2 is bendyness(not having to do with fish morphometrics, just having to do with how the fish was placed for photo)


#pink odd
#odd pink (p2)
ggMorpho_pinkodd<- ggplot(data=df.pink2021.orig)+
  aes(x=Comp3, y=Comp1, color=origin)+
  geom_point(size=2)+scale_color_manual(values=wildhatch) + 
  stat_ellipse() + theme_bw() + 
  labs(x= "RW 3", y= "RW 1")+
  guides(color="none")+
  labs(x=NULL, y=NULL)

#coho (all)
ggMorpho_coho <- ggplot(df.coho) + aes(x=Comp3, y=Comp1, color=Wild.or.hatch) + #changed from df.coho.long  to df.coho
  geom_point(size=2)+ scale_color_manual(values=mycolors.coho) + 
  stat_ellipse() + 
  labs(y= "RW 1", x= "RW 3") + 
  theme_bw()+
  guides(color="none")+
  labs(x=NULL, y=NULL)
#"RW 1 (snout, roughly)", "RW 3 (depth, roughly)"


#combine all 3
#ggMorpho_pinkeven + ggMorpho_pinkodd + ggMorpho_coho

morpho_base <- plot_grid(ggMorpho_pinkeven, ggMorpho_pinkodd, ggMorpho_coho, nrow=1, ncol=3, scale=0.95)
morpho_base2 <- plot_grid(NULL, morpho_base, ncol = 2, rel_widths= c(0.3,9.7))
morpho_base3 <- plot_grid(NULL, morpho_base2, NULL, nrow=3, rel_heights = c(0.5, 9.0, 0.5))

#now axes and labels... 
plot_morpho_male <- ggdraw(morpho_base3) +
  draw_label("Even-year pink", x = 0.21, y = 0.965, fontfamily = "Arial", fontface="bold", size = 15)+
  draw_label("Odd-year pink", x=0.53, y=0.965, fontfamily = "Arial", fontface="bold", size = 15)+
  draw_label("Coho", x=0.85, y=0.965, fontfamily = "Arial", fontface="bold", size = 15) +
  draw_label ("RW3", x = 0.56, y = 0.05, size = 15) + 
  draw_label (("RW1"), angle= 90, x = 0.03, y = 0.50, size = 15)
plot_morpho_male 

dev.new (width = 10, height =3.5, unit = "in", noRStudioGD = T); last_plot()
#ggsave ("Plots/Male_morpho_toppart.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()
dev.off()

#I'll need to combine this morpho graph with teh fish graph, which I did on ppt. Get the fish graph

par(mfrow = c(1, 3), mar=c(1,1,1,0), family  = "Arial")
plotRefToTarget(ref.w.p, ref.h.p, method="vector", mag=5, mar=c(1,1,1,0))
plotRefToTarget(ref.w.p2, ref.h.p2, method="vector", mag=5, mar=c(1,1,1,0)) #wild dot, hatchery arrows
plotRefToTarget(ref.w.c, ref.h.c, method="vector", mag=5, mar=c(1,1,1,0))
dev.new (width = 10, height =3.5, unit = "in", noRStudioGD = T); last_plot()
#ggsave ("Plots/Male_morpho_bottompart.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()
##hmm. Had to save last plot manually
dev.off()




##########################################################################################
##########################################################################################

#post-hoc power analysis

###############################################################################################
#############################################################################################


library(pwr)

###########################################################################
#MALE 
##these were general linear models with R, so use: pwr.f2.test for all snout and depth for males
?pwr.f2.test
#############################################################################
#P1#############################################
##dataframe df.pink2
#snout
##should u be one? u=1 in all of Charlie power analyses
##v: df based on sample size. sample size-2?
###what's the sample size here?? 59. So my df is... 57? v=57? Check results tables for df
##df for this model (sum_p_snout) os 2 and 56. Idk if to use 1 and 57 or two and 56.
length(df.pink2$snoutmm)
##f2 is the effect size. Where's my model?
sum_p_snout #factor(Wild.or.hatch)W  0.127469

#the test
pwr_p1_s <- pwr.f2.test(u=1, v=57, f2=0.127469, sig.level=0.05, power=NULL) #power is 0.77. That means, it can detect a medium size?
##an interpretation? This means that, under the specified conditions, the statistical power of your test to detect the given 
##effect size is estimated to be around 76.91%. In other words, if the true effect size is 0.127469 (as specified by f2), and 
##the test is conducted with the specified degrees of freedom and significance level, you have a 76.91% chance of correctly rejecting 
##the null hypothesis (assuming it is false).


#depth
sum_p_depth #effect size for wild or hatch: 1.80562
pwr_p1_d<-pwr.f2.test(u=1, v=57, f2=1.80562, sig.level=0.05, power=NULL) #power is 1.
##so fish origin (Wild or hatch) definitely has an impact on depth?

#P2#####################################
##dataframe: p2.males1
length(p2.males1$Length.mm.) #sample size is 98

#snout

sum_depth_p2 #factor(Otolith.reading)PORT ARMSTRONG -10.78867  #f2 must be postive so... just make it positive?
pwr_p2_d <- pwr.f2.test(u=1, v=98, f2=10.78867, sig.level=0.05, power=NULL) #power is 1. 

#depth

sum_snout_p2 #factor(Otolith.reading)PORT ARMSTRONG  -2.59889
pwr_p2_s <- pwr.f2.test(u=1, v=98, f2=2.59889, sig.level=0.05, power=NULL) #power is 1.

#C
##dataframe: df.coho
length(df.coho$snoutmm) #length is 60

#snout
sum_c_5 #factor(Wild.or.hatch)W -1.325819
pwr_c_s <- pwr.f2.test(u=1, v=60, f2=1.325819, sig.level=0.05, power=NULL) #power is 1

#depth 
sum_c_depth #factor(Wild.or.hatch)W  -5.61080
pwr_c_d <- pwr.f2.test(u=1, v=60, f2=5.61080, sig.level=0.05, power=NULL) #power is 1

######
#create csv
power_males <- data.frame(snout_power=c(pwr_p1_s$power, pwr_p2_s$power ,pwr_c_s$power), depth_power=c(pwr_p1_d$power, pwr_p2_d$power, pwr_c_d$power))
  
rownames(power_males) = c("Pink even", "Pink odd", "Coho")

write.csv(power_males, "Results/Male post-hoc power analysis.csv")


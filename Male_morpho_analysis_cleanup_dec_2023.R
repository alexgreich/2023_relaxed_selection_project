#######################################
###male morpho analysis cleanup - dec 2023
####################################

#year one data will be first, then year two.

#year one (pink 2020, coho) is copy-pasted from "Morpho_yr1_pinkcoho.Rmd" 
#and year two (pink 2021) is copy-pasted from "Morpho_yr2_Rmd"
##so if any issues, look there
#I may also add on the simply graphing.R file at the end, if it makes sense


#################################################3
#load packages
library(ggplot2)
library(geomorph)
library(ggfortify)
library(dplyr)
library(tidyr)
library(ggforce)
library(Morpho)
library(RColorBrewer)
library(cowplot)
library(patchwork)
######################################################

##########################
########################
#pink 2020
#########################
########################

#################################
#next, load TPS file
specID <- c("HM3", "HM19", "WM7", "HM13", "WM5", "HM29", "HM11", "WM15", "HM14", "WM16", 
            "HM30", "HM7", "WM25", "HM2", "HM22", "HM8", "HM16", "WM19", "WM28", "WM23",
            "WM9", "HM18", "WM29", "HM25", "HM1", "WM14", "WM2", "HM12", "WM6",
            "WM13", "WM18", "HM9", "WM10", "WM1", "HM27", "WM22", "HM20", "WM30", "HM5",
            "HM21", "HM24", "WM21", "HM17", "WM26", "HM4", "WM11", "WM3", "HM28", "WM24",
            "HM23", "HM6", "WM20", "WM8", "WM27", "HM26", "WM4", "HM10", "WM17", "HM15", "WM12")
#specID is the Id.
#QC comments: this looks good

pinkdat.lands <- readland.tps("Data/pinks3.tps", readcurves=TRUE)
###########################################

################################################
#remove outlier
pinkdat.landsno56 <- pinkdat.lands
pinkdat.landsno56 <- pinkdat.landsno56[,,-55]
#pinkdat.landsnoHM26 <- pinkdat.lands[,,-56]
dim(pinkdat.landsno56)

plotAllSpecimens(pinkdat.lands)
plotAllSpecimens(pinkdat.landsno56)
#plotAllSpecimens(pinkdat.landsnoHM26)
######################################


#############################################data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAkCAYAAAD7PHgWAAAA00lEQVR42mNgGAWjYBQMUxAauorZLL452TyhZQUtMMhs47QGLrIdaJ7QmtSyeP+5fTc//N98+e1/agGQWSvOvPqfNGHnRbO4lnjyHRjfvHzvzff/Zx5+/r9x60OqORBkFgg3bHnw1yy+ZQkFIdiyAuRbmIHUdiAIg+wYdeCoA0cdOOrAUQdSyYG0wKMOHHUgOQ6kNGOMOhCXpaMOHHXgiHTgSmDva9A6ENRvTejfcYFWDkzs33kBZAfZDvTMncQO6huDup+06rhbhvZxjg6RjILBDgAZYqbbTdtPRgAAAABJRU5ErkJggg==
#this is the code that was used to create the sliding semi-landmarks. It was needed to create the file for semilandmarks
##and does not need to be run again.

#pinks
#define.sliders(pinkdat.landsno56[,,1], nsliders=3)
#needs to be done not in a chuk

#coho
#define.sliders(coho_lands.no_ad[,,1], nsliders=3)

########################################### #line 74


##################################################################
#now transform into Procrustes shape space

#pinkcurves <- as.matrix(read.csv("/Users/alexandrareich/Desktop/THESIS TIME!/Compilation of GITHUB project code/Relaxed_selection_proj_reproducible_science/DATA/curveslide.csv", header=T))
#getwd()
#library(here)
pinkcurves <- as.matrix(read.csv("Data/curveslide.csv", header=T))

gpa.no56 <- gpagen(pinkdat.landsno56, curves = pinkcurves)
plot(gpa.no56)

gpa_withoutlier <-  gpagen(pinkdat.lands, curves = pinkcurves)


plot(gpa_withoutlier)

##############################################################


####################################################################### #line 95 of Morpho_yr1_pinkcoho
#Now run a PCA/RWA (principal component analysis aka Relative Warp Analysis)
pca.no56 <- gm.prcomp(gpa.no56$coords)
plot(pca.no56)
plot(pca.no56, axis1=1, axis2=3)
#how to add categorical variables though?
#wait, I did that.

pca_withoutlier <- gm.prcomp(gpa_withoutlier$coords)

plot(pca_withoutlier)
###########################################################################


###########################################3
#convert to dataframe

summary(pca.no56)
#convert array to dataframe
#do I need to convert to dataframe?
df.pink2 <- as.data.frame(pca.no56$x[,1:3])
class(df.pink2) 


#ok, need to add river column based on origin
#probaly an easier way to do this, but:
#trying to make a vector of length 180 first...
#or make smaller vectors that fit the river, and then combine, and then add to dataframe

#create wild/hatch vector based on specID
Wild.or.hatch2 <- c("H", "H", "W", "H", "W", "H", "H", "W", "H", "W", 
                    "H", "H", "W", "H", "H", "H", "H", "W", "W", "W",
                    "W", "H","W", "H", "H", "W", "W", "H", "W",
                    "W", "W", "H", "W", "W", "H", "W", "H", "W", "H",
                    "H", "H", "W", "H", "W", "H", "W", "W", "H", "W",
                    "H", "H", "W", "W", "W", "W", "H", "W", "H", "W")

#and add SpecID into there somewhere?
specID2 <- specID[-55]

head(df.pink2)
df.pink2$Wild.or.hatch <- Wild.or.hatch2
df.pink2$ID <- specID2 #GO BACK AND LOOK AT FISH PHOTOTS
head(df.pink2)

##############################################################################

###############################################################################
#exploratory plot

wildhatch <- c("orange", "blue")

(plotp<-ggplot(data=df.pink2)+aes(x=Comp1, y=Comp2, color=Wild.or.hatch)+geom_point()+scale_color_manual(values=wildhatch)+
    stat_ellipse() + theme_cowplot())

#10/12/21 edit, tinkering with graph

ggplot(data=df.pink2)+aes(x=Comp3, y=Comp1, color=Wild.or.hatch)+geom_point()+scale_color_manual(values=wildhatch) + stat_ellipse() + theme_bw() + labs(x= "RW 3 (hump)", y= "RW 1 (hump and snout)")


##################################################################################


##################################################################################
#OK, NOW fix this by placing wild or hatch category within the Procrustes array
#My long term goal: get groups... gpagen groups?
 # To color code wild/hatch

summary(gpa.no56)
summary(pinkdat.landsno56)

##################################################################################


####################################################################################
#data analysis
##coresponds to chunk on line 220 of Morpho_yr1_pinkcoho.Rmd code
Y.gpa <- gpa.no56
gdf <- geomorph.data.frame(Y.gpa, origin = Wild.or.hatch2)

fit.p <- procD.lm(coords ~ log(Csize) + origin, data = gdf,  #all coords taken into account
                  iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.p)  #NOT SIGNIFICANT...well origin is not sig, but csize is...

#are interaction effects a thing?
fit.p.int <- procD.lm(coords ~ log(Csize) + origin + log(Csize)*origin, data = gdf,  #all coords taken into account
                      iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.p.int)  #NOT SIGNIFICANT, and interaction effects are not significant

#ok. got the fit . now what? test it? plot it?...I think that's it. Maybe just plot the relative warps too.

#log of csize is significant, but not origin categorical variabble
#A MANCOVA


#making the results csv, edits 12/7/23 for increased reproducibility
sum_p <- summary(fit.p)
F_origin_p1 <- 1.6659 #from fit.p 
F_centroidsize_p1 <- 12.4231 #log of centroid size F
p_val_origin <- 0.1108
p_val_centroid_size <- 1e-04

names_p1_geomorph <- c("F origin", "F log centroid size", "P value origin", "P value centroid size")
results_p1_geomorph <- c(F_origin_p1,F_centroidsize_p1, p_val_origin, p_val_centroid_size)

p1_geomorph <- cbind(names_p1_geomorph, results_p1_geomorph)

#instead of saving what I did above, I'm just going to save the output table from sum_p as my results csv, dec 23

write.csv(sum_p$table, "Results/morpho_p1_table.csv")

#write.csv(p1_geomorph, "Even-year pink geomorph results.csv")
###############################################################################################

#################################################################################################
#more (morpho) plots for pinks
#these plots were important for seeing what the difference in relative warp was
plotp

#plotReftoTarget
#pinks are in Y.gpa
ref <- mshape(Y.gpa$coords)
w.p <- numeric()
h.p <- numeric()
i=1

for (i in 1:length(Wild.or.hatch2)){
  if(Wild.or.hatch2[i]=="W"){
    w.p <- c(w.p, i)
  }else{
    h.p <- c(h.p, i)
  }
}

ref.w.p <-mshape(Y.gpa$coords[,,w.p]) #wild mean
ref.h.p <-mshape(Y.gpa$coords[,,h.p]) #hatchery mean

#wild to hatch, hatch to wild examine
plotRefToTarget(ref.w.p, ref.h.p, method="vector", mag=5)
plotRefToTarget(ref.h.p, ref.w.p, method="vector", mag=5)


#?mshape()
#plotRefToTarget()

#pca plot ref to target
ref <- mshape(gpa.no56$coords)
ref_p1 <- mshape(gpa.no56$coords)

plotRefToTarget(pca.no56$shapes$shapes.comp1$min, ref, method="TPS")
plotRefToTarget(pca.no56$shapes$shapes.comp1$max, ref, method="TPS")

plotRefToTarget(pca.no56$shapes$shapes.comp2$min, ref, method="TPS")
plotRefToTarget(pca.no56$shapes$shapes.comp2$max, ref, method="TPS")


plotRefToTarget(pca.no56$shapes$shapes.comp1$min, ref, method="vector")
plotRefToTarget(pca.no56$shapes$shapes.comp1$max, ref, method="vector")

plotRefToTarget(pca.no56$shapes$shapes.comp2$min, ref, method="vector")
plotRefToTarget(pca.no56$shapes$shapes.comp2$max, ref, method="vector")

plotRefToTarget(pca.no56$shapes$shapes.comp3$min, ref, method="vector")
plotRefToTarget(pca.no56$shapes$shapes.comp3$max, ref, method="vector") 
#THIS IS IT! COMP 3 is the hump!!!

plotRefToTarget(pca.no56$shapes$shapes.comp4$min, ref, method="vector")
plotRefToTarget(pca.no56$shapes$shapes.comp4$max, ref, method="vector") 

#can you also plot the ID of the important PCs???
pvar <- (pca.no56$sdev^2)/(sum(pca.no56$sdev^2))
names(pvar) <- seq(1:length(pvar))
barplot(pvar, xlab= "Principal Components", ylab = "% Variance")

#also plot the CVA, please
vari <- gpa.no56$coords
facto <- Wild.or.hatch2

CVA.pink <- CVA(vari, groups=facto)
#############################################################################################

#####################################################################################################
######################################################################################################
#############################               COHO               #######################################
#####################################################################################################3
######################################################################################################
#corresponds to line 322 in Morpho_yr1_pinkcoho.Rmd code


####################################################################################################
#load in coho data

#coho_lands <- readland.tps("/Users/alexandrareich/Desktop/THESIS TIME!/Compilation of GITHUB project code/Relaxed_selection_proj_reproducible_science/DATA/Coho091021_rand_fixed062221.tps", readcurves=TRUE)
coho_lands <- readland.tps("Data/Coho091021_rand_fixed062221.tps", readcurves=TRUE)

coho_lands<-coho_lands[c(-5,-6),,]
#the photo names should tell you the ID, probably need to go make a vector of the coho ID's tho

#need to reorder my length vector to match the coho ID
#or reorder my coho ID to match the length vector....
the_order <- c(45,43,41, 49, 3, 16, 44, 13, 54, 47, 17, 55, 18, 21, 34, 57, 39, 2, 4, 53, 24, 28, 22, 31, 35, 15, 60,
               52, 56, 37, 33, 51, 20, 19, 29, 42, 27, 48, 30, 7, 26, 40, 59, 25, 11, 12, 14, 8, 38, 36, 32, 46, 5, 50,
               23, 9, 58, 1, 6, 10)
#WHAT IS length_c....
#######3
#Coho length/depth stuff
#coho.dat <- read.csv("/Users/alexandrareich/Desktop/THESIS TIME!/Compilation of GITHUB project code/Relaxed_selection_proj_reproducible_science/DATA/MASTERMaleCohoQC copy.csv")
coho.dat <- read.csv("Data/MASTERMaleCohoQC copy.csv")
names(coho.dat)
length_c <- coho.dat$Length..mm.

########3


length_ordered <- c(602,583,length_c[41], length_c[49], length_c[3], length_c[16], length_c[44], length_c[13], length_c[54], length_c[47], length_c[17], length_c[55], length_c[18], length_c[21], length_c[34], length_c[57], length_c[39], length_c[2], length_c[4], length_c[53], length_c[24], length_c[28], length_c[22], length_c[31], length_c[35], length_c[15], length_c[60],
                    length_c[52], length_c[56], length_c[37], length_c[33], length_c[51], length_c[20], length_c[19], length_c[29], length_c[42], length_c[27], length_c[48], length_c[30], length_c[7], length_c[26], length_c[40], length_c[59], length_c[25], length_c[11], length_c[12], length_c[14], length_c[8], length_c[38], length_c[36], length_c[32], length_c[46], length_c[5], length_c[50],
                    length_c[23], length_c[9], length_c[58], length_c[1], length_c[6], length_c[10])

#l.c <- numeric()
#i=1
#for (i in 1:length(length_ordered)){
 # if (length_ordered[i] > 525){
  #  l.c <- c(l.c,i)
  #}
#}

Wild.or.hatch.d <- c("H","H", "H", "H", "W", "W", "H", "W", "H", "H", #47
                     "W", "H", "W", "W", "H", "H", "H", "W", "W", "H",#53
                     "W", "W", "W", "H", "H", "W", "H", "H", "H", "H", #37
                     "H", "H", "W", "W", "W", "H", "W", "H", "W", "W", "W", #26
                     "H", "H", "W", "W", "W", "W", "W", "H", "H", "H", #32
                     "H", "W", "H", "W", "W", "H", "W", "W", "W")
#######################################################################################################


###########################################################################################################
#Plot all speciments. hashtagging out long specimens, vecause we are not subsetting in this version

#Wild.or.hatch.long<- Wild.or.hatch.d[l.c]
#coho_lands_long <- coho_lands[,,l.c]
plotAllSpecimens(coho_lands)
#plotAllSpecimens(coho_lands_long)
###############################################################################################################


###############################################################################################################3
#Add in the SEMILANDMARKS (special to this iteration of data)
#should be the same as the otehr matrix, so not adding anything new here - line 389
################################################################################################################


#################################################################################################################
#procrustes shape space transformation #line 393

gpa_coho <- gpagen(coho_lands, curves=pinkcurves) #should be the same matrix as pinkcurves, nothing special here. but be careful
#plotAllSpecimens(coho_lands)
plot(gpa_coho)

#gpa_coho_long <- gpagen(coho_lands_long, curves=pinkcurves) #pink sliders same as coho for id purposes, so OK to use here.
#plotAllSpecimens(coho_lands_long)
#plot(gpa_coho_long)
####################################################################################################################


#####################################################################################################################
#PCA/rel warps
#run a PCA
pca_coho <- gm.prcomp(gpa_coho$coords)
plot(pca_coho)
plot(pca_coho, axis1=1, axis2=3)
#summary(pca_coho)
#names(pca_coho)
#pca 1 : 53.14% of variance explained
#pca 2: 16.53
#pca 3: 6.81% var explained


#pca_coho_long <- gm.prcomp(gpa_coho_long$coords)
#plot(pca_coho_long)
#summary(pca_coho_long)
#names(pca_coho_long)

#can you also plot the ID of the important PCs???
pvar_cl <- (pca_coho$sdev^2)/(sum(pca_coho$sdev^2))
names(pvar_cl) <- seq(1:length(pvar_cl))
barplot(pvar_cl, xlab= "Principal Components", ylab = "% Variance")

######################################################################################################################


######################################################################################################################
#PCA ggplot

df.coho <- as.data.frame(pca_coho$x[,1:3])
df.coho$Wild.or.hatch <- Wild.or.hatch.d 

ggplot(df.coho) + aes(x=Comp1, y=Comp2, color=Wild.or.hatch) + geom_point()+
  scale_color_brewer(palette="Set2")+
  stat_ellipse()


#df.coho.long <- as.data.frame(pca_coho_long$x[,1:3])
#df.coho.long$Wild.or.hatch <- Wild.or.hatch.long

#ggplot(df.coho.long) + aes(x=Comp1, y=Comp2, color=Wild.or.hatch) + geom_point()+
 # scale_color_brewer(palette="Set2")+
#  stat_ellipse()

#ggplot(df.coho.long) + aes(x=Comp1, y=Comp3, color=Wild.or.hatch) + geom_point()+
 # scale_color_brewer(palette="Set2")+
#  stat_ellipse()

######################################################################
#need to do a RW1 vs RW 3 graph 10/12/21 for coho and coho long #SHIT #delete??
mycolors.coho=c("orange", "blue")
#coho long
#ggplot(df.coho.long) + aes(x=Comp3, y=Comp1, color=Wild.or.hatch) + geom_point(size=2)+ scale_color_manual(values=mycolors.coho) + stat_ellipse() + labs(y= "RW 1 (snout, roughly)", x= "RW 3 (depth, roughly)") + theme_bw()  #RW1 is roughtly hump.RW3depth is 

#COHO all
ggplot(df.coho) + aes(x=Comp3, y=Comp1, color=Wild.or.hatch) + geom_point(size=2)+ scale_color_manual(values=mycolors.coho) + stat_ellipse() + labs(y= "RW 1 (snout and hump, roughly)", x= "RW 3 (depth, roughly)") + theme_bw() 

##########################################################################################################################


##################################################################################################
#coho models mancova 461

Y.gpa.c <- gpa_coho
gdf.c <- geomorph.data.frame(Y.gpa.c, origin = Wild.or.hatch.d)

fit.c <- procD.lm(coords ~ log(Csize) + origin, data = gdf.c,  #all coords taken into account
                  iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.c)
#sig, yes

#INT EFFECTS CHECK!!
fit.c.int <- procD.lm(coords ~ log(Csize) + origin + log(Csize)*origin, data = gdf.c,  #all coords taken into account
                      iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.c.int) #ooh, some interaction effects. That's expected though

#12/07/23
#ok, here's my main results
sum_c <- summary(fit.c)
#write.csv(sum_c$table, "Results/morpho_c_table.csv")
#note that this model can be swapped for the one with interaction effects

#AIC(fit.c, fit.c.int) #takes forever





###############################################################################################


#################################################################################################
#more plots for coho (the shapey fish plots) line 494

#plotReftoTarget
#pinks are in Y.gpa
ref <- mshape(Y.gpa.c$coords)
w.c <- numeric()
h.c <- numeric()
i=1

for (i in 1:length(Wild.or.hatch.d)){
  if(Wild.or.hatch.d[i]=="W"){
    w.c <- c(w.c, i)
  }else{
    h.c <- c(h.c, i)
  }
}

ref.w.c <-mshape(Y.gpa.c$coords[,,w.c]) #wild mean
ref.h.c <-mshape(Y.gpa.c$coords[,,h.c]) #hatchery mean

plotRefToTarget(ref.w.c, ref.h.c, method="vector", mag=5) #wild dots to hatchery arrows

#?mshape()
#plotRefToTarget()

#pca plot ref to target

plotRefToTarget(pca_coho$shapes$shapes.comp1$min, ref, method="TPS")
plotRefToTarget(pca_coho$shapes$shapes.comp1$max, ref, method="TPS")

plotRefToTarget(pca_coho$shapes$shapes.comp2$min, ref, method="TPS")
plotRefToTarget(pca_coho$shapes$shapes.comp2$max, ref, method="TPS")


plotRefToTarget(pca_coho$shapes$shapes.comp1$min, ref, method="vector")
plotRefToTarget(pca_coho$shapes$shapes.comp1$max, ref, method="vector")

plotRefToTarget(pca_coho$shapes$shapes.comp2$min, ref, method="vector")
plotRefToTarget(pca_coho$shapes$shapes.comp2$max, ref, method="vector")

plotRefToTarget(pca_coho$shapes$shapes.comp3$min, ref, method="vector")
plotRefToTarget(pca_coho$shapes$shapes.comp3$max, ref, method="vector") 
#THIS IS IT! COMP 3 is the hump!!!

plotRefToTarget(pca_coho$shapes$shapes.comp4$min, ref, method="vector")
plotRefToTarget(pca_coho$shapes$shapes.comp4$max, ref, method="vector") 
#LOOKS LIKE THE DORSAL.....

#can you also plot the ID of the important PCs???
pvar <- (pca_coho$sdev^2)/(sum(pca_coho$sdev^2))
names(pvar) <- seq(1:length(pvar))
barplot(pvar, xlab= "Principal Components", ylab = "% Variance")

#also plot the CVA, please
vari <- gpa_coho$coords
facto <- Wild.or.hatch.d

CVA.pink <- CVA(vari, groups=facto)

#10/12/21
#RW1 vs. RW3 PCA graph
##################################################################################################



################################################################################################################################
###############################################################################################################################
############################################     Yr1 linear morph analysis (p1, coho; males) line 743     ###############################
##################################################################################################################################
#####################################################################################################################################


###########################################
#######snout length using geomorph
#?interlmkdist
#landmakrs of interest for the coho:
#14 is the eye - #2 is end of snout
#1 is bottom of snout -  #3 is above nare
#interlmkdist(coho_lands) #using the un-tranformed data
lmks <- data.frame(snoutL=c(2,14), snoutH=c(1,3))
snoutlineardists <- interlmkdist(coho_lands, lmks)
#snoutlineardists
snout_df <- as.data.frame(snoutlineardists)
df.coho$snoutL <- snout_df$snoutL
df.coho$length <- length_ordered
df.coho$ID <- the_order #cleaning up my dataframe

#depth
depth_coho_ord <- numeric(60)
depth_c <- coho.dat$Body.Depth.mm.

for (i in (1:length(the_order))){
  depth_coho_ord[i] <- depth_c[the_order[i]]
}

df.coho$depth <- depth_coho_ord



#NOW linear regression on that
#snout_df$Wild.or.hatch <- Wild.or.hatch.long #woohoo, right #'s, no doubles
#need snout height and snout length as covariates, wild or hatch as a factor
#or, just use snout lenght and wild or hatch, one sided t-test


#revisit 10/25/21
#t.test(x=hatch_SL, y=wild_SL, alternative= "two.sided")

#revisit 11/23/21
##want to check error assumptions (especially constant variance)

#05/04/22

#1-0.002575/2

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

################################################


###############################
#DEPTH COHO #line 944

#revisit 10/25/21
#t.test(x=H2$depth, y=W2$depth, alternative= "two.sided")

###########3
#revisit 11/23/21 #what was the df for coho.all?
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


#see line 1007 in Morpho_yr1_pinkcoho.Rmd code for possible alternative: depth^3


#revisit 05/04/22 : to make ANVOVAS where needed to incorporate length (and..other stuff too?)
##lets do an ANCOVA test, one -way t test

aov.depth.c <- lm(depth ~ length + factor(Wild.or.hatch), df.coho)
sum_c_depth <- summary(aov.depth.c)
#plot(aov.depth.c) #looks fine ish


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





#PINK 1 SNOUTS AND DEPTHS


#line 1020 goes into pink snouts

########################################
#pink snouts and depth
lmks.p <- data.frame(snoutL=c(1,12))
snoutlineardists.p <- interlmkdist(pinkdat.landsno56, lmks.p)
#snoutlineardists
snout_df.no56 <- as.data.frame(snoutlineardists.p)
#ah! so easy!
#snout_df.no56 <- as.data.frame(snout_df.p[-56,]) #why dont delete this one?
#depth



df.pink2$snout <- snout_df.no56$snoutL
df.pink2


#wild_SL.p <- wild_S.p$`snout_df.p[-56, ]`
#hatch_SL.p<- hatch_S.p$`snout_df.p[-56, ]`



#LENGTHS PINK
#get the lengths for pink
#p.dat <- read.csv("/Users/alexandrareich/Desktop/THESIS TIME!/Compilation of GITHUB project code/Relaxed_selection_proj_reproducible_science/DATA/MASTERMalePinksQC.reorder.csv")
p.dat <- read.csv("Data/MASTERMalePinksQC.reorder.csv")
#names(p.dat)
length_p <- p.dat$Length..mm.

order_p <- c(33, 49, 7, 43, 5, 59, 41, 15, 44, 16, 
             60, 37, 25, 32, 52, 38, 46, 19, 28, 23,
             9, 48, 29, 55, 31, 14, 2, 42, 6,
             13, 18, 39, 10, 1, 57, 22, 50, 30, 35,
             51, 54, 21, 47, 26, 34, 11, 3, 58, 24,
             53, 36, 20, 8, 27, 56, 4, 40, 17, 45, 12)
anyDuplicated(order_p) #yay, no duplicates!


length_p_ord <- numeric(60)

i=1
for (i in (1:length(order_p))){
  length_p_ord[i] <- length_p[order_p[i]]
}
#YES!!!
length_no_56 <- length_p_ord[-55]

df.pink2$length <- length_no_56
#and need to restructure PC data to fit in ggplot
#Comp1_p <- as.data.frame(pca_coho$x[,1])
#names(Comp1_c)

#order depth
depth_p <- p.dat$Body.Depth.mm.

depth_p_ord <- numeric(60)

i=1
for (i in (1:length(order_p))){
  depth_p_ord[i] <- depth_p[order_p[i]]
}
#YES!!!
depth_no_56 <- depth_p_ord[-55]

df.pink2$depth <- depth_no_56

write.csv(df.pink2, "df.pink2.csv") #do I need this?



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
base_d <- lm(depth~ factor(Wild.or.hatch), data=df.pink2)

AIC(mod.p.d.l.int, p1_depth_mod, base_d ) #depth
AIC(mod.p.s.l.int,p1_snout_mod, base_s) #snout

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


###############################################################
#write csv's for the simply graphing R script (obsolete?)
#write.csv(df.coho.long2, "df.coho.long2.csv")
write.csv(df.coho, "df.coho.csv")
####################################################################




################################################################################################################################
##################################################################################################################################3
###################################################Morpho yr 2 pinks aka. pink2 aka pink odd##########################################
#####################################################################################################################################
#######################################################################################################################################
#corresponds to R file Morpho_yr2_pinks.Rmd
##I tried digital unbending (both cubic and quadratic) in the above file. I'm just using basic here. Not unbending was fine.


###################################################################################################################################3
#load tps file line 35

p2.lands <- readland.tps("Data/pink2021_landmarks.fixed73.deleteextras.appendhere_010522.tps", readcurves=T, specID= "imageID")


#plotall
plotAllSpecimens(p2.lands)


#on semilandmarks: all the pink semilandmark curve references are the same, so I'll use the same file (curveslide in the MORPHO IN ACTION folder for the pinks). line 116
pinkcurves <- as.matrix(read.csv("Data/curveslide.csv", header=T))


#transform into procrustes shape space
gpa.p2.orig <- gpagen(p2.lands, curves = pinkcurves)
plot(gpa.p2.orig)

#run a PCA #line 144
pca.p2.orig <- gm.prcomp(gpa.p2.orig$coords)
plot(pca.p2.orig)
plot(pca.p2.orig, axis1=1, axis2=3)

#add the wild/hatchery factor
wild.or.hatch <- read.csv("Data/pink2021.wildorhatch.andorder.csv")

#convert to data frame
df.pink2021.orig <- as.data.frame(pca.p2.orig$x[,1:3])

df.pink2021.orig$origin <- wild.or.hatch$Oto
df.pink2021.orig$ID <- wild.or.hatch$Fish.ID

#plot PCA/RWA (Principal component analysis/ relative warp analysis. They are the same thing here.)
##line 190
wildhatch <- c("orange", "blue")

ggplot(data=df.pink2021.orig)+aes(x=Comp3, y=Comp1, color=origin)+geom_point()+scale_color_manual(values=wildhatch) + stat_ellipse() + theme_bw() + labs(x= "RW 3", y= "RW 1")

ggplot(data=df.pink2021.orig)+aes(x=Comp2, y=Comp1, color=origin)+geom_point()+scale_color_manual(values=wildhatch) + stat_ellipse() + theme_bw() + labs(x= "RW 2", y= "RW 1")


#What do these RW's mean???
##Comp 1: hump, and a little bit of snout
##Comp 2: bendy fish
##comp 3: overall fish largeness, some depth involved
##Comp4: snout, strech of hump linearly
#pca plot ref to target
ref <- mshape(gpa.p2.orig$coords) 
ref_p2 <- mshape(gpa.p2.orig$coords) 
#comp1
plotRefToTarget(pca.p2.orig$shapes$shapes.comp1$min, ref, method="vector")
plotRefToTarget(pca.p2.orig$shapes$shapes.comp1$max, ref, method="vector")
#comp2
plotRefToTarget(pca.p2.orig$shapes$shapes.comp2$min, ref, method="vector")
plotRefToTarget(pca.p2.orig$shapes$shapes.comp2$max, ref, method="vector")
#comp3
plotRefToTarget(pca.p2.orig$shapes$shapes.comp3$min, ref, method="vector")
plotRefToTarget(pca.p2.orig$shapes$shapes.comp3$max, ref, method="vector")
#comp4
plotRefToTarget(pca.p2.orig$shapes$shapes.comp4$min, ref, method="vector")
plotRefToTarget(pca.p2.orig$shapes$shapes.comp4$max, ref, method="vector")


#plot wild to hatch comparison #line 289
#pinks are in Y.gpa ??
Y.gpa <- gpa.p2.orig
ref4 <- mshape(gpa.p2.orig$coords)
w.p <- numeric()
h.p <- numeric()
i=1

for (i in 1:length(df.pink2021.orig$origin)){
  if(df.pink2021.orig$origin[i]=="W"){
    w.p <- c(w.p, i)
  }else{
    h.p <- c(h.p, i)
  }
}

ref.w.p2 <-mshape(Y.gpa$coords[,,w.p]) #wild mean
ref.h.p2 <-mshape(Y.gpa$coords[,,h.p]) #hatchery mean

#wild to hatch, hatch to wild examine
plotRefToTarget(ref.w.p2, ref.h.p2, method="vector", mag=5) #what is happening here? is wild plotted to hatch or hatch plotted to wild?
plotRefToTarget(ref.h.p2, ref.w.p2, method="vector", mag=5) #hatchery dots, wild arrows
#HYPOTHESIS SUPPORTED! DAMN!


plotRefToTarget(ref.w.p2, ref.h.p2, method="vector", mag=5) #wild dot, hatchery arrows
plotRefToTarget(ref.w.p2, ref4, method="vector", mag=5) #wild to mean fish
plotRefToTarget(ref.h.p2, ref4, method="vector", mag=5) #hatchery to mean fish


#PCA/RWA influence plots
#can you also plot the ID of the important PCs???
pvar2 <- (pca.p2.orig$sdev^2)/(sum(pca.p2.orig$sdev^2)) #sooo...does this work?
names(pvar2) <- seq(1:length(pvar2))
barplot(pvar2, xlab= "Principal Components", ylab = "% Variance")

#also plot the CVA, please
vari <- gpa.p2.orig$coords
facto <- df.pink2021.orig$origin

CVA.pink <- CVA(vari, groups=facto)


#significance tests: MANCOVA
Y.gpa <- gpa.p2.orig
gdf <- geomorph.data.frame(Y.gpa, origin = df.pink2021.orig$origin)

fit.p <- procD.lm(coords ~ log(Csize) + origin, data = gdf,  #all coords taken into account
                  iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.p)  #says YES, significant difference between wild and hatch

#are interaction effects a thing?
fit.p.int <- procD.lm(coords ~ log(Csize) + origin + log(Csize)*origin, data = gdf,  #all coords taken into account
                      iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.p.int)


#we test for date below

#RESULTS TABLE CREATION (from fit.p)
sum_p2 <- summary(fit.p)
names(sum_p2)
sum_p2$table

write.csv(sum_p2$table, "Results/morpho_p2_table.csv")


#####################################################################################################################################



###################################################################################################################
#testing for date: an exploration (but no results to take away from this block, I believe)
#what dataframe is needed for this combination??
p_2_old <- read.csv("Data/Male.p2.Rdata.2_alt_ID.csv") #not working...why?
names(p_2_old)
p_2_old$ID
df.pink2021.orig$ID
##sigh. They're not identical. I'll need to data wrangle - 05/03/22
df_pink2021_orig_05 <- left_join(df.pink2021.orig, p_2_old, by="ID")
names(df_pink2021_orig_05)
df_pink2021_orig_05$Date


#05/27/22: I'll make a Jdate adjustment, it will make things work betteer
library(lubridate)
names(df_pink2021_orig_05)
Dates_toJ <- mdy(df_pink2021_orig_05$Date)
Julian_morph <- julian(Dates_toJ, origin=as.Date("2021-01-01"))
class(Julian_morph)
str(Julian_morph)
df_pink2021_orig_05$Julian <- as.vector(Julian_morph)


Y.gpa <- gpa.p2.orig
gdf <- geomorph.data.frame(Y.gpa, origin = df_pink2021_orig_05$origin, date=df_pink2021_orig_05$Julian )
gdf_Julian <- gdf
df_p2_reversed_Julian <- as.data.frame(gdf$coords, gdf$origin, gdf$date, gdf$Csize)
names(gdf)

fit.p_04 <- procD.lm(coords ~ log(Csize) + origin, data = gdf,  #all coords taken into account
                     iter = 9999, print.progress = FALSE, RRPP=TRUE)
summary(fit.p_04) 

fit.p_05 <- procD.lm(coords ~ log(Csize) + origin + date, data = gdf,  #all coords taken into account
                     iter = 9999, print.progress = FALSE, RRPP=TRUE)

fit.p_06 <- procD.lm(coords ~ log(Csize) + date + origin, data = gdf,  #all coords taken into account
                     iter = 9999, print.progress = FALSE, RRPP=TRUE)

fit.p_07 <- procD.lm(coords ~ log(Csize) +origin+ date, data = gdf,  #all coords taken into account
                     iter = 9999, print.progress = FALSE, RRPP=TRUE)

fit.p_08 <- procD.lm(coords ~ log(Csize) * date * origin, data = gdf,  #all coords taken into account
                     iter = 9999, print.progress = FALSE, RRPP=TRUE) #no interactions

#in the case of date not sig, we go back to original


summary(fit.p_04)#this one if date is not sig 
summary(fit.p_05)
summary(fit.p_06) #this one? #I think this one... #with origin last...?
summary(fit.p_07)
summary(fit.p_08)

#date is significant predictor of shape, until I tested J-date, now its not.


#gotta graph to figure out date relationship on coords

#AIC(fit.p_04)
#AIC(fit.p_05)
anova(fit.p_04, fit.p_05) #fit 5 is better than 4. Keep the mod with date. Graph shape vs. date? How?
names(fit.p_05)
fit.p_05$aov.table

#AIC(fit.p_04, fit.p_05, fit.p_06, fit.p_07, fit.p_08) #argh I forgot this takes forever

#date was not included in the model because turns out it is more driven by location than... by date


#################################################################################################################################

###############################################################################################################################
#################################################################################################################################
############################################p2 (pink odd, pink year 2, pink 2021) linear morpho analysis#########################
########################################## AND the creation of linear morpho results dataframes##################################3
###################################################################################################################################

#######################################################################################################
#pink 2021 linear morph data
p2.males <- read.csv("Data/Male.p2.Rdata.2.csv")


#exploratory graphs
#use graphs you already made, dummy!
#fuck, can't find them
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

write.csv(SNOUT_RESULTS, "Results/morpho_snout_results.csv")
write.csv(DEPTH_RESULTS, "Results/morpho_depth_results.csv")

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
S_6 <- lm(log(Snout.length.mm.) ~ Length.mm. + Otolith.reading, data=p2.males_dateadj)
S_7 <- lm(log(Snout.length.mm.) ~ Otolith.reading, data=p2.males_dateadj)

B_glob <- lm(log(Body.depth.mm.) ~ Length.mm. * Otolith.reading * Julian, data=p2.males_dateadj)
B_2 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Length.mm.:Julian + Otolith.reading:Julian, data=p2.males_dateadj)
B_3 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Length.mm.:Otolith.reading + Otolith.reading:Julian, data=p2.males_dateadj)
B_4 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian + Otolith.reading:Julian, data=p2.males_dateadj)
B_5 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading + Julian, data=p2.males_dateadj)
B_6 <- lm(log(Body.depth.mm.) ~ Length.mm. + Otolith.reading, data=p2.males_dateadj)
B_7 <- lm(log(Body.depth.mm.) ~ Otolith.reading, data=p2.males_dateadj)


AIC(B_glob, B_2, B_3, B_4, B_5, B_6, B_7)
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


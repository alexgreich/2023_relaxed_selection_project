#######################################
###male morpho analysis cleanup - dec 2023
####################################

#year one data will be first, then year two.

#year one (pink 2020, coho) is copy-pasted from "Morpho_yr1_pinkcoho.Rmd" 
#and year two (pink 2021) is copy-pasted from "Morpho_yr2_Rmd"
##so if any issues, look there
#I may also add on the simply graphing.R file at the end, if it makes sense


#################################################3
#load packages- do I need ALL of these?
library(ggplot2)
library(geomorph)
#library(shapes)
#library(MASS)
#library(ape)
#library(vegan)
#library(rgl)
#library(scatterplot3d)
library(ggfortify)
library(dplyr)
#library(reshape2)
#library(reshape)
library(tidyr)
library(ggforce)
library(Morpho)
library(RColorBrewer)
library(cowplot)
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
#summary(pca_coho)
#names(pca_coho)


#pca_coho_long <- gm.prcomp(gpa_coho_long$coords)
#plot(pca_coho_long)
#summary(pca_coho_long)
#names(pca_coho_long)

#can you also plot the ID of the important PCs???
#pvar_cl <- (pca_coho_long$sdev^2)/(sum(pca_coho_long$sdev^2))
#names(pvar_cl) <- seq(1:length(pvar_cl))
#barplot(pvar_cl, xlab= "Principal Components", ylab = "% Variance")

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
write.csv(sum_c$table, "Results/morpho_c_table.csv")
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
sd_s_w <- 10*sd #qc'ed, checks out.

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
plot(aov.depth.c) #looks fine ish


#prepare the table
df_coho_depth_results <- data.frame(t=sum_c_depth$coefficients[3,3], df= sum_c_depth$df[2], p_one_sided=1-sum_c_depth$coefficients[3,4]/2,
                                    length_t = sum_c_depth$coefficients[2,3], length_p_two_Sided = sum_c_depth$coefficients[2,4], 
                                    hatchery_mean= mean_d_hc, hatchery_sd = sd_d_hc,
                                    wild_mean=mean_d_wc, wild_sd=sd_d_wc)

#some density plots
ggplot(df.coho) + aes(x=snoutL) + geom_density()
################################


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
#NOW linear regression on that... how to do again?


df.pink2$snout <- snout_df.no56$snoutL
df.pink2


#first separate wild and hatchery
wild_S.p <-  df.pink2 %>% filter(Wild.or.hatch=="W")
hatch_S.p <-  df.pink2 %>% filter(Wild.or.hatch=="H")

#wild_SL.p <- wild_S.p$`snout_df.p[-56, ]`
#hatch_SL.p<- hatch_S.p$`snout_df.p[-56, ]`


#pink 2020 sd
sd(hatch_S.p$snout)
sd(wild_S.p$snout)

names(wild_S.p)
sd(hatch_S.p$depth)
sd(wild_S.p$depth)

#find sample size:
length(hatch_S.p$depth)
length(wild_S.p$depth)

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



#exploratory
qqnorm(df.pink2$snout)
qqline(df.pink2$snout)
qqnorm(log(df.pink2$snout))
qqline(log(df.pink2$snout))
ggplot(df.pink2) + aes(x=log(snout)) + geom_density()
ggplot(df.pink2) + aes(x=log(depth)) + geom_density()
ggplot(df.pink2) + aes(x=depth) + geom_density()

#actual analysis - no more logging
p1_snout_mod <- lm(snout~ length + factor(Wild.or.hatch), data=df.pink2)
plot(p1_snout_mod)

p1_depth_mod <- lm(depth~ length + factor(Wild.or.hatch), data=df.pink2)
plot(p2_depth_mod)

#interaction effects?
mod.p.d.l.int <- aov(depth ~ length + Wild.or.hatch + length*Wild.or.hatch, df.pink2)
summary.lm(mod.p.d.l.int) #not sig

mod.p.s.l.int <- aov(snout ~ Wild.or.hatch+length + length*Wild.or.hatch, df.pink2)
summary.lm(mod.p.s.l.int) #no int effects

#results


###########################################

#line 1495- formal length tests, hatch vs wild, both for coho and pink 2020


#problem: the snout analysis above looks incomplete. I also need the: depth analysis. Where is that? Is it on a different R file, that I didn't copy over
##something is fishy here, and not in a good way...
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



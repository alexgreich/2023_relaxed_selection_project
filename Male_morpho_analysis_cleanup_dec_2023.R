#######################################
###male morpho analysis cleanup - dec 2023
####################################

#year one data will be first, then year two.

#year one (pink 2020, coho) is copy-pasted from "Morpho_yr1_pinkcoho.Rmd" 
#and year two (pink 2021) is copy-pasted from "Morpho_yr2_Rmd"
##so if any issues, look there


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


#############################################
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


#######################################################################
#Now run a PCA/RWA (principal component analysis aka Relative Warp Analysis)

###########################################################################




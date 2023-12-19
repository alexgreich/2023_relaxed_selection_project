#####################################################
#####################################################
#Alexandra Reich
#Relaxed selection project
#Post-hoc power analyses
#Creation date: 12/19/23
#################################################
##################################################

#OBJECTIVE: do post-hoc power analyses on... the main analyses
###then create a table for co-authors
#There are a lot of analyses... We'll go through one by one

#maybe do indivudually in their own R data files, that way I don't have to read in data again?


library(pwr)

###########################################################################
#MALE 
##these were general linear models with R, so use: pwr.f2.test for all snout and depth for males
?pwr.f2.test
#############################################################################
#P1

#snout

#depth


#P2

#snout

#depth


#C

#snout

#depth


################################################################################
#FEMALE
##GSI - two sample, one tailed t-test, sometimes with same # of samples, sometimes with diff
?pwr.t.test #one and two samples for equal sample sizes
?pwr.t2n.test #two samples of differnt sizes
################################################################################
#P1

#GSI

#eggs


#P2

#GSI

#eggs


#C

#GSI

#eggs
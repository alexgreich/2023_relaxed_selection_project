##Figure 3 graphing##
##all lengths compared figure

library(ggplot2)
library(cowplot)
library(dplyr)


#############################3
#female coho (all)
##make sure the non-sub version has run to give you the correct c.GSI.clean
##or just load it here:
c.GSI.clean <- read.csv("Data/c.GSI.clean.csv")

coho_f_boxplot<-ggplot(c.GSI.clean) + aes(y=Length..mm., x=Wild.or.Hatch, color=Wild.or.Hatch) +
  geom_boxplot(color="black") +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y=element_blank(), x=element_blank()) +
  coord_cartesian(ylim=c(473,650))+
  scale_y_continuous(breaks=c(500,550,600,650), expand=c(0,0)) +
  theme(plot.margin = margin(t=7,r=12,l=5,b=5))+
  scale_x_discrete(labels=c("H", "W"))


#dev.new (width = 7, height = 4, unit = "in", noRStudioGD = T); last_plot() #perfect

#ggsave("Length_compared_femalecoho.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()
#####################################

##############################
#male coho (all)
##LOAD IN coho.dat!!!
coho.dat <- read.csv("Data/MASTERMaleCohoQC copy.csv")
#write.csv(coho.dat, "coho.dat")
#coho.dat <- read.csv("coho.dat")

coho_m_boxplot<-ggplot(coho.dat) + aes(y=Length..mm., x=Wild.or.Hatch, color=Wild.or.Hatch) +
  geom_boxplot(color="black", outlier.alpha = 0) +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y="Length (mm)", x=element_blank()) +
  scale_x_discrete(labels=c("H", "W")) +
  scale_y_continuous(breaks=c(450, 550, 650), expand=c(0,0))+
  coord_cartesian(ylim=c(400,650))#+
  #ylab("Coho salmon") #hmm. I'll have to add the external titles in ppt or something. or just do it the coding long way
 


#updated to here dec 2023

#####
even_male_pink <- p.dat <- read.csv("Data/MASTERMalePinksQC.reorder.csv")
## dec 2023 edit: we remove an outlier from the male dataset?
##yes, do it. and re-generate graph
even_male_pink <- even_male_pink %>% filter(Body.Depth.mm. < 183) ### dec 2023 edit here

even_female_pink <- read.csv("Data/FemalePinkGSIotodataadd1_without_gaps copy.csv")
odd_male_pink <-  read.csv("Data/Male.p2.Rdata.2.csv")
odd_female_pink <- read.csv("Data/Female.p2.Rdata.3.csv")
##eventually smooth this part out
#even-year pinks
##Wait -> I think we need to ezclude the big hatchery fish.
#males - use wild or hatch to differentiate, we don't have otos
even_male_boxplot <-ggplot(even_male_pink) + aes(y=Length..mm., x=Wild.or.Hatch, color=Wild.or.Hatch) +
  geom_boxplot(color="black", outlier.alpha = 0) +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y=element_blank(), x=element_blank()) +
  scale_x_discrete(labels=c("H", "W"))+
  scale_y_continuous(breaks=c(400, 440, 480, 520), expand=c(0,0))+
  coord_cartesian(ylim=c(366, 520))+
  #xlab("Even-year pink salmon")+
  ylab("Length (mm)")+
  labs(title = "Male")+
  theme(plot.title = element_text(hjust = 0.5))
#females
library(forcats) 
?fct_rev
even_female_pink_clean <- even_female_pink %>% filter(Weird == "n", Otolith.results != "unknown") 
even_female_boxplot<-ggplot(even_female_pink_clean) + aes(y=Length..mm., x=fct_rev(Otolith.results), color=fct_rev(Otolith.results)) +
  geom_boxplot(color="black", outlier.alpha = 0) +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y=element_blank(), x=element_blank()) +
  scale_x_discrete(labels=c("H", "W"))+
  scale_y_continuous(breaks=c(400, 440, 480), expand=c(0,0))+
  coord_cartesian(ylim=c(395,480))+
  labs(title="Female")+
  theme(plot.title = element_text(hjust = 0.5))
##looks good there

#############33
#pink odd
#males
odd_male_pink_clean <- odd_male_pink %>% filter(Otolith.reading != "Overground")
odd_male_boxplot <- ggplot(odd_male_pink_clean) + aes(y=Length.mm., x=fct_rev(Otolith.reading), color=fct_rev(Otolith.reading)) +
  geom_boxplot(color="black", outlier.alpha = 0) +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y=element_blank(), x=element_blank()) +
  scale_x_discrete(labels=c("H", "W"))+
  scale_y_continuous(breaks=c(350, 400, 450, 500), expand=c(0,0))+
  coord_cartesian(ylim=c(342,500))+ #looks good
  #xlab("Odd-year pink salmon")
  ylab("Length (mm)")

#odd_male_boxplot + 
  #annotate("Odd-year male pink salmon", hjust=0, vjust =1, angle=90)+
 # theme( plot.title=element_text(angle=90, hjust=0, vjust=0.5), plot.title.position = "plot")

#females
odd_female_pink_clean <- odd_female_pink %>% filter(Weird == "n", Oto.reading != "No Oto", Oto.reading != "Overground") 
odd_female_boxplot <- ggplot(odd_female_pink_clean) + aes(y=Length.mm., x=fct_rev(Oto.reading), color=fct_rev(Oto.reading)) +
  geom_boxplot(color="black", outlier.alpha = 0) +
  geom_jitter(size=2) +
  scale_color_manual(values =c("orange", "blue"))+
  theme_cowplot()+
  guides(color="none") +
  labs(y=element_blank(), x=element_blank()) +
  scale_x_discrete(labels=c("H", "W")) +
  scale_y_continuous(breaks=c(375, 425, 475), expand=c(0,0))+
  coord_cartesian(ylim=c(340,475))
  
##Qced positions, it's correct (min length is 344, which is a hatchery fish)




#link everything together with patchwork
##I'll need to adjust axes later
###I'm thinking 6 panel, 3 rows by 2 columns. Even year pink, odd year pink, coho by row
library(patchwork)
#(even_male_boxplot/odd_male_boxplot/coho_m_boxplot)+(even_female_boxplot/odd_female_boxplot/coho_f_boxplot)

length_boxplot_base <-(even_male_boxplot+even_female_boxplot)/(odd_male_boxplot+odd_female_boxplot)/(coho_m_boxplot+coho_f_boxplot)
#I ENDED UP USING THIS ONE ALL IS WELL
#SHIT, HOW DO I SAVE PLOTS AGAIN?
# 7 BY 8 DIMENSIONS

length_boxplot_base

#ggdraw(length_boxplot_base) + draw_label("Even-year pink salmon", x = 0, y = 0.88, size = 15, angle=90) +  #giving up and doing this on ppt.
 # draw_label ((ylab_GSI), angle= 90, x = 0.02, y = 0.6, size = 15) + #gi
  #draw_label ("Even-year pink", x = 0.85, y = 0.98, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Odd-year pink", x = 0.85, y = 0.67, fontfamily = "Arial", fontface="bold", size = 13) + draw_label ("Coho", fontfamily = "Arial", fontface="bold", size = 13, x = 0.85, y = 0.37) 


dev.new (width = 8, height = 10, unit = "in", noRStudioGD = T); last_plot() #perfect
#ggsave ("Plots/FIG3_LENGTH.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()
dev.off()

###ALL BELOW IS NOT SO RELEVANT I THINK

##this one looks better. I'll have to add labels and make dots smaller
###How to add the labels tho?
###I'd like to add labels: MALE, FEMALE to colums and EP, OP, and C (or of the sort) to the rows
#01/11/23.FIGURED IT OUT!
#(even_male_boxplot+even_female_boxplot)/(odd_male_boxplot+odd_female_boxplot)/(coho_m_boxplot+coho_f_boxplot)+
  #plot_layout(tag_level = 'new') +
  #plot_annotation(tag_levels = list(c('CM', 'CF', 'EPM', 'EPF', 'OPM', 'OPF' )))
#on the rght track, but I'm not stoked
##
#(even_male_boxplot+even_female_boxplot)/(odd_male_boxplot+odd_female_boxplot)/(coho_m_boxplot+coho_f_boxplot)+



####how to add them labels.... Look at old stuff? Code from data viz class
###get rid of how everything says length? In the first two plots

#nice, we've fixed the axes. Good job
#now let;s add some labels
#length_boxplot_base + annotate("text", x=0.100, y=0.100, hjust=0, size=10, label= "Female", family="Times New Roman")#+
  #annotate("text", x=1, y=5600, hjust=0, size=5, label="Male", family="Times New Roman")
#well it's esentially done, just need to add the labels

#might give up and just add labels in ppt. Sometimes doing things in R makes sense, sometiems it does not.



###I'll do some length comparison here##
#even vs odd pink length (megan wants to know)
t.test(odd_male_pink_clean$Length.mm., even_male_pink$Length..mm.)
#sig that even is larger
t.test(odd_female_pink_clean$Length.mm., even_female_pink_clean$Length..mm.)
#t.test that even is larger


#sample sizes
W <- c.GSI.clean %>% filter(Wild.or.Hatch=="wild")
H <- c.GSI.clean %>% filter(Wild.or.Hatch=="hatchery")
length(W$ID)
length(H$ID)


#################################3
#extracting male length data restuls
#p1
##even_male_pink
h_p1_m<-even_male_pink %>% filter(Wild.or.Hatch== "hatchery")
w_p1_m<-even_male_pink %>% filter(Wild.or.Hatch== "wild")

mean_length_p1_hm <- mean(h_p1_m$Length..mm.)
mean_length_p1_wm <- mean(w_p1_m$Length..mm.)

sd_length_p1_hm <- sd(h_p1_m$Length..mm.)
sd_length_p1_wm <- sd(w_p1_m$Length..mm.)

p1_m_L_ttest <- t.test(h_p1_m$Length..mm., w_p1_m$Length..mm., var.equal=T)

t_p1ml <- p1_m_L_ttest$statistic
p_p1ml <- p1_m_L_ttest$p.value
df_p1ml <- p1_m_L_ttest$parameter


#p2
h_p2_m<-odd_male_pink_clean %>% filter(Otolith.reading== "PORT ARMSTRONG")
w_p2_m<-odd_male_pink_clean %>% filter(Otolith.reading== "No Mark")

mean_length_p2_hm <- mean(h_p2_m$Length.mm.)
mean_length_p2_wm <- mean(w_p2_m$Length.mm.)

sd_length_p2_hm <- sd(h_p2_m$Length.mm.)
sd_length_p2_wm <- sd(w_p2_m$Length.mm.)

p2_m_L_ttest <- t.test(h_p2_m$Length.mm., w_p2_m$Length.mm., var.equal=T)

t_p2ml <- p2_m_L_ttest$statistic
p_p2ml <- p2_m_L_ttest$p.value
df_p2ml <- p2_m_L_ttest$parameter

#c
##coho.dat

h_c_m<-coho.dat %>% filter(Wild.or.Hatch== "hatch")
w_c_m<-coho.dat %>% filter(Wild.or.Hatch== "wild")

mean_length_c_hm <- mean(h_c_m$Length..mm.)
mean_length_c_wm <- mean(w_c_m$Length..mm.)

sd_length_c_hm <- sd(h_c_m$Length..mm.)
sd_length_c_wm <- sd(w_c_m$Length..mm.)

c_m_L_ttest <- t.test(h_c_m$Length..mm., w_c_m$Length..mm., var.equal=T)

t_cml <- c_m_L_ttest$statistic
p_cml <- c_m_L_ttest$p.value
df_cml <- c_m_L_ttest$parameter


#combine and write csv for length restuls
male_length_results <- data.frame(
  t= c(t_p1ml, t_p2ml, t_cml),
  df = c(df_p1ml, df_p2ml, df_cml),
  p = c(p_p1ml, p_p2ml, p_cml),
  hatch_mean = c(mean_length_p1_hm, mean_length_p2_hm, mean_length_c_hm),
  hatch_sd = c(sd_length_p1_hm, sd_length_p2_hm, sd_length_c_hm),
  wild_mean = c(mean_length_p1_wm, mean_length_p2_wm, mean_length_c_wm),
  wild_sd = c(sd_length_p1_wm, sd_length_p2_wm, sd_length_c_wm)
)

rownames(male_length_results) <- c("pink 2020", "pink 2021", "coho")

write.csv(male_length_results, "Results/Male length results.csv")




#01/26/24
#the graph that sam requested: rationalizing combining sashin and lovers, male and female. Make it like your length graph

##p1.clean and p2.GSI.clean. Cleaning up plots L11+L33
L111 <- ggplot(even_female_pink_clean) + aes(y=Length..mm., x=Location, color=Otolith.results) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Even-year pink females")+
  labs(y=element_blank(),x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(390, 480), breaks=c(400, 440, 480), expand=c(0,0))
#?position_dodge()

L333 <- ggplot(odd_female_pink_clean) + aes(x=Location, y=Length.mm., color=Oto.reading) +
  geom_boxplot() + geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Odd-year pink females")+
  labs(y=element_blank(), x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(340, 480), breaks=c(360, 400, 440, 480), expand=c(0,0))

#male even
M111 <- ggplot(even_male_pink) + aes(y=Length..mm., x=Location, color=factor(Wild.or.Hatch, levels=c("wild", "hatchery"))) +
  geom_boxplot(outlier.color="white", outlier.fill="white") +geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Even-year pink males")+
  labs(y="Length (mm)", x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(368, 520), breaks=c(400, 440, 480, 520), expand=c(0,0))


#male odd
L444 <- ggplot(odd_male_pink_clean) + aes(x=Location, y=Length.mm., color=Otolith.reading) +
  geom_boxplot() + geom_jitter(position=position_dodge(width=0.75)) + scale_color_manual(values=c("blue", "orange"), labels=c("Wild", "Hatchery"))+
  theme_cowplot()+
  ggtitle("Odd-year pink males")+
  labs(y="Length (mm)", x=element_blank(), color=NULL)+
  scale_x_discrete(labels=c("Port Armstrong", "Lovers Cove", "Sashin Creek"))+
  scale_y_continuous(limits=c(340, 500), breaks=c(350, 400, 450, 500), expand=c(0,0))

(M111 + L111)/(L444 + L333) + plot_layout(guides="collect")

dev.new (width = 10, height = 6, unit = "in", noRStudioGD = T); last_plot() #perfect

ggsave("Plots/Length_compared_loverssashin.jpg", width = dev.size()[1], height = dev.size()[2]); dev.off()

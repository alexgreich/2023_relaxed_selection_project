#snag
#Alexandra Reich
#12/2/24

#snag graph for reviewers


#libraries
library(tidyverse)
library(ggplot2)
library(cowplot)

#data
P_f_o <- read.csv("Data/Snag/p2_female_snag.csv")
P_m_o <- read.csv("Data/Snag/p2_male_snag.csv")

#graph
P_f_o_w <- P_f_o %>% filter(Oto.reading == "No mark")
ggplot(P_f_o) + aes(x=Catch.method, y=Length.mm.) +
  geom_boxplot() + geom_point()

P_m_o_w <- P_m_o %>% filter(Otolith.reading == "No Mark")
plot_m <-ggplot(P_m_o) + aes(x=Catch.method, y=Length.mm.) +
  geom_boxplot() + geom_point() + theme_cowplot() +
  labs(y="Length (mm)", x="Catch method")

ggsave(filename = "Plots/male geartype plot.jpg",plot_m)

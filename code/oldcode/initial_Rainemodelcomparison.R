# Initial draft code to assess Raine's model output

library(tidyverse)

d1 <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSites.csv") %>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_kelpbiomass")

d2 <- read.csv("data/Large_files/Raine_modeloutput/EndModelRunsAllSitesObsKelp.csv")%>%
  pivot_longer(cols = c(ABUR1:WOOD3), names_to = "id", values_to = "predicted_kelpbiomass")


d2 %>%
  separate(id, into = c("site", "transect"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  ggplot(aes(x = tmpt, y = predicted_kelpbiomass))+
  geom_line(aes(color = group, linetype = as.factor(transect)))+
  facet_wrap(~site)







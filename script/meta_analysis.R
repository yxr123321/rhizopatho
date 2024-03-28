library(metafor)
library(magrittr)
library(tidyr)
library(dplyr)
make_pct <- function(x) (exp(x) - 1) * 100

richness=read.csv("data/all_shannon.csv",header=T)
d2<-escalc(measure="ROM",data=richness,m1i=average.y,sd1i=sd.y,n1i=n.y,m2i=average.x,sd2i=sd.x,n2i=n.x)
###Total####
total <- rma(yi,vi, data=d2,method = "DL")
total
total.n <- d2  %>% drop_na(.,yi) %>% summarise(n = n()) 
total.df <- coef(summary(total)) %>% mutate(type="", 
                                          factor="Total",
                                          size=total.n$n)
####Land use#####
Land.use<-rma(yi, vi, mods=~Crop_Managment -1, data=d2,method = "DL")
Land.use.n <- d2 %>% drop_na(.,yi) %>% group_by(Crop_Managment) %>% summarise(n = n()) %>% drop_na()
Land.use.df <- coef(summary(Land.use)) %>% mutate(type="Land use", 
                                                  factor=levels(as.factor(d2$Crop_Managment)),
                                                  size=Land.use.n$n)
##Crop#####
#d2 <- subset(d2,Crop_Managment!='Wetland')
Crop<-rma(yi, vi, mods=~Crop_Family -1, data=d2,method = "DL")
Crop.n <- d2 %>%  drop_na(.,yi) %>% group_by(Crop_Family) %>% summarise(n = n()) %>% drop_na()
Crop.df <- coef(summary(Crop)) %>% mutate(type="Crops", 
                                          factor=levels(as.factor(Crop.n$Crop_Family)),
                                          size=Crop.n$n)
meta.df <- bind_rows(total.df,Land.use.df,Crop.df)


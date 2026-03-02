library(tidyverse)
library(writexl)
library(readxl)
library(stringi)
library(stringr)
library(sf)
library(terra)
library(raster)
library(ape)
library(spatstat)

path <- "//files.ugent.be/modhondt/shares/endure/SUSANA/plantenbakken Zeebrugge/ROOTS"
setwd(path)
getwd()

roots <- read_xlsx("./Roots.xlsx")
summary(roots)
roots$layer <- as.factor(roots$layer)
roots$layer <- ordered(roots$layer, levels = c("shoots", "top", "1", "2", "3", "4"))
roots$substrate <- as.factor(roots$substrate)
#roots$substrate <- ordered(roots$substrate, levels = c("D", "87.5", "75.0", "62.5", "50.0", "37.5", "25.0", "12.5", "W", "R"))
roots$substrate <- ordered(roots$substrate, levels = c("D", "87.5", "75", "62.5", "50", "37.5", "25", "12.5", "W", "R"))
summary(roots)

ggplot(roots , aes(x = as.factor(layer), y = biomass_per_L)) +
  geom_boxplot() +
  ggtitle("Marram grass biomass allocation according to depth") + labs(y="Total biomass (g) per L of substrate", x="Layer") 

ggplot(roots , aes(x = as.factor(layer), y = total_biomass)) +
  geom_point(aes(color=substrate)) +
  ggtitle("Marram grass biomass allocation on alternative substrates") + labs(y="Total biomass (g)", x="Layer") +
  facet_wrap(~as.factor(substrate))
 
ggplot(roots , aes(x = as.factor(substrate), y = total_biomass)) +
  geom_point(aes(color=substrate)) +
  ggtitle("Marram grass biomass allocation on alternative substrates") + labs(y="Total biomass (g)", x="Substrate") +
  facet_wrap(~as.factor(layer))     

ggplot(roots , aes(x = as.factor(substrate), y = biomass_per_L)) +
  geom_point(aes(color=substrate)) +
  ggtitle("Marram grass biomass allocation on alternative substrates") + labs(y="Total biomass (g) per L substrate", x="Substrate") +
  facet_wrap(~as.factor(layer)) 

ggplot(roots , aes(x = as.factor(layer), y = biomass_per_L)) +
  geom_point(aes(color=substrate)) +
  ggtitle("Marram grass biomass allocation on alternative substrates") + labs(y="Total biomass (g) per L substrate", x="Substrate") +
  facet_wrap(~as.factor(substrate))   


################################################################################

# SHOOT BIOMASS
shoots <- filter(roots, layer == 'shoots')

bakken_all <- read_excel("//files.ugent.be/modhondt/shares/endure/SUSANA/plantenbakken Zeebrugge/opvolging plantenbakken Zeebrugge.xlsx")
str(bakken_all)
bakken_all$row <- as.character(bakken_all$row)
bakken_all$substrate <- as.factor(bakken_all$substrate)
bakken_all$substrate <- ordered(bakken_all$substrate, levels = c("D", "87.5", "75", "62.5", "50", "37.5", "25", "12.5", "W", "R"))
bakken_all$plant <- as.factor(bakken_all$plant)
bakken_all$length <- as.numeric(bakken_all$length)
bakken_all$nr <- as.numeric(bakken_all$nr)

bakken_all$t <- as.factor(bakken_all$t)
bakken_all$t <- ordered(bakken_all$t, levels = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12", "t13"))

control       <- bakken_all[which(bakken_all$plant=='controle'),]
bakken        <- bakken_all[which(bakken_all$plant=='plant' | bakken_all$plant=='replant'),]

# bakken tem t13
bakken        <- bakken[which(bakken$t %in% c('t1', 
                                              't2', # no data for length
                                              't3', 
                                              't4', 
                                              't5', 
                                              't6',
                                              't7', 
                                              't8',
                                              't9',
                                              't10', 
                                              't11',
                                              't12',
                                              't13')),]

bakken_shoots <- merge(bakken,shoots,by=c("column", "row"))
bakken_shoots = subset(bakken_shoots, select = c(column, row, ID.x, substrate.x, plant, length, nr, t, time, date, layer, shoot_biomass) )
bakken_shoots$ID = bakken_shoots$ID.x
bakken_shoots$substrate = bakken_shoots$substrate.x
bakken_shoots = subset(bakken_shoots, select = -c(ID.x, substrate.x) )

bakken_shoots_t13 <- bakken_shoots[which(bakken_shoots$t == 't13'),]
bakken_shoots_t12 <- bakken_shoots[which(bakken_shoots$t == 't12'),]


ggplot(bakken_shoots, aes(x=shoot_biomass, y=nr)) + 
  geom_point(aes(colour=substrate)) +
  geom_smooth(method = "loess", colour="darkred", se=TRUE) +
  ggtitle("Marram nr of leaves vs shoots biomass") + labs(y="leaf nr", x="total shoot biomass") +
  facet_wrap(~as.factor(t), scale="free")

ggplot(bakken_shoots, aes(x=shoot_biomass, y=length)) + 
  geom_point(aes(colour=substrate)) +
  geom_smooth(method = "loess", colour="darkred", se=TRUE) +
  ggtitle("Marram height vs shoots biomass") + labs(y="Marran height (cm)", x="total shoot biomass") +
  facet_wrap(~as.factor(t), scale="free")

# nr of leaves doesnt seem correlated with shoot biomass, 
# the max height of the grass does


################################################################################


# MEAN ROOT BIOMASS
roots_layers234 <- roots[which(roots$layer == '2' | roots$layer == '3' | roots$layer == '4' ),]
roots_layers234$ID <- paste(roots_layers234$column, roots_layers234$row, sep="_")
roots_meanbiomass <- roots_layers234 %>% group_by(ID) %>%
  summarise( mean_biomass_per_L = mean(biomass_per_L)) 
bakken_roots <- merge(bakken,roots_meanbiomass,by=c("ID"))

ggplot(bakken_roots, aes(x=mean_biomass_per_L, y=nr)) + 
  geom_point(aes(colour=substrate)) + #xlim(0.25, 0.7) +
  geom_smooth(method = "loess", colour="darkred", se=TRUE) +
  ggtitle("Marram nr of leaves vs mean root biomass per L of substrate") + labs(y="leaf nr", x="mean root biomass per L of substrate") +
  facet_wrap(~as.factor(t), scale="free")

ggplot(bakken_roots, aes(x=mean_biomass_per_L, y=length)) + 
  geom_point(aes(colour=substrate)) +  #xlim(0.25, 0.7) +
  geom_smooth(method = "loess", colour="darkred", se=TRUE) +
  ggtitle("Marram height vs mean root biomass per L of substrate") + labs(y="Marran height (cm)", x="mean root biomass per L of substrate") +
  facet_wrap(~as.factor(t), scale="free")


################################################################################


# ROOT/SHOOT ratio
box_shoot_biomass <- aggregate(roots["shoot_biomass"], roots["box_ID"], sum, na.rm = TRUE)
box_root_biomass <- aggregate(roots["root_biomass"], roots["box_ID"], sum, na.rm = TRUE)
root_shoot_ratio <- merge(box_shoot_biomass, box_root_biomass, by=c("box_ID"))
roots <- merge(roots, root_shoot_ratio, by=c("box_ID"))
roots$root_shoot_ratio <- roots$root_biomass.y/roots$shoot_biomass.y

ggplot(roots , aes(x = shoot_biomass.y, y = root_biomass.y)) +
  geom_smooth(method = "lm", colour="darkred", se=T) +
  geom_point(aes(color=substrate, size=root_shoot_ratio)) +
  ggtitle("Root/shoot ratio") + labs(y="root biomass (g)", x="shoot biomass (g)")

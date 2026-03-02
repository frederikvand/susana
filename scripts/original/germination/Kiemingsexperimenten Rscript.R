
library("readxl")
library("ggplot2")
library("magrittr")
library("dplyr")
library("lattice")
library("car") 


# Read data
df_kiem <- read_excel("N:/shares/endure/SUSANA/Kiemingsexperimenten/kiemingsexperimenten data.xlsx")
#df_kiem <- read_excel("//files.ugent.be/modhondt/shares/endure/SUSANA/Kiemingsexperimenten/kiemingsexperimenten data.xlsx")

df_kiem$substrate <- as.factor(df_kiem$substrate)
df_kiem$substrate <- ordered(df_kiem$substrate, levels = c("D", "D75", "D50", "D25", "W", "R"))
df_kiem$t <- as.factor(df_kiem$t)
df_kiem$t <- ordered(df_kiem$t, levels = c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12"))
df_kiem$nr <- as.numeric(df_kiem$nr)
df_kiem$length <- as.numeric(df_kiem$length)
k = rep(96, 12)
df_kiem$time <- c(rep(1, 96),
                  rep(2, 96),
                  rep(3, 96),
                  rep(4, 96),
                  rep(5, 96),
                  rep(6, 96),
                  rep(7, 96),
                  rep(8, 96),
                  rep(9, 96),
                  rep(10,96),
                  rep(11,96),
                  rep(12,96))

str(df_kiem)


## plots
# number of germinated seeds (t1 and t2) or number of shoots (t6 and t12)
ggplot(df_kiem, aes(x=substrate, y=nr, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Number of germinated seeds") + labs(y="number of germinated seeds", x="substrate") +
  facet_grid(. ~ t)

# percentage of germinated seeds (t1 and t2) or number of shoots (t6 and t12)
df_kiem_t1 <- df_kiem[which(df_kiem$t=='t1'),]
df_kiem_t1$percentage_germinated <- df_kiem_t1$nr*10
df_kiem <- df_kiem[which(df_kiem$t!='t1'),]
ggplot(df_kiem_t1, aes(x=substrate, y=percentage_germinated, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Percentage of germinated seeds") + labs(y="% germinated seeds", x="substrate") +
  facet_grid(. ~ t)

df_kiem_t2 <- df_kiem[which(df_kiem$t=='t2'),]
df_kiem_t2$percentage_germinated <- df_kiem_t2$nr*10
ggplot(df_kiem_t2, aes(x=substrate, y=percentage_germinated, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Percentage of germinated seeds") + labs(y="% germinated seeds", x="substrate") +
  facet_grid(. ~ t)

df_kiem_t6 <- df_kiem[which(df_kiem$t=='t6'),]
ggplot(df_kiem_t6, aes(x=substrate, y=nr, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Number of shoots") + labs(y="Number of shoots", x="substrate") +
  facet_grid(. ~ t)

df_kiem_t11 <- df_kiem[which(df_kiem$t=='t11'),]
ggplot(df_kiem_t11, aes(x=substrate, y=nr, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Number of shoots") + labs(y="number of shoots", x="substrate") +
  facet_grid(. ~ t)

# length of seedlings per time
ggplot(df_kiem, aes(x=substrate, y=length, group=substrate)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("length of seedlings") + labs(y="length of seedlings", x="substrate") +
  facet_grid(. ~ t)

# length of seedlings per substrate
ggplot(df_kiem, aes(x=time, y=length, group=t)) + 
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("length of seedlings") + labs(y="length of seedlings", x="time") +
  facet_grid(. ~ substrate) 

# growth rate curves per substrate
ggplot(df_kiem, aes(y = length, x = time, colour = substrate)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05, jitter.height = 0.5, dodge.width = 0.4)) +
  geom_smooth(method = "loess", se = FALSE, inherit.aes = TRUE) +
  ggtitle("Length of seedlings") + labs(y="length of seedlings (cm)", x="time") +
  scale_x_continuous(breaks=c(1:11))



## ANOVA
#  1. is er een verschil in het aantal gekiemde zaden per substraat?
#  2. is er een verschil in de lengte van de bladeren per substraat?
#  3. is er een verschil in groeisterkte tussen de substraten? 
#  (groeisterkte = r van groeicurve)


# ANOVA 1. % gekiemde zaden
g = aov(percentage_germinated ~ substrate, data = df_kiem_t2)
Anova(g, type="III") # ***
TukeyHSD(g)

# ANOVA 2. length
l = aov(length ~ substrate * t, data = df_kiem)
Anova(l, type="III") # ***
TukeyHSD(l)

##  LOGGER DATA
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# read all logger data 
# T1 soil temperature
# T2 surface temperature
# T3 air temperature
logger_IDs <- c( "l_Dc", "l_75c", "l_50c", "l_25c", "l_Wc", "l_Rc")
#for (logger_ID in logger_IDs) {
#  #logger_ID = "l_Wc"
#  logger = read.csv(paste0("N:/shares/endure/SUSANA/Kiemingsexperimenten/logger data/", logger_ID, ".csv"), sep=";", header=F)
#  colnames(logger) = c("pkey", "datetime", "T1", "T2", "T3", "Moisture", "Vol")
#  logger$Time = as.POSIXct(logger$datetime, format= "%Y.%m.%d %H:%M") 
#  logger$Day = as.Date(logger$datetime, format = "%Y.%m.%d")
#  logger_exp = subset(logger, Day > "2024-04-06") # opzet was 5 april 
#  logger_exp = logger_exp %>% select(c("T1", "T2", "T3", "Moisture", "Vol", "Time", "Day"))
#  logger_exp = as_tibble(logger_exp)
#  logger_exp_pivot = logger_exp %>% tidyr::pivot_longer(cols = c("T1", "T2", "T3", "Moisture"), 
#                                                        names_to = "Clim_var", 
#                                                        values_to = "Clim_val")
#  # all variables together
#  ggplot(logger_exp_pivot, aes(y = Clim_val, x = Time, color = Clim_var)) +
#    geom_line(size = 0.5) +
#    scale_color_manual(values = c("turquoise4", "goldenrod", "darkgoldenrod", "darkorange4")) +
#    labs(title = logger_ID, y = "Value") +
#    facet_wrap(~Clim_var, scales = "free")
#  ggsave(paste0("./Lolly exported data/Logger", logger_ID, ".jpg"), width = 7, height = 7)
#  
#  ggplot(logger_exp, aes(x=Time)) +
#    geom_line(size=0.2, aes(y=T1), color="darkorange4") +
#    geom_line(size=0.2, aes(y=T2), color="chartreuse4") +
#    geom_line(size=0.2, aes(y=T3), color="steelblue") +
#    labs(title = logger_ID, y = "Soil temperatire")
#  
#  assign(paste0(logger_ID, "_", "df"), logger_exp)
#  assign(paste0(logger_ID, "_", "df_pivot"), logger_exp_pivot)
#}

l_Dc_df$substrate  <- rep("D", nrow(l_Dc_df))
l_Dc_df$color_code <- rep("#440154", nrow(l_Dc_df))
l_75c_df$substrate <- rep("D75", nrow(l_75c_df))
l_75c_df$color_code <- rep("#414487", nrow(l_75c_df))
l_50c_df$substrate <- rep("D50", nrow(l_50c_df))
l_50c_df$color_code <- rep("#2a788e", nrow(l_50c_df))
l_25c_df$substrate <- rep("D25", nrow(l_25c_df))
l_25c_df$color_code <- rep("#22a884", nrow(l_25c_df))
l_Wc_df$substrate  <- rep("W", nrow(l_Wc_df))
l_Wc_df$color_code <- rep("#7ad151", nrow(l_Wc_df))
l_Rc_df$substrate  <- rep("R", nrow(l_Rc_df))
l_Rc_df$color_code <- rep("#fde725", nrow(l_Rc_df))

l_Dc_df_pivot$substrate  <- rep("D", nrow(l_Dc_df_pivot))
l_Dc_df_pivot$color_code <- rep("#440154", nrow(l_Dc_df_pivot))
l_75c_df_pivot$substrate <- rep("D75", nrow(l_75c_df_pivot))
l_75c_df_pivot$color_code <- rep("#414487", nrow(l_75c_df_pivot))
l_50c_df_pivot$substrate <- rep("D50", nrow(l_50c_df_pivot))
l_50c_df_pivot$color_code <- rep("#2a788e", nrow(l_50c_df_pivot))
l_25c_df_pivot$substrate <- rep("D25", nrow(l_25c_df_pivot))
l_25c_df_pivot$color_code <- rep("#22a884", nrow(l_25c_df_pivot))
l_Wc_df_pivot$substrate  <- rep("W", nrow(l_Wc_df_pivot))
l_Wc_df_pivot$color_code <- rep("#7ad151", nrow(l_Wc_df_pivot))
l_Rc_df_pivot$substrate  <- rep("R", nrow(l_Rc_df_pivot))
l_Rc_df_pivot$color_code <- rep("#fde725", nrow(l_Rc_df_pivot))

loggers_all <- rbind(l_Dc_df, l_75c_df, l_50c_df, l_25c_df, l_Wc_df, l_Rc_df)
loggers_all$color_code <- ordered(loggers_all$color_code, levels = c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"))
loggers_all$substrate <- ordered(loggers_all$substrate, levels = c("D", "D75", "D50", "D25", "W", "R"))
loggers_all$splittime <- loggers_all$Time
loggers_all <- separate(loggers_all, splittime, c("date", "time"), sep = " ") # NA's are 00:00:00
loggers_all_4 <- loggers_all[which(loggers_all$time=='04:00:00'| loggers_all$time=='16:00:00'),]
  
loggers_all_pivot <- rbind(l_Dc_df_pivot, l_75c_df_pivot, l_50c_df_pivot, l_25c_df_pivot, l_Wc_df_pivot, l_Rc_df_pivot)
loggers_all_pivot$color_code <- ordered(loggers_all_pivot$color_code, levels = c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"))
loggers_all_pivot$substrate <- ordered(loggers_all_pivot$substrate, levels = c("D", "D75", "D50", "D25", "W", "R"))
loggers_all_pivot_T <- loggers_all_pivot[which(loggers_all_pivot$Clim_var=='T1' | loggers_all_pivot$Clim_var=='T2' | loggers_all_pivot$Clim_var=='T3'),]
loggers_all_pivot_T$splittime <- loggers_all_pivot_T$Time
loggers_all_pivot_T <- separate(loggers_all_pivot_T, splittime, c("date", "time"), sep = " ") # NA's are 00:00:00
loggers_all_pivot_T_4 <- loggers_all_pivot_T[which(loggers_all_pivot_T$time=='04:00:00'| loggers_all_pivot_T$time=='16:00:00'),]


loggers_all_pivot_Moisture <- loggers_all_pivot[which(loggers_all_pivot$Clim_var=='Moisture'),]

c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")


# Temperature graphs
ggplot(loggers_all, aes(x=Time)) +
  geom_line(size=0.2, aes(y=T1, color=substrate)) +
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Soil temperature" , y = "Soil temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)
min(loggers_all$T1) #2.31
max(loggers_all$T1) #45.5

ggplot(loggers_all, aes(x=Time)) +
  geom_line(size=0.2, aes(y=T2, color=substrate)) +
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Surface temperature", y = "Surface temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)
min(loggers_all$T2) #2.5
max(loggers_all$T2) #41.62

ggplot(loggers_all, aes(x=Time)) +
  geom_line(size=0.2, aes(y=T3, color=substrate)) +
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Air temperature", y = "Air temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)
min(loggers_all$T3) #3.31
max(loggers_all$T3) #38.12

ggplot(loggers_all, aes(x=Time)) + # Air temperature and moisture color scale
  geom_line(size=0.2, aes(y=T1, color=Vol)) +
  #scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Soil temperature and moisture content over time", y = "Air temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)

ggplot(loggers_all_pivot_T, aes(x=Clim_var, y=Clim_val)) +
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Climate variables") + labs(y="climate values", x="climate variable") +
  facet_grid(.~substrate)

ggplot(loggers_all_pivot_T, aes(x=substrate, y=Clim_val)) +
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Climate variables") + labs(y="climate values", x="climate variable") +
  facet_grid(.~Clim_var)

# Temperature 4am vs 4pm 
ggplot(loggers_all_4, aes(x=Time)) +
  geom_point(size=0.7, aes(y=T1, color=time)) +
  scale_color_manual(values=c("#33ceff", "#ff5733")) +
  labs(title = "Soil temperature", y = "Soil temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)

ggplot(loggers_all_4, aes(x=Time)) +
  geom_point(size=0.7, aes(y=T2, color=time)) +
  scale_color_manual(values=c("#33ceff", "#ff5733")) +
  labs(title = "Surface temperature", y = "Surface temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)

ggplot(loggers_all_4, aes(x=Time)) +
  geom_point(size=1, aes(y=T1, color=time)) + #soil
  geom_point(size=2, aes(y=T2, color=time)) + #surface
  geom_point(size=3, aes(y=T3, color=time)) + #air
  scale_color_manual(values=c("#33ceff", "#ff5733")) +
  labs(title = "Temperature", y = "Temperature (°C)", x = "Date") +
  facet_wrap(.~substrate)

ggplot(loggers_all_pivot_T_4, aes(x=Time)) +
  geom_point(size=1.2, aes(y=Clim_val, color=Clim_var)) + #soil
  geom_smooth(aes(y = Clim_val, group = Clim_var, color=Clim_var)) +
  scale_color_manual(values=c("#33ceff", "#ff5733", "yellow"),
                     labels = c("T1" = "soil", "T2" = "surface", "T3" = "air")) +
  labs(title = "Temperature", y = "Temperature (°C)", x = "Date") +
  guides(color = guide_legend(title = "Temperature")) +
  facet_wrap(.~substrate+time)

ggplot(loggers_all_pivot_T_4, aes(x=Time)) +
  geom_point(size=1.2, aes(y=Clim_val, color=substrate)) + #soil
  geom_smooth(aes(y = Clim_val, group = substrate, color=substrate)) +
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Temperature", y = "Temperature (°C)", x = "Date") +
  facet_wrap(.~Clim_var+time)

ggplot(loggers_all_pivot_T_4, aes(x=substrate, y=Clim_val)) +
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Temperature") + labs(y="substrate", x="temperature") +
  facet_grid(.~time+Clim_var)


# Moisture graphs
ggplot(loggers_all, aes(x=Time)) +
  geom_line(size=0.2, aes(y=Vol, color=substrate)) + #volumetric soil moisture (0-100% vol.)
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Soil moisture", y = "volumetric soil moisture (% vol.)")

ggplot(loggers_all, aes(x=Time)) +
  geom_line(size=0.2, aes(y=Moisture, color=substrate)) + #raw soil moisture values (~500-3600)
  scale_color_manual(values=c("#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725")) +
  labs(title = "Soil moisture", y = "Soil Moisture")

ggplot(loggers_all_pivot_Moisture, aes(x=substrate, y=Clim_val)) +
  geom_boxplot(aes(fill=substrate)) +
  ggtitle("Soil moisture") + labs(y="Soil moisture", x="substrate")

ggplot(loggers_all, aes(x=Vol, y=T1)) +
  geom_point(size=1, aes(color=Time)) +
  ggtitle("Soil moisture in function of soil temperature") + labs(y="Soil temperature (°C)", x="Vol") +
  facet_wrap(.~substrate, scale="free")

ggplot(loggers_all_4, aes(x=T1, y=Vol)) + # water uptake is higher at lower temperature
  geom_point(size=1, aes(color=Time)) +
  geom_smooth(aes(x = T1, y=Vol)) +
  ggtitle("Soil moisture in function of soil temperature") + labs(y="Vol", x="Soil temperature (°C)") +
  facet_wrap(.~substrate+time, scale="free")


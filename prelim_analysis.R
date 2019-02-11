# preliminary analysis

library(tidyverse)
library(lubridate)
library(viridis)
library(readxl)
library(lme4)
library(ggdag)
theme_set(theme_bw(base_size = 20)) 

photoperiod <- function(lat, doy, p = 1.5){
  theta <- 0.2163108 + 
    2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi <- asin(0.39795 * cos(theta))
  D <- 24 - (24 / pi) * acos(
    (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi))/
      (cos(lat * pi / 180) * cos(phi))
  )
}


# from uspest.org 2018 weather data
gdd <- read.csv("data/gdd2018.txt", row.names=NULL, sep="")
names(gdd) <- c("month", "day", "tmax", "tmin", "precip", "degdays", "accumdd")
gdd$date <- lubridate::mdy(paste(gdd$month, gdd$day, 2018, sep = "-"))
gdd$daylength <- photoperiod(44.5646, lubridate::yday(gdd$date), p = 1.5)

dat <- readxl::read_xlsx("data/CommonGarden_Flowers.xlsx", na = "NA")
flowermass <- dat %>% 
  group_by(Population, Pool) %>% 
  summarise(season_total = sum(FlowerMass, na.rm = TRUE))

rootmass <- readxl::read_xlsx("data/CommonGarden_Roots.xlsx", na = "NA")


# basic plot of mean flower mass by population
boxplt <- ggplot(flowermass, aes(x = Population, y = season_total)) +
  geom_boxplot()
boxplt

boxplt <- ggplot(rootmass, aes(x = Population, y = Root_mass)) +
  geom_boxplot()
boxplt

# extract explanatory variables from other datasets
dat <- readxl::read_xlsx("data/CommonGarden_Setup.xlsx", na = "NA")
controlvars <- dat %>% 
  group_by(Population, Pool) %>% 
  mutate(Aphids = (Aphid_July15 + Aphid_July2) / 2,
         Senescence = (PlantSenescence_Aug14 + PlantSenescence_July29) / 2,
         Regrowth = (PlantRegrowth_Oct16 + PlantRegrowth_Sept28) / 2)

dat <- readxl::read_xlsx("data/CommonGarden_F1_AdultCounts.xlsx", na = "NA")
names(dat)[4:15] <- paste0("x", names(dat[4:15]))
counts <- dat %>% 
  gather(key = date, value = counts, x43276:x43304) %>% 
  mutate(date = as.Date(as.numeric(sub(pattern = "x", replacement = "", x = date, fixed = TRUE)),
                        origin = "1899-12-30")) %>% 
  filter(Pool != "ALL") %>% 
  mutate(Pool = as.numeric(as.character(Pool)))

# notice spikes/drops after the redistribution when beetles collected for diapause tests
# different pools mixed at that point
tsplt <- ggplot(counts, aes(x = date, y = counts, group = Pool, color = Pool)) +
  geom_path() +
  facet_wrap(~Population)
tsplt

redistrib <- data.frame(Population = unique(counts$Population),
                        date_redistrib = c("2018-07-10", "2018-07-09",
                                           NA, "2018-07-07", "2018-07-05",
                                           "2018-07-09", "2018-07-06")) %>% 
  mutate(date_redistrib = as.Date(date_redistrib))

f1counts <- left_join(counts, redistrib) %>% 
  filter(date <= date_redistrib) %>% 
  group_by(Population, Pool) %>% 
  summarise(F1counts = sum(counts, na.rm = TRUE),
            mF1counts = max(counts, na.rm = TRUE))

# F2 adults removed as counted
dat <- readxl::read_xlsx("data/CommonGarden_F2_AdultRemoval.xlsx", na = "NA")
names(dat)[3:9] <- paste0("x", names(dat[3:9]))
removal <- dat %>% 
  gather(key = date, value = counts, x43326:x43371) %>% 
  mutate(date = as.Date(as.numeric(sub(pattern = "x", replacement = "", x = date, fixed = TRUE)),
                        origin = "1899-12-30")) %>% 
  filter(Pool != "ALL") %>% 
  mutate(Pool = as.numeric(as.character(Pool))) %>% 
  group_by(Population, Pool) %>% 
  summarise(F2counts = sum(counts, na.rm = TRUE))

# F1 diapause proportion for each population
dat <- readxl::read_xlsx("data/CommonGarden_F1_Diapause.xlsx", na = "NA")
diapause <- dat %>% 
  filter(Pool != "ALL") %>% 
  mutate(Pool = as.numeric(as.character(Pool)),
         Total = Diapause + Feeder + Repro_female,
         Perc_diap = Diapause / Total,
         PotID = paste(Population, Pool, sep = "_"),
         Date = as.Date(Date)) %>% 
  filter(Total > 0) %>% 
  left_join(f1counts) %>% 
  left_join(controlvars) %>% 
  left_join(gdd, by = c("Date" = "date")) %>% 
  mutate(elapsed_degdays = accumdd - StartAccumDD,
         zF1counts = scale(F1counts),
         zPlantSize = scale(TotalStemHeight),
         zAphids = scale(Aphids))

# Growth chamber photoperiod
cp <- data.frame(Population = unique(diapause$Population),
                 CP = c(15.56, 16.07, 15.81, 15.12, 15.94, 16.2)) %>% 
  mutate(zCP = scale(CP))
diapause <- left_join(diapause, cp)

moddiap <- glmer(Perc_diap ~ zF1counts 
                 + zPlantSize
                 + zAphids 
                 + zCP 
                   + (1 | Population)
                 + (1|Pool), weights = Total, 
                 family = binomial(link = "logit"), data = diapause)
summary(moddiap)

# table of model output F1 diapause
write.csv(broom.mixed::tidy(moddiap), file = "F1diapause_glmer.csv", row.names = FALSE)


popdiap <- ranef(moddiap)$Population
popdiap$Population <- row.names(popdiap)
names(popdiap)[1] <- "Diapause"

preddf <- cp %>% 
  mutate(zF1counts = 0,
         zPlantSize = 0,
         zAphids = 0,
         Pool = NA)
  
preddf$pred <- predict(moddiap, newdata = preddf, re.form = ~(1 | Population), type = "response")

boxplt <- ggplot(diapause, aes(x = Population, y = Perc_diap)) +
  geom_point(aes(size = log(Total)), alpha = 0.3) + 
  geom_point(data = preddf, aes(x = Population, y = pred), color = "red", shape = 17, size = 3) +
  ylab("Percent entering diapause") +
  ggtitle("F1 generation: diapause by population", subtitle = "Observed data and modeled population averages")
boxplt
# F2 diapause proportion for each population
dat <- readxl::read_xlsx("data/CommonGarden_F2_Diapause.xlsx", na = "NA")
diapause <- dat %>% 
  filter(Pool != "ALL") %>% 
  mutate(Pool = as.numeric(as.character(Pool)),
         Total = Diapause + Feeder + Repro_female,
         Perc_diap = Diapause / Total,
         PotID = paste(Population, Pool, sep = "_"),
         Date = as.Date(Date_test)) %>% 
  filter(Total > 0) %>% 
  left_join(removal) %>% 
  left_join(controlvars) %>% 
  left_join(gdd, by = c("Date" = "date")) %>% 
  mutate(elapsed_degdays = accumdd - StartAccumDD,
         zF2counts = scale(F2counts),
         zPlantSize = scale(TotalStemHeight),
         zAphids = scale(Aphids),
         zSenesc = scale(Senescence))
diapause <- left_join(diapause, cp)

photo_exp <- diapause %>% 
  dplyr::select(Population:Collected_Aug22) %>% 
  group_by(Population, Pool, Date_test) %>% 
  summarise(photo = weighted.mean(c(14.29, 14.10, 13.92), w = c(Collected_Aug14, Collected_Aug18, Collected_Aug22)),
            density = Collected_Aug14 + Collected_Aug18 + Collected_Aug22) %>% 
  ungroup() %>% 
  mutate(   zphoto = scale(photo)[,1],
            zdens = scale(density)[,1])

diapause <- left_join(diapause, photo_exp) %>% 
  left_join(f1counts) %>% 
  mutate(zF1counts = scale(F1counts))


moddiap <- glmer(Perc_diap ~ zF2counts
                 + zPlantSize
                 + zCP
                 + zphoto
                 # + zdens
                 + zAphids 
                 + zSenesc
                 + (1 | Population)
                 + (1|Pool), weights = Total, 
                 family = binomial(link = "logit"), data = diapause)
summary(moddiap)

write.csv(broom.mixed::tidy(moddiap), file = "F2diapause_glmer.csv", row.names = FALSE)



popdiap <- ranef(moddiap)$Population
popdiap$Population <- row.names(popdiap)
names(popdiap)[1] <- "Diapause"

preddf <- cp %>% 
  mutate(zF2counts = 0,
         zPlantSize = 0,
         zAphids = 0,
         zphoto = 0,
         zSenesc = 0,
         Pool = NA)

preddf <- data.frame(zF2counts = 0,
                     zPlantSize = 0,
                     zCP = 0,
                     zphoto = 0,
                     zSenesc = 0,
                     zAphids = 0,
                     Population = c("BL", "BS", "M", "S", "V", "Y"),
                     Pool = NA)
preddf$pred <- predict(moddiap, newdata = preddf, re.form = ~(1 | Population), type = "response")

boxplt <- ggplot(diapause, aes(x = Population, y = Perc_diap)) +
  geom_point(aes(size = log(Total)), alpha = 0.3) + 
  geom_point(data = preddf, aes(x = Population, y = pred), color = "red", shape = 17, size = 3) +
  ylab("Percent entering diapause") +
  ggtitle("F2 generation: diapause by population", subtitle = "Observed data and modeled population averages")
boxplt


# can diapause differences be partially explained by temperature?
diaptemp <- gdd %>% 
  filter(month == 7, day <= 14) %>% 
  summarise(meantemp = mean((tmax + tmin)/2))

# 
latef1counts <- left_join(counts, redistrib) %>% 
  filter(date > date_redistrib) %>% 
  group_by(Population, Pool) %>% 
  summarise(lateF1counts = sum(counts, na.rm = TRUE))



moddat <- controlvars %>% 
  left_join(flowermass) %>% 
  left_join(f1counts) %>%
  left_join(latef1counts) %>% 
  left_join(removal) %>% 
  # left_join(popdiap) %>% 
  left_join(rootmass) %>% 
  left_join(cp) %>% 
  ungroup() %>% 
  mutate(F1counts = replace_na(F1counts, 0),
         F2counts = replace_na(F2counts, 0)) %>% 
  mutate(PotID = paste(Population, Pool, sep = "_"),
         TotalCounts = F1counts + F2counts,
         HerbivTiming = F1counts / TotalCounts,
         zF1 = scale(F1counts),
         zF2 = scale(F2counts),
         logF1 = log(F1counts + 1),
         logF2 = log(F2counts + 1),
         Pool = as.factor(Pool),
         zPlantSize = scale(TotalStemHeight),
         zAphids = scale(Aphids),
         zSenesc = scale(Senescence),
         zRegrowth = scale(Regrowth),
         zBeetles = log(TotalCounts + 1))
moddat$season_total[is.na(moddat$season_total)] <- 0
moddat$Control <- ifelse(moddat$Population == "C", "yes", "no")
moddat$logflower <- log(moddat$season_total + 1)

pairs(moddat[, c("season_total", "zF1", "zF2", "Aphids", "Senescence", "Regrowth")])

# population growth rates
# Do these make sense with the redistribution we did after F1???
# /3 is for f1counts calculated by sum instead of max (3 times more)
moddat <- moddat %>% filter(Population != "C", F1counts > 0, F2counts > 0)
lam1 <- log((moddat$F1counts/3) / 8)
lam2 <- log(moddat$F2counts / (moddat$lateF1counts/3))
ggplot(moddat, aes(x = log(F1counts / 8 / 3), y = log(F2counts / (lateF1counts/3)), 
                   group = Population, color = Population)) + 
  geom_point(size = 2.5, alpha = 0.5) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  xlab('P to F1 population growth rate (log)') +
  ylab('F1 to F2 population growth rate (log)') +
  ggtitle('Galerucella population growth rates')

cg_dag <- dagify(Fitness ~ Flowers + Roots,
                 Roots ~ Herbivory + PlantSize,
                 Flowers ~ Herbivory + PlantSize,
                 Herbivory ~ Aphids + F1 + F2,
                 Aphids ~~ F1,
                 F1 ~ Population + PlantSize,
                 F2 ~ F1 + Diapause,
                 Diapause ~ Population + F1 + Aphids + PlantSize,
                 labels = c("Herbivory" = "Herbivory",
                            "PlantSize" = "Plant Size",
                            "Aphids" = "Aphids",
                            "Beetles" = "Beetles",
                            "F1" = "F1",
                            "F2" = "F2",
                            "Population" = "Beetle\nPopulation",
                            "Diapause" = "Diapause",
                            "Fitness" = "Fitness",
                            "Reproduction" = "Reproduction",
                            "Roots" = "Roots",
                            "Flowers" = "Reproduction"),
                 # latent = "unhealthy",
                 # exposure = "smoking",
                 outcome = "Fitness")

ggdag(cg_dag, text = FALSE, use_labels = "label")

cg_dag <- dagify(Fitness ~ Herbivory,
                 Herbivory ~ F1 + F2,
                 F2 ~ F1 + Diapause,
                 F1 ~ Population,
                 Diapause ~ Population,
                 labels = c("Herbivory" = "Herbivory",
                            "PlantSize" = "Plant Size",
                            "Aphids" = "Aphids",
                            "Beetles" = "Beetles",
                            "F1" = "F1",
                            "F2" = "F2",
                            "Population" = "Beetle\nPopulation",
                            "Diapause" = "Diapause",
                            "Fitness" = "Plant\nFitness",
                            "Reproduction" = "Reproduction",
                            "Roots" = "Roots",
                            "Flowers" = "Reproduction"),
                 # latent = "unhealthy",
                 # exposure = "smoking",
                 outcome = "Fitness")

ggdag_classic(cg_dag, text_label = "label")+ theme_dag_blank()



# control have more flowers
mod <- glmer(round(season_total) ~ 
              # Control +
              zPlantSize +
              zAphids +
               zBeetles +
              (1|Pool) + (1|Population), 
            family = poisson(link = "log"),
            data = moddat)
summary(mod)


write.csv(broom.mixed::tidy(mod), file = "flower_beetles_glmer.csv", row.names = FALSE)
MuMIn::r.squaredGLMM(mod)

# and roots?
mod <- glmer(round(Root_mass) ~ 
               # Control +
               zPlantSize +
               zAphids +
               zBeetles +
               # zCP +
              zSenesc +
              zRegrowth +
               (1|Population) +
              (1|Pool), 
             family = poisson(link = "log"),
             data = moddat)
summary(mod)

write.csv(broom.mixed::tidy(mod), file = "root_allvars_glmer.csv", row.names = FALSE)
MuMIn::r.squaredGLMM(mod)


mod <- glm(round(Root_mass) ~ Population +
             log(TotalStemHeight) +
             Aphids +
             zF1 + zF2 + Senescence + Regrowth, 
           family = poisson(link = "log"),
           data = moddat)
summary(mod)


boxplt <- ggplot(moddat, aes(x = Population, y = Senescence)) +
  geom_boxplot()
boxplt


# alternatively, 



trtdat <- moddat %>% filter(Population != 'C') %>% filter(Population != "V")
mod <- lmer(logflower ~ 
              (TotalCounts + HerbivTiming)^2 +
              (1|Pool), 
            # family = poisson(link = "log"), 
            data = trtdat)
summary(mod)
mod <- lmer(logflower ~ 
              (zF1 + zF2)^2 +
              (1|Population), 
            # family = poisson(link = "log"), 
            data = trtdat)
mod <- lmer(logflower ~ 
              (logF1 + logF2)^2 +
              (1|Population), 
            # family = poisson(link = "log"), 
            data = trtdat)
summary(mod)

# what about quasipoisson?
mod <- glm(season_total ~ 
             Population - 1 +
             TotalCounts * HerbivTiming, 
           family = quasipoisson(link = "log"),
           data = trtdat)
summary(mod)

ggplot(trtdat, aes(x = Diapause, y = HerbivTiming)) +
  geom_point()

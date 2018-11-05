# preliminary analysis

library(tidyverse)
library(lubridate)
library(viridis)
library(readxl)
library(lme4)
theme_set(theme_bw(base_size = 20)) 


dat <- readxl::read_xlsx("data/CommonGarden_Flowers.xlsx", na = "NA")
flowermass <- dat %>% 
  group_by(Population, Pool) %>% 
  summarise(season_total = sum(FlowerMass, na.rm = TRUE))

# basic plot of mean flower mass by population
boxplt <- ggplot(flowermass, aes(x = Population, y = season_total)) +
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
  summarise(F1counts = sum(counts, na.rm = TRUE))

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
         PotID = paste(Population, Pool, sep = "_")) %>% 
  filter(Total > 0)

moddiap <- glmer(Perc_diap ~ 1 + (1|Population) + (1|Pool) + (1|PotID), weights = Total, 
                 family = binomial(link = "logit"), data = diapause)

popdiap <- ranef(moddiap)$Population
popdiap$Population <- row.names(popdiap)
names(popdiap)[1] <- "Diapause"

moddat <- controlvars %>% 
  left_join(flowermass) %>% 
  left_join(f1counts) %>% 
  left_join(removal) %>% 
  left_join(popdiap) %>% 
  ungroup() %>% 
  mutate(PotID = paste(Population, Pool, sep = "_"),
         TotalCounts = F1counts + F2counts,
         HerbivTiming = F1counts / TotalCounts,
         zF1 = scale(F1counts),
         zF2 = scale(F2counts),
         logF1 = log(F1counts + 1),
         logF2 = log(F2counts + 1),
         Pool = as.factor(Pool))
moddat$season_total[is.na(moddat$season_total)] <- 0
moddat$Control <- ifelse(moddat$Population == "C", "yes", "no")
moddat$logflower <- log(moddat$season_total + 1)

# control have more flowers
mod <- lmer(logflower ~ Control + 
              TotalStemHeight +
              Aphids +
              (1|Pool) + (1|Population), 
            # family = poisson(link = "log"), 
            data = moddat)
summary(mod)

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


ggplot(trtdat, aes(x = Diapause, y = HerbivTiming)) +
  geom_point()
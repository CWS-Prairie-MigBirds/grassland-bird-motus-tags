#### GBMP Script 02 - Data cleaning ####
# GBMP = Grassland Bird Monitoring Program
# ECCC contract 3000792304
# adapted from the Motus R book by Sarah Endenburg
# This script contributes to deliverable 1

#### 1 - Prepare working environment ####

library(motus)
library(tidyverse)
library(lubridate)
library(ggforce)
library(stringr)
library(rnaturalearth)
library(geosphere)
library(osfr) # for interacting with the OSF repository I created to store the GBMP Motus data

# Motus data is stored in UTC, set system time to match
Sys.setenv(tz = "UTC")

#### 2 - loading data ####

# you can load my file from OSF or you can load the data from your working directory after running through script 01

# Authenticate OSF account to download/upload data
# paste your OSF PAT within the quotation marks below to authenticate 
osf_auth(token = "")

project <- osf_retrieve_node("ge8bd")
# the code in the quotation marks is the project id on OSF

files <- osf_ls_files(project)

file_to_download <- files %>%
  filter(name == "df.alltags.912.filter.nodups.2025Jan20.rds")
# retrieving the file object (just using file name doesn't work, there might be multiples of the same name)

osf_download(file_to_download, conflicts = "rename")

# now you can load it into R
df.alltags.912.filter.nodups <- readRDS("df.alltags.912.filter.nodups.2025Jan20.rds")

#### 3 - initial track map ####

# to see what we're working with

world <- ne_countries(scale = "medium", returnclass = "sf")

states_CanUS <- ne_states(country = c("canada", "united states of america"), returnclass = "sf")

lakes <- ne_download(scale = "large", type = 'lakes', category = 'physical',
returnclass = "sf", destdir = "mapData") # only need this first time downloading

#lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
                 #returnclass = "sf",
                 #destdir = "mapData") # use this if already downloaded shapefiles

fun.getpath <- function(df) {
  df %>%
    filter(!is.na(recvDeployLat) | !(recvDeployLat == 0)) %>% 
    group_by(motusTagID, mfgID, runID, recvDeployName, ambigID, 
             tagDepLon, tagDepLat, recvDeployLat, recvDeployLon, speciesEN) %>%
    summarize(max.runLen = max(runLen), time = mean(ts), .groups = "drop") %>%
    arrange(motusTagID, time) %>%
    data.frame()
}

path912 <- fun.getpath(df.alltags.912.filter.nodups)

xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1

ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = states_CanUS, colour = "black") +
  geom_sf(data = lakes, colour = NA, fill = "white") +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  labs(x = "", y = "", colour = "Species", fill = "Species") +
  geom_path(data = filter(path912, is.na(speciesEN) == FALSE), 
            aes(x = recvDeployLon, y = recvDeployLat, 
                colour = as.factor(speciesEN), group = motusTagID)) +
  geom_point(data = filter(path912, is.na(speciesEN) == FALSE), aes(x = recvDeployLon, y = recvDeployLat), fill = "black", shape = 21) +
  geom_point(data = filter(path912, is.na(speciesEN) == FALSE), 
             aes(x = tagDepLon, y = tagDepLat), fill = "red", colour = "black", shape = 23, size = 2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  facet_wrap(~speciesEN)
# yeah this looks crazy.

#### 4 - noisy tower filter for lotek tags ####

# this is a filter that I developed for my BANS tags - it calculates the proportion of hits, per individual Motus tag, per receiver, per day, that have a freqsd >0.1. High freqsd is an indicator of false detections, and I found that a cutoff of a proportion of noisy hits >0.25 removes the vast majority of detections that are suspicious (as confirmed by looking at detections for individual tags and signal strength plots, etc.)
# it only works for Lotek tags though, there's no freqsd for CTT tags. but lotek tags are much more prone to false detections anyways

lotekTags <- df.alltags.912.filter.nodups |>
  filter(mfg == "Lotek")

df.prop.fail.freq.filt.912 <- lotekTags |>
  mutate(freqFilter = ifelse(freqsd > 0.1, 0, 1),
         date = as.Date(ts)) |>
  group_by(motusTagID, recvDeployName, date, freqFilter) %>% # how many detections were flagged as false for that tag at that receiver on that day
  count() %>% 
  pivot_wider(names_from = freqFilter, values_from = n, values_fill = 0) %>%
  rename(passedFreqFilter = `1`,
         failedFreqFilter = `0`) %>%
  mutate(propFailedFreqFilter = round(failedFreqFilter/
                                        (failedFreqFilter + passedFreqFilter), 
                                      digits = 2))

# to view by receiver (16 per page):

df.prop.fail.freq.filt.912 %>% 
  ggplot(aes(x = date, y = propFailedFreqFilter)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  facet_wrap_paginate(~recvDeployName, # plot subset of data per page
                      scales = "free_x",
                      ncol = 4, nrow = 4,
                      page = 1)
# change the page = # from 1-8 in the last line to view all pages

# this filter catches a LOT of hits that I think are unlikely - many are in Eastern Canada (Maritimes, Atlantic coast of the USA), South America, receivers that I know are noisy from my BANS work (FDSHQ, Reserva de la Biosfera Ria Celestun)
# will verify below before removing

# need a unique code to ID individual tag-receiver-date combinations in each dataframe. they look like this: 51898-Nashlyn-2022-06-13
df.prop.fail.freq.filt.912 <- df.prop.fail.freq.filt.912 |>
  mutate(dateTagRecvID = as.character(paste(motusTagID, recvDeployName, date, sep = "-")))

df.alltags.912.filter.nodups <- df.alltags.912.filter.nodups |>
  mutate(date = as.Date(ts)) |>
  mutate(dateTagRecvID = as.character(paste(motusTagID, recvDeployName, date, sep = "-")))

failFreqFilter912 <- df.prop.fail.freq.filt.912 |>
  filter(propFailedFreqFilter >= 0.25)

# and adding the filter column
df.alltags.912.filter.nodups <- df.alltags.912.filter.nodups |>
  mutate(noisyTowerFilter = ifelse(dateTagRecvID %in% failFreqFilter912$dateTagRecvID, 0, 1))

table(df.alltags.912.filter.nodups$noisyTowerFilter)
# 5578 hits fail, 419041 pass. so it's a small proportion, only 1.3% of hits, but it makes a big difference to the track maps.

noisyTowerHits <- df.alltags.912.filter.nodups |>
  filter(noisyTowerFilter == 0) |>
  mutate(reasonRemoved = "noisyTowerHits")

#### 5 - Baird's Sparrow ####

# I like to take a look at the "unfiltered" detections to get an idea of where each individual/species was detected so that I can more easily tell if something looks suspicious/needs a more in-depth inspection

BAIS <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Baird's Sparrow")

BAIStags <- unique(BAIS$motusTagID) # 17

for(i in 1:length(BAIStags)){
  print(ggplot(data = filter(BAIS, 
                             motusTagID == BAIStags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName, shape = as.factor(noisyTowerFilter))) +
          labs(title = paste(BAIStags[i], "BAIS", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# detections that I think look suspicious/require further investigation:

# 74912 was only detected at kennekuk3 and alaksen, check those. 
# 61608 also looks suspicious - clear lake and fletcher pass the noisy tower filter
# 67371 looks mostly good, check FDSHQ and Finca Cristina - short run lengths, towers that I know are noisy
# 64536 - Reserva de la Biosfera and white rock are the only detections, white rock might be ok
# 64538 also possibly has some false detections at shuswap (in BC)
# 64539 at reserva de la biosfera
# 64537 only at FDSHQ and reserva de la biosfera, all runs of three
# 74904 - APR - schmoekel is probably good. short runs at audubon center, hacienda la esperanza
# 85295 only at kennekuk 4&6
# 74905 only at Bowdoin NWR and Alaksen
# 74911 at facultad de ciencias agropecuarias (noisy tower), Miles city, alaksen
# 74906 only a run of three at alaksen
# 85292 - Kennekuk 6
# 85291 only a run of 3 at kennekuk 6

# one BAIS was detected in michigan, which is quite far east

eastBais <- df.alltags.912.clean |>
  filter(speciesEN == "Baird's Sparrow" & recvDeployLon > -90)
# clear lake & fletcher, two runs of three. 61608
# lotek tag and the freq, freqsd and sig look odd. all perfectly round numbers, freqsd is zero 

#### 6 - Sprague's pipit ####

# SPPI has a mix of CTT and Lotek tags

SPPI <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Sprague's Pipit")

SPPItags <- unique(SPPI$motusTagID) # 56

for(i in 1:length(SPPItags)){
  print(ggplot(data = filter(SPPI, 
                             motusTagID == SPPItags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(SPPItags[i], "SPPI", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

table(SPPI$tagDepYear)
# detections to verify for SPPI:

# 67509 has what appears to be some noise - st-jean-baptiste? passes the noisy tower filter but its in Quebec...
# 67506 at para la tierra - fails noisy tower filter
# 61006 - FDSHQ and kejimkujik are sketchy.
# 61014 at reserva de la biosfera, kejimkujik, tadoussac will all probably be removed by noisy tower filter
# 61012 at kejimkujik too
# 51879 only detected in April, spring migration? its a CTT tag. all look plausible
# 75410 at summit fire tower - will hopefully get removed by noisy tower filter. kennekuk 6 - noisy filter looks fine, check freq. kennekuk recv are all fairly far east
# 61020 at kennekuk 6 - timing looks off, was detected on breeding grounds and then back south again in July?
# 75401 at quilchena and nashlyn will probably get removed, kennekuk 6?
# 75402 apparently at Quilchena in BC at the same time as at ellice archie in SK/MB
# 75398 - mcFarlane park, Tadoussac, RPBO-Donnecke, FERS, Campbellville, Bose Elementary School, Scout Islands
# 85080 at kennekuk 6, 2024_D4_Hou
# 85765 at kennekuk 6, before detected at hanleywest for almost a month

# there's one SPPI that was detected in southern QC, which is quite far east

eastSPPI <- df.alltags.912.clean |>
  filter(speciesEN == "Sprague's Pipit" & recvDeployLon > -90)
# 67509 - run of three at St-Jean-Baptiste on June 10

#### 7 - Horned Lark ####

HOLA <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Horned Lark")

HOLAtags <- unique(HOLA$motusTagID) # 21

for(i in 1:length(HOLAtags)){
  print(ggplot(data = filter(HOLA, 
                             motusTagID == HOLAtags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(HOLAtags[i], "HOLA", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# hits to verify/that look suspicious:

# 67510 at FL panther, tadoussac, reserva de la biosfera - could remove pre-departure hits and that might work. or flight speed. or restrict longitude to get rid of tadoussac
# 67510 looks insane. I think it was "detected" on every noisy tower possible.
# 75411 at kennekuk 6
# 64534 all look false - high park, kejimkujik and mary's garden. all east coast
# 85089, 85086 & 85085 only have runs of 3 at Kennekuk 6
# 75827 too?
# 67510 at GFAFB - Fire station (location ok?), florida panther

# one HOLA looks like it went down to FL and then back east...67510 goes to FL panther and then GFAFB?

#### 8 - Chestnut-collared Longspur ####

CCLO <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Chestnut-collared Longspur")

CCLOtags <- unique(CCLO$motusTagID) # 67

for(i in 1:length(CCLOtags)){
  print(ggplot(data = filter(CCLO, 
                             motusTagID == CCLOtags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(CCLOtags[i], "CCLO", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

ggplot(data = filter(df.alltags.912.clean, 
                     motusTagID == 67388), 
       aes(x = ts, y = recvDeployLat)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  geom_point(aes(size = runLen, color = recvDeployName)) +
  labs(title = paste(67388, "CCLO", sep = "-"))

# sketchy detections to look into:

# 74909 at kennekuk 3&5&6, facultad de ciencias agropecuarias, alaksen?? most alaksen hits fail the noisy tower filter
# 67508 at kejimkujik
# 64600 at FDSHQ
# 61609 is very messy - detected at a lot of noisy stations, I think they'll mostly be filtered out by the noisy tower filter
# 67505 looks like it was dropped - continuously detected at battlecreek. unless the bird stayed there. check sig strength
# 74907 at kent island? far east
# 61606 at kennekuk 3
# 61008 at magnetic hill zoo
# 74916 at alaksen in december?
# 61018 has short runs at kejimkujik...
# 61021 at kennekuk 6
# 67384 only detected at FDSHQ, runs of 3
# 74910 at alaksen in November? 74914 too, only detected there. 74913 too. 
# 85822 was at nashlyn until october - seems ok according to birds of the world
# 85813 at kennekuk 6 - range shape file would help with that
# 85752 at kennekuk 6
# 85821 only detected at GPCD SPLT. looks good though...

#### 9 - Thick-billed longspur ####

TBLO <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Thick-billed Longspur")

TBLOtags <- unique(TBLO$motusTagID) # 8

for(i in 1:length(TBLOtags)){
  print(ggplot(data = filter(TBLO, 
                             motusTagID == TBLOtags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(TBLOtags[i], "TBLO", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# detections to verify:

# 61022 & 85082 at kennekuk 6

#### 10 - Grasshopper sparrow ####

GRSP <- df.alltags.912.filter.nodups |>
  filter(speciesEN == "Grasshopper Sparrow")

GRSPtags <- unique(GRSP$motusTagID) # 2

for(i in 1:length(GRSPtags)){
  print(ggplot(data = filter(GRSP, 
                             motusTagID == GRSPtags[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(GRSPtags[i], "GRSP", sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# all runs of three, kennekuk6, point loma, battlecreek and snake butte (the latter two are probably fine) 

# 11 - looking at individual receiver signal strength plots ####

reservaLaBiosfera <- df.alltags.912.clean |>
  filter(grepl("Reserva", recvDeployName))

ggplot(data = filter(df.alltags.912.filter.nodups,
                     recvDeployName == "Battlecreek"), 
       aes(x = ts, y = sig)) +
  theme_bw() +
  geom_point(aes(size = runLen, color = as.factor(motusTagID))) +
  labs(title = "Battlecreek", x = "ts", y = "sig")

# 67505 does look like it was dropped near battlecreek - looks good until maybe late August and then the signal strength stays pretty constant, detected throughout the winter.

kennekuk <- df.alltags.912.filter.nodups %>%
  filter(str_detect(recvDeployName, regex("kennekuk", ignore_case = TRUE)))
# kennekuk max run length is 4 (at any kennekuk receiver)
# yeah I don't like these ones - the sigsd and freqsd are 0 exactly for every hit, and the freq are all perfectly round numbers. it's sketchy.

FLpanther <- df.alltags.912.clean |>
  filter(recvDeployName == "Florida Panther NWR, FL") |>
  select(ts, motusTagID, runLen, freq, freqsd, recvDeployName, speciesEN)

table(FLpanther$runLen) # mostly runs of three, some runs of 4 & 5, one 6, 7, 8
# could get rid of short runs since it's a noisy tower. but runs of 7-8? they're not for HOLA, they're for CCLO which I don't think should be in FL...

#### 12 - applying the noisy tower filter ####

df.alltags.912.clean <- df.alltags.912.filter.nodups |>
  filter(noisyTowerFilter == 1)

#### 13 - filtering out detections that are outside of the expected range of each species ####

# I wanted to download an ebird raster or something to do this, but I don't think it would work super well - eg. for BAIS, which has discrete breeding and wintering distributions, it would've filtered out migratory hits... 
# instead, I'm just applying min and max lat and long based roughly on the range maps from Birds of the World
# eBird does show some detections outside of those ranges, so it is theoretically possible. Most of these receivers I've had issues with before anyways though so I am skeptical. will set aside for future investigation

df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(minRecvLon = case_when(
    speciesEN == "Baird's Sparrow" ~ as.numeric("-115"),
    speciesEN == "Chestnut-collared Longspur" ~ as.numeric("-115"),
    speciesEN == "Horned Lark" ~ as.numeric("-160"),
    speciesEN == "Thick-billed Longspur" ~ as.numeric("-115"),
    speciesEN == "Sprague's Pipit" ~ as.numeric("-115"),
    speciesEN == "Grasshopper Sparrow" ~ as.numeric("-160"),
    TRUE ~ as.numeric("0")),
    
  maxRecvLon = case_when(
    speciesEN == "Baird's Sparrow" ~ as.numeric("-90"),
    speciesEN == "Chestnut-collared Longspur" ~ as.numeric("-92"),
    speciesEN == "Horned Lark" ~ as.numeric("-87"), # technically, HOLA range extends further east, but this just removed a FL detection that seemed unlikely anyways
    speciesEN == "Thick-billed Longspur" ~ as.numeric("-100"),
    speciesEN == "Sprague's Pipit" ~ as.numeric("-90"),
    speciesEN == "Grasshopper Sparrow" ~ as.numeric("-60"),
    TRUE ~ as.numeric("0")),
  
  minRecvLat = case_when(
    speciesEN == "Baird's Sparrow" ~ as.numeric("18"),
    speciesEN == "Chestnut-collared Longspur" ~ as.numeric("18"),
    speciesEN == "Horned Lark" ~ as.numeric("18"),
    speciesEN == "Thick-billed Longspur" ~ as.numeric("23"),
    speciesEN == "Sprague's Pipit" ~ as.numeric("14"),
    speciesEN == "Grasshopper Sparrow" ~ as.numeric("9"),
    TRUE ~ as.numeric("0")),

  maxRecvLat = case_when(
    speciesEN == "Baird's Sparrow" ~ as.numeric("55"),
    speciesEN == "Chestnut-collared Longspur" ~ as.numeric("55"),
    speciesEN == "Horned Lark" ~ as.numeric("60"),
    speciesEN == "Thick-billed Longspur" ~ as.numeric("53"),
    speciesEN == "Sprague's Pipit" ~ as.numeric("56"),
    speciesEN == "Grasshopper Sparrow" ~ as.numeric("55"),
    TRUE ~ as.numeric("0")))

hitsOutOfRange <- df.alltags.912.clean |>
  filter(recvDeployLon < minRecvLon | recvDeployLon > maxRecvLon | recvDeployLat < minRecvLat | recvDeployLat > maxRecvLat) |>
  mutate(reasonRemoved = "outOfSpeciesRange")
# 534

table(hitsOutOfRange$recvDeployName)
# all recv I've had issues with

df.alltags.912.clean <- df.alltags.912.clean |>
  filter(!(recvDeployLon < minRecvLon | recvDeployLon > maxRecvLon | recvDeployLat < minRecvLat | recvDeployLat > maxRecvLat))

#### 14 - removing individual tags, detections with strange freq and freqsd ####

tagsToRemove <- c(67505, 75839)
# 67505 looks like it was dropped near the Battlecreek receiver - detections are pretty well continuous at the same receiver for extended periods of time
# 75839-HOLA does look like it was dropped in mid-late august, hence the detections at Govenlock until December

droppedTags <- df.alltags.912.clean |>
  filter(motusTagID %in% c(67505, 75839)) |>
  mutate(reasonRemoved = "droppedTag")

df.alltags.912.clean <- df.alltags.912.clean |>
  filter(!motusTagID %in% c(67505, 75839))

# and some towers have very odd freq and freqsd - they'll be perfectly even numbers and the freqsd is exactly zero. that combination never happens at any receiver that I trust (sometimes freqsd will be zero but freq should always have some decimals)
weirdTowers <- df.alltags.912.clean |>
  filter((freqsd == 0 & freq %% 1 == 0)) |>
  mutate(reasonRemoved = "weirdTower")

table(weirdTowers$recvDeployName) # 33 hits
# kennekuk 6 and mary's garden. I've never had any good hits there, they're both more eastern, especially mary's garden

df.alltags.912.clean <- df.alltags.912.clean |>
  filter(!(freqsd == 0 & freq %% 1 == 0 & tagModel %in% c("NTQB2-2-M", "NTQB2-3-2-M", "NTQB2-5-1-M"))) # have to specify that this is only for lotek tags or it'll remove hits for ctt tags too

#### 15 - removing short runs at receivers with a high proportion of noisy hits ####

# even after applying the noisy tower filter, some hits at very noisy towers still passed through.

noisyTowers <- df.alltags.912.filter.nodups %>%
  group_by(recvDeployName) %>% 
  summarise(failure_percentage = mean(noisyTowerFilter == 0), .groups = "drop") %>%  # Calculate the percentage of failures (noisyTowerFilter == 0)
  filter(failure_percentage >= 0.5)

moreNoisyHits <- df.alltags.912.clean |>
  filter(recvDeployName %in% noisyTowers$recvDeployName & runLen <5 & tagModel %in% c("NTQB2-2-M", "NTQB2-3-2-M", "NTQB2-5-1-M")) |>
  mutate(reasonRemoved = "shortRunNoisyRecv")
#233, only Lotek tags detected there which is interesting.

df.alltags.912.clean <- df.alltags.912.clean |>
  filter(!(recvDeployName %in% noisyTowers$recvDeployName & runLen <5 & tagModel %in% c("NTQB2-2-M", "NTQB2-3-2-M", "NTQB2-5-1-M")))

#### 15 - burst interval filter ####

# calculates the time elapsed between consecutive hits in a run and flags any that deviate substantially from a multiple of the tag burst interval. Motus does allow for some missed hits before starting a new run, and some slight deviation is to be expected, which is why I don't expect hits to be exact multiples of the burst interval

# I do this in two separate steps because it's a bit computationally expensive
df.alltags.912.clean <- df.alltags.912.clean %>%
  arrange(runID, ts) %>%
  group_by(runID) %>%
  mutate(timeDiff = ts - lag(ts))

df.alltags.912.clean <- df.alltags.912.clean %>%
  arrange(runID, ts) %>%
  group_by(runID) %>%
  mutate(timeDiff_divBy_BI = as.numeric(timeDiff/tagBI))

hist(df.alltags.912.clean$timeDiff_divBy_BI)
# a very small number appear to be >120 times the BI, which seems high

highTimeDiff <- df.alltags.912.clean |>
  filter(timeDiff_divBy_BI > 20)
# 1641 hits

table(highTimeDiff$recvDeployName) # all seem like plausible locations

table(highTimeDiff$tagModel)
# they're all for CTT tags. I've messaged Motus to try and better understand how runs work for CTT tags, apparently it's different than for Lotek tags since burst interval isn't required to verify the tag ID of CTT tags. with lotek tags, the max # missing hits before a new run is started seems to be about 60, which is still high. but it seems to be higher for CTT tags. I'll leave them in until I find out more information

hist(highTimeDiff$timeDiff_divBy_BI)

BI_flag <- df.alltags.912.clean %>%
  filter(abs(timeDiff_divBy_BI - round(timeDiff_divBy_BI)) > 0.2)

table(BI_flag$tagModel) # all CTT tags... maybe there is more wiggle room for BI with CTT tags than Lotek? will investigate and come back to this if it seems suspicious

table(BI_flag$recvDeployName) # the majority are at nashlyn, which is a very probable location for real detections. they all seem plausible really

#### 16 - travel speed filter ####

# this filter flags detections that would require biologically unfeasible travel speeds for the individual to actually move between those receivers

fun_distTime <- function(df) {
  df %>% 
    ungroup() |>
    mutate(distBtwRecvM = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName)),
                                  c(NA, distVincentyEllipsoid(df[,c("recvDeployLon", "recvDeployLat")])),
                                  NA_real_),
           distBtwRecvKM = distBtwRecvM/1000,
           timeH = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName) & distBtwRecvKM >30 ),
                           as.numeric(difftime(ts, lag(ts), units = "hours")),
                           NA_real_),
           timeD = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName) & distBtwRecvKM > 30),
                           as.numeric(difftime(ts, lag(ts), units = "days")),
                           NA_real_),
           timeS = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName) & distBtwRecvKM > 30),
                           as.numeric(difftime(ts, lag(ts), units = "secs")),
                           NA_real_),
           gspeedMS = distBtwRecvM/as.numeric(timeS),
           gspeedKMH = 3.6*gspeedMS)
}

df.alltags.912.clean <- df.alltags.912.clean |>
  ungroup()|>
  arrange(motusTagID, ts)

df.alltags.912.clean <- fun_distTime(df.alltags.912.clean)

# Filter high-speed detections (gspeedMS > 100)
high_speed_individuals <- df.alltags.912.clean %>%
  filter(gspeedMS > 100 | is.infinite(gspeedMS)) %>%  # Detections with high flight speeds
  select(motusTagID, runID, hitID, recvDeployName, gspeedMS) %>%  # Keep only necessary columns
  distinct()  # Get distinct combinations of tag, run, and receiver

# Join this filtered data back to the original dataframe to get all detections for those individuals, so that I can look at them in order of ts
all_high_speed_detections <- df.alltags.912.clean %>%
  filter(motusTagID %in% high_speed_individuals$motusTagID)  # Keep all detections for these individuals

# For each individual, run, and receiver, select the first detection (`ts`)
all_high_speed_detections <- all_high_speed_detections %>%
  #filter(is.infinite(gspeedMS) | !is.infinite(gspeedMS)) %>%
  group_by(runID) %>%
  arrange(ts) %>%  # Ensure we get the first timestamp for each receiver
  filter(ts == min(ts) | is.infinite(gspeedMS) | gspeedMS > 100) %>%  # Select the first row per group OR any rows where the gspeed appears infinite because the detections occur at the exact same time (but the hits are not the first ones in each run) OR hits that resulted in very high flight speeds
  ungroup() |> # Ungroup after selection
  select(runID, hitID, ts, motusTagID, runLen, tagModel, tagBI, speciesEN, recvDeployName, distBtwRecvM, distBtwRecvKM, timeH, timeD, timeS, gspeedMS, gspeedKMH) |>
  arrange(motusTagID, ts)

# checking that I've caught all of the high or infinite speed cases

check <- high_speed_individuals |>
  filter(!hitID %in% all_high_speed_detections$hitID)
# none, perfect

# looking at the df, most high (>100m/s or Inf) flight speeds are for receivers <65km apart
# reserva de la biosfera gets flagged, 3570 apart and time diff of under an hour

# 51896  has some infinite flight speeds

df.51896 <- all_high_speed_detections |>
  filter(motusTagID == 51896) |>
  select(runID, ts, recvDeployName, distBtwRecvM, distBtwRecvKM, timeH, timeS, gspeedMS, gspeedKMH)
# all between APR-Schmoekel and Blue Ridge. they're 65km apart, so it's possible. but there are two closer receivers that the tag wasn't picked up on simultaneously
# both receivers have lots of good looking hits for multiple individuals. so I think they're real, just a bit odd. 
# APR has yagi antennas and so does blue ridge
# the ones in between are also yagis though so I'd expect hits there too

# checking out the hits at Reserva de la Biosfera Ria Celestun that were flagged
df.67510.reserva <- df.alltags.912.clean |>
  filter(motusTagID == 67510 & grepl("Reserva", recvDeployName)) |>
  mutate(reasonRemoved = "flightSpeed")
# a single run of five that occurred only 6 days after the tag was deployed. I will remove it (runID 749628849)

df.alltags.912.clean <- df.alltags.912.clean |>
  filter(!runID %in% df.67510.reserva$runID)

#### 17 - saving data files ####

saveRDS(df.alltags.912.clean, "df.alltags.912.clean.2025Jan20.rds")

# data that I removed
saveRDS(noisyTowerHits, "noisyTowerHits.2025Jan20.rds")
saveRDS(hitsOutOfRange, "hitsOutOfRange.2025Jan20.rds") 
saveRDS(droppedTags, "droppedTag.2025Jan20.rds")
saveRDS(weirdTowers, "weirdRecv.2025Jan20.rds")
saveRDS(moreNoisyHits, "shortRunsNoisyRecv.2025Jan20.rds")
saveRDS(df.67510.reserva, "flightSpeedHitsRemoved.2025Jan20.rds")

#### 18 - compiling removed hits ####
# to make them easier to look at later and compare to the clean data, I will merge together all of the dataframes for the data I've removed. they all have the "reasonRemoved" column so they can be filtered by issue

# don't need to load them if you already have them in your R environment from running script 01

shortRuns <- readRDS("shortRuns912.2025Jan20.rds")
undeployedTags <- readRDS("detsUndepTags.2025Jan20.rds")
dets2021 <- readRDS("2021Detections.2025Jan15.rds")
preDepDets <- readRDS("preDeploymentTagDets.2025Jan15.rds")
noTagData <- readRDS("postDepDets_missingTagInfo.2025Jan15.rds")
dups <- readRDS("duplicateDataRemoved.2025Jan20.rds")
noRecvLatLon <- readRDS("noRecvInfo.2025Jan20.rds")
noisyTowers <- readRDS("noisyTowerHits.2025Jan20.rds")
range <- readRDS("hitsOutOfRange.2025Feb11.rds")
droppedTag <- readRDS("droppedTag.2025Feb11.rds")
weirdRecv <- readRDS("weirdRecv.2025Feb11.rds")
shortRunsNoisyRecv <- readRDS("shortRunsNoisyRecv.2025Feb11.rds")
flightSpeed <- readRDS("flightSpeedHitsRemoved.2025Jan20.rds")

allHitsRemoved <- shortRuns |>
  bind_rows(undeployedTags, dets2021, preDepDets, noTagData, dups, noRecvLatLon, noisyTowers, range, droppedTag, weirdRecv, shortRunsNoisyRecv,flightSpeed)

saveRDS(allHitsRemoved, "allHitsRemoved.2025Feb11.rds")

# and then can upload to OSF - note that if the file already exists there (if it's one of the ones I uploaded there), you may have to specify conflict = "rename" to upload. I just put the upload option here in case you want it. 
local_files <- c("noisyTowerHits.2025Jan20.rds", "hitsOutOfRange.2025Feb11.rds", "droppedTag.2025Feb11.rds",  "weirdRecv.2025Feb11.rds", "shortRunsNoisyRecv.2025Feb11.rds", "flightSpeedHitsRemoved.2025Jan20.rds", "allHitsRemoved.2025Feb11.rds")

# Upload the file to the project (path argument specifies folder in the project, leaving blank uploads it to the main directory)
osf_upload(project, local_files)

# 19 - summary table for detections removed ###

summary_df <- allHitsRemoved %>%
  group_by(reasonRemoved) %>%
  summarize(nHitsRemoved = n()) %>%
  arrange(factor(reasonRemoved, levels = unique(allHitsRemoved$reasonRemoved)))

df.alltags.912 <- readRDS("df.alltags.912.2025Jan8.rds")

total_number_detections <- as.numeric(nrow(df.alltags.912))

# Add a total row at the top
total_row <- data.frame(
  reasonRemoved = "Total number of hits",  # Row label for total
  nHitsRemoved = 0,  # No hits removed for the total row
  total_detections = total_number_detections  # Start with total detections
)

# Add the total row at the top
summary_df <- bind_rows(total_row, summary_df)

# Initialize `remaining_detections` as a column of NAs
summary_df$remaining_detections <- NA_real_

# Use `lapply` to apply the calculation row by row
for (i in 1:nrow(summary_df)) {
  if (i == 1) {
    # For the first row, set remaining_detections to total_number_detections
    summary_df$remaining_detections[i] <- total_number_detections
  } else {
    # For subsequent rows, subtract nHitsRemoved from previous remaining_detections
    summary_df$remaining_detections[i] <- summary_df$remaining_detections[i - 1] - summary_df$nHitsRemoved[i]
  }
}

summary_df <- summary_df |>
  select(-total_detections)

print(summary_df)
# this table shows the total # unfiltered detections in the top row, the number removed with each filtering step, and the remaining detections after each step. I've also provided an Excel version of this table on GitHub

saveRDS(summary_df, "tableOfHitsRemoved.2025Feb11.rds")





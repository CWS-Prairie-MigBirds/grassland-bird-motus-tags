#### GBMP Script 03 - Summary statistics, departure dates, migration pace and maps ####
# GBMP = Grassland Bird Monitoring Program
# ECCC contract 3000792304
# adapted from the Motus R book by Sarah Endenburg
# This script contributes to deliverable 2 and constitutes prep work for the R Markdown report

#### 1 - Prepare working environment ####

library(motus)
library(tidyverse)
library(lubridate)
library(ggforce)
library(stringr)
library(rnaturalearth)
library(geosphere)
library(ggrepel)
library(sf)

# Motus data is stored in UTC, set system time to match
Sys.setenv(tz = "UTC")

# turning off scientific notation
options(scipen = 9999)

#### 2 - load data ####
df.alltags.912.clean <- readRDS("df.alltags.912.clean.2025Feb11.rds")

# note that you can download the files from OSF using R - but it just downloads it into your working directory for you, whereas you could just go to the OSF page and download all the relevant files without going through the steps here
# if you want to do it through R, I included code in scripts 1&2 for how to do that. It takes more space in the scripts and I didn't find it saved time or effort so I didn't include that code in scripts 3&4

#### 3 - number of detections by year and species ####

# re-doing the number of days detected now that detections have been filtered

df.alltags.912.clean <- df.alltags.912.clean |>
  group_by(motusTagID) |>
  mutate(numDaysDetected = as.numeric(difftime(max(ts), min(ts), units = "days")))

hist(df.alltags.912.clean$numDaysDetected)
# 0 to 500 days

# summarizing detections by species and year
detections_by_Species_Year <- df.alltags.912.clean |>
  ungroup() |>
  group_by(tagDepYear, speciesEN) |>
  summarise(nTagsDetected = n_distinct(motusTagID),
            nHits = n(),
            avgHitsPerTag = nHits/nTagsDetected,
            avgDaysDetected = mean(numDaysDetected),
            nRecv = n_distinct(recvDeployName),
            .groups = 'drop') |>
  arrange(speciesEN, tagDepYear)

# adding in number of tags deployed per year and species, to compare with number detected
tags912 <- read_csv("MotusTagsTemplate-912.csv")

tags912 <- tags912 |>
  mutate(speciesEN = case_when(
    speciesName == "Centronyx bairdii" ~ "Baird's Sparrow",
    speciesName == "Calcarius ornatus" ~ "Chestnut-collared Longspur",
    speciesName == "Eremophila alpestris" ~ "Horned Lark",
    speciesName == "Anthus spragueii" ~ "Sprague's Pipit",
    speciesName == "Ammodramus savannarum" ~ "Grasshopper Sparrow",
    speciesName == "Rhynchophanes mccownii" ~ "Thick-billed Longspur",
    TRUE ~ "NA"
  ))

tagDeps <- tags912 |>
  filter(is.na(decimalLatitude) == FALSE) |>
  rename(tagDepYear = utcYearStart) |>
  group_by(tagDepYear, speciesEN) |>
  summarise(nTagsDeployed = n())

detections_by_Species_Year <- detections_by_Species_Year |>
  left_join(tagDeps, by = c("tagDepYear", "speciesEN")) |>
  mutate(percentTagsDetected = (nTagsDetected/nTagsDeployed)*100) |>
  select(tagDepYear, speciesEN, nTagsDeployed, nTagsDetected, percentTagsDetected, nHits, avgHitsPerTag, avgDaysDetected, nRecv) |>
  mutate(across(where(is.numeric), ~ round(., 1))) |>
  arrange(speciesEN, tagDepYear) 

#### 4 - number of detections/detection days as a function of deployment distance from tower ####

## 4a) using the closest receiver on the exact date of tag deployment ####
# need to find closest tower to each tag deployment location
# note that this can give slightly misleading results because some tags were deployed a few weeks before a tower was deployed close by. 4b addresses that issue

depTags912 <- tags912 |>
  filter(is.na(utcYearStart) == FALSE) |>
  select(tagID, mfgID, period, utcYearStart, utcMonthStart, utcDayStart, utcHourStart, utcMinuteStart, decimalLatitude, decimalLongitude) |>
  mutate(tagDeployDate = make_datetime(
    utcYearStart,
    utcMonthStart,
    utcDayStart
  )) |>
  rename(motusTagID = tagID,
         tagBI = period,
         tagDepLat = decimalLatitude,
         tagDepLon = decimalLongitude) |>
  filter(is.na(tagDepLat) == FALSE) # there are two tags with no deployment lat/long, need those from nancy still

# I need the closest receiver to each tag site, but not just the closest receiver with a detection because for many tags they were only detected on very distant receivers. instead, I need the closest deployed receiver, regardless of whether there was a detection there

# loading a file that should contain all receivers for all time, downloaded from the Motus site

allRecv <- read_csv("receiver-deployments.csv")

allRecv$dtStart <- as.POSIXct(allRecv$dtStart)
allRecv$dtEnd <- as.POSIXct(allRecv$dtEnd)

# Step 1: Filter active receivers at the time of each tag deployment

allRecv <- allRecv |>
  rename(recvDeployName = deploymentName,
         recvDeployLat = latitude,
         recvDeployLon = longitude)

# A function to check if a receiver was active at the time of tag deployment
get_active_receivers <- function(tagDeployDate) {
  active_receivers <- allRecv %>% # have to have recvDeployName in there, it's not in allRecv
    filter(dtStart <= tagDeployDate & (dtEnd >= tagDeployDate | is.na(dtEnd))) %>% # recv that are still active have an NA for their end date!
    select(recvDeployName, recvDeployLat, recvDeployLon, dtStart, dtEnd)
  return(active_receivers)
}

# Step 2: Calculate the distance between tag deployments and active receivers
# A function to calculate the closest receiver using the Vincenty formula

calculate_closest_receiver <- function(tag_row) {
  # Get the active receivers at the time of tag deployment
  active_receivers <- get_active_receivers(tag_row$tagDeployDate)
  
  # Calculate distances to all active receivers
  distances <- distVincentyEllipsoid(cbind(active_receivers$recvDeployLon, active_receivers$recvDeployLat),
                                  c(tag_row$tagDepLon, tag_row$tagDepLat))
  
  # Find the closest receiver
  closest_idx <- which.min(distances)
  
  # Return the closest receiver information
  return(data.frame(tag_id = tag_row$motusTagID, 
                    closest_receiver = active_receivers$recvDeployName[closest_idx],
                    distance_km = distances[closest_idx] / 1000))  # convert meters to km
}

# Step 3: Apply the calculation for each tag deployment

closest_receivers <- do.call(rbind, lapply(1:nrow(depTags912), function(i) {
  tag_row <- depTags912[i, ] # get the current tag's row
  calculate_closest_receiver(tag_row) # pass the tag row to the function
}))

# Combine the results
closest_receivers

closest_receivers <- closest_receivers |>
  rename(motusTagID = tag_id,
         distToClosestRecv = distance_km)

detections_deploy_dist_from_recv <- df.alltags.912.clean |>
  ungroup() |>
  group_by(speciesEN, tagModel, motusTagID) |>
  summarise(nHits = n(),
            numDaysDetected = mean(numDaysDetected),
            nRecv = n_distinct(recvDeployName),
            .groups = 'drop')

detections_deploy_dist_from_recv <- detections_deploy_dist_from_recv |>
  left_join(closest_receivers, by = "motusTagID")

# plotting the results
ggplot(detections_deploy_dist_from_recv) +
  geom_point(aes(x = distToClosestRecv, y = nRecv, fill = speciesEN), size = 4, shape = 21, alpha = 0.5, position = position_jitter(width = 0.2, height = 0.2)) +
  theme_classic() +
  labs(x = "Distance from tag deployment site to closest active receiver (km)", y = "Number of detections") +
  facet_wrap(~speciesEN) + 
  scale_fill_viridis_d() + 
  theme(axis.title.x = element_text(size = 11), 
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)) +
  guides(fill = "none")

## 4b) using distance to any receiver that was active during breeding season ####
# we'll say that the breeding season ends at the end of July - I didn't have any reasonable looking departures until at least August

# Step 1: Filter active receivers at the time of each tag deployment
# A function to check if a receiver was active at the time of tag deployment

get_active_receivers <- function(tagDeployDate) {
  
  tagDepYear <- format(tagDeployDate, "%Y")
  tagDepMonth <- format(tagDeployDate, "%m")
  
  end_of_july <- as.Date(paste(tagDepYear, "07", "31", sep = "-"))
  
  active_receivers <- allRecv %>% 
    filter(dtStart <= end_of_july & (dtEnd >= end_of_july | is.na(dtEnd))) %>% # recv that are still active have an NA for their end date!
    select(recvDeployName, recvDeployLat, recvDeployLon, dtStart, dtEnd)
  return(active_receivers)
}

# Step 2: Calculate the distance between tag deployments and active receivers
# A function to calculate the closest receiver using the Haversine formula

calculate_closest_receiver <- function(tag_row) {
  
  # Get the active receivers at the time of tag deployment
  active_receivers <- get_active_receivers(tag_row$tagDeployDate)
  
  # Calculate distances to all active receivers
  distances <- distVincentyEllipsoid(cbind(active_receivers$recvDeployLon, active_receivers$recvDeployLat),
                                     c(tag_row$tagDepLon, tag_row$tagDepLat))
  
  # Find the closest receiver
  closest_idx <- which.min(distances)
  
  # Return the closest receiver information
  return(data.frame(tag_id = tag_row$motusTagID, 
                    closest_receiver = active_receivers$recvDeployName[closest_idx],
                    distance_km = distances[closest_idx] / 1000))  # convert meters to km
}

closest_receivers <- do.call(rbind, lapply(1:nrow(depTags912), function(i) {
  tag_row <- depTags912[i, ] # get the current tag's row
  calculate_closest_receiver(tag_row) # pass the tag row to the function
}))

# Combine the results
closest_receivers

closest_receivers <- closest_receivers |>
  rename(motusTagID = tag_id,
         distToClosestRecv = distance_km)

detections_deploy_dist_from_recv <- df.alltags.912.clean |>
  ungroup() |>
  group_by(speciesEN, tagModel, motusTagID) |>
  summarise(nHits = n(),
            numDaysDetected = mean(numDaysDetected),
            nRecv = n_distinct(recvDeployName),
            .groups = 'drop')

detections_deploy_dist_from_recv <- detections_deploy_dist_from_recv |>
  left_join(closest_receivers, by = "motusTagID")

# plotting the results. I used this plot (or a version of it) for the R markdown
ggplot(detections_deploy_dist_from_recv) +
  geom_point(aes(x = distToClosestRecv, y = numDaysDetected, fill = speciesEN), size = 4, shape = 21, alpha = 0.5, position = position_jitter(width = 0.2, height = 0.2)) +
  theme_classic() +
  labs(x = "Distance from tag deployment site to closest active receiver (km)", y = "Number of detections") +
  facet_wrap(~speciesEN) + 
  scale_fill_viridis_d() + 
  theme(axis.title.x = element_text(size = 11), 
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  guides(fill = "none")
# can use numDaysDetected or nRecv as the y variable

#### 5 - number of detections/detection days as a function of tag type ####

# this is what tags were detected
detections_by_tagModel <- df.alltags.912.clean |>
  ungroup() |>
  group_by(tagModel) |>
  summarise(nTagsDetected = n_distinct(motusTagID),
            nHits = n(),
            avgHitsPerTag = nHits/nTagsDetected,
            avgDaysDetected = mean(numDaysDetected),
            nRecv = n_distinct(recvDeployName),
            .groups = 'drop')

# to compare with what was deployed:
# there are 4 types of tags in the data but only 3 in the banding sheet - HybridTag and LifeTag just say CTT solar. it shows on the motus site but for some reason not in the MotusTagTemplate912 that can be downloaded...

# I downloaded a different file that has that info
tags <- read_csv("tags.csv")

tagModels <- tags |>
  select(tagID, model) |>
  rename(motusTagID = tagID,
         tagModel = model) 

table(tagModels$tagModel)
# there are 5 types, two are CTT (hybrid and life), three are lotek. 

depTagModels <- tagModels |>
  group_by(tagModel) |>
  summarise(nTagsDeployed = n())

detections_by_tagModel <- depTagModels |>
  left_join(detections_by_tagModel, by = "tagModel") |>
  mutate(manufacturer = case_when(
    tagModel %in% c("LifeTag", "HybridTag") ~ "CTT",
    TRUE ~ "Lotek"
  )) |>
  mutate(percentTagsDetected = (nTagsDetected/nTagsDeployed)*100) |>
  select(manufacturer, tagModel, nTagsDeployed, nTagsDetected, percentTagsDetected, nHits, avgHitsPerTag, avgDaysDetected, nRecv)

# if we ignore confounding factors like species, deployment location, etc, a higher percentage of lotek tags were detected. that could be because CTT tag digital IDs apparently degrade more easily than Lotek signals, and might not get interpreted as a tag signal. there also are probably more false detections for Lotek tags because of the way the coding works. and 19/51 of the CTT tags were solar, so not emitting signals if they haven't been exposed to sun (or maybe the signal is weaker if battery is low)

#### 6 - departure dates by species and year ####

# there are both fall and spring migration detections in these data, so those will need to be separate departure dates

# fall departure dates ####

df.alltags.912.clean <- df.alltags.912.clean |> 
  ungroup() |>
  mutate(distFromTagSiteM = distVincentyEllipsoid(df.alltags.912.clean[,c("recvDeployLon", "recvDeployLat")],
                                                  df.alltags.912.clean[,c("tagDepLon", "tagDepLat")]),
         distFromTagSiteKM = distFromTagSiteM/1000)

# fall migration would be same calendar year as the deployment, and obviously starting from where the bird was tagged (or nearby)
# actually, some birds (eg. 75838-HOLA) had fall departures the following year too - I will look at those later on

# Ellison et al. 2017 found that CCLO tagged with geolocators moved up to ~120 km post-breeding but pre-migration. I think the error on geolocators can be up to ~100km...but that seems to match roughly with the detections that I see anyways 
# I have not seen any similar data for the other species so I'm going to look at 100km, 150km and 200km for all species - I also tried 50km before (what I used for BANS) and got almost the exact same results for departure dates

fallDepDates100km <- df.alltags.912.clean |>
  filter(distFromTagSiteKM <= 100 & tagDepYear == year & numDaysDetected > 1) |>
  group_by(motusTagID, speciesEN, tagDepYear, tagDepLat, tagDepLon, numDaysDetected) |>
  summarize(depDate100km = max(ts)) |>
  arrange(speciesEN, tagDepYear)

fallDepDates150km <- df.alltags.912.clean |>
  filter(distFromTagSiteKM <= 150 & tagDepYear == year & numDaysDetected > 1) |>
  group_by(motusTagID, speciesEN, tagDepYear, tagDepLat, tagDepLon, numDaysDetected) |>
  summarize(depDate150km = max(ts)) |>
  arrange(speciesEN, tagDepYear)

fallDepDates200km <- df.alltags.912.clean |>
  filter(distFromTagSiteKM <= 200 & tagDepYear == year & numDaysDetected > 1) |>
  group_by(motusTagID, speciesEN, tagDepYear, tagDepLat, tagDepLon, numDaysDetected) |>
  summarize(depDate200km = max(ts)) |>
  arrange(speciesEN, tagDepYear)

# adding departure dates to the main dataframe
deps100 <- fallDepDates100km |>
  ungroup()|>
  select(motusTagID, depDate100km)

deps150 <- fallDepDates150km |>
  ungroup()|>
  select(motusTagID, depDate150km)

deps200 <- fallDepDates200km |>
  ungroup() |>
  select(motusTagID, depDate200km, speciesEN)

df.alltags.912.clean <- df.alltags.912.clean |>
  left_join(deps100, by = "motusTagID") 

df.alltags.912.clean <- df.alltags.912.clean |>
  left_join(deps150, by = "motusTagID")

df.alltags.912.clean <- df.alltags.912.clean |>
  left_join(deps200, by = "motusTagID") |>
  select(-speciesEN.y) |>
  rename(speciesEN = speciesEN.x)

# and combining all dep dates together 
allDepDates <- deps200 |>
  left_join(deps150, by = "motusTagID") |>
  left_join(deps100, by = "motusTagID") |>
  select(motusTagID, speciesEN, depDate100km, depDate150km, depDate200km) |>
  arrange(speciesEN, depDate100km)

# CCLO - Ellison et al. 2017 earliest departure (from SK) was 30 Sept, Birds of the World says fall migration starts around mid-Sept to early-Oct
# according to Birds of the World, SPPI fall migration starts in September in the prairies
# according to Birds of the World, Baird's sparrow has left the breeding grounds by mid-late Sept...no indication when they start leaving though
# COSEWIC assessment for BAIS suggests they start leaving in Sept-Oct (Green et al. 2002 ref, which is Birds of North America and no longer available...)

# comparing the departure dates that I pulled above with my detection plots to verify 

df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(bandingCode = case_when(
    speciesEN == "Horned Lark" ~ "HOLA",
    speciesEN == "Chestnut-collared Longspur" ~ "CCLO",
    speciesEN == "Sprague's Pipit" ~ "SPPI",
    speciesEN == "Baird's Sparrow" ~ "BAIS",
    speciesEN == "Thick-billed Longspur" ~ "TBLO",
    speciesEN == "Grasshopper Sparrow" ~ "GRSP",
    TRUE ~ "NA"
  ))

allDepDates <- allDepDates |>
  mutate(bandingCode = case_when(
    speciesEN == "Horned Lark" ~ "HOLA",
    speciesEN == "Chestnut-collared Longspur" ~ "CCLO",
    speciesEN == "Sprague's Pipit" ~ "SPPI",
    speciesEN == "Baird's Sparrow" ~ "BAIS",
    speciesEN == "Thick-billed Longspur" ~ "TBLO",
    speciesEN == "Grasshopper Sparrow" ~ "GRSP",
    TRUE ~ "NA"
  ))

for(i in 1:nrow(allDepDates)){
  print(ggplot(data = filter(df.alltags.912.clean, 
                             motusTagID == allDepDates$motusTagID[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(allDepDates$motusTagID[i], allDepDates$bandingCode[i], sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# sometimes the bird was detected later on in migration (indicating the tag didn't drop or completely die), but with such a large time gap between the breeding detections to migration that it's unlikely to have been a direct movement. rather, there must have been long periods without detections

# 67391-BAIS has an august departure but with good tracks down into Montana. seems legit.
# 61577-BAIS probably not a real departure date (July) - only detected at Nashlyn and Snake Butte which are only about 100km from one another. snake butte is about 110km from the tag site which is why that isn't the departure date itself

# going to look at departures along with the next detection for each individual that has migratory detections
df_sorted <- df.alltags.912.clean %>%
  filter(!is.na(depDate100km) & ts >= depDate100km) |>
  arrange(motusTagID, recvDeployName, ts)

# pulling just the first post-departure detection
df_next_detection <- df_sorted |>
  filter(ts > depDate100km) |>
  group_by(motusTagID) |>
  arrange(ts) |>
  slice_head(n = 1)

# and just the departures themselves
df_deps <- df_sorted |>
  filter(ts == depDate100km)

# binding together
df_deps_next_det <- df_deps |>
  bind_rows(df_next_detection) |>
  select(motusTagID, speciesEN, ts, depDate100km, recvDeployName, recvDeployLat, distBtwRecvKM, timeD, numDaysDetected, tagDepYear) |>
  mutate(depMonth100 = month(depDate100km))|>
  arrange(speciesEN, motusTagID, depDate100km, ts)
# not all birds have a next detection

# from looking at the data in df_deps_next_det and comparing with the tracks on the figures generated above, I came up with the following classification scheme to rank how confident I am that the departure dates I calculated are true migratory departures

# birds that have detections at over 100km from the tag site (migratory detections, since I'm using the last detection within 100km as the departure date)
migBirds <- df.alltags.912.clean |>
  filter(distFromTagSiteKM >100)

newDepDates <- df_deps_next_det |>
  mutate(typeOfDep = case_when(
    depMonth100 %in% c(8, 9, 10) & motusTagID %in% migBirds$motusTagID & timeD < 16 ~ "migDep", # departure month is reasonable and bird was detected not too long after departing
    depMonth100 %in% c(6, 7) & motusTagID %in% migBirds$motusTagID & !is.na(timeD) ~ "regionalMovement", # departure month is super early but they were detected again later on so the tag wasn't lost
    depMonth100 %in% c(6, 7) & !motusTagID %in% migBirds$motusTagID & is.na(timeD) ~ "possibleTagLoss", # departure month is super early and were not detected again
    depMonth100 %in% c(8, 9, 10) & motusTagID %in% migBirds$motusTagID & timeD > 16 ~ "migDepOrRegional", # departure month is realistic but the long gap between detections makes me less confident that the date is correct
    depMonth100 %in% c(8, 9, 10) & !motusTagID %in% migBirds$motusTagID ~ "possibleMigDep", # departure month is reasonable but were not detected again so can't be sure if they were migrating, tag dropped or regional movement
    TRUE ~ "unknown" # this would be everything else, aka anything with an early departure and that wasn't detected again - tag loss?
  )) |>
  group_by(motusTagID) |>
  filter(ts == max(ts)) |>
  select(-ts) |> 
  rename(nextRecv = recvDeployName,
         nextRecvLat = recvDeployLat, 
         distToNextRecvKM = distBtwRecvKM,
         timeToNextRecvDays = timeD) |>
  mutate(bandingCode = case_when(
    speciesEN == "Horned Lark" ~ "HOLA",
    speciesEN == "Chestnut-collared Longspur" ~ "CCLO",
    speciesEN == "Sprague's Pipit" ~ "SPPI",
    speciesEN == "Baird's Sparrow" ~ "BAIS",
    speciesEN == "Thick-billed Longspur" ~ "TBLO",
    speciesEN == "Grasshopper Sparrow" ~ "GRSP",
    TRUE ~ "NA"
  ))

for(i in 1:nrow(newDepDates)){
  print(ggplot(data = filter(df.alltags.912.clean, 
                             motusTagID == newDepDates$motusTagID[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(newDepDates$motusTagID[i], newDepDates$bandingCode[i], newDepDates$typeOfDep[i], sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

regional <- newDepDates |>
  filter(typeOfDep %in% c("regionalMovement", "migDepOrRegional"))

for(i in 1:nrow(regional)){
  print(ggplot(data = filter(df.alltags.912.clean, 
                             motusTagID == regional$motusTagID[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(regional$motusTagID[i], regional$bandingCode[i], regional$typeOfDep[i], sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# 75839-HOLA does look like it was dropped in mid-late august, hence the detections at Govenlock until December
# 75838 at escape coulee - looks legit, stops in october and is apparently back in late march. sig strength is lower in march but probably cause tag is dying. good run lengths and still variable signal strength and was gone in the winter. was at sage creek in July and then october, not in the spring but sage creek was only up in 2023. looks good though.
# same with 75836 at nashlyn, only at battlecreek in october. but all looks good.
# and 67390-CCLO at Ellice Archie

# fall departures the following year
fallDepDates100km_nextYear <- df.alltags.912.clean |>
  filter(distFromTagSiteKM <= 100 & year == as.numeric(tagDepYear + 1) & numDaysDetected > 1) |>
  group_by(motusTagID, speciesEN, tagDepYear, tagDepLat, tagDepLon, numDaysDetected) |>
  summarize(depDate100km = max(ts)) |>
  arrange(speciesEN, tagDepYear)

deps100_nextYear <- fallDepDates100km_nextYear |>
  ungroup()|>
  select(motusTagID, depDate100km) |>
  rename(depDate100km_nextYear = depDate100km)

allDepDates <- allDepDates |>
  left_join(deps100_nextYear, by = "motusTagID")

df.alltags.912.clean <- df.alltags.912.clean |>
  left_join(deps100_nextYear, by = "motusTagID")

# going to look at departures along with the next detection for each individual that has migratory detections
df_sorted <- df.alltags.912.clean %>%
  filter(!is.na(depDate100km_nextYear) & ts >= depDate100km_nextYear) |>
  arrange(motusTagID, recvDeployName, ts)

# pulling just the first post-departure detection
df_next_detection <- df_sorted |>
  filter(ts > depDate100km_nextYear) |>
  group_by(motusTagID) |>
  arrange(ts) |>
  slice_head(n = 1)
# ok none of them have any further detections

# and just the departures themselves
df_deps <- df_sorted |>
  filter(ts == depDate100km_nextYear)

# binding together
df_deps_next_det <- df_deps |>
  bind_rows(df_next_detection) |>
  select(motusTagID, speciesEN, ts, depDate100km_nextYear, recvDeployName, recvDeployLat, distBtwRecvKM, timeD, numDaysDetected, tagDepYear) |>
  mutate(depMonth100_nextYear = month(depDate100km_nextYear))|>
  arrange(speciesEN, motusTagID, depDate100km_nextYear, ts)
# none of the birds have a next detection
# only 75836 and 75838 HOLA are likely to be actual departure dates, the rest are all in July or earlier

# plotting fall departure dates (same year as tag deployments) ####

# creating a month-day column just for plotting
newDepDates$Month_Day100 <- format(as.Date(newDepDates$depDate100km), "%m-%d")
newDepDates$Month_Day100 <- as.Date(newDepDates$Month_Day100, "%m-%d")

ggplot() +
  geom_boxplot(data = filter(newDepDates, typeOfDep %in% c("migDep", 
                                                           "possibleMigDep") & speciesEN %in% c("Sprague's Pipit", "Chestnut-collared Longspur")), aes(x = as.factor(tagDepYear), y = Month_Day100, fill = as.factor(speciesEN)), colour = "#505050", outlier.shape = NA, alpha = 0.3) +
  geom_point(data = filter(newDepDates, typeOfDep %in% c("migDep", 
                                                         "possibleMigDep")), aes(x = as.factor(tagDepYear), y = Month_Day100, fill = as.factor(speciesEN), shape = as.factor(typeOfDep)), size = 3, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.1), colour = "black") +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b %d")+
  coord_flip()+
  facet_wrap(~speciesEN) +
  theme_classic() +
  scale_colour_manual(values = c("Sprague's Pipit" = "#440154", "Horned Lark" = "#35b779", "Chestnut-collared Longspur" = "#31688e", "Baird's Sparrow" = "#fde725", "Thick-billed Longspur" = "#21918c"))+
  scale_fill_manual(values = c("Sprague's Pipit" = "#440154", "Horned Lark" = "#35b779", "Chestnut-collared Longspur" = "#31688e", "Baird's Sparrow" = "#fde725", "Thick-billed Longspur" = "#21918c")) +
  scale_shape_manual(values = c("migDep" = 21, "possibleMigDep" = 24), # Default shapes for each type
                     labels = c("Migratory departure", "Possible migratory departure")) +  # Custom labels for legend
  labs(y = "Departure Date", x = "Year") +
  guides(colour = "none", fill = "none", shape = guide_legend(title = "Type of Departure")) +
  theme(legend.position = "bottom",   # Place legend at the bottom
        legend.direction = "horizontal",  # Arrange legend items horizontally
        legend.justification = c("center", "top"),  # Center the legend horizontally
        legend.box = "horizontal",  # Ensure the items are in a horizontal box
        legend.box.margin = margin(t = 10)) +  # Optional: space above the legend +
  theme(text = element_text(size = 16))

# spring departure dates ####

# spring migration would be the following calendar year
# I'm going to find the southernmost point and use the last detection within ~100km that occurs the year after tag deployment
# since I can't currently determine wintering areas, this would be the earliest northward movement essentially...

springDepDates <- df.alltags.912.clean |>
  filter((year == tagDepYear + 1)) |>
  group_by(motusTagID, speciesEN) |>
  filter(recvDeployLat == min(recvDeployLat)) |>
  filter(ts == max(ts)) |>
  select(ts, motusTagID, bandingCode, recvDeployLat, recvDeployName) |>
  filter(recvDeployLat < 49) |>
  arrange(speciesEN, recvDeployLat)
# 25 birds with detections the following year
# some of these are right back on the breeding grounds though
# changed it to recvLat < 49 --> 17 birds. otherwise I was getting some breeding ground detections in there

# checking that they make sense
for(i in 1:nrow(winterLocs)){
  print(ggplot(data = filter(df.alltags.912.clean, 
                             motusTagID == winterLocs$motusTagID[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(winterLocs$motusTagID[i], winterLocs$bandingCode[i], sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# plotting

springDepDates$Month_Day <- format(as.Date(springDepDates$ts), "%m-%d") # creates a string
springDepDates$Month_Day <- as.Date(springDepDates$Month_Day, "%m-%d")

ggplot(springDepDates) +
  geom_point(aes(x = Month_Day, y = as.factor(year(ts)), fill = as.factor(speciesEN)), colour = "black", shape = 21, size = 4) +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 20)) +
  labs(x = "Date of first spring detection", y = "Year", fill = "Species")

#### 7 - Spring arrival dates ####
# this is the first detection within 100km of the tagging location in the year following tag deployment

springArrivals <- df.alltags.912.clean |>
  filter(year == (tagDepYear + 1) & distFromTagSiteKM <= 100) |>
  group_by(motusTagID) |>
  arrange(ts) |>
  slice_head(n = 1) |>
  select(speciesEN, bandingCode, ts, recvDeployName, recvDeployLat, distFromTagSiteKM, tagDepYear) |>
  arrange(speciesEN, recvDeployLat) 

springArrivals$Month_Day <- format(as.Date(springArrivals$ts), "%m-%d") # creates a string
springArrivals$Month_Day <- as.Date(springArrivals$Month_Day, "%m-%d") # turns it into a date for plotting

# 75838-HOLA arrived in march? so did 75836-HOLA. the rest are in April except for one SPPI in may
# yeah seems to be the case. interesting that both HOLA with arrival dates arrived in March

springArrivals |>
  group_by(speciesEN, tagDepYear) |>
  summarize(nTags = n_distinct(motusTagID))
# none for 2024 tags yet because they're not back yet

ggplot(springArrivals) +
  geom_point(aes(x = Month_Day, y = as.factor(year(ts)), fill = as.factor(speciesEN)), colour = "black", shape = 21, size = 4) +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 20)) +
  labs(x = "Spring breeding ground arrival date", y = "Year", fill = "Species")
  
#### 8 - Migration pace by species and year ####

## fall migration pace ####

# I'm going to calculate an average migration pace by calculating the distance between departure and the southernmost detection and divide by time elapsed

# I want the departure dates (for birds with valid departures) and furthest south detection in the same calendar year as tag deployment for each bird
# also taking the detection with the highest signal strength at each receiver (on each day) for each individual - represents the point at which the bird was probably the closest to that receiver
# also only want birds with detections at least 200km away - this is an arbitrary number, but it's double what I used for departure date cutoffs, so just selects for birds that had actual migratory detections
migBirds2 <- df.alltags.912.clean |>
  filter(distFromTagSiteKM > 200)

migPaceDets <- df.alltags.912.clean |>
  mutate(ordinalDate = yday(ts)) |>
  filter(!is.na(depDate100km) & ts >= depDate100km & year == tagDepYear & motusTagID %in% migBirds2$motusTagID) |>
  group_by(motusTagID, recvDeployName) |>
  filter(sig == max(sig)) # some birds have multiple detections with the same signal strength at an individual receiver

# now selecting just the departure date and southernmost point (for birds that have them)
migPaceDets2 <- migPaceDets |>
  group_by(motusTagID) |>
  arrange(ts) |>
  filter(ts == depDate100km | recvDeployLat == min(recvDeployLat)) |>
  select(motusTagID, speciesEN, tagDepYear, ts, recvDeployName, recvDeployLat, recvDeployLon, sig)

migPaceDets2 <- migPaceDets2 |>
  ungroup() |>
  arrange(motusTagID, ts) 

migPace1 <- migPaceDets2 |>
  mutate(timeDiff = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName )),
                               difftime(ts, lag(ts), units = "days"),
                               NA),
            totalDistM = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName)),
                                                 c(NA, distVincentyEllipsoid(migPaceDets2[,c("recvDeployLon", "recvDeployLat")])),
                                                 NA_real_),
            totalDistKM = totalDistM/1000,
            pace = totalDistKM/as.numeric(timeDiff)) |>
  filter(is.na(pace) == FALSE)

# 85822-CCLO has the highest pace by far - most are 200km/day or under, this one is over 1000km/day because it traveled 500km in a half day

# fall migration pace plot/table - note that I've slightly improved on these for the R markdown report
migPace1 |>
  group_by(speciesEN, tagDepYear) |>
  summarize(nTags = n_distinct(motusTagID))

ggplot(filter(migPace1, pace < 500)) +
  geom_point(aes(x = pace, y = as.factor(speciesEN), fill = as.factor(speciesEN), shape = as.factor(year(ts))), colour = "black", size = 4, position = position_jitter(height = 0.1, width = 0.1), alpha = 0.7) +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 14)) +
  labs(x = "Migration pace (km/day)", y = "Species", shape = "Year") +
  scale_colour_viridis_d() +
  scale_shape_manual(values = c(21, 22, 24)) +
  guides(fill = "none")

paceTable <- migPace1 |>
  filter(pace < 500) |>
  group_by(speciesEN) |>
  summarize(nTags = n_distinct(motusTagID),
            maxPace = max(pace),
            minPace = min(pace),
            avgPace = mean(pace))

## spring migration pace ####
# calculating average pace between arrival on breeding grounds and the first spring migration detection
springArrivalDates <- springArrivals |>
  select(motusTagID, ts) |>
  rename(springArrivalDate = ts)

springMigPaceDets <- df.alltags.912.clean |>
  filter(year == (tagDepYear + 1))

colnames(springMigPaceDets)
   
springMigPaceDets <- springMigPaceDets |>
  left_join(springArrivalDates, by = "motusTagID") |>
  filter(is.na(springArrivalDate) == FALSE) |>
  group_by(motusTagID, recvDeployName) |>
  filter(sig == max(sig) | ts == springArrivalDate) |> # selecting max signal strength at each receiver, unless it's the arrival date which is just the first detection back on the breeding grounds
  select(motusTagID, speciesEN, tagDepYear, ts, recvDeployName, recvDeployLat, recvDeployLon, sig, springArrivalDate)

# 88 detections for 15 birds
# some birds ONLY have a spring arrival date because they weren't detected during spring migration at all, so I can't calculate pace for them

springMigPaceDets2 <- springMigPaceDets |>
  group_by(motusTagID) |>
  arrange(ts) |>
  filter(recvDeployLat == min(recvDeployLat) | ts == springArrivalDate) |> # each bird should have two points here, their lowest latitude point and their arrival point back on the breeding grounds
  select(motusTagID, speciesEN, tagDepYear, ts, recvDeployName, recvDeployLat, recvDeployLon, sig, springArrivalDate)

# 75837 somehow has 36, 75836 has three
# the recv names are all the same though, so my code below won't use those to calculate pace
# I think maybe because the bird was only detected on the breeding grounds, so the min recvDeployLat pulls all of those hits

springMigPaceDets2 <- springMigPaceDets2 |>
  ungroup() |>
  arrange(motusTagID, ts)

springMigPace1 <- springMigPaceDets2 |>
  mutate(timeDiff = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName )),
                            difftime(ts, lag(ts), units = "days"),
                            NA),
         totalDistM = if_else((motusTagID == lag(motusTagID) & recvDeployName != lag(recvDeployName)),
                              c(NA, distVincentyEllipsoid(springMigPaceDets2[,c("recvDeployLon", "recvDeployLat")])),
                              NA_real_),
         totalDistKM = totalDistM/1000,
         pace = totalDistKM/as.numeric(timeDiff)) |>
  filter(is.na(pace) == FALSE)

# all reasonable except 67509-SPPI is over 1000 km/day. only traveled ~180km though (Montana into Saskatchewan). the very low pace (0.6 km/day) is for 75406-SPPI that only traveled 7km which is not far enough between receivers - was basically only detected on the breeding grounds and doesn't have a spring "departure" date...I will limit this plot to only birds with spring departure dates

# spring mig pace plot

ggplot(filter(springMigPace1, pace < 500 & motusTagID %in% springArrivals$motusTagID)) +
  geom_point(aes(x = pace, y = as.factor(speciesEN), fill = as.factor(speciesEN), shape = as.factor(year(ts))), colour = "black", size = 4, position = position_jitter(height = 0.1, width = 0.1), alpha = 0.7) +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 14)) +
  labs(x = "Migration pace (km/day)", y = "Species", shape = "Year") +
  scale_colour_viridis_d() +
  scale_shape_manual(values = c(21, 22, 24)) +
  guides(fill = "none")

#### 9 - Migration maps by species ####

# loading winter range maps from eBird

CCLO_range <- st_read("chclon_range_2022/chclon_range_2022.gpkg")
CCLO_winter_range <- CCLO_range |>
  filter(season == "nonbreeding")

SPPI_range <- st_read("sprpip_range_2022/sprpip_range_2022.gpkg")
SPPI_winter_range <- SPPI_range |>
  filter(season == "nonbreeding")

BAIS_range <- st_read("baispa_range_2022/baispa_range_2022.gpkg")
BAIS_winter_range <- BAIS_range |>
  filter(season == "nonbreeding")

HOLA_range <- st_read("horlar_range_2022/horlar_range_2022.gpkg")
HOLA_winter_range <- HOLA_range |>
  filter(season == "nonbreeding")

TBLO_range <- st_read("mcclon_range_2022/mcclon_range_2022.gpkg")
TBLO_winter_range <- TBLO_range |>
  filter(season == "nonbreeding")

GRSP_range <- st_read("graspa_range_2022/graspa_range_2022.gpkg")
GRSP_winter_range <- GRSP_range |>
  filter(season == "nonbreeding")

allWinterRanges <- CCLO_winter_range |>
  bind_rows(SPPI_winter_range, BAIS_winter_range, HOLA_winter_range, TBLO_winter_range, GRSP_winter_range) |>
  rename(speciesEN = common_name)

# loading basemap layers
world <- ne_countries(scale = "medium", returnclass = "sf")

states_CanUS <- ne_states(country = c("canada", "united states of america"), returnclass = "sf")

lakes <- ne_download(scale = "large", type = 'lakes', category = 'physical',
                     returnclass = "sf", destdir = "mapData") # only need this first time downloading

#lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
                 #returnclass = "sf",
                 #destdir = "mapData") # use this if already downloaded shapefiles

# function to connect consecutive receiver detections for each individual
fun.getpath <- function(df) {
  df %>%
    filter(!is.na(recvDeployLat) | !(recvDeployLat == 0)) %>% 
    group_by(motusTagID, mfgID, runID, recvDeployName, ambigID, 
             tagDepLon, tagDepLat, recvDeployLat, recvDeployLon, speciesEN) %>%
    summarize(max.runLen = max(runLen), time = mean(ts), .groups = "drop") %>%
    arrange(motusTagID, time) %>%
    data.frame()
}

# generating the path data
path912 <- fun.getpath(df.alltags.912.clean)

tagDepYears <- tags912 |>
  rename(motusTagID = tagID) |>
  select(motusTagID, utcYearStart) |>
  rename(tagDepYear = utcYearStart) |>
  filter(is.na(tagDepYear) == FALSE)

path912 <- path912 |>
  left_join(tagDepYears, by = "motusTagID")

# Identifying the first receiver for each tag so that I can connect the tag deployment location to the first receiver (fun.getpath just connects receivers)
first_receivers <- path912 %>%
  group_by(motusTagID) %>%
  slice(1) %>%
  ungroup()

# adding these columns just for plotting purposes
path912$tagSite <- "Tag Site"
allWinterRanges$Season <- "Nonbreeding Season Range"
path912$recvWithDets <- "Receiver With Detections"
path912$lineType <- "Apparent Migration Route"
first_receivers$lineType <- "Apparent Migration Route"

# defining the boundaries of the map
xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 10
ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1

migTracks <- ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = allWinterRanges, aes(fill = as.factor(Season)), color = NA, alpha = 0.3) + 
  geom_sf(data = states_CanUS, colour = "black", fill = NA) +
  geom_sf(data = lakes, colour = NA, fill = "white") +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  labs(x = "", y = "") +
  geom_path(data = filter(path912, is.na(speciesEN) == FALSE), 
            aes(x = recvDeployLon, y = recvDeployLat, 
                colour = as.factor(speciesEN), group = motusTagID), show.legend = FALSE) +
  geom_segment(data = first_receivers, aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat, colour = speciesEN), show.legend = FALSE) +
  geom_point(data = filter(path912, is.na(speciesEN) == FALSE ), aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
  geom_point(data = filter(path912, is.na(speciesEN) == FALSE), # plotting tag deployment sites
             aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
  scale_x_continuous(
    breaks = c(-110, -100, -90, -80)
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~speciesEN) +
  # Manual color scales for different layers
  scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
  scale_colour_manual(values = c("Baird's Sparrow" = "#440154", "Chestnut-collared Longspur" = "#443983", "Grasshopper Sparrow" = "#31688e", "Horned Lark" = "#21918c", "Sprague's Pipit" = "#35b779", "Thick-billed Longspur" = "#90d743")) +
  
  # Create a combined, unified legend using guides()
  guides(
    fill = guide_legend(order = 2, title = "Legend", ncol = 1))

# saving the plot
ggsave("All species migration tracks 2025Feb24.tiff", plot = migTracks, width = 10, height = 6, dpi = 600)

#### 10 - Migration maps for individuals ####

## HOLA maps ####
HOLA <- df.alltags.912.clean |>
  filter(speciesEN == "Horned Lark")

HOLAtags <- unique(HOLA$motusTagID) # 15

first_receivers <- path912 %>%
  group_by(motusTagID) %>%
  slice(1) %>%
  ungroup()

activeRecv <- allRecv |>
  filter(dtStart <= "2024-12-01 00:00:00 UTC" & (dtEnd >= "2022-06-01 00:00:00 UTC" | is.na(dtEnd)))
activeRecv$recvType <- "Active Receivers" # for plotting

# note - running this for loop will save a file of each plot in your working directory, you can disable that by # out that line of code (line 947). at the moment, I've created separate files for each species and the maps are saving there, you can modify the code as required
for(i in 1:length(HOLAtags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == HOLAtags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == HOLAtags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
    p <- ggplot(data = world) + 
      geom_sf(colour = "black") +
      geom_sf(data = states_CanUS, colour = "black") +
      geom_sf(data = filter(allWinterRanges, speciesEN == "Horned Lark"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
      geom_sf(data = lakes, colour = NA, fill = "white") +
             coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
             theme_bw() + 
             labs(x = "", y = "") +
      geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
      geom_path(data = path, 
                aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, colour = as.factor(lineType)), size = 1) +
      geom_segment(data = filter(first_receivers, motusTagID == HOLAtags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat, colour = as.factor(lineType)), show.legend = FALSE, size = 1) + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
      geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
      geom_point(data = path, # plotting tag deployment sites
                 aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
      geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
      geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
    labs(title = paste("HOLA", HOLAtags[i], sep = "-")) +
      
      # Manual color scales for different layers
      scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
      scale_colour_manual(values = c("Active Receivers" = "gray45", "Apparent Migration Route" = "#21918c")) +
    
    # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1)) +  # This ensures that shapes are displayed in the legend
      theme(legend.spacing.y = unit(0, "cm"))
    
    print(p)
    
    # saving within a folder called "HOLA maps" in my working directory
    ggsave(filename = file.path("HOLA maps", paste0("HOLA-", HOLAtags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
    
  readline(prompt = 'Press return/enter to continue')}

# HOLA tags to print on the markdown report: 75828, 75830, 85084?
df.75828 <- HOLA |>
  filter(motusTagID == 75828)

df.75830 <- HOLA |>
  filter(motusTagID == 75830)

## SPPI ####
SPPI <- df.alltags.912.clean |>
  filter(speciesEN == "Sprague's Pipit")

SPPItags <- unique(SPPI$motusTagID) # 55

path912 <- path912 |>
  mutate(migrationSeason = case_when(
    as.numeric(year(time)) == tagDepYear ~ "Apparent Fall Migration Route",
    as.numeric(year(time)) == (tagDepYear + 1) ~ "Apparent Spring Migration Route",
    TRUE ~ "NA"
  ))

# as with the HOLA code, this for loop will save a copy of each map in a folder I created in my working directory. edit line 1027 to change that
for(i in 1:length(SPPItags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == SPPItags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == SPPItags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
  p <- ggplot(data = world) + 
    geom_sf(colour = "black") +
    geom_sf(data = filter(allWinterRanges, speciesEN == "Sprague's Pipit"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
    geom_sf(data = lakes, colour = NA, fill = "white") +
    geom_sf(data = states_CanUS, colour = "black", fill = NA) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_bw() + 
    labs(x = "", y = "") +
    geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
    geom_path(data = filter(path, as.numeric(year(time)) == tagDepYear), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#35b779") +
    geom_path(data = filter(path, as.numeric(year(time)) == (tagDepYear + 1)), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#35b779") +
    geom_segment(data = filter(first_receivers, motusTagID == SPPItags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat), show.legend = FALSE, size = 1, colour = "#35b779") + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
    geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
    geom_point(data = path, # plotting tag deployment sites
               aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
    geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
    geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("SPPI", SPPItags[i], sep = "-")) +
    
    # Manual color scales for different layers
    scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
    scale_colour_manual(values = c("Active Receivers" = "gray45")) +
    scale_linetype_manual(values = c("Apparent Fall Migration Route" = "solid", "Apparent Spring Migration Route" = "dotdash")) +
    
  # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1),
      linetype = guide_legend(order = 4, title = NULL, ncol = 1)
      ) +  # This ensures that shapes are displayed in the legend
    theme(legend.spacing.y = unit(0, "cm"))
  
  print(p)
  
  # saving within a folder called "SPPI maps" in my working directory
  ggsave(filename = file.path("SPPI maps", paste0("SPPI-", SPPItags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
  
  readline(prompt = 'Press return/enter to continue')}

# SPPI tags to print on the markdown report: 51879, 51896, 61006, 61020, 61024, 67386, 67507 (connecting the two southern points is probably misleading), 75401, 75402? 75404? 75405? 75406? (same thing with connecting line), 75410, 75412? 75413 (got down to TX), 75781? 75837 (same thing with connecting line), 85757, 85767

df.61024 <- SPPI |>
  filter(motusTagID == 61024)
# Mexico Nov 4

df.75413 <- SPPI |>
  filter(motusTagID == 75413)
# TX Nov 9

df.51896 <- SPPI |>
  filter(motusTagID == 51896)
# NM/TX Oct 27

df.61010 <- SPPI |>
  filter(motusTagID == 61010)
# KS Oct 20

df.61020 <- SPPI |>
  filter(motusTagID == 61020)
# oct 17 in Nebraska 

df.61014 <- SPPI |>
  filter(motusTagID == 61014)
# oct 17 in KS 

## CCLO ####

CCLO <- df.alltags.912.clean |>
  filter(speciesEN == "Chestnut-collared Longspur")

CCLOtags <- unique(CCLO$motusTagID) # 61

# this for loop will save a copy of each map in a folder I created in my working directory. edit line 1120 to change that
for(i in 1:length(CCLOtags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == CCLOtags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == CCLOtags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
  p <- ggplot(data = world) + 
    geom_sf(colour = "black") +
    geom_sf(data = filter(allWinterRanges, speciesEN == "Chestnut-collared Longspur"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
    geom_sf(data = lakes, colour = NA, fill = "white") +
    geom_sf(data = states_CanUS, colour = "black", fill = NA) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_bw() + 
    labs(x = "", y = "") +
    geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
    geom_path(data = filter(path, as.numeric(year(time)) == tagDepYear), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#443983") +
    geom_path(data = filter(path, as.numeric(year(time)) == (tagDepYear + 1)), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#443983") +
    geom_segment(data = filter(first_receivers, motusTagID == CCLOtags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat), show.legend = FALSE, size = 1, colour = "#443983") + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
    geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
    geom_point(data = path, # plotting tag deployment sites
               aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
    geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
    geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("CCLO", CCLOtags[i], sep = "-")) +
    
    # Manual color scales for different layers
    scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
    scale_colour_manual(values = c("Active Receivers" = "gray45")) +
    scale_linetype_manual(values = c("Apparent Fall Migration Route" = "solid", "Apparent Spring Migration Route" = "dotdash")) +
    
    # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1),
      linetype = guide_legend(order = 4, title = NULL, ncol = 1)
    ) +  # This ensures that shapes are displayed in the legend
    theme(legend.spacing.y = unit(0, "cm"))
  
  print(p)
  
  # saving within a folder called "CCLO maps" in my working directory
  ggsave(filename = file.path("CCLO maps", paste0("CCLO-", CCLOtags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
  
  readline(prompt = 'Press return/enter to continue')}

# 61028, 61569, 61587, 64603, 85818

df.61028 <- CCLO |>
  filter(motusTagID == 61028)
# Nov 8 in CO

df.85818 <- CCLO |>
  filter(motusTagID == 85818)
# Dec 4 in NM

df.67387 <- CCLO |>
  filter(motusTagID == 67387)
# oct 19 in NM

## BAIS ####

BAIS <- df.alltags.912.clean |>
  filter(speciesEN == "Baird's Sparrow")

BAIStags <- unique(BAIS$motusTagID) # 61

# this for loop will save a copy of each map in a folder I created in my working directory. edit line 1201 to change that
for(i in 1:length(BAIStags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == BAIStags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == BAIStags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
  p <- ggplot(data = world) + 
    geom_sf(colour = "black") +
    geom_sf(data = filter(allWinterRanges, speciesEN == "Baird's Sparrow"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
    geom_sf(data = lakes, colour = NA, fill = "white") +
    geom_sf(data = states_CanUS, colour = "black", fill = NA) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_bw() + 
    labs(x = "", y = "") +
    geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
    geom_path(data = filter(path, as.numeric(year(time)) == tagDepYear), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#440154") +
    geom_path(data = filter(path, as.numeric(year(time)) == (tagDepYear + 1)), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#440154") +
    geom_segment(data = filter(first_receivers, motusTagID == BAIStags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat), show.legend = FALSE, size = 1, colour = "#440154") + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
    geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
    geom_point(data = path, # plotting tag deployment sites
               aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
    geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
    geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("BAIS", BAIStags[i], sep = "-")) +
    
    # Manual color scales for different layers
    scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
    scale_colour_manual(values = c("Active Receivers" = "gray45")) +
    scale_linetype_manual(values = c("Apparent Fall Migration Route" = "solid", "Apparent Spring Migration Route" = "dotdash")) +
    
    # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1),
      linetype = guide_legend(order = 4, title = NULL, ncol = 1)
    ) +  # This ensures that shapes are displayed in the legend
    theme(legend.spacing.y = unit(0, "cm"))
  
  print(p)
  
  # saving within a folder called "BAIS maps" in my working directory
  ggsave(filename = file.path("BAIS maps", paste0("BAIS-", BAIStags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
  
  readline(prompt = 'Press return/enter to continue')}

## GRSP ####

GRSP <- df.alltags.912.clean |>
  filter(speciesEN == "Grasshopper Sparrow")

GRSPtags <- unique(GRSP$motusTagID) # 61

# this for loop will save a copy of each map in a folder I created in my working directory. edit line 1268 to change that
for(i in 1:length(GRSPtags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == GRSPtags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == GRSPtags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
  p <- ggplot(data = world) + 
    geom_sf(colour = "black") +
    geom_sf(data = filter(allWinterRanges, speciesEN == "Grasshopper Sparrow"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
    geom_sf(data = lakes, colour = NA, fill = "white") +
    geom_sf(data = states_CanUS, colour = "black", fill = NA) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_bw() + 
    labs(x = "", y = "") +
    geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
    geom_path(data = filter(path, as.numeric(year(time)) == tagDepYear), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#31688e") +
    geom_path(data = filter(path, as.numeric(year(time)) == (tagDepYear + 1)), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#31688e") +
    geom_segment(data = filter(first_receivers, motusTagID == GRSPtags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat), show.legend = FALSE, size = 1, colour = "#31688e") + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
    geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
    geom_point(data = path, # plotting tag deployment sites
               aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
    geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
    geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("GRSP", GRSPtags[i], sep = "-")) +
    
    # Manual color scales for different layers
    scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
    scale_colour_manual(values = c("Active Receivers" = "gray45")) +
    scale_linetype_manual(values = c("Apparent Fall Migration Route" = "solid", "Apparent Spring Migration Route" = "dotdash")) +
    
    # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1),
      linetype = guide_legend(order = 4, title = NULL, ncol = 1)
    ) +  # This ensures that shapes are displayed in the legend
    theme(legend.spacing.y = unit(0, "cm"))
  
  print(p)
  
  # saving within a folder called "GRSP maps" in my working directory
  ggsave(filename = file.path("GRSP maps", paste0("GRSP-", GRSPtags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
  
  readline(prompt = 'Press return/enter to continue')}


## TBLO ####

TBLO <- df.alltags.912.clean |>
  filter(speciesEN == "Thick-billed Longspur")

TBLOtags <- unique(TBLO$motusTagID) # 61

# this for loop will save a copy of each map in a folder I created in my working directory. edit line 1336 to change that
for(i in 1:length(TBLOtags)){
  
  # filtering the path data for each individual tag
  path <- filter(path912, motusTagID == TBLOtags[i])
  
  bird <- filter(df.alltags.912.clean, motusTagID == TBLOtags[i])
  
  # defining the limits of the plot - some birds were not detected near the colony so have no recv around the tag site in their path data, so I am selecting the max/min based on tag sites as well
  
  xmin <- min(path912$recvDeployLon, na.rm = TRUE) - 2
  xmax <- max(path912$recvDeployLon, na.rm = TRUE) + 2
  ymin <- min(path912$recvDeployLat, na.rm = TRUE) - 1
  ymax <- max(path912$recvDeployLat, na.rm = TRUE) + 1
  # get a warning about fractional recycling, which doesn't matter in this instance - it's doing what I want it to
  
  p <- ggplot(data = world) + 
    geom_sf(colour = "black") +
    geom_sf(data = filter(allWinterRanges, speciesEN == "Thick-billed Longspur"), color = NA, alpha = 0.3, aes(fill = as.factor(Season))) + 
    geom_sf(data = lakes, colour = NA, fill = "white") +
    geom_sf(data = states_CanUS, colour = "black", fill = NA) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_bw() + 
    labs(x = "", y = "") +
    geom_point(data = activeRecv, aes(x = recvDeployLon, y = recvDeployLat, colour = as.factor(recvType)), shape = 1) +
    geom_path(data = filter(path, as.numeric(year(time)) == tagDepYear), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#90d743") +
    geom_path(data = filter(path, as.numeric(year(time)) == (tagDepYear + 1)), 
              aes(x = recvDeployLon, y = recvDeployLat, group = motusTagID, linetype = as.factor(migrationSeason)), size = 1, colour = "#90d743") +
    geom_segment(data = filter(first_receivers, motusTagID == TBLOtags[i]), aes(x = tagDepLon, y = tagDepLat, xend = recvDeployLon, yend = recvDeployLat), show.legend = FALSE, size = 1, colour = "#90d743") + # this connects tag sites to the first recv, which is helpful for birds that were not detected around the tag site
    geom_point(data = path, aes(x = recvDeployLon, y = recvDeployLat, fill = as.factor(recvWithDets)), shape = 21) + # plotting receivers
    geom_point(data = path, # plotting tag deployment sites
               aes(x = tagDepLon, y = tagDepLat, fill = as.factor(tagSite)), colour = "black", shape = 23, size = 2) +
    geom_text(data = filter(states_CanUS, postal %in% c("MT", "ND", "WY", "SD", "NE", "KS", "OK", "NM", "TX", "CO", "WY", "FL", "ID", "UT", "AZ")), aes(x = longitude, y = latitude, label = postal)) +
    geom_text(data = filter(states_CanUS, postal %in% c("AB", "SK", "MB")), aes(label = postal), y = 51, x = c(-112.5, -105, -97.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("TBLO", TBLOtags[i], sep = "-")) +
    
    # Manual color scales for different layers
    scale_fill_manual(values = c("Nonbreeding Season Range" = "purple4", "Tag Site" = "red", "Receiver With Detections" = "black")) +
    scale_colour_manual(values = c("Active Receivers" = "gray45")) +
    scale_linetype_manual(values = c("Apparent Fall Migration Route" = "solid", "Apparent Spring Migration Route" = "dotdash")) +
    
    # Create a combined, unified legend using guides()
    guides(
      color = guide_legend(order = 1, title = "Legend", ncol = 1),
      fill = guide_legend(order = 2, title = NULL, ncol = 1),
      shape = guide_legend(order = 3, ncol = 1),
      linetype = guide_legend(order = 4, title = NULL, ncol = 1)
    ) +  # This ensures that shapes are displayed in the legend
    theme(legend.spacing.y = unit(0, "cm"))
  
  print(p)
  
  # saving within a folder called "TBLO maps" in my working directory
  ggsave(filename = file.path("TBLO maps", paste0("TBLO-", TBLOtags[i], ".tiff")), plot = p, width = 10, height = 8, dpi = 600)
  
  readline(prompt = 'Press return/enter to continue')}

#### 11 - Possible overwintering areas ####

# southernmost detections for each individual

winterLocs <- df.alltags.912.clean |>
  group_by(motusTagID, speciesEN) |>
  filter(recvDeployLat == min(recvDeployLat)) |>
  filter(ts == min(ts)) |> # first arrival on wintering grounds (or lowest latitude detection)
  select(ts, motusTagID, speciesEN, tagDepYear, bandingCode, recvDeployLat, recvDeployName) |>
  arrange(speciesEN, recvDeployLat)

# this will be investigated in further detail for deliverable 3

#### 12 - Saving files ####

# change dates in file names as desired
saveRDS(df.alltags.912.clean, "df.alltags.912.clean.depDates.2025Mar25.rds")
saveRDS(allDepDates, "allDepDates.2025Mar25.rds")
saveRDS(migPace1, "fallMigPaceData.2025Feb25.rds")
saveRDS(springDepDates, "springDepDates.2025Feb25.rds")
saveRDS(springArrivals, "springArrivalDates.2025Feb25.rds")
saveRDS(springMigPace1, "springMigPace.2025Feb25.rds")
saveRDS(closest_receivers, "closestRecv.2025Mar6.rds")

#### 13 - References ####

# Beason, R. C. (2020). Horned Lark (Eremophila alpestris), version 1.0. In Birds of the World (S. M. Billerman, Editor). Cornell Lab of Ornithology, Ithaca, NY, USA. https://doi-org.proxy.library.carleton.ca/10.2173/bow.horlar.01

# Bleho, B., K. Ellison, D. P. Hill, and L. K. Gould (2020). Chestnut-collared Longspur (Calcarius ornatus), version 1.0. In Birds of the World (P. G. Rodewald, Editor). Cornell Lab of Ornithology, Ithaca, NY, USA. https://doi-org.proxy.library.carleton.ca/10.2173/bow.chclon.01

# Davis, S. K., M. B. Robbins, and B. C. Dale (2020). Sprague's Pipit (Anthus spragueii), version 1.0. In Birds of the World (A. F. Poole, Editor). Cornell Lab of Ornithology, Ithaca, NY, USA. https://doi-org.proxy.library.carleton.ca/10.2173/bow.sprpip.01

# Ellison, K., E. McKinnon, S. Zack, S. Olimb, R. Sparks, and E. Strasser (2017). Migration and winter distribution of the Chestnut-collared Longspur. Animal Migration (4): 37-40. https://doi.org/10.1515/ami-2017-0005 

# Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson, W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, I. Davies, A. Rodewald, V. Ruiz-Gutierrez, C. Wood. 2023. eBird Status and Trends, Data Version: 2022; Released: 2023. Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.2022



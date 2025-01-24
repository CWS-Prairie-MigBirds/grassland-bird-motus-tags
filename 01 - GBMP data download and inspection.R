#### GBMP Script 01 - Data download and initial inspection ####
# GBMP = Grassland Bird Monitoring Program
# ECCC contract 3000792304
# adapted from the Motus R book by Sarah Endenburg
# This script contributes to deliverable 1 

#### 1 - Prepare working environment ####

library(motus)
library(tidyverse)
library(lubridate)
library(RSQLite)
library(rnaturalearthdata)
library(rnaturalearth)
library(osfr) # for interacting with the OSF repository I created to store the GBMP Motus data

# This code works with data loaded from the OSF repository and saves new data files there
# That requires an OSF account and personal access token (PAT)
# I just operate out of a single working directory on my local computer, you may need to specify file paths if you want to create separate folders

# Motus data is stored in UTC, set system time to match
Sys.setenv(tz = "UTC")

# Authenticate OSF account to download/upload data
# paste your OSF PAT within the quotation marks below to authenticate 
osf_auth(token = "")

project <- osf_retrieve_node("ge8bd")
# the code in the quotation marks is the project id on OSF

print(osf_ls_files(project)) # view all of the files

#### 2 - Checking tag deployment data ####
# tag csv file downloaded from motus.org for project 912 (Prairie Landbird Monitoring - Canada)
# saved on the GBMP OSF repository
files <- osf_ls_files(project)

file_to_download <- files[files$name == "MotusTagsTemplate-912.csv"]
# retrieving the file object (just using file name doesn't work, there might be multiples of the same name)

osf_download(file_to_download, conflicts = "rename") # renames the file if there is already one with the same name in the destination folder

# now have to read the downloaded file into R
tags912 <- read_csv("MotusTagsTemplate-912.csv")

# filtering for tags that have deployment dates 
depTags912 <- tags912 |>
  filter(is.na(utcYearStart) == FALSE)
# 241 tags with deployment dates

table(depTags912$speciesName)
# six species - Ammodramus savannarum (grasshopper sparrow, GRSP), Anthus spragueii (Sprague's pipit, SPPI), Calcarius ornatus (Chestnut-collared longspur, CCLO), Centronyx bairdii (Baird's sparrow, BAIS), Eremophila alpestris (Horned lark, HOLA), Rhynchophanes mccownii (Thick-billed longspur, TBLO)

#### 3 - downloading SQL detection database ####
# this section and section 4 can be skipped unless you want to re-download the data. under step 5, I will load the unfiltered data from the OSF repository

sql.motus.912 <- tagme(projRecv = 912, forceMeta = TRUE)

# logout
motusLogout()

#### 4 - converting SQL alltags table to a flat dataframe for easier use/storage ####

tbl.alltags <- tbl(sql.motus.912, "alltags")

df.alltags.912 <- tbl.alltags %>% 
  collect() %>% 
  as.data.frame() %>% 
  mutate(
    # adjust time variables to yyyy-mm-dd format, in a new column to keep the original ts
    ts.orig = ts,
    ts = as_datetime(ts, tz = "UTC", origin = "1970-01-01"),
    tagDeployStart = as_datetime(tagDeployStart, origin = "1970-01-01"),
    
    # create variable for tag deployment year
    tagDepYear = year(tagDeployStart),
    
    # round lat/lon coordinates to 4 decimal places
    recvDeployLat = plyr::round_any(recvDeployLat, 0.0001),
    recvDeployLon = plyr::round_any(recvDeployLon, 0.0001),
    
  )

saveRDS(df.alltags.912, "df.alltags.912.2025Jan8.rds")

# 7 828 574 detections in the unprocessed data

#### 5 - removing short runs ####

# first, I will load the unfiltered data from OSF

file_to_download <- files %>%
  filter(name == "df.alltags.912.2025Jan8.rds")

# retrieving the file object (just using file name doesn't work, there might be multiples of the same name)

osf_download(file_to_download, conflicts = "rename") # renames the file if there is already one with the same name in the destination folder

# now have to read the downloaded file into R
df.alltags.912 <- readRDS("df.alltags.912.2025Jan8.rds")

# runs of <3 hits are more likely to be junk data for Lotek tags, and "runs" of 1 are likely to be false for CTT tags. short runs are also difficult to verify. the motus filter considers any runs <5 to be false, but many of these lotek tags have long burst intervals, making it more likely that short runs will occur, especially if the bird is at the limit of the receiver detection range. that's why I kept runs of 3-4, but anything smaller gets hard to tell if it's real or not. They can be revisited later if desired

df.alltags.912.filter <- df.alltags.912 |>
  filter(!(tagModel %in% c("NTQB2-2-M", "NTQB2-3-2-M", "NTQB2-5-1-M") & runLen < 3) & !(tagModel %in% c("HybridTag", "LifeTag") & runLen == 1))

nrow(df.alltags.912) - nrow(df.alltags.912.filter)
# removed 127 068 hits

# pulling the short runs to save separately for inspection later if desired
short_runs_912 <- df.alltags.912 |>
  filter((tagModel %in% c("NTQB2-2-M", "NTQB2-3-2-M", "NTQB2-5-1-M") & runLen < 3) | (tagModel %in% c("HybridTag", "LifeTag") & runLen == 1)) |>
  mutate(reasonRemoved = "shortRun")

#### 6 - removing data for undeployed tags, pre-deployment detections & detections past expected tag lifespan ####
# for some reason, the majority of this dataset seems to be for undeployed tags or detections outside of deployments. I will focus on detections for tags that were truly deployed

dets_undep_tags <- df.alltags.912.filter |>
  filter(!motusTagID %in% depTags912$tagID) |>
  mutate(reasonRemoved = "undeployedTag")

nrow(dets_undep_tags) # 4 843 670

df.alltags.912.filter <- df.alltags.912.filter |>
  filter(motusTagID %in% depTags912$tagID)
# 2 859 064 detections for deployed tags. 

unique(df.alltags.912.filter$motusTagID)
# 206 motusTagIDs in there at the moment

# there were some 2021 detections despite all of the tags being deployed in 2022 or later (I think some were registered in 2021 and that's why they were "detected"...)
df.alltags.912.filter <- df.alltags.912.filter |>
  mutate(year = year(ts))

table(df.alltags.912.filter$year) # 241 286 hits from 2021. I will save those and remove

detections2021 <- df.alltags.912.filter |>
  filter(year == 2021) |>
  mutate(reasonRemoved = as.character("2021detection"))

df.alltags.912.filter <- df.alltags.912.filter |>
  filter(!year == 2021)

# I think there are still some pre-deployment detections in there that need to be removed though
# I think most of them show up without tag deployment info, some of them have the info for whatever reason

tagDeployInfo <- tags912 |>
  select(tagID, utcYearStart, utcMonthStart, utcDayStart, utcHourStart, utcMinuteStart) |>
  rename(motusTagID = tagID)

df.alltags.912.filter <- df.alltags.912.filter |>
  left_join(tagDeployInfo, by = "motusTagID") |>
  mutate(tagDeployDate = make_datetime(
    utcYearStart,
    utcMonthStart,
    utcDayStart,
    utcHourStart,
    utcMinuteStart
  ))

allDetsPreDeployment <- df.alltags.912.filter |>
  filter(ts < tagDeployDate) |>
  mutate(reasonRemoved = "preDeployment")
# 2 192 905 are pre-deployment

preDepDets_withInfo <- df.alltags.912.filter |>
  filter(ts < tagDeployDate & !is.na(tagDeployStart)) |>
  select(motusTagID, ts, hitID, runID, recvDeployName, tagDeployDate)
# 135 pre-deployment detections do have information
# they're all garbage though - no recv info, seem to all have multiple duplicates
# they're all for 75830 and all on 2023-06-12 when tag was deployed on 2023-06-13. weird, will put aside to review later if desired. it's a CTT tag too, which is supposed to have fewer false detections, so maybe it was a testing situation.

# and removing from the main dataset

df.alltags.912.filter2 <- df.alltags.912.filter |>
  filter(ts >= tagDeployDate)
# 424 873 detections left

# now checking to see if tag info is there for post-deployment detections

postDepDets_missingTagInfo <- df.alltags.912.filter2 |>
  filter(is.na(tagDeployStart) | is.na(tagDepLat)) |>
  mutate(reasonRemoved = "missingTagInfo")
# there are 72 detections that are missing tag deployment info but are post-deployment dates.

table(postDepDets_missingTagInfo$recvDeployName) 
# some might be legit, at Battlecreek/Nashlyn/Snake Butte/Lake Seventeen.
# the rest are towers I've had issues with before - Audubon Center at Riverlands, Facultad de Ciencias Agropecuarias, Summit Fire Tower, Alaksen

table(postDepDets_missingTagInfo$tagModel)
# how long do these tags last? they are all NTQB2-5-1-M and NTQB2-2-M

table(postDepDets_missingTagInfo$tagBI) # 27 - 30.7s BIs

table(postDepDets_missingTagInfo$runLen) # mostly runs of three but some are longer

table(postDepDets_missingTagInfo$motusTagID) 

# 67390 - dets are sept 29 2023, deployed in july 2022...oh, was expected only to last until Sept 28 2023
# 75408 dets in 2023 july, aug, nov, all post-deployment (dep june 2023). was expected to last until sept 2024 so I'm not sure why theyre not showing the info. 
# 75410 dets in sept 2024, deployed in july 2023. tag was only expected to last until sept 25 and hits are on the 26
# 74916 - hits are post-deployment
# I'm going to remove them and put them aside for later inspection if desired

df.alltags.912.filter2 <- df.alltags.912.filter2 |>
  filter(!is.na(tagDeployStart) & !is.na(tagDepLat))

# tag expected lifespan depends on the model and burst interval 

df.alltags.912.filter2$tagBIround <- round(df.alltags.912.filter2$tagBI, digits = 1)

table(df.alltags.912.filter2$tagBIround)
# 5 (CTT tags), 27.1, 27.7, 20.7, 34.9

df.alltags.912.filter2 |>
  group_by(tagModel, tagBIround) |>
  summarize(nTags = n_distinct(motusTagID)) # 9 combinations of BI and tagModel

# I got the expected lifespans (in days) for each tag model and BI combo from the Motus portal. I used the "estimated buffered total lifespan" which is pretty generous, there's also a lower number (estimated total lifespan). I think it's hard to tell how long tags actually last since we don't get many back and we don't know if a lack of detections is because of the tag or where the bird is relative to receivers, etc. I'll put a column for the lower number in case we want the cutoff to be more stringent

df.alltags.912.filter2 <- df.alltags.912.filter2 |>
  mutate(bufferedTagLife = case_when(
    tagModel == "HybridTag" ~ as.numeric("1825"),
    tagModel == "LifeTag" ~ as.numeric("1825"),
    (tagModel == "NTQB2-2-M" & tagBIround == 27.1) ~ as.numeric("230"),
    (tagModel == "NTQB2-2-M" & tagBIround == 27.7) ~ as.numeric("232"),
    (tagModel == "NTQB2-2-M" & tagBIround == 34.9) ~ as.numeric("261"),
    tagModel == "NTQB2-3-2-M" ~ as.numeric("537"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 27.1) ~ as.numeric("421"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 27.7) ~ as.numeric("426"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 30.7) ~ as.numeric("449"),
    TRUE ~ as.numeric("0")
  ))

# the two lifespan estimates are identical for the CTT tags
df.alltags.912.filter2 <- df.alltags.912.filter2 |>
  mutate(tagLife = case_when(
    tagModel == "HybridTag" ~ as.numeric("1825"),
    tagModel == "LifeTag" ~ as.numeric("1825"),
    (tagModel == "NTQB2-2-M" & tagBIround == 27.1) ~ as.numeric("153"),
    (tagModel == "NTQB2-2-M" & tagBIround == 27.7) ~ as.numeric("155"),
    (tagModel == "NTQB2-2-M" & tagBIround == 34.9) ~ as.numeric("174"),
    tagModel == "NTQB2-3-2-M" ~ as.numeric("358"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 27.1) ~ as.numeric("280"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 27.7) ~ as.numeric("284"),
    (tagModel == "NTQB2-5-1-M" & tagBIround == 30.7) ~ as.numeric("299"),
    TRUE ~ as.numeric("0")
  ))

table(df.alltags.912.filter2$tagLife)

# now compare to how long the tags were detected for

df.alltags.912.filter2 <- df.alltags.912.filter2 |>
  group_by(motusTagID) |>
  mutate(numDaysDetected = as.numeric(difftime(max(ts), min(ts), units = "days")), # should this be since deployment? not all tags were detected right away when deployed
         daysSinceDep = as.numeric(difftime(ts, tagDeployDate, units = "days")))

lateHits <- df.alltags.912.filter2 |>
  filter(daysSinceDep > bufferedTagLife) |>
  mutate(reasonRemoved = "pastBufferedTagLife")
# there are way more late hits when you use the smaller estimate for tag lifespan, as expected.
# actually there are none that are past the buffered lifespan estimate

# maybe can do something with this later on - like for noisy towers (<50% hits noisy?), remove hits outside of the shorter estimated tag lifespan. 

#### 7 - removing duplicates ####

# Motus data often has duplicate detections that are identical except for the hitID, I think it's an artifact of the data upload or download process...
# First, I am pulling both the "original" and duplicated rows to make sure they're actually identical before removing the duplicates

# Identify duplicates based on the specified columns
duplicate_rows <- df.alltags.912.filter2[duplicated(df.alltags.912.filter2[, c("mfgID", "ts", "sig", "runLen", "recvDeployName", "motusTagID", "port")]) | 
                                          duplicated(df.alltags.912.filter2[, c("mfgID", "ts", "sig", "runLen", "recvDeployName", "motusTagID", "port")], fromLast = TRUE), ] |>
  arrange(ts, runID)
# I include port here because it is theoretically possible for a tag to have simultaneous hits on two ports of a receiver, which would result in seemingly identical or duplicated detections

view(duplicate_rows)
# there are 340 duplicated rows 
# i ordered the df by ts and runID and did a spot check. as expected, the duplicates seem to all be identical except for the hitID, runID and bootNum. So I think the duplicates got created when data were uploaded with different boot numbers, and Motus assigned them different hitID and runIDs but the data are otherwise identical
# so, I will remove the duplicates and save them separately in case we want to

# and remove one of each duplicate pair
dups <- df.alltags.912.filter2 %>% 
  arrange(batchID) %>%
  dplyr::select(mfgID, ts, sig, runLen, recvDeployName, motusTagID, port) %>%
  duplicated(fromLast = TRUE)

df.alltags.912.filter.nodups <- df.alltags.912.filter2[!dups,]

duplicateData <- df.alltags.912.filter2[dups,] |>
  mutate(reasonRemoved = "duplicate")

nrow(df.alltags.912.filter2) - nrow(df.alltags.912.filter.nodups)
# removed 170 detections, which is 1/2 the size of the duplicate pair file

#### 8 - checking for ambiguous detections ####

table(df.alltags.912.filter.nodups$ambigID) # none. that was easy.

#### 9 - checking recv deploy lat and long ####
# if receivers are missing this info, it's not useful in analysis/can indicate a missing receiver deployment in the system

sum(is.na(df.alltags.912.filter.nodups$recvDeployLat))
# 12  have NA for recvDeployLat

noRecvInfo <- df.alltags.912.filter.nodups |>
  filter(is.na(recvDeployLat)) |>
  mutate(reasonRemoved = "noRecvInfo")
# some have absolutely no recv info but most have at least the recv ID (antenna serial number)

table(noRecvInfo$recv)

# CTT-V3023D0DEA8D - 6 hits
# escape coulee...deployed by CWS - Prairie Region Stations
# the hits without info are from just a few minutes before the official deployment date/time of the receiver.

# the rest don't even have a recv serial number
# I'll remove all - not many hits lost and if the recv was not actually deployed, we don't necessarily know if the location is correct

df.alltags.912.filter.nodups <- df.alltags.912.filter.nodups |>
  filter(is.na(recvDeployLat) == FALSE)

#### 10 - map of tag deployment sites ####

tagDepSites <- df.alltags.912.filter.nodups |>
  ungroup() |>
  select(tagDepLat, tagDepLon, speciesEN) |>
  mutate(tagDepLat = round(tagDepLat, digits = 1),
         tagDepLon = round(tagDepLon, digits = 1)) |>
  distinct() |>
  filter(is.na(tagDepLat) == FALSE)

world <- ne_countries(scale = "medium", returnclass = "sf")

states_Can <- ne_states(country = c("canada"), returnclass = "sf")

lakes <- ne_download(scale = "large", type = 'lakes', category = 'physical',
returnclass = "sf", destdir = "mapData") # only need this first time downloading

#lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
                 #returnclass = "sf",
                 #destdir = "mapData") # use this if already downloaded shapefiles

xmax <- -90
xmin <- -120
ymax <- 60
ymin <- 45

ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = lakes, fill = "lightblue3", colour = "gray44") +
  geom_sf(data = states_Can, colour = "gray44")+
  geom_sf(data = world, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  labs(x = "", y = "") +
  geom_jitter(data = tagDepSites, aes(x = tagDepLon, y = tagDepLat, fill = as.factor(speciesEN)), colour = "black", shape = 21, size = 3)+
  geom_text(data = states_Can, aes(label = postal, x = longitude, y = latitude), size = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue3")) +
  facet_wrap(~speciesEN) +
  theme(legend.position = "none") +
  scale_x_continuous(
    breaks = c(-115, -110, -105, -100, -95),  # Manually specify the breaks (in degrees)
    labels = c("115°W", "110°W", "105°W", "100°W", "95°W")  # Custom labels for the breaks
  ) +
  # Rotate x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### 11 - pulling recv & tag site lat and long to plot on Google Maps ####

receivers <- df.alltags.912 |>
  dplyr::select(recvDeployLat, recvDeployLon, recvDeployName) |>
  distinct()
# 1113 recv in unfiltered data
  
# already pulled tagDepSites earlier

# the map can be viewed here (anyone with the link can view but the map cannot be searched on the internet)
# https://www.google.com/maps/d/u/0/edit?mid=1YUYdqCOwcfgHutPt7zcVWtkaobgEgWw&usp=sharing

#### 12 -  saving files ####

# main dataframes, first has duplicates, the second doesn't. there are on OSF already
saveRDS(df.alltags.912.filter2, "df.alltags.filter.2025Jan20.rds")
saveRDS(df.alltags.912.filter.nodups, "df.alltags.filter.nodups.2025Jan20.rds")

# tag deployment locations & recv with detections
write_csv(tagDepSites, "GBMP_tagDeploymentSites.2025Jan7.csv")
write_csv(receivers, "recv.LatLon.2025Jan9.csv")

# data that I removed
saveRDS(shortRuns, "shortRuns912.2025Jan20.rds")
saveRDS(dets_undep_tags, "detsUndepTags.2025Jan20.rds")
saveRDS(detections2021, "2021Detections.2025Jan15.rds")
saveRDS(allDetsPreDeployment, "preDeploymentTagDets.2025Jan15.rds")
saveRDS(postDepDets_missingTagInfo, "postDepDets_missingTagInfo.2025Jan15.rds")
saveRDS(duplicateData, "duplicateDataRemoved.2025Jan20.rds")
saveRDS(noRecvInfo, "noRecvInfo.2025Jan20.rds")

# and then can upload to OSF - note that if the file already exists there (if it's one of the ones I uploaded there), you may have to specify conflict = "rename" to upload. I just put the upload option here in case you want it 
local_files <- c("df.alltags.filter.2025Jan20.rds", "df.alltags.filter.nodups.2025Jan20.rds", "shortRuns912.2024Jan15.rds","detsUndepTags.2025Jan20.rds", "2021Detections.2025Jan15.rds", "preDeploymentTagDets.2025Jan15.rds", "postDepDets_missingTagInfo.2025Jan15.rds", "duplicateDataRemoved.2025Jan20.rds", "noRecvInfo.2025Jan20.rds", "GBMP_tagDeploymentSites.2025Jan7.csv", "recv.LatLon.2025Jan9.csv")

# Upload the file to the project (path argument specifies folder in the project, leaving blank uploads it to the main directory)
osf_upload(project, local_files)


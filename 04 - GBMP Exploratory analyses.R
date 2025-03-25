#### GBMP Script 04 - Exploratory analyses ####
# GBMP = Grassland Bird Monitoring Program
# ECCC contract 3000792304
# adapted from the Motus R book by Sarah Endenburg
# This script contributes to deliverable 3

#### 1 - Prepare working environment ####
library(tidyverse)
library(lubridate)
library(rnaturalearthdata) # plotting maps
library(rnaturalearth)
library(geosphere) # calculating distance
library(suncalc) # to calculate solar elevation/altitude
library(lutz) # for determining time zones
library(DHARMa) # for checking model residuals
library(MuMIn) # for comparing model AIC 

# Motus data is stored in UTC, set system time to match
Sys.setenv(tz = "UTC")

# turning off scientific notation
options(scipen = 9999)

#### 2 - load data ####
# ensure that correct file version is used or columns will get duplicated
df.alltags.912.clean <- readRDS("df.alltags.912.clean.depDates.2025Mar25.rds")

#### 3 - possible overwintering areas ####

# were there any detections between November-February? Ellison et al. 2017 found that CCLO arrived on wintering grounds starting in November

df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(month = month(ts))

table(df.alltags.912.clean$month) # there are no Jan or Feb hits, only November and December

winterDets <- df.alltags.912.clean |>
  filter(month %in% c(11, 12)) 

winterBirds <- winterDets |>
  select(motusTagID, speciesEN) |>
  distinct() # 12 bird when I include Nov, only three when I just look at December

# this shows all detections for each of the 12 birds with winter detections
for(i in 1:nrow(winterBirds)){
  print(ggplot(data = filter(df.alltags.912.clean, 
                             motusTagID == winterBirds$motusTagID[i]), 
               aes(x = ts, y = recvDeployLat)) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
          geom_point(aes(size = runLen, color = recvDeployName)) +
          labs(title = paste(winterBirds$motusTagID[i], winterBirds$speciesEN[i], sep = "-")))
  readline(prompt = 'Press return/enter to continue')}

# this loops through the winter detections and plots them (only the winter detections, plotted by receiver and motusTagID to look at signal strength)
for(tag_id in unique(winterBirds$motusTagID)) {
  
  # Filter the data based on the current motusTagID
  filtered_data <- filter(winterDets, motusTagID == tag_id)
  
  # Loop through each unique recvDeployName
  for(deploy_name in unique(filtered_data$recvDeployName)) {
    
    # Filter the data for the current recvDeployName
    deploy_data <- filter(filtered_data, recvDeployName == deploy_name)
    
    # Extract species name (assuming speciesEN is a column in the data)
    species <- unique(deploy_data$speciesEN)
    
    # Create the plot
    p <- ggplot(deploy_data, aes(x = ts, y = sig)) + 
      theme_bw() + 
      geom_point(aes(size = runLen)) + 
      labs(
        title = paste(tag_id, species, deploy_name, sep = "-"), 
        x = "ts", 
        y = "sig"
      )
    
    # Print the plot
    print(p)
    
    # Wait for the user to press Enter before continuing to the next plot
    readline(prompt = "Press Enter to see the next plot...")
  }
}

# 61569 at GPCD SPLT looks like a flyby - all within about a 10 minute span
# 75828 at Sterling Wildlife Management Area-ID - detections between Nov 10 and Dec 2
# 75828 at Minidoka National Wildlife Refuge, ID Nov 20 & Dec 1, Cold Water Rest Area-ID Dec 1 & 16, Curlew National Grassland (space after the name) looks like a flyby all within two minutes, strling wildlife management area between Nov 10 - Dec 2, market lake WMA all within a 5min period, deer parks also looks like a flyby, 
# 85818 at GPCD Bosque del Apache NWR - single run of three
# 61024 at rancho el uno - looks like a flyby
# 61587 at GPCD SPLT looks like a flyby, same with 61606 at the same recv
# 61617 at ruth bader ginsburg is a run of 4
# 75412 at cheyenne bottoms - a few runs but all within about 20min of one another
# 75413 at mimms ranch - flyby
# 85088 at nashyn in Nov?? and battlecreek? battlecreek looks good, nashlyn is just a run of three. late departure? also can I filter these out of this analysis using eBird wintering areas? because this is unlikely the wintering area (HOLA). I had filtered this out of depDate analyses because the tag was detected for <24h
# 85753 at salt plains is a single run of 8

# plotting out these winter detections

unique(winterDets$speciesEN) # SPPI, CCLO and HOLA only

world <- ne_countries(scale = "medium", returnclass = "sf")

states_CanUS <- ne_states(country = c("canada", "united states of america"), returnclass = "sf")

lakes <- ne_download(scale = "large", type = 'lakes', category = 'physical',
                     returnclass = "sf", destdir = "mapData")

CCLO_range <- st_read("chclon_range_2022/chclon_range_2022.gpkg")
CCLO_winter_range <- CCLO_range |>
  filter(season == "nonbreeding")

SPPI_range <- st_read("sprpip_range_2022/sprpip_range_2022.gpkg")
SPPI_winter_range <- SPPI_range |>
  filter(season == "nonbreeding")

HOLA_range <- st_read("horlar_range_2022/horlar_range_2022.gpkg")
HOLA_winter_range <- HOLA_range |>
  filter(season == "nonbreeding")

allWinterRanges <- CCLO_winter_range |>
  bind_rows(SPPI_winter_range, HOLA_winter_range) |>
  rename(speciesEN = common_name) |>
  mutate(season2 = "Nonbreeding season range")


xmin <- min(winterDets$recvDeployLon, na.rm = TRUE) - 10
xmax <- max(winterDets$recvDeployLon, na.rm = TRUE) + 10
ymin <- min(winterDets$recvDeployLat, na.rm = TRUE) - 10
ymax <- max(winterDets$recvDeployLat, na.rm = TRUE) + 10

winterDets <- winterDets |>
  mutate(newTagID = paste(motusTagID, bandingCode, sep = "-"))

ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = allWinterRanges, color = NA, alpha = 0.3, aes(fill = as.factor(season2))) +
  geom_sf(data = lakes, colour = NA, fill = "white") +
  geom_sf(data = states_CanUS, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  geom_text(data = states_CanUS, aes(x = longitude, y = latitude, label = postal)) + 
  geom_point(data = winterDets, aes(x = tagDepLon, y = tagDepLat), colour = "black", fill = "red", shape = 23) +
  geom_point(data = winterDets, aes(x = recvDeployLon, y = recvDeployLat, size = as.numeric(runLen), fill = as.factor(newTagID)), shape = 21) +
  labs(x = "", y = "", title = "Winter (November & December) detections", size = "Run length") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(order = 2, title = "Legend", ncol = 1)) +
  facet_wrap(~speciesEN)
# not sure that this provides much additional information from the track maps, aside from knowing that these detections were in November or December


#### 4 - possible stopover areas ####

# this code identifies possible stopover areas based on length of time spent there (post-departure from the breeding grounds)
# Hagelin et al used 2+ days as length cutoff, akesson et al 2012 too
# Hill et al 2018 used 2+ consecutive twilights
# Delmore, Stutchbury et al also use 2 days
# I shortened the length of time to 24h to see if that identified more sites, I don't think it changed much from using 48h

# let's look at any post-departure detections where the bird was detected there for 1+ days
# there are sometimes receivers within ~70km or so of one another that result in clusters of detections, so I want to account for that (I tried 50km first but some individuals had clusters ~70km so I expanded the distance range)
# BAIS, TBLO and GRSP didn't have anything (when I checked) so I will just do CCLO, SPPI and HOLA

# I do this by species since the code can take a long time to run. just change the speciesEN below and save the valid_detections_df below before running the next species
df_filtered2 <- df.alltags.912.clean %>%
  filter(speciesEN == "Chestnut-collared Longspur" & ts > depDate100km)

# Function to calculate distance between two points (Vincenty formula)
calculate_distance <- function(lat1, lon1, lat2, lon2) {
  distVincentyEllipsoid(c(lon1, lat1), c(lon2, lat2))} # Returns distance in meters

# create an empty list to store detections in
valid_detections <- list()

for(tag_id in unique(df_filtered2$motusTagID)) {
  
  # Filter data for the current bird (motusTagID)
  bird_data <- filter(df_filtered2, motusTagID == tag_id) %>%
    group_by(runID) %>%
    slice_head(n = 1) %>%  # just keeping the first hit in each run to reduce the amount of data I'm using
    ungroup() %>%
    group_by(motusTagID, recvDeployName, ts) %>%
    filter(sig == max(sig)) %>%
    ungroup() %>%
    arrange(ts)
  
  # Create an empty list to track groupings
  groupings <- list()
  assigned_detections <- c()  # Track detections already assigned to a group
  
  # Loop through each detection and check its distance from others
  for(i in 1:nrow(bird_data)) {
    # Skip if this detection has already been assigned to a group
    if (i %in% assigned_detections) next
    
    current_row <- bird_data[i, ]
    
    # Initialize a group with the current detection
    group <- list(current_row)
    
    # Check if the current detection is within 50km of any other detection
    for(j in 1:nrow(bird_data)) {
      if(i != j && !(j %in% assigned_detections)) { # each run should only be used in a single group
        other_row <- bird_data[j, ]
        
        # Calculate the distance between the two detections if the receiver names are different
        if (current_row$year == other_row$year) {  # current_row$recvDeployName != other_row$recvDeployName & 
          dist <- calculate_distance(
            current_row$recvDeployLat, current_row$recvDeployLon,
            other_row$recvDeployLat, other_row$recvDeployLon
          )
          
          # If within 70 km, add to the same group
          if(dist <= 70000) {
            group <- append(group, list(other_row))
            assigned_detections <- c(assigned_detections, j)  # Mark detection as assigned
          }
        }
      }
    }
    
    # Only keep groups that have detections spanning at least 24 hours
    group_df <- do.call(rbind, group)
    min_ts <- min(group_df$ts)
    max_ts <- max(group_df$ts)
    
    # Calculate the time difference between the first and last detection
    time_diff_hours <- as.numeric(difftime(max_ts, min_ts, units = "hours"))
    
    # If time span is at least 24 hours, store the group (can adjust to any desired timespan)
    if(time_diff_hours >= 24) {
      # Create a unique group identifier (can use a combination of tag_id and min_ts)
      group_id <- paste(tag_id, min_ts, sep = "-")
      
      # Add the group_id to each detection in the group
      group_df$group_id <- group_id
      
      # Store the group with the unique group_id
      groupings[[group_id]] <- group_df
    }
  }
  
  # After grouping, add valid groupings to final results
  for(group_name in names(groupings)) {
    # If valid_detections is a data frame, we need to bind the results properly
    if (exists("valid_detections")) {
      valid_detections <- rbind(valid_detections, groupings[[group_name]])
    } else {
      # If it's the first iteration, initialize the valid_detections data frame
      valid_detections <- groupings[[group_name]]
    }
  }
}

# Combine the valid detections into one data frame
valid_detections_df <- bind_rows(valid_detections) %>%
  group_by(motusTagID, group_id) %>%
  arrange(motusTagID, group_id, ts)

# Calculate time difference from the first detection for each motusTagID
valid_detections_df <- valid_detections_df %>%
  group_by(group_id) %>%
  mutate(
    first_detection_time = first(ts),
    time_diff_hours = as.numeric(difftime(ts, first_detection_time, units = "hours"))
  )

# Calculate distance to previous receiver only if the receiver names are different
valid_detections_df <- valid_detections_df %>%
  ungroup() %>%
  arrange(motusTagID, ts) %>%
  mutate(
    distance_from_previous_recv = c(NA, sapply(2:nrow(valid_detections_df), function(i) {
      if(valid_detections_df$motusTagID[i] == valid_detections_df$motusTagID[i - 1] &
         valid_detections_df$recvDeployName[i] != valid_detections_df$recvDeployName[i - 1] ) {
        calculate_distance(
          valid_detections_df$recvDeployLat[i], valid_detections_df$recvDeployLon[i],
          valid_detections_df$recvDeployLat[i - 1], valid_detections_df$recvDeployLon[i - 1]
        )
      } else {
        NA  # No distance calculation if receiver names are the same
      }
    }))
  )

# just selecting the relevant columns to make the data easier to look at
valid_detections_df <- valid_detections_df |>
  dplyr::select(ts, runID, motusTagID, tagDepLat, tagDepLon, group_id, recvDeployLat, recvDeployLon, recvDeployName, speciesEN, year, daysSinceDep, date, depDate100km, distBtwRecvKM, distance_from_previous_recv, time_diff_hours, timeDiff)

# save each one as a different object so that I can re-run the code for the other species without losing the data
valid_detections_df_CCLO <- valid_detections_df
valid_detections_df_HOLA <- valid_detections_df
valid_detections_df_SPPI <- valid_detections_df

## plotting stopover sites ####

# basemap layers already loaded above

### HOLA ####
# adding some columns for plotting purposes
plotDataHOLA <- valid_detections_df_HOLA |>
  group_by(group_id) |>
  mutate(stopoverLengthDays = difftime(max(ts), min(ts), units = "days"),
         stopoverLengthRound = round(stopoverLengthDays, digits = 0),
         meanLat = mean(recvDeployLat),
         meanLon = mean(recvDeployLon),
         min_ts = as.Date(min(ts)), # adding the start and end points for plotting
         max_ts = as.Date(max(ts))) |>
  filter(recvDeployLat < 49) |> # can filter any way you'd like, I'm just trying to remove breeding ground detections from the following spring
  mutate(siteType = "Tag site")

xmin <- min(plotDataHOLA$meanLon, na.rm = TRUE) - 10
xmax <- max(plotDataHOLA$meanLon, na.rm = TRUE) + 10
ymin <- min(plotDataHOLA$meanLat, na.rm = TRUE) - 10
ymax <- max(plotDataHOLA$meanLat, na.rm = TRUE) + 10

ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = filter(allWinterRanges, speciesEN == "Horned Lark"), color = NA, alpha = 0.3, aes(fill = as.factor(season2))) +
  geom_sf(data = lakes, colour = NA, fill = "white") +
  geom_sf(data = states_CanUS, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  geom_point(data = plotDataHOLA, aes(x = tagDepLon, y = tagDepLat), colour = "black", fill = "red", shape = 23) +
  geom_point(data = plotDataHOLA, aes(x = meanLon, y = meanLat, size = as.numeric(stopoverLengthRound), fill = as.factor(motusTagID)), shape = 21, alpha = 0.5) +
  labs(x = "", y = "", title = "Horned Lark Stopover Sites") +
  geom_text(data = states_CanUS, aes(x = longitude, y = latitude, label = postal)) + 
  geom_label(data = plotDataHOLA, aes(x = tagDepLon, y = tagDepLat), label = "Tag Site", colour = "black", size = 3, hjust = 0.5, vjust = -0.5) +
  geom_label(data = filter(plotDataHOLA, motusTagID == 75828), aes(label = paste("Date range:", min_ts, "-", max_ts)), color = "black", size = 3, vjust = -1, y = 40, x = -113) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Adjust the legend for stopoverLengthDays to range from 17 (shortest stopover) to 25 (longest stopover) by 2
  scale_size_continuous(name = "Stopover Length (days)", 
                        breaks = seq(min(as.numeric(plotDataHOLA$stopoverLengthRound)), max(as.numeric(plotDataHOLA$stopoverLengthRound)), by = 2),  # Define the breaks for the legend
                        labels = seq(min(as.numeric(plotDataHOLA$stopoverLengthRound)), max(as.numeric(plotDataHOLA$stopoverLengthRound)), by = 2),  # Define the labels for the legend
                        range = c(2, 10)) +  # Set the size range
  scale_fill_manual(values = c("Nonbreeding season range" = "purple4", "75828" = "#21918c")) +
  guides(fill = guide_legend(order = 2, title = "Legend", ncol = 1))

### CCLO ####
plotDataCCLO <- valid_detections_df_CCLO |>
  group_by(group_id) |>
  mutate(stopoverLengthDays = difftime(max(ts), min(ts), units = "days"),
         stopoverLengthRound = round(stopoverLengthDays, digits = 0),
         meanLat = mean(recvDeployLat),
         meanLon = mean(recvDeployLon),
         min_ts = as.Date(min(ts)), # adding the start and end points for plotting
         max_ts = as.Date(max(ts))) |>
  filter(recvDeployLat < 49) |>
  mutate(siteType = "Tag site")

xmin <- min(plotDataCCLO$meanLon, na.rm = TRUE) - 10
xmax <- max(plotDataCCLO$meanLon, na.rm = TRUE) + 10
ymin <- min(plotDataCCLO$meanLat, na.rm = TRUE) - 10
ymax <- max(plotDataCCLO$meanLat, na.rm = TRUE) + 10

ggplot(data = world) + 
  geom_sf(colour = "black") +
  geom_sf(data = filter(allWinterRanges, speciesEN == "Chestnut-collared Longspur"), color = NA, alpha = 0.3, aes(fill = as.factor(season2))) + 
  geom_sf(data = lakes, colour = NA, fill = "white") +
  geom_sf(data = states_CanUS, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  geom_point(data = plotDataCCLO, aes(x = tagDepLon, y = tagDepLat), fill = "red", colour = "black", shape = 23) +
  geom_point(data = filter(plotDataCCLO, !motusTagID == 61596), aes(x = meanLon, y = meanLat, size = as.numeric(stopoverLengthRound), fill = as.factor(motusTagID)), shape = 21, alpha = 0.5) +
  geom_point(data = filter(plotDataCCLO, motusTagID == 61596), aes(x = meanLon, y = meanLat, size = as.numeric(stopoverLengthRound), fill = as.factor(motusTagID)), shape = 21, alpha = 0.5) +
  labs(x = "", y = "", title = "Chestnut-collared Longspur Stopover and/or Wintering Sites", fill = "Tag ID") +
  geom_text(data = states_CanUS, aes(x = longitude, y = latitude, label = postal)) + 
  geom_label(data = filter(plotDataCCLO, motusTagID == 61596), aes(y = tagDepLat), x = -108.5, label = "Tag Sites", colour = "black", size = 3, hjust = 0.75, vjust = -0.5) +
  geom_label(data = filter(plotDataCCLO, motusTagID == 61587), aes(x = meanLon, label = paste("Date range:", min_ts, "-", max_ts)), color = "black", size = 3, vjust = -1, y = 39) + 
  geom_segment(data = filter(plotDataCCLO, motusTagID == 61587), aes(x = meanLon, y = 38.7, 
                                                                     xend = meanLon, yend = 39.5),
               color = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_label(data = filter(plotDataCCLO, motusTagID == 67508), aes(y = meanLat, label = paste("Date range:", min_ts, "-", max_ts)), color = "black", size = 3, vjust = -1, x = -100) +
  geom_segment(data = filter(plotDataCCLO, motusTagID == 67508), x = -102, y = 48.6, 
                                                                     xend = -107.5, yend = 48.1,
               color = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_label(data = filter(plotDataCCLO, motusTagID == 61596), aes(x = meanLon, label = paste("Date range:", min_ts, "-", max_ts)), color = "black", size = 3, vjust = -1, y = 45) +
  geom_segment(data = filter(plotDataCCLO, motusTagID == 61596), aes(x = meanLon, y = 46.75, 
                                                                     xend = meanLon, yend = 47.6),
               color = "black", size = 0.5, linetype = "solid", arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Adjust the legend for stopoverLengthDays to range from 17 (shortest stopover) to 25 (longest stopover) by 2
  scale_size_continuous(name = "Stopover Length (days)", 
                        breaks = seq(min(as.numeric(plotDataCCLO$stopoverLengthRound)), max(as.numeric(plotDataCCLO$stopoverLengthRound)), by = 2),  # Define the breaks for the legend
                        labels = seq(min(as.numeric(plotDataCCLO$stopoverLengthRound)), max(as.numeric(plotDataCCLO$stopoverLengthRound)), by = 2),  # Define the labels for the legend
                        range = c(2, 10)) +
  scale_fill_manual(values = c("Nonbreeding season range" = "purple4", "61587" = "#440154", "67508" = "#21918c", "61596" = "#fde725")) +
  guides(
    fill = guide_legend(order = 2, title = "Legend", ncol = 1))

### SPPI ####

plotDataSPPI <- valid_detections_df_SPPI |>
  group_by(group_id) |>
  mutate(stopoverLengthDays = difftime(max(ts), min(ts), units = "days"),
         stopoverLengthRound = round(stopoverLengthDays, digits = 0),
         meanLat = mean(recvDeployLat),
         meanLon = mean(recvDeployLon),
         min_ts = as.Date(min(ts)), # adding the start and end points for plotting
         max_ts = as.Date(max(ts))) |>
  filter(recvDeployLat < 49) |>
  mutate(siteType = "Tag site")

xmin <- min(plotDataSPPI$meanLon, na.rm = TRUE) - 10
xmax <- max(plotDataSPPI$meanLon, na.rm = TRUE) + 10
ymin <- min(plotDataSPPI$meanLat, na.rm = TRUE) - 10
ymax <- max(plotDataSPPI$meanLat, na.rm = TRUE) + 10

ggplot(data = world) + 
  geom_sf(colour = "black") +
  #geom_sf(data = filter(allWinterRanges, speciesEN == "Sprague's Pipit"), aes(fill = as.factor(season2)), color = NA, alpha = 0.3) + ignoring this layer since SPPI didn't have any stopovers >24h as far south as the wintering range
  geom_sf(data = lakes, colour = NA, fill = "white") +
  geom_sf(data = states_CanUS, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_bw() + 
  geom_point(data = plotDataSPPI, aes(x = tagDepLon, y = tagDepLat), colour = "black", fill = "red", shape = 23) +
  geom_point(data = filter(plotDataSPPI, stopoverLengthRound > 1), aes(x = meanLon, y = meanLat, size = as.numeric(stopoverLengthRound), fill = as.factor(motusTagID)), shape = 21, alpha = 0.5) +
  labs(x = "", y = "", title = "Sprague's Pipit Stopover Sites", fill = "Tag ID") +
  geom_text(data = states_CanUS, aes(x = longitude, y = latitude, label = postal)) + 
  geom_label(data = plotDataSPPI, aes(x = tagDepLon, y = tagDepLat), label = "Tag Site", colour = "black", size = 3, hjust = 0.75, vjust = -0.5) +
  geom_label(data = plotDataSPPI, aes(x = meanLon, y = meanLat, label = paste("Date range:", min_ts, "-", max_ts)), color = "black", size = 3, vjust = -1) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_size_continuous(name = "Stopover Length (days)", 
                        breaks = c(6),  # Define the breaks for the legend
                        labels = c(6),  # Define the labels for the legend
                        range = c(2, 10)) +
  scale_fill_manual(values = c("61024" = "#440154"))

#### 5 - regression of distance to recv vs. number of detections ####
# exploring whether we can use the tag deployment distance from a receiver to predict whether a tag will be detected
# I limited this to detections within 30 days of tag deployment since deployment distance should have no bearing on whether a tag was detected post-departure from the breeding grounds
# I verified visually that within 30 days, the majority of tags were only detected at local receivers. a few started moving around towards the end of that period but since this is a binary response and they were also detected locally, I don't think it should impact the model

# loading the file of closest receivers that I generated in script 3
closestRecv <- readRDS("closestRecv.2025Mar6.rds")

# loading tag file and filtering for deployed tags only
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

depTags912 <- tags912 |>
  filter(is.na(utcYearStart) == FALSE) |>
  select(tagID, mfgID, speciesEN, period, utcYearStart, utcMonthStart, utcDayStart, utcHourStart, utcMinuteStart, decimalLatitude, decimalLongitude) |>
  mutate(tagDeployDate = make_datetime(
    utcYearStart,
    utcMonthStart,
    utcDayStart
  )) |>
  rename(motusTagID = tagID,
         tagBI = period,
         tagDepLat = decimalLatitude,
         tagDepLon = decimalLongitude) |>
  filter(is.na(tagDepLat) == FALSE)

# combining deployed tags with closest receivers
depTags912 <- depTags912 |>
  left_join(closestRecv, by = "motusTagID")

# adding a column for whether the tag has BREEDING SEASON detections (within 30 days of deployment) in df.alltags.912.clean 

# Calculate the date difference from deployment (in days)
df.alltags.912.clean$daysSinceDep <- as.numeric(difftime(df.alltags.912.clean$ts, df.alltags.912.clean$tagDeployDate, units = "days"))

breedingDets <- df.alltags.912.clean |>
  filter(daysSinceDep <= 30)

depTags912 <- depTags912 |>
  mutate(tagDetected = case_when(
    motusTagID %in% breedingDets$motusTagID ~ "y",
    !motusTagID %in% breedingDets$motusTagID ~ "n",
    TRUE ~ "NA"
  ))

# adding a column for tag model, which is only available in another file that I downloaded from the Motus site
tags <- read_csv("tags.csv")

tagModels <- tags |>
  select(tagID, model) |>
  rename(motusTagID = tagID,
         tagModel = model)

depTags912 <- depTags912 |>
  left_join(tagModels, by = "motusTagID")

# negative binomial regression for whether tag was detected or not as a function of deployment distance from closest receiver
# Ensure the variables are in the right format (factors where necessary)
depTags912 <- depTags912 %>%
  mutate(
    tagDetected = factor(tagDetected, levels = c("n", "y")), # Ensure tagDetected is a factor with "n" and "y" levels
    tagModel = factor(tagModel),
    speciesEN = factor(speciesEN),
    utcYearStart = factor(utcYearStart)
  )

# running negative binomial models separately for speciesEN and tagModel because they are confounding factors
binom_model_speciesEN <- glm(tagDetected ~ distToClosestRecv + speciesEN, 
                            data = depTags912, 
                            family = binomial(link = "logit"))
summary(binom_model_speciesEN)

binom_model_tagModel <- glm(tagDetected ~ distToClosestRecv + tagModel, 
                             data = depTags912, 
                             family = binomial(link = "logit"))
summary(binom_model_tagModel)

AICc(binom_model_speciesEN, binom_model_tagModel) # species model is slightly better apparently

# checking model fit using residuals (DHARMa package)
residuals <- simulateResiduals(binom_model_speciesEN)
plot(residuals)
# DHARMa doesn't seem to find any significant issues with the model. the top line on the right plot is a bit wonky but everything else looks fine

# Check relationship between categorical predictors and tagDetected
table(depTags912$tagModel, depTags912$tagDetected)
table(depTags912$speciesEN, depTags912$tagDetected)

# can I model each species/tag type?
table(depTags912$speciesEN, depTags912$tagModel)
table(depTags912$tagModel, depTags912$tagDetected)
# one issue with running the model with all of the data is that certain tag types and species will be highly correlated with one another. so let's try breaking the data down and running subsets of the global model to get at differences in detectability for tag models and species

depTags912 <- depTags912 |>
  mutate(tagManufacturer = case_when(
    tagModel %in% c("LifeTag", "HybridTag") ~ "CTT",
    tagModel %in% c("NTQB2-5-1-M", "NTQB2-2-M","NTQB2-3-2-M") ~ "Lotek", # SPPI only has 2-5-1s for lotek tags
    TRUE ~ "NA" # to catch anything that might be missed
  ))

# CCLO only model to pull apart tag model effects - CCLO had mostly NTQB2-2-M and NTQB2-5-1-M plus 6 LifeTags. had to remove the lifeTags though because none of the 6 were detected within 30 days of deployment
binom_model_CCLO <- glm(tagDetected ~ distToClosestRecv + tagModel, 
                   data = filter(depTags912, speciesEN == "Chestnut-collared Longspur" & tagModel %in% c("NTQB2-2-M", "NTQB2-5-1-M")), 
                   family = binomial(link = "logit"))
# Warning message:glm.fit: fitted probabilities numerically 0 or 1 occurred. I was able to avoid that by filtering only for distToClosestRecv <10 because nothing over that had any breeding season detections

summary(binom_model_CCLO)

hist(depTags912$distToClosestRecv)

CCLO <- depTags912 |>
  filter(speciesEN == "Chestnut-collared Longspur" & tagModel %in% c("NTQB2-2-M", "NTQB2-5-1-M"))
# to get the model to run, I had to get rid of anything with a distance to closest recv > 10km because above that NOTHING was detected within 30 days of deployment, so the model wouldn't fit properly. now dist is not significant in the model, it was before but with fit issues... seems like 10km is roughly a cutoff for CCLO

# HOLA - going to combine both Lotek tag types into one because there were only 4 NTQB2-3-2-M, and I'll combine CTT tags because there were just two lifeTags

HOLA <- depTags912 |>
  filter(speciesEN == "Horned Lark") |>
  mutate(tagManufacturer = case_when(
    tagModel %in% c("LifeTag", "HybridTag") ~ "CTT",
    tagModel %in% c("NTQB2-3-2-M", "NTQB2-5-1-M") ~ "Lotek",
    TRUE ~ "NA" # to catch anything that might be missed
  ))

binom_model_HOLA <- glm(tagDetected ~ distToClosestRecv + tagManufacturer, 
                        data = filter(depTags912, speciesEN == "Horned Lark"), 
                        family = binomial(link = "logit"))

summary(binom_model_HOLA)
# detection probability does go down with increasing distance, no impact of tag type

resids <- simulateResiduals(binom_model_HOLA)
plot(resids)
# no significant issues it seems

# SPPI - going to combine both CTT tag types into one because there were only 5 + 9 of each
SPPI <- depTags912 |>
  filter(speciesEN == "Sprague's Pipit") |>
  mutate(tagManufacturer = case_when(
    tagModel %in% c("LifeTag", "HybridTag") ~ "CTT",
    tagModel == "NTQB2-5-1-M" ~ "Lotek", # SPPI only has 2-5-1s for lotek tags
    TRUE ~ "NA" # to catch anything that might be missed
  ))

binom_model_SPPI <- glm(tagDetected ~ distToClosestRecv + tagManufacturer, 
                        data = filter(depTags912, speciesEN == "Sprague's Pipit"), 
                        family = binomial(link = "logit"))

summary(binom_model_SPPI)
# detection probability does go down with increasing distance, no significant impact of tag type although lotek tags were slightly more likely to be detected (p = 0.09)

resids <- simulateResiduals(binom_model_SPPI)
plot(resids)
# no significant issues it seems

# species effect - will pull individual tag types that had decent sample sizes across species

binom_model_251 <- glm(tagDetected ~ distToClosestRecv + speciesEN, 
                        data = filter(depTags912, tagModel == "NTQB2-5-1-M"), 
                        family = binomial(link = "logit"))

summary(binom_model_251) # distance is still the key factor

resids <- simulateResiduals(binom_model_251)
plot(resids)
# some issues detected but it doesn't look that bad to me?

binom_model_22 <- glm(tagDetected ~ distToClosestRecv + speciesEN, 
                       data = filter(depTags912, tagModel == "NTQB2-2-M"), 
                       family = binomial(link = "logit"))

summary(binom_model_22) # distance is still the key factor

resids <- simulateResiduals(binom_model_22)
plot(resids)
# no issues

# now I will try to plot all of the distance relationships
# note that I don't think this makes a lot of sense, but I kept it in here in case it's useful at some point
# the plot makes it look like all of the predictions are from the same model, and covers distances that were not in the data...

# plotting the model predicted probabilities

# Create a grid of 'distToClosestRecv' values (100 values between min and max)
dist_values <- seq(min(depTags912$distToClosestRecv), 
                   max(depTags912$distToClosestRecv), 
                   length.out = 100)

# Create a data frame that contains all combinations of distToClosestRecv and the model labels
plot_data <- expand.grid(distToClosestRecv = dist_values,
                         model = c("SPPI", "HOLA", "CCLO", "22", "251"),
                         tagManufacturer = "Lotek", # no tag differences so it doesn't matter which i pick
                         speciesEN = "Chestnut-collared Longspur", # this species was present in both tag models, and there were no species differences 
                         tagModel = "NTQB2-5-1-M") # same thing here

# Add predicted probabilities for each model
plot_data <- plot_data %>%
  mutate(
    predicted_prob = case_when(
      model == "SPPI" ~ predict(binom_model_SPPI, newdata = plot_data, type = "response"),
      model == "HOLA" ~ predict(binom_model_HOLA, newdata = plot_data, type = "response"),
      model == "CCLO" ~ predict(binom_model_CCLO, newdata = plot_data, type = "response"),
      model == "22"   ~ predict(binom_model_22, newdata = plot_data, type = "response"),
      model == "251"  ~ predict(binom_model_251, newdata = plot_data, type = "response")
    )
  )

# Plot the predicted probabilities
ggplot(plot_data, aes(x = distToClosestRecv, y = predicted_prob, color = model)) +
  geom_line() +
  labs(
    title = "Predicted Probability of Tag Detection vs. Distance to Closest Receiver",
    x = "Distance to Closest Receiver (km)",
    y = "Predicted Probability of Detection"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange"))  # Optional: customize colors


#### 6 - detection probability with harness tube or no tube ####

# reading in banding data to get the info on which tags had a tube vs no tube
bandingData <- read_csv("GBMP_and_Suffield_banding_data_Motusonly.csv")

# selecting the columns I need
harnessData <- bandingData |>
  select(Year, Month, Day, Species_Code, Tag_Type, Tag_ID, tube_notube)

# now need to match to motusTagID based on deployment date, species and tagID (there is no motusTagID in the banding data)
harnessData <- harnessData |>
  mutate(tagDeployDate = as.Date(paste(Year, Month, Day, sep = "-")),
         tube = case_when(tube_notube == "no_tube" ~ "n",
                          tube_notube == "tube" ~ "y",
                          TRUE ~ "NA")) |>
  select(-Year, - Month, -Day, -Tag_Type, -tube_notube)

# need to round mfgID because some of them have decimals but the IDs in the harness data do not
df.alltags.912.clean <- df.alltags.912.clean %>%
  mutate(mfgID_round = ifelse(grepl("^[-+]?[0-9]*\\.?[0-9]+$", mfgID), 
                             round(as.numeric(mfgID), 0), 
                             mfgID))
# get a warning about NAs introduced by coercion but it appears to have worked
table(df.alltags.912.clean$mfgID_round)

harnessData <- harnessData |>
  rename(mfgID_round = Tag_ID,
         bandingCode = Species_Code)

# adding to the detections data
df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(tagDeployDate = as.Date(tagDeployDate)) |>
  left_join(harnessData, by = c("mfgID_round", "bandingCode", "tagDeployDate"))

# adding in the number of hits and receivers each tag was detected on
df.alltags.912.clean <- df.alltags.912.clean |>
  group_by(motusTagID) |>
  mutate(nHits = n(),
         nRecv = n_distinct(recvDeployName))

# summary of detected tags, by harness type
detections_by_harness_type <- df.alltags.912.clean |>
  ungroup() |>
  filter(!tube == "NA" & !is.na(tube)) |>
  group_by(tube, speciesEN) |>
  summarise(nTagsDetected = n_distinct(motusTagID),
            nHitsAllTags = n(),
            avgHitsPerTag = mean(nHits),
            avgHitsPerTag_SD = sd(nHits, na.rm = TRUE),
            avgDaysDetected = mean(numDaysDetected),
            avgDaysDetected_SD = sd(numDaysDetected, na.rm = TRUE),
            meanNRecv = mean(nRecv),
            nRecv_sd = sd(nRecv),
            nTagsInNextYear = n_distinct(motusTagID[year == as.numeric(tagDepYear + 1)]),
            .groups = 'drop') |>
  mutate(
    avgHitsPerTag = paste0(round(avgHitsPerTag, 0), " ± ", round(avgHitsPerTag_SD, 0)),  # Combine mean and SD
    avgDaysDetected = paste0(round(avgDaysDetected, 0), " ± ", round(avgDaysDetected_SD, 0)),  # Combine mean and SD
    meanNRecv = paste0(round(meanNRecv, 0), " ± ", round(nRecv_sd, 0))
  ) |>
  select(-avgHitsPerTag_SD, -avgDaysDetected_SD, -nRecv_sd) |>
  arrange(speciesEN, tube)

# just for Lotek 2-5-1 nanotags
detections_by_harness_type_251 <- df.alltags.912.clean |>
  ungroup() |>
  filter(!tube == "NA" & !is.na(tube) & tagModel == "NTQB2-5-1-M") |>
  group_by(tube, speciesEN) |>
  summarise(nTagsDetected = n_distinct(motusTagID),
            nHitsAllTags = n(),
            avgHitsPerTag = mean(nHits),
            avgHitsPerTag_SD = sd(nHits, na.rm = TRUE),
            avgDaysDetected = mean(numDaysDetected),
            avgDaysDetected_SD = sd(numDaysDetected, na.rm = TRUE),
            meanNRecv = mean(nRecv),
            nRecv_sd = sd(nRecv),
            nTagsInNextYear = n_distinct(motusTagID[year == as.numeric(tagDepYear + 1)]),
            .groups = 'drop') |>
  mutate(
    avgHitsPerTag = paste0(round(avgHitsPerTag, 0), " ± ", round(avgHitsPerTag_SD, 0)),  # Combine mean and SD
    avgDaysDetected = paste0(round(avgDaysDetected, 0), " ± ", round(avgDaysDetected_SD, 0)),  # Combine mean and SD
    meanNRecv = paste0(round(meanNRecv, 0), " ± ", round(nRecv_sd, 0))
  ) |>
  select(-avgHitsPerTag_SD, -avgDaysDetected_SD, -nRecv_sd) 

# now, how many of each harness type were deployed
tube_deps <- harnessData |>
  filter(!tube == "NA" & !is.na(tube)) |>
  mutate(speciesEN = case_when(bandingCode == "SPPI" ~ "Sprague's Pipit",
                               bandingCode == "CCLO" ~ "Chestnut-collared Longspur",
                               bandingCode == "BAIS" ~ "Baird's Sparrow",
                               bandingCode == "HOLA" ~ "Horned Lark", 
                               bandingCode == "GRSP" ~ "Grasshopper Sparrow",
                               bandingCode == "TBLO" ~ "Thick-billed Longspur",
                               TRUE ~ "NA")) |>
  group_by(tube, speciesEN) |>
  summarise(nTagsDeployed = n())
  
# combining deployments with detections for each harness type
detections_by_harness_type <- tube_deps |>
  left_join(detections_by_harness_type, by = c("tube", "speciesEN")) |>
  mutate(percentTagsDetected = round((nTagsDetected/nTagsDeployed)*100, digits = 0),
         percentTagsDetectedNextYear = round((nTagsInNextYear/nTagsDeployed)*100, digits = 0)) |>
  select(tube, speciesEN, nTagsDeployed, nTagsDetected, percentTagsDetected, nHitsAllTags, avgHitsPerTag, avgDaysDetected, meanNRecv, nTagsInNextYear, percentTagsDetectedNextYear) |>
  arrange(speciesEN, tube)

# banding data doesn't have tag model in there, so need to use the tag file that I loaded previously that has that info to look just at # 2-5-1 tags deployed with and without harnesses...
deployed_tags <- tags |>
  select(mfgID, dtStart, motusScientificName, model) |>
  filter(is.na(dtStart) == FALSE) |>
  mutate(mfgID_round = ifelse(grepl("^[-+]?[0-9]*\\.?[0-9]+$", mfgID), 
                              round(as.numeric(mfgID), 0), 
                              mfgID),
         tagDeployDate = as.Date(dtStart),
         bandingCode = case_when(
           motusScientificName == "Anthus spragueii" ~ "SPPI",
           motusScientificName == "Calcarius ornatus" ~ "CCLO",
           motusScientificName == "Rhynchophanes mccownii" ~ "TBLO",
           motusScientificName == "Eremophila alpestris" ~ "HOLA",
           motusScientificName == "Centronyx bairdii" ~ "BAIS",
           motusScientificName == "Ammodramus savannarum" ~ "GRSP",
           TRUE ~ NA),
         speciesEN = case_when(
           motusScientificName == "Anthus spragueii" ~ "Sprague's Pipit",
           motusScientificName == "Calcarius ornatus" ~ "Chestnut-collared Longspur",
           motusScientificName == "Rhynchophanes mccownii" ~ "Thick-billed Longspur",
           motusScientificName == "Eremophila alpestris" ~ "Horned Lark",
           motusScientificName == "Centronyx bairdii" ~ "Baird's Sparrow",
           motusScientificName == "Ammodramus savannarum" ~ "Grasshopper Sparrow",
           TRUE ~ NA
         )) |>
  select(-mfgID, -dtStart, -motusScientificName)

tube_deps_251 <- harnessData |>
  left_join(deployed_tags, by = c("bandingCode", "tagDeployDate", "mfgID_round")) |>
  filter(model == "NTQB2-5-1-M") |>
  group_by(tube, speciesEN) |>
  summarise(nTagsDeployed = n())

detections_by_harness_type_251 <- tube_deps_251 |>
  left_join(detections_by_harness_type_251, by = c("tube", "speciesEN")) |>
  mutate(percentTagsDetected = round((nTagsDetected/nTagsDeployed)*100, digits = 0),
         percentTagsDetectedNextYear = round((nTagsInNextYear/nTagsDeployed)*100, digits = 0)) |>
  select(tube, speciesEN, nTagsDeployed, nTagsDetected, percentTagsDetected, nHitsAllTags, avgHitsPerTag, avgDaysDetected, meanNRecv, nTagsInNextYear, percentTagsDetectedNextYear) |>
  arrange(speciesEN, tube)

#### 7 - Time of day ####
# what time of day were birds detected during migration?
# my code for this is based on solar elevation at the lat/long and time of detection

df.alltags.912.clean <- df.alltags.912.clean %>%
  ungroup() |>
  mutate(
    solar_elevation = mapply(function(dt, lat, lon) {
      # Get the solar position for each datetime, latitude, and longitude
      solar_pos <- getSunlightPosition(dt, lat, lon)
      
      # Return the solar elevation (altitude) - THIS IS IN RADIANS, need to convert to degrees below before determining time of day
      solar_pos$altitude
    }, ts, recvDeployLat, recvDeployLon)
  )

# converting radians to degrees - that's what NOAA uses (see link below) and what my existing code for time of day (also below) was based on
df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(solar_elevation_deg = solar_elevation * (180/pi))

# I spot-checked some of these using the NOAA calculator. just NOTE - there are sometimes 1h discrepancies between the local time that I calculated and the time difference that NOAA uses, I think due to daylight savings time. so just make sure that the NOAA UTC offset matches the difference between ts and localDatetime, otherwise you'll sometimes get incorrect values
# https://gml.noaa.gov/grad/solcalc/ 

df.alltags.912.clean <- df.alltags.912.clean |>
  mutate(timeOfDay = case_when(
    (solar_elevation < -18) ~ "night",
    (solar_elevation > -18 & solar_elevation < -12) ~ "astronomical twilight",
    (solar_elevation > -12 & solar_elevation < -6) ~ "nautical twilight",
    (solar_elevation > -6 & solar_elevation <= 0) ~ "civil twilight",
    solar_elevation > 0 ~ "day",
    TRUE ~ "NA"))

table(df.alltags.912.clean$timeOfDay) # all day and civil twilight

migDets <- df.alltags.912.clean |>
  filter(ts > depDate100km & year == as.numeric(tagDepYear))

table(migDets$timeOfDay)
# day and civil twilight about equally

# Count the number of detections (rows) for each species and timeOfDay combination
timeOfDay_counts <- df.alltags.912.clean %>%
  group_by(speciesEN, timeOfDay) %>%
  summarise(detections = n()) %>%
  ungroup()

# Create the paired bar plot
ggplot(timeOfDay_counts, aes(x = timeOfDay, y = detections, fill = timeOfDay, group = timeOfDay)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Species", y = "Number of Detections", title = "Time of Day of Species Detections") +
  scale_fill_manual(values = c("civil twilight" = "skyblue", "day" = "orange")) +  # Customize colors
  theme_classic() +
  facet_wrap(~ speciesEN, scales = "free_y") +
  theme(legend.position = "none")

# 8 - saving files ####

saveRDS(df.alltags.912.clean, "df.alltags.912.clean.depDates.timeOfDay.2025Mar25.rds")

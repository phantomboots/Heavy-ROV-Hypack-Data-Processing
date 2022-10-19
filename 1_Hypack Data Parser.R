#===============================================================================
# Script Name: 1_Hypack Data Parser_Phantom.R
#
# Script Function: This script is designed to unpack the Hypack .LOG files, and 
# to extract various pieces of the data that are contained within the logs. The 
# script locates files with the .LOG extension, reads them in and extract the 
# timestamps, device (sensor) names, and sensor data for positions, depths, 
# heading, altitude, speed and slant range (from camera). The script searches 
# for duplicate timestamp records for data sources that update faster than 1 Hz, 
# and then removes extra records to preserved a common 1 Hz time series for all 
# sensors.
# 
# The script design at this point is to merge all records from the .LOG files in 
# the working directory into one 'master file', and to then use transect start/
# end time to 'trim' the master file time series to the periods of interest. 
# Specifically, data is trimmed to transect start/end times based on divelog. 
#
# The script will search for the preferred data source first (i.e. CTD depth, 
# rather than onboard depth sensor) and will fall back to extracting the 
# secondary source as required, while also writing a data flag.
# 
# Note: All datetimes in dive log must be in yyyy-mm-dd hh:mm:ss format
#
# Script Author: Ben Snow, adapted by Jessica Nephin
# Script Date: Aug 27, 2019, adapted in Jan 2022
# R version 4.2.1


#===============================================================================
# Packages and session options

# Check for the presence of packages shown below, install if needed
packages <- c("lubridate", "purrr", "dplyr", "zoo", "future.apply", "terra")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required packages
lapply(packages, require, character.only = TRUE)

# Set the number of sig figs high, to facilitate inspection of UTM data records
options(digits = 12)

# Set multisession so future_lapply runs in parallel
plan(multisession)


#===============================================================================
# Select options

# Choose the number of minutes for padding start and end transect times
# If onlyTransects is FALSE, padded time can be set to 0
padtime <- 2

# Clip to transect or dive?
# TRUE for on-transect data only
# FALSE for entire dive including descent and ascent 
# FALSE requires Launch and Recovery times in the dive log
onlyTransects <- TRUE

# Should the secondary depth source be converted to meters?
# Multiply by 3.28084
convert_depth <- FALSE


#===============================================================================
# STEP 1 - SET PATHS AND MAKE EXPORT DIRECTORY

# Set working directory
# Use getwd() if you are already in the right directory
# The project folder and needs to be in the working directory

wdir <- getwd()

# Enter Project folder
project_folder <- "Anchorages_2022"

# Directory where Hypack .RAW files are stored
hypack_path <- file.path(wdir, project_folder, "Raw")

# Directory where dive log csv file is stored
divelog_path <- file.path(wdir, project_folder, "Dive_Logs/Dive_Log.csv")

# Create directory for saving .CSV files
save_dir <- file.path(wdir, project_folder, "1.Hypack_Processed_Data")
dir.create(save_dir, recursive = TRUE) # Will warn if already exists



#===============================================================================
# STEP 2 - SET HYPACK RAW FILE DATA STREAM SOURCES

# Device types to extract data from
device_types <- c("POS","EC1","HCP","GYR","TID")

# Set column names for position, depth, heading, tide (altitude) and heave (speed) 
# Must match the names as listed in hardware devices. If a device is not present, 
# write NULL. MAKE SURE DEVICE NAMES MATCH .RAW FILES 
# No roll, pitch or slant range through hypack with BlueROV
ship_heading_pref <- "Hemisphere_GPS"
GPS_pref <- "Hemisphere_GPS"
depth_pref <- "BlueROV_Depth_Heading" 
pos_pref <- "ROV_USBL_200"
phantom_heading_pref <- "BlueROV_Depth_Heading" 
speed_pref <- "Waterlinked_DVL" 
altitude_pref <- "Waterlinked_DVL" 

# Set names for secondary hardware devices, for cases were primary device may 
# be malfunctioning, NULL if there is no secondary
pos_secondary <- NULL
depth_secondary <- NULL



#===============================================================================
# STEP 3 - START LOG FILE

# Sink output to file
rout <- file( file.path(wdir,project_folder,"1.Hypack_Data_Parser.log"), 
              open="wt" )
sink( rout, split = TRUE ) # display output on console and send to log file
sink( rout, type = "message" ) # send warnings to log file
options(warn=1) # print warnings as they occur

# Start the timer
sTime <- Sys.time( )

# Messages
message("Parsing ", project_folder, " project data on ", Sys.Date(), "\n")
message("Padding start and end of transects by ", padtime, " minutes\n")
message("Only export on-transect data? ", onlyTransects, "\n")

# Print device types and names
cat("Device types:\n")
cat(paste0(device_types, collapse = ", "), "\n")
cat("\nDevice names:\n")
cat("ship_heading_pref = '", ship_heading_pref, "'\n", sep="")
cat("GPS_pref = '", GPS_pref, "'\n", sep="")
cat("depth_pref = '", depth_pref, "'\n", sep="")
cat("pos_pref = '", pos_pref, "'\n", sep="")
cat("phantom_heading_pref = '", phantom_heading_pref, "'\n", sep="")
cat("speed_pref = '", speed_pref, "'\n", sep="")
cat("altitude_pref = '", altitude_pref, "'\n", sep="")
cat("pos_secondary = '", pos_secondary, "'\n", sep="")
cat("depth_secondary = '", depth_secondary, "'\n\n", sep="")



#===============================================================================
# STEP 4 - READ IN RAW DATA

# Lists Hypack files, apply function to extract data on list of files, combine
# extracted data into a single dataframe for entire cruise

# Function to extract data from hypack file
extractHypack <- function( hfile ){
  # Packages
  packages <- c("lubridate","dplyr","purrr")
  lapply(packages, require, character.only = TRUE)
  # Read in lines from hypack files
  hlines <- readLines(hfile)
  # Extract device info 
  # Question: Could grep OFF offset info as well, could that be important?
  dev <- hlines[grepl("DEV", hlines)]
  devices <- sub('.*"(.*)".*', "\\1", dev)
  ID <-  sub('.*DEV ', "", dev) %>% strsplit(., " ") %>% sapply(.,"[[",1)
  dev <- data.frame(ID=ID, Device=devices)
  # Extract UTM Zone
  # From the 3rd column of line 5
  longitude <- hlines[5] %>% strsplit(., " ") %>% 
    unlist() %>% .[3] %>% as.numeric()
  # Get UTM zone based on longitude
  zone <- NA
  zone[longitude == -123] <- 10
  zone[longitude == -129] <- 9
  zone[longitude == -135] <- 8
  # Check
  if(is.na(zone)) warning("Longitude does not match UTM zone 8, 9  or 10", 
                          call.= FALSE)
  # Extract start day
  # from the 3rd column of line 9
  day <- hlines[9] %>% strsplit(., " ") %>% 
    unlist() %>% .[3] 
  # Extract data stream
  # Start after EOH line
  startline <- which(grepl("EOH", hlines)) + 1
  ds_lines <- hlines[startline:length(hlines)]
  # Only include lines with device_types set in STEP 2
  ds_lines <- ds_lines[grepl(paste0(device_types, collapse = "|"), ds_lines)]
  # rbind lines together in dataframe, fills in blanks with NA
  ds <- ds_lines %>% strsplit(., " ") %>% map(., ~ c(X=.)) %>% 
    bind_rows(.) %>% as.data.frame()
  # Rename field 1
  names(ds)[1] <- "Device_type"
  # Set all columns except X1 and X2 as numeric
  if( ncol(ds) > 3 ){
    ds[,-c(1:2)] <- apply(ds[,-c(1:2)], 2, as.numeric)
  } else if (ncol(ds) == 3) {
    ds[,3] <- as.numeric(ds[,3])
  } else {
    stop("Too few columns in RAW file ", hfile, call.=F)
  }
  # Add datetime field
  # Date formated as m/d/y in raw files
  ds$Datetime <- ymd_hms(mdy(day) + seconds(ds$X3))
  # Adjust datetime for second day if needed
  if ( any(ds$X3 < ds$X3[1]) ){
    ds$Datetime[ds$X3 < ds$X3[1]] <- ymd_hms((mdy(day)+1) + 
                                               seconds(ds$X3[ds$X3 < ds$X3[1]]))
  }
  # Merge device and data stream tables
  alldat <- merge(dev, ds, by.x = "ID", by.y = "X2")
  # Add UTM zone field
  alldat$Zone <- zone
  # Return
  return(alldat)
}

# List of hypack files
hypack_files <- list.files(pattern = ".RAW", path = hypack_path, full.names = T)

# Apply function across list of input files
all_list <- future_lapply(hypack_files, FUN=extractHypack, future.seed=42)

# Normally data field are columns X4 and X5, but HCP and TID can have column X6
# If not all files from a survey have an HCP device active, add dummy values to 
# column X6, so that the rbind process below can complete successfully.
for(i in 1:length(all_list)){
  if(is.null(all_list[[i]]$X6)){
    all_list[[i]]$X6 <- rep(NA, length(all_list[[i]]$ID))  
    }
}
  
# Bind list of tables into single dataframe
dat <- do.call("rbind", all_list)



#===============================================================================
# STEP 5 - PROCESS SENSOR DATA

# Split sensor data into individual dataframes, process, then merge together, so
# there is a field for each data type (e.g. depth, heading, etc.) in the final
# data frame. Probably don't want to replace beacon with ship at this point,
# maybe later in #3 if needed. Or if there if there are no positions for a transect.
# Double check what was done in original code.


#=============#
#    Depth    #
#=============#
# Message
message("\nExtracting depth")
# device type == 'EC1' device type
# primary device == depth_pref
# secondary device == depth_secondary
depth_data <- dat[dat$Device_type == "EC1" & dat$Device == depth_pref, 
                  c("X4", "Datetime")]
names(depth_data)[1] <- "Depth_m"
depth_data2 <- dat[dat$Device_type == "EC1" & dat$Device == depth_secondary, 
                   c("X4", "Datetime")]
names(depth_data2)[1] <- "Depth2"
# Merge depths together with datetime
depths <- merge(depth_data, depth_data2, by="Datetime", all=T)
# Assign NA to all negative values
depths$Depth_m[ depths$Depth_m < 0 ] <- NA
depths$Depth2[ depths$Depth2 < 0 ] <- NA
# Convert Depth to meters
# Question: What units should depth_secondary be in?  
if( convert_depth ) depths$Depth2 <- depths$Depth2 * 3.28084
# Remove duplicated
depths <- depths[!duplicated(depths$Datetime),]
# Plot, look for tight relationship between depth sensors
if ( any(!is.na(depths$Depth2)) ) {
  plot(depths$Depth_m, depths$Depth2)
  abline(a=0, b=1, col="red")
}
# Use depth_secondary to fill in gaps in depth_pref where possible
depths$Depth <- ifelse( is.na(depths$Depth_m), depths$Depth2, depths$Depth_m )
# Depth source
depths$Depth_Source <- paste0("Hypack_", depth_pref)
depths$Depth_Source[depths$Depth_m == depths$Depth2] <- 
  paste0("Hypack_", depth_secondary)
table(depths$Depth_Source)
# Subset
depths <- depths[c("Datetime","Depth_m", "Depth_Source")]
# Summary
print(summary(depths))


#==============#
#   Position   #
#==============#
# Message
message("\nExtracting the ROV and ship position")
# Ships Position
# device type == 'POS' device type
# primary device == GPS_pref
ship_GPS_data <- dat[dat$Device_type == "POS" & dat$Device == GPS_pref, 
                     c("X4", "X5", "Datetime", "Zone")]
names(ship_GPS_data)[1:2] <- c("Ship_Easting","Ship_Northing")
# Remove duplicated
ship_gps <- ship_GPS_data[!duplicated(ship_GPS_data$Datetime),]

# ROV position
# device type == 'POS' device type
# primary device == pos_pref
# secondary device == pos_secondary
position_data <- dat[dat$Device_type == "POS" & dat$Device == pos_pref, 
                     c("X4", "X5", "Datetime", "Zone")]
names(position_data)[1:2] <- c("Beacon_Easting","Beacon_Northing")
position_data2 <- dat[dat$Device_type == "POS" & dat$Device == pos_secondary, 
                      c("X4", "X5", "Datetime", "Zone")]
names(position_data2)[1:2] <- c("Beacon_Easting2","Beacon_Northing2")
# Sort primary by Datetime
position_data <- position_data[order(position_data$Datetime),]
# Calculate gap time in seconds between primary position records 
position_data$Beacon_Gaps <- as.numeric(
  c( difftime(tail(position_data$Datetime, -1), 
              head(position_data$Datetime, -1)), 0 ))

# Merge primary, secondary and ship positions together using datetime and zone
pos_list <- list(position_data, position_data2, ship_gps)
positions <- Reduce(function(x, y) merge(x, y, by=c("Datetime","Zone"), all=TRUE), 
                    pos_list, accumulate=FALSE)
# Remove duplicated and order
positions <- positions[!duplicated(positions$Datetime),]
positions <- positions[order(positions$Datetime),]
# Plot, look for tight relationship between position sensor data
# Primary v secondary
if( !is.null(pos_secondary) ){
  plot(positions$Beacon_Easting, positions$Beacon_Easting2)
  abline(a=0, b=1, col="red")
}
# Primary v ship
plot(positions$Beacon_Easting, positions$Ship_Easting)
abline(a=0, b=1, col="red")

# Check for NA in first Beacon_Easting (a gap at the start in the primary)
# If there is a gap at the start, measure it in minutes
if( is.na(positions$Beacon_Easting[1]) ){
  # First Beacon_Easting position
  f <- min(which(!is.na(positions$Beacon_Easting)))
  positions$Beacon_Gaps[1] <- as.numeric(difftime(positions$Datetime[f], 
                                                 positions$Datetime[1]))
}
# Expand gap values
positions$Beacon_Gaps <- na.locf(positions$Beacon_Gaps, fromLast = FALSE)
# Warn: filling in gap values with non-primary data sources
if( any(positions$Beacon_Gaps > 60 & is.na(positions$Beacon_Northing)) ){
  warning("Gaps greater than 60 seconds exist in the primary position sensor.\n", 
          "Attempting to fill using secondary position source.")
}
# Replace primary with secondary when primary is NA for more than 60 seconds
positions$Beacon_Northing <-ifelse( positions$Beacon_Gaps > 60 & 
                                      is.na(positions$Beacon_Northing), 
                                     positions$Beacon_Northing2, 
                                     positions$Beacon_Northing )
positions$Beacon_Easting <- ifelse( positions$Beacon_Gaps > 60 & 
                                      is.na(positions$Beacon_Easting), 
                                    positions$Beacon_Easting2, 
                                    positions$Beacon_Easting )
# Add source field
positions$Beacon_Source <- "Primary"
positions$Beacon_Source[is.na(positions$Beacon_Easting)] <- NA
positions$Beacon_Source[!is.na(positions$Beacon_Easting2)] <- "Secondary"
table(positions$Beacon_Source)

# Convert UTM to lat/lon
# Empty lon and lat columns to fill
positions$Beacon_Longitude <- NA
positions$Beacon_Latitude <- NA
positions$Beacon_Longitude2 <- NA
positions$Beacon_Latitude2 <- NA
positions$Ship_Longitude <- NA
positions$Ship_Latitude <- NA
# Loop through UTM zones
for (z in unique(positions$Zone) ){
  # Subset by zone
  tmp <- filter(positions, Zone == z)
  # Convert to vector layer
  vb <- terra::vect(tmp, geom=c("Beacon_Easting", "Beacon_Northing"),
                    crs=paste0("+proj=utm +zone=", z," +datum=WGS84 +units=m"))
  vb2 <- terra::vect(tmp, geom=c("Beacon_Easting2", "Beacon_Northing2"),
                    crs=paste0("+proj=utm +zone=", z," +datum=WGS84 +units=m"))
  vs <- terra::vect(tmp, geom=c("Ship_Easting", "Ship_Northing"),
                    crs=paste0("+proj=utm +zone=", z," +datum=WGS84 +units=m"))
  # Project to lat/lon
  pb <- terra::project(vb, "+proj=longlat +datum=WGS84")
  pb2 <- terra::project(vb2, "+proj=longlat +datum=WGS84")
  ps <- terra::project(vs, "+proj=longlat +datum=WGS84")
  # Replace with lon and lat values for zone == z rows
  positions[positions$Zone == z,"Beacon_Longitude"] <- geom(pb)[, "x"]
  positions[positions$Zone == z,"Beacon_Latitude"] <- geom(pb)[, "y"]
  positions[positions$Zone == z,"Beacon_Longitude2"] <- geom(pb2)[, "x"]
  positions[positions$Zone == z,"Beacon_Latitude2"] <- geom(pb2)[, "y"]  
  positions[positions$Zone == z,"Ship_Longitude"] <- geom(ps)[, "x"]
  positions[positions$Zone == z,"Ship_Latitude"] <- geom(ps)[, "y"]
}
# Replace NaN generate during lat/long conversion with NA values
positions[is.na(positions)] <- NA

#If  "Beacon_Latitude" is missing, try to get get data for the missing row from Beacon2.
#Only replace when the primary Beacon Lat/Long is missing, and there is an available datum for the secondary
#Beacon (Beacon_Latitude2). Same process for longitude.
positions$Beacon_Longitude <- ifelse(is.na(positions$Beacon_Longitude) & !is.na(positions$Beacon_Longitude2),
                                     positions$Beacon_Longitude2, 
                                     positions$Beacon_Longitude)

positions$Beacon_Latitude <- ifelse(is.na(positions$Beacon_Latitude) & !is.na(positions$Beacon_Latitude2),
                                     positions$Beacon_Latitude2, 
                                     positions$Beacon_Latitude)

# Subset to columns of interest
pos <- positions[c("Datetime", "Beacon_Source", "Beacon_Gaps",
                   "Beacon_Longitude","Beacon_Latitude", 
                   "Ship_Longitude","Ship_Latitude")]

# Summary
print(summary(pos))


#=============#
#   Heading   #
#=============#
# Message
message("\nExtracting the ROV and ship heading")
# device type == 'GYR' device type
# primary device == phantom_heading_pref
# secondary device == ship_heading_pref
rov_heading_data <- dat[dat$Device_type == "GYR" & dat$Device == phantom_heading_pref, 
                        c("X4", "Datetime")]
names(rov_heading_data)[1] <- "ROV_heading"
ship_heading_data <- dat[dat$Device_type == "GYR" & dat$Device == ship_heading_pref, 
                         c("X4", "Datetime")]
names(ship_heading_data)[1] <- "Ship_heading"
# Merge headings together with datetime
headings <- merge(rov_heading_data, ship_heading_data, by="Datetime", all=T)
# Assign NA to all negative values
headings$ROV_heading[ headings$ROV_heading < 0 ] <- NA
headings$Ship_heading[ headings$Ship_heading < 0 ] <- NA
# Remove duplicated
headings <- headings[!duplicated(headings$Datetime),]
# Summary
print(summary(headings))


#=============#
#   Altitude  #
#=============#
# Message
message("\nExtracting altitude")
# device type == 'TID' device type
# primary device == altitude_pref
altitude_data <- dat[dat$Device_type == "TID" & dat$Device == altitude_pref, 
                     c("X4", "Datetime")]
names(altitude_data)[1] <- "Altitude_m"
# Assign NA to all negative values
altitude_data$Altitude_m[ altitude_data$Altitude_m < 0 ] <- NA
#altitude_data$Altitude_m[ altitude_data$Altitude_m > 9.99 ] <- NA
# Remove duplicated
altitude <- altitude_data[!duplicated(altitude_data$Datetime),]
# Summary
print(summary(altitude))

#=============#
#    Speed    #
#=============#
# Message
message("\nExtracting ROV speed")
# device type == 'HCP' device type
# primary device == speed_pref
speed_data <- dat[dat$Device_type == "HCP" & dat$Device == speed_pref, 
                  c("X4", "Datetime")]
names(speed_data)[1] <- "Speed_kts"
# Assign NA to all negative values
# todo maybe filter out large values of speed, over 5 knots?
speed_data$Speed_kts[ speed_data$Speed_kts < 0 ] <- NA
#speed_data$Speed_kts[ speed_data$Speed_kts > 5 ] <- NA
# Remove duplicated
speed <- speed_data[!duplicated(speed_data$Datetime),]
# Check for high speeds (> 30 knots), reject these, as they are obviously false.
hist(speed$Speed_kts, breaks=30)
# Summary
print(summary(speed))



#=============#
#   Combine   #
#=============#
# Message
message("\nCombining all sensor data")
# Merge all together based on datetime
df_list <- list(pos, depths, headings, altitude, speed )
sdat <- Reduce(function(x, y) merge(x, y, by="Datetime", all=TRUE), 
               df_list, accumulate=FALSE)

# Sort by Datetime
sdat <- sdat[order(sdat$Datetime),]



#===============================================================================
# STEP 6 - EXTRACT PLANNED LINES FROM RAW HYPACK FILES

# Planned line file coordinates are logged in each Hypack .RAW file, and have 
# the line designations 'PTS' and 'LNN'. The exact line number will vary 
# depending on how many devices are present in the Hypack Hardware file.

# Function to extract planned lines data from hypack file
extractPlannedLines <- function( hfile ){
  # Read in lines from hypack files
  hlines <- readLines(hfile)
  # Extract UTM Zone
  # From the 3rd column of line 5
  longitude <- hlines[5] %>% strsplit(., " ") %>% 
    unlist() %>% .[3] %>% as.numeric()
  # Get UTM zone based on longitude
  zone <- NA
  zone[longitude == -123] <- 10
  zone[longitude == -129] <- 9
  zone[longitude == -135] <- 8
  # Check
  if(is.na(zone)) warning("Longitude does not match UTM zone 8, 9  or 10", 
                          call.= FALSE)
  # Extract planned lines info
  plines <- hlines[grepl("PTS", hlines)]
  PlannedStart <- as.numeric(strsplit(plines[1], " ")[[1]][-1])
  PlannedEnd <- as.numeric(strsplit(plines[2], " ")[[1]][-1])
  lname <- sub("LNN", "", hlines[grepl("LNN", hlines)])
  pldf <- data.frame(lname, t(PlannedStart), t(PlannedEnd), zone)
  names(pldf) <- c("Name", "Start_x", "Start_y", "End_x", "End_y", "Zone")
  # Return
  return(pldf)
}

# Apply function across list of input files
pl_list <- future_lapply(hypack_files, FUN=extractPlannedLines, future.seed=42)
pl <- do.call("rbind", pl_list)

# Remove any duplicates
pl <- pl[!duplicated(pl),]
#Remove any lines that NA values for their names
pl <- pl[!is.na(pl$Name),]

# Convert UTM to lat/lon
# Empty lon and lat columns to fill
pl$Start_Longitude <- NA
pl$Start_Latitude <- NA
pl$End_Longitude <- NA
pl$End_Latitude <- NA
# Loop through UTM zones
for (z in unique(pl$Zone) ){
  # Subset by zone
  tmp <- filter(pl, Zone == z)
  # Convert to vector layer
  sv <- terra::vect(tmp, geom=c("Start_x", "Start_y"),
                    crs=paste0("+proj=utm +zone=", z," +datum=WGS84 +units=m"))
  ev <- terra::vect(tmp, geom=c("End_x", "End_y"),
                    crs=paste0("+proj=utm +zone=", z," +datum=WGS84 +units=m"))
  # Project to lat/lon
  sll <- terra::project(sv, "+proj=longlat +datum=WGS84")
  ell <- terra::project(ev, "+proj=longlat +datum=WGS84")
  # Replace with lon and lat values for zone == z rows
  pl[pl$Zone == z,"Start_Longitude"] <- geom(sll)[, "x"]
  pl[pl$Zone == z,"Start_Latitude"] <- geom(sll)[, "y"]
  pl[pl$Zone == z,"End_Longitude"] <- geom(ell)[, "x"]
  pl[pl$Zone == z,"End_Latitude"] <- geom(ell)[, "y"]
}
# Replace NaN with NA
pl[is.na(pl)] <- NA
# Subset
plannedTransects <- pl[c("Name", "Start_Longitude", "Start_Latitude",
                         "End_Longitude","End_Latitude")]

# Message
message( "\n", nrow(pl), " planned transects found in Hypack for ", 
         project_folder, "\n")


#===============================================================================
# STEP 7 - READ IN DIVE LOG AND MERGE WITH SENSOR DATA

# Read in the start and end time from the Dive Log
dlog <- read.csv(divelog_path)
message( "\n", "Start of dive log... \n")
print(head(dlog))

# Check for all required fields
dnames <- c("Dive_Name", "Transect_Name")
# For only transects data processing
if( onlyTransects ){
  tnames <- c("Start_UTC", "End_UTC")
  cnames <- c(dnames, tnames)
  if( any(!cnames %in% names(dlog)) ) {
    stop( "Missing fields in dive log -> ", 
          paste(cnames[!cnames %in% names(dlog)], collapse = ", "))
  } else {
    # Remove all but required fields
    dlog <- dlog[cnames]
  }
# For launch to recovery data processing (onlyTransects=FALSE)
} else {
  tnames <- c("Start_UTC", "End_UTC","Launch_UTC", "Recovery_UTC")
  cnames <- c(dnames, tnames)
  if( any(!cnames %in% names(dlog)) ) {
    stop( "Missing fields in dive log -> ", 
          paste(cnames[!cnames %in% names(dlog)], collapse = ", "))
  } else {
    # Remove all but required fields
    dlog <- dlog[cnames]
  }
}
# Set datetime format
dlog[tnames] <- lapply( dlog[tnames], FUN=function(x) ymd_hms(x) )
# Pad out the transect start/end times by 'padtime' minutes each to ensure overlap with 
# transect annotations
dlog$Start_UTC_pad <- dlog$Start_UTC - minutes(padtime)
dlog$End_UTC_pad <- dlog$End_UTC + minutes(padtime)
# Crop padded times so there is no overlap between transects
for (i in 2:nrow(dlog)){
 if( dlog$End_UTC_pad[i-1] > dlog$Start_UTC_pad[i] ){
   # Replace end times with non-padded times
   dlog$End_UTC_pad[i-1] <- dlog$End_UTC[i-1]
   # Set start time to previous endtime
   dlog$Start_UTC_pad[i] <- dlog$End_UTC[i-1]
 }
}
# Generate by second sequence of datetimes from start to end of transects 
slog <- NULL
for(i in 1:nrow(dlog)){
  tmp <- data.frame(
    Dive_Name = dlog$Dive_Name[i], 
    Transect_Name = dlog$Transect_Name[i],
    Datetime = seq(dlog$Start_UTC_pad[i],dlog$End_UTC_pad[i], 1))
  # Bind each transect to the previous
  slog <- rbind(slog, tmp)
}
# Add sequence from launch to recovery when onlyTransects=FALSE
if( !onlyTransects ){
  # Dives only
  dives <- dlog[!duplicated(dlog$Dive_Name),]
  lrlog <- NULL
  for(i in 1:nrow(dives)){
    ltmp <- data.frame(
      Dive_Name = dives$Dive_Name[i], 
      Transect_Name = "",
      Datetime = seq(dives$Launch_UTC[i],dives$Recovery_UTC[i], 1) )
    # Bind each transect to the previous
    lrlog <- rbind(lrlog, ltmp)
  }
  # Removes datetimes already in slog
  lrlog <- lrlog[!(lrlog$Datetime %in% slog$Datetime),]
  # Bind to slog
  slog <- rbind(lrlog, slog)
}
# Order
slog <- slog[order(slog$Datetime),]
# Add dive phases
slog$Dive_Phase <- ""
# On transect
for (n in unique(dlog$Transect_Name)){
  p_ind <- which(slog$Datetime >= dlog$Start_UTC_pad[dlog$Transect_Name == n] & 
                 slog$Datetime <= dlog$End_UTC_pad[dlog$Transect_Name == n])
  slog$Dive_Phase[p_ind] <- "Padded_transect"
  t_ind <- which(slog$Datetime >= dlog$Start_UTC[dlog$Transect_Name == n] & 
                 slog$Datetime <= dlog$End_UTC[dlog$Transect_Name == n])
  slog$Dive_Phase[t_ind] <- "On_transect"
}
# Add off_transect phases if onlyTransect=FALSE
if( !onlyTransects ){
  # If empty, add off_transect dive phase
  slog$Dive_Phase[slog$Dive_Phase == ""] <- "Off_transect"
}

# Merge the hypack processed data with the 1 Hz  dive log sequence
dat <- merge(slog, sdat, by = "Datetime", all.x=T)

# Message
cat("\n", length(unique(dat$Dive_Name)), "dives: \n")
cat(paste0(unique(dat$Dive_Name), collapse = "\n"), "\n")
# Check
# Map each dive
for (i in unique(dat$Dive_Name)){
  tmp <- dat[dat$Dive_Name == i,]
  if( !all(is.na(tmp$Ship_Longitude))){
    plot(tmp$Ship_Longitude,tmp$Ship_Latitude, asp=1, main=i, pch=16, cex=.5)
    points(tmp$Beacon_Longitude,tmp$Beacon_Latitude, col="blue", pch=16, cex=.5)
  }
}

# Summary
message("\nSummary of expanded 1Hz hypack data", "\n")
print(summary(dat))


#===============================================================================
# STEP 8 - EXPORT DATA

# Planned transects
# Save as RData for use in later processing scripts
save(plannedTransects, file=file.path(save_dir,"Hypack_PlannedTransects.RData"))

# All hypack sensor data by second
# Save as RData for use in later processing scripts
save(dat, file=file.path(save_dir, "HypackData_bySecond.RData"))
# Save as CSV
write.csv(dat, file = file.path(save_dir, "HypackData_bySecond.csv"), 
          quote = F, row.names = F)



#===============================================================================
# Print end of file message and elapsed time
cat( "\nFinished: ", sep="" )
print( Sys.time( ) - sTime )

# Stop sinking output
sink( type = "message" )
sink( )
closeAllConnections()





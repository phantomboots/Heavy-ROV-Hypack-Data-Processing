#===============================================================================
# Script Name: 3_QAQC_Interpolate__Smooth.R
# Function: This script reads in the .CSV files created from "1_Hypack Data 
#           Parser_****.R', "2_ASDL Data Processing.R' and"2b_Manual Beacon 
#           Position Calculator.R". It fills in the best position and depth data 
#           source for each dive, and interpolates Lat/Longs for Ship GPS 
#           position and vehicle beacon. Next, distance to GPS track is 
#           calculated for each beacon point, and points outside of median+2.5 *
#           IQRare identified as outliers, and beacon fixes at these time points 
#           are removed. Beacon positions are interpolated a second time (with 
#           outliers removed). Beacon positions are then smoothed using a 
#           rolling median (overlapping medians), using stats::runmed(). 
#           Bandwidth is set in the 'STEP 1' portion of these scripts. 
#           Final data for each dive are written to .CSV.
#
# Script Author: Ben Snow, adapted by Jessica Nephin
# Script Date: Sep 9, 2019 adapted in Feb 2022
# R Version: 3.5.1, version 4.0.2


################################################################################
#                                    CHANGE LOG                                #
################################################################################
#
# May 13, 2020: Based on conversations with J.Nephin and S. Jeffery, decided to 
#               add linear interpolation of the depth data stream. Depth record 
#               completion now checks ASDL to see if additional data records can 
#               be found.
# May 24, 2020: Added a check to see if there were additional records that could 
#               be recovered from ASDL for Altitude and Slant Range. FOR 2019 
#               ONLY, added loop to pull Bottom_Velocity_Y (forward velocity). 
#               For some dives in 2019, the 3D dimensional velocity had been 
#               recorded by Hypack, rather then simply forward velocity. Added 
#               linear interpolation of Slant Range and Altitude, but only for 
#               cases where there is just one NA value in a row. This is done by 
#               setting zoo::na.approx(..., maxgap = 1).
# June 3, 2020: Reworked the initial loop that checks ASDL to replace NA values 
#               from the Hypack read-in. This was previously only working 
#               properly for the positional data streams (it checked that GAPS 
#               == T for all data, not just position)
# June 3, 2020: In some cases, interpolation of a -9999 and a normal Altitude or
#               Slant Range values produced values in the range of -4000 or 
#               -5000. These is now a check against this, and values are reset 
#               to -9999 after interpolation.
# July 2, 2020: Removed the portion of this script that was previously used to 
#               generate .KML and .SHP files, and put it in a new script called
#               "4b_Generate_KML_and_Shapefiles.R"
# Aug 20, 2021: Added if() statement to check if the manual beacon tracking 
#               script has been run, before trying to interpolate from this data 
#               source. If this script has not been run, skip this portion.
# Jan 2022: - Wrote function for filling gaps
#           - Applied function to each ASDL source used for filling gaps
#           - Returns number of gaps detected and filled each time
#           - Fixed speed unit error, used to replace knots with m/s
#           - Check the relationship between variables in tofill and forfilling
#             before filling gaps
#           - Question: what to do when relationship is bad, don't fill?
#             eg. speed_kts
#           - attempted to make code more explicit, removed use of column order
#           - rewrote offset section to fix issue with if statements, need to
#             confirm it is working as intended
#           - Using distGeo instead of distm function for step 7
#           - Only removing upper qauntile outliers (not lower)
#           - added log file
#           - saves transect maps to file comparing rov and ship positions
#           - exports all transects and each transect NAV files to csv
#           - removed ship_name as it wasn't used
#           - added original script 4 to this script to reduce # of steps
#           - uses RBR CTD data from excel files first, then ASDL backup second
#           - longlat rounded to 6th decimal place
#           - renamed some attributes
#           - doesn't use ship pos to fill gaps in beacon
#           - check ships pos for outliers (greater than 25m between pos)
#           - outlier cutoff too sensitive, changed to 2.5 IQR
################################################################################



#===============================================================================
# Packages and session options

# Check if necessary packages are present, install as required.
packages <- c("lubridate","dplyr","stringr","imputeTS","geosphere","readxl")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required packages
lapply(packages, require, character.only = TRUE)

# Force display of long digit numbers
options(digits = 12)



#===============================================================================
# STEP 1 - SET PATHS AND SELECT OPTIONS, MAKE EXPORT DIRECTORY

# Set working directory
# Use getwd() if you are already in the right directory
# The project folder and needs to be in the working directory
wdir <- getwd() 

# Project folder
project_folder <- "Feb2022_ROV_Sponge_Coral_MPA"

# Directory where Hypack processed data are stored
hypack_path <- file.path(wdir, project_folder, "Data/1.Hypack_Processed_Data")

# Directory where RBR CTD excel files are stored
RBR_path <- file.path(wdir, project_folder, "Data/RBR_CTD_Data")

# Directory with ASDL master files
ASDL_path <- file.path(wdir, project_folder, "Data/2.ASDL_Processed_Data")

# Export directories 
final_dir <- file.path(wdir, project_folder, "Data/3.Final_Processed_Data")
dir.create(final_dir, recursive = TRUE) # Will warn if already exists
fig_dir <- file.path(wdir, project_folder, "Data/3.Figures")
dir.create(fig_dir, recursive = TRUE) # Will warn if already exists

# # Specify offsets for ship GPS source. If more than one GPS is used, specify 
# # both sources independently. Offset to the port side are positive values for 
# # 'GPS_abeam' and offset towards the bow are positive for 'GPS_along'.
# # Question: are these the distances from the center of the ship, or from the 
# # ROV transponder? 
# GPS_abeam <- -4.1
# GPS_along <- -12.57

# Do you want to use manual trackman ROV position as a backup ROV position?
# TRUE for yes, FALSE for no
trackman <- FALSE

# Set the value to use for the 'window' size (in seconds) of the running median 
# smoothing of the beacon position. Must be odd number!
smooth_window <- 31

# Set the LOESS span values. This is a parameter that described the proportion 
# of the total data set to use when weighting (e.g. span = 0.05 would use 5% 
# of the total data set as the local weighting window).
loess_span = 0.05

# Set the maximum distance (in meters) that can occur between the ROV and ship
# Will differ depending on the ROV
# Question - what is this distance for phantom (100m) and boots (1000m)?
max_dist <- 100


# Fields to include in final processed files
flds <- c("Datetime","Dive_Name", "Transect_Name", "Dive_Phase", 
          "ROV_Longitude_loess", "ROV_Latitude_loess", "ROV_Longitude_smoothed", 
          "ROV_Latitude_smoothed", "ROV_Longitude_unsmoothed", 
          "ROV_Latitude_unsmoothed", "ROV_Source", "Ship_Longitude", 
          "Ship_Latitude", "Depth_m", "Depth_Source", "ROV_heading", 
          "Ship_heading", "Speed_kts", "Altitude_m", "Slant_range_m", 
          "Rogue_pitch", "Rogue_roll", "MiniZeus_zoom_percent",
          "MiniZeus_focus_percent", "MiniZeus_aperture_percent", 
          "MiniZeus_pitch","MiniZeus_roll","ROV_pitch","ROV_roll",
          "Conductivity_mS_cm", "Temperature_C", "Pressure_dbar", 
          "DO_Sat_percent", "Sea_Pressure_dbar", "Salinity_PSU", 
          "Sound_Speed_m_s", "Specific_Cond_uS_cm", "Density_kg_m3", 
          "DO_conc_mgL")


#===============================================================================
# STEP 2 - START LOG FILE

# Sink output to file
rout <- file( file.path(wdir,project_folder,"Data","3.QAQC_Processing.log" ), 
              open="wt" )
sink( rout, split = TRUE ) # display output on console and send to log file
sink( rout, type = "message" ) # send warnings to log file
options(warn=1) # print warnings as they occur

# Start the timer
sTime <- Sys.time( )

# Messages
message("QAQC ", project_folder, " project data on ", Sys.Date(), "\n\n")
# message( "GPS abeam distance = ", GPS_abeam)
# message( "GPS along distance = ", GPS_along)
message( "Smoothing window = ", smooth_window)
message( "loess span = ", loess_span)
message( "Maximum allowable distance between ship and ROV = ", max_dist, "\n\n")


#===============================================================================
# STEP 3 - READ IN THE HYPACK AND ASDL PROCESSED DATA, AND RBR EXCEL FILES

# Read all Master Log files
Slant_Range_Master <- read.csv(file.path(ASDL_path,"Tritech_SlantRange_MasterLog.csv"))
Slant_Range_Master$Datetime <- ymd_hms(Slant_Range_Master$Datetime)
Manual_Track_Master <- read.csv(file.path(ASDL_path,"Manual_Beacon_Tracking_MasterLog.csv"))
Manual_Track_Master$Datetime <- ymd_hms(Manual_Track_Master$Datetime)
Manual_Track_Master$ID <- "Manual tracking backup"
Hemisphere_Master <- read.csv(file.path(ASDL_path,"Hemisphere_GPS_MasterLog.csv"))
Hemisphere_Master$Datetime <- ymd_hms(Hemisphere_Master$Datetime)
Hemisphere_Master$ID <- "Ship GPS backup"
Ship_Heading_Master <- read.csv(file.path(ASDL_path,"Hemisphere_Heading_MasterLog.csv"))
Ship_Heading_Master$Datetime <- ymd_hms(Ship_Heading_Master$Datetime)
DVL_Master <- read.csv(file.path(ASDL_path,"ROWETECH_DVL_MasterLog.csv"))
DVL_Master$Datetime <- ymd_hms(DVL_Master$Datetime)
ROV_Heading_Depth_Master <- read.csv(file.path(ASDL_path,"ROV_Heading_Depth_MasterLog.csv"))
ROV_Heading_Depth_Master$Datetime <- ymd_hms(ROV_Heading_Depth_Master$Datetime)
ROV_Heading_Depth_Master$ID <- "ROV Heading Depth ASDL"
RBR_Master <- read.csv(file.path(ASDL_path,"RBR_CTD_MasterLog.csv"))
RBR_Master$Datetime <- ymd_hms(RBR_Master$Datetime)
RBR_Master$ID <- "RBR CTD ASDL"
MiniZeus_ZFA_Master <- read.csv(file.path(ASDL_path,"MiniZeus_ZFA_MasterLog.csv"))
MiniZeus_ZFA_Master$Datetime <- ymd_hms(MiniZeus_ZFA_Master$Datetime)
ROV_MiniZeus_IMUS_Master <- read.csv(file.path(ASDL_path,"Zeus_ROV_IMU_MasterLog.csv"))
ROV_MiniZeus_IMUS_Master$Datetime <- ymd_hms(ROV_MiniZeus_IMUS_Master$Datetime)

# Read in RBR data from .xlsx files
# Create blank data frame to fill in the loop below.
RBR_Data <- data.frame()
# Read in all RBR .xlsx files in the directory, merge into one larger dataframe
RBR_files <- list.files(path=RBR_path, pattern = ".xlsx", full.names = T)
for(i in RBR_files){
  # High max guess range, sometimes there are many rows missing at beginning
  tmp <- read_xlsx(i, sheet = "Data", skip = 1, guess_max = 10000)
  RBR_Data <- bind_rows(RBR_Data, tmp)
}
# Rename
names(RBR_Data) <- c("Datetime", "Conductivity_mS_cm", "Temperature_C",
                     "Pressure_dbar", "DO_Sat_percent", "Sea_Pressure_dbar",
                     "Depth_m", "Salinity_PSU", "Sound_Speed_m_s",
                     "Specific_Cond_uS_cm", "Density_kg_m3", "DO_conc_mgL")
# Add source ID
RBR_Data$ID <- "RBR CTD XLS"
# Summary 
message("RBR data summary from merged excel files:")
print(summary(RBR_Data))
# Check to see if the RBR data is from year 2000 (default in system)
if ( any(year(RBR_Data$Datetime) == 2000) ){
  warning("RBR year is incorrect!! No RBR CTD data with match up.")
}

# Read in transect data processed at a 1Hz from Hypack (named: dat)
load(file=file.path(hypack_path, "HypackData_bySecond.RData"))

# Is hypack data by transect or by dive?
# Look for 'off-transect' values
if( any(dat$Dive_Phase == 'Off_transect') ) {
  grp <- "Dive_Name"
} else {
  grp <- "Transect_Name"
}

#===============================================================================
# STEP 4 - USE ASDL AND RBR TO FILL IN HYPACK DATA GAPS 

# Function to fill gaps
# Arguments:
# 'tofill' is the dataframe with gaps that need filling
# 'forfilling' is the dataframe to attempt to fill the gaps with
# 'sourcefields' are the names of the columns that you want filled
# 'fillfields' are the names of the columns to fill sourcefields with
# 'type' differentiates between position and other data types
fillgaps <- function(tofill, forfilling, type="", sourcefields, fillfields ){
  # Check if fields are present in the data
  if( any(!sourcefields %in% names(tofill))  ){
    stop("'sourcefields' were not found in 'tofill' dataframe", call.=F)
  }
  if( any(!fillfields %in% names(forfilling)) ){
    stop("'fillfields' were not found in 'forfilling' dataframe", call.=F)
  }
  # Use the first field in sourcefields to find gaps
  field <- sourcefields[1]
  # Find indices with remaining ROV position gaps in 'tofill' DF
  # for position data only fill large gaps, greater than 60 seconds
  if( type == "pos"){
    gapstofill <- which(is.na(tofill[[field]]) & tofill$Beacon_Gaps > 60)
  } else {
    gapstofill <- which(is.na(tofill[[field]]))
  }
  # Find matching indices in 'forfilling' DF by datetime
  fillingrows <- match(tofill$Datetime[gapstofill], forfilling$Datetime)
  # If there are gaps to fill and matches found
  if( length(gapstofill) > 0 && any(!is.na(fillingrows)) ){
    # Fill gaps
    tofill[gapstofill, sourcefields] <- forfilling[fillingrows, fillfields]
  }
  # Remaining gaps
  if( type == "pos"){
    stillgaps <- which(is.na(tofill[[field]]) & tofill$Beacon_Gaps > 60)
  } else {
    stillgaps <- which(is.na(tofill[[field]]))
  }
  # Message
  message( "\n", length(gapstofill), " gaps detected in ", field, "\n",
           length(stillgaps), " gaps remain after filling with ", 
           fillfields[1], "\n")
  # Return
  return(tofill)
}


#===============#
#    POSITION   #
#===============#

# Compare hypack and manual tracking backup beacon positions by dive
tmp <- merge(dat, Manual_Track_Master, by="Datetime")
# Map each dive
for (i in unique(dat$Dive_Name)){
  ot <- dat[dat$Dive_Name == i,]
  mt <- tmp[tmp$Dive_Name == i,]
  plot(ot$Ship_Longitude, ot$Ship_Latitude, 
       asp=1, main=i, pch=16, cex=.5, col="#009E73", 
       xlab = "Longitude", ylab="Latitude")
  points(ot$Beacon_Longitude, ot$Beacon_Latitude, 
         pch=16, cex=.5, col="#0072B2")
  points(mt$Beacon_Longitude.y, mt$Beacon_Latitude.y, 
         pch=16, cex=.5, col="#D55E00")
  legend("bottom", horiz=T, bty = "n", legend = c("Ship", "Hypack", "Manual"),
         col = c("#009E73","#0072B2","#D55E00"), pch=16)
}

# Fill gaps in Beacon long and lat with TrackMan manual ROV GPS
if (trackman ){
  message( "Filling position with TrackMan manual ROV GPS:")
  dat <- fillgaps(tofill=dat,
                  forfilling=Manual_Track_Master,
                  type = "pos",
                  sourcefields=c("Beacon_Longitude", "Beacon_Latitude", 
                                 "Beacon_Source"),
                  fillfields=c("Beacon_Longitude", "Beacon_Latitude", "ID") )
}

#==============#
#    HEADING   #
#==============#
# Check for relationship between tofill and forfilling
tmp <- merge(dat, ROV_Heading_Depth_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$ROV_heading.x, tmp$ROV_heading.y)
# Fill gaps in ROV heading with 'ROV_Heading_Depth_MasterLog.csv'
message( "Filling ROV heading with 'ROV_Heading_Depth_MasterLog.csv':")
dat <- fillgaps(tofill=dat,
                forfilling=ROV_Heading_Depth_Master,
                sourcefields="ROV_heading",
                fillfields="ROV_heading" )
# Check for relationship between tofill and forfilling
tmp <- merge(dat, Ship_Heading_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Ship_heading.x, tmp$Ship_heading.y)
# Fill gaps in ships heading with 'Hemisphere_Heading_MasterLog.csv'
message( "Filling ship heading with 'Hemisphere_Heading_MasterLog.csv':")
dat <- fillgaps(tofill=dat,
                forfilling=Ship_Heading_Master,
                sourcefields="Ship_heading",
                fillfields="Ship_heading" )


#==============#
#     DEPTH    #
#==============#
# Check for relationship between tofill and forfilling
tmp <- merge(dat, RBR_Data, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Depth_m.x, tmp$Depth_m.y)
# Fill gaps in depth with 'ROV_Heading_Depth_MasterLog.csv'
message( "Filling depth with RBR CTD data:")
dat <- fillgaps(tofill=dat,
                forfilling=RBR_Data,
                sourcefields=c("Depth_m","Depth_Source"),
                fillfields=c("Depth_m","ID"))
# Check for relationship between tofill and forfilling
tmp <- merge(dat, RBR_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Depth_m.x, tmp$Depth_m.y)
# Fill gaps in depth with 'ROV_Heading_Depth_MasterLog.csv'
message( "Filling depth with RBR CTD backup:")
dat <- fillgaps(tofill=dat,
                forfilling=RBR_Master,
                sourcefields=c("Depth_m","Depth_Source"),
                fillfields=c("Depth_m","ID"))
# Check for relationship between tofill and forfilling
tmp <- merge(dat, ROV_Heading_Depth_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Depth_m.x, tmp$Depth_m.y)
# Fill remaining gaps in depth with 'ROV_Heading_Depth_MasterLog.csv'
message( "Filling depth with 'ROV_Heading_Depth_MasterLog.csv':")
dat <- fillgaps(tofill=dat,
                forfilling=ROV_Heading_Depth_Master,
                sourcefields=c("Depth_m","Depth_Source"),
                fillfields=c("Depth_m","ID"))

#==================#
#    SLANT RANGE   #
#==================#
# Check for relationship between tofill and forfilling
tmp <- merge(dat, Slant_Range_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Slant_range_m.x, tmp$Slant_range_m.y)
# Fill gaps in slant range with 'Tritech_SlantRange_MasterLog.csv'
message( "Filling slant range with Tritech backup:")
dat <- fillgaps(tofill=dat,
                forfilling=Slant_Range_Master,
                sourcefields="Slant_range_m",
                fillfields="Slant_range_m" )


#===============#
#    ALTITUDE   #
#===============#
# Check for relationship between tofill and forfilling
tmp <- merge(dat, DVL_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Altitude_m.x, tmp$Altitude_m.y)
# Question why leave out of range values as -9999, instead of NA? 
# Fill gaps in altitude with 'ROWETECH_DVL_MasterLog.csv'
message( "Filling altitude with DVL backup:")
dat <- fillgaps(tofill=dat,
                forfilling=DVL_Master,
                sourcefields="Altitude_m",
                fillfields="Altitude_m" )


#==============#
#     SPEED    #
#==============#
# Check for relationship between tofill and forfilling
tmp <- merge(dat, DVL_Master, by="Datetime")
if(nrow(tmp) > 0) plot(tmp$Speed_kts.x, tmp$Speed_kts.y)
# Fill gaps in speed with 'ROWETECH_DVL_MasterLog.csv'
message( "Filling speed with DVL backup:")
dat <- fillgaps(tofill=dat,
                forfilling=DVL_Master,
                sourcefields="Speed_kts",
                fillfields="Speed_kts" )



#===============================================================================
# STEP 5 - ADD ASDL and RBR CTD DATA NOT IN HYPACK DATA 

# Merge MiniZeus fields
dat <- left_join(dat, MiniZeus_ZFA_Master, by = "Datetime")
dat <- left_join(dat, ROV_MiniZeus_IMUS_Master, by = "Datetime")
dat <- left_join(dat, RBR_Data[!names(RBR_Data) %in% c("Depth_m","ID")], 
                 by = "Datetime")

# Use ASDL data to fill RBR CTD data gaps
# Looks for gaps using the first field => Conductivity_mS_cm   
message( "Filling RBR CTD gaps with ASDL RBR master backup:")
dat <- fillgaps(tofill=dat,
                forfilling=RBR_Master,
                sourcefields=c("Conductivity_mS_cm", "Temperature_C",
                               "Pressure_dbar", "DO_Sat_percent", 
                               "Sea_Pressure_dbar", "Salinity_PSU", 
                               "Sound_Speed_m_s", "Specific_Cond_uS_cm"),
                fillfields=c("Conductivity_mS_cm", "Temperature_C",
                             "Pressure_dbar", "DO_Sat_percent", 
                             "Sea_Pressure_dbar", "Salinity_PSU", 
                             "Sound_Speed_m_s", "Specific_Cond_uS_cm"))
# Summary
summary(dat)


#===============================================================================
# STEP 6 - INTERPOLATE TO FILL GAPS

# Interpolate function
# Only interpolates if there are more than 2 non-NA values
interpGaps <- function( x, variable ){
  # Check total not missing values
  total_not_missing <- sum(!is.na(x[[variable]]))
  # check there is sufficient data for na_interpolation 
  if(total_not_missing >= 2) {
    # For altitude and slant range over fill 1 row gaps
    if( grepl("altitude|slant", variable, ignore.case = T) ){
      # Interpolate, max gap = 1
      x[[variable]] <- na_interpolation(x[[variable]], 
                                        option = "linear", maxgap=1)
    } else {
      # Interpolate, max gap = INF
      x[[variable]] <- na_interpolation(x[[variable]], option = "linear")
    }
  } else {
    # Don't interpolate if there aren't enough non-NA values
    x[[variable]] <- x[[variable]]
  }
  # Return
  return(x)
}

# Variables to interpolate
# ROV heading, speed, rogue or minizeus camera variables not interpolated
variables <- c("Beacon_Longitude", "Beacon_Latitude", "Ship_Longitude",
               "Ship_Latitude", "Depth_m", "Ship_heading",
               "Altitude_m", "Slant_range_m")

# Interpolate within each dive or transect by variable
for (v in variables){
  dat <- dat %>% group_by(.dots=grp) %>% 
    group_modify(~interpGaps(.x, variable={{v}})) %>% as.data.frame()
}
# Summary
message("\nSummary after variable gaps were filled with linear interpolation", "\n")
print(summary(dat))

# Set out of bound alitude and slant range to NA
# > 20m for altitude and > 10m for slant range
dat$Altitude_m[dat$Altitude_m > 20] <- NA
dat$Slant_range_m[dat$Slant_range_m > 10] <- NA

# Add "interp" label to beacon source
dat$Beacon_Source[is.na(dat$Beacon_Source)] <- "Interpolation"
# Check
message("\nBeacon position sources:")
print(table(dat$Beacon_Source))



#===============================================================================
# STEP 7 - APPLY OFFSETS TO POSITIONS DATA

# Calculate the bearing and distance from the center of the ship to the GPS 
# antenna. The bearing and distance are used to calculate new lat and long
# coordinates at the GPS antenna location (offset from the center of the ship).
# Offset for the antenna to the port side are (+) for 'GPS_abeam' and offsets 
# towards the bow are (+) for 'GPS_along'. For trig purposes, abeam = opposite 
# and along = adjacent.

# # Compute length of hypotenuse to determine offset distance, in meters.
# # Offset dist will be 0 if GPS_abeam and GPS_along are 0
# offset_dist <- sqrt((GPS_abeam^2) + (GPS_along^2))
# 
# # Calculate the angle from the center of the ship to the antenna location
# # Absolute value
# offset_angle <- atan(GPS_abeam/GPS_along)
# offset_angle <- abs(offset_angle * (180/pi)) #Convert to degrees.
# 
# # Set bearing when the GPS either along=0 or abeam == 0 or both
# # GPS antenna dead center on the ship
# if(GPS_abeam == 0 & GPS_along == 0) {
#   bearing <- 0
# # GPS antenna along keel line, forward of the center of ship
# } else if (GPS_abeam == 0 & GPS_along > 0) {
#   bearing <- dat$Ship_heading
#   # GPS antenna along keel line, aft of the center of ship
# } else if (GPS_abeam == 0 & GPS_along < 0) {
#   bearing <- dat$Ship_heading - 180
#   # GPS antenna centered fore/aft, but port of the keel line
# } else if (GPS_abeam > 0 & GPS_along == 0) {
#   bearing <- dat$Ship_heading - 90
#   # GPS antenna centered fore/aft, but starboard of the keel line
# } else if (GPS_abeam < 0 & GPS_along == 0) {
#   bearing <- dat$Ship_heading + 90
# }
# 
# # Calculate bearing with offset angle when along and abeam are not == 0
# # When GPS is starboard and aft, abeam (-) and along (-)
# if (GPS_abeam < 0 & GPS_along < 0){
#   # Subtract offset angle
#   bearing <- dat$Ship_heading - offset_angle - 180
# # When GPS is port and aft, abeam (+) and along (-)
# } else if (GPS_abeam > 0 & GPS_along < 0){
#   # Add offset angle
#   bearing <- dat$Ship_heading + offset_angle + 180
# # When GPS is starboard and forward, abeam (-) and along (+)
# } else if (GPS_abeam < 0 & GPS_along > 0){
#   # Add offset angle
#   bearing <- dat$Ship_heading + offset_angle
# # When GPS is port and forward, abeam (+) and along (+)
# } else if (GPS_abeam > 0 & GPS_along > 0){
#   # Subtract offset angle
#   bearing <- dat$Ship_heading - offset_angle
# }
# 
# # If bearing is negative, add to 360
# bearing[bearing < 0] <- bearing[bearing < 0] + 360
# # If bearing is greater than 360, subtract 360 from bearing
# bearing[bearing > 360] <- bearing[bearing > 360] - 360
# 
# # Apply offsets from GPS to center of ship using bearing and offset_dist
# # Question: should this offset only be applied to specific beacon sources?
# # e.g. primary but not ship?
# offsets <- destPoint(p=dat[c("Beacon_Longitude",
#                                "Beacon_Latitude")],
#                      b=bearing,
#                      d=offset_dist)
# # Add new offset long/lat to dat
# dat$ROV_Longitude <- offsets[,1]
# dat$ROV_Latitude <- offsets[,2]
dat$ROV_Longitude <- dat$Beacon_Longitude
dat$ROV_Latitude <- dat$Beacon_Latitude


#===============================================================================
# STEP 8 - REMOVE SHIP AND ROV POSITION OUTLIERS

# Order dataframe by datetime
dat <- dat[order(dat$Datetime),]

# Start while loop to remove outliers from ship positions
alongdist_ship <- 100
counter <- 1
while( any(alongdist_ship > 5) ){
  # Calculate distance between adjacent points in ship track
  alongdist_ship <- c(0, geosphere::distGeo(
    as.matrix(dat[1:(nrow(dat)-1), c("Ship_Longitude","Ship_Latitude")]),
    as.matrix(dat[2:nrow(dat), c("Ship_Longitude","Ship_Latitude")])))
  # Set distance to zero if between dives or transects
  for(i in 2:nrow(dat)){
    if( dat$Dive_Name[i-1] != dat$Dive_Name[i] ) alongdist_ship[i] <- 0
    if( dat$Transect_Name[i-1] != dat$Transect_Name[i] ) alongdist_ship[i] <- 0
  }
  # Check
  message("\nWhile loop #", counter)
  message("Summary of distance between adjacent ship positions", "\n")
  print(summary(alongdist_ship))
  # Check
  message("\nRemoved ", length(which(alongdist_ship > 5)) ,
          " outliers greater than 5 m between along track ship positions")
  # Set long/lat values outside of range to NA
  dat$Ship_Longitude[alongdist_ship > 5] <- NA
  dat$Ship_Latitude[alongdist_ship > 5] <- NA
  # Re-interpolate ship long/lat values now that outliers have been removed
  # Variables to interpolate
  variables <- c("Ship_Longitude", "Ship_Latitude")
  # Interpolate within each dive or transect by variable
  for (v in variables){
    dat <- dat %>% group_by(.dots=grp) %>% 
      group_modify(~interpGaps(.x, variable={{v}})) %>% as.data.frame()
  }
  # Number of loops counter
  counter <- counter + 1
} # End of while loop


# Calculate cross track distance between ship pos and beacon positions,
# after along track outliers have been removed and interpolated
crossdist <- geosphere::distGeo(
  as.matrix(dat[c("Ship_Longitude","Ship_Latitude")]), 
  as.matrix(dat[c("ROV_Longitude","ROV_Latitude")]))
# Check
message("\nSummary of distance between ship and ROV position", "\n")
print(summary(crossdist))

# Calculate outlier distance for each dive or transect
# Loop through each dive or transect
for (g in unique(dat[[grp]])){
  # Indices for the group
  ind <- which(dat[[grp]] == g)
  # subset by grp
  distance <- crossdist[ind]
  # Look for outliers (+/- 2.5 * Interquartile range)
  outlier_limit <- median(distance) + (2.5*IQR(distance))
  # Check
  message("\n", grp, " ", g, ": Removed ",
          length(which(distance > max_dist | distance > outlier_limit)),
          " outliers greater than ", round(min(c(outlier_limit, max_dist)),1),
          " m between ship and ROV")
  # Set long/lat values outside of range to NA
  rem <- ind[distance > max_dist | distance > outlier_limit]
  dat$ROV_Longitude[rem] <- NA
  dat$ROV_Latitude[rem] <- NA
}

# Re-interpolate ROV long/lat values now that outliers have been removed
# Variables to interpolate
variables <- c("ROV_Longitude", "ROV_Latitude")
# Interpolate within each dive or transect by variable
for (v in variables){
  dat <- dat %>% group_by(.dots=grp) %>% 
    group_modify(~interpGaps(.x, variable={{v}})) %>% as.data.frame()
}



#===============================================================================
# STEP 9 - APPLY LOESS AND RUNNING MEDIAN SMOOTHING TO ROV POSITION DATA 

# Smoothing function
# Round values to the 6th decimal
smoothCoords <- function( x, variable ){
  # loess smooth
  lname <- paste(variable, "loess", sep="_")
  x[[lname]] <- round(loess(x[[variable]] ~ as.numeric(x[['Datetime']]), 
                            span = loess_span)$fitted, 6)
  # window smooth
  sname <- paste(variable, "smoothed", sep="_")
  x[[sname]] <- round(runmed(x[[variable]], smooth_window), 6)
  # Return
  return(x)
}

# Variables to smooth
variables <- c("ROV_Longitude", "ROV_Latitude")
# Smooth using loess and window method within each dive or transect by variable
for (v in variables){
  dat <- dat %>% group_by(.dots=grp) %>% 
    group_modify(~smoothCoords(.x, variable={{v}})) %>% as.data.frame()
}

# Round unsmoothed to 6 decimals places
dat[variables] <- round(dat[variables], 6)
# Rename unsmoothed long/lat
names(dat)[names(dat) %in% variables] <- c("ROV_Longitude_unsmoothed", 
                                           "ROV_Latitude_unsmoothed")
# Rename beacon source
names(dat)[names(dat) == "Beacon_Source"] <- "ROV_Source"


# Function to add a scalebar to a base-graphics plot
# Adds a scalebar in meters roughly one fifth of plot width
myScalebar <- function(){
  # Get plot coordinates
  pc <- par("usr") 
  # Convert to mercator to get units in meters
  p1 <- geosphere::mercator(pc[c(1,3)])
  p2 <- geosphere::mercator(pc[c(2,4)])
  # Get 1/5th of width of plot as scale bar distance, 2 sig digits
  # Calculate scale bar in meractor
  dist_scale <- signif((p2[1] - p1[1]) * 0.2, digits=2)
  units_label <- paste(dist_scale, "Meters")
  x1 <- p2[1] - (dist_scale * 1.15)
  x2 <- x1 + dist_scale
  y <- p1[2] + ((p2[2] - p1[2]) * 0.05)
  # Convert scale coorinates back to lat/long
  l1 <- geosphere::mercator(c(x1, y), inverse = T)
  l2 <- geosphere::mercator(c(x2, y), inverse = T)
  # Position scale line between last two major x-axis tick marks
  # and 1/10th of the total y-range above the lower y-axis coordinate
  lines(x=c(l1[1],l2[1]), y=c(l1[2],l2[2]), lwd=4, lend=2)
  # Place the units label at the midpoint of and just below the scale line
  text(x=mean(c(l1[1],l2[1])), y=l1[2], label=units_label, adj=c(0.5, -1))
}

# Check
# Map each transect
for (i in unique(dat$Dive_Name)){
  tmp <- dat[dat$Dive_Name == i,]
  png(filename=file.path(fig_dir, paste0("Dive_", i, ".png")),
      width=12, height=10, units = "in", res = 120)
  par(mar = c(5, 4, 4, 8), xpd = TRUE)
  plot(tmp$Ship_Longitude,tmp$Ship_Latitude, 
       asp=1, main=i, pch=16, cex=.5, col="#009E73", 
       xlab = "Longitude", ylab="Latitude")
  points(tmp$Beacon_Longitude,tmp$Beacon_Latitude, 
         pch=16, cex=.5, col="#0072B2")
  points(tmp$ROV_Longitude_unsmoothed, tmp$ROV_Latitude_unsmoothed, 
         pch=16, cex=.5, col="#D55E00")
  points(tmp$ROV_Longitude_smoothed, tmp$ROV_Latitude_smoothed, 
         pch=16, cex=.4, col="grey20")
  points(tmp$ROV_Longitude_loess, tmp$ROV_Latitude_loess, 
         pch=16, cex=.4, col="grey50")
  myScalebar()
  legend("right", inset = c(-0.17, 0), horiz=F, bty = "n",
         legend = c("Ship", "ROV OG", "ROV unsmooth", "ROV smooth", "ROV loess"),
         col = c("#009E73","#0072B2","#D55E00","grey20","grey50"), pch=16)
  dev.off()
}


#===============================================================================
# STEP 10 - WRITE FINAL PROCESSED DATA TO FILE

# Final dataset
fdat <- dat[flds]

# Save as Rdata
save(fdat, file=file.path(final_dir, paste0(project_folder, 
                                            "_SensorData_Georeferenced.RData")))
# Save all transects in one csv
write.csv(fdat, quote = F, row.names = F,
          file = file.path(final_dir, paste0(project_folder, 
                                             "_SensorData_Georeferenced.csv")))

# Export by dive
for (i in unique(fdat$Dive_Name)){
  tmp <- fdat[fdat$Dive_Name == i,]
  write.csv(tmp, quote = F, row.names = F,
            file = file.path(final_dir, 
                             paste0(project_folder, "_", i, 
                                    "_SensorData_Georeferenced_All_Transects.csv")))
}



#===============================================================================
# Print end of file message and elapsed time
cat( "\nFinished: ", sep="" )
print( Sys.time( ) - sTime )

# Stop sinking output
sink( type = "message" )
sink( )
closeAllConnections()



#===============================================================================
# Script Name: 2_ASDL Data Processing_Phantom
#
# Script Function: Reads in the data files created by Advanced Serial Data 
#   logger, parses each file and fills in missing date_time values. Writes one 
#   master file for the whole cruise. Create files with the suffix '_MasterLog', 
#   for use in later scripts Processes records for:
#   - Hemisphere GPS position and heading
#   - TrackMan data 
#   - Waterlinked DVL
#
# Script Author: Ben Snow, adapted by Jessica Nephin
# Script Date: Oct 2022
# R Version: 4.2.1
#
#
# Notes: 
# - depth heading backup come from csv file
# - BlueROV has hemisphere, waterlinked and trackman log files
# - adapted from phantom code


#===============================================================================
# Packages and session options

# Check for the presence of packages shown below, install missing
packages <- c("lubridate","dplyr","stringr","imputeTS","purrr", 
              "future.apply","geosphere")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required packages
lapply(packages, require, character.only = TRUE)

# Force display of long digit numbers
options(digits = 12)

# Set multisession so future_lapply runs in parallel
plan(multisession)


#===============================================================================
# STEP 1 - SET PATHS AND OFFSETS, MAKE EXPORT DIRECTORY

# Set working directory
# Use getwd() if you are already in the right directory
# The project folder and needs to be in the working directory
wdir <- getwd() 

# Enter Project folder name
project_folder <- "Anchorages_2022"

# Directory where the ASDL files are stored
# Path must start from your working directory, check with getwd(), or full paths
ASDL_dir <- file.path(wdir, project_folder, "ASDL")

# Set the directory for saving of the master files.
# Path must start from your working directory, check with getwd(), or full paths
save_dir <- file.path(wdir, project_folder, "2.ASDL_Processed_Data")
dir.create(save_dir, recursive = TRUE) # Will warn if already exists

# Beacon ID number for ROV beacon? (not clump beacon)
ROV_beacon <- 4

# Should the ROV depth source be converted to meters?
# Multiply by 3.28084
convert_depth <- FALSE

#===============================================================================
# STEP 2 - START LOG FILE

# Sink output to file
rout <- file( file.path(wdir,project_folder,"Data","2.ASDL_Data_Parser.log" ), 
              open="wt" )
sink( rout, split = TRUE ) # display output on console and send to log file
sink( rout, type = "message" ) # send warnings to log file
options(warn=1) # print warnings as they occur

# Start the timer
sTime <- Sys.time( )

# Message
message("Parsing ", project_folder, " project data on ", Sys.Date(), "\n")


#===============================================================================
# STEP 3 - READ AND PROCESS ASDL DATA BY TYPE

# Function to read in and process ASDL
# If a file fails to be read or processed, continue to next file
# Return error messages
readASDL <- function( afile, type ){
  tryCatch({
    # Packages
    packages <- c("lubridate","dplyr","stringr","imputeTS","purrr")
    lapply(packages, require, character.only = TRUE)
    # Read in lines from ASDL files
    alines <- readLines(afile, skipNul=FALSE)
    # Remove <COM> messages that can interfere with data columns
    if (type != "track") alines <- sub("<COM.*>","",alines)
    # If 'pos' crop to first GPZDA row, first row with datetime and 
    # subset GPZDA and GPGGA rows only
    if ( type == "pos" | type == "head" ){
      alines <- alines[min(grep("GPZDA",alines)):length(alines)]
      if (type == "pos") alines <- alines[grep("GPZDA|GPGGA",alines)]
      if (type == "head") alines <- alines[grep("GPZDA|HEHDT|GPHDT",alines)]
    } else {
      # Count the number of columns (commas) in each row
      colcount <- str_count(alines,",")
      # Remove lines without the correct number of columns
      # Filters out rows with nulls or errors
      alines <- alines[colcount == median(colcount)]
    }
    # Bind lines together in dataframe, fills in blanks with NA
    ds <-  alines %>% strsplit(., ",")  %>%
      map(., ~ c(X=.)) %>% bind_rows(.) %>% as.data.frame()
    # Make a copy of X1 prior to datetime formatting
    dtype <- ds$X1 # Copy first column
    # Keep only relevant columns based on type
    if( type == "pos" | type == "head" ) ds <- ds[,c("X1","X2","X3","X4","X5")]
    if( type == "track" ) ds <- ds[,paste0("X", c(1:2,4:18))]
    if( type == "dvl" ) ds <- ds[,paste0("X", c(2:4))]
    # Extract datetime in X1
    if( type == "dvl" ){
      names(ds)[1] <- 'X1'
      ds$X1 <- ymd_hms(ds$X1)
    } else if( type == "pos" | type == "head" ){
      date_str <- paste0(ds$X5, ds$X4, ds$X3)
      date_str[dtype == "$GPGGA"] <- "19990909" # place holder to be removed
      time_str <- sub("[.].*","",ds$X2)
      ds$X1 <- suppressWarnings(ymd_hms(paste(date_str, time_str, sep = " ")))
      ds$X1[dtype == "$GPGGA"] <- NA # remove place holder datetimes
    } else {
      date_str <- str_extract(ds$X1, "\\d{8}")
      time_str <- str_extract(ds$X1, "\\d{2}\\:\\d{2}\\.\\d{2}")
      time_str <- gsub("\\.",":",time_str)
      ds$X1 <- suppressWarnings(ymd_hms(paste(date_str, time_str, sep = " ")))
    }
    # Interpolate the datetime series to fill gaps
    ds$X1 <- floor_date(as_datetime(na_interpolation(as.numeric(ds$X1))), "second")
    # For position, reformat lat/lon
    if( type == "pos" ){
      # Subset rows and columns GPS position
      ds <- ds[dtype == "$GPGGA", c("X1","X3","X5")]
      # Reformat to decimal degrees with 7 decimal points
      ds$X3 <- round(as.numeric(substring(ds$X3, 1, 2)) + 
                       (as.numeric(substring(ds$X3, 3, nchar(ds$X3)))/60), 7)
      ds$X5 <- -1 * round(as.numeric(substring(ds$X5, 1, 3)) + 
                            (as.numeric(substring(ds$X5, 4, nchar(ds$X5)))/60), 7)
    } 
    # If heading, subset GPHDT and HEHDT rows only
    if ( type == "head" ){
      ds <- ds[grep("GPHDT|HEHDT",dtype),]
      ds <- ds[,1:2]
    }
    # Remove duplicate time stamps
    ds <- ds[!duplicated(ds$X1),]
    # Set all columns except X1 as numeric
    if( ncol(ds) > 2 ){
      ds[,-1] <- apply(ds[,-1], 2, as.numeric)
    } else {
      ds[,2] <- as.numeric(ds[,2])
    }
    # Remove NA values by checking second row
    ds <- ds[!is.na(ds[,2]),]
    # Return
    return(ds)
  },
  error=function(cond) {
    message(paste("File caused an error:", sub(".*[/]","",afile)))
    message(conditionMessage(cond))
    return(NULL)
  })    
}



#=========================#
#   HEMISPHERE POSITION   #
#=========================#

# Ship GPS backup
# The heading data from this source is not used in subsequently, so was removed 
# from processing to simplify

# List files
GPS_files <- list.files(pattern = "Hemisphere_GPS", path = ASDL_dir, full.names = T)

# Run if files exist
if( length(GPS_files) > 0 ){
  # Message
  message("\nCreating 'Hemisphere_GPS_MasterLog.csv'", "\n")
  # Apply function in parallel
  gpslist <- future_lapply(GPS_files, FUN=readASDL, type="pos")
  # Bind into data frame
  GPS_all <- do.call("rbind", gpslist)
  # Rename
  names(GPS_all) <- c("Datetime","Latitude","Longitude")
  # Summary
  print(summary(GPS_all))
  # Check
  plot(GPS_all$Longitude, GPS_all$Latitude, asp=1)
  # Write
  write.csv(GPS_all, file.path(save_dir,"Hemisphere_GPS_MasterLog.csv"),
            quote = F, row.names = F)
} else {
  stop("\nNo hemisphere position files were found!\n")
}


#========================#
#   HEMISPHERE HEADING   #
#========================#

# Run if files exist
if( length(GPS_files) > 0 ){
  # Message
  message("\nCreating 'Hemisphere_Heading_MasterLog.csv'", "\n")
  # Apply function in parallel
  headlist <- future_lapply(GPS_files, FUN=readASDL, type="head")
  # Bind into data frame
  heading_all <- do.call("rbind", headlist)
  # Rename
  names(heading_all) <- c("Datetime","Ship_heading")
  # Summary
  print(summary(heading_all))
  # Write
  write.csv(heading_all, file.path(save_dir,"Hemisphere_Heading_MasterLog.csv"),
            quote = F, row.names = F)
} else {
  stop("\nNo hemisphere heading files were found!\n")
}

#===============#
#   TRACKMAN   #
#===============#

# List files
track_files <- list.files(pattern = "^TrackMan", path = ASDL_dir, full.names = T)
# Run if files exist
if( length(track_files) > 0 ){
  # Message
  message("\nCreating 'TrackMan_Beacons_MasterLog.csv'", "\n")
  # Apply function in parallel
  tracklist <- future_lapply(track_files, FUN=readASDL, type="track")
  # Bind into data frame
  Track_all <- do.call("rbind", tracklist)
  # Rename
  names(Track_all) <- c("Datetime","Beacon_ID","Phase_Counts_A",
                        "Phase_Counts_B","Phase_Counts_C","Quality_Factor",
                        "Error_Code", "Target_Slant_Range_m","Depression_Angle",
                        "Target_Bearing","DistanceX_m","DistanceY_m","DistanceZ_m",
                        "Ship_Heading","TSS_Pitch","TSS_Roll","Temp_C")
  # Summary
  print(summary(Track_all))
  # Write
  write.csv(Track_all, file.path(save_dir,"TrackMan_Beacons_MasterLog.csv"),
            quote = F, row.names = F)
} else {
  stop("\nNo TrackMan files were found!\n")
}


#======================#
#     Waterlinked DVL  #
#======================#

# List files
DVL_files <- list.files(pattern = "DVL", path = ASDL_dir, full.names = T)
# Run if files exist
if( length(DVL_files) > 0 ){
  # Message
  message("\nCreating 'DVL_MasterLog.csv'", "\n")
  # Apply function in parallel
  dvllist <- future_lapply(DVL_files, FUN=readASDL, type="dvl")
  # Bind into data frame
  DVL_all <- do.call("rbind", dvllist)
  # Set -99.999 values to NA
  DVL_all[,-1] <- apply(DVL_all[,-1], 2, function(x) ifelse(x == -99.999, NA, x))
  # Rename
  names(DVL_all) <- c("Datetime","Altitude_m", "Speed_ms")
  # Convert -1 to NA
  DVL_all$Altitude_m[DVL_all$Altitude_m == -1] <- NA
   # Convert speed to knots
  DVL_all$Speed_kts <- DVL_all$Speed_ms * 1.94384 # m/s to knots
  # todo maybe filter out large values of speed, over 5 knots?
  # Summary
  print(summary(DVL_all))
  # Check for high speeds
  hist(DVL_all$Speed_kts, breaks=30)
  # Write
  write.csv(DVL_all, file.path(save_dir,"DVL_MasterLog.csv"),
            quote = F, row.names = F)
} else {
  stop("\nNo DVL files were found!\n")
}



#===============================================================================
# STEP 4 - CALCULATE ROV POSITION FROM TRACKMAN 

# Computes the coordinates of each beacon using the X, Y distance and bearing 
# from the hydrophone specified in the TrackMan master file and converts these 
# to new Lat/Long positions for the vehicle. No offsets are applied to ships 
# position so calculations are from GPS. Trackman bearing is from center of ship
# so these manually calculated ROV positions could be off if the distance between
# the GPS and center of the ship is large.

# Message
message("\nCreating 'Manual_Beacon_Tracking_MasterLog.csv'", "\n")

# Read in data if not already loaded in workspace
if( !exists("Track_all") ){
  Track_all <- read.csv(file.path(save_dir,"TrackMan_Beacons_MasterLog.csv"))
  Track_all$Datetime <- ymd_hms(Track_all$Datetime)
}
if( !exists("GPS_all") ){
  GPS_all <- read.csv(file.path(save_dir,"Hemisphere_GPS_MasterLog.csv"))
  GPS_all$Datetime <- ymd_hms(GPS_all$Datetime)
}

# Select only ROV beacon
Track_all <- filter(Track_all, Beacon_ID == ROV_beacon)

# Filter out TrackMan readings where no Distance or bad error codes
TrackMan <- filter(Track_all, DistanceX_m != 0 &  DistanceY_m != 0 &
                     Error_Code == 0 | 
                     (Error_Code > 20 & Error_Code != 23 & Error_Code != 73))
# Join the Lat/Long positions to the TrackMan DF
New_Tracking <- left_join(TrackMan, GPS_all, by = "Datetime")
# Calculate distance to ship from DistanceX and DistanceY variables
New_Tracking$Distance <- sqrt(New_Tracking$DistanceX_m^2 + 
                                New_Tracking$DistanceY_m^2)
# Remove records with NA coordinates 
New_Tracking <- New_Tracking[!is.na(New_Tracking$Latitude),]
# Check
hist(New_Tracking$Distance, breaks = 30)

# Generate new points with X,Y distance and bearing from existing Lat/Longs
Beacon_Coords <- destPoint(p=New_Tracking[c("Longitude","Latitude")], 
                           b=New_Tracking$Target_Bearing, 
                           d=New_Tracking$Distance)
# Slot the new Beacon tracking points back into New_Tracking DF
New_Tracking$Beacon_Longitude <- Beacon_Coords[,1]
New_Tracking$Beacon_Latitude <- Beacon_Coords[,2]
# Check
plot(New_Tracking$Beacon_Longitude, New_Tracking$Beacon_Latitude, asp=1)
# Drop unnecessary columns
New_Tracking <- New_Tracking[c("Datetime","Beacon_Longitude",
                               "Beacon_Latitude")]
# Summary
print(summary(New_Tracking))
# Write
write.csv(New_Tracking, quote = F, row.names = F,
          file.path(save_dir,"Manual_Beacon_Tracking_MasterLog.csv"))




#===============================================================================
# Print end of file message and elapsed time
cat( "\nFinished: ", sep="" )
print( Sys.time( ) - sTime )

# Stop sinking output
sink( type = "message" )
closeAllConnections()
  
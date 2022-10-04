#===============================================================================
# Script Name: 4_Generate_KML_and_Shapefiles_Phantom.R
# Script Function: This script reads in the .CSV files created from 
#                 "3_QAQC_Interpolation_Offset_and_Data_Smoothing.R". Creates 
#                 Hypack Planned lines, GPS tracks, and unsmoothed/smoothed 
#                 transponder fixes are written to .KML and .SHP files. 
#                 Shapefiles include a spatial lines file, as well as Spatial 
#                 Points file for each dive. The Spatial Points file contains
#                 all collected data, which is written to the attribute table 
#                 for the file. Spatial lines .SHP are generated from both the 
#                 running median and LOESS smoothing of the Lat/Longs. Spatial 
#                 Points .SHP are generated from the LOESS data only (i.e. 
#                 generally the more accurate of the two smoothing functions).
#
# Script Author: Ben Snow, adapted by Jessica Nephin
# Script Date: Jul 2, 2020 adapted in Feb 2022
# R Version: 3.5.1, version 4.0.2


################################################################################
#                                 CHANGE LOG                                   #
################################################################################
#
# Jan 2022: - Moved processing of hypack raw files for planned transects to #1
#           - Wrote function for creating sp lines
#           - Reduced the number of ind files by exporting all points together
#           - Wrote function for spatial file export
#           - Can export line as dives or transects
################################################################################




#===============================================================================
# Packages 

# Check if necessary packages are present, install as required
packages <- c("rgdal","dplyr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load required packages
lapply(packages, require, character.only = TRUE)



#===============================================================================
# STEP 1 - SET PATHS AND SELECT OPTIONS, MAKE EXPORT DIRECTORY

# Set working directory
# Use getwd() if you are already in the right directory
# The project folder and needs to be in the working directory
wdir <- getwd() 

# Project folder
project_folder <- "Pac2021-054_phantom"

# Name of Ship used in the survey
ship_name <- "Manyberries"

# Directories where final processed data are stored
imports_dir <- file.path(wdir, project_folder, "Data/3.Final_Processed_Data")
hypack_dir <- file.path(wdir, project_folder, "Data/1.Hypack_Processed_Data")

# Export directories 
KML_path <- file.path(wdir, project_folder, "Data/4.Spatial_Data/KML")
dir.create(KML_path, recursive = TRUE) # Will warn if already exists
SHP_path <- file.path(wdir, project_folder, "Data/4.Spatial_Data/Shapefiles")
dir.create(SHP_path, recursive = TRUE) # Will warn if already exists



#===============================================================================
# STEP 2 - LOAD PROCESSED DATA

# Load all sensor data (fdat)
load(file=file.path(imports_dir, paste0(project_folder, 
                                        "_SensorData_Georeferenced.RData")))
# Load planned lines (plannedTransects)
load(file=file.path(hypack_dir, "Hypack_PlannedTransects.RData"))


# Is sensor data by transect or by dive?
# Look for 'off-transect' values
if( any(fdat$Dive_Phase == 'Off_transect') ) {
  type <- "Dive_Name"
} else {
  type <- "Transect_Name"
}



#===============================================================================
# STEP 3 - CREATE SPATIAL LINES OBJECTS

# Create sp lines for planned transects
# Empty list to fill
slines <- list()
for (i in 1:nrow(plannedTransects)) {
  slines[[i]] <- Lines(Line(rbind(
    as.numeric(plannedTransects[i, c("Start_Longitude","Start_Latitude")]), 
    as.numeric(plannedTransects[i, c("End_Longitude", "End_Latitude")]))), 
    ID = as.character(i))
}
# Create spatial lines data frame from list of lines
Planned_Lines <- SpatialLinesDataFrame(
  SpatialLines(slines, proj4string = CRS("+proj=longlat +datum=WGS84")), 
  plannedTransects[,"Name", drop=F], match.ID = FALSE)

# Check
plot(Planned_Lines)


# Create sp lines from actual surveyed positions
# Smoothed and unsmoothed for the ROV and ship

# Create lines function 
createLines <- function(dat, x, y){
  # Empty list to fill
  slines <- list()
  # Create lines for each transect or dive
  for (i in unique(dat[[type]])) {
    tmp <- dat[dat[[type]] == i,]
    slines[[i]] <- Lines(list(Line(tmp[c(x,y)])), ID = as.character(i))
  }
  # Create spatial lines data frame from list of lines
  sldf <- SpatialLinesDataFrame(
    SpatialLines(slines, proj4string = CRS("+proj=longlat +datum=WGS84")), 
    data.frame(Name=unique(dat[[type]])), match.ID = FALSE)  
  # Return
  return(sldf)
}

# Create ship lines
Ship_Lines <- createLines(dat=fdat,
                            x="Ship_Longitude",
                            y="Ship_Latitude")
# Create ROV lines
Unsmoothed_Lines <- createLines(dat=fdat,
                                x="ROV_Longitude_unsmoothed",
                                y="ROV_Latitude_unsmoothed")
# Create smoothed ROV lines
Smooth_Lines <- createLines(dat=fdat,
                            x="ROV_Longitude_smoothed",
                            y="ROV_Latitude_smoothed")
# Create loess ROV lines
Loess_Lines <- createLines(dat=fdat,
                            x="ROV_Longitude_loess",
                            y="ROV_Latitude_loess")
# Check
plot(Ship_Lines)
plot(Unsmoothed_Lines, add=T, col="red")
plot(Smooth_Lines, add=T, col="blue")
plot(Loess_Lines, add=T, col="green")



#===============================================================================
# STEP 4 - CREATE A SPATIAL POINTS DATA FRAME 

# Create a SpatialPointsDataFrame for all transect data. This will contain all 
# sensor data in the attribute table. Uses the LOESS smoothing coordinates
Loess_Points <- SpatialPointsDataFrame(
  coords= fdat[c("ROV_Longitude_loess","ROV_Latitude_loess")],
  data=fdat, proj4string = CRS("+proj=longlat +datum=WGS84"), match.ID=FALSE)

# Shorten names for shapefile
names(Loess_Points) <- sub("Latitude", "y", names(Loess_Points))
names(Loess_Points) <- sub("Longitude", "x", names(Loess_Points))

# Check
plot(Loess_Points)



#===============================================================================
# STEP 5 - EXPORT FILES AS KML AND SHAPEFILES

# Export spatial layers function
exportSpatial <- function( layer, name, driver ){
  if( driver == "KML"){
    writeOGR(layer, file.path(KML_path, paste0(name,".kml")), layer="Name", 
             driver=driver, overwrite_layer = TRUE)
  } else if (driver == "ESRI Shapefile"){
    writeOGR(layer, dsn= SHP_path, layer=name, driver=driver, 
             overwrite_layer = TRUE)
  }
}

# suffix
suf <- sub("_.*", "", type)

# Write lines as KML files
kdriver <- "KML"
exportSpatial(Ship_Lines, paste0(ship_name,"_Track"), kdriver)
exportSpatial(Unsmoothed_Lines, paste0("ROV_Unsmoothed_", suf), kdriver)
exportSpatial(Smooth_Lines, paste0("ROV_Smoothed_", suf), kdriver)
exportSpatial(Loess_Lines, paste0("ROV_Loess_", suf), kdriver)
exportSpatial(Planned_Lines, "Planned_Transects", kdriver)

# Write line as shapefiles
edriver <- "ESRI Shapefile"
exportSpatial(Ship_Lines, paste0(ship_name,"_Track"), edriver)
exportSpatial(Unsmoothed_Lines, paste0("ROV_Unsmoothed_", suf), edriver)
exportSpatial(Smooth_Lines, paste0("ROV_Smoothed_", suf), edriver)
exportSpatial(Loess_Lines, paste0("ROV_Loess_", suf), edriver)
exportSpatial(Planned_Lines, "Planned_Transects", edriver)

# Write points as shapefile
exportSpatial(Loess_Points, "ROV_Loess_Points", edriver)

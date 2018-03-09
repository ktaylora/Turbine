#
# Automated Web-scraping Thingy
# Script will automatically fetch digital obstruction data from the FAA and process
# wind turbines into a WGS84 projected shapefile.  Create a CRON job to run ever 2 months and go crazy.
#
# Author: Kyle Taylor (kyle.taylor@pljv.org)
#

require(raster)
require(rgdal)

# ' scrape a static faa.gov URL for digital obstruction zipfiles and download the most recent zipfile.
# '
# @value returns the filename of the download FAA zipfile
web_scrape_faa <- function(){
  FAA_OBSTRUCTIONS <- "https://www.faa.gov/air_traffic/flight_info/aeronav/digital_products/dof/"
  download.file(FAA_OBSTRUCTIONS,"index.html",quiet=T)
  # parse lines for zipfile path
  lines <- readLines("index.html")
    lines <- lines[grepl(lines,pattern="a href")] # line is a link
      lines <- lines[grepl(lines,pattern="zip")]  # line references a zip
  # parse link to something we can combine to a base URL
  pos <- gregexpr(pattern="\"",text=lines)
  for(i in 1:length(lines)){
    lines[[i]] <- substr(lines[[i]],pos[[i]][1]+1,pos[[i]][2]-1)
  }
  # extract the DOY numeric string from the zipfile name and find the max date
  dates <- as.numeric(na.omit(suppressWarnings(as.numeric(unlist(strsplit(unlist(strsplit(lines,split="[.]")),split="_"))))))
    dates <- dates ==  max(dates)
  # concat the string into something we can fetch -- then download it
  lines <- lines[dates]
    file <- strsplit(lines,split="/")[[1]][2]
  download.file(paste(FAA_OBSTRUCTIONS,lines,sep=""),destfile=file,quiet=T)
  return(file)
}
# ' unpack an FAA zipfile and return as SpatialPointsDataFrame to the user.
# '
# @param x zipfile name containing FAA obstruction data
# @param write optionally write SpatialPointsDataFrame to current working directory.
# @value returns wind turbine data to the user as a SpatialPointsDataFrame
unpack_faa_zip <- function(x,write=T){
  unzip(x)
  coords <- data.frame()
  t <- read.csv("DOF.DAT",t=" ",header=F,skip=4,stringsAsFactors=F)
    t <- t[apply(t,1,function(x) sum(grepl(x,pattern="WINDMILL"))) > 0,] # remove everything that doesn't a WINDMILL field somewhere
  cat(" -- processing ", length(t), " turbine points:", sep="")
  for(i in 1:length(t)){
    ln <- strsplit(t[i],split=" ")[[1]]
    # parse northing and westing, then convert to decimal degrees (WGS84)
    # NORTHING
    northing <- ln[nchar(ln)==6]
      northing <- northing[grepl(northing[nchar(northing)==6],pattern=".*.[.].*N")]
    degrees_n <- as.numeric(ln[which(grepl(ln,pattern=northing))-2])
      minutes_n <- as.numeric(ln[which(grepl(ln,pattern=northing))-1])
        seconds_n <- as.numeric(substr(northing,1,nchar(northing)-1))
    northing <- degrees_n+(minutes_n/60)+(seconds_n/3600)
    # WESTING
    westing <- ln[nchar(ln)==6]
    if(sum(grepl(westing[nchar(westing)==6],pattern=".*.[.].*W"))>0){ # are we using a true westing coordinate?
      westing <- westing[grepl(westing[nchar(westing)==6],pattern=".*.[.].*W")]
      degrees_w <- as.numeric(ln[which(grepl(ln,pattern=westing))-2])
        minutes_w <- as.numeric(ln[which(grepl(ln,pattern=westing))-1])
          seconds_w <- as.numeric(substr(westing,1,nchar(westing)-1))
      westing <- -1*(degrees_w+(minutes_w/60)+(seconds_w/3600))
    } else if(sum(grepl(westing[nchar(westing)==6],pattern=".*.[.].*E"))>0){ # this must be an easting coordinate
      westing <- westing[grepl(westing[nchar(westing)==6],pattern=".*.[.].*E")]
      degrees_w <- as.numeric(ln[which(grepl(ln,pattern=westing))-2])
        minutes_w <- as.numeric(ln[which(grepl(ln,pattern=westing))-1])
          seconds_w <- as.numeric(substr(westing,1,nchar(westing)-1))
      westing <- degrees_w+(minutes_w/60)+(seconds_w/3600)
    } else {
      stop(paste("unknown westing coordinate at:",westing))
    }
    # YEAR
    year <- as.numeric(substr(ln[length(ln)-1],1,4))
    # BIND to master table
    coords <- rbind(coords,data.frame(x=westing,y=northing,year=year))
    if(i%%10==0){ cat(".") }
  }; cat("\n");
  # re-format as a SpatialPointsDataFrame and return to user
  pts <- sp::SpatialPointsDataFrame(coords=data.frame(x=coords$x,y=coords$y),data=data.frame(year=coords$year))
    raster::projection(pts) <- raster::projection("+init=epsg:4326")
  if(write){
    rgdal::writeOGR(pts,".",paste(unlist(strsplit(x,split="[.]"))[[1]],"_pts",sep=""),driver="ESRI Shapefile",overwrite=T)
  }
  # clean-up
  unlink(list.files(pattern=".DAT$|.Dat$"))
  return(pts)
}

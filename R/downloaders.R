# Title     : Web Scrapers and Downloaders
# Objective : Tasks related to scraping, downloading, unpacking, and generating points from FAA obstruction data
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

# ' hidden function that will scrape a static faa.gov URL for digital obstruction zipfiles and
# ' download the most recent zipfile.
# @value returns the filename of the download FAA zipfile
web_scrape_faa_digital_obstructions <- function(FAA_OBSTRUCTIONS="https://www.faa.gov/air_traffic/flight_info/aeronav/digital_products/dof/", write=T){
  download.file(FAA_OBSTRUCTIONS, destfile="index.html", quiet=T)
  # parse lines for zipfile path
  lines <- readLines("index.html")
  lines <- lines[grepl(lines,pattern="a href")] # line is a link
  lines <- lines[grepl(lines,pattern="zip")]    # line references a zip
  # parse link to something we can combine with a base URL
  pos <- gregexpr(pattern="\"",text=lines)
  urls <- sapply(
    X=1:length(lines),
    FUN=function(i){
      substr(lines[[i]], pos[[i]][1]+1, pos[[i]][2]-1)
  })
  # extract the DOY numeric string from the zipfile name and find the max date
  dates <- as.numeric(na.omit(suppressWarnings(
    as.numeric(unlist(strsplit(unlist(strsplit(urls,split="[.]")),
    split="_")))
  )))
  # of all of the FAA obstruction zipfiles found, subset the most recent
  urls <- urls[which.max(dates)]
  # parse a filename from our full download URL
  file <- strsplit(urls,split="/")[[1]]
  file <- file[length(file)] # filename is the last element in the vector
  # clean-up our original index.html file
  file.remove("index.html")
  # should we download the file?
  if(write){
    # download the file
    download.file(urls, destfile=file, quiet=T)
  }
  # return the file path to the user
  return(file)
}
#' check to see if the proposed zip file is the most recent
check_fetch_most_recent_obstruction_file <- function(dir=".", proposed_zip=NULL){
  # do we even have existing zips to compare against? if not, fetch the most recent
  existing_zips <- list.files(dir, pattern="^DOF.*.[.]zip$")
  if(length(existing_zips) == 0) return(TRUE)
  # grab the date strings for our existing zip files
  dates <- as.numeric(na.omit(suppressWarnings(
    as.numeric(unlist(strsplit(unlist(strsplit(existing_zips,split="[.]")),
    split="_")))
  )))
  # parse the date string for our proposed zipfile
  proposed <- as.numeric(na.omit(suppressWarnings(
    as.numeric(unlist(strsplit(unlist(strsplit(proposed_zip,split="[.]")),
    split="_")))
  )))
  # is our proposed zip's date greater-than any existing dates?
  if( sum(proposed <= dates) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# ' hidden function that will unpack an FAA zipfile and return as SpatialPointsDataFrame to the user.
# '
# @param x zipfile name containing FAA obstruction data
# @param write optionally write SpatialPointsDataFrame to current working directory.
# @value returns wind turbine data to the user as a SpatialPointsDataFrame
unpack_faa_zip <- function(x,write=T){
  unzip(x)

  t <- read.csv("DOF.DAT",t=" ",header=F,skip=4,stringsAsFactors=F)
  t <- t[apply(t,1,function(x) sum(grepl(x,pattern="WINDMILL"))) > 0,] # remove everything that doesn't a WINDMILL field somewhere

  cat(" -- processing ", length(t), " turbine points:", sep="")

  coords <- do.call(rbind, lapply(
      X=1:length(t),
      FUN=function(i){
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
        # return
        if(i%%100==0){ cat(".") }
        return(data.frame(x=westing,y=northing,year=year))
      }
  )); cat("\n");
  # re-format as a SpatialPointsDataFrame and return to user
  pts <- sp::SpatialPointsDataFrame(coords=data.frame(x=coords$x,y=coords$y),data=data.frame(year=coords$year))
    raster::projection(pts) <- raster::projection("+init=epsg:4326")
  if(write){
    rgdal::writeOGR(pts,".",paste(unlist(strsplit(x,split="[.]"))[[1]],"_pts",sep=""),driver="ESRI Shapefile",overwrite=T)
  }
  # clean-up
  file.remove(list.files(".", pattern=".DAT$|.Dat$|[.]exe$|[.]EXE$|[.]pdf$"))
  return(pts)
}

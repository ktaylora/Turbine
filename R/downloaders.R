# Title     : Web Scrapers and Downloaders
# Objective : Tasks related to scraping, downloading, unpacking, and generating points from FAA obstruction data
# Created by: Kyle Taylor (kyle.taylor@pljv.org)
# Created on: 4/19/18

# ' hidden function that will scrape a static faa.gov URL for digital obstruction zipfiles and
# ' download the most recent zipfile.
# @value returns the filename of the download FAA zipfile
web_scrape_faa_digital_obstructions <- function(FAA_OBSTRUCTIONS="https://www.faa.gov/air_traffic/flight_info/aeronav/digital_products/dof/"){
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
  # download the file
  download.file(urls, destfile=file, quiet=T)
  # return the file path to the user
  return(file)
}

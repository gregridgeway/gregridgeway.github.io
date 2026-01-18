##############################################################################
##############################################################################
# This script is not essential. It just gives the background on how I
# assembled the individual datasources
##############################################################################
##############################################################################


setwd("z:/presentations/ASW2017/")

library(RSocrata)
library(lubridate)

# load Connecticut traffic stop data
#    details at https://data.ct.gov/Public-Safety/Racial-Profiling-Prohibition-Project-Traffic-Stop-/baka-5j97
#    takes several minutes
dataCT <- read.socrata("https://data.ct.gov/resource/uzvs-9qkn.json")
dataCT$day_of_week <- NULL
dataCT$identification_category_description_text <- NULL
dataCT$intervention_location_description_text <- NULL
dataCT$intervention_identification_id <- NULL
dataCT$organization_activity_text <- NULL
dataCT$statutatory_citation <- NULL
dataCT$statute_code_identification_id <- NULL
dataCT$towed_indicator <- NULL
dataCT$town_resident_indicator <- NULL
dataCT$resident_indicator <- NULL
dataCT$organization_identification_id <- NULL
dataCT$reporting_office_ridentification_id <- NULL

# remove columns we don't need
save(dataCT,file="dataCT.RData",compress=TRUE)


# get darkness data
# 2013
a <- scan("http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2013&task=2&state=CT&place=middletown",
          what="",sep="\n")
a <- a[29:59]
civt.tab <- read.fwf(textConnection(a),
                     widths=c(2,rep(c(6,5),12)),
                     colClasses = "character")
civt <- cbind(expand.grid(day=1:31,month=1:12,year=2013,
                          stringsAsFactors=FALSE),
              end=unlist(civt.tab[,2*(1:12)+1]),
              stringsAsFactors=FALSE)
# 2014
a <- scan("http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2014&task=2&state=CT&place=middletown",
          what="",sep="\n")
a <- a[29:59]
civt.tab <- read.fwf(textConnection(a),
                     widths=c(2,rep(c(6,5),12)),
                     colClasses = "character")
a <- cbind(expand.grid(day=1:31,month=1:12,year=2014,
                       stringsAsFactors=FALSE),
           end=unlist(civt.tab[,2*(1:12)+1]),
           stringsAsFactors=FALSE)
civt <- rbind(civt, a)

# 2015
a <- scan("http://aa.usno.navy.mil/cgi-bin/aa_rstablew.pl?ID=AA&year=2015&task=2&state=CT&place=middletown",
          what="",sep="\n")
a <- a[29:59]
civt.tab <- read.fwf(textConnection(a),
                     widths=c(2,rep(c(6,5),12)),
                     colClasses = "character")
a <- cbind(expand.grid(day=1:31,month=1:12,year=2015,
                       stringsAsFactors=FALSE),
           end=unlist(civt.tab[,2*(1:12)+1]),
           stringsAsFactors=FALSE)
civt <- rbind(civt, a)

civt$end   <- gsub(" ","",civt$end)
civt <- subset(civt, end!="")
civt$end <- gsub("([0-9][0-9])([0-9][0-9])","\\1:\\2",civt$end)
civt$date <- with(civt, ymd(paste0(year,"/",month,"/",day)))
civt$datetime <- with(civt, ymd_hm(paste0(year,"/",month,"/",day," ",end)))
save(civt,file="civt.RData")


###################################
# NYPD data

dataNYPD <- read.csv("2011.csv.gz")
dataNYPD$addrnum <- NULL
dataNYPD$stname  <- NULL
dataNYPD$stinter <- NULL
dataNYPD$crossst <- NULL
dataNYPD$aptnum  <- NULL
dataNYPD$city    <- NULL
dataNYPD$state   <- NULL
dataNYPD$zip     <- NULL

save(dataNYPD,file="dataNYPD.RData",compress=TRUE)

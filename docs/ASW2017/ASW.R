#####################
# Set to your working directory where you downloaded the datasets
#   Or in RStudio click Session->Set Working Directory->Choose Directory
#####################
setwd("z:/presentations/ASW2017/")

# for working with dates and times
library(lubridate)
library(splines)

# just one time
library(devtools)
install_github("gbm-developers/gbm")
install_github("gregridgeway/fastDR")
# for getting propensity score/DR estimates
library(fastDR)

# for making contour plots
library(ks)

##############################################################################
##############################################################################
##  Veil of Darkness exercise
##############################################################################
##############################################################################

# Full details about creating the dataset are in dataprep.R
# load the Connecticut dataset
load("dataCT.RData")
# load the table with end of civil twilight data
load("civt.RData")

# adjust for DST
i <- with(civt, which((date>="2013-03-10" & date<"2013-11-03") |
                      (date>="2014-03-09" & date<"2014-11-02") |
                      (date>="2015-03-08" & date<"2015-11-01")))
civt$datetime[i] <- civt$datetime[i] + hours(1)

# convert to date and time objects
dataCT$date <- with(dataCT, ymd_hm(paste(intervention_date,intervention_time)))
range(dataCT$date)

# subset stops to be within 1 month of DST
i <- with(dataCT, which(abs((date%--%ymd("2013-11-03"))/ddays(1)) < 30 |
                        abs((date%--%ymd("2014-03-09"))/ddays(1)) < 30 |
                        abs((date%--%ymd("2014-11-02"))/ddays(1)) < 30 |
                        abs((date%--%ymd("2015-03-08"))/ddays(1)) < 30))
dataCT <- dataCT[i,]

# was it dark?
i <- match(date(dataCT$date), civt$date)
dataCT$dark <- as.numeric(dataCT$date > civt$datetime[i])

# focus on times near darkness
with(subset(dataCT,dark==1), min(hour(date)+minute(date)/60))
with(subset(dataCT,dark==0), max(hour(date)+minute(date)/60))

dataCT <- subset(dataCT, (hour(date) >= 17) & (hour(date) < 20))

# create minority/non--white outcome
dataCT$minority <- with(dataCT, as.numeric(subjec_trace_code!="W" |
                                           subjec_tethnicit_ycode!="N"))
table(dataCT$minority)

# create a clocktime covariate
dataCT$t <- with(dataCT, hour(date)+minute(date)/60)
hist(dataCT$t)

plot(t~decimal_date(dataCT$date),
     data=dataCT,col=2-dataCT$dark,pch=".")


# estimate VoD
glmVOD <- glm(minority~dark+ns(t,df=15),
              family=binomial,
              data=dataCT)
summary(glmVOD)

a <- coef(summary(glmVOD))["dark",]
exp(-a["Estimate"] + c(0,-2,2)*a["Std. Error"])


##############################################################################
##############################################################################
##  Propensity score/doubly robust with stop data
##############################################################################
##############################################################################

load("dataNYPD.RData")

# compare black and white pedestrians
dataNYPD <- subset(dataNYPD, race %in% c("B","W"))

# create outcomes
dataNYPD$frisk  <- as.numeric(dataNYPD$frisked=="Y")
dataNYPD$search <- as.numeric(dataNYPD$searched=="Y")
dataNYPD$force <- with(dataNYPD,
                       as.numeric(pf_hands=="Y" |
                                  pf_wall =="Y" |
                                  pf_grnd =="Y" |
                                  pf_drwep=="Y" |
                                  pf_ptwep=="Y" |
                                  pf_baton=="Y" |
                                  pf_hcuff=="Y" |
                                  pf_pepsp=="Y" |
                                  pf_other=="Y"))

# define "treatment" indicator
dataNYPD$black <- as.numeric(dataNYPD$race=="B")

# rescale so that X,Y are not so large (1000000+)
dataNYPD$x <- scale(dataNYPD$xcoord)
dataNYPD$y <- scale(dataNYPD$ycoord)
plot(y~x, data=dataNYPD, pch=".")

dataNYPD$precinct <- as.factor(dataNYPD$pct)

# convert datestop and timestop to properly formatted dates and times
a <- gsub("(.*)([0-9]{2})([0-9]{4})$","\\3-\\1-\\2",
          as.character(dataNYPD$datestop))
b <- sprintf("%04d",dataNYPD$timestop)
b <- gsub("([0-9]{2})([0-9]{2})","\\1:\\2:00",b)
dataNYPD$date <- ymd_hms(paste(a,b))

# compute some useful stop features
dataNYPD$month <- month(dataNYPD$date,label=TRUE)
dataNYPD$wday  <- wday(dataNYPD$date,label=TRUE)
dataNYPD$hour  <- hour(dataNYPD$date)+minute(dataNYPD$date)/60
hist(dataNYPD$hour)

# code major categories of crime suspected
dataNYPD$crim <- as.character(dataNYPD$crimsusp)
rev(sort(table(dataNYPD$crim)))[1:50]
dataNYPD$crim[dataNYPD$crim %in% c("F","FEL")] <- "FELONY"
dataNYPD$crim[dataNYPD$crim %in% c("M","MIS","MISD.","CRIM MIS","MISDEMEANOR")] <- "MISD"
dataNYPD$crim[grep("ROB", dataNYPD$crim)] <- "ROBBERY"
dataNYPD$crim[grep("BURG",dataNYPD$crim)] <- "BURGLARY"
dataNYPD$crim[dataNYPD$crim %in% c("GRAND LARC","G/L","FEL/G/L","G.L.","G.LARCENY","GL","GRAN LARC","GR LARC")] <- "GRAND LARCENY"
dataNYPD$crim[dataNYPD$crim %in% c("P/L","PETIT LARC","PL","PET LARC")] <- "PETIT LARCENY"
dataNYPD$crim[dataNYPD$crim %in% c("GLA","FEL/GLA","FELONY/GLA","GRAND LARCENY AUTO","GLA (FEL)","GLA/FEL","G.L.A.","FELONY-GLA")] <- "AUTO THEFT"
dataNYPD$crim[dataNYPD$crim %in% c("CRIM TRES","CRIM TRESPASS","MIS/CRIM TRES","CRIM TRESS","CRIM TRES (MISD)","CRIMINAL TRESSPASS","CRIM-TRES","CRIMINAL TRES","CRIMINAL TRESS","CRIM TRESP","CRIMINAL TRESPASS")] <- "CRIMINAL TRESPASS"
dataNYPD$crim[dataNYPD$crim %in% c("CPM","CRIMINAL POSSESSION OF MARIHUA","MIS/CPM","MISD /CPM","CPM (MISD)","MISD/CPM","MISD/ CPM","CPM/MISD")] <- "MARIJUANA"
dataNYPD$crim[dataNYPD$crim %in% c("CRIM MISC","CRIMINAL MIS","CRIM MISCHIEF")] <- "CRIMINAL MISCHIEF"
dataNYPD$crim[dataNYPD$crim %in% c("PL 490.10","CT")] <- "TERRORISM"
dataNYPD$crim[grep("TOS",dataNYPD$crim)] <- "THEFT OF SERVICE"
dataNYPD$crim[grep("GRAFFITI",dataNYPD$crim)] <- "GRAFFITI"
dataNYPD$crim[grep("GRAFFITTI",dataNYPD$crim)] <- "GRAFFITI"
dataNYPD$crim[grep("GRAFITTI",dataNYPD$crim)] <- "GRAFFITI"
dataNYPD$crim[grep("CPW",dataNYPD$crim)] <- "CPW"
dataNYPD$crim[grep("ASS",dataNYPD$crim)] <- "ASSAULT"
dataNYPD$crim[dataNYPD$crim %in% c("C.P.W.","C.P.W","CRIMINAL POSSESSION WEAPON")] <- "CPW"

a <- table(dataNYPD$crim)
# crime types with at least 200 cases get most of the stops
sum(a[a>200])/sum(a)

dataNYPD <- subset(dataNYPD, crim %in% names(a[a>200]))
dataNYPD$crim <- factor(dataNYPD$crim)

# generate unique case identifier
dataNYPD$key <- 1:nrow(dataNYPD)

# just to make it run a little faster for the demo
# i <- sample(which(dataNYPD$black==1),size=339758-100000,replace=FALSE)
# dataNYPD <- dataNYPD[-i,]

####################
# fit propensity score/DR estimates
# return to slides while this is running
####################
dr1 <- fastDR(list(y.form=~frisk+search+force,
                   t.form=~black,
                   x.form=~x+y+precinct+    # location
                           inout+           # inside or outside
                           trhsloc+         # transit or housing
                           month+wday+hour+ # when
                           radio+           # radio run/call for service
                           crim,            # crime suspected
                   key.form=~key),
              data=dataNYPD,
              # all outcomes are 0/1
              y.dist=c("quasibinomial","quasibinomial","quasibinomial"),
              interaction.depth=3,
              # smaller shrinkage always better, but requires more trees
              shrinkage=0.5,
              # pick n.trees large, but still computable
              n.trees=1000,    # number of boosting iterations
              verbose=TRUE)

# number of treated cases
dr1$n1
# effective sample size of control cases
dr1$ESS

# attach the weights to the dataset
dataNYPD$w <- 0
i <- match(names(dr1$w),dataNYPD$key)
dataNYPD$w[i] <- dr1$w

####################
# check balance
####################
bal <- dr1$balance.tab

i <- grep("wday",rownames(bal))
barplot(100*t(bal[i,c("treatment","control")]),
        beside=TRUE,
        ylim=c(0,20))
i <- grep("month",rownames(bal))
barplot(100*t(bal[i,c("treatment","control")]),
        beside=TRUE,
        ylim=c(0,10))

kdeHour1 <- with(subset(dataNYPD,black==1), kde(hour,h=0.3))
plot(kdeHour1)
kdeHour0 <- with(subset(dataNYPD,black==0), kde(hour,h=0.4,w=w))
plot(kdeHour0,add=TRUE,col="red")

i <- grep("precinct",rownames(bal))
barplot(100*t(bal[i,c("treatment","control")]),
        beside=TRUE,
        ylim=c(0,7))

# gbm also matches on interactions
plot(jitter(y,128)~jitter(x,128),
     data=subset(dataNYPD,black==1 & !is.na(x) & !is.na(y) & precinct==75),
     pch=".")
H <- with(subset(dataNYPD,black==1 & !is.na(x) & !is.na(y) & precinct==75),
             kde(cbind(x,y),w=w))$H
kde1 <- with(subset(dataNYPD,black==1 & !is.na(x) & !is.na(y) & precinct==75),
             kde(cbind(x,y),w=w,H=4*H))
plot(kde1,add=TRUE)
kde0 <- with(subset(dataNYPD,black==0 & !is.na(x) & !is.na(y) & precinct==75),
             kde(cbind(x,y),w=w))
plot(kde0,add=TRUE,col="red")


kdeHour1 <- with(subset(dataNYPD,black==1 & wday=="Fri" & precinct==75),
                 kde(hour))
plot(kdeHour1)
kdeHour0 <- with(subset(dataNYPD,black==0 & wday=="Fri" & precinct==75),
                 kde(hour,h=1,w=w))
plot(kdeHour0,add=TRUE,col="red")

####################
# check outcomes
####################
dr1
dr1$effects

par(mfrow=c(3,1))
with(dr1$effects,
     barplot(100*c(frisk["un","E.y1"],frisk[,"E.y0"]),
             ylab="Frisk rate",
             names.arg=c("Black","White (unadjusted)",
                         "White (PS adjusted)",
                         "White (DR adjusted)")))
with(dr1$effects,
     barplot(100*c(search["un","E.y1"],search[,"E.y0"]),
             ylab="Search rate",
             names.arg=c("Black","White (unadjusted)",
                         "White (PS adjusted)",
                         "White (DR adjusted)")))
with(dr1$effects,
     barplot(100*c(force["un","E.y1"],force[,"E.y0"]),
             ylab="Force rate",
             names.arg=c("Black","White (unadjusted)",
                         "White (PS adjusted)",
                         "White (DR adjusted)")))

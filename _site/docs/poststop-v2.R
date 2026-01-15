
# poststop.R
# Greg Ridgeway
# RAND, Santa Monica, CA
# gregr@rand.org
#
# Description: R script for implementing the propensity score analysis for 
#    post-stop outcomes. An example dataset is included with this script to
#    demonstrate the process. Utilizing your own data will take a little bit
#    of work to read in the data, select the data that interests you, select
#    the important variables and outcomes of interest.
#    For information on obtaining and using the R statistical environment visit
#    www.r-project.org. In particular review 
#    "Introduction to R" at http://cran.r-project.org/doc/manuals/R-intro.pdf
#
# For a detailed discussion of this methodology for studies of racial bias see
#   G. Ridgeway (2006). "Assessing the effect of race bias in 
#   post-traffic stop outcomes using propensity scores." Journal of 
#   Quantitative Criminology 22(1). 
#
# History:
#    07/25/2005 Created
#    03/07/2006 Modified to use the twang package



###############################################################################
# configure some needed packages and helper functions
###############################################################################


# check whether necessary packages are installed
needed.packages <- c("foreign","chron","gbm","survey","twang","xtable")
i <- is.element(needed.packages,
                installed.packages()[,1])
if(!all(i))
{
    a <- paste(needed.packages[!i],collapse=" and ")
    
    cat("Need to install: ",a,". Only needs to be done once. Installing now...\n")
    install.packages(needed.packages[!i])
}


# load the necessary packages
library(foreign)  # methods for importing data from dbf, stata, spss, sas, etc.
library(chron)    # methods for dealing with times and dates
library(survey)   # methods for correct analyses of weighted data
library(twang)    # methods for fitting and evaluating gbm propensity score models
library(xtable)   # for making prettier tables


###############################################################################
# setup the data
###############################################################################

# set the working directory
#    change this to wherever the data file is and where the output will go
#    Note that R uses "/" in file names
setwd("d:/poststop")

# read in data
# read in the example data already prepared in R format
#    This dataset is simulated and does not represent any data from a real police
#    department. I created this dataset so that there is some relationship between
#    race and citation and between race and the other variables.
load("example.RData")
# ...or read in from a dbf file that is easily created in Excel
data <- read.dbf("example.dbf")

# subset the data at this point if needed, for example
#    only black vs. white
#       data <- subset(data, (race=="W") | (race=="B"))
#    drop certain neighborhoods
#       data <- subset(data, !(nhood %in% c(1,3,5,7)))


# make sure that R knows which variables are categorical
# R should do this by default but categorical variables that are coded as integers
#    (like neighborhood codes) will get interpreted as numeric values. factor()
#    will correct this, and can add better labels (e.g. see reason)
data$nhood  <- factor(data$nhood)
data$age    <- factor(data$age)
data$race   <- factor(data$race)
data$reason <- factor(data$reason,
                      labels=c("mech/reg","MV danger","MV nodanger"))

# format the dates and times. This tells R that 03:45:00 is 3:45am and allows us
#   to calculate when the stop took place. Make sure that the dates and times are 
#   formatted correctly. For example, times need the ":00" seconds pasted on the end.
#   If you are missing the seconds just run the following:
# data$time <- paste(as.character(data$time),":00",sep="")
data$date <- chron(dates=as.character(data$date),format="m/d/y")
data$time <- chron(times=as.character(data$time),format="h:m:s")

# create 4 hour blocks of time
data$time.block <- cut(24*data$time,seq(0,24,4))

###############################################################################
# estimate a propensity score model, target=black, comparison=non-black drivers
###############################################################################

# create an indicator variable for the "target" group, here set to black
data$target <- as.numeric(data$race=="B")
# to do a target=Hispanic and comparison=White, for example do
#   data <- subset(data, (race=="H") | (race=="W"))
#   data$target <- as.numeric(data$race=="H")

# fit the propensity score model. The model parameters are set to defaults that 
#   are likely reasonable. For the gory details on the method and parameters 
#   type at the R prompt
# vignette("twang")
# To use your own data you will need to replace "nhood + male + ..." with the 
#    list of your own variables separated by "+". The "+" does *not* mean that
#    these variables are being added together *nor* does it mean that model
#    is linear. This is just R's notation for variables in the model 
ps1 <- ps(target~nhood + male + age + resident + reason + time.block,
          data=data,
          stop.method = stop.methods$es.stat.max,
          plots = FALSE,     
          n.trees=3000,
          shrinkage=0.005,
          interaction.depth=3,
          verbose=TRUE)

# plot the relative influence of each variable on the propensity score
#   and save the variable names and ordering
vars <- as.character(summary(ps1$gbm.obj,
                             n.trees=ps1$desc$es.stat.max$n.trees,
                             plotit=TRUE)[,1])

# show the balance table. 
# Be sure to stare at this table before proceeding. Note that at this point you
#   have not looked at any outcomes like citation or search. So the risk of 
#   the analyst biasing the results is limited. Make sure now that the column
#   labeled target and the PSadjusted column are very close (they should be).
bal.table(ps1)$es.stat.max

# create in LaTeX format. To make a table in Word, cut and paste the result 
#   into Word, highlight it, Table->Convert->Text to Table, 
#   under separate text at select other, and set it equal to &
xtable(bal.table(ps1)$es.stat.max)

###############################################################################
# propensity score analysis of citation rates
###############################################################################

# store the propensity score weights
data$w <- get.weights(ps1,stop.method="es.stat.max")

# create design object that will deal with the propensity score weights
design.ps <- svydesign(ids=~1,weights=~w,data=data)
# create a design object for the unweighted analysis for comparison
design.un <- svydesign(ids=~1,data=data)

# citation rates
# create the propensity score weighted outcomes table 
#    shows the citation rates for the comparison and target groups
temp <- svytable(~target+citation, design.ps)
temp <- (temp/rowSums(temp))[,2]
names(temp) <- c("Comparison group (adjusted)","Target group")
print(round(100*temp,1))

# Test whether the citation rates are significantly different 
svychisq(~citation+target,design=design.ps,statistic="F")

# If the balance table still had some differences remaining (e.g. black drivers 
#    were 50% while non-black drivers are 45% male) this can be adjusted with
#    covariates. Since the groups are very similar *this* weighted regression
#    is very robust to violations of the additive model structure. Most likely
#    The p-value will not change very much.
glm1.ps <- svyglm(citation ~ target + male, design=design.ps, family=quasibinomial)
summary(glm1.ps)


# Note: You can replace citation with other categorical outcome variable that 
#    you might have including type of search (none, consent, pat, probable cause)
#    or a discretization of the length of the search (less than 10 minutes, 
#    greater than 10 minutes).
# For continuous outcomes (e.g. duration of stop) svytable() does not make 
#    sense and svychisq() will not work. Simply replace all of this with the
#    following replacing "y" with the name of the variable
#  svymean(~y, design=subset(design.ps,target==1))
#  svymean(~y, design=subset(design.ps,target==0))
#  summary(svyglm(yy~target,design=design.ps))
#  glm1.ps <- svyglm(citation ~ target + male, design=design.ps, family=gaussian)


# sensitivity analysis
sensitivity(ps1,data,"citation")


###############################################################################
# Common alternatives to the propensity score analysis
###############################################################################

# The most common in racial profiling studies is the direct comparison of 
#    outcomes across races (e.g. citation rate for black drivers vs. citation
#    rate for white drivers) without any adjustment.
#    To show the dangers of the naive analysis, here's the unadjusted analysis
#    shows the citation rates for the comparison and target groups.
#    The target group citation rate will be unchanged from before.
temp <- svytable(~target+citation, design.un)
temp <- (temp/rowSums(temp))[,2]
names(temp) <- c("Comparison group (unadjusted)","Target group")
print(round(100*temp,1))


# Logistic regression comparison
#    Those studies that correctly acknowledge that some adjustment is necessary
#    often turn to logistic regression, but the results are very sensitive to 
#    the modeling assumptions (linearity on the log-odds scale) when the target 
#    and comparison groups do not overlap.
glm1 <- glm(citation~target + nhood + male + age + resident + reason + time.block,
            data=data,
            family=binomial)
summary(glm1)
# Here's an estimate of the citation rate for similarly situated white drivers 
#    relying on the logistic regression model using "recycled predictions" to 
#    put the estimate on the same scale as the propensity score analysis.
data.temp <- subset(data,target==1)
data.temp$target <- 0
mean(predict(glm1,newdata=data.temp,type="response"),na.rm=TRUE)

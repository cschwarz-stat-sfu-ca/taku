# CODE CREATED BY CARL SCHWARZ 2019
# THESE ARE CUSTOM FUNCTIONS TO IMPLEMENT BTSPAS FOR TAKU SOCKEYE
# BTSPAS package information is available at https://github.com/cschwarz-stat-sfu-ca/BTSPAS
# This function library lives at https://github.com/cschwarz-stat-sfu-ca/taku

# Source directly from there using:
#     library(devtools)
#     devtools::source_url("https://raw.githubusercontent.com/cschwarz-stat-sfu-ca/taku/master/FUNCTIONS_BTSPAS_Wrappers.R")

# Note: Running these functions first requires installing and loading BTSPAS:
# 			library(devtools)
#			devtools::install_github("cschwarz-stat-sfu-ca/BTSPAS", dependencies = TRUE,
#                        build_vignettes = TRUE, force=TRUE)
#			library(BTSPAS)
#			library(ggplot2)
#			library(plyr)

# Version Log
# 2019-05-01 GP: added some minor tweaks (default label for catch column in BTSPAS_input function, target folder in fit.BTSPAS)
# 2019-04-29 CS: Added check for n1=0 but ones added. Changed the n1=0 to n1=1.
# 2019-04-24 CS: Added the "add.ones.at.start" flat to fit.BTSPAS and fit.BTSPAS.dropout functions.
##########################################
# Create the data structures required for BTSPAS for releases
# This is stratum index, n1, m2, and u2
#   n1 is the number of tagged releases
#   m2 is a matrix with columns representing recoveries in the same stratum as release, the next stratum of release etc
#   u2 is the number of recoveries

BTSPAS_input <- function(recovery, commercial, 
                         rel.stratum, rec.stratum, 
                         stratum.index, catch.var="CdnCommCt_All"){  #GP edit: changed default to "_All"
    # recovery = data frame of releases and recoveries. One line per tagged fish
    # commercial = data frame of commercial recoveries. One line per stratum with the total count given
    #              missing strata imply a commercial catch of 0
    # rel.stratum = variable names defining the release strata
    # rec.stratum = variable names defining the recovery strata
    #             We assume that character version of stratum numbers wil properly sort, i.e. store as 01 etc
    # stratum.index has index number and stratum label.
    # only those strata that belong to stratum index will be used (to enable you truncate at front/end of study)

    # check that data frame as the required variables
    if(length(recovery[,rel.stratum])==0 | length(recovery[,rec.stratum])==0){
        stop("Input strata variables don't exists in recovery data frame")
    }
    if(length(commercial[,rec.stratum])==0){
        stop("Recovery strata variables don't exists in commercial data frame")
    }  
    #browser()
    # only want recoveries from commercial fishery. If the RecoveryType is anything else
    # then convert the rec.stratum to missing
    recovery[,rec.stratum][!recovery$RecoveryType %in% "Commercial"] <- NA
  
    # convert release stratum to character form and check for missing values
    # if the stratum has an "NA" in the stratum value, it is also considered to be missing
    recovery[, rel.stratum] <- as.character(recovery[, rel.stratum])
    # Check that no missing values in release statum
    if(sum(grepl("NA",recovery[,rel.stratum]) | is.na(recovery[,rel.stratum]))){
        stop("Missing values not allowed in release stratum")
    }
    # Check for no missing values in commercial recovery strata
        if(sum(grepl("NA",commercial[,rec.stratum]) | is.na(commercial[,rec.stratum]))){
        stop("Missing values not allowed in commercial recoveries")
        }
    
    #browser()
    # convert release stratum to character form and set to missing any stratum value that contains NA 
    recovery[, rec.stratum] <- as.character(recovery[, rec.stratum])
    recovery[, rec.stratum][ grepl("NA", recovery[, rec.stratum])] <- NA

    # only keep those release strata that appears in statum.index$labels
    recovery <- recovery[ recovery[,rel.stratum] %in% stratum.index$stratum.label,]
    # only keep those recovery strata that appears in statum.index$labels
    recovery <- recovery  [ recovery  [,rec.stratum] %in% stratum.index$stratum.label | is.na(recovery[,rec.stratum]),]
    commercial<-commercial[ commercial[,rec.stratum] %in% stratum.index$stratum.label,]

    # add the stratum index values to original data frame for releases and recoveries
    recovery <- merge(recovery, plyr::rename(stratum.index, c("stratum.index"="rel.index")), 
                   by.x=rel.stratum, by.y="stratum.label", all.x=TRUE)
    recovery <- merge(recovery, plyr::rename(stratum.index, c("stratum.index"="rec.index")), 
                   by.x=rec.stratum, by.y="stratum.label", all.x=TRUE)
 
    # compute statistics
    
    # total number of releases by the stratum index
    n1.df <- plyr::ddply(recovery, "rel.index", plyr::summarize, n1=length(rel.index))
    # impute any zeroes for no releases in a stratum. We only impute up to larges stratum index
    max.rel.index <- max(n1.df$rel.index)
    n1.df <- merge(n1.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rel.index")), all.y=TRUE)
    n1.df <- n1.df[ n1.df$rel.index <= max.rel.index,]
    n1.df$n1[ is.na(n1.df$n1)] <- 0
    n1.df <- n1.df[ order(n1.df$rel.index),]
    #browser()
    # get the cross classification of releases x recoveries. We first need to insert 0's in missing combinations
    m2.large <- plyr::ddply(recovery, "rel.index", function(x){
       m2.df <- plyr::ddply(x, "rec.index", plyr::summarize, m2=length(rec.index))
       # impute any zeroes for no recoveries in a stratum. 
       m2.df <- merge(m2.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rec.index")),
                      all.y=TRUE)
       m2.df$m2[ is.na(m2.df$m2)] <- 0
       m2.df
    })
    # insert 0 for any release strata with no tags put out
    m2.all.comb <- expand.grid(rel.index= n1.df$rel.index, rec.index=stratum.index$stratum.index)
    m2.large <- merge(m2.large, m2.all.comb, all.y=TRUE)
    m2.large$m2[ is.na(m2.large$m2)] <- 0
    
    # now we cross classify
    m2.full <- as.matrix(xtabs(m2~rel.index+rec.index, data=m2.large))
    
    # rotate each row to the left because you can't recover before releases
    # we put the amount of (shift+1) in the first column
    m2.red <- plyr::aaply(cbind(n1.df$rel.index,m2.full), 1, function(x){
        #browser()
        rot.vec <- c((x[1]+1):length(x), 2:(x[1]+1))
        x <- x[rot.vec]
        x[-length(x)]
    })
    #browser()
    # remove columns at the right that are all zero
    all.zero <- apply(m2.red==0, 2, all)
    remove.right <- rev(cumprod(rev(all.zero)))
    m2.red <- m2.red[, !remove.right]
    
    # remove columns at the left that are all zero
    all.zero <- apply(m2.red==0, 2, all)
    remove.left <- cumprod(all.zero)
    m2.red <- m2.red[, !remove.left]
  
    # Find the commercial recoveries by rec.stratum
    # Merge the stratum index
    # browser()
    commercial <- merge(commercial, plyr::rename(stratum.index, c("stratum.index"="rec.index")), 
                   by.x=rec.stratum, by.y="stratum.label", all.x=TRUE)
    #browser()
    u2.df <- plyr::ddply(commercial, "rec.index", function(x){
        u2=sum(x[,catch.var])
        data.frame(u2=u2)
    })
    # impute zeros for missing recovery strata
    u2.df <-merge(u2.df,  plyr::rename(stratum.index[,"stratum.index", drop=FALSE], c("stratum.index"="rec.index")), all.y=TRUE)
    u2.df$u2[ is.na(u2.df$u2)] <- 0
    u2.df <- u2.df[ order(u2.df$rec.index),]
    #browser()
    commercial.catch.df <- u2.df 
    # subtract the recoveries of marked fish
    u2.df.before.correction <- u2.df$u2 - apply(m2.full,2,sum)
    u2.df$u2 <- pmax(0,u2.df.before.correction)
   
    #browser()
    # round everything to integer values
    n1.df <- round(n1.df)
    m2.full <- round(m2.full)
    m2.red  <- round(m2.red)
    commercial.catch.df <- round(commercial.catch.df)
    u2.df.before.correction <- round(u2.df.before.correction)
    u2.df <- round(u2.df)
    
    list(Year=commercial$Year[1],
         recovery=recovery,
         commercial=commercial, 
         stratum.index=stratum.index,
         n1.df=n1.df,
         m2.full=m2.full,
         m2.red=m2.red,
         commercial.catch.df=commercial.catch.df,
         u2.df.before.correction=u2.df.before.correction,
         u2.df=u2.df)
}





#####################################################
#####################################################
######################################
# Make the call to BTSPAS

fit.BTSPAS <- function( input.data, prefix="BTSPAS", folder=NULL , #GP edit: added folder argument
						debug=FALSE, increase.iterations.factor=1,
                        add.ones.at.start=FALSE, add.ones.at.end=add.ones.at.start){
   # many.more.iterations = increase interations by factor of x
   # add.ones.at.start = TRUE. Assumes that first few years of releases had no recoveries
   #           because commercial fishery started later due to small numbers of fish.
   #           It will add 1 to the diagonal elements for the first few releases to make
   #           the full m2.full matrix diagonal with at least 1's along the diagonal
  
   # takes the input data created earlier and fits the BTSPAS model
   #browser()
   cat("Starting to fit year ", input.data$Year, "\n")
   dir.name <- paste(prefix,"-", input.data$Year, sep="") # GP Modified
   if(!is.null(folder)){dir.name <- paste(folder,dir.name,sep="/")} # GP Modified

   cat("Directory name is ", dir.name, "\n")
   if(file.access(dir.name)!=0){
     dir.create(dir.name, showWarnings=TRUE)
   }
   setwd(dir.name)
   taku.prefix <- paste(prefix,"-",input.data$Year, sep="")
   taku.title  <- paste(taku.prefix," TSPND NP")

   # see if we need to add 1 to the first few rows to make the recoveries
   # diagonal
   # browser()
   # add a 1 to diagonal for every row with no recoveries from any releases at start of recoveries
   if(add.ones.at.start){
      total.recaps <- apply(input.data$m2.full,2,sum)
      index.to.add <- (1:nrow(input.data$m2.red))[as.logical(cumprod(total.recaps==0))]
      input.data$m2.red[index.to.add,1] <- 1
      input.data$m2.full[ cbind(index.to.add, index.to.add)] <- 1
      input.data$n1.df$n1[index.to.add]<- ifelse(input.data$n1.df$n1[index.to.add]==0,1,input.data$n1.df$n1[index.to.add])
   }
   
   # add a 1 to the diagonal for every row with no recoveries from any releases at END of recoveries
   # The default action is the same as add.ones.at.start, but can be changed
   if(add.ones.at.end){
      #browser()
      total.recaps <- apply(input.data$m2.full,2,sum)[1:nrow(input.data$n1.df)]  # only diagonal elements have 1's added
      index.to.add <- rev((1:nrow(input.data$m2.red)))[as.logical(cumprod(rev(total.recaps==0)))]
      input.data$m2.red[index.to.add,1] <- 1
      input.data$m2.full[ cbind(index.to.add, index.to.add)] <- 1
      input.data$n1.df$n1[index.to.add]<- ifelse(input.data$n1.df$n1[index.to.add]==0,1,input.data$n1.df$n1[index.to.add])
   }
   
   #write out a csv file with the releases and recoveries in one giant data frame
   #browser()
   n1.df.copy <- input.data$n1.df
   m2.full.copy <- as.data.frame.matrix(input.data$m2.full)
   colnames(m2.full.copy) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
   
   out.df <- cbind(n1.df.copy, m2.full.copy)
   u2.copy <- as.data.frame.matrix(t(input.data$u2.df$u2))
   colnames(u2.copy)<- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
   u2.copy$n1 <- NA
   u2.copy$rel.index <- NA
   
   total.cap <- as.data.frame(t(unlist(apply(input.data$m2.full,2,sum) + input.data$u2.df$u2)))
   colnames(total.cap) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
   total.cap$n1 <- NA
   total.cap$rel.index <- NA
   
   out.df <- rbind(out.df, u2.copy, total.cap)
   row.names(out.df) <- c(input.data$stratum.index$stratum.label[1:nrow(n1.df.copy)],"untagged","total.captures")
    
   write.csv(out.df, paste(taku.prefix,"-SummaryStatisticsInMatrixForm.csv",sep=""), row.names=TRUE)
   

   # what is the strata identification number (statistical week from start of year)?
   taku.sweek <- input.data$stratum.index$stratum.index
   
   # releases, recoveries, and commercial catch
   taku.n1 <- input.data$n1.df$n1
   taku.m2 <- input.data$m2.red  # we want the reduced matrix
   taku.u2 <- input.data$u2.df$u2

   # are there any jumps in the abundance?
   taku.jump.after <- NULL    # list sample times after which jump in number occurs

   # are there any bad values that need to be adjusted?
   taku.bad.n1     <- c()     # list sample times of bad n1 values
   taku.bad.m2     <- c()     # list sample times of bad m2 values
   taku.bad.u2     <- c()     # list sample times of bad u2 values

   # are there any days where the capture probability is fixed in advance?, i.e. because no commercial fishery
#  taku.logitP.fixed        <- taku.sweek[which(taku.u2 ==0)]
   taku.logitP.fixed        <- taku.sweek[which(total.cap ==0)]
   taku.logitP.fixed.values <- rep(-10, length(taku.logitP.fixed))
   #browser()
   ## Run TSPNDE model
   taku.fit.tspndenp <- TimeStratPetersenNonDiagErrorNP_fit(
                  title=      taku.title,
                  prefix=     taku.prefix,
                  time=       taku.sweek,
                  n1=         taku.n1,
                  m2=         taku.m2,
                  u2=         taku.u2,
                  jump.after= taku.jump.after,
                  bad.n1=     taku.bad.n1,
                  bad.m2=     taku.bad.m2,
                  bad.u2=     taku.bad.u2,
                  logitP.fixed=taku.logitP.fixed,
                  logitP.fixed.values=taku.logitP.fixed.values,
                  n.iter=ifelse(debug,1000,20000*increase.iterations.factor), 
                  n.burnin=ifelse(debug,100,1000), n.sims=ifelse(debug,30,500),
                  debug=FALSE
                  )
   # Rename files that were created.

   file.rename("data.txt",       paste(taku.prefix,".data.txt",sep=""))
   file.rename("CODAindex.txt",  paste(taku.prefix,".CODAindex.txt",sep=""))
   file.rename("CODAchain1.txt", paste(taku.prefix,".CODAchain1.txt",sep=""))
   file.rename("CODAchain2.txt", paste(taku.prefix,".CODAchain2.txt",sep=""))
   file.rename("CODAchain3.txt", paste(taku.prefix,".CODAchain3.txt",sep=""))
   file.rename("inits1.txt",     paste(taku.prefix,".inits1.txt",sep=""))
   file.rename("inits2.txt",     paste(taku.prefix,".inits2.txt",sep=""))
   file.rename("inits3.txt",     paste(taku.prefix,".inits3.txt",sep=""))


   # Save the information for later retreival if needed
   save(list=c("taku.fit.tspndenp"), file="taku-fit-tspndenp-saved.Rdata")  # save the results from this run
   setwd("..")# return back
}

#####################################################
#####################################################
######################################
# Make the call to BTSPAS allowing for dropout/fallback

fit.BTSPAS.dropout <- function( input.data, prefix="BTSPAS", folder=NULL , #GP edit: added folder argument
								debug=FALSE, n, dropout,
                                add.ones.at.start=FALSE, add.ones.at.end=add.ones.at.start){
   # takes the input data created earlier and fits the BTSPAS model
   # add.ones.at.start = TRUE. Assumes that first few years of releases had no recoveries
   #           because commercial fishery started later due to small numbers of fish.
   #           It will add 1 to the diagonal elements for the first few releases to make
   #           the full m2.full matrix diagonal with at least 1's along the diagonal
  #browser()
  cat("Starting to fit year ", input.data$Year, "\n")
  dir.name <- paste(prefix,"-", input.data$Year, sep="") # GP Modified
  if(!is.null(folder)){dir.name <- paste(folder,dir.name,sep="/")} # GP Modified
  
  
  cat("Directory name is ", dir.name, "\n")
  if(file.access(dir.name)!=0){
    dir.create(dir.name, showWarnings=TRUE)
  }
  setwd(dir.name)
  taku.prefix <- paste(prefix,"-",input.data$Year, sep="")
  taku.title  <- paste(taku.prefix," TSPND NP")
  
   # see if we need to add 1 to the first few rows to make the recoveries
   # diagonal
   # browser()
   # add a 1 to diagonal for every row with u2=2 at start of recoveries
   if(add.ones.at.start){
      total.recaps <- apply(input.data$m2.full,2,sum)
      index.to.add <- (1:nrow(input.data$m2.red))[as.logical(cumprod(total.recaps==0))]
      input.data$m2.red[index.to.add,1] <- 1
      input.data$m2.full[ cbind(index.to.add, index.to.add)] <- 1
      input.data$n1.df$n1[index.to.add]<- ifelse(input.data$n1.df$n1[index.to.add]==0,1,input.data$n1.df$n1[index.to.add])
   }
  
   # add a 1 to the diagonal for every row with no recoveries from any releases at END of recoveries
   # The default action is the same as add.ones.at.start, but can be changed
   if(add.ones.at.end){
      #browser()
      total.recaps <- apply(input.data$m2.full,2,sum)[1:nrow(input.data$n1.df)]  # only diagonal elements have 1's added
      index.to.add <- rev((1:nrow(input.data$m2.red)))[as.logical(cumprod(rev(total.recaps==0)))]
      input.data$m2.red[index.to.add,1] <- 1
      input.data$m2.full[ cbind(index.to.add, index.to.add)] <- 1
      input.data$n1.df$n1[index.to.add]<- ifelse(input.data$n1.df$n1[index.to.add]==0,1,input.data$n1.df$n1[index.to.add])
   }
  #write out a csv file with the releases and recoveries in one giant data frame
  #browser()
  n1.df.copy <- input.data$n1.df
  m2.full.copy <- as.data.frame.matrix(input.data$m2.full)
  colnames(m2.full.copy) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  
  out.df <- cbind(n1.df.copy, m2.full.copy)
  u2.copy <- as.data.frame.matrix(t(input.data$u2.df$u2))
  colnames(u2.copy)<- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  u2.copy$n1 <- NA
  u2.copy$rel.index <- NA
  
  total.cap <- as.data.frame(t(unlist(apply(input.data$m2.full,2,sum) + input.data$u2.df$u2)))
  colnames(total.cap) <- input.data$stratum.index$stratum.label[1:ncol(m2.full.copy)]
  total.cap$n1 <- NA
  total.cap$rel.index <- NA
  
  out.df <- rbind(out.df, u2.copy, total.cap)
  row.names(out.df) <- c(input.data$stratum.index$stratum.label[1:nrow(n1.df.copy)],"untagged","total.captures")
  
  write.csv(out.df, paste(taku.prefix,"-SummaryStatisticsInMatrixForm.csv",sep=""), row.names=TRUE)
  
  
  # what is the strata identification number (statistical week from start of year)?
  taku.sweek <- input.data$stratum.index$stratum.index
  
  # releases, recoveries, and commercial catch
  taku.n1 <- input.data$n1.df$n1
  taku.m2 <- input.data$m2.red  # we want the reduced matrix
  taku.u2 <- input.data$u2.df$u2
  
  # are there any jumps in the abundance?
  taku.jump.after <- NULL    # list sample times after which jump in number occurs
  
  # are there any bad values that need to be adjusted?
  taku.bad.n1     <- c()     # list sample times of bad n1 values
  taku.bad.m2     <- c()     # list sample times of bad m2 values
  taku.bad.u2     <- c()     # list sample times of bad u2 values
  
  # are there any days where the capture probability is fixed in advance?, i.e. because no commercial fishery
  #taku.logitP.fixed        <- taku.sweek[which(taku.u2 ==0)]
  taku.logitP.fixed        <- taku.sweek[which(total.cap ==0)]
  taku.logitP.fixed.values <- rep(-10, length(taku.logitP.fixed))
  #browser()
  ## Run TSPNDE model
  taku.fit.tspndenp <- TimeStratPetersenNonDiagErrorNPMarkAvail_fit(
    title=      taku.title,
    prefix=     taku.prefix,
    time=       taku.sweek,
    n1=         taku.n1,
    m2=         taku.m2,
    u2=         taku.u2,
    jump.after= taku.jump.after,
    bad.n1=     taku.bad.n1,
    bad.m2=     taku.bad.m2,
    bad.u2=     taku.bad.u2,
    logitP.fixed=taku.logitP.fixed,
    logitP.fixed.values=taku.logitP.fixed.values,
    n.iter=ifelse(debug,1000,30000), n.burnin=ifelse(debug,100,1000), n.sims=ifelse(debug,30,500),
    marked_available_n=n, marked_available_x=n-dropout,
    debug=FALSE
  )
  # Rename files that were created.
  
  file.rename("data.txt",       paste(taku.prefix,".data.txt",sep=""))
  file.rename("CODAindex.txt",  paste(taku.prefix,".CODAindex.txt",sep=""))
  file.rename("CODAchain1.txt", paste(taku.prefix,".CODAchain1.txt",sep=""))
  file.rename("CODAchain2.txt", paste(taku.prefix,".CODAchain2.txt",sep=""))
  file.rename("CODAchain3.txt", paste(taku.prefix,".CODAchain3.txt",sep=""))
  file.rename("inits1.txt",     paste(taku.prefix,".inits1.txt",sep=""))
  file.rename("inits2.txt",     paste(taku.prefix,".inits2.txt",sep=""))
  file.rename("inits3.txt",     paste(taku.prefix,".inits3.txt",sep=""))
  
  
  # Save the information for later retreival if needed
  save(list=c("taku.fit.tspndenp"), file="taku-fit-tspndenp-saved.Rdata")  # save the results from this run
  setwd("..")# return back
}


#####################################################
#####################################################
BTSPAS_input_from_matrix <- function(sheet, workbook){
  # read the worksheet from the workbook and creates the data structures for BTSPAS input from
  # summary data at the stat week for prior years
  # 
  # Format of the matrix
  #    Row 1. Stat week of recovery (starting second column)
  #    Column 1. Stat week of release (numeric)
  #    intersection of stat week of release and recover= # of tags captured
  #    Then there can be several rows for the total recoveries, test fisheries, etc until
  #        final row which is TOTAL harvest/test fisher INCLUDING all tags. This MUST BE LABELLED as "CatchTotal" or "Total catch" in column A
  #    final two columns are total recoveries (never used) and total releases.
  # 
  # Stat week of recovery must start at te first stat week of release (yes, many of the first columns will be zero)
  # The last statweek of recovery must be at least as large as the last stat week of release.
  
  # Blanks indicate 0 
  #browser()
  cat("Extracting matrix from ", workbook, "; sheet:",sheet, "\n")
  # get the stat weeks of releases from the first column
  stat.week.rel <- readxl::read_excel(workbook, 
                                      sheet=sheet,
                                      col_names=FALSE,
                                      .name_repair="universal",
                                      range=cellranger::cell_cols("A"))
  stat.week.rel <- as.numeric(stat.week.rel[,1,drop=TRUE])
  stat.week.rel <- stat.week.rel[!is.na(stat.week.rel)]
  
  # get the stat week of recoveries from the first row
  stat.week.rec <- readxl::read_excel(workbook,
                                      sheet=sheet,
                                      col_names=FALSE,
                                      .name_repair="universal",
                                      range=cellranger::cell_rows(1))
  stat.week.rec <- as.numeric(stat.week.rec[1,,drop=TRUE])
  stat.week.rec <- stat.week.rec[ !is.na(stat.week.rec)]
  
  # get the number of releases
  colA.entries <- readxl::read_excel(workbook, 
                                     sheet=sheet,
                                     col_names=FALSE,
                                     .name_repair="universal",
                                     range=cellranger::cell_cols("A"))
  #browser()
  colA.entries <- colA.entries[!is.na(as.numeric(colA.entries[,1,drop=TRUE])),,drop=TRUE]

  n1 <- readxl::read_excel(workbook,
                           sheet=sheet,
                           col_names=FALSE,
                           .name_repair="universal",
                           range=cellranger::cell_cols(3+length(stat.week.rec)))
  n1 <- as.numeric(n1[,1,drop=TRUE])
  n1 <- n1[ -1]  # drop the first entry
  n1 <- n1[ 1:length(colA.entries)] # only those rows with numeric statweek
  n1[ is.na(n1)] <- 0  # missing values are assumed to be zero
  n1.df <- data.frame(n1=n1[1:length(stat.week.rel)], rel.index=1:length(stat.week.rel))
  
  # get the full recovery matrix
  m2.full <- readxl::read_excel(workbook,
              sheet=sheet, col_names=FALSE,
              .name_repair="universal",
              range=cellranger::cell_limits(ul=c(2,2),
                                            lr=c(1+length(n1),
                                                 1+length(stat.week.rec))) )
  m2.full[ is.na(m2.full)] <- 0
  m2.full <- as.matrix(m2.full)
  #browser()
  # Now to compute the reduced recovery matrix  by shifting rows to the left 
  m2.red <- plyr::aaply(cbind(n1.df$rel.index,m2.full), 1, function(x){
        #browser()
        rot.vec <- c((x[1]+1):length(x), 2:(x[1]+1))
        x <- x[rot.vec]
        x[-length(x)]
  })
  # remove columns at the right that are all zero
  all.zero <- apply(m2.red==0, 2, all)
  remove.right <- rev(cumprod(rev(all.zero)))
  m2.red <- m2.red[, !remove.right]
    
  # remove columns at the left that are all zero
  all.zero <- apply(m2.red==0, 2, all)
  remove.left <- cumprod(all.zero)
  m2.red <- m2.red[, !remove.left]
  
  # get the total number of recoveries
  # This is the last row in the spreadsheet. We wil ignore any rows between the
  # end of the statweeks and the last row
  colA.entries <- readxl::read_excel(workbook, 
                                      sheet=sheet,
                                      col_names=FALSE,
                                      .name_repair="universal",
                                      range=cellranger::cell_cols("A"))
  #browser()
  colA.entries <- colA.entries[!is.na(colA.entries[,1,drop=TRUE]),,drop=TRUE]
  last.row <- length(colA.entries)
  n2 <- readxl::read_excel(workbook,
                           sheet=sheet,
                            col_names=FALSE,
                            .name_repair="universal",
                            range=cellranger::cell_rows(last.row))
  n2 <- as.numeric(n2[1,,drop=TRUE])
  n2 <- n2[-1] # drop the first column
  n2 <- n2[ 1:length(stat.week.rec)]
  
  u2.df <- data.frame(u2=n2-apply(m2.full,2,sum ))
  
  
  #browser()

 
  list(Year=sheet,
       stratum.index=data.frame(stratum.index=1:length(stat.week.rec), stratum.label=stat.week.rec),
       n1.df=as.data.frame(n1.df),
       m2.full=m2.full,
       m2.red=as.matrix(m2.red),
       commercial.catch.df=data.frame(commercial.catch=n2),
       u2.df=as.data.frame(u2.df)) 
  
}


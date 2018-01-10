#' @export
Get.Station.Data.9var <- function(stnid, wdir) {

  #wdir = obsdir   ####"F:/SForecast-TestRun/Database/obs_samples/asos-daily"
  #stnid = "ID090"

  setwd(wdir)

  srchstr = paste("*", stnid, "*.csv", sep="")
  obsfile <- list.files(wdir, pattern = glob2rx(srchstr), full.names=F)

  obs = read.csv(obsfile, header=T, na.strings = c("-99.00", "NA"))
  obs = obs[which(!is.na(obs$Year)),]
  colcnt = ncol(obs)
  datestr <- paste(obs[[1]],obs[[2]],obs[[3]],sep="-")
  date <- as.Date(datestr,"%Y-%m-%d")

  if(colcnt == 9){
    obs = cbind(date, obs[,c(4:9)])
    obs[,8:10] = NA
  }

  if(colcnt == 12){
    obs = cbind(date, obs[,c(4:12)])
  }
  year=as.numeric (format(obs[,1],"%Y"))
  month=as.numeric (format(obs[,1],"%m"))
  yearmon=as.character (format(obs[,1],"%Y-%m"))
  data <- cbind(year, month, yearmon, obs)

  return(data)
}

#' @export
Cal.Monthly.Station.Mean <- function(EnvList, varnms){


  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  obsdir = EnvList$obsdir
  syear_obs = as.numeric(EnvList$syear_obs)

  VarDFile = paste(vardir, "/", varfile, sep="")

  ###### Get Station ID, lat, and Lon information
  setwd(stndir)
  stninfo = read.csv(stnfile, header=T)
  colnames(stninfo) = c("Lon","Lat", "Elev", "ID", "Ename", "SYear")
  stninfo = stninfo[which(stninfo$SYear <= syear_obs),c("ID", "Lon", "Lat", "SYear")]
  stnnms = matrix(stninfo$ID)

  stncnt = length(stnnms)
  cat("process now running...Station...");cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  for(j in 1:stncnt){
    setTxtProgressBar(pb,j)
    stnid = stnnms[j]

    obsdaydata = Get.Station.Data.9var(stnid, obsdir)
    colnames(obsdaydata) = c("year", "month", "yearmon", "date", "prec", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")
    obsdaydata$t2m = (obsdaydata$tmax + obsdaydata$tmin)/2
    obsdaydata$stnid = stnid

    # select common period
    #obsdaydata = obsdaydata[which(obsdaydata$year >= syear_obs & obsdaydata$year <= eyear_obs),]
    obsdaydata = obsdaydata[which(obsdaydata$year >= syear_obs),]

    if(j == 1){
      data = obsdaydata
    } else {
      temp = obsdaydata
      data = rbind(data, temp)
    }
  }
  close(pb)
  mdata = melt(data, id=c("stnid", "yearmon", "year", "month", "date"))

  outdata = summaryBy(value~variable+yearmon, data=mdata, FUN = mean, na.rm=T)
  meandata = reshape(outdata, idvar = c("yearmon"), timevar = "variable", direction = "wide")
  colnames(meandata) = c("yearmon", "prec", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg", "t2m")

  ##### just 2 varialbes #########################
  #meandata = meandata[c("yearmon", "prec", "t2m")]
  meandata = meandata[c("yearmon", varnms)]
  setwd(vardir)
  write.csv(meandata, varfile, row.names=F)

}

#' @export
Fill.Missing.KMA.Daily.9var <- function(EnvList){

  #obsdir = "F:/StationDB/weather/kma_asos_9var"
  #smpldir = "F:/SForecast-RT/Database/kma-asos_9var"
  #nrstfile = "kma_asos_57stns_nearest_10stns.csv"
  #syear_obs = 1976
  #eyear = 2013
  obsdir = EnvList$obsdir
  smpldir = EnvList$smpldir
  stndir = EnvList$stndir
  syear = as.numeric(EnvList$syear_obs)
  nrstfile = EnvList$nrstfile
  eyear = as.numeric(EnvList$eyear_sim)

  NrstDFile = paste(stndir, "/", nrstfile, sep="")
  nrstdata = read.csv(NrstDFile, header=T)

  # Get Station List
  stnnms = unique(nrstdata$InputID)
  stncnt = length(stnnms)

  # Read original file, fill the missing values, save the adjusted file
  cat("process now running...");cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  for(i in 1:stncnt){

    stnnm = stnnms[i]
    stnid = sprintf("ID%03d", stnnm)
    SrcDFile = paste(obsdir, "/", stnid, ".csv", sep="")

    obs = read.csv(SrcDFile, header=T, na.strings = "-99.00")
    orgheader = colnames(obs)

    datestr <- paste(obs[[1]],obs[[2]],obs[[3]],sep="-")
    date <- as.Date(datestr,"%Y-%m-%d")
    obs = cbind(date, obs[,c(4:12)])

    year=as.numeric (format(obs[,1],"%Y"))
    month=as.numeric (format(obs[,1],"%m"))
    yearmon=as.character (format(obs[,1],"%Y-%m"))
    data <- cbind(year, month, yearmon, obs)

    firstdate = as.Date(paste(syear, "-1-1", sep=""))
    lastdate = as.Date(paste(eyear, "-12-31", sep=""))
    #lastdate = Sys.Date()
    date = as.data.frame(seq(firstdate, lastdate, by="day")); colnames(date) = c("date")

    data = merge(date, data, by="date", all=T)
    data = data[which(as.numeric(substr(data$date, 1, 4)) >= syear & as.numeric(substr(data$date, 1, 4)) <= eyear), ]

    colnames(data) = c("date", "year", "month", "yearmon", "prec", "tmax", "tmin", "wind", "rhum", "srad", "sshine", "cloud", "tavg")

    varnms = names(data[, c(5:13)])

    nstnnms = nrstdata[which(nrstdata$InputID == stnnm),]

    #data = Fill.Missing.Daily (data, varnms, nstnnms, obsdir)
    data = Fill.Missing.Daily (data, varnms, nstnnms, obsdir, syear, eyear)

    data$Day = as.numeric(substr(data$date, 9, 10))

    data = data[c("year", "month", "Day", "prec", "tmax", "tmin", "wind", "rhum", "srad", "sshine", "cloud", "tavg")]
    colnames(data) = orgheader

    outdir = paste(smpldir, "/asos-daily", sep="")
    SetWorkingDir(outdir)
    OutDFile = paste(outdir, "/", stnid, ".csv", sep="")
    write.csv(data, OutDFile, row.names=F)
    setTxtProgressBar(pb,i)

  }
  close(pb)
}

#' @export
Fill.Missing.Daily <- function(data, varnms, nstnnms, srcdir, syear, eyear){

  #varnms = varnms[!varnms %in% "srad"]

  varcnt = length(varnms)

  for(k in 1:varcnt){
    varnm = varnms[k]

    firstdate = as.Date(paste(syear, "-1-1", sep=""))
    lastdate = as.Date(paste(eyear, "-12-31", sep=""))
    #firstdate = as.Date(data[1, "date"])
    #lastdate = as.Date(data[nrow(data), "date"])

    if(sum(is.na(data[,c(varnm)])) > 0){
      nrcnt = 1

      repeat {
        nstnnm = nstnnms[nrcnt, c("TargetID")]
        nstnid = sprintf("ID%03d", nstnnm)

        NstnDFile = paste(srcdir, "/", nstnid, ".csv", sep="")

        if(file.exists(NstnDFile)){

          obs = read.csv(NstnDFile, header=T, na.strings = "-99.00")
          datestr <- paste(obs[[1]],obs[[2]],obs[[3]],sep="-")
          date <- as.Date(datestr,"%Y-%m-%d")
          obs = cbind(date, obs[,c(4:12)])

          year=as.numeric (format(obs[,1],"%Y"))
          month=as.numeric (format(obs[,1],"%m"))
          yearmon=as.character (format(obs[,1],"%Y-%m"))
          ndata <- cbind(year, month, yearmon, obs)

          colnames(ndata) = c("year", "month", "yearmon", "date", "prec", "tmax", "tmin", "wind", "rhum", "srad", "sshine", "cloud", "tavg")

          date = as.data.frame(seq(firstdate, lastdate, by="day")); colnames(date) = c("date")

          ndata = merge(date, ndata, by="date", all=T)
          ndata = ndata[which(as.numeric(substr(ndata$date, 1, 4)) >= syear & as.numeric(substr(ndata$date, 1, 4)) <= eyear), ]
          #ndata = ndata[which(as.numeric(substr(ndata$date, 1, 4)) >= syear), ]

          if(nrow(data) == nrow(ndata)){
            data[is.na(data[, c(varnm)]), c("date", varnm)] = ndata[is.na(data[, c(varnm)]), c("date", varnm)]
          } else {
            cat(sprintf("     OBS: Number of rows are different: nearest = %s!\n", nstnid))
          }

        }

        nrcnt = nrcnt + 1

        if((sum(is.na(data[,c(varnm)])) == 0) | nrcnt == nrow(nstnnms)){

          # if there is sill missing after considering nearest 10 stations, fill the missing using interpolation
          if(sum(is.na(data[,c(varnm)])) <= (length(data[,c(varnm)]-2))) {data[,c(varnm)] = na.approx(data[,c(varnm)], na.rm=F)}

          break

        }

      }

    }


  }

  return(data)

}

#' @export
Fill.Missing.KMA.Hourly.9var <- function(EnvList){

  #srcdir = "F:/StationDB/weather/kma_asos_hourly/station"
  #outdir = "F:/SForecast-RT/Database/kma-asos-hourly"
  #nrstfile = "kma_asos_63stns_nearest_10stns.csv"
  #syear = 1976
  #eyear = 2013
  obsdir = EnvList$obsdir
  srcdir = paste(obsdir, "/asos_hourly_download/station", sep="")
  smpldir = EnvList$smpldir
  stndir = EnvList$stndir
  nrstfile = EnvList$nrstfile
  syear = as.numeric(EnvList$syear_obs)
  eyear = as.numeric(EnvList$eyear_obs)


  options(warn=-1)
  options(stringsAsFactors = FALSE)

  outdir = paste(smpldir, "/asos-hourly", sep="")
  DSumDir = paste(outdir, "/daily_summary", sep="")
  SetWorkingDir(outdir)
  SetWorkingDir(DSumDir)

  NrstDFile = paste(stndir, "/", nrstfile, sep="")
  nrstdata = read.csv(NrstDFile, header=T)

  # Get Station List
  stnnms = unique(nrstdata$InputID)
  stncnt = length(stnnms)

  # Read original file, fill the missing values, save the adjusted file
  for(i in 1:stncnt){

    stnnm = stnnms[i]
    stnid = sprintf("ID%03d", stnnm)
    SrcDFile = paste(srcdir, "/", stnid, "-Hourly.csv", sep="")

    if(file.exists(SrcDFile)){
      # hourly data contains 30min interval data, so remove them
      obs = read.csv(SrcDFile, header=T, na.strings = "-999")
      ostart = as.POSIXct(paste(obs$datetime[1], ":00", sep=""))
      oend = as.POSIXct(paste(obs$datetime[nrow(obs)], ":00", sep=""))
      odatetime = seq(ostart, oend, by="hour")
      odatetime = as.data.frame(substr(odatetime, 1, 16))
      colnames(odatetime) = c("datetime")
      obs = merge(odatetime, obs, by ="datetime")
      obs = obs[which(as.numeric(substr(obs$datetime, 1, 4)) >= syear & as.numeric(substr(obs$datetime, 1, 4)) <= eyear), ]
      obs = unique(obs)

      fdatetime = as.POSIXct(sprintf("%4d-01-01 00:00:00", syear))
      ldatetime = as.POSIXct(sprintf("%4d-12-31 23:00:00", eyear))
      datetime = seq(fdatetime, ldatetime, by="hour")
      datetime = as.data.frame(substr(datetime, 1, 16))
      colnames(datetime) = c("datetime")

      data = merge(datetime, obs, by="datetime", all=T)
      data$stnid = stnnm

      varnms = names(data[, c(3:11)])

      nstnnms = nrstdata[which(nrstdata$InputID == stnnm),]

      data = Fill.Missing.Hourly (data, varnms, nstnnms, srcdir, syear, eyear)

      data = unique(data)

      # Change columns and header based on NIER modeling system
      hdata = data[c("datetime", "stnid", "wind", "pres", "atem", "dewp", "rhum", "prec", "clou", "solr")]
      colnames(hdata) = c("TM", "STNID", "WIND", "PRES", "T.m", "Tdew.m", "RH", "P.m", "CC", "RAD")

      OutDFile = paste(outdir, "/", stnid, ".csv", sep="")
      write.csv(hdata, OutDFile, row.names=F)

      # Create daily summary
      data$date = substr(data$datetime, 1, 10)
      ddata_avg = aggregate(cbind(stnid, wind, pres, atem, dewp, rhum, clou) ~ date, data = data, FUN = mean)
      ddata_sum = aggregate(cbind(prec, solr, sshine) ~ date, data = data, FUN = sum)
      ddata_max = aggregate(atem ~ date, data = data, FUN = max); colnames(ddata_max) = c("date", "tmax")
      ddata_min = aggregate(atem ~ date, data = data, FUN = min); colnames(ddata_min) = c("date", "tmin")
      ddata = merge(ddata_sum, ddata_avg, by="date")
      ddata = merge(ddata, ddata_max, by="date")
      ddata = merge(ddata, ddata_min, by="date")
      ddata$year = substr(ddata$date,1,4)
      ddata$Mon = substr(ddata$date,6,7)
      ddata$Day = substr(ddata$date,9,10)
      ddata$rhum = ddata$rhum/100
      ddata = ddata[c("year", "Mon", "Day", "prec", "tmax", "tmin", "wind", "rhum", "solr", "sshine", "clou", "atem")]

      ddata$prec = sprintf("%5.1f", ddata$prec)
      ddata$tmax = sprintf("%5.1f", ddata$tmax)
      ddata$tmin = sprintf("%5.1f", ddata$tmin)
      ddata$wind = sprintf("%5.1f", ddata$wind)
      ddata$solr = sprintf("%5.1f", ddata$solr)
      ddata$sshine = sprintf("%5.1f", ddata$sshine)
      ddata$clou = sprintf("%5.1f", ddata$clou)
      ddata$atem = sprintf("%5.1f", ddata$atem)
      ddata$rhum = sprintf("%5.3f", ddata$rhum)

      colnames(ddata) = c("Year", "Mon", "Day", "Pcp(mm)", "Tmax(C)", "Tmin(C)", "WSpeed(m/s)", "Humidity(fr)", "SRad(MJ/m2)", "SShine(hr)", "Cloud", "Tavg(C)")

      OutDFile = paste(DSumDir, "/", stnid, ".csv", sep="")
      write.csv(ddata, OutDFile, row.names=F)

      cat(sprintf("     OBS: Filling missing values has been finished: Station = %s!\n", stnid))

    } else {
      cat(sprintf("     OBS: Observation file does not exist: Station = %s!\n", stnid))
    }

  }

}

#' @export
Fill.Missing.Hourly <- function(data, varnms, nstnnms, srcdir, syear, eyear){

  options(warn=-1)
  options(stringsAsFactors = FALSE)


  varcnt = length(varnms)
  for(k in 1:varcnt){
    varnm = varnms[k]



    if(sum(is.na(data[,c(varnm)])) > 0){
      nrcnt = 1

      repeat {
        nstnnm = nstnnms[nrcnt, c("TargetID")]
        nstnid = sprintf("ID%03d", nstnnm)

        NstnDFile = paste(srcdir, "/", nstnid, "-Hourly.csv", sep="")

        if(file.exists(NstnDFile)){

          obs = read.csv(NstnDFile, header=T, na.strings = "-999")
          ostart = as.POSIXct(paste(obs$datetime[1], ":00", sep=""))
          oend = as.POSIXct(paste(obs$datetime[nrow(obs)], ":00", sep=""))
          odatetime = seq(ostart, oend, by="hour")
          odatetime = as.data.frame(substr(odatetime, 1, 16))
          colnames(odatetime) = c("datetime")

          obs = merge(odatetime, obs, by ="datetime")
          obs = obs[which(as.numeric(substr(obs$datetime, 1, 4)) >= syear & as.numeric(substr(obs$datetime, 1, 4)) <= eyear), ]
          obs = unique(obs)

          fdatetime = as.POSIXct(sprintf("%4d-01-01 00:00:00", syear))
          ldatetime = as.POSIXct(sprintf("%4d-12-31 23:00:00", eyear))
          datetime = seq(fdatetime, ldatetime, by="hour")
          datetime = as.data.frame(substr(datetime, 1, 16))
          colnames(datetime) = c("datetime")

          ndata = merge(datetime, obs, by="datetime", all=T)

          if(nrow(data) == nrow(ndata)){
            data[is.na(data[, c(varnm)]), c(varnm)] = ndata[is.na(data[, c(varnm)]), c(varnm)]
          } else {
            cat(sprintf("     OBS: Number of rows are different: nearest = %s!\n", nstnid))
          }
        }

        nrcnt = nrcnt + 1

        if((sum(is.na(data[,c(varnm)])) == 0) | nrcnt == nrow(nstnnms)){

          # if there is sill missing after considering nearest 10 stations, fill the missing using interpolation
          if(sum(is.na(data[,c(varnm)])) <= (length(data[,c(varnm)]-2))) {data[,c(varnm)] = na.approx(data[,c(varnm)], na.rm=F)}

          break

        }

      }

    }


  }

  return(data)

}

#' @export
kma.asos.hourly.station.update <- function(EnvList) {

  obsdir = EnvList$obsdir
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  rawdir = paste(obsdir, "/asos_hourly_download", sep="")
  outdir = paste(rawdir, "/station", sep="")
  SetWorkingDir(outdir)

  # Check available stations
  #chk = read.csv(paste(rawdir, "/2014.csv", sep=""), header=T)
  #stnnms = unique(chk[, 1])
  #stncnt = length(stnnms)
  StnDFile = paste(stndir, "/", stnfile, sep="")
  stninfo = read.csv(StnDFile, header=T)
  stnnms = as.numeric(substr(matrix(stninfo$ID), 3, 5))
  stncnt = length(stnnms)

  # Extract staion data from yearly downloaded data
  for(i in 1:stncnt){
    stnid= stnnms[i]

    OutDFile = paste(outdir, "/", sprintf("ID%03d", stnid), ".csv", sep="")

    srchstr = paste("????.csv", sep="")
    flist = list.files(rawdir, pattern = glob2rx(srchstr), full.names = F)
    fcnt = length(flist)

    for(j in 1:fcnt){
      fname = flist[j]
      InDFile  = paste(rawdir, "/", flist[j], sep="")
      data = read.csv(InDFile, header=T)
      data = data[which(data[, 1] == stnid), c(2, 1, 5, 10, 3, 9, 7, 4, 15, 13, 12)]
      if(j == 1){
        out = data
      } else {
        tmp = data
        out = rbind(out, tmp)
      }

    }

    colnames(out) = c("datetime", "stnid", "wind", "pres", "atem",  "dewp", "rhum", "prec", "clou", "solr", "sshine")

    write.csv(out, OutDFile, row.names=F)

    # Interpolate within a day
    if(sum(is.na(out$wind)) <= (length(out$wind)-2)){out[, c("wind")] = na.approx(out[, c("wind")], na.rm=F, maxgap=24)}
    if(sum(is.na(out$pres)) <= (length(out$pres)-2)){out[, c("pres")] = na.approx(out[, c("pres")], na.rm=F, maxgap=24)}
    if(sum(is.na(out$atem)) <= (length(out$atem)-2)){out[, c("atem")] = na.approx(out[, c("atem")], na.rm=F, maxgap=24)}
    if(sum(is.na(out$dewp)) <= (length(out$dewp)-2)){out[, c("dewp")] = na.approx(out[, c("dewp")], na.rm=F, maxgap=24)}
    if(sum(is.na(out$rhum)) <= (length(out$rhum)-2)){out[, c("rhum")] = na.approx(out[, c("rhum")], na.rm=F, maxgap=24)}
    out$prec[is.na(out$prec)] = 0
    if(sum(is.na(out$clou)) <= (length(out$clou)-2)){out[, c("clou")] = na.approx(out[, c("clou")], na.rm=F, maxgap=24)}
    if(sum(is.na(out$solr)) <= (length(out$solr)-2)){out[which(as.numeric(substr(out$datetime,12, 13)) <= 5 | as.numeric(substr(out$datetime,12, 13)) >= 21), c("solr")] = 0}
    if(sum(is.na(out$solr)) <= (length(out$solr)-2)){out[, c("solr")] = na.approx(out[, c("solr")], na.rm=F, maxgap=13)}
    if(sum(is.na(out$sshine)) <= (length(out$sshine)-2)){out[which(as.numeric(substr(out$datetime,12, 13)) <= 5 | as.numeric(substr(out$datetime,12, 13)) >= 21), c("sshine")] = 0}
    if(sum(is.na(out$sshine)) <= (length(out$sshine)-2)){out[, c("sshine")] = na.approx(out[, c("sshine")], na.rm=F, maxgap=13)}

    if(any(is.na(out$wind))){
      out$wind[is.na(out$wind)] = -999
      #cat(sprintf("Station=%03d Variable=Wind Speed contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$pres))){
      out$pres[is.na(out$pres)] = -999
      #cat(sprintf("Station=%03d Variable=Pressure contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$atem))){
      out$atem[is.na(out$atem)] = -999
      #cat(sprintf("Station=%03d Variable=Temperature contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$dewp))){
      out$dewp[is.na(out$dewp)] = -999
      #cat(sprintf("Station=%03d Variable=dewp contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$rhum))){
      out$rhum[is.na(out$rhum)] = -999
      #cat(sprintf("Station=%03d Variable=rhum contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$clou))){
      out$clou[is.na(out$clou)] = -999
      #cat(sprintf("Station=%03d Variable=clou contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$solr))){
      out$solr[is.na(out$solr)] = -999
      #cat(sprintf("Station=%03d Variable=solr contains NA and replaced by -999!\n", stnid))
    }

    if(any(is.na(out$sshine))){
      out$sshine[is.na(out$sshine)] = -999
      #cat(sprintf("Station=%03d Variable=sshine contains NA and replaced by -999!\n", stnid))
    }


    OutDFile = paste(outdir, "/", sprintf("ID%03d", stnid), "-Hourly.csv", sep="")
    write.csv(out, OutDFile, row.names=F)

    cat(sprintf("KMA hourly data has been extracted: Station=%s\n", stnid))

  }

}


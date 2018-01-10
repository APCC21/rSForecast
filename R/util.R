#' @export
Set.Working.Environment <- function(envfile, override=list()) {

  options(stringsAsFactors=FALSE)

  data <- yaml::yaml.load_file(envfile)
  data <- lapply(data, function(x) if (is.character(x)) gsubfn::gsubfn("\\$\\((.*?)\\)", data, x) else x)

  if(data$stndir == "User") data$stndir = data$obsdir
  if(data$stndir == "GHCN") data$stndir = data$ghcndir

  if(!file.exists(data$prjdir)) dir.create(data$prjdir, showWarnings=F, recursive=T)
  if(!file.exists(data$mmedir)) dir.create(data$mmedir, showWarnings=F, recursive=T)
  if(!file.exists(data$vardir)) dir.create(data$vardir, showWarnings=F, recursive=T)
  if(!file.exists(data$idxdir)) dir.create(data$idxdir, showWarnings=F, recursive=T)
  if(!file.exists(data$rnldir)) dir.create(data$rnldir, showWarnings=F, recursive=T)
  if(!file.exists(data$bnddir)) dir.create(data$bnddir, showWarnings=F, recursive=T)
  if(!file.exists(data$aphrodir)) dir.create(data$aphrodir, showWarnings=F, recursive=T)
  if(!file.exists(data$stndir)) dir.create(data$stndir, showWarnings=F, recursive=T)
  if(!file.exists(data$smpldir)) dir.create(data$smpldir, showWarnings=F, recursive=T)

#  if (missing(envfile)) {
#    envfile = file.path(basedir, sprintf("%s.txt", prjname))
#  } else {
#    envfile = file.path(basedir, envfile)
#  }

  if(!file.exists(envfile)){
    stop(sprintf("Environment file %s does not exist!", envfile))
  }

#  tbl = read.table(envfile, header=F, sep="=", stringsAsFactors=FALSE)
#  colnames(tbl) = c("varnm", "varval")
#  for(i in 1:nrow(tbl)) {
#    assign(tbl$varnm[i], unlist(strsplit(tbl$varval[i], split=",")))
#  }

  outList = list("prjdir"=data$prjdir,
                 "dbdir"=data$dbdir,
                 "mmedir"=data$mmedir,
                 "bnddir"=data$bnddir,
                 "aphrodir"=data$aphrodir,
                 "vardir"=data$vardir,
                 "idxdir"=data$idxdir,
                 "rnldir"=data$rnldir,
                 "obsdir"=data$obsdir,
                 "asosdir"=data$asosdir,
                 "ghcndir"=data$ghcndir,
                 "stndir"=data$stndir,
                 "smpldir"=data$smpldir,
                 "mdlnms_3mon"=data$mdlnms_3mon,
                 "mdlnms_6mon"=data$mdlnms_6mon,
                 "ptrnms"=data$ptrnms,
                 "NBest"=data$NBest,
                 "cpcidxs"=data$cpcidxs,
                 "apccidxs"=data$apccidxs,
                 "nrange"=data$nrange,
                 "minlat"=data$minlat,
                 "maxlat"=data$maxlat,
                 "MinGrdCnt"=data$MinGrdCnt,
                 "MaxGrdCnt"=data$MaxGrdCnt,
                 "smonth"=data$smonth,
                 "emonth"=data$emonth,
                 "syear_obs"=data$syear_obs,
                 "eyear_obs"=data$eyear_obs,
                 "syear_mme"=data$syear_mme,
                 "eyear_mme"=data$eyear_mme,
                 "eyear_sim"=data$eyear_sim,
                 "stnfile"=data$stnfile,
                 "varfile"=data$varfile,
                 "idxfile"=data$idxfile,
                 "ptrfile"=data$ptrfile,
                 "bndfile"=data$bndfile,
                 "nrstfile"=data$nrstfile,
                 "BCPointOpt"=data$BCPointOpt,
                 "BCAreaOpt"=data$BCAreaOpt,
                 "CIRegOpt"=data$CIRegOpt,
                 "MWRegOpt"=data$MWRegOpt,
                 "MWRObsOpt"=data$MWRObsOpt,
                 "tscale" =data$tscale,
                 "fcstmode" =data$fcstmode,
                 "combnmode" =data$combnmode,
                 "fiyearmode" =data$fiyearmode,
                 "AcuMonths"=data$AcuMonths,
                 "precopt"=data$precopt,
                 "CRAdj"=data$CRAdj)

  # override
  for (varname in names(override)) {
    outList[[varname]] = override[[varname]]
  }

  return(outList)
}

#' @export
Update.Climate.Information <- function(EnvList, updatemode) {

  if(EnvList$CIRegOpt == "On"){
    CIReg.Update.Climate.Index (EnvList)
  }
  if(EnvList$MWRegOpt == "On" | EnvList$BCPointOpt == "On" | EnvList$BCAreaOpt == "On"){
    Update.APCC.forecast.data (EnvList, updatemode)
  }
  if(EnvList$MWRObsOpt == "On"){
    MWRObs.Download.Reanalysis1 (EnvList)
  }

}

#' @export
Daily.Hourly.Sampling <- function(EnvList, fiyearmon, smplmoncnt) {

  tscale = EnvList$tscale
  fcstmode = EnvList$fcstmode
  combnmode = as.logical(EnvList$combnmode)
  fiyearmode = as.logical(EnvList$fiyearmode)

  # Update Integrated Table, Monthly Summary Tables, Error statistics
  RTFcst.Create.Model.Integration.Table (EnvList, AcuMonths=1)
  RTFcst.Create.Summary.Table(EnvList, precopt=F)
  RTFcst.Calculate.Error.Statistics(EnvList)

  Decide.Range.Values(EnvList)
  RTFcst.Mahalanobis.Sampling.Observed.9var(EnvList, fiyearmon, smplmoncnt)
  RTFcst.Create.Summary.Graph.Table(EnvList, fiyearmon)

}

#' @export
Run.RealTime.Forecast.Model <- function(EnvList, fiyearmon) {

  BCPointOpt = EnvList$BCPointOpt
  BCAreaOpt = EnvList$BCAreaOpt
  CIRegOpt = EnvList$CIRegOpt
  MWRegOpt = EnvList$MWRegOpt
  MWRObsOpt = EnvList$MWRObsOpt
  prjdir = EnvList$prjdir
  ltcnt = 6

  syear = as.numeric(substr(fiyearmon, 1, 4))
  eyear = as.numeric(substr(seq(as.Date(paste(fiyearmon, "-01", sep="")), by = "month", length = (ltcnt+1))[ltcnt + 1],1,4))

  # Delete existing folders
  for (yr in syear:(as.numeric(substr(Sys.Date(), 1, 4))+1)){
    if(CIRegOpt == "On"){
      cirdir = paste(prjdir, "CIReg", yr, sep="/")
      unlink(cirdir, recursive=T)
    }

    if(MWRegOpt == "On"){
      for(mmetype in c("3MON", "6MON")) {
        mwrdir = paste(prjdir, "MWReg", mmetype, yr, sep="/")
        unlink(mwrdir, recursive=T)
      }
    }

    if(MWRObsOpt == "On"){
      mwobsdir = paste(prjdir, "MWRObs", yr, sep="/")
      unlink(mwobsdir, recursive=T)
    }
  }


  if(BCPointOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCPoint.Extract.MME.Stations(EnvList, mmetype, varnms = c("prec", "t2m"))
      BCPoint.Bias.Correction (EnvList, mmetype, varnms = c("prec", "t2m"))
      BCPoint.Create.Summary.Table(EnvList, mmetype)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="BCPoint", mmetype)
    }
  }

  if(BCAreaOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCArea.Calculate.SME.Area.Averages(EnvList, mmetype, varnms=c("prec"), pointopt=F)
      BCArea.SME.Bias.Correction(EnvList, mmetype)
      BCArea.Create.Summary.Table(EnvList, mmetype, AcuMonths=1, precopt=F, CRAdj=1.0)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="BCArea", mmetype)
    }
  }

  if(CIRegOpt == "On"){
    for(tyear in syear:eyear){
      CIReg.Copy.Hindcast.Results(EnvList, tyear)
      CIReg.Select.BestModel.and.Run.Best.Regression(EnvList, tyear)
      CIReg.Create.Summary.Table(EnvList)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="CIReg")
    }
  }

  if(MWRegOpt == "On"){
    for(tyear in syear:eyear){
      for(mmetype in c("3MON", "6MON")) {
        MWReg.Copy.Hindcast.Results(EnvList, mmetype, tyear)
        MWReg.Extract.Best.Predictor.TSeries(EnvList, mmetype, tyear)
        MWReg.Run.Best.Regression(EnvList, mmetype, tyear)
        MWReg.Draw.Best.Predictor.Map(EnvList, mmetype, tyear)
        MWReg.Create.Summary.Table(EnvList, mmetype)
        Calculate.Error.Statistics.Selected(EnvList, dsmethod="MWReg", mmetype)
      }
    }
  }

  if(MWRObsOpt == "On"){
    for(tyear in syear:eyear){
      MWRObs.Copy.Hindcast.Results (EnvList, tyear)
      MWRObs.Extract.Best.Predictor.TSeries (EnvList, tyear)
      MWRObs.Run.Best.Regression (EnvList, tyear)
      MWRObs.Draw.Best.Predictor.Map (EnvList, tyear)
      MWRObs.Create.Summary.Table (EnvList)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="MWRObs")
    }
  }

  # Create Table for Model Selection and Create Summary Table
  RTFcst.Create.Model.Integration.Table (EnvList, AcuMonths=1)
  RTFcst.Create.Summary.Table(EnvList, precopt=F)
  RTFcst.Calculate.Error.Statistics(EnvList)

}

#' @export
Integrate.Individual.Forecast.Model <- function(EnvList) {

  BCPointOpt = EnvList$BCPointOpt
  BCAreaOpt = EnvList$BCAreaOpt
  CIRegOpt = EnvList$CIRegOpt
  MWRegOpt = EnvList$MWRegOpt
  MWRObsOpt = EnvList$MWRObsOpt

  if(BCPointOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCPoint.Create.Summary.Table(EnvList, mmetype)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="BCPoint", mmetype)
    }
  }

  if(BCAreaOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCArea.Create.Summary.Table(EnvList, mmetype, AcuMonths=1, precopt=F, CRAdj=1.0)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="BCArea", mmetype)

    }
  }

  if(CIRegOpt == "On"){
    CIReg.Create.Summary.Table(EnvList)
    Calculate.Error.Statistics.Selected(EnvList, dsmethod="CIReg")
  }

  if(MWRObsOpt == "On"){
    MWRObs.Create.Summary.Table (EnvList)
    Calculate.Error.Statistics.Selected(EnvList, dsmethod="MWRObs")
  }

  if(MWRegOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      MWReg.Create.Summary.Table(EnvList, mmetype)
      Calculate.Error.Statistics.Selected(EnvList, dsmethod="MWReg", mmetype)
    }
  }


}

#' @export
Hindcast.Forecast.Model.Construction <- function(EnvList) {

  BCPointOpt = EnvList$BCPointOpt
  BCAreaOpt = EnvList$BCAreaOpt
  CIRegOpt = EnvList$CIRegOpt
  MWRegOpt = EnvList$MWRegOpt
  MWRObsOpt = EnvList$MWRObsOpt
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)

  if(BCPointOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCPoint.Extract.MME.Stations(EnvList, mmetype, varnms = c("prec", "t2m"))
      BCPoint.Bias.Correction (EnvList, mmetype, varnms = c("prec", "t2m"))
    }
  }

  if(BCAreaOpt == "On"){
    for(mmetype in c("3MON", "6MON")) {
      BCArea.Calculate.SME.Area.Averages(EnvList, mmetype, varnms=c("prec"), pointopt=F)
      BCArea.SME.Bias.Correction(EnvList, mmetype)
    }
  }

  if(CIRegOpt == "On"){
    for(tyear in syear_mme:eyear_obs){
      CIReg.Calculate.All.Regression.Split.Validation (EnvList, tyear)
      CIReg.Select.BestModel.and.Run.Best.Regression(EnvList, tyear)
    }
  }

  if(MWRObsOpt == "On"){
    for(tyear in syear_mme:eyear_obs){
      MWRObs.Calculate.All.Regression (EnvList, tyear)
      MWRObs.Extract.Best.Predictor.TSeries (EnvList, tyear)
      MWRObs.Run.Best.Regression (EnvList, tyear)
      MWRObs.Draw.Best.Predictor.Map (EnvList, tyear)
    }
  }

  if(MWRegOpt == "On"){
    for(tyear in syear_mme:eyear_obs){
      for(mmetype in c("3MON", "6MON")) {
        MWReg.Calculate.All.Regression(EnvList, mmetype, tyear)
        MWReg.Extract.Best.Predictor.TSeries(EnvList, mmetype, tyear)
        MWReg.Run.Best.Regression(EnvList, mmetype, tyear)
        MWReg.Draw.Best.Predictor.Map(EnvList, mmetype, tyear)
      }
    }
  }

}

#' @export
ReadReanalysis1 <- function(wdir, ncfile) {

  #ncfile = "F:/2016-SDS_Training/SForecast/Database/reanalysis1/omega.mon.mean.nc"

  #setwd(wdir)

  ncfile = paste(wdir, "/", ncfile, sep="")

  fin = nc_open(ncfile)

  # Get the file information into array
  finfo = capture.output(print(fin))

  # Get the dimension, variable name, dimension names
  nrow=length(finfo[])
  for(i in 1:nrow){
    str = finfo[i]
    str = str_split_fixed(str,"float ",2)[2]
    if(!(str == "")){
      varnm = str_split_fixed(str, fixed("["), 2)[1]
      str = str_split_fixed(str, fixed("["), 2)[2]
      str = str_split_fixed(str, fixed("]"), 2)[1]
      dimsz = str_count(str,fixed(",")) + 1
      dimnms = str_split_fixed(str, fixed(","), dimsz)
    }
  }

  lon = ncvar_get(fin, dimnms[1])
  lat = ncvar_get(fin, dimnms[2])
  if(dimsz == 3) {
    lvl = NA
    tim = ncvar_get(fin, dimnms[3])
  }
  if(dimsz == 4) {
    lvl = ncvar_get(fin, dimnms[3])
    tim = ncvar_get(fin, dimnms[4])
  }

  variable = ncvar_get(fin, varnm)

  timestr = ncatt_get(fin, "time", "units")$value
  tstr = str_split_fixed(timestr, fixed(" "), 4)[3]
  tunit = str_split_fixed(timestr, fixed(" "), 4)[1]
  tyear = as.integer(str_split_fixed(tstr, fixed("-"), 3)[1])
  tmonth = as.integer(str_split_fixed(tstr, fixed("-"), 3)[2])
  tday = as.integer(str_split_fixed(tstr, fixed("-"), 3)[3])

  # http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm
  if(tunit == "hours") tim = as.Date(chron::chron(tim/24, origin = c(tmonth, tday, tyear)))
  if(tunit == "days") tim = as.Date(chron::chron(tim, origin = c(tmonth, tday, tyear)))

  outList = list("x"=lon, "y"=lat, "l"=lvl, "t"=tim, "var"=variable, "fin"=finfo)

  nc_close(fin)

  return(outList)

}

#' @export
ReadNetCDF4 <- function(wdir, ncfile) {


  setwd(wdir)

  fin = ncdf4::nc_open(ncfile)

  # Get the file information into array
  finfo = capture.output(print(fin))

  # Get the dimension, variable name, dimension names
  nrow=length(finfo[])
  for(i in 1:nrow){
    str = finfo[i]
    str = str_split_fixed(str,"float ",2)[2]
    if(!(str == "")){
      varnm = str_split_fixed(str, fixed("["), 2)[1]
      str = str_split_fixed(str, fixed("["), 2)[2]
      str = str_split_fixed(str, fixed("]"), 2)[1]
      dimsz = str_count(str,fixed(",")) + 1
      dimnms = str_split_fixed(str, fixed(","), dimsz)
    }
  }

  lon = ncvar_get(fin, dimnms[1])
  lat = ncvar_get(fin, dimnms[2])
  variable = ncvar_get(fin,varnm)

  # In the case of multiple layers, get the first layer (surface ) values
  if(varnm == 'hur'){
    variable = variable[,,1,]
  }

  outList = list("x"=lon, "y"=lat, "var"=variable, "fin"=finfo)

  nc_close(fin)

  return(outList)

}

#' @export
SetWorkingDir <- function(wdir) {

  dir.create(wdir, showWarnings=F,recursive=T)
  setwd(wdir)

}

#' @export
APCCMME.Get.MMEfile.Info <- function(mmedir){

  #mmedir = "E:/SForecast-Test/Database/apcc-mme/3MON/CWB/APR"

  dlist = list.dirs(mmedir)

  sydir = dlist[2]
  srchstr = paste("*.nc", sep="")
  flist = list.files(sydir, pattern = glob2rx(srchstr), full.names = T)
  ncfile = flist[1]
  fin = nc_open(ncfile)
  finfo = capture.output(print(fin))

  nrow=length(finfo[])
  for(i in 1:nrow){
    str = finfo[i]

    str1 = str_split_fixed(str,"float ",2)[2]
    if(!(str1 == "")){
      varnm = str_split_fixed(str1, fixed("["), 2)[1]
      str1 = str_split_fixed(str1, fixed("["), 2)[2]
      str1 = str_split_fixed(str1, fixed("]"), 2)[1]
      dimsz = str_count(str1,fixed(",")) + 1
      dimnms = str_split_fixed(str1, fixed(","), dimsz)
    }

    str2 = str_split_fixed(str,"time  Size:",2)[2]
    if(!(str2 == "")){
      MaxLTime = as.numeric(substr(str2, 1, 3))
    }

    str3 = str_split_fixed(str,"level  Size:",2)[2]
    if(!(str3 == "")){
      esmcnt = as.numeric(substr(str3, 1, 3))
    }
    str4 = str_split_fixed(str,"lat  Size:",2)[2]
    if(!(str4 == "")){
      latcnt = as.numeric(substr(str4, 1, 3))
    }
    str5 = str_split_fixed(str,"lon  Size:",2)[2]
    if(!(str5 == "")){
      loncnt = as.numeric(substr(str5, 1, 3))
    }

  }

  outList = list("MaxLTime"=MaxLTime, "esmcnt"=esmcnt, "latcnt"=latcnt, "loncnt"=loncnt, "d"=d)

  nc_close(fin)

  return(outList)

}

#' @export
Decide.Range.Values <- function(EnvList){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  #smonth = as.numeric(EnvList$smonth)
  #emonth = as.numeric(EnvList$emonth)
  nrange = as.numeric(EnvList$nrange)

  #smonth = 1
  #emonth = 12
  #nrange = 3
  outdir = paste(prjdir, "/0_RTForecast/ProbabilityRanges", sep="")
  SetWorkingDir(outdir)

  # Read variable files
  VarDFile = paste(vardir, "/", varfile, sep="")
  data = read.csv(VarDFile, header=T)

  varnms = colnames(data)[2:ncol(data)]
  varcnt = length(varnms)

  for(i in 1:varcnt){
    varnm = varnms[i]

    vardata = data[c("yearmon", varnm)]
    if(varnm == "prec"){
      vardata = prcp.mmday2mmmon(vardata)
    }

    for(j in 1:12){

      curmon = j
      mondata = vardata[which(as.numeric(substr(vardata$yearmon, 6, 7)) == curmon), c(varnm)]

      if(nrange == 3){
        quant = as.data.frame(quantile(mondata, c(0.33, 0.67)))
        curmon = as.data.frame(curmon)
        cvals = cbind(curmon, t(quant))
        colnames(cvals) = c("month", "1st", "2nd")
      } else if(nrange == 5) {
        quant = as.data.frame(quantile(mondata, c(0.2, 0.4, 0.6, 0.8)))
        curmon = as.data.frame(curmon)
        cvals = cbind(curmon, t(quant))
        colnames(cvals) = c("month", "1st", "2nd", "3rd", "4th")
      }

      if(j == 1){
        out = cvals
      } else {
        tmp = cvals
        out = rbind(out, tmp)
      }

    }


    OutDFile = paste(outdir, "/", varnm, "-probability_range.csv", sep="")
    write.csv(out, OutDFile, row.names=F)

  }

}

#' @export
Sampling.Daily.Data.9var <- function(outdir, sampledir, stninfo, mdlnm) {

  #outdir = "F:/SForecast-RT/Korea/0_RTForecast/6MON/Daily_Sampling/2010-01"
  #outdir = scndir
  #sampledir = dsmpldir

  stnnms = stninfo$ID
  stncnt = length(stnnms)
  fname = paste(outdir, "/", mdlnm, "/", stnnms[stncnt], ".csv", sep="")

  if(!file.exists(fname))  {

    # read station_average for MME data
    bestfile = paste(outdir, "/BestMon-", mdlnm, ".csv", sep="")

    bestarray = read.csv(bestfile, header=T)

    if(any(is.na(bestarray$bestmon))){
      cat(sprintf("     Util: Daily sampling has been skipped: Model=%s\n", mdlnm))
    } else {

      firstdate = as.Date(paste(bestarray$yearmon[1], "-1", sep=""))
      lastimsi = as.Date(paste(bestarray$yearmon[nrow(bestarray)], "-1", sep="")) #First day of the last month
      lastimsi = seq(lastimsi,length=2,by="months")-1   #End days of the last two months
      lastdate = lastimsi[2]
      #lastdate = as.Date(paste(bestarray$yearmon[nrow(bestarray)], "-31", sep=""))
      date = seq(firstdate, lastdate, by="day")

      moncnt = nrow(bestarray)
      for(ii in 1:stncnt){

        stnid = stnnms[ii]
        # Lat and Elev info for calculating sunshine hours
        lat = stninfo$Lat[ii]
        elev = stninfo$Elev[ii]

        #sampledir ="F:/SForecast-RT/Database/obs_samples/kma-asos_9var"
        srchstr = paste("*", stnid, "*.csv", sep="")
        SmplDFile = list.files(sampledir, pattern = glob2rx(srchstr), full.names = T)
        obs = read.csv(SmplDFile, header=T)

        #obs = Get.Station.Data.9var(stnid, obsdir)
        if(ncol(obs) == 12){
          colnames(obs) = c("year", "month", "day", "prcp", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")
        }
        if(ncol(obs) == 9){
          colnames(obs) = c("year", "month", "day", "prcp", "tmax", "tmin", "wspd", "rhum", "rsds")
        }

        for(jj in 1:moncnt){

          if(!is.na(bestarray$bestmon[jj])){
            byear = as.numeric(substr(bestarray$bestmon[jj], 1, 4))
            bmon = as.numeric(substr(bestarray$bestmon[jj], 6, 7))
            pratio = bestarray$prec[jj]/bestarray$prec_obs[jj]
            tdiff = bestarray$t2m[jj]-bestarray$t2m_obs[jj]

            # Select daily data for the best month from the observed time series
            bestobs = obs[which(obs$year == byear & obs$month == bmon),]
            bestobs[, c("prcp")]= bestobs[, c("prcp")] * pratio
            bestobs[, c("tmax")] = bestobs[, c("tmax")] + tdiff
            bestobs[, c("tmin")] = bestobs[, c("tmin")] + tdiff
            if(ncol(obs) == 12){ bestobs[, c("tavg")] = bestobs[, c("tavg")] + tdiff }
            bestobs = bestobs[, 4:ncol(obs)]

            # Check the number of dates and adjust for the February
            if(bmon == 2){
              nbestday = nrow(bestobs)
              nesmday = numberOfDays(as.Date(paste(bestarray$yearmon[jj], "-1", sep="")))
              if(nbestday > nesmday){ bestobs = bestobs[1:nesmday,] }
              if(nbestday < nesmday){
                dummy = bestobs[nbestday,]
                bestobs = rbind(bestobs, dummy)
              }
            }

          } else { ### If bestmon is NA, use observed data
            byear = as.numeric(substr(bestarray$yearmon[jj], 1, 4))
            bmon = as.numeric(substr(bestarray$yearmon[jj], 6, 7))

            # Select daily data for the best month from the observed time series
            bestobs = obs[which(obs$year == byear & obs$month == bmon),]
            bestobs = bestobs[,4:ncol(obs)]
          }
          if(jj == 1){
            out = bestobs
          } else {
            tmp = bestobs
            out = rbind(out, tmp)
          }

        }

        outdata = cbind(date, out)
        if(ncol(obs) == 12){
          colnames(outdata) = c("date", "prcp", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")
        }
        if(ncol(obs) == 9){
          colnames(outdata) = c("date", "prcp", "tmax", "tmin", "wspd", "rhum", "rsds")
        }


        # Fill solar radiation
        if(!any(is.na(outdata$tmin)) & !any(is.na(outdata$tmax))){
          outdata$tmin[outdata$tmin > outdata$tmax] = outdata$tmax[outdata$tmin > outdata$tmax]
        }

        jday = strptime(outdata$date, "%Y-%m-%d")$yday+1

        isNA = which(is.na(outdata[,"rsds"]) & !is.na(outdata[,"tmax"]) & !is.na(outdata[,"tmin"]))
        if(length(isNA[]) >= 1) {
          outdata[isNA,"rsds"] = Solar(lat,jday[isNA], outdata[isNA,"tmax"], outdata[isNA,"tmin"], latUnits='degrees')/1000.0
          # Unit conversion kj/m2/day --> MJ/m2 (Default unit of Solar() is Kj/m2)
        }

        wdir = paste(outdir, "/", mdlnm, sep="")
        SetWorkingDir(wdir)
        outfile = paste(stnid, ".csv", sep="")
        write.csv(outdata, outfile, row.names=F)


      }

      cat(sprintf("     Util: Daily or Hourly data has been created: Model=%s\n", mdlnm))

    } # End If


  }

  return(fname)

}

#' @export
Sampling.Hourly.Data.9var <- function(outdir, sampledir, stninfo, mdlnm) {

  #outdir = scndir
  #mdlnm = "MIN"
  #sampledir = "F:/SForecast-RT/Database/obs_samples/kma-asos-hourly"

  stnnms = stninfo$ID
  stncnt = length(stnnms)
  fname = paste(outdir, "/", mdlnm, "/", stnnms[stncnt], ".csv", sep="")

  if(!file.exists(fname))  {

    # read station_average for MME data
    bestfile = paste(outdir, "/BestMon-", mdlnm, ".csv", sep="")

    bestarray = read.csv(bestfile, header=T)

    #hourlydata = read.csv("F:/SForecast-RT/Database/kma-asos-hourly/hourly-data.csv", header=T)

    if(any(is.na(bestarray$bestmon))){
      cat(sprintf("     Util: Hourly sampling has been skipped: Model=%s\n", mdlnm))
    } else {

      firstdate = as.POSIXct(paste(bestarray$yearmon[1], "-1 00:00", sep=""))
      lastimsi = as.Date(paste(bestarray$yearmon[nrow(bestarray)], "-1", sep="")) #First day of the last month
      lastimsi = seq(lastimsi,length=2,by="months")-1   #End days of the last two months
      lastdate = as.POSIXct(paste(lastimsi[2], "23:00", sep=""))
      date = seq(firstdate, lastdate, by="hour")

      moncnt = nrow(bestarray)
      for(ii in 1:stncnt){
        stnnm = stnnms[ii]
        stnid = as.numeric(substr(stnnm,3,5))
        # Lat and Elev info for calculating sunshine hours
        lat = stninfo$Lat[ii]
        elev = stninfo$Elev[ii]

        #obs = Get.Station.Data.Hourly(hourlydata, stnid)
        srchstr = paste("*", stnid, "*.csv", sep="")
        SmplDFile = list.files(sampledir, pattern = glob2rx(srchstr), full.names = T)
        obs = read.csv(SmplDFile, header=T)

        #seq(from = as.POSIXct("2012-05-15 00:00"), to = as.POSIXct("2012-05-17 23:00"), by = "hour")
        #colnames(obs) = c("year", "month", "yearmon", "date", "prcp", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")

        for(jj in 1:moncnt){

          if(!is.na(bestarray$bestmon[jj])){
            byear = as.numeric(substr(bestarray$bestmon[jj], 1, 4))
            bmon = as.numeric(substr(bestarray$bestmon[jj], 6, 7))
            pratio = bestarray$prec[jj]/bestarray$prec_obs[jj]
            tdiff = bestarray$t2m[jj]-bestarray$t2m_obs[jj]

            # Select daily data for the best month from the observed time series
            bestobs = obs[which(as.numeric(substr(obs$TM, 1, 4)) == byear & as.numeric(substr(obs$TM, 6, 7)) == bmon),]
            bestobs = unique(bestobs)
            bestobs[, c("P.m")]= bestobs[, c("P.m")] * pratio
            bestobs[, c("T.m")] = bestobs[, c("T.m")] + tdiff

            # Check the number of dates and adjust for the February
            if(bmon == 2){
              nbestday = as.numeric(substr(bestobs[nrow(bestobs), c("TM")],9,10))
              nesmday = numberOfDays(as.Date(paste(bestarray$yearmon[jj], "-1", sep="")))
              if(nbestday > nesmday){ bestobs = bestobs[1:(nesmday*24),] }
              if(nbestday < nesmday){
                dummy = bestobs[(nrow(bestobs)-23):nrow(bestobs),]
                bestobs = rbind(bestobs, dummy)
              }
            }

            bestobs = bestobs[,2:ncol(bestobs)]


          } else { ### If bestmon is NA, use observed data
            byear = as.numeric(substr(bestarray$yearmon[jj], 1, 4))
            bmon = as.numeric(substr(bestarray$yearmon[jj], 6, 7))

            # Select daily data for the best month from the observed time series
            bestobs = obs[which(as.numeric(substr(obs$TM, 1, 4)) == byear & as.numeric(substr(obs$TM, 5, 6)) == bmon),]
            bestobs = bestobs[,2:ncol(bestobs)]

          }
          if(jj == 1){
            out = bestobs
          } else {
            tmp = bestobs
            out = rbind(out, tmp)
          }

        }

        #write.csv(out, "F:/SForecast-RT/Korea/0_RTForecast/3MON/Hourly_Sampling/2012-07/test.csv", row.names=F)

        outdata = cbind(date, out)
        colnames(outdata) = c("TM", "STNID", "WIND", "PRES", "T.m", "Tdew.m", "RH", "P.m", "CC", "RAD")

        wdir = paste(outdir, "/", mdlnm, sep="")
        SetWorkingDir(wdir)
        outfile = paste(stnnm, ".csv", sep="")
        write.csv(outdata, outfile, row.names=F)


      }

      cat(sprintf("     Util: Daily data has been created: Model=%s\n", mdlnm))

    } # End If


  }

  return(fname)

}

#' @export
Select.BestFit.Month <- function(wdir, meanfile) {

  #wdir = "F:/SForecast-Dev/Korea/0_RTForecast/Daily_Sampling"
  #meanfile = "StnMean-All-Simulations.csv"

  # read station_average for MME data
  MeanDFile = paste(wdir, "/", meanfile, sep="")

  BestDFile = paste(wdir, "/BestMon-All-Simulations.csv", sep="")

  # read station_average for MME data
  meanarray = read.csv(MeanDFile, header=T)
  mcolnms = colnames(meanarray)
  moncnt = nrow(meanarray)
  colnum = ncol(meanarray)
  colcnt = colnum + 3
  bestarray = as.data.frame(array(NA, dim=c(moncnt, colcnt)))  #Empty array
  colnames(bestarray) = c("yearmon", "model", "LTime", "ESM", "prec", "t2m", "bestmon", "prec_obs", "t2m_obs")
  colnames(bestarray) = c(mcolnms, "bestmon", "prec_obs", "t2m_obs")

  for(i in 1:moncnt){

    bestarray[i, 1:colnum] = meanarray[i, 1:colnum]

    if(is.na(meanarray[i,c("prec")]) | is.na(meanarray[i,c("t2m")])){
      bestarray[i,c("bestmon", "prec_obs", "t2m_obs")] = NA
    } else {
      curmon = as.numeric(substr(meanarray[i,1], 6, 7))
      curdata = t(as.vector(meanarray[i,c("prec", "t2m")]))

      obsname = paste(wdir, "/ObsMean-", sprintf("M%02d", curmon), ".csv", sep="" )
      obsarray = read.csv(obsname, header=T)
      obsdata = as.matrix(obsarray[,c("prec", "temp")])

      covar = cov(obsarray[,c("prec", "temp")])

      MD = mahalanobis(obsarray[,c("prec", "temp")], curdata, covar)

      bestrnum = which(MD == min(MD))
      bestarray[i,c("bestmon")] = as.character(obsarray[bestrnum,c("yearmon")])
      bestarray[i,c("prec_obs", "t2m_obs")] = obsarray[bestrnum,c("prec", "temp")]
    }

  }

  write.csv(bestarray, BestDFile, row.names=F)

}

#' @export
numberOfDays <- function(date) {
  m <- format(date, format="%m")

  while (format(date, format="%m") == m) {
    date <- date + 1
  }

  return(as.integer(format(date - 1, format="%d")))
}

#' @export
prcp.mmday2mmmon <- function(data) {

  moncnt = nrow(data)
  for(i in 1:moncnt){
    daycnt  = numberOfDays(as.Date(paste(data$yearmon[i], "-1", sep="")))
    data[i,2:ncol(data)] = data[i,2:ncol(data)] * daycnt
  }
  return(data)

}

#' @export
prcp.mmmon2mmday <- function(data) {

  moncnt = nrow(data)
  for(i in 1:moncnt){
    daycnt  = numberOfDays(as.Date(paste(data$yearmon[i], "-1", sep="")))
    data[i,2:ncol(data)] = data[i,2:ncol(data)] / daycnt
  }
  return(data)

}

#' @export
critical.r <- function( n, alpha = 0.05 ) {
  df <- n - 2
  critical.t <- qt( alpha/2, df, lower.tail = F )
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return( critical.r )
}

#' @export
Extract.Point.Value <- function(xin, yin, ncfile, wdir) {

  # get x and y information
  nc = ReadNetCDF4(wdir, ncfile)
  x = nc$x
  # In the case, lon is west side (- value) convert it into 0-360 ranges
  if(x < 0){x = x + 360}
  y = nc$y
  var = nc$var

  ncols = length(x)
  nrows = length(y)

  xcoord = array(0, dim = c(ncols, 2))
  ycoord = array(0, dim = c(nrows, 2))

  for(i in 1:ncols){
    if(i == ncols){
      xcoord[i,1] = x[i] - (x[i]-x[i-1])/2.0
      xcoord[i,2] = x[i] + (x[i]-x[i-1])/2.0
    }else{
      xcoord[i,1] = x[i] - (x[i+1]-x[i])/2.0
      xcoord[i,2] = x[i] + (x[i+1]-x[i])/2.0
    }
  }

  for(i in 1:nrows){
    if(i == nrows){
      ycoord[i,1] = y[i] - (y[i]-y[i-1])/2.0
      ycoord[i,2] = y[i] + (y[i]-y[i-1])/2.0
    }else{
      ycoord[i,1] = y[i] - (y[i+1]-y[i])/2.0
      ycoord[i,2] = y[i] + (y[i+1]-y[i])/2.0
    }
  }

  colnum = which(xin >= xcoord[,1] & xin < xcoord[,2])
  rownum = which(yin >= ycoord[,1] & yin < ycoord[,2])

  if(length(dim(var)) == 3){
    var = var[colnum, rownum, ]
  } else if (length(dim(var)) == 4){
    var = var[colnum, rownum, , ]
  } else {
    print("Dimensions of the variable is greather than 4!")
  }

  outList = list("rown"=rownum, "coln"=colnum, "var"=var)

  return(outList)


}

#' @export
GIS.ncVar2raster <- function(x, y, var, ext=NULL) {


  rvar = t(var)[ncol(var):1,]
  r = raster(rvar)
  prj = "+proj=longlat +datum=WGS84 +no_defs"
  projection(r) = prj
  xdiff = (x[2]-x[1])/2.0
  tdiff = (y[length(y[])]-y[(length(y[]) - 1)])/2.0
  bdiff = (y[2]-y[1])/2.0
  xmin(r) = min(x) - xdiff
  xmax(r) = max(x) + xdiff
  ymin(r) = min(y) - bdiff
  ymax(r) = max(y) + tdiff

  if(!missing(ext)){
    r = crop(r, ext, snap='out')
  }

  return(r)

}

#' @export
Calculate.Probability <- function(nrange, vals, rngvals){

  tot = length(vals)
  if(nrange == 3){
    rng1 = rngvals[1]; rng2 = rngvals[2]
    p1 = length(which(vals < rng1))/tot * 100
    p2 = length(which(vals >= rng1 & vals < rng2))/tot * 100
    p3 = length(which(vals >= rng2))/tot * 100

    prob = cbind(p1, p2, p3)

  }

  if(nrange == 4){
    rng1 = rngvals[1]; rng2 = rngvals[2]; rng3 = rngvals[3]
    p1 = length(which(vals < rng1))/tot * 100
    p2 = length(which(vals >= rng1 & vals < rng2))/tot * 100
    p3 = length(which(vals >= rng2 & vals < rng3))/tot * 100
    p4 = length(which(vals >= rng3))/tot * 100

    prob = cbind(p1, p2, p3, p4)

  }

  if(nrange == 5) {
    rng1 = rngvals[1]; rng2 = rngvals[2]; rng3 = rngvals[3]; rng4 = rngvals[4]
    p1 = length(which(vals < rng1))/tot * 100
    p2 = length(which(vals >= rng1 & vals < rng2))/tot * 100
    p3 = length(which(vals >= rng2 & vals < rng3))/tot * 100
    p4 = length(which(vals >= rng3 & vals < rng4))/tot * 100
    p5 = length(which(vals >= rng4))/tot * 100

    prob = cbind(p1, p2, p3, p4, p5)
  }

  return(prob)

}

#' @export
Create.Contingency.Table <- function(rngvals, ctdata){

  nrange = length(rngvals) + 1
  maxval = max(ctdata)
  rngvals[nrange] = maxval*2
  dcnt = nrow(ctdata)

  CTable = array(0, dim=c(nrange, nrange))  #Empty array
  colnames(CTable) = paste("Obs", seq(1,nrange,1), sep="")
  rownames(CTable) = paste("Sim", seq(1,nrange,1), sep="")
  for(k in 1:dcnt){
    sim = ctdata[k,1]
    obs = ctdata[k,2]

    for(ii in 1:nrange){

      if(sim < rngvals[ii]) {
        for(jj in 1:nrange){

          if(obs < rngvals[jj]){
            CTable[ii, jj] = CTable[ii, jj] + 1
            break
          }

        } # OBS Loop
        break
      }

    } # sim Loop

  }

  return(CTable)

}

#' @export
Calculate.Error.Statistics.Selected <- function(EnvList, dsmethod, mmetype=NULL) {

  prjdir = EnvList$prjdir

  if(dsmethod == "BCArea" | dsmethod == "BCPoint" | dsmethod == "MWReg"){
    indir = paste(prjdir, "/0_Analysis/", dsmethod, "/", mmetype, "/Monthly-TSeries", sep="")
    outdir = paste(prjdir, "/0_Analysis/", dsmethod, "/", mmetype, sep="")
  } else {
    indir = paste(prjdir, "/0_Analysis/", dsmethod, "/Monthly-TSeries", sep="")
    outdir = paste(prjdir, "/0_Analysis/", dsmethod, sep="")

  }

  if(dsmethod == "BCArea"){
    flist = list.files(indir, pattern = glob2rx("BCArea*.csv"), full.names = F)
  } else if(dsmethod == "BCPoint") {
    flist = list.files(indir, pattern = glob2rx("BCPoint*.csv"), full.names = F)
  } else if(dsmethod == "MWReg") {
    flist = list.files(indir, pattern = glob2rx("MWReg*.csv"), full.names = F)
  } else if(dsmethod == "CIReg") {
    flist = list.files(indir, pattern = glob2rx("CIReg*.csv"), full.names = F)
  } else if(dsmethod == "MWRObs") {
    flist = list.files(indir, pattern = glob2rx("MWRObs*.csv"), full.names = F)
  }

  if(length(flist) > 0){
    varnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,2])
    monnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,3])
    monnms =  monnms[order(monnms)]
    acunms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,4])
    acunms = substr(acunms, 1, nchar(acunms)-4)

    varcnt = length(varnms); acucnt = length(acunms)

    for(i in 1:varcnt){
      varnm  = varnms[i]

      for(j in 1:acucnt){
        acunm = acunms[j]

        #for(k in 1:ltcnt){
        #  ltnm = ltnms[k]
        cnt = 1
        for(ii in 1:12){

          #monnm = monnms[ii]
          monnm = sprintf("M%02d", ii)

          if(dsmethod == "BCArea"){
            srchstr = paste("BCArea-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          } else if(dsmethod == "BCPoint") {
            srchstr = paste("BCPoint-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          } else if(dsmethod == "MWReg") {
            srchstr = paste("MWReg-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          } else if(dsmethod == "CIReg") {
            srchstr = paste("CIReg-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          } else if(dsmethod == "ITReg") {
            srchstr = paste("ITReg-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          } else if(dsmethod == "MWRObs") {
            srchstr = paste("MWRObs-", varnm, "-", monnm, "-", acunm, ".csv", sep="")
          }
          InDFile = list.files(indir, pattern = glob2rx(srchstr), full.names = T)

          if(length(InDFile) > 0){

            data = read.csv(InDFile, header=T)

            colcnt = ncol(data) - 2

            for(jj in 2:colcnt){
              colnm = names(data)[jj]
              errdata = data[,c(colnm,"OBS")]
              errdata = na.omit(errdata)
              colnames(errdata) = c("SIM", "OBS")

              if(nrow(errdata) == 0){
                cor = NA; nrmse = NA
              } else {
                cor = format(cor(errdata$SIM, errdata$OBS, method="pearson"), digits=2)
                nrmse = format(nrmse(errdata$SIM, errdata$OBS, norm="sd")/100, digits=2)
              } # End if

              if(jj == 2){
                cortbl = c(monnm, cor)
              } else {
                cortbl = c(cortbl, cor)
              }

            } # Column Loop

            cortbl = t(data.frame(cortbl))
            colnames(cortbl) = c("month", names(data)[2:colcnt])

            if(cnt == 1){
              corerr = as.data.frame(cortbl)
            } else {
              cortmp = as.data.frame(cortbl)
              corerr = rbind.fill(corerr, cortmp)
            }
            cnt = cnt + 1

          } else{
            if(cnt == 1){
              month = monnm; MME = NA; corerr = data.frame(month, MME)
            } else {
              month = monnm; MME = NA; cortmp = data.frame(month, MME)
              corerr = rbind.fill(corerr, cortmp)
            }
            cnt = cnt + 1
          } # IF file.exists

        } # Month Loop

        # Change column order
        mmedata = corerr[c("MME")]
        mondata = corerr[c("month")]
        valdata = corerr[-which(names(corerr) %in% c("month", "MME"))]
        if(ncol(valdata) > 1){valdata = valdata[, order(substr(names(valdata), (nchar(names(valdata))-1), nchar(names(valdata))))]}
        corerr = cbind(mondata, valdata, mmedata)

        if(dsmethod == "BCArea"){
          OutDFile = paste(outdir, "/BCArea-", varnm, "-", acunm, "-TCC.csv", sep="")
        } else if(dsmethod == "BCPoint") {
          OutDFile = paste(outdir, "/BCPoint-", varnm, "-", acunm, "-TCC.csv", sep="")
        } else if(dsmethod == "MWReg") {
          OutDFile = paste(outdir, "/MWReg-", varnm, "-", acunm, "-TCC.csv", sep="")
        } else if(dsmethod == "CIReg") {
          OutDFile = paste(outdir, "/CIReg-", varnm, "-", acunm, "-TCC.csv", sep="")
        } else if(dsmethod == "ITReg") {
          OutDFile = paste(outdir, "/ITReg-", varnm, "-", acunm, "-TCC.csv", sep="")
        } else if(dsmethod == "MWRObs") {
          OutDFile = paste(outdir, "/MWRObs-", varnm, "-", acunm, "-TCC.csv", sep="")
        }
        write.csv(corerr, OutDFile, row.names=F)

        #} # Lead Time
      } # Accumulation Months
    } # Variable


  } # If selected model exists


}

#' @export
Calculate.Covariance.Matrix <- function(obsdir, outdir, stnnms, syear_obs, eyear_obs) {

  varnms = c("prec", "t2m")

  SetWorkingDir(outdir)

  covfile = paste(outdir, "/covariance_matrix_monthly.csv", sep="")

  if(file.exists(covfile)) {

    covdata = read.csv(covfile, header=T)

  } else {

    varcnt = length(varnms)
    stncnt = length(stnnms)

    #### define empty array for covariance ************************
    covdata = array(NA, dim=c(12, varcnt*varcnt))  #Empty array

    for(i in 1:12){

      #### define empty array
      nmonths = (eyear_obs - syear_obs + 1) #Number of rows of array
      #Empty array ************************************************
      dat = array(NA, dim=c(nmonths, stncnt, varcnt))

      for(j in 1:stncnt){
        stnid = stnnms[j]

        #obsdaydata = Get.Station.Data(stnid, obsdir)
        obsdaydata = Get.Station.Data.9var(stnid, obsdir)

        colnames(obsdaydata) = c("year", "month", "yearmon", "date", "prec", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")
        obsdaydata$t2m = (as.numeric(obsdaydata$tmax) + as.numeric(obsdaydata$tmin))/2
        obsdaydata$prec = as.numeric(obsdaydata$prec)

        # select common period
        obsdaydata = obsdaydata[which(obsdaydata$year >= syear_obs & obsdaydata$year <= eyear_obs & obsdaydata$month == i),]

        # calculate monthly mean
        obsmondata = aggregate(cbind(prec, t2m) ~ yearmon, data = obsdaydata, FUN = mean)

        yearmon = as.matrix(obsmondata[,1])

        for(k in 1:varcnt){
          varnm = varnms[k]
          if(dim(obsmondata)[1] == dim(dat)[1]){
            dat[,j,k] = obsmondata[,k+1]
          } else {
            dat[,j,k] = NA
            cat(sprintf("     Util: Month=%s Stn=%s  Var=%s has been removed from the list due to missing data\n", i, stnid, varnm))
          }
        }
      }

      # calculate station average
      prec = rowMeans(dat[,,1], na.rm=T)
      temp = rowMeans(dat[,,2], na.rm=T)
      mondata = data.frame(cbind(prec, temp))

      obsmean = cbind(yearmon, prec, temp)
      colnames(obsmean) = c("yearmon", "prec", "temp")
      outname = paste(outdir, "/ObsMean-", sprintf("M%02d", i), ".csv", sep="")
      write.csv(obsmean, outname, row.names=F)

      covar = as.vector(cov(mondata))

      covdata[i,] = covar
    }


    write.csv(covdata, covfile, row.names=F)

  }

  return(covdata)

}

#' @export
Solar2Sunshine <- function(Rs_MJ, julday, lat, elev) {

  # It is coded based on http://biomet.ucdavis.edu/biomet/Radiation/Wton.htm
  # n : hours of bright sunshine
  # N : potential bright sunshine
  # Rs_mm : surface solar radiation (mm/d)
  # Ra_mm : extraterrestrial solar radiation (mm/d)
  # n/N = 2*(Rs_mm/Ra_mm - 0.25)

  Gsc = 0.082
  PI = lat * pi / 180
  Rs_cal = Rs_MJ * 23.88
  Rs_mm = Rs_MJ/2.45
  df = 1+ 0.033 * cos(2*pi/365*julday)
  del = 0.409 * sin(2*pi/365*julday - 1.39)
  ws = acos(-tan(PI)*tan(del))
  N = ws*24/pi
  Ra_mm = ((24*60/pi)*Gsc*df*(ws*sin(PI)*sin(del)+cos(PI)*cos(del)*sin(ws)))/2.45
  R_ratio = Rs_mm/Ra_mm
  n_ratio = 2*(R_ratio - 0.25)
  n = n_ratio * N

  n[n <= 0] = 0

  n_N = n/N

  out = cbind(n, n_N)
  colnames(out) = c("n", "n_N")

  return(out)

}

#' @export
world.map<- function(database,...){

  center = 200

  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])

  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})

  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1] + 160
  Obj[[2]] <- polygons[,2]

  map(Obj,...)
}




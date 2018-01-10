#' @export
Update.APCC.forecast.data <- function(EnvList, updatemode){

  #library(RCurl)
  #library(XML)
  #library(stringr)

  mdlnms_3mon = EnvList$mdlnms_3mon
  mdlnms_6mon = EnvList$mdlnms_6mon
  mmedir = EnvList$mmedir

  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  odapurl = "http://jpcho96:Lauren02@cis.apcc21.org/opendap"
  downurl = "ftp://Data_RESDAT:12345@adss.apcc21.org"
  maxyear = as.numeric(substr(Sys.Date(),1,4))+1

  if(updatemode == F){
    syear = 1983
  } else {
    syear = as.numeric(substr(Sys.Date(),1,4))
  }
  eyear = as.numeric(substr(Sys.Date(),1,4))+1

  for(ltnm in c("3-MON", "6-MON")){

    if(ltnm == "3-MON") {
      mdlcnt = length(mdlnms_3mon)
      mdlnms = mdlnms_3mon
      ltdir = "3MON"
    } else  {
      mdlcnt = length(mdlnms_6mon)
      mdlnms = mdlnms_6mon
      ltdir = "6MON"
    }

    for(j in 1:mdlcnt){
      mdlnm = mdlnms[j]

      moncnt = length(monnms)
      cat(paste("process now running...",mdlnm," ",j,"/",mdlcnt,sep=""));cat("\n")
      pb <- txtProgressBar(min=1,max=moncnt,style=3)
      for(k in 1:moncnt){
        monnm = monnms[k]
        #pb_yr <- winProgressBar(title=paste("Process now running.. Year period:",syear,"-",eyear,sep=""),
        #                        label="0% done", min=0, max=100, initial=0)
        #mm <- 1
        for(m in syear:eyear){
          #info <- sprintf("%d %d%% done", m,round(mm/(eyear-syear+1)*100),sep="")
          #setWinProgressBar(pb_yr, mm/(eyear-syear+1)*100, label=info)
          #mm <- mm+1
          for(fcstmode in c("HINDCAST", "FORECAST")){

            url = paste(odapurl, "/INDV_MODEL/", ltnm, "/", fcstmode, "/", mdlnm, "/", monnm, "/", m, "/contents.html", sep="")

            if(url.exists(url)) {
              page_parse <- htmlParse(url, encoding = "utf-8")
              links <- getHTMLLinks(url)
              fnames <- links[str_detect(links, ".nc$")]
              tmp  <- matrix(unlist(strsplit(fnames, "/")), nrow = length(fnames), byrow=T)
              fnames <- tmp[, ncol(tmp)]
              #fnames <- fnames[nchar(fnames) < 10]

              fcnt = length(fnames)
              for(n in 1:fcnt){
                fname = fnames[n]
                srcnm = paste(downurl, "/INDV_MODEL/", ltnm, "/", fcstmode, "/", mdlnm, "/", monnm, "/", m, "/", fname, sep ="")
                trgdir = paste(mmedir, "/", ltdir, "/", mdlnm, "/", monnm, "/", m, sep="")
                if(!dir.exists(trgdir)){dir.create(trgdir, showWarnings=F,recursive=T)}
                #SetWorkingDir(trgdir)
                trgnm = paste(trgdir, "/", fname, sep="")
                if(!file.exists(trgnm)){
                  download.file(srcnm, trgnm, method = "auto", quiet = TRUE, mode="wb", cacheOK = TRUE)
                }
              }
            } # EndIF
          } # Forecast mode

        } # Year
        #close(pb_yr)
        setTxtProgressBar(pb,k)
      } # Month
    close(pb)
    } # Model

    cat(sprintf("     DBUpdate: Individual model data has been downloaded: LTime=%s\n", ltdir))

  } # Lead-Time

}

#' @export
MWRObs.Download.Reanalysis1 <- function(EnvList) {

  rnldir = EnvList$rnldir

  setwd(rnldir)

  # Monthly mean air temperature
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/air.mon.mean.nc"
  download.file(url, "air.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean geopotential height
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc"
  download.file(url, "hgt.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean specific humidity
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/shum.mon.mean.nc"
  download.file(url, "shum.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean vertical velocity (Omega)
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/omega.mon.mean.nc"
  download.file(url, "omega.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean U-Wind
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/uwnd.mon.mean.nc"
  download.file(url, "uwnd.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean V-Wind
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc"
  download.file(url, "vwnd.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean sea level pressure
  url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/slp.mon.mean.nc"
  download.file(url, "slp.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean precipitable water
  #url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/pr_wtr.mon.mean.nc"
  #download.file(url, "pr_wtr.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly mean relative humidity
  #url = "ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/rhum.mon.mean.nc"
  #download.file(url, "rhum.mon.mean.nc", quiet = TRUE, mode = "wb")

  # Monthly Sea Surface Temperature
  #url = "ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc"
  #download.file(url, "sst.mnmean.nc", quiet = TRUE, mode = "wb")

  cat(sprintf("     DBUpdate: Reanalysis1 data was successfully updated!\n"))

}

#' @export
CIReg.Get.APCC.Climate.Index <- function(EnvList) {

  idxdir = EnvList$idxdir
  apccidxs = EnvList$apccidxs

  destdir = paste(idxdir, "/APCC", sep="")
  SetWorkingDir(destdir)

  fnames = c("AAO_2D.txt", "AO_2D.txt", "ATLTRI_2D.txt", "EOFPAC_2D.txt", "NAO_2D.txt", "NINO3_2D.txt", "NINO4_2D.txt", "NINO12_2D.txt", "NINO34_2D.txt", "NOI_2D.txt", "NP_2D.txt", "ONI_2D.txt", "PACWARM_2D.txt", "PNA_2D.txt", "QBO_2D.txt", "SOI_2D.txt", "TNA_2D.txt", "TNI_2D.txt", "TSA_2D.txt", "WP_2D.txt")
  fcnt = length(fnames)

  for(i in 1:fcnt){
    fname = fnames[i]
    url = paste("http://www.apcc21.org/cmm/fms/FileDown2.do?atchFileId=", fname, "&fileSn=1", sep="")
    download.file(url, fname, quiet = TRUE, mode = "wb")

    idxnm = substr(fname, 1, nchar(fname)-7)
    #data = read.table(fname, header = T, skip = 1)
    data = read.table(fname, header = T)
    syear = data[1, 1]; eyear = data[nrow(data), 1]
    val = as.vector(t(data[,2:13]))
    yearmon = substr(seq(as.Date(paste(syear,"-01-01", sep="")), as.Date(paste(eyear,"-12-01", sep="")), by="month"), 1, 7)

    df = cbind(yearmon, val)
    colnames(df) = c("yearmon", idxnm)

    if(i == 1){
      outdata = df
    } else {
      imsi = df
      outdata = merge(outdata, imsi, by="yearmon", all=T)
    }
  }

  outdata[outdata == "-999.999"] = NA
  # Select only data included in the index list
  outdata = outdata[c("yearmon", apccidxs)]

  SetWorkingDir(idxdir)
  write.csv(outdata, "CIndex-APCCC.csv", row.names=F)

  return(outdata)

}

#' @export
CIReg.Get.ESRL.Climate.Index <- function(EnvList) {

  cpcidxs = EnvList$cpcidxs
  idxdir = EnvList$idxdir

  cat("process now running...");cat("\n")
  pb <- txtProgressBar(min=1,max=40,style=3)
  # Loop for 40 Climate Index
  for(i in 1:40){

    ######### Define URL for each index
    if(i == 1){ # PNA: Pacific North American Index ############################
      IName = "PNA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/pna.data"
    }
    if(i == 2){ # EP: East Pacific/North Pacific Oscillation
      IName = "EP"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/epo.data"
    }
    if(i == 3){ # WP: West Pacific Oscillation ################################
      # WP: daily Western Pacific Index(WP)
      IName = "WP"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/wp.data"
    }
    if(i == 4){ # EA: Eastern Asia/Western Russia
      IName = "EA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/ea.data"
    }
    if(i == 5){ # NAO: North Atlantic Oscillation ##############################
      IName = "NAO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/nao.data"
    }
    if(i == 6){ # SOI: Southern Oscillation Index ##############################
      IName = "SOI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/soi.data"
    }
    if(i == 7){ # NINO3: Eastern Tropical Pacific SST
      IName = "NINO3"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/nina3.data"
    }
    if(i == 8){ # BEST: Bivariate ENSO Timeseries
      IName = "BEST"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/censo.data"
    }
    if(i == 9){ # TNA: Tropical Northern Atlantic Index ########################
      IName = "TNA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/tna.data"
    }
    if(i == 10){ # TSA: Tropical Southern Atlantic Index #######################
      IName = "TSA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/tsa.data"
    }
    if(i == 11){ # WHWP: Western Hemisphere warm pool ###########################
      IName = "WHWP"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/whwp.data"
    }
    if(i == 12){ # ONI: Oceanic Nino Index #####################################
      IName = "ONI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/oni.data"
    }
    if(i == 13){ # MEI: Multivariate ENSO Index ################################
      IName = "MEI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/mei.data"
    }
    if(i == 14){ # NINO12: Extreme Eastern Tropical Pacific SST
      IName = "NINO12"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/nina1.data"
    }
    if(i == 15){ # NINO4: Central Tropical Pacific SST
      IName = "NINO4"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/nina4.data"
    }
    if(i == 16){ # NINO34: East Central Tropical Pacific SST
      IName = "NINO34"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/nina34.data"
    }
    if(i == 17){ # PDO: Pacific Decadal Oscillation ############################
      IName = "PDO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/pdo.data"
    }
    if(i == 18){ # NOI: Northern Oscillation Index #############################
      IName = "NOI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/noi.data"
    }
    if(i == 19){ # NP: North Pacific pattern
      IName = "NP"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/np.data"
    }
    if(i == 20){ # TNI: Indices of El Nino evolution
      IName = "TNI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/tni.data"
    }
    if(i == 21){ # AO: Antarctic Oscillation ###################################
      IName = "AO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/ao.data"
    }
    if(i == 22){ # AAO: Antarctic Oscillation
      IName = "AAO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/aao.data"
    }
    if(i == 23){ # PACWARM: Pacific Warmpool (1st EOF of SST (60e-170E, 15S-15N) SST EOF)
      IName = "PACWARM"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/pacwarm.data"
    }
    if(i == 24){ # EOFPAC: Tropical Pacific SST EOF (1st EOF of SST 20N-20S, 120E-60W)
      IName = "EOFPAC"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/eofpac.data"
    }
    if(i == 25){ # ATLTRI: Atlantic Tripole SST EOF (1st EOF of SST 10N-70N, 0-80W)
      IName = "ATLTRI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/atltri.data"
    }
    if(i == 26){ # AMO: Atlantic  Multidecadal Oscillation #####################
      IName = "AMO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/amon.us.data"
    }
    if(i == 27){ # AMM: Atlantic Meridional Mode
      IName = "AMM"
      URL = "http://www.esrl.noaa.gov/psd/data/timeseries/monthly/AMM/ammsst.data"
    }
    if(i == 28){ # NTA: North Tropical Atlantic SST Index
      IName = "NTA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/NTA.data"
    }
    if(i == 29){ # CAR: Caribbean SST Index
      IName = "CAR"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/CAR.data"
    }
    if(i == 30){ # AMOSM, smoothed: Atlantic Multidecadal Osillation
      IName = "AMOSM"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/amon.sm.data"
    }
    if(i == 31){ # QBO: Quasi-Biennial Oscillation
      IName = "QBO"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/qbo.data"
    }
    if(i == 32){ # GIAM: Globally Integrated Angular Momentum
      IName = "GIAM"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/glaam.data.scaled"
    }
    if(i == 33){ # ESPI: ENSO precipitation index
      IName = "ESPI"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/espi.data"
    }
    if(i == 34){ # CIP: Central Indian Precipitation
      IName = "CIP"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/indiamon.data"
    }
    if(i == 35){ # SahelRain: Sahel Standardized Rainfall (20-8N, 20W-10E)
      IName = "SahelRain"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/sahelrain.data"
    }
    if(i == 36){ # SWM: SahelArea averaged precipitation for Arizona and New Mexico
      # SWM: SW Monsoon Region Rainfall #############################
      IName = "SWM"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/swmonsoon.data"
    }
    if(i == 37){ # NBRA: Northeast Brazil Rainfall Anomaly
      IName = "NBRA"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/brazilrain.data"
    }
    if(i == 38){ # GML: Global Mean Lan/Ocean Temperature ######################
      IName = "GML"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/gmsst.data"
    }
    if(i == 39){ # Solar: Solar Flux (10.7cm)
      IName = "Solar"
      URL = "http://www.esrl.noaa.gov/psd/data/correlation/solar.data"
    }
    if(i == 40){ # ESL: Equatorial Eastern Pacific SLP #########################
      IName = "ESL"
      URL = "http://www.cpc.ncep.noaa.gov/data/indices/repac_slpa.for"
    }


    ######## Extract table from the URL
    tbl = readLines(URL)
    tbl = strsplit(tbl, " +")


    cnt = 1
    syear = 10000
    eyear = 100
    nline = length(tbl)

    for(j in 1:nline){

      colnum = length(unlist(tbl[j]))
      if(colnum >= 13 & colnum <= 15){

        scolnum = colnum - 11
        yearcol = colnum - 12
        if(cnt == 1){
          curyear = as.numeric(unlist(tbl[j])[yearcol])
          var = data.frame(unlist(tbl[j])[scolnum:colnum])
          if(curyear <= syear) {syear = curyear}
          if(curyear >= eyear) {eyear = curyear}
          cnt = cnt + 1
        } else {
          curyear = as.numeric(unlist(tbl[j])[yearcol])
          temp = data.frame(unlist(tbl[j])[scolnum:colnum])
          var = rbind(var, temp)
          if(curyear <= syear) {syear = curyear}
          if(curyear >= eyear) {eyear = curyear}
          cnt = cnt + 1
        }
      }
    }

    nyear = eyear - syear + 1
    year = rep(seq(syear, eyear), each=12)
    mon = rep(seq(1,12), nyear)
    yearmon = paste(year, "-", formatC(mon, width=2, flag="0"), sep="")

    outtbl = cbind(yearmon, var)
    colnames(outtbl) = c("yearmon",  IName)

    # Combine(merge) each indice into one table, outdata
    if(i == 1){
      outdata = outtbl
    } else {
      imsi = outtbl
      outdata = merge(outdata, imsi, by="yearmon", all=T)
    }
  setTxtProgressBar(pb,i)
  }
  close(pb)
  # Extract indices which are real-time updated
  outdata = outdata[,c("yearmon", cpcidxs)]


  # each index has different NA values
  outdata[outdata == "-99"] = NA
  outdata[outdata == "-99.90"] = NA
  outdata[outdata == "-99.900"] = NA
  outdata[outdata == "-99.99"] = NA
  outdata[outdata == "-99.990"] = NA
  outdata[outdata == "-99.000"] = NA
  outdata[outdata == "-9.99"] = NA
  outdata[outdata == "-9.90"] = NA
  outdata[outdata == "-9.000"] = NA
  outdata[outdata == "-999"] = NA
  outdata[outdata == "-999.0"] = NA
  outdata[outdata == "-999.00"] = NA
  outdata[outdata == "-999.000"] = NA
  outdata[outdata == "9999"] = NA
  outdata[outdata == "999.9"] = NA

  #idxdir = SetWorkingDir(idxdir)
  #SetWorkingDir(idxdir)
  write.csv(outdata, paste(idxdir,"/CIndex-CPC.csv",sep=""), row.names=F)

  return(outdata)

}

#' @export
CIReg.Update.Climate.Index <- function(EnvList){

  idxdir = EnvList$idxdir
  idxfile = EnvList$idxfile

  options(warn=-1)

  # Add one more years to consider lagtime
  sdate = as.Date("1948-01-01")
  eyear = as.numeric(format(Sys.Date(), format="%Y")) + 1
  edate = as.Date(paste(eyear, "-12-31", sep=""))
  yearmon = as.data.frame(substr(seq(sdate, edate, by="mon"),1,7))
  colnames(yearmon) = c("yearmon")

  cpc = CIReg.Get.ESRL.Climate.Index(EnvList)
  apcc = CIReg.Get.APCC.Climate.Index (EnvList)

  #idx = merge(yearmon, cpc, by="yearmon", all=T)
  idx = merge(cpc, apcc, by="yearmon", all=T)

  IdxDFile = paste(idxdir, "/", idxfile, sep="")
  write.csv(idx, IdxDFile, row.names=F)

  cat(sprintf("     DBUpdate: Climate Indices provided by CPC and APCC are successfully combined!\n"))
  #return(idx)

}

#' @export
kma.asos.daily.update <- function(EnvList) {

  # 1:지점, 2:일시,	3:평균기온(°C),	4:최저기온(°C),	5:최저기온시각(hhmi),	6:최고기온(°C), 7:최고기온시각(hhmi), 8:강수계속시간(hr),	9:10분최다강수량(mm),	10:10분최다강수량시각(hhmi)
  # 11:1시간최다강수량(mm), 12:1시간최다강수량시각(hhmi), 13:일강수량(mm), 14:최대순간풍속(m/s), 15:최대순간풍속풍향(16방위), 16:최대순간풍속시각(hhmi)	17:최대풍속(m/s), 18:최대풍속풍향(16방위), 19:최대풍속시각(hhmi), 20:평균풍속(m/s)
  #	21:풍정합(100m), 22:평균이슬점온도(°C), 23:최소상대습도(%), 24:최소상대습도시각(hhmi), 25:평균상대습도(%), 26:평균증기압(hPa), 27:평균현지기압(hPa), 28:최고해면기압(hPa), 29:최고해면기압시각(hhmi), 30:최저 해면기압(hPa),
  # 31:최저해면기압시각(hhmi), 32:평균해면기압(hPa), 33:가조시간(hr), 34:합계일조시간(hr), 35:1시간최다일사시각(hhmi), 36:1시간최다일사량(MJ/m2), 37:합계일사(MJ/m2), 38:일최심신적설(cm), 39:일최심신적설시각(hhmi), 40:일최심적설(cm)
  # 41:일최심적설시각(hhmi)	42:합계3시간신적설(cm)	43:평균전운량(1/10)	44:평균중하층운량(1/10)	45:평균지면온도(°C)	46:최저초상온도(°C)	47:평균5cm지중온도(°C)	48:평균10cm지중온도(°C)	49:평균20cm지중온도(°C)	50:평균30cm지중온도(°C)
  # 51:0.5m지중온도(°C)	52:1.0m지중온도(°C)	53:1.5m지중온도(°C)	54:3.0m지중온도(°C)	55:5.0m지중온도(°C)	56:합계대형증발량(mm)	57:합계소형증발량(mm)	58:9-9강수(mm)	59:안개계속시간(hr)
  #==> Selected:

  obsdir = EnvList$obsdir
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  downdir = paste(obsdir, "/asos_daily_download", sep="")

  flist = list.files(downdir, pattern = "*.csv", full.names = T)
  fcnt = length(flist)

  for(i in 1:fcnt){

    fname = flist[i]
    if(i == 1){
      data = read.csv(fname, header = T)
    } else {
      imsi = read.csv(fname, header = T)
      data = rbind(data, imsi)
    }
  }

  data = data[ , c(1, 2, 13, 6, 4, 20, 25, 37, 34, 43, 3)]
  data[is.na(data[,3]), 3] = 0
  data[, 7] = data[, 7] / 100
  colnames(data) = c("Stn", "Date", "Pcp(mm)", "Tmax(c)", "Tmin(c)", "WSpeed(m/s)", "RHumidity(fr)", "SRad(MJ/m2)", "SShine(hr)", "Cloud(1/10)", "Tavg(c)")

  #stnids = unique(data$Stn)
  StnDFile = paste(stndir, "/", stnfile, sep="")
  stninfo = read.csv(StnDFile, header=T)
  stnids = matrix(stninfo$ID)
  stncnt = length(stnids)

  # Get the latest date from all stations
  edate = max(as.Date(data$Date))

  cat("process now running... Station data update");cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  for(j in 1:stncnt){
    setTxtProgressBar(pb,j)
    stnid = as.numeric(substr(stnids[j], 3, 5))
    stndata = data[which(data$Stn == stnid),]
    stndata$Date = as.Date(stndata$Date)

    stndata = stndata[order(stndata$Date),]
    sdate = as.Date(stndata[1, "Date"])
    #edate = as.Date(stndata[nrow(stndata), "Date"])
    daytime = as.data.frame(seq(sdate, edate, by="day")); colnames(daytime) = c("Date")


    stndata = merge(daytime, stndata, by="Date", all=T)

    stndata$Year = substr(stndata$Date, 1, 4)
    stndata$Mon = substr(stndata$Date, 6, 7)
    stndata$Day = substr(stndata$Date, 9, 10)

    stndata[is.na(stndata)] = "-99.00"

    stndata = stndata[, c("Year", "Mon", "Day", "Pcp(mm)", "Tmax(c)", "Tmin(c)", "WSpeed(m/s)", "RHumidity(fr)", "SRad(MJ/m2)", "SShine(hr)", "Cloud(1/10)", "Tavg(c)")]
    ofname = sprintf("ID%03d.csv", stnid)
    setwd(obsdir)
    write.csv(stndata, ofname, row.names = F)
  }
  close(pb)

}

#' @export
ghcn.daily.update <- function(EnvList, cntry) {

  ghcndir = EnvList$ghcndir
  obsdir = EnvList$obsdir
  stndir = EnvList$stndir
  syear_obs = as.numeric(EnvList$syear_obs)
  #syear_obs = 1984
  eyear_obs = as.numeric(EnvList$eyear_obs)

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  setwd(ghcndir)

  stnfile = paste(ghcndir, "/ghcnd-stations.txt", sep="")
  ivntfile = paste(ghcndir, "/ghcnd-inventory.txt", sep="")

  if(!file.exists(stnfile)){
    url = "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
    download.file(url, "ghcnd-stations.txt", quiet = TRUE, mode = "wb")
  }

  url = "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
  download.file(url, "ghcnd-inventory.txt", quiet = TRUE, mode = "wb")

  stndata = read.fwf(stnfile, header = F, widths = c(11, 10, 10, 7, 80))
  colnames(stndata) = c("ID", "Lon", "Lat", "Elev", "Ename")
  srchstr = paste("^", cntry, sep ="")
  stndata = stndata %>% filter(str_detect(ID, srchstr))
  stnnms = stndata$ID
  stncnt = length(stnnms)

  ivntdata = read.table(ivntfile, header = F, sep = "")
  colnames(ivntdata) = c("ID", "Lon", "Lat", "var", "syear", "eyear")
  cnt = 1
  for(i in 1:stncnt){
    stnnm = stnnms[i]
    ivnt = ivntdata[which(ivntdata$ID == stnnm & (ivntdata$var == "PRCP" | ivntdata$var == "TMAX" | ivntdata$var == "TMIN")), ]

    if(nrow(ivnt) > 0){

      syear = max(ivnt$syear)
      eyear = min(ivnt$eyear)

      if(syear <= syear_obs & eyear >= eyear_obs){
        stninfo = stndata[which(stndata$ID == stnnm), c("Lon", "Lat", "Elev", "ID", "Ename")]
        stninfo$SYear = syear

        url <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/", stnnm, ".dly", sep = "")
        download.file(url, paste(stnnm, ".dly", sep = ""), quiet = TRUE, mode = "wb")

        dlyfile = paste(stnnm, ".dly", sep="")
        csvfile = paste(stnnm, ".csv", sep="")
        convert.ghcn.daily (ghcndir, obsdir, dlyfile, csvfile)

        if(cnt == 1) {
          stnout = stninfo
        } else {
          imsi = stninfo
          stnout = rbind(stnout, imsi)
        }
        cnt = cnt + 1
      }

    }
  }

  stnoutfile = paste(stndir, "/", cntry , "-GHCN_Stations.csv", sep ="")
  write.csv(stnout, stnoutfile, row.names=F)

}

#' @export
convert.ghcn.daily <- function(dlydir, csvdir, dlyfile, csvfile) {

  #dlydir = "E:/SForecast-Test/Database/observed/ghcn-daily_download"
  #csvdir = "E:/SForecast-Test/Database/observed"
  #dlyfile = "KS000047112.dly"
  #csvfile = "KS000047112.csv"

  DlyDFile = paste(dlydir, "/", dlyfile, sep ="")
  CsvDFile = paste(csvdir, "/", csvfile, sep ="")

  data = read.fwf(DlyDFile, header = F, widths = c(11, 6, 4, rep(c(5, 3), 31)))
  data = data[, c(1:3, seq(4, 64, by=2))]
  colnames(data) = c("ID", "yearmon", "varnm", seq(1, 31))
  data[data == "-9999"] = NA

  varnms = c("PRCP", "TMAX", "TMIN", "AWND", "PSUN", "TAVG")
  varcnt = length(varnms)

  sdate = as.Date(sprintf("%s-%s-01", substr(data[1, "yearmon"], 1, 4), substr(data[1, "yearmon"], 5, 6)))
  daycnt  = numberOfDays(as.Date(paste(substr(data[nrow(data), "yearmon"], 1, 4), "-", substr(data[nrow(data), "yearmon"], 5, 6), "-01", sep="")))
  edate = as.Date(sprintf("%s-%s-%02d", substr(data[nrow(data), "yearmon"], 1, 4), substr(data[nrow(data), "yearmon"], 5, 6), daycnt))
  dates = as.data.frame(seq(sdate, edate, by = "day"))
  colnames(dates) = "date"
  dates$date = as.character(dates$date)

  cnt = 1
  for(i in 1:varcnt){

    varnm = varnms[i]
    val = as.vector(t(data[which(data$varnm == varnm), 4:34]))

    if(length(val) > 0){

      if(varnm == "PSUN") {
        val = as.numeric(val) / 100.0
      } else {
        val = as.numeric(val) / 10.0
      }

      yearmon = data[which(data$varnm == varnm), "yearmon"]
      yearmonstr = rep(yearmon, each=31)
      daystr = sprintf("%02d",rep(seq(1:31), length(yearmon)))
      datestr = paste(substr(yearmonstr, 1, 4), "-", substr(yearmonstr, 5, 6), "-", daystr, sep="")

      if(cnt == 1){
        varout = cbind(datestr, val)
        colnames(varout) = c("date", varnm)

      } else {
        imsi = cbind(datestr, val)
        colnames(imsi) = c("date", varnm)
        varout = merge(varout, imsi, all=T)
      }
      cnt = cnt +1
    }

  }

  out = merge(dates, varout)

  out[is.na(out)] = "-99.00"
  out$Year = substr(out$date, 1, 4)
  out$Mon = substr(out$date, 6, 7)
  out$Day = substr(out$date, 9, 10)

  if(!("PRCP" %in% colnames(out))) out$PRCP = "-99.00"
  if(!("TMAX" %in% colnames(out))) out$TMAX = "-99.00"
  if(!("TMIN" %in% colnames(out))) out$TMIN = "-99.00"
  if(!("TAVG" %in% colnames(out))) out$TAVG = "-99.00"
  if(!("AWND" %in% colnames(out))) out$AWND = "-99.00"
  if(!("PSUN" %in% colnames(out))) out$PSUN = "-99.00"
  if(!("SRAD" %in% colnames(out))) out$SRAD = "-99.00"
  if(!("RHUM" %in% colnames(out))) out$RHUM = "-99.00"
  if(!("CLOD" %in% colnames(out))) out$CLOD = "-99.00"

  out = out[c("Year", "Mon", "Day", "PRCP", "TMAX", "TMIN", "AWND", "RHUM", "SRAD", "PSUN", "CLOD", "TAVG")]
  colnames(out) = c("Year", "Mon", "Day", "Pcp(mm)", "Tmax(c)", "Tmin(c)", "WSpeed(m/s)", "RHumidity(fr)", "SRad(MJ/m2)", "SShine(hr)", "Cloud(1/10)", "Tavg(c)")

  write.csv(out, CsvDFile, row.names = F)

}

#' @export
KMA.ASOS.Observation.Update <- function(EnvList, fillingmode) {

  tscale = EnvList$tscale

  kma.asos.daily.update (EnvList)
  Cal.Monthly.Station.Mean (EnvList, varnms = c("prec", "t2m"))

  if(tscale=="daily" & fillingmode == T){
    Fill.Missing.KMA.Daily.9var (EnvList)
  }

  if(tscale=="hourly" & fillingmode == T){
    kma.asos.hourly.station.update (EnvList)
    Fill.Missing.KMA.Hourly.9var (EnvList)
  }

}


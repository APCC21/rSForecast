#' @export
BCPoint.Get.Data.Periods <- function(mdldir) {

  #mdldir = "F:/SForecast-RT/Database/apcc-mme/3MON/NASA"

  syear = 20000
  eyear = 0

  # Get list of all sub folders
  dirnms = list.dirs(mdldir, recursive = F)
  strs = matrix(unlist(strsplit(dirnms, "/")), nrow=length(dirnms), byrow=T)
  dirnms = dirnms[which(nchar(strs[,ncol(strs)]) == 3)] #3MON-NASA contains 2009 folder --> Exclude nchar == 4
  monnms = matrix(unlist(strsplit(dirnms, "/")), nrow=length(dirnms), byrow=T)
  monnms = monnms[,ncol(monnms)]

  moncnt = length(monnms)
  for(i in 1:moncnt){
    yeardirnms = list.dirs(dirnms[i], recursive = F)
    yearnms = as.numeric(substr(yeardirnms, nchar(yeardirnms)-3, nchar(yeardirnms)))
    temyear = min(yearnms, na.rm=T)
    if(i==1 & temyear < syear) syear = temyear
    if(i>1 & temyear > syear) syear = temyear
    temyear = max(yearnms, na.rm=T)
    if(temyear > eyear) eyear = temyear

  }

  outList = list("syear"=syear, "eyear"=eyear)

  return(outList)

}

#' @export
BCPoint.Extract.MME.Stations <- function(EnvList, mmetype, varnms) {

  prjdir = EnvList$prjdir
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  mmedir = EnvList$mmedir
  syear = as.numeric(EnvList$syear_mme)
  eyear = as.numeric(substr(Sys.Date(), 1, 4))

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon


  ##### Define folder names
  cmmedir = paste(mmedir, "/", mmetype, sep="")
  outdir = paste(prjdir, "/BCPoint/", mmetype, "/01_extracted", sep="")
  comdir = paste(outdir, "/common", sep="")
  SetWorkingDir(comdir)

  # remove all files from the common dir
  flist = list.files(comdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  ### Get Station ID, lat, and Lon information
  setwd(stndir)
  stns = read.csv(stnfile, header=T)
  stns = stns[,c("ID", "Lon", "Lat")]

  ltcnt = as.numeric(substr(mmetype,1,1))

  # Count for Loops
  stncnt = nrow(stns); varcnt = length(varnms); mdlcnt = length(mdlnms)
  cat(paste("process now running ",mmetype,"... extract data to common directory",sep=""));cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  #### Station Loop
  for(i in 1:stncnt){
    setTxtProgressBar(pb,i)
    lon = stns[i,2]; lat = stns[i,3]; StnID = stns[i,1]
    ### Variable Loop
    for(j in 1:varcnt){

      varnm = varnms[j]
      ## Model Loop
      for(k in 1:mdlcnt){
        mdlnm = mdlnms[k]
        wdir = paste(cmmedir, "/", mdlnm, sep="")
        FName = list.files(wdir, pattern = glob2rx("*.nc"), recursive = T)[1]
        dir1 = strsplit(FName, "/")[[1]][1]
        dir2 = strsplit(FName, "/")[[1]][2]
        wdir = file.path(cmmedir, mdlnm, dir1, dir2)
        ncfile = strsplit(FName, "/")[[1]][3]
        out = Extract.Point.Value(lon, lat, ncfile, wdir)
        rown = out$rown; coln = out$coln; var = out$var   # var(ensembles, lead time)

        comname = paste(comdir, "/", sprintf("R%03d", rown), sprintf("C%03d", coln), "-", varnm, "-", mdlnm, "-01.csv", sep="")

        if(!file.exists(comname)){

          # define empty array
          esmcnt = max(nrow(var), 100)  #Expected Max. Ensambel numbers (JMA-2007-DEC: esms = 51)
          nmonths = (eyear - syear + 1) * 12 + (ltcnt-1)  #Number of rows of array
          dat = array(NA, dim=c(nmonths, esmcnt, ltcnt))  #Empty array


          ######### Save data into array ############
          #Dummy Loop for year folders
          for(ii in syear:eyear){
            curyear = ii

            #Dummy Loop for month folders
            for(jj in 1:12){
              monnm = monnms[jj]

              wdir = paste(cmmedir, "/", mdlnm, "/", monnm, "/", curyear, sep="")
              ncfile = paste(varnm, ".nc", sep="")
              fname = paste(wdir, "/", ncfile, sep="")

              if(file.exists(fname) & file.size(fname) > 0) {
                out = Extract.Point.Value(lon, lat, ncfile, wdir)
                var = out$var
                var[which(var < -1000)] <- NA


                ### Unit conversion
                if(varnm == "t2m" & mdlnm != "SCM") {var = var - 272.15}   # Kelvin to celsius
                if(varnm == "prec" & mdlnm == "BCC"){var[which(var > 30)] <- NA}

                if(mdlnm == "SCM"){
                  tvar = data.frame(var)
                } else {
                  tvar = t(var)
                }

                esmcol = ncol(tvar)

                # Lead-Time Loop
                for(kk in 1:ltcnt){

                  currow = (curyear - syear) * 12 + jj + (kk-1)
                  dat[currow,1:esmcol,kk] = tvar[kk,1:esmcol]

                } # End of Lead-Time Loop

              } else {

                for(kk in 1:ltcnt){

                  currow = (curyear - syear) * 12 + jj + (kk-1)
                  dat[currow, ,kk] = NA

                } # End of Lead-Time Loop

              } # End of IF

            } # Dummy Loop for Month
          } # Dummy Loop for Year

          ########## Create output file ##############
          for(kk in 1:ltcnt){

            val = as.data.frame(dat[,,kk])
            MME = rowMeans(val, na.rm=T)
            sdate = as.Date(paste(syear, "/1/1", sep=""))
            edate = as.Date(paste(eyear+1, "/", ltcnt-1, "/1", sep=""))
            dates = seq(sdate, edate, "month")
            yearmon = strtrim(dates, 7)

            if(mdlnm == "SCM"){
              output = cbind(yearmon, val)
              colnames(output) = c("yearmon", "SCM")
            } else {
              output = cbind(yearmon, val, MME)
            }

            # exclude columns which have all NA values
            output = output[,-which(apply(output,2,function(x)all(is.na(x))))]

            LTime = sprintf("%02d",kk)
            comname = paste(comdir, "/", sprintf("R%03d", rown), sprintf("C%03d", coln), "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
            write.csv(output, comname, row.names=F)
          }

          #cat(sprintf("     BCPoint: Values are extracted: Station=%s  Variable=%s  Model=%s\n", StnID, varnm, mdlnm))

        }

      } ## End of Model Loop

    }  ### End of Variable Loop
    #close(pb_yr)
  }   #### End of Station Loop
  close(pb)
  ################################################################
  ####### Copy files from Common dir to 01_extracted dir
  cat(paste("copy files from common directory to 01_extracted directory now ",mmetype,sep=""));cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  for(i in 1:stncnt){
    setTxtProgressBar(pb,i)
    lon = stns[i,2]; lat = stns[i,3]; StnID = stns[i,1]

    ### Variable Loop
    for(j in 1:varcnt){
      varnm = varnms[j]

      ## Model Loop
      for(k in 1:mdlcnt){

        mdlnm = mdlnms[k]

        # get start and end year for selected model
        wdir = paste(cmmedir, "/", mdlnm, "/JAN", sep="")
        FileName = list.files(wdir, pattern = glob2rx("*.nc"), recursive = T)[1]
        syear = as.numeric(strsplit(FileName, "/")[[1]][1])
        fname = strsplit(FileName, "/")[[1]][2]
        wdir = paste(cmmedir, "/", mdlnm, "/JAN/", syear, sep="")
        ncfile = paste(wdir, "/prec.nc", sep="")
        out = Extract.Point.Value(lon, lat, ncfile, wdir)
        rown = out$rown; coln = out$coln; var = out$var   # var(ensembles, lead time)

        for(kk in 1:ltcnt){

          LTime = sprintf("%02d",kk)
          fromname = paste(comdir, "/", sprintf("R%03d", rown), sprintf("C%03d", coln), "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
          toname = paste(outdir, "/", StnID, "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
          file.copy(fromname, toname, overwrite=T)
        }

      } ## End of Model Loop
    }  ### End of Variable Loop
  }   #### End of Station Loop
  close(pb)

}

#' @export
BCPoint.Bias.Correction <- function(EnvList, mmetype, varnms){

  #mmetype = "3MON"
  #varnms = c("prec", "t2m")


  prjdir = EnvList$prjdir
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  obsdir = EnvList$obsdir
  mmedir = EnvList$mmedir

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon

  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  ##### Define folder names
  cmmedir = paste(mmedir, "/", mmetype, sep="")
  indir = paste(prjdir, "/BCPoint/", mmetype, "/01_extracted", sep="")
  outdir = paste(prjdir, "/BCPoint/", mmetype, "/02_BiasCorrected", sep="")
  SetWorkingDir(outdir)


  ###### Get Station ID, lat, and Lon information
  setwd(stndir)
  stninfo = read.csv(stnfile, header=T)
  stninfo = stninfo[,c("ID", "Lon", "Lat")]
  stnnms = matrix(stninfo$ID)

  ltcnt = as.numeric(substr(mmetype,1,1))

  # Count for Loops
  stncnt = nrow(stnnms); varcnt = length(varnms); mdlcnt = length(mdlnms)
  #BCPoint: Bias Correction
  cat("process now running...BCPoint: Bias Correction");cat("\n")
  pb <- txtProgressBar(min=1,max=stncnt,style=3)
  #### Station Loop
  for(i in 1:stncnt){
    stnid = stnnms[i]

    ### Variable Loop
    for(j in 1:varcnt){
      varnm = varnms[j]

      ## Model Loop
      for(k in 1:mdlcnt){
        mdlnm = mdlnms[k]

        # get start and end years from MME data
        fname = paste(indir, "/", stnid, "-", varnm, "-", mdlnm, "-01.csv", sep="")
        dat = read.csv(fname, header=T)
        years = as.numeric(substr(dat$yearmon, 1, 4))
        syear = min(years)
        eyear = max(years)

        # read obs data
        obsdaydata = Get.Station.Data.9var(stnid, obsdir)
        colnames(obsdaydata) = c("year", "month", "yearmon", "date", "prec", "tmax", "tmin", "wspd", "rhum", "rsds", "sshine", "cloud", "tavg")
        obsdaydata$t2m = (obsdaydata$tmax + obsdaydata$tmin)/2

        # select common period
        obsdaydata = obsdaydata[which(obsdaydata$year >= syear & obsdaydata$year <= eyear),]

        # calculate monthly mean
        obsmondata = aggregate(cbind(prec, t2m) ~ yearmon, data = obsdaydata, FUN = mean)
        obsdata = obsmondata[,c("yearmon", varnm)]
        colnames(obsdata) = c("yearmon", "obs")
        obsdata$mon = as.numeric(substr(obsdata$yearmon, 6, 7))

        # Lead-Time Loop
        for(ii in 1:ltcnt){

          LTime = sprintf("%02d",ii)

          outname = paste(outdir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, "-BC.csv", sep="")

          mmename = paste(indir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
          mmedata = read.csv(mmename, header=T)
          #if(which(mmedata[,2] > 30000)){print(1)}
          #if(which(mmedata[,2] < -30000)){print(1)}
          #Esemble Loop
          esmcnt = ncol(mmedata) - 1
          for(jj in 1:esmcnt){

            mmeesm = mmedata[, c(1, jj+1)]
            colnames(mmeesm) = c("yearmon", "mme")
            mmeesm$mon = as.numeric(substr(mmeesm$yearmon, 6, 7))

            ############ Calculae QMapping fit
            for(kk in 1:12){

              obsmon = obsdata[which(obsdata$mon == kk),c("yearmon", "obs")]
              #obsmon = as.data.frame(obsmon)

              mmemon = mmeesm[which(mmeesm$mon == kk),c("yearmon", "mme")]
              #mmemon = as.data.frame(mmemon)

              if(all(is.na(mmemon$mme))){
                mmetmp = mmemon
                #mmetmp = as.data.frame(mmeesm[which(mmeesm$mon == kk),c("yearmon")])
                #mmetmp$val = NA
                colnames(mmetmp) = c("yearmon", varnm)
              } else {

                comdata = merge(obsmon, mmemon)
                comdata = na.omit(comdata)

                # Use the overall avearage if there is no common period
                if(nrow(comdata) == 0){
                  obsmean = mean(obsmon$obs, na.rm=T)
                  mmemean = mean(mmemon$mme, na.rm=T)
                }else{
                  obsmean = mean(comdata$obs)
                  mmemean = mean(comdata$mme)
                }

                if(varnm == "t2m"){
                  mmetmp = mmemon
                  mmetmp$mme = obsmean + (mmemon$mme - mmemean)
                  colnames(mmetmp) = c("yearmon", varnm)
                }else{
                  mmetmp = mmemon
                  # in case greater than average, add
                  mmetmp$mme[which((mmetmp$mme - mmemean) >= 0)] = obsmean + (mmetmp$mme[which((mmetmp$mme - mmemean) >= 0)] - mmemean)
                  # in case less than avearage, use ratio
                  mmetmp$mme[which((mmemon$mme - mmemean) < 0)] = obsmean / mmemean * mmetmp$mme[which((mmemon$mme - mmemean) < 0)]
                  colnames(mmetmp) = c("yearmon", varnm)

                }

              }

              if(kk == 1) {
                mmeadj = mmetmp
              }else {
                mmeadj = rbind(mmeadj, mmetmp)
              }
            }  ## End of KK Loop

            #colnames(mmeadj) = c("yearmon",varnm)
            mmeadj = mmeadj[order(mmeadj$yearmon),]
            if(jj == esmcnt){
              colnm = "MME"
            } else {
              colnm = sprintf("E%02d",jj)
            }
            colnames(mmeadj) = c("yearmon", colnm)

            if(jj == 1){
              mmeout = mmeadj
            } else {
              mmeout = merge(mmeout, mmeadj, all=T)
            }

          } ## End of jj Loop


          outname = paste(outdir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, "-BC.csv", sep="")
          write.csv(mmeout, outname, row.names=F)

        }
        #cat(sprintf("     BCPoint: Bias Correction has been finished: Station=%s  Variable=%s  Model=%s\n", stnid, varnm, mdlnm))

      } ## End of Model Loop
    }  ### End of Variable Loop
    setTxtProgressBar(pb,i)
  }   #### End of Station Loop
  close(pb)
  #   ################### Do Bias Correction separately for SCM ####################
  #   for(i in 1:stncnt){
  #     stnid = stnnms[i]
  #
  #     for(j in 1:varcnt){
  #       varnm = varnms[j]
  #
  #       for(ii in 1:ltcnt){
  #         LTime = sprintf("%02d",ii)
  #
  #         for(k in 1:mdlcnt){
  #           mdlnm = mdlnms[k]
  #           inname = paste(outdir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, "-BC.csv", sep="")
  #
  #           if(k == 1){
  #             data = read.csv(inname, header=T)
  #             data = data[c("yearmon", "MME")]
  #             colnames(data) = c("yearmon", mdlnm)
  #           } else {
  #             temp = read.csv(inname, header=T)
  #             temp = temp[c("yearmon", "MME")]
  #             colnames(temp) = c("yearmon", mdlnm)
  #
  #             data = merge(data, temp, by="yearmon")
  #           }
  #         }
  #         yearmon = data$yearmon
  #         imsi = data[,2:ncol(data)]
  #         mmedata = rowMeans(imsi)
  #
  #         outdata = cbind(as.character(yearmon), mmedata)
  #         colnames(outdata) = c("yearmon", "MME")
  #
  #         outname = paste(outdir, "/", stnid, "-", varnm, "-SCM-", LTime, "-BC.csv", sep="")
  #         write.csv(outdata, outname, row.names=F)
  #
  #       }
  #     }
  #   }

}

#' @export
BCPoint.Combine.Forecast.Output <- function(EnvList, mmetype, varnm, lagtime){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon

  syear = syear_mme
  eyear = as.numeric(substr(Sys.Date(),1,4))+1
  #eyear = eyear_obs

  # define working directory
  bcadir = paste(prjdir, "/BCPoint/", mmetype, "/02_BiasCorrected", sep="")

  mdlcnt = length(mdlnms)

  # Max LeadTime from mmetype
  ltcnt = as.numeric(substr(mmetype,1,1))
  cnt = 1
  out = NA
  for(i in lagtime:ltcnt){

    for(j in 1:mdlcnt){

      mdlnm = mdlnms[j]

      srchstr = paste("*-", varnm, "-*", mdlnm, sprintf("-%02d",i), "-BC.csv", sep="")
      flist = list.files(bcadir, pattern = glob2rx(srchstr), full.names = T)
      fnames = list.files(bcadir, pattern = glob2rx(srchstr), full.names = F)
      stncnt = length(flist)

      for(k in 1:stncnt){
        InDFile = flist[k]
        fname = fnames[k]

        stnnm = unique(matrix(unlist(strsplit(fname[], "-")), nrow=length(fname), byrow=T)[,1])

        data = read.csv(InDFile, header=T)
        data = data[c("yearmon", "MME")]
        colnames(data) = c("yearmon", stnnm)

        if(k == 1){
          StnData = data
        } else {
          imsi = data
          StnData = merge(StnData, imsi, by="yearmon", all=T)
        }
      }

      # calculate station meam
      rcnt = nrow(StnData)
      StnData$MME = NA
      for(ii in 1:rcnt){
        if(!all(is.na(StnData[ii, 2:(ncol(StnData)-1)]))){
          StnData[ii, c("MME")] = mean(as.numeric(StnData[ii, 2:(ncol(StnData)-1)]), na.rm=T)
        }
      }

      StnMME = StnData[c("yearmon", "MME")]
      colnames(StnMME) = c("yearmon", sprintf("%s%02d",mdlnm,i))

      if(cnt == 1){
        out = StnMME
        cnt = cnt + 1
      } else {
        tmp = StnMME
        out = merge(out, tmp, by="yearmon", all=T)
        cnt = cnt + 1
      }
    } # Model Loop

  } # LeadTime Loop

  if(!is.na(out)){
    # Overall MME by considering different models and lead times
    rcnt = nrow(out)
    out$MME = NA
    for(i in 1:rcnt){
      if(!all(is.na(out[i, 2:(ncol(out)-1)]))){
        out[i, c("MME")] = mean(as.numeric(out[i, 2:(ncol(out)-1)]), na.rm=T)
      }
    }

    ### Get observed data and merge to output df
    VarDFile = paste(vardir, "/", varfile, sep="")
    obs = read.csv(VarDFile, header=T)
    obs = obs[,c("yearmon", varnm)]
    colnames(obs) = c("yearmon", "OBS")

    out = merge(out, obs, by="yearmon", all=T)

    ### Calculate Climatology and merge to output df
    obs$month = as.numeric(substr(obs$yearmon, 6, 7))
    clim = aggregate(OBS ~ month, data = obs, FUN = mean)
    out$CLIM = NA
    for(i in 1:12){
      out[which(as.numeric(substr(out$yearmon,6,7)) == i), c("CLIM")] = clim$OBS[i]
    }

    out = out[which(as.numeric(substr(out$yearmon,1,4)) >= syear & as.numeric(substr(out$yearmon,1,4)) <= eyear), ]
    out = out[order(as.numeric(substr(out$yearmon, 1, 4)), as.numeric(substr(out$yearmon, 6, 7))),]
    out = out[!(is.na(out$MME)) | !(is.na(out$OBS)),]

  }

  return(out)

}

#' @export
BCPoint.Create.Summary.Table <- function(EnvList, mmetype){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  AcuMonths = as.numeric(EnvList$AcuMonths)
  precopt = as.logical(EnvList$precopt)
  CRAdj = as.double(EnvList$CRAdj)

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon

  #prjdir, mmetype, mdlnms, vardir, varfile, syear_mme, AcuMonths, precopt, CRAdj

  #AcuMonths = 1; dsmethod="BCPoint"; precopt=F

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  ltcnt = as.numeric(substr(mmetype, 1, 1))

  # Empty existing summary files
  outdir = paste(prjdir, "/0_Analysis/BCPoint/", mmetype, sep="")
  SetWorkingDir(outdir)
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  monoutdir = paste(prjdir, "/0_Analysis/BCPoint/", mmetype, "/Monthly-TSeries", sep="")
  SetWorkingDir(monoutdir)
  flist = list.files(monoutdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)


  if(AcuMonths > ltcnt){
    AcuMonths = ltcnt
    cat(sprintf("     BCPoint: Accumulation month is greater than max. lead time!"))
  }

  cat("process now running...Create Summary Table",mmetype);cat("\n")
  pb <- txtProgressBar(min=1,max=varcnt,style=3)
  for(i in 1:varcnt){

    varnm = varnms[i]

    ### Combine forecasted output into table and convert unit
    smry = BCPoint.Combine.Forecast.Output(EnvList, mmetype, varnm, lagtime=1)

    if(!is.na(smry)){
      if(varnm == "prec" | precopt == T){ smry = prcp.mmday2mmmon(smry) }

      for(ii in 1:AcuMonths){

        acumon = ii

        colcnt = ncol(smry); rowcnt = nrow(smry)
        acusmry = as.data.frame(array(NA, dim=c((rowcnt - acumon + 1), colcnt)))  #Empty array
        colnames(acusmry) = colnames(smry)
        acusmry$yearmon = smry$yearmon[1:(rowcnt - acumon + 1)]
        for(jj in 1:(rowcnt - acumon + 1)){
          acusmry[jj, 2:colcnt] = colMeans(smry[jj:(jj+acumon-1), 2:colcnt])
        }

        # Define output dir and file
        fname = paste("BCPoint-", varnm, sprintf("-AM%02d",acumon), "-Summary.csv", sep="")
        OutDFile = paste(outdir, "/", fname, sep="")
        write.csv(acusmry, OutDFile, row.names=F)

        for(j in 1:12){

          #ctdata = smry[which(as.numeric(substr(smry$yearmon, 6,7)) == j),c("yearmon", "MME", "OBS")]
          ctdata = acusmry[which(as.numeric(substr(acusmry$yearmon, 6,7)) == j), ]

          colcnt = ncol(ctdata) - 3
          outdata = as.data.frame(ctdata$yearmon); colnames(outdata) = "yearmon"
          for(jj in 2:colcnt){
            colnm = names(ctdata)[jj]
            mdlnm = substr(colnm, 1, (nchar(colnm)-2))
            errdata = ctdata[,c(colnm,"OBS")]
            rsltthold = nrow(errdata[!is.na(errdata$OBS),]) * 0.8
            varthold = nrow(errdata[!is.na(errdata$OBS),]) * 0.5
            errdata = na.omit(errdata)
            colnames(errdata) = c("SIM", "OBS")
            monnm = sprintf("M%02d", j)
            #ltnm = sprintf("LT%02d", k)
            CR = critical.r(nrow(errdata))

            if(nrow(errdata) >= rsltthold){


              cor = format(cor(errdata$SIM, errdata$OBS, method="pearson"), digits=2)

              if(cor > CR * CRAdj){
                outtmp = ctdata[, c("yearmon", colnm)]
                outdata = merge(outdata, outtmp, by="yearmon", all=T)

                #Post-process: Here
              }

            } # End if
          } # Column Loop

          if(ncol(outdata) > 1){

            rcnt = nrow(outdata)
            outdata$MME = NA
            for(kk in 1:rcnt){
              if(!all(is.na(outdata[kk, 2:(ncol(outdata)-1)]))){
                outdata[kk, c("MME")] = mean(as.numeric(outdata[kk, 2:(ncol(outdata)-1)]), na.rm=T)
              }
            }

            obsdata = ctdata[c("yearmon", "OBS")]
            outdata = merge(outdata, obsdata, by="yearmon", all=T)

            climdata = ctdata[c("yearmon", "CLIM")]
            outdata = merge(outdata, climdata, by="yearmon", all=T)

            outfile = paste("BCPoint-", varnm, sprintf("-M%02d", j), sprintf("-AM%02d",acumon), ".csv", sep="")
            MonDFile = paste(monoutdir, "/", outfile, sep="")
            write.csv(outdata, MonDFile, row.names=F)

          } # End IF

        } # Month Loop

      } # Accumulation Loop

    }
    setTxtProgressBar(pb,i)

  } # Variable Loop
  close(pb)
}


#' @export
BCPoint.Get.Common.DataPeriods <- function(mdldir) {

  monnms_all = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  syear = 0
  eyear = 20000

  # Get list of all sub folders
  dirnms = list.dirs(mdldir, recursive = F)
  monnms = matrix(unlist(strsplit(dirnms, "/")), nrow=length(dirnms), byrow=T)
  monnms = monnms[,ncol(monnms)]

  moncnt = length(monnms)
  for(i in 1:moncnt){
    if(monnms[i] %in% monnms_all){
      yeardirnms = list.dirs(dirnms[i], recursive = F)
      yearnms = as.numeric(substr(yeardirnms, nchar(yeardirnms)-3, nchar(yeardirnms)))
      minyear = min(yearnms, na.rm=T)
      if(minyear > syear) syear = minyear
      #if(i>1 & minyear > syear) syear = minyear
      maxyear = max(yearnms, na.rm=T)
      if(maxyear < eyear) eyear = maxyear + 1
    }
  }

  outList = list("syear"=syear, "eyear"=eyear)

  return(outList)

}

#' @export
BCArea.Calculate.SME.Area.Averages <- function(EnvList, mmetype, varnms, pointopt){

  #prjdir, mmedir, mmetype, mdlnms, bnddir, bndfile, varnms, pointopt, weightopt
  prjdir = EnvList$prjdir
  mmedir = EnvList$mmedir
  bnddir = EnvList$bnddir
  bndfile = EnvList$bndfile
  syear = as.numeric(EnvList$syear_mme)
  eyear = as.numeric(substr(Sys.Date(), 1, 4))

  weightopt = F

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon


  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  cmmedir = paste(mmedir, "/", mmetype, sep="")
  outdir = paste(prjdir, "/BCArea/", mmetype, "/01_Area-Average", sep="")
  SetWorkingDir(outdir)

  ### Get model related information
  ltcnt = as.numeric(substr(mmetype,1,1))

  # read boundary shape file (gis file should have "ID" column)
  #GisDFile = paste(bnddir, "/", substr(bndfile, 1, (nchar(bndfile)-4)) , sep="")
  GisDFile = paste(bnddir, "/", bndfile , sep="")
  if(pointopt == 'T') {
    bond = readShapePoints(GisDFile, IDVar="ID")
  } else {
    bond = readShapePoly(GisDFile, IDvar="ID")
  }
  # station count
  stncnt = length(bond$ID)
  stnnms = bond$ID

  varcnt = length(varnms); mdlcnt = length(mdlnms)
  ### Variable Loop (prec, t2m)
  for(i in 1:varcnt){
    varnm = varnms[i]

    ## Model Loop
    for(j in 1:mdlcnt){

      mdlnm = mdlnms[j]

      # get start and end year for selected model
      mdldir = paste(cmmedir, "/", mdlnm, sep="")
      #years = BCPoint.Get.Common.DataPeriods(mdldir)
      #syear = years$syear; eyear = years$eyear

      nmonths = (eyear - syear + 1) * 12
      dat = array(NA, dim=c(stncnt, ltcnt, nmonths))  #Empty array

      ######### Save data into array ############
      #Dummy Loop for year folders
      cnt = 1
      for(k in syear:eyear){
        curyear = k

        # Loop for month folders
        for(ii in 1:12){
          monnm = monnms[ii]

          wdir = paste(cmmedir, "/", mdlnm, "/", monnm, "/", curyear, sep="")
          ncfile = paste(varnm, ".nc", sep="")

          NcDFile = paste(wdir, "/", ncfile, sep="")
          if(file.exists(NcDFile)){
            nc = ReadNetCDF4(wdir, ncfile)
            var = nc$var; x = nc$x; y = nc$y # [lat, lon, ESM, LTime]

            for(jj in 1:ltcnt){
              varg = apply(var[,,, jj], 1:2, mean)
              r = GIS.ncVar2raster(x, y, varg)
              v = extract(r, bond, fun=mean, weights=weightopt, small=T)

              for(m in 1:stncnt){
                dat[m,jj,cnt] = v[m] # dat[stncnt, ltcnt, nmonths]
              }

            }
          } else {
            #dat[,jj,cnt] = NA
            dat[ , ,cnt] = NA
          }

          cnt = cnt+1

          cat(sprintf("     BCArea: Area average are calculated: Variable=%s Model=%s  Year=%s  Month=%s\n", varnm, mdlnm, curyear, ii))
        } # Loop for month
      } # Loop for year


      for(i in 1:stncnt){
        stnnm = stnnms[i]

        for(j in 1:ltcnt){

          #### Define yearmon information
          smon = j; emon = 11+j
          if(emon >= 13){
            emon2 = emon  -12
            eyear2 = eyear + 1
          } else {
            emon2 = emon
            eyear2 = eyear
          }

          sdate = as.Date(paste(syear, "-", smon, "-01", sep=""))
          edate = as.Date(paste(eyear2, "-", emon2, "-01", sep=""))
          # Get yearmon infomration
          yearmon = substr(seq(sdate, edate, by="mon"),1,7)

          # dat[stncnt, ltcnt, nmonths]
          val = as.data.frame(dat[i,j,])

          out = cbind(yearmon, val)
          colnames(out) = c("yearmon", varnm)

          #### Write file
          LTime = sprintf("LT%02d",j)
          ofname = paste(outdir, "/", stnnm, "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
          write.csv(out, ofname, row.names=F)
        }

      }
      cat(sprintf("     BCArea: Values are extracted: Variable=%s  Model=%s\n", varnm, mdlnm))

    } ## End of Model Loop
  }  ### End of Variable Loop

}

#' @export
BCArea.SME.Bias.Correction <- function(EnvList, mmetype){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile

  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  ##### Define folder names
  indir = paste(prjdir, "/BCArea/", mmetype, "/01_Area-Average", sep="")
  outdir = paste(prjdir, "/BCArea/", mmetype, "/02_BiasCorrected", sep="")
  SetWorkingDir(outdir)


  ###### Get Station ID, lat, and Lon information
  flist = list.files(indir, pattern = glob2rx("*.csv"), full.names = F)
  stnnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,1])
  varnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,2])
  mdlnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,3])
  ltnms = substr(unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,4]),1,4)
  ltcnt = length(ltnms)

  # Count for Loops
  stncnt = length(stnnms); varcnt = length(varnms); mdlcnt = length(mdlnms)

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
        fname = paste(indir, "/", stnid, "-", varnm, "-", mdlnm, "-LT01.csv", sep="")
        dat = read.csv(fname, header=T)
        years = as.numeric(substr(dat$yearmon, 1, 4))
        syear = min(years)
        eyear = max(years)

        # read obs data
        ObsDFile = paste(vardir, "/", varfile, sep="")
        obsdata = read.csv(ObsDFile, header=T)

        # select common period
        obsdata = obsdata[which(substr(obsdata$yearmon, 1,4) >= syear & substr(obsdata$yearmon, 1,4) <= eyear),]

        obsdata = obsdata[,c("yearmon", stnid)]
        colnames(obsdata) = c("yearmon", "obs")
        obsdata$mon = as.numeric(substr(obsdata$yearmon, 6, 7))

        # Lead-Time Loop
        for(ii in 1:ltcnt){

          LTime = sprintf("LT%02d",ii)

          outname = paste(outdir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, ".-BCcsv", sep="")
          if(!file.exists(outname)){

            mmename = paste(indir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, ".csv", sep="")
            mmedata = read.csv(mmename, header=T)


            colnames(mmedata) = c("yearmon", "mme")
            mmedata$mon = as.numeric(substr(mmedata$yearmon, 6, 7))

            ############ Calculae QMapping fit
            for(kk in 1:12){

              obsmon = obsdata[which(obsdata$mon == kk),c("yearmon", "obs")]
              mmemon = mmedata[which(mmedata$mon == kk),c("yearmon", "mme")]

              if(all(is.na(mmemon$mme))){
                mmetmp = mmemon
                colnames(mmetmp) = c("yearmon", varnm)

              } else {

                comdata = merge(obsmon, mmemon)
                comdata = na.omit(comdata)

                # Use the overall avearage if there is no common period
                if(nrow(comdata) < 20){
                  obsmean = mean(obsmon$obs, na.rm=T)
                  mmemean = mean(mmemon$mme, na.rm=T)
                }else{
                  obsmean = mean(comdata$obs)
                  mmemean = mean(comdata$mme)
                }

                # Do bias correction
                if(varnm == "t2m"){
                  mmetmp = mmemon
                  mmetmp$mme = obsmean + (mmemon$mme - mmemean)
                  colnames(mmetmp) = c("yearmon", varnm)
                }else {
                  mmetmp = mmemon
                  # in case greater than average, add
                  mmetmp$mme[which((mmetmp$mme - mmemean) >= 0)] = obsmean + (mmetmp$mme[which((mmetmp$mme - mmemean) >= 0)] - mmemean)
                  # in case less than avearage, use ratio
                  mmetmp$mme[which((mmemon$mme - mmemean) < 0)] = obsmean / mmemean * mmetmp$mme[which((mmemon$mme - mmemean) < 0)]
                  colnames(mmetmp) = c("yearmon", varnm)
                }

              }  ## End if

              if(kk == 1) {
                mmeadj = mmetmp
              }else {
                mmeadj = rbind(mmeadj, mmetmp)
              }

            } ## End of kk Loop (month)

            #colnames(mmeadj) = c("yearmon",varnm)
            mmeadj = mmeadj[order(mmeadj$yearmon),]
            colnames(mmeadj) = c("yearmon", varnm)


            outname = paste(outdir, "/", stnid, "-", varnm, "-", mdlnm, "-", LTime, "-BC.csv", sep="")
            write.csv(mmeadj, outname, row.names=F)

          } else {

            cat(sprintf("     BCArea: File already exists! : File=%s\n", outname))
          }

        } # End of Lead Time Loop

        cat(sprintf("     BCArea: Bias Correction has been finished: Station=%s  Variable=%s  Model=%s\n", stnid, varnm, mdlnm))

      } ## End of Model Loop
    }  ### End of Variable Loop
  }   #### End of Station Loop

}

#' @export
BCArea.Create.Summary.Table <- function(EnvList, mmetype, AcuMonths, precopt, CRAdj){

  #AcuMonths = 1; precopt=F; CRAdj=1.0
  #prjdir, mmetype, mdlnms, vardir, varfile, syear_mme, AcuMonths, precopt, CRAdj
  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon


  options(warn=-1)
  options(stringsAsFactors = FALSE)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  ltcnt = as.numeric(substr(mmetype, 1, 1))

  # Empty existing summary files
  outdir = paste(prjdir, "/0_Analysis/BCArea/", mmetype, sep="")
  SetWorkingDir(outdir)
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  monoutdir = paste(prjdir, "/0_Analysis/BCArea/", mmetype, "/Monthly-TSeries", sep="")
  SetWorkingDir(monoutdir)
  flist = list.files(monoutdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)


  if(AcuMonths > ltcnt){
    AcuMonths = ltcnt
    cat(sprintf("     BCArea: Accumulation month is greater than max. lead time!"))
  }

  for(i in 1:varcnt){

    varnm = varnms[i]

    ### Combine forecasted output into table and convert unit
    smry = BCArea.Combine.Forecast.Output(EnvList, mmetype, varnm, lagtime=1)

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
        fname = paste("BCArea-", varnm, sprintf("-AM%02d",acumon), "-Summary.csv", sep="")
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

            outfile = paste("BCArea-", varnm, sprintf("-M%02d", j), sprintf("-AM%02d",acumon), ".csv", sep="")
            MonDFile = paste(monoutdir, "/", outfile, sep="")
            write.csv(outdata, MonDFile, row.names=F)

          } # End IF

        } # Month Loop

      } # Accumulation Loop

    }

  } # Variable Loop

}

#' @export
BCArea.Combine.Forecast.Output <- function(EnvList, mmetype, varnm, lagtime){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  mdlnms_3mon = EnvList$mdlnms_3mon
  mdlnms_6mon = EnvList$mdlnms_6mon

  if(mmetype == "3MON") mdlnms = EnvList$mdlnms_3mon
  if(mmetype == "6MON") mdlnms = EnvList$mdlnms_6mon

  syear = syear_mme
  eyear = as.numeric(substr(Sys.Date(),1,4))+1
  #eyear = eyear_obs

  # define working directory
  bcadir = paste(prjdir, "/BCArea/", mmetype, "/02_BiasCorrected", sep="")

  mdlcnt = length(mdlnms)

  # Max LeadTime from mmetype
  ltcnt = as.numeric(substr(mmetype,1,1))
  cnt = 1
  out = NA
  for(i in lagtime:ltcnt){

    for(j in 1:mdlcnt){

      mdlnm = mdlnms[j]

      srchstr = paste(varnm, "-*", mdlnm, sprintf("-LT%02d",i), "-BC.csv", sep="")
      OutDFile = list.files(bcadir, pattern = glob2rx(srchstr), full.names = T)

      data = read.csv(OutDFile, header=T)
      colnames(data) = c("yearmon", sprintf("%s%02d",mdlnm,i))

      if(cnt == 1){
        out = data
        cnt = cnt + 1
      } else {
        tmp = data
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


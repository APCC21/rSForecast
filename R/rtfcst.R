#' @export
RTFcst.Create.Model.Integration.Table <- function(EnvList, AcuMonths){

  prjdir = EnvList$prjdir
  mdlnms_3mon = EnvList$mdlnms_3mon
  mdlnms_6mon = EnvList$mdlnms_6mon
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  BCAreaOpt = EnvList$BCAreaOpt
  BCPointOpt = EnvList$BCPointOpt
  CIRegOpt = EnvList$CIRegOpt
  MWRegOpt = EnvList$MWRegOpt
  MWRObsOpt = EnvList$MWRObsOpt

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  ################## Create combined table using 3MON and 6MON simulations
  if(AcuMonths > 3){
    AcuMonths = 3
    cat(sprintf("     RTFcst: Accumulation month is greater than max. lead time!"))
  }

  for(i in 1:varcnt){
    varnm = varnms[i]

    for(ii in 1:AcuMonths){
      acumon = ii

      cnt = 1

      #BCPoint
      if(BCPointOpt == "On") {
        for(mmetype in c("3MON", "6MON")){
          ltcnt = as.numeric(substr(mmetype, 1, 1))
          sbcdir = paste(prjdir, "/0_Analysis/BCPoint/", mmetype, sep="")

          SbcDFile = paste(sbcdir, "/BCPoint-", varnm, sprintf("-AM%02d",acumon), "-TCC.csv", sep="")

          if(file.exists(SbcDFile)){
            sbcdata = read.csv(SbcDFile, header=T)
            sbcdata = sbcdata[!(colnames(sbcdata) %in% c("MME"))]
            colnames(sbcdata) = c("month", paste("B_", colnames(sbcdata[2:ncol(sbcdata)]), sep=""))
          } else {
            month = "M01"; sbcdata = data.frame(month)
          }

          if(mmetype == "3MON"){
            sbcout = sbcdata
          } else {
            sbctmp = sbcdata
            sbctmp = sbctmp[!(colnames(sbctmp) %in% colnames(sbcout[2:ncol(sbcout)]))]
            sbcout = merge(sbcout, sbctmp, by="month", all=T )
          }
        }
        OutTbl = sbcout
        cnt = cnt + 1
      }

      #BCArea
      if(BCAreaOpt == "On"){
        for(mmetype in c("3MON", "6MON")){
          ltcnt = as.numeric(substr(mmetype, 1, 1))
          sbcdir = paste(prjdir, "/0_Analysis/BCArea/", mmetype, sep="")

          SbcDFile = paste(sbcdir, "/BCArea-", varnm, sprintf("-AM%02d",acumon), "-TCC.csv", sep="")

          if(file.exists(SbcDFile)){
            sbcdata = read.csv(SbcDFile, header=T)
            sbcdata = sbcdata[!(colnames(sbcdata) %in% c("MME"))]
            colnames(sbcdata) = c("month", paste("B_", colnames(sbcdata[2:ncol(sbcdata)]), sep=""))
          } else {
            month = "M01"; sbcdata = data.frame(month)
          }

          if(mmetype == "3MON"){
            sbcout = sbcdata
          } else {
            sbctmp = sbcdata
            sbctmp = sbctmp[!(colnames(sbctmp) %in% colnames(sbcout[2:ncol(sbcout)]))]
            sbcout = merge(sbcout, sbctmp, by="month", all=T )
          }
        }
        if(cnt == 1) {
          OutTbl = sbcout
        } else {
          cat("     RTFcst: Both BCPoint and BCArea are turned on\n")
        }
        cnt = cnt + 1
      }

      # CIReg
      if(CIRegOpt == "On"){
        cirdir = paste(prjdir, "/0_Analysis/CIReg/", sep="")
        CirDFile = paste(cirdir, "/CIReg-", varnm, sprintf("-AM%02d",acumon), "-TCC.csv", sep="")

        if(file.exists(CirDFile)){
          cirdata = read.csv(CirDFile, header=T)
          cirdata = cirdata[!(colnames(cirdata) %in% c("MME"))]
          colnames(cirdata) = c("month", paste("C_", colnames(cirdata[2:ncol(cirdata)]), sep=""))
        } else {
          month = "M01"; cirdata = data.frame(month)
        }
        if(cnt == 1){
          OutTbl = cirdata
        } else {
          OutTbl = merge(OutTbl, cirdata, by="month", all=T)
        }

        cnt = cnt + 1
      }

      # MWReg
      if(MWRegOpt == "On"){
        mwrcnt = 0; month = "M01"; mwrout = data.frame(month)
        for(mmetype in c("3MON", "6MON")){
          ltcnt = as.numeric(substr(mmetype, 1, 1))
          mwrdir = paste(prjdir, "/0_Analysis/MWReg/", mmetype, sep="")
          MwrDFile = paste(mwrdir, "/MWReg-", varnm, sprintf("-AM%02d",acumon), "-TCC.csv", sep="")

          if(file.exists(MwrDFile)){
            mwrdata = read.csv(MwrDFile, header=T)
            mwrdata = mwrdata[!(colnames(mwrdata) %in% c("MME"))]
            colnames(mwrdata) = c("month", paste("M_", colnames(mwrdata[2:ncol(mwrdata)]), sep=""))
          } else {
            month = "M01"; mwrdata = data.frame(month)
          }

          #if(mmetype == "3MON" | (mmetype == "6MON" & nrow(mwrdata) == 1)){
          if(ncol(mwrdata) > 1 & mwrcnt == 0){
            mwrout = mwrdata
            mwrcnt = mwrcnt + 1
          } else if (mwrcnt > 0) {
            mwrtmp = mwrdata
            mwrtmp = mwrtmp[!(colnames(mwrtmp) %in% colnames(mwrout[2:ncol(mwrout)]))]
            mwrout = merge(mwrout, mwrtmp, by="month", all=T )
          }
        }

        if(cnt == 1){
          OutTbl = mwrout
        } else {
          OutTbl = merge(OutTbl, mwrout, by="month", all=T)
        }
        cnt = cnt + 1
      }

      # MWRObs
      if(MWRObsOpt == "On"){
        mwodir = paste(prjdir, "/0_Analysis/MWRObs/", sep="")
        MwoDFile = paste(mwodir, "/MWRObs-", varnm, sprintf("-AM%02d",acumon), "-TCC.csv", sep="")

        if(file.exists(MwoDFile)){
          mwodata = read.csv(MwoDFile, header=T)
          mwodata = mwodata[!(colnames(mwodata) %in% c("MME"))]
          colnames(mwodata) = c("month", paste("O_", colnames(mwodata[2:ncol(mwodata)]), sep=""))
        } else {
          month = "M01"; mwodata = data.frame(month)
        }
        mwoout = mwodata

        if(cnt == 1 & nrow(mwoout) > 1) {
          OutTbl = mwoout
        } else if(cnt > 1 & nrow(mwoout) > 1) {
          OutTbl = merge(OutTbl, mwoout, by="month", all=T)
        }
        cnt = cnt + 1
      }

      #outdir = paste(prjdir, "/0_RTForecast/ALL", sep="")
      outdir = paste(prjdir, "/0_RTForecast", sep="")
      SetWorkingDir(outdir)
      OutDFile = paste(outdir, "/RTFcst-Table-", varnm, sprintf("-AM%02d",acumon), ".csv", sep="")
      write.csv(OutTbl, OutDFile, row.names=F)

    }

  }


}

#' @export
RTFcst.Create.Summary.Table <- function(EnvList, precopt){

  prjdir = EnvList$prjdir
  mdlnms_3mon = EnvList$mdlnms_3mon
  mdlnms_6mon = EnvList$mdlnms_6mon
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  BCAreaOpt = EnvList$BCAreaOpt
  BCPointOpt = EnvList$BCPointOpt
  CIRegOpt = EnvList$CIRegOpt
  MWRegOpt = EnvList$MWRegOpt
  MWRObsOpt = EnvList$MWRObsOpt

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  outdir = paste(prjdir, "/0_RTForecast/Monthly-TSeries/Original", sep="")
  SetWorkingDir(outdir)
  # reset files
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  outdir2 = paste(prjdir, "/0_RTForecast/Monthly-TSeries", sep="")
  SetWorkingDir(outdir2)
  # reset files
  flist = list.files(outdir2, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  ############## Combined 3MON and 6MON
  ltcnt = 6
  for(i in 1:varcnt){
    varnm = varnms[i]

    for(m in 1:ltcnt){

      lagtime = m

      #for(ii in 1:2){ # Forecasting Mode (1: Continous, 2: one-time)
      for(ii in 1:1){

        acumon = 1
        FcstDFile = paste(prjdir, "/0_RTForecast/RTFcst-Table-", varnm, sprintf("-AM%02d",acumon), ".csv", sep="")
        fcstbl_all = read.csv(FcstDFile, header=T)

        fcstbl = NA
        if(ii == 1){

          for(j in m:ltcnt){
            tmp = as.data.frame(fcstbl_all[, c(1, grep(sprintf("%02d", j), colnames(fcstbl_all)))])
            if(ncol(tmp) > 1){
              if(j == m | is.na(fcstbl)){
                fcstbl = tmp
              } else if(!is.na(fcstbl)) {
                fcstbl = merge(fcstbl, tmp, by="month")
              }
            } else if(ncol(tmp) == 1 & is.na(fcstbl)) {
              fcstbl = NA
            }
          }

        } else {

          fcstbl = fcstbl_all[, c(1, grep(sprintf("%02d", lagtime), colnames(fcstbl_all)))]

        }

        if(!is.na(fcstbl)){
          fcstcnt = nrow(fcstbl)
          for(k in 1:fcstcnt){

            # current month
            curmon = as.numeric(substr(fcstbl[k, c("month")], 2, 3))

            tmp = as.data.frame(t(fcstbl[k, 2:ncol(fcstbl)])); colnames(tmp) = "TCC"
            if(nrow(tmp) == 1){
              tmp$colnm = colnames(fcstbl)[2]
            } else {
              tmp$colnm = row.names(tmp)
            }
            tmp = na.omit(tmp)
            #tmp$colnm = row.names(tmp)
            rcnt = nrow(tmp)

            if(rcnt > 0){
              cnt = 1
              for(kk in 1:rcnt){

                orgcolnm = tmp$colnm[kk]
                colnm = substr(orgcolnm, 3, nchar(orgcolnm))
                curltime = as.numeric(substr(colnm, (nchar(colnm)-1), nchar(colnm)))
                dsm = substr(orgcolnm, 1, 2)

                if(dsm == "B_" & BCPointOpt == "On"){
                  dsmethod="BCPoint"
                  dir3mon = paste(prjdir, "/0_Analysis/BCPoint/3MON/Monthly-TSeries", sep="")
                  dir6mon = paste(prjdir, "/0_Analysis/BCPoint/6MON/Monthly-TSeries", sep="")
                }

                if(dsm == "B_" & BCAreaOpt ==  "On"){
                  dsmethod="BCArea"
                  dir3mon = paste(prjdir, "/0_Analysis/BCArea/3MON/Monthly-TSeries", sep="")
                  dir6mon = paste(prjdir, "/0_Analysis/BCArea/6MON/Monthly-TSeries", sep="")
                }

                if(dsm == "C_" & CIRegOpt ==  "On"){
                  dsmethod="CIReg"
                  dir3mon = paste(prjdir, "/0_Analysis/CIReg/Monthly-TSeries", sep="")
                  dir6mon = paste(prjdir, "/0_Analysis/CIReg/Monthly-TSeries", sep="")
                }

                if(dsm == "M_" & MWRegOpt ==  "On"){
                  dsmethod="MWReg"
                  dir3mon = paste(prjdir, "/0_Analysis/MWReg/3MON/Monthly-TSeries", sep="")
                  dir6mon = paste(prjdir, "/0_Analysis/MWReg/6MON/Monthly-TSeries", sep="")
                }

                if(dsm == "O_" & MWRObsOpt ==  "On"){
                  dsmethod="MWRObs"
                  dir3mon = paste(prjdir, "/0_Analysis/MWRObs/Monthly-TSeries", sep="")
                  dir6mon = paste(prjdir, "/0_Analysis/MWRObs/Monthly-TSeries", sep="")
                }


                InDFile = paste(dir3mon, "/", dsmethod, "-", varnm, "-", sprintf("M%02d",curmon), sprintf("-AM%02d",acumon), ".csv", sep="")
                if(file.exists(InDFile) & curltime <= 3){
                  data = read.csv(InDFile, header=T)
                  if(length(grep(colnm, colnames(data))) == 0 ){
                    InDFile = paste(dir6mon, "/", dsmethod, "-", varnm, "-", sprintf("M%02d",curmon), sprintf("-AM%02d",acumon), ".csv", sep="")
                    if(file.exists(InDFile)){
                      data = read.csv(InDFile, header=T)
                    }
                  }
                } else {
                  InDFile = paste(dir6mon, "/", dsmethod, "-", varnm, "-", sprintf("M%02d",curmon), sprintf("-AM%02d",acumon), ".csv", sep="")
                  if(file.exists(InDFile)){
                    data = read.csv(InDFile, header=T)
                  }
                }

                if(length(grep(colnm, colnames(data))) > 0 ){

                  if(cnt == 1){
                    outdata = data[, c("yearmon", colnm)]
                    colnames(outdata) = c("yearmon", orgcolnm)
                  } else {
                    outtmp = data[, c("yearmon", colnm)]
                    colnames(outtmp) = c("yearmon", orgcolnm)
                    outdata = merge(outdata, outtmp, by="yearmon", all=T)
                  }
                  cnt = cnt + 1

                }

              }

              # There can be missing year, and merging using all=T
              enddate = seq(Sys.Date(), by = "month", length = ltcnt)[ltcnt]
              yearmon = substr(seq(as.Date(paste(syear_mme, "-01-01", sep="")), enddate, by = "month") ,1,7)
              yearmon = data.frame(yearmon[which(as.numeric(substr(yearmon, 6, 7)) == curmon)])
              colnames(yearmon) = c("yearmon")
              outdata = merge(yearmon, outdata, by="yearmon", all=T)

              if(ncol(outdata) > 1){

                rcnt2 = nrow(outdata)
                outdata$MME = NA
                for(kk in 1:rcnt2){
                  if(!all(is.na(outdata[kk, 2:(ncol(outdata)-1)]))){
                    outdata[kk, c("MME")] = mean(as.numeric(outdata[kk, 2:(ncol(outdata)-1)]), na.rm=T)
                  }
                }

                obsdata = data[c("yearmon", "OBS")]
                outdata = merge(outdata, obsdata, by="yearmon", all=T)

                climdata = data[c("yearmon", "CLIM")]
                outdata = merge(outdata, climdata, by="yearmon", all=T)

                if(ii == 1){
                  outfile = paste("RTFcst-", varnm, sprintf("-M%02d", curmon), "-Cont", sprintf("-LT%02d", m), ".csv", sep="")

                } else {
                  outfile = paste("RTFcst-", varnm, sprintf("-M%02d", curmon), "-Once", sprintf("-LT%02d", m), ".csv", sep="")

                }

                #Post-process: Insert function Here
                outdata2 = RTFcst.Adjust.Mean.SD (outdata, varnm)

                MonDFile = paste(outdir, "/", outfile, sep="")
                write.csv(outdata, MonDFile, row.names=F)

                MonDFile2 = paste(outdir2, "/", outfile, sep="")
                write.csv(outdata2, MonDFile2, row.names=F)

              } # End IF

            } # If Selected models is more than

          } # Month Loop

        }

      } # Accumulation Loop

    } # Lead-time loop

  } # Variable Loop


}

#' @export
RTFcst.Calculate.Error.Statistics <- function(EnvList) {

  prjdir = EnvList$prjdir

  outdir = paste(prjdir, "/0_RTForecast",  sep="")
  # reset files
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  indir = paste(prjdir, "/0_RTForecast/Monthly-TSeries", sep="")
  flist = list.files(indir, pattern = glob2rx("RTFcst*.csv"), full.names = F)

  if(length(flist) > 0) {
    varnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,2])
    monnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,3])
    monnms =  monnms[order(monnms)]
    acunms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,4])
    ltnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,5])
    ltnms = substr(ltnms, 1, nchar(ltnms)-4)
    ltcnt = length(ltnms)

    varcnt = length(varnms); acucnt = length(acunms)

    for(i in 1:varcnt){
      varnm  = varnms[i]

      for(j in 1:acucnt){
        acunm = acunms[j]

        for(k in 1:ltcnt){
          ltnm = ltnms[k]
          cnt = 1
          for(ii in 1:12){

            #monnm = monnms[ii]
            monnm = sprintf("M%02d", ii)

            srchstr = paste("RTFcst-", varnm, "-", monnm, "-", acunm, "-", ltnm, ".csv", sep="")
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
                  nrmsetbl = c(monnm, nrmse)
                } else {
                  cortbl = c(cortbl, cor)
                  nrmsetbl = c(nrmsetbl, nrmse)
                }

              } # Column Loop

              cortbl = t(data.frame(cortbl))
              colnames(cortbl) = c("month", names(data)[2:colcnt])

              nrmsetbl = t(data.frame(nrmsetbl))
              colnames(nrmsetbl) = c("month", names(data)[2:colcnt])

              if(cnt == 1){
                corerr = as.data.frame(cortbl)
                nrmseerr = as.data.frame(nrmsetbl)
              } else {
                cortmp = as.data.frame(cortbl)
                corerr = rbind.fill(corerr, cortmp)

                nrmsetmp = as.data.frame(nrmsetbl)
                nrmseerr = rbind.fill(nrmseerr, nrmsetmp)
              }
              cnt = cnt + 1

            } else{
              if(cnt == 1){
                month = monnm; MME = NA; corerr = data.frame(month, MME)
                month = monnm; MME = NA; nrmseerr = data.frame(month, MME)
              } else {
                month = monnm; MME = NA; cortmp = data.frame(month, MME)
                corerr = rbind.fill(corerr, cortmp)

                month = monnm; MME = NA; nrmsetmp = data.frame(month, MME)
                nrmseerr = rbind.fill(nrmseerr, nrmsetmp)
              }
              cnt = cnt + 1
            } # IF file.exists

          } # Month Loop

          # Change column order
          valdata = corerr[2:ncol(corerr)]
          valdata = valdata[, order(names(valdata))]
          mondata = corerr[1]
          corerr = cbind(mondata, valdata)

          OutDFile = paste(outdir, "/RTFcst-", varnm, "-", acunm, "-", ltnm, "-TCC.csv", sep="")
          write.csv(corerr, OutDFile, row.names=F)


          valdata = nrmseerr[2:ncol(nrmseerr)] # cor --> nrmse
          valdata = valdata[, order(names(valdata))]
          mondata = nrmseerr[1]
          nrmseerr = cbind(mondata, valdata)

          OutDFile = paste(outdir, "/RTFcst-", varnm, "-", acunm, "-", ltnm, "-NRMSE.csv", sep="")
          write.csv(nrmseerr, OutDFile, row.names=F)

        } # Lead Time
      } # Accumulation Months
    } # Variable

    ###### Climatology
    for(i in 1:varcnt){
      varnm  = varnms[i]
      for(ii in 1:12){

        #monnm = monnms[ii]
        monnm = sprintf("M%02d", ii)

        srchstr = paste("RTFcst-", varnm, "-", monnm, "-*.csv", sep="")
        InDFile = list.files(indir, pattern = glob2rx(srchstr), full.names = T)[1]

        if(!is.na(InDFile)){
          data = read.csv(InDFile, header=T)

          errdata = data[,c("CLIM","OBS")]
          errdata = na.omit(errdata)
          colnames(errdata) = c("SIM", "OBS")

          if(nrow(errdata) == 0){
            cor = NA; nrmse = NA
          } else {
            cor = format(cor(errdata$SIM, errdata$OBS, method="pearson"), digits=2)
            nrmse = format(nrmse(errdata$SIM, errdata$OBS, norm="sd")/100, digits=2)
          } # End if

          if(ii == 1){
            corerr = as.data.frame(cor)
            nrmseerr = as.data.frame(nrmse)
          } else {
            cortmp = as.data.frame(cor)
            corerr = rbind(corerr, cortmp)

            nrmsetmp = as.data.frame(nrmse)
            nrmseerr = rbind(nrmseerr, nrmsetmp)
          }
        } else {
          if(ii == 1){
            month = monnm; corerr = as.data.frame(NA); colnames(corerr) =c("cor")
            month = monnm; nrmseerr = as.data.frame(NA); colnames(nrmseerr) =c("nrmse")
          } else {
            month = monnm; cortmp = as.data.frame(NA); colnames(cortmp) =c("cor")
            corerr = rbind(corerr, cortmp)

            month = monnm; nrmsetmp = as.data.frame(NA); colnames(nrmsetmp) =c("nrmse")
            nrmseerr = rbind(nrmseerr, nrmsetmp)
          }
        }

      } # Month Loop

      # Change column order

      OutDFile = paste(outdir, "/RTFcst-", varnm, "-CLIM-TCC.csv", sep="")
      write.csv(corerr, OutDFile, row.names=F)


      OutDFile = paste(outdir, "/RTFcst-", varnm, "-CLIM-NRMSe.csv", sep="")
      write.csv(nrmseerr, OutDFile, row.names=F)


    } # Variable


  }

}

#' @export
RTFcst.Create.Summary.Graph.Table <- function(EnvList, fiyearmon){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)
  nrange = as.numeric(EnvList$nrange)

  options(warn=-1)

  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  fsyearmon = substr(seq(as.Date(paste(fiyearmon, "-01", sep="")), by = "month", length = 2)[2],1,7)
  fsyear = as.numeric(substr(fsyearmon, 1, 4))
  fsmon = monnms[as.numeric(substr(fsyearmon, 6, 7))]

  nmonth = 6

  feyearmon = substr(seq(as.Date(paste(fsyearmon, "-01", sep="")), by = "month", length = nmonth)[nmonth],1,7)
  feyear = as.numeric(substr(feyearmon, 1, 4))
  femon = monnms[as.numeric(substr(feyearmon, 6, 7))]


  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  for(i in 1:varcnt){

    varnm = varnms[i]

    smry = RTFcst.Combine.Forecast.Output (EnvList, varnm, lagtime=1)

    if(!is.na(smry)){

      ###############################################################
      ### Probability forecast for future period
      temp = smry[which(smry$yearmon == fsyearmon):nrow(smry), ]
      yearmon = as.data.frame(substr(seq(as.Date(paste(fsyearmon, "-01", sep="")), as.Date(paste(feyearmon, "-01", sep="")), by= "months"), 1,7))
      colnames(yearmon) = "yearmon"
      fdata = merge(yearmon, temp, all = T)
      fdata = fdata[1:nmonth, !(colnames(fdata) %in% c("MME", "OBS", "CLIM"))]


      for(j in 1:nrow(fdata)){

        cmonth = as.numeric(substr(fdata[j,"yearmon"], 6, 7))

        # get range values
        rngdir = paste(prjdir, "/0_RTForecast/ProbabilityRanges", sep="")
        RngDFile = paste(rngdir,  "/", varnm, "-probability_range.csv", sep="")
        rng = read.csv(RngDFile, header=T)
        rngvals = as.vector(t(rng[which(rng$month == cmonth), 2:ncol(rng)]))


        vals = as.vector(na.omit(t(fdata[j,2:ncol(fdata)])))

        if(j == 1){
          prob = t(Calculate.Probability(nrange, vals, rngvals))
        } else {
          temp = t(Calculate.Probability(nrange, vals, rngvals))
          prob = cbind(prob, temp)
        }
      }

      colnames(prob) = fdata$yearmon
      if(nrange == 3){
        rownames(prob) = c("Low", "Normal", "High")
      }
      if(nrange == 5){
        rownames(prob) = c("Very Low", "Low", "Normal", "High", "Very High")
      }

      prob = format(prob, digits=1)
      colnames(prob) = monnms[as.numeric(substr(colnames(prob),6,7))]

      for(k in 1:nmonth){

        smry = RTFcst.Combine.Forecast.Output (EnvList, varnm, lagtime=k)

        if(!is.na(smry)){
          sdata = na.omit(smry[, c("yearmon", "MME", "OBS")])
          if(k == 1){
            ssyear = as.numeric(substr(sdata$yearmon[1], 1, 4))
            ssmon = monnms[as.numeric(substr(sdata$yearmon[1], 6, 7))]
            seyear = as.numeric(substr(sdata$yearmon[nrow(sdata)], 1, 4))
            semon = monnms[as.numeric(substr(sdata$yearmon[nrow(sdata)], 6, 7))]
          }

          for(j in 1:12){

            rngvals = as.vector(t(rng[which(rng$month == j), 2:ncol(rng)]))

            ctdata = sdata[which(as.numeric(substr(sdata$yearmon, 6,7)) == j),c("yearmon", "MME", "OBS")]

            if(nrow(ctdata) == 0){
              cor = "-"; nrmse = "-"; pc = "-"; hss = "-"
            } else {

              ctdata = ctdata[c("MME", "OBS")]
              colnames(ctdata) = c("SIM", "OBS")

              # Create contingency table for each month
              CTable = Create.Contingency.Table(rngvals, ctdata)

              cor = format(cor(ctdata$SIM, ctdata$OBS), digits=2)
              nrmse = format(nrmse(ctdata$SIM, ctdata$OBS, norm="sd")/100, digits=2)


              SScore = multi.cont(CTable, baseline = NULL)
              pc = format(SScore$pc, digits=2)
              hss = format(SScore$hss, digits=2)

            }

            if(j == 1){
              #sc = rbind(cor, nrmse, pc, hss)
              sc = rbind(cor, hss)
              colnames(sc) = monnms[j]
            } else {
              #tmp = rbind(cor, nrmse, pc, hss)
              tmp = rbind(cor, hss)
              colnames(tmp) = monnms[j]
              sc = cbind(sc, tmp)
            }

          }
          #rnms = c("TCC", "NRMSE", "Accuracy", "HSS")
          rnms = c("TCC", "HSS")

          if(k == 1){
            sctbl = sc
            rownames(sctbl) = paste("LT", k, "-", rnms, sep="")
          } else {
            tmptbl = sc
            rownames(tmptbl) = paste("LT", k, "-", rnms, sep="")
            sctbl = rbind(sctbl, tmptbl)
          }
        }



      }



      ################################################################
      ### subset of df for graph

      smry = RTFcst.Combine.Forecast.Output (EnvList, varnm, lagtime=1)

      gsyear = feyear - 1
      temp = smry[which(as.numeric(substr(smry$yearmon, 1,4)) >= gsyear & as.numeric(substr(smry$yearmon, 1,4)) <= feyear), ]
      yearmon = as.data.frame(substr(seq(as.Date(paste(gsyear, "-01-01", sep="")), as.Date(paste(feyearmon, "-01", sep="")), by= "months"), 1,7))
      colnames(yearmon) = "yearmon"
      gdata = merge(yearmon, temp, all = T)
      gdata[which(as.numeric(substr(gdata$yearmon, 1,4)) == feyear & as.numeric(substr(gdata$yearmon, 6,7)) > as.numeric(substr(feyearmon, 6, 7))), 2:ncol(gdata)] = NA
      gdata$yearmon = factor(gdata$yearmon)

      ymax = max(gdata[,2:ncol(gdata)], na.rm=T)*1.2
      ymin = min(gdata[,2:ncol(gdata)], na.rm=T)*1.2

      if(varnm == "t2m"){
        gdata[,2:ncol(gdata)] = gdata[,2:ncol(gdata)] - gdata[,c("CLIM")]
        ymin = -4; ymax = 4 # Fix min and max values for temperature
      }

      ldata = gdata[,c("yearmon", "MME", "OBS", "CLIM")]
      bdata = gdata[-which(names(gdata) %in% c("MME", "OBS", "CLIM"))]

      bd = na.omit(melt(bdata, id=c("yearmon")))


      #pngfile = paste(prjdir, "/0_RTForecast/", mmetype, "/RTFcst-", varnm, "-", fiyearmon, ".png", sep="")
      pngfile = paste(prjdir, "/0_RTForecast/RTFcst-", varnm, "-", fiyearmon, ".png", sep="")
      png(pngfile, width = 11, height = 8, units = 'in', res = 400)

      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),  widths=c(3,2), heights=c(3,2))

      # 1st panel
      #par(mar = c(2,2,4,2))
      if(varnm == "t2m"){
        #mtitle = paste(substr(mmetype,1,1), "-month Anomaly Temperature (2m) Forecast for ", fsmon, ". ", fsyear, " - ", femon, ". ", feyear, " (Issued: ", fiyearmon, ")", sep="")
        mtitle = paste("6-month Anomaly Temperature (2m) Forecast for ", fsmon, ". ", fsyear, " - ", femon, ". ", feyear, " (Issued: ", fiyearmon, ")", sep="")
        ylabel = "Anomaly temperature (C)"
      } else {
        #mtitle = paste(substr(mmetype,1,1), "-month precipitation Forecast for ", fsmon, ". ", fsyear, " - ", femon, ". ", feyear, " (Issued: ", fiyearmon, ")", sep="")
        mtitle = paste("6-month precipitation Forecast for ", fsmon, ". ", fsyear, " - ", femon, ". ", feyear, " (Issued: ", fiyearmon, ")", sep="")
        ylabel = "Precipitation (mm/month)"
      }

      boxplot(value~yearmon,data=bd, boxwex=0.3, xlab="Year-Month", ylab=ylabel, ylim = c(ymin, ymax), main= mtitle)
      lines(ldata$yearmon, ldata$CLIM, lwd=1, lty=2, col="blue")
      points(ldata$yearmon, ldata$CLIM, pch=2, col="blue")
      lines(ldata$yearmon, ldata$MME, lwd=2, col="black")
      lines(ldata$yearmon, ldata$OBS, lwd=2, lty=3, col="red")
      points(ldata$yearmon, ldata$OBS, pch=0, col="red")
      legend("top", inset=.02, bty="n", y.intersp=-0.5, horiz=T, c("Climatology", "Observed", "Forecasted MME"),
             cex=0.9, col=c("blue", "red", "black"), pch=c(2,0,-1), lwd=c(1,2,2), lty=c(2,3,1))

      # 2nd panel
      par(mar = c(2,1,0,1))
      plot(c(0:8), c(0:8), axes = F, xlab = "", ylab = "", type = "n")
      addtable2plot(0, 1, sctbl, bty="o", cex=1.0, display.colnames=T, display.rownames=T, hlines=T, vlines=T,title=paste("Monthly skil score for ", ssmon, ". ", ssyear, " - ", semon, ". ", seyear, sep=""))

      # 3rd panel
      plot(c(0:8), c(0:8), axes = F, xlab = "", ylab = "", type = "n")
      addtable2plot(0, 1, prob, bty="o", cex=1.0, display.colnames=T, display.rownames=T, hlines=T, vlines=T,title= paste("Probapbility forecast for ", feyear, sep=""))

      dev.off()

    }


  } # Variable Loop

}

#' @export
RTFcst.Combine.Forecast.Output <- function(EnvList, varnm, lagtime){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)

  syear = syear_mme

  indir = paste(prjdir, "/0_RTForecast/Monthly-TSeries", sep="")

  cnt = 1
  for(ii in 1:12){

    monnm = sprintf("M%02d", ii)

    ### Continous forecast
    srchstr = paste("RTFcst-", varnm, "-", monnm, "-Cont-", sprintf("LT%02d", lagtime), ".csv", sep="")
    InDFile = list.files(indir, pattern = glob2rx(srchstr), full.names = T)

    if(length(InDFile) > 0){

      data = read.csv(InDFile, header=T)

      if(cnt == 1){
        out = data
      } else {
        tmp = data
        out = rbind.fill(out, tmp)
      }
      cnt = cnt + 1

    }

  } # Month Loop

  #setwd("T:/")
  #write.csv(out, "test-rbindfill.csv", row.names=F)
  enddates = seq(Sys.Date(), length.out = lagtime, by="month")
  yearmon = as.data.frame(substr(seq(as.Date(paste(syear, "-01-01", sep="")), enddates[length(enddates)], by="month"), 1, 7))
  colnames(yearmon) = "yearmon"

  if(cnt > 1){
    # Change column order
    mondata = out[c("yearmon")]
    valdata = out[-which(names(out) %in% c("yearmon", "MME", "OBS", "CLIM"))]
    if(ncol(valdata) > 1){valdata = valdata[, order(substr(names(valdata), (nchar(names(valdata))-1), nchar(names(valdata))))]}
    out = cbind(mondata, valdata)

    # Overall MME by considering different models and lead times
    rcnt = nrow(out)
    out$MME = NA
    for(i in 1:rcnt){
      if(!all(is.na(out[i, 2:(ncol(out)-1)]))){
        out[i, c("MME")] = mean(as.numeric(out[i, 2:(ncol(out)-1)]), na.rm=T)
      }
    }

    out = merge(yearmon, out, all=T)

    ### Get observed data and merge to output df
    VarDFile = paste(vardir, "/", varfile, sep="")
    obs = read.csv(VarDFile, header=T)
    obs = obs[,c("yearmon", varnm)]
    colnames(obs) = c("yearmon", "OBS")
    if(varnm == "prec"){ obs = prcp.mmday2mmmon(obs) }

    out = merge(out, obs, by="yearmon", all=T)

    ### Calculate Climatology and merge to output df
    obs$month = as.numeric(substr(obs$yearmon, 6, 7))
    clim = aggregate(OBS ~ month, data = obs, FUN = mean)
    out$CLIM = NA
    for(i in 1:12){
      out[which(as.numeric(substr(out$yearmon,6,7)) == i), c("CLIM")] = clim$OBS[i]
    }

    out = out[which(as.numeric(substr(out$yearmon,1,4)) >= syear), ]
    out = out[order(as.numeric(substr(out$yearmon, 1, 4)), as.numeric(substr(out$yearmon, 6, 7))),]
    #out = out[!(is.na(out$MME)) | !(is.na(out$OBS)),]

  } else {
    out = NA
  }

  return(out)
}

#' @export
RTFcst.Adjust.Mean.SD <- function(indata, varnm){

  #setwd("E:/SForecast-TestRun/Korea/0_RTForecast/Monthly-TSeries/Original")
  #infile = "RTFcst-prec-M02-Cont-LT01.csv"
  #indata = read.csv(infile, header=T)

  yearmon = indata$yearmon
  OBS = indata$OBS
  CLIM = indata$CLIM

  omean = mean(indata$OBS, na.rm=T)
  osd = sd(indata$OBS, na.rm=T)

  nc = ncol(indata)
  for(i in 2:(nc-2)){

    simdata = indata[, i]
    colnm = colnames(indata[c(i)])

    smean = mean(simdata, na.rm=T)
    ssd = sd(simdata, na.rm=T)

    adjsd <- function(AdjCoef) { sd(omean + AdjCoef * (simdata-smean), na.rm=T) }
    diffsd <- function(AdjCoef) { abs(osd - adjsd(AdjCoef)) }
    opt = optimize(diffsd, interval = c(0, 10), maximum = F)

    adjsim = as.data.frame(omean + opt$minimum * (simdata-smean))
    colnames(adjsim) = colnm

    if(varnm == "prec") {adjsim[which(adjsim[1] < 0),] = 0}

    if(i == 2){
      outdata = cbind(yearmon, adjsim)
    } else {
      outdata = cbind(outdata, adjsim)
    }

  }

  outdata = cbind(outdata, OBS)
  outdata = cbind(outdata, CLIM)

  return(outdata)

}

#' @export
RTFcst.Combine.Forecaste.prec.t2m <- function(EnvList){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  fcstmode = EnvList$fcstmode
  combnmode = as.logical(EnvList$combnmode)
  tscale = EnvList$tscale

  if(tscale=="daily"){outdir = paste(prjdir, "/0_RTForecast/Daily_Sampling", sep="")}
  if(tscale=="hourly"){outdir = paste(prjdir, "/0_RTForecast/Hourly_Sampling", sep="")}
  SetWorkingDir(outdir)

  # There can be missing year, and merging using all=T
  ltcnt = 6
  enddate = seq(Sys.Date(), by = "month", length = ltcnt+1)[ltcnt+1]
  yearmon = data.frame(substr(seq(as.Date(paste(syear_mme, "-01-01", sep="")), enddate, by = "month") ,1,7))
  colnames(yearmon) = c("yearmon")

  OutDFile = paste(outdir, "/StnMean-All-Simulations.csv", sep="")
  rtfdir = paste(prjdir, "/0_RTForecast/Monthly-TSeries", sep="")

  syear = syear_mme
  eyear = as.numeric(substr(enddate, 1, 4))

  if(combnmode == T){
    cnt = 1
    for(j in 1:ltcnt){
      ltnm = sprintf("LT%02d",j)

      rtfdir = paste(prjdir, "/0_RTForecast/Monthly-TSeries", sep="")

      for(ii in 1:12){

        monnm = sprintf("M%02d", ii)

        ##### Precipitation######################
        PrcpDFile = paste(rtfdir, "/RTFcst-prec-", monnm, "-", fcstmode, "-",  ltnm, ".csv", sep="")
        if(file.exists(PrcpDFile)) {

          prcp = read.csv(PrcpDFile, header=T)
          prcp = prcp[which(substr(prcp$yearmon,1,4) >= syear &  substr(prcp$yearmon,1,4) <= eyear),]

          prcp = prcp.mmmon2mmday(prcp)

          rcnt = nrow(prcp)
          prcp$RType = NA
          for(i in 1:rcnt){
            if(ncol(prcp) == 5) {
              if(!is.na(prcp[i, c("MME")])) {
                prcp[i, c("RType")] = "ONE"
              } else {
                prcp[i, c("RType")] = "CLIM"
              }
            } else {
              if(is.na(prcp[i, c("MME")])) {
                prcp[i, c("RType")] = "CLIM"
              } else if (sum(!is.na(prcp[i, 2:(ncol(prcp)-4)])) == 1) {
                prcp[i, c("RType")] = "ONE"
                prcp[i, c("MME")] = NA
              } else {
                prcp[i, c("RType")] = "MMM"
              }
            }
          }

          ## if MME is NA, fill using as CLIM
          #for(i in 1:rcnt){
          #  if(is.na(prcp[i, c("MME")])){prcp[i, c("MME")] = prcp[i, c("CLIM")]}
          #}

          # Add mininum column
          prcp$MIN = NA
          for(i in 1:rcnt){
            if(sum(!is.na(prcp[i, 2:(ncol(prcp)-5)])) >= 2){
              prcp[i, c("MIN")] = min(as.numeric(prcp[i, 2:(ncol(prcp)-5)]), na.rm=T)
              #if(is.na(prcp[i, c("MIN")])){prcp[i, c("MIN")] = prcp[i, c("CLIM")]}
            }
          }

          # Add maximum column
          prcp$MAX = NA
          for(i in 1:rcnt){
            if(sum(!is.na(prcp[i, 2:(ncol(prcp)-6)])) >= 2 ){
              prcp[i, c("MAX")] = max(as.numeric(prcp[i, 2:(ncol(prcp)-6)]), na.rm=T)
              #if(is.na(prcp[i, c("MAX")])){prcp[i, c("MAX")] = prcp[i, c("CLIM")]}
            }
          }

          prcp$ONE = NA
          for(i in 1:rcnt){
            if(prcp[i, "RType"] == "ONE"){
              prcp[i, c("ONE")] = max(as.numeric(prcp[i, 2:(ncol(prcp)-7)]), na.rm=T)
            }
            if(prcp[i, "RType"] == "CLIM"){
              prcp[i, c("ONE")] = prcp[i, c("CLIM")]
            }

          }

          prcp = prcp[,c("yearmon", "RType", "MME", "MIN", "MAX", "ONE")]

        } else {

          VarDFile = paste(vardir, "/", varfile, sep="")
          prcp = read.csv(VarDFile, header=T)
          prcp = merge(yearmon, prcp, by="yearmon", all=T)
          prcp = prcp[which(substr(prcp$yearmon,1,4) >= syear & substr(prcp$yearmon,1,4) <= eyear),c("yearmon", "prec")]
          prcp$month = as.numeric(substr(prcp$yearmon, 6, 7))
          clim_prcp = aggregate(prec ~ month, data = prcp, FUN = mean, na.rm=T)
          prcp$RType = "CLIM"
          prcp$MME = NA
          prcp$MIN = NA
          prcp$MAX = NA
          for(i in 1:12){
            prcp[which(as.numeric(substr(prcp$yearmon,6,7)) == i), c("ONE")] = clim_prcp$prec[i]
          }
          prcp = prcp[which(as.numeric(substr(prcp$yearmon,6,7)) == ii),c("yearmon", "RType", "MME", "MIN", "MAX", "ONE")]

        }

        ##### Temperature ######################
        T2mDFile = paste(rtfdir, "/RTFcst-t2m-", monnm, "-", fcstmode, "-",  ltnm, ".csv", sep="")
        if(file.exists(T2mDFile)){
          t2m = read.csv(T2mDFile, header=T)

          t2m = t2m[which(substr(t2m$yearmon,1,4) >= syear &  substr(t2m$yearmon,1,4) <= eyear),]

          rcnt = nrow(t2m)
          t2m$TType = NA
          for(i in 1:rcnt){
            if(is.na(t2m[i, c("MME")])) {
              t2m[i, c("TType")] = "CLIM"
            } else {
              t2m[i, c("TType")] = "MME"
            }
          }

          # if MME is NA, fill using as CLIM
          for(i in 1:rcnt){
            if(is.na(t2m[i, c("MME")])){t2m[i, c("MME")] = t2m[i, c("CLIM")]}
          }

          t2m = t2m[,c("yearmon", "TType", "MME")]

        } else {

          ### Get observed data and merge to output df
          VarDFile = paste(vardir, "/", varfile, sep="")
          t2m = read.csv(VarDFile, header=T)
          t2m = merge(yearmon, t2m, by="yearmon", all=T)
          t2m = t2m[which(substr(t2m$yearmon,1,4) >= syear & substr(t2m$yearmon,1,4) <= eyear),c("yearmon", "t2m")]
          t2m$month = as.numeric(substr(t2m$yearmon, 6, 7))
          clim_t2m = aggregate(t2m ~ month, data = t2m, FUN = mean)
          t2m$MME = NA
          t2m$TType = "CLIM"
          for(i in 1:12){
            t2m[which(as.numeric(substr(t2m$yearmon,6,7)) == i), c("MME")] = clim_t2m$t2m[i]
          }
          t2m = t2m[which(as.numeric(substr(t2m$yearmon,6,7)) == ii),c("yearmon", "TType", "MME")]

        }

        colstrs = c("MME", "MIN", "MAX", "ONE")
        for(k in 1:length(colstrs)){

          out = prcp
          out$LTime = j
          out$ESM = fcstmode
          out$model = colstrs[k]
          outdat = out[c("yearmon", "model", "LTime", "ESM", "RType", colstrs[k])]
          outdat = merge(outdat, t2m, by="yearmon", all=T)
          colnames(outdat) = c("yearmon", "model", "LTime", "ESM","RType",  "prec", "TType", "t2m")

          if(cnt == 1){
            data = outdat
            cnt = cnt + 1
          } else {
            tmp = outdat
            data = rbind(data, tmp)
            cnt = cnt + 1
          }


        } # End k

      } # Month Loop


    } # LeadTime Loop

    data = na.omit(data)

  # combnmode == F
  } else {

    colstrs = c("MME", "MIN", "MAX")

    cnt = 1
    for(j in 1:ltcnt){
      ltnm = sprintf("LT%02d",j)

      for(ii in 1:12){

        monnm = sprintf("M%02d", ii)

        ##### Precipitation######################
        PrcpDFile = paste(rtfdir, "/RTFcst-prec-", monnm, "-", fcstmode, "-",  ltnm, ".csv", sep="")
        if(file.exists(PrcpDFile)){

          prcp = read.csv(PrcpDFile, header=T)
          prcp = prcp[which(substr(prcp$yearmon,1,4) >= syear &  substr(prcp$yearmon,1,4) <= eyear),]

          prcp = prcp.mmmon2mmday(prcp)

          rcnt = nrow(prcp)
          # if MME is NA, fill using as CLIM
          for(i in 1:rcnt){
            if(is.na(prcp[i, c("MME")])){prcp[i, c("MME")] = prcp[i, c("CLIM")]}
          }

          # Add mininum column
          prcp$MIN = NA
          for(i in 1:rcnt){
            if(!all(is.na(prcp[i, 2:(ncol(prcp)-3)]))){
              prcp[i, c("MIN")] = min(as.numeric(prcp[i, 2:(ncol(prcp)-3)]), na.rm=T)
              if(is.na(prcp[i, c("MIN")])){prcp[i, c("MIN")] = prcp[i, c("CLIM")]}
            }
          }

          # Add maximum column
          prcp$MAX = NA
          for(i in 1:rcnt){
            if(!all(is.na(prcp[i, 2:(ncol(prcp)-3)]))){
              prcp[i, c("MAX")] = max(as.numeric(prcp[i, 2:(ncol(prcp)-3)]), na.rm=T)
              if(is.na(prcp[i, c("MAX")])){prcp[i, c("MAX")] = prcp[i, c("CLIM")]}
            }
          }

          prcp = prcp[,c("yearmon", "MME", "MIN", "MAX")]

        } else {

          ### Get observed data and merge to output df
          VarDFile = paste(vardir, "/", varfile, sep="")
          prcp = read.csv(VarDFile, header=T)
          prcp = merge(yearmon, prcp, by="yearmon", all=T)
          prcp = prcp[which(substr(prcp$yearmon,1,4) >= syear & substr(prcp$yearmon,1,4) <= eyear),c("yearmon", "prec")]
          prcp$month = as.numeric(substr(prcp$yearmon, 6, 7))
          clim_prcp = aggregate(prec ~ month, data = prcp, FUN = mean, na.rm=T)
          prcp$MME = NA
          for(i in 1:12){
            prcp[which(as.numeric(substr(prcp$yearmon,6,7)) == i), c("MME")] = clim_prcp$prec[i]
          }
          prcp$MIN = prcp$MME; prcp$MAX = prcp$MME
          prcp = prcp[which(as.numeric(substr(prcp$yearmon,6,7)) == ii),c("yearmon", "MME", "MIN", "MAX")]

        }

        ##### Temperature ######################
        T2mDFile = paste(rtfdir, "/RTFcst-t2m-", monnm, "-", fcstmode, "-",  ltnm, ".csv", sep="")
        if(file.exists(T2mDFile)){
          t2m = read.csv(T2mDFile, header=T)

          t2m = t2m[which(substr(t2m$yearmon,1,4) >= syear &  substr(t2m$yearmon,1,4) <= eyear),]

          rcnt = nrow(t2m)

          # if MME is NA, fill using as CLIM
          for(i in 1:rcnt){
            if(is.na(t2m[i, c("MME")])){t2m[i, c("MME")] = t2m[i, c("CLIM")]}
          }

          t2m = t2m[,c("yearmon", "MME")]

        } else {

          ### Get observed data and merge to output df
          VarDFile = paste(vardir, "/", varfile, sep="")
          t2m = read.csv(VarDFile, header=T)
          t2m = merge(yearmon, t2m, by="yearmon", all=T)
          t2m = t2m[which(substr(t2m$yearmon,1,4) >= syear & substr(t2m$yearmon,1,4) <= eyear),c("yearmon", "t2m")]
          t2m$month = as.numeric(substr(t2m$yearmon, 6, 7))
          clim_t2m = aggregate(t2m ~ month, data = t2m, FUN = mean)
          t2m$MME = NA
          for(i in 1:12){
            t2m[which(as.numeric(substr(t2m$yearmon,6,7)) == i), c("MME")] = clim_t2m$t2m[i]
          }
          t2m = t2m[which(as.numeric(substr(t2m$yearmon,6,7)) == ii),c("yearmon", "MME")]

        }

        for(k in 1:3){
          out = prcp
          out$LTime = j
          out$ESM = fcstmode
          out$model = colstrs[k]
          outdat = out[c("yearmon", "model", "LTime", "ESM", colstrs[k])]
          outdat = merge(outdat, t2m, by="yearmon", all=T)
          colnames(outdat) = c("yearmon", "model", "LTime", "ESM", "prec", "t2m")

          if(cnt == 1){
            data = outdat
            cnt = cnt + 1
          } else {
            tmp = outdat
            data = rbind(data, tmp)
            cnt = cnt + 1
          }


        } # End k

      } # Month Loop

      ### Check Below is necessary
      # VarDFile = paste(vardir, "/", varfile, sep="")
      # vardata = read.csv(VarDFile, header=T)
      # prcp = vardata[which(substr(vardata$yearmon,1,4) >= syear & substr(vardata$yearmon,1,4) <= eyear),c("yearmon", "prec")]
      # prcp$month = as.numeric(substr(prcp$yearmon, 6, 7))
      # clim_prcp = aggregate(prec ~ month, data = prcp, FUN = mean)
      #
      # t2m = vardata[which(substr(vardata$yearmon,1,4) >= syear & substr(vardata$yearmon,1,4) <= eyear),c("yearmon", "t2m")]
      # t2m$month = as.numeric(substr(t2m$yearmon, 6, 7))
      # clim_t2m = aggregate(t2m ~ month, data = t2m, FUN = mean)
      #
      # dcnt = nrow(data)
      # for(i in 1:dcnt){
      #   if(is.na(data[i, c("prec")])){data[i, c("prec")] = clim_prcp$prec[as.numeric(substr(data$yearmon[i],6,7))]}
      #   if(is.na(data[i, c("t2m")])){data[i, c("t2m")] = clim_t2m$t2m[as.numeric(substr(data$yearmon[i],6,7))]}
      # }

    } # LeadTime Loop


  }

  write.csv(data, OutDFile, row.names=F)
  cat(sprintf("     RTFcst: Combining of prec and t2m has been compleated!\n"))


}

#' @export
RTFcst.Mahalanobis.Sampling.Observed.9var <- function(EnvList, fiyearmon, smplmoncnt){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  smpldir = EnvList$smpldir
  obsdir = EnvList$obsdir
  syear_obs = as.numeric(EnvList$syear_obs)
  eyear_obs = as.numeric(EnvList$eyear_obs)

  tscale = EnvList$tscale
  fcstmode = EnvList$fcstmode
  combnmode = as.logical(EnvList$combnmode)
  fiyearmode = as.logical(EnvList$fiyearmode)

  mdlnms = c("MIN", "MME", "MAX")

  options(warn=-1)
  options(stringsAsFactors = FALSE)
  ltcnt = 6

  if(tscale=="daily"){outdir = paste(prjdir, "/0_RTForecast/Daily_Sampling", sep="")}
  if(tscale=="hourly"){outdir = paste(prjdir, "/0_RTForecast/Hourly_Sampling", sep="")}

  scndir = paste(outdir, "/", fiyearmon, sep="")
  # reset files
  unlink(scndir, recursive=TRUE)
  SetWorkingDir(scndir)

  ###### Get Station ID, lat, and Lon information
  #setwd(stndir)
  StnDFile = paste(stndir, stnfile, sep="/")
  stninfo = read.csv(StnDFile, header=T)
  stninfo = stninfo[which(stninfo$SYear <= syear_obs),c("ID", "Lon", "Lat", "Elev", "SYear")]
  colnames(stninfo) = c("ID", "Lon", "Lat", "Elev", "SYear")
  stnnms = matrix(stninfo$ID)
  stncnt = length(stnnms)

  ### save station means based on KMA observed data(monthly)
  #if(tscale=="daily"){ obsdir = paste(smpldir, "/asos-daily", sep="") }
  #if(tscale=="hourly"){ obsdir = paste(smpldir, "/asos-hourly/daily_summary", sep="") }
  #obsdir = smpldir

  covdata = Calculate.Covariance.Matrix(obsdir, outdir, stnnms, syear_obs, eyear_obs)

  RTFcst.Combine.Forecaste.prec.t2m(EnvList)

  ### Save best month information
  Select.BestFit.Month(outdir, "StnMean-All-Simulations.csv")

  ################# Select bet-fit month based on given scenario
  RTFcst.Select.Forecast.Ensemble (EnvList, fiyearmon, smplmoncnt)

  #if(tscale=="daily"){scndir = paste(prjdir, "/0_RTForecast/Daily_Sampling", "/", fiyearmon, sep="")}
  #if(tscale=="hourly"){scndir = paste(prjdir, "/0_RTForecast/Hourly_Sampling", "/", fiyearmon, sep="")}

  ################# Daily or Hourly data sampling ################################################
  srchstr = paste("BestMon-*.csv", sep="")
  flist = list.files(scndir, pattern = glob2rx(srchstr), full.names = F)
  mdlnms = unique(matrix(unlist(strsplit(flist[], "-")), nrow=length(flist), byrow=T)[,2])
  mdlnms = substr(mdlnms, 1, nchar(mdlnms)-4)
  mdlcnt = length(mdlnms)

  for(i in 1:mdlcnt){
    mdlnm = mdlnms[i]

    if(tscale == "daily"){
      #dsmpldir = smpldir
      dsmpldir = paste(smpldir, "/asos-daily", sep="")
      outdata = Sampling.Daily.Data.9var(scndir, dsmpldir, stninfo, mdlnm)
    }

    if(tscale == "hourly"){
      #hsmpldir = smpldir
      hsmpldir = paste(smpldir, "/asos-hourly", sep="")
      outdata = Sampling.Hourly.Data.9var(scndir, hsmpldir, stninfo, mdlnm)
    }
  }

}

#' @export
RTFcst.Select.Forecast.Ensemble <- function(EnvList, fiyearmon, smplmoncnt){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  stndir = EnvList$stndir
  stnfile = EnvList$stnfile
  smpldir = EnvList$smpldir
  syear_obs = as.numeric(EnvList$syear_obs)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  fiyearmode = as.logical(EnvList$fiyearmode)
  tscale = EnvList$tscale
  combnmode = as.logical(EnvList$combnmode)
  fiyearmode = as.logical(EnvList$fiyearmode)


  if(tscale=="daily"){outdir = paste(prjdir, "/0_RTForecast/Daily_Sampling", sep="")}
  if(tscale=="hourly"){outdir = paste(prjdir, "/0_RTForecast/Hourly_Sampling", sep="")}

  scndir = paste(outdir, "/", fiyearmon, sep="")
  SetWorkingDir(scndir)
  ltcnt = 6

  # Function to get model name
  Get.Model.Name <- function(ltdata, rowno){
    if(ltdata$model[rowno] == "MAX"){ mdnm = "X"}
    if(ltdata$model[rowno] == "MIN"){ mdnm = "N"}
    if(ltdata$model[rowno] == "MME"){ mdnm = "M"}
    if(ltdata$RType[rowno] == "ONE"){ mdnm = "O"}
    if(ltdata$RType[rowno] == "CLIM"){ mdnm = "C"}
    return(mdnm)
  }

  # Function to combine row data
  Row.Combine.Data <- function(lt0d, lt1d, lt2d, lt3d, lt4d, lt5d, lt6d, fiyearmode){
    if(!any(is.na(lt0d)) & fiyearmode == T){ mondata = lt0d}
    if(!any(is.na(lt1d))){
      if(fiyearmode == T) {
        mondata = rbind(mondata, lt1d)
      } else {
        mondata = lt1d
      }
    }
    if(!any(is.na(lt2d))){ mondata = rbind(mondata, lt2d)}
    if(!any(is.na(lt3d))){ mondata = rbind(mondata, lt3d)}
    if(!any(is.na(lt4d))){ mondata = rbind(mondata, lt4d)}
    if(!any(is.na(lt5d))){ mondata = rbind(mondata, lt5d)}
    if(!any(is.na(lt6d))){ mondata = rbind(mondata, lt6d)}

    mondata = mondata[order(mondata$yearmon),]

    return(mondata)
  }

  data = read.csv(paste(outdir, "/BestMon-All-Simulations.csv", sep=""), header=T)

  if(combnmode == T){

    fyearmons = substr(seq(as.Date(paste(fiyearmon, "-01", sep="")), by = "month", length = (ltcnt+1)),1,7)
    lt1 = data[which(data$yearmon == fyearmons[2] & data$LTime == 1), ]; lt1cnt = nrow(lt1)
    lt2 = data[which(data$yearmon == fyearmons[3] & data$LTime == 2), ]; lt2cnt = nrow(lt2)
    lt3 = data[which(data$yearmon == fyearmons[4] & data$LTime == 3), ]; lt3cnt = nrow(lt3)
    lt4 = data[which(data$yearmon == fyearmons[5] & data$LTime == 4), ]; lt4cnt = nrow(lt4)
    lt5 = data[which(data$yearmon == fyearmons[6] & data$LTime == 5), ]; lt5cnt = nrow(lt5)
    lt6 = data[which(data$yearmon == fyearmons[7] & data$LTime == 6), ]; lt6cnt = nrow(lt6)

    lt0d = data[which(data$yearmon == fyearmons[1] & data$LTime == 1 & (data$model == "MME" | data$model == "ONE")), ]
    for(i1 in 1:lt1cnt){
      if(smplmoncnt >= 1){
        md1 = Get.Model.Name (lt1, i1)
        lt1d = lt1[i1, ]
      } else {
        md1 = ""; lt1d = NA
      }

      for(i2 in 1:lt2cnt){
        if(smplmoncnt >= 2){
          md2 = Get.Model.Name (lt2, i2)
          lt2d = lt2[i2, ]
        } else {
          md2 = ""; lt2d = NA
        }

        for(i3 in 1:lt3cnt){
          if(smplmoncnt >= 3){
            md3 = Get.Model.Name (lt3, i3)
            lt3d = lt3[i3, ]
          } else{
            md3 = ""; lt3d = NA
          }

          for(i4 in 1:lt4cnt){
            if(smplmoncnt >= 4){
              md4 = Get.Model.Name (lt4, i4)
              lt4d = lt4[i4, ]
            } else{
              md4 = ""; lt4d = NA
            }

            for(i5 in 1:lt5cnt){
              if(smplmoncnt >= 5){
                md5 = Get.Model.Name (lt5, i5)
                lt5d = lt5[i5, ]
              } else{
                md5 = ""; lt5d = NA
              }

              for(i6 in 1:lt6cnt){
                if(smplmoncnt >= 6){
                  md6 = Get.Model.Name (lt6, i6)
                  lt6d = lt6[i6, ]
                } else{
                  md6 = ""; lt6d = NA
                }

                mondata = Row.Combine.Data (lt0d, lt1d, lt2d, lt3d, lt4d, lt5d, lt6d, fiyearmode)
                mdlnm = paste(md1, md2, md3, md4, md5, md6, sep ="")

                OutDFile = paste(scndir, "/BestMon-", mdlnm, ".csv",  sep="")
                write.csv(mondata, OutDFile, row.names=F)
                cat(sprintf("     RTFcst: Best-fit month was selected based on scenario: Model=%s\n" ,mdlnm))

              }
            }
          }
        }
      }
    }

  } else {

    mdlnms = c("MME", "MAX", "MIN")

    mdlcnt = length(mdlnms)
    for(i in 1:mdlcnt){
      mdlnm = mdlnms[i]

      # fiyearmon: year-month issuing the forecast
      fyearmons = substr(seq(as.Date(paste(fiyearmon, "-01", sep="")), by = "month", length = (smplmoncnt+1)),1,7)

      if(fiyearmode == T){
        startno = 1
      } else {
        startno = 2
      }

      for(k in startno:(smplmoncnt+1)){

        cyearmon = fyearmons[k]

        if(k == startno){
          if(k==1) {
            mondata = data[which(data$model == mdlnm & data$yearmon == cyearmon & data$LTime == (k)), ]
          } else {
            mondata = data[which(data$model == mdlnm & data$yearmon == cyearmon & data$LTime == (k-1)), ]
          }
          if(nrow(mondata) == 0){
            mondata = as.data.frame(cbind(cyearmon, mdlnm, k, "SME", NA, NA, NA, NA, NA))
            colnames(mondata) = c("yearmon", "model", "LTime", "ESM", "prec", "t2m", "bestmon", "prec_obs", "t2m_obs")
          }
        } else {
          imsi = data[which(data$model == mdlnm & data$yearmon == cyearmon & data$LTime == (k-1)), ]
          if(nrow(imsi) == 0){
            imsi = cbind(cyearmon, mdlnm, (k-1), "SME", NA, NA, NA, NA, NA)
            colnames(imsi) = c("yearmon", "model", "LTime", "ESM", "prec", "t2m", "bestmon", "prec_obs", "t2m_obs")
          }
          mondata = rbind(mondata, imsi)
        }
      }
      mondata = mondata[order(mondata$yearmon),]

      OutDFile = paste(scndir, "/BestMon-", mdlnm, ".csv",  sep="")
      write.csv(mondata, OutDFile, row.names=F)
      cat(sprintf("     RTFcst: Best-fit month was selected based on scenario: Model=%s\n" ,mdlnm))

    }

  }

}

#' @export
MWRObs.Copy.Hindcast.Results <- function(EnvList, tyear){

  prjdir = EnvList$prjdir
  eyear_obs = as.numeric(EnvList$eyear_obs)

  if(tyear > eyear_obs){

    srcdir = paste(prjdir, "/MWRObs/", eyear_obs, "/01_Regression-all", sep="")
    dstdir = paste(prjdir, "/MWRObs/", tyear, sep="")
    SetWorkingDir(dstdir)

    file.copy(srcdir, dstdir, recursive = T)

  }

}

#' @export
MWRObs.Extract.Predictor.Grid <- function(rnldir, ptrfile, cmonth, LTime, syear, eyear){

  loncnt = 144; latcnt = 73

  ptrinfo = paste(rnldir, "/", ptrfile, sep="")
  ptrlist = read.csv(ptrinfo, header=T)
  ptrcnt = nrow(ptrlist)

  yrcnt = eyear - syear + 1

  # define empty array
  rnl = array(NA,dim=list(loncnt,latcnt,yrcnt,ptrcnt))

  for(i in 1:ptrcnt) {

    ptrnm = ptrlist[i, "ptrnm"]; ncfile = ptrlist[i, "fname"]; ptrlvl = ptrlist[i, "level"]

    #rnldir = "F:/SForecast-ADSS/Database/reanalysis1"
    #ncfile = "vwnd.mon.mean.nc"
    out = ReadReanalysis1(rnldir, ncfile)
    lvls = out$l; times = out$t; var = out$var

    # Select given level
    if(length(lvls) > 1){
      loclvl = which(lvls == ptrlvl)
      var = var[ , , loclvl, ]
    }

    # Consider LagTime: if LTime = 3month, predictor for current(eg. 2016-04) comes from 2016-01
    nmon = dim(var)[3]
    var[ , , (1+LTime):nmon] = var[ , , 1:(nmon-LTime)]
    var[ , , 1:LTime] = NA

    # Select current month
    imsi = var[ , , which(as.numeric(substr(times, 1, 4)) >= syear &  as.numeric(substr(times, 1, 4)) <= eyear & as.numeric(substr(times, 6, 7)) == cmonth)]
    if(dim(imsi) ==  3) {
      yrimsi = dim(imsi)[3]
      if(yrimsi == yrcnt){
        rnl[ , , , i] = imsi
      } else {
        rnl[ , , 1:yrimsi, i] = imsi
      }
    }

    if(dim(imsi) ==  2) {
      rnl[ , , 1, i] = imsi
    }

  }

  return(rnl)

}

#' @export
MWRObs.Calculate.All.Regression <- function(EnvList, tyear){

  prjdir = EnvList$prjdir
  rnldir = EnvList$rnldir
  ptrfile = EnvList$ptrfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_mme = as.numeric(EnvList$eyear_mme)
  eyear_sim = as.numeric(EnvList$eyear_sim)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)
  minlat = as.numeric(EnvList$minlat)
  maxlat = as.numeric(EnvList$maxlat)
  MinGrdCnt = as.numeric(EnvList$MinGrdCnt)
  MaxGrdCnt = as.numeric(EnvList$MaxGrdCnt)

  lagcnt = 8

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  ptrinfo = paste(rnldir, "/", ptrfile, sep="")
  ptrlist = read.csv(ptrinfo, header=T)
  ptrnms = ptrlist$ptrnm

  minrow = round((minlat - (-90))/2.5 + 1, digits=0)
  maxrow = round(72- (90 - maxlat)/2.5 + 1, digits=0)

  # define and create output folder
  outdir = paste(prjdir, "/MWRObs/", tyear, "/01_Regression-all", sep="")
  SetWorkingDir(outdir)

  ### Get variable(predictand) data
  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms_mwr = names(var)[2:ncol(var)]

  ######## Regression Analysis for all combinations on Cross Validation mode
  varcnt = length(varnms_mwr)
  for(i in 1:varcnt){
    varnm = varnms_mwr[i]

    for(j in smonth:emonth){

      cmonth = j

      # Check the eyear_mme by comparing to the max. year of observed data(myear_obs)
      myear_obs = max(as.numeric(substr(var[which(as.numeric(substr(var$yearmon,6,7)) == j), c("yearmon")],1,4)))
      if(eyear_mme > myear_obs){
        eyear_mme = myear_obs
        cat(sprintf("     MWRObs: End of training period exceeded the end of observed data: Var=%s,  Mon=%s\n", varnm, cmonth))
      }

      # Check eyear_sim by comparing to current yearand month
      curyear = as.numeric(substr(Sys.Date(), 1, 4))
      curmon = as.numeric(substr(Sys.Date(), 6, 7))
      if(cmonth >= curmon & eyear_sim >= curyear){
        eyear_sim = eyear_sim - 1
        cat(sprintf("     MWRObs: End year of simulation has been changed!\n"))
      }

      varmon = var[which(as.numeric(substr(var$yearmon,6,7)) == j & as.numeric(substr(var$yearmon,1,4)) >= syear_mme & as.numeric(substr(var$yearmon,1,4)) <= eyear_mme),]
      # Exclude target year
      varmon = varmon[which(!(as.numeric(substr(varmon$yearmon,1,4)) == tyear)),]

      for(ii in 1:lagcnt){
        LTime = ii

        outfile = paste(outdir, "/MWRObs-", varnm, "-M", sprintf("%02d", cmonth), "-LT", sprintf("%02d", LTime), ".csv", sep="")

        if(!file.exists(outfile)){

          cnt = 1

          # Extract Single Model Ensemble (SME) data mme[nlon, nlat, nyear, npredictors]
          mme = MWRObs.Extract.Predictor.Grid(rnldir, ptrfile, cmonth, LTime, syear_mme, eyear_sim)
          # Exclude target year and consider min & max lat
          mme = mme[ , minrow:maxrow, , ]

          skipyr = tyear - syear_mme + 1
          cvmme = mme[ , ,-skipyr , ]


          ##############################################################
          # 2/3 of data (First 2/3, ending 2/3)
          cend = round(nrow(varmon)*2/3, digits=0)
          vend = nrow(varmon)
          vstart = vend - cend + 1

          cvardata = varmon[1:cend, c(varnm)]
          vvardata = varmon[vstart:vend, c(varnm)]

          # define number of random smapling sets
          nrsmpls = (vend-cend) * 3

          for(tt in 1:nrsmpls){
            yearnm = paste("yearnos", tt, sep="")
            vardatanm = paste("rvardata", tt, sep="")
            assign(yearnm, sample(1:vend, cend, replace=F))
            assign(vardatanm, varmon[get(yearnm), c(varnm)])
          }

          ptrcnt = length(ptrnms)
          for(m in 1:ptrcnt){

            ntsmpls = nrsmpls + 2
            nlon = length(mme[,1,1,1])
            nlat = length(mme[1,,1,1])
            cormme = array(NA, dim=c(nlon, nlat, ntsmpls))

            ptrnm = ptrnms[m]

            cmmedata = cvmme[,,1:cend, m]
            vmmedata = cvmme[,,vstart:vend, m]

            cormme[,,1] = apply(cmmedata[,,],1:2,function(x) {cor(cvardata,x, use = "na.or.complete")})
            cormme[,,2] = apply(vmmedata[,,],1:2,function(x) {cor(vvardata,x, use = "na.or.complete")})


            for(tt in 1:nrsmpls){
              mmedatanm = paste("rmmedata", tt, sep="")
              yearnm = paste("yearnos", tt, sep="")
              vardatanm = paste("rvardata", tt, sep="")
              assign(mmedatanm, cvmme[,,get(yearnm), m])
              cormme[,,tt+2] = apply(get(mmedatanm),1:2,function(x) {cor(get(vardatanm),x, use = "na.or.complete")})
            }

            # Critical TCC based on data period (syear_mme ~ eyear_mme)
            CR = critical.r(vend)

            #### Check positive correlation
            #MaxCR = 1
            MinCR = CR
            grdcnt = 0
            repeat{

              corplus = cormme
              corplus[which(corplus < MinCR)] = 0

              for(tt in 1:ntsmpls){
                if(tt == 1){
                  cbndplus = corplus[,,1]
                } else {
                  cbndplus = cbndplus * corplus[,,tt]
                }
              }

              rowcnt = nrow(cbndplus); colcnt = ncol(cbndplus)

              coor = which(cbndplus > 0, arr.ind=T)
              coorcnt = nrow(coor)

              if(coorcnt > 0){
                # Nearest neighbor distance and comparison
                nb = spdep::dnearneigh(coor, 0, 2)
                comp = n.comp.nb(nb)

                cluster = as.data.frame(cbind(coor, comp$comp.id))
                colnames(cluster) = c("row", "col", "rgns")
                cnttbl = as.data.frame(table(cluster$rgns))
                cnttbl = cnttbl[order(cnttbl$Freq, decreasing=T), ]

                #ngrids = cnttbl[1, "Freq"]
                rgnnm = cnttbl[1, "Var1"]
                coornm = cluster[which(cluster$rgns == rgnnm), c("row", "col")]
                grdcnt = nrow(coornm)
                bndnm = array(0, dim=c(rowcnt, colcnt))
                for(mm in 1:grdcnt){
                  rown = coornm[mm, "row"]
                  coln = coornm[mm, "col"]
                  bndnm[rown, coln] = 100
                }
              }

              #if(cntplus > MinGrdCnt | MaxCR <= CR){
              if(grdcnt < MaxGrdCnt | MinCR > 1 | coorcnt == 0){

                break
              }
              #MaxCR = MaxCR - 0.01
              MinCR = MinCR + 0.01
            } # End repeat

            #if(grdcnt > MinGrdCnt){
            if(grdcnt < MaxGrdCnt & grdcnt > MinGrdCnt){
              if(cnt == 1){
                pptrout = cbind(varnm, cmonth, LTime, ptrnm, "Positive", grdcnt, MinCR)
                colnames(pptrout) = c("variable", "month", "LTime", "predictor", "Slope", "Cells", "CC")
                cnt = cnt + 1
              } else {
                pptrtmp = cbind(varnm, cmonth, LTime, ptrnm, "Positive", grdcnt, MinCR)
                colnames(pptrtmp) = c("variable", "month", "LTime", "predictor", "Slope", "Cells", "CC")
                pptrout = rbind(pptrout, pptrtmp)
                cnt = cnt + 1
              }

              #cbndplus[which(cbndplus > 0)] = 100
              BndDFile = paste(outdir, "/MWRObs-", varnm, "-M", sprintf("%02d", cmonth), "-LT", sprintf("%02d", LTime), "-", ptrnm, "-Positive.csv", sep="")
              write.csv(bndnm, BndDFile, row.names=F)
            }

            #### Check negative correlation
            #MaxCR = 1
            MinCR = CR
            grdcnt = 0
            repeat{
              corminus = cormme
              corminus[which(corminus > -MinCR)] = 0

              for(tt in 1:ntsmpls){
                if(tt == 1){
                  cbndminus = corminus[,,1]
                } else {
                  cbndminus = cbndminus * corminus[,,tt]
                }
              }

              rowcnt = nrow(cbndminus); colcnt = ncol(cbndminus)

              #if(ntsmpls %% 2 == 0){coor = which(cbndminus > 0, arr.ind=T)}
              #if(ntsmpls %% 2 != 0){coor = which(cbndminus < 0, arr.ind=T)}
              coor = which(!(cbndminus == 0), arr.ind=T)
              coorcnt = nrow(coor)

              if(coorcnt > 0){
                # Nearest neighbor distance and comparison
                nb = dnearneigh(coor, 0, 2)
                comp = n.comp.nb(nb)

                cluster = as.data.frame(cbind(coor, comp$comp.id))
                colnames(cluster) = c("row", "col", "rgns")
                cnttbl = as.data.frame(table(cluster$rgns))
                cnttbl = cnttbl[order(cnttbl$Freq, decreasing=T), ]

                #ngrids = cnttbl[1, "Freq"]
                rgnnm = cnttbl[1, "Var1"]
                coornm = cluster[which(cluster$rgns == rgnnm), c("row", "col")]
                grdcnt = nrow(coornm)
                bndnm = array(0, dim=c(rowcnt, colcnt))
                for(mm in 1:grdcnt){
                  rown = coornm[mm, "row"]
                  coln = coornm[mm, "col"]
                  bndnm[rown, coln] = 100
                }
              }

              #cntminus = length(which(cbndminus != 0))
              #if(cntminus > MinGrdCnt | MaxCR <= CR){
              if(grdcnt < MaxGrdCnt | MinCR > 1 | coorcnt == 0){
                break
              }
              MinCR = MinCR + 0.01
            }

            #if(cntminus >= MinGrdCnt){
            if(grdcnt < MaxGrdCnt & grdcnt > MinGrdCnt){
              if(cnt == 1){
                pptrout = cbind(varnm, cmonth, LTime, ptrnm, "Negative", grdcnt, MinCR)
                colnames(pptrout) = c("variable", "month", "LTime", "predictor", "Slope", "Cells", "CC")
                cnt = cnt + 1
              } else {
                pptrtmp = cbind(varnm, cmonth, LTime, ptrnm, "Negative", grdcnt, MinCR)
                colnames(pptrtmp) = c("variable", "month", "LTime", "predictor", "Slope", "Cells", "CC")
                pptrout = rbind(pptrout, pptrtmp)
                cnt = cnt + 1
              }

              cbndminus[which(cbndminus != 0)] = 100
              BndDFile = paste(outdir, "/MWRObs-", varnm, "-M", sprintf("%02d", cmonth), "-LT", sprintf("%02d", LTime), "-", ptrnm, "-Negative.csv", sep="")
              #write.csv(cbndminus, BndDFile, row.names=F)
              write.csv(bndnm, BndDFile, row.names=F)
            }

          } # Loop for predictors

          ############# Time series
          if(cnt > 1){

            pptrout = as.data.frame(pptrout)
            pptrout = pptrout[order(pptrout$CC, decreasing = T),]
            write.csv(pptrout, outfile, row.names=F)

            pptrcnt = nrow(pptrout)
            for(jj in 1:pptrcnt){

              # read boundary file
              bptrnm = pptrout$predictor[jj]
              bptrno = which(ptrnms == bptrnm)
              bptrtype = pptrout$Slope[jj]

              BndDFile = paste(outdir, "/MWRObs-", varnm, "-M", sprintf("%02d", cmonth), "-", sprintf("LT%02d",LTime), "-", bptrnm, "-", bptrtype, ".csv", sep="")
              bbnd = as.data.frame(read.csv(BndDFile, header=T)/100)

              yearmon = as.data.frame(paste(seq(syear_mme, eyear_sim), "-", sprintf("%02d", cmonth), sep=""))
              colnames(yearmon) = "yearmon"
              nyr = length(mme[1,1,,1])
              ptrdata = array(NA,dim=list(nyr))
              for(tt in 1:nyr){
                bmme = mme[,,tt,bptrno] * bbnd
                bmme[bmme == 0] = NA
                ptrdata[tt] = mean(unlist(bmme), na.rm=T)
              }
              ptrout = cbind(yearmon, ptrdata)
              ptrout$varnm = varnm
              #ptrout$mdlnm = mdlnm
              ptrout$LTime = LTime
              ptrout$ptrnm = bptrnm
              ptrout$cmonth = cmonth
              colnames(ptrout) = c("yearmon", "APCC", "Variable", "LTime", "predictor", "month")
              OutDFile = paste(outdir, "/MWRObs-", varnm, "-", sprintf("M%02d",cmonth), "-", sprintf("LT%02d",LTime), "-", bptrnm, "-", bptrtype, "-TSeries.csv",  sep="")
              write.csv(ptrout, OutDFile, row.names=F)

            }

          } # IF (cnt > 1)

        }  # IF(file.exist)

        cat(sprintf("     MWRObs: Analysis has been finished Var=%s Month=%s, LTime=%s\n", varnm, cmonth, LTime))

      } # LTime loop(ii)

    } # Month loop(j)
  } # Variable loop(i)

}

#' @export
MWRObs.Extract.Best.Predictor.TSeries <- function(EnvList, tyear){

  #prjdir, rnldir, ptrfile, vardir, varfile, syear_mme, eyear_mme, eyear_obs, eyear_sim, tyear, minlat, maxlat, precopt
  prjdir = EnvList$prjdir
  rnldir = EnvList$rnldir
  ptrfile = EnvList$ptrfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_mme = as.numeric(EnvList$eyear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  eyear_sim = as.numeric(EnvList$eyear_sim)
  minlat = as.numeric(EnvList$minlat)
  maxlat = as.numeric(EnvList$maxlat)
  precopt = as.logical(EnvList$precopt)

  lagcnt = 8

  minrow = round((minlat - (-90))/2.5 + 1, digits=0)
  maxrow = round(72- (90 - maxlat)/2.5 + 1, digits=0)

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  indir = paste(prjdir, "/MWRObs/", tyear, "/01_Regression-all", sep="")
  outdir = paste(prjdir, "/MWRObs/", tyear, "/02_BPredictor_TSeries", sep="")

  SetWorkingDir(outdir)

  ### Get variable(predictand) data
  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  srchstr = "MWRObs-*TSeries.csv"
  flist = list.files(indir, pattern = glob2rx(srchstr), full.names = F)
  fcnt = length(flist)

  ChkDFile1 = paste(outdir, "/selected_with_same_sign.csv", sep="")
  if(!file.exists(ChkDFile1) & fcnt > 0){

    bcnt = 0
    for(i in 1:fcnt){

      fname = flist[i]
      varnm = unlist(strsplit(fname, "-"))[2]
      curmon = unlist(strsplit(fname, "-"))[3]
      #mdlnm = unlist(strsplit(fname, "-"))[4]
      ltime = unlist(strsplit(fname, "-"))[4]
      ptrnm = unlist(strsplit(fname, "-"))[5]
      dirnm = unlist(strsplit(fname, "-"))[6]

      PtrDFile = paste(indir, "/", fname, sep="")
      ptrdata = read.csv(PtrDFile, header=T)
      ptrdata = ptrdata[c("yearmon", "APCC")]
      colnames(ptrdata) = c("yearmon", "x")

      obsdata = var[c("yearmon", varnm)]
      colnames(obsdata) = c("yearmon", "y")

      regdata = merge(ptrdata, obsdata, by="yearmon")

      data = regdata[which(as.numeric(substr(regdata$yearmon,1,4)) >= syear_mme & as.numeric(substr(regdata$yearmon,1,4)) <= eyear_mme), ]

      # Training period (TCC between predictor and Observed)
      cnt = 1
      for(jj in syear_mme:eyear_mme){
        curx = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c("x")]
        # if curx is nothing
        if(length(curx) == 0){curx = NA}

        cyearmon = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c("yearmon")]
        if(!is.na(curx)){
          # Exclude current year
          cvaldata = data[which(!as.numeric(substr(data$yearmon,1,4)) == jj), c("yearmon", "y", "x")]
          # Exclude target year
          cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]

          fit = lm(y ~ x, data=cvaldata)

          summary = summary(fit)
          coeff = coefficients(fit)
          R2 = summary$r.squared
          sigma = summary$sigma
          pVal = anova(fit)$'Pr(>F)'[1]
          cury = coeff[1] + coeff[2]*curx
          if((precopt==T | varnm == "prec") & (!is.na(cury)) & (cury <= 0)) { cury = 0}
        }

        if(cnt == 1){
          sumout = cbind(cyearmon, t(coeff), R2, sigma, pVal, cury)
          cnt = cnt + 1
        } else {
          sumtmp = cbind(cyearmon, t(coeff), R2, sigma, pVal, cury)
          sumout = rbind(sumout, sumtmp)
          cnt = cnt + 1
        }
      }

      # Verification period (TCC between predictor and Observed)
      cvaldata = data
      # Exclude target year
      cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]

      colnames(cvaldata) = c("yearmon", "y", "x")
      fit = lm(y ~ x, data=cvaldata)
      summary = summary(fit)
      coeff = coefficients(fit)
      R2 = summary$r.squared
      sigma = summary$sigma
      pVal = anova(fit)$'Pr(>F)'[1]

      data = regdata[which(as.numeric(substr(regdata$yearmon,1,4)) >= (eyear_mme + 1) & as.numeric(substr(regdata$yearmon,1,4)) <= eyear_obs), ]
      for(jj in (eyear_mme + 1):eyear_obs){
        curx = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c("x")]
        # if curx is nothing
        if(length(curx) == 0){curx = NA}

        cyearmon = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c("yearmon")]
        if(!is.na(curx)){
          #coeff = coefficients(fit)
          #R2 = summary$r.squared
          #sigma = summary$sigma
          #pVal = anova(fit)$'Pr(>F)'[1]
          cury = coeff[1] + coeff[2]*curx
          if((precopt==T | varnm == "prec") & (!is.na(cury)) & (cury <= 0)) { cury = 0}
        }

        if(cnt == 1){
          sumout = cbind(cyearmon, t(coeff), R2, sigma, pVal, cury)
          cnt = cnt + 1
        } else {
          sumtmp = cbind(cyearmon, t(coeff), R2, sigma, pVal, cury)
          sumout = rbind(sumout, sumtmp)
          cnt = cnt + 1
        }
      } # End of Year Loop for Cross Validation

      ####### Compared TCC between forecasted and observed precipitation or temperature
      simdata = as.data.frame(sumout[,c("cyearmon", "cury")])
      colnames(simdata) = c("yearmon", "sim")

      obsdata = var[c("yearmon", varnm)]
      colnames(obsdata) = c("yearmon", "obs")

      outdata = merge(simdata, obsdata, by="yearmon")
      # Exclude target year
      outdata = outdata[which(!(as.numeric(substr(outdata$yearmon,1,4)) == tyear)), ]

      # Training and Verification period
      trndata = outdata[which(as.numeric(substr(outdata$yearmon,1,4)) >= syear_mme & as.numeric(substr(outdata$yearmon,1,4)) <= eyear_mme), ]
      veridata = outdata[which(as.numeric(substr(outdata$yearmon,1,4)) >= (eyear_mme + 1) & as.numeric(substr(outdata$yearmon,1,4)) <= eyear_sim), ]

      TrnCR = critical.r(nrow(trndata))
      VeriCR = critical.r(nrow(veridata))

      TrnTcc = cor(as.numeric(trndata$sim), as.numeric(trndata$obs), method="pearson")
      VeriTcc = cor(as.numeric(veridata$sim), as.numeric(veridata$obs), method="pearson")

      #if(!is.na(TrnTcc) & !is.na(VeriTcc) & TrnTcc>TrnCR & VeriTcc>VeriCR & (TrnTcc * VeriTcc) > 0){
      # If TCC signs are same for training and verification periods
      if(!is.na(TrnTcc) & !is.na(VeriTcc) & (TrnTcc * VeriTcc) > 0){
        bcnt = bcnt + 1
        if(bcnt == 1){
          bestptr = cbind(varnm, curmon, ltime, ptrnm, dirnm, TrnTcc, VeriTcc)
          bcnt = bcnt + 1
        } else {
          besttmp = cbind(varnm, curmon, ltime, ptrnm, dirnm, TrnTcc, VeriTcc)
          bestptr = rbind(bestptr, besttmp)
          bcnt = bcnt + 1
        }
      }

      cat(sprintf("     MWRObs: Sign of TCC for training and verification periods has been checked: %s\n", fname))

    } # Loop for File Count

    if(bcnt > 0){
      bestptr = as.data.frame(bestptr)
      write.csv(bestptr, paste(outdir, "/selected_with_same_sign.csv", sep=""), row.names=F)
    }
  }


  ############ Select best with highest TCC for training period
  ChkDFile2 = paste(outdir, "/selected_with_highest_TCC.csv", sep="")
  if(file.exists(ChkDFile1) & !file.exists(ChkDFile2)){

    varnms = unique(bestptr$varnm); varcnt = length(varnms)
    monthnms = unique(bestptr$curmon); moncnt = length(monthnms)
    ltnms = unique(bestptr$ltime); ltcnt = length(ltnms)
    ptrnms = unique(bestptr$ptrnm); ptrcnt = length(ptrnms)
    dirnms = unique(bestptr$dirrnm); dircnt = length(dirnms)

    bcnt = 1
    for(i in 1:varcnt){
      for(j in 1:moncnt){
        for(ii in 1:ltcnt){

          varnm = varnms[i]; curmon = monthnms[j]; ltime = ltnms[ii]
          pptrs = bestptr[which(bestptr$varnm == varnm & bestptr$curmon == curmon  & bestptr$ltime == ltime),]

          if(nrow(pptrs) == 0){
            bestmodel = cbind(varnm, curmon, ltime, NA, NA, NA, NA)
            colnames(bestmodel) = c("varnm", "curmon", "ltime", "ptrnm", "dirnm", "TrnTcc", "VeriTcc")
          } else {
            pptrs = pptrs[order(pptrs$TrnTcc, decreasing=T),]
            bestmodel = pptrs[1,]
          }

          if(bcnt == 1){
            bestout = bestmodel
            bcnt = bcnt + 1
          } else {
            besttmp = bestmodel
            bestout = rbind(bestout, besttmp)
            bcnt = bcnt + 1
          }
        }
      }
    }

    write.csv(bestout, paste(outdir, "/selected_with_highest_TCC.csv", sep=""), row.names=F)
    cat(sprintf("     MWRObs: Predictor with highest TCC for training period has been selected!\n"))

  }

  ######## Update time-series of the best predictors
  if(file.exists(ChkDFile2)){
    bestout = read.csv(paste(outdir, "/selected_with_highest_TCC.csv", sep=""), header=T)
    rowcnt = nrow(bestout)
    for(i in 1:rowcnt){

      varnm = bestout[i,c("varnm")]; cmonnm = bestout[i,c("curmon")]
      ltime = bestout[i,c("ltime")]; ptrnm = bestout[i,c("ptrnm")]; dirnm = bestout[i,c("dirnm")]

      if(!is.na(ptrnm)){

        # Copy time-series from 01 to 02_ folder
        fname = paste("MWRObs-", varnm, "-", cmonnm, "-", ltime, "-", ptrnm, "-", dirnm, "-TSeries.csv", sep="")
        fromname = paste(indir, "/", fname, sep="")
        toname = paste(outdir, "/", fname, sep="")
        file.copy(fromname, toname, overwrite=T)

        # Read old time-series from 02_BPredictor_TSeries folder
        TsDFile = paste(outdir, "/", fname, sep="")
        tsdata = read.csv(TsDFile, header=T)
        tsdata = tsdata[1:(nrow(tsdata)-2),]
        syear_mme = as.numeric(substr(tsdata[nrow(tsdata), c("yearmon")],1,4)) + 1

        #mdldir = paste(mmedir, "/", mmetype, "/", mdlnm, sep="")
        cmonth = as.numeric(substr(cmonnm,2,3))
        #curmon = monnms[cmonth]
        LTime = as.numeric(substr(ltime, 3, 4))
        #ptrnms = ptrnm
        ptrinfo = paste(rnldir, "/", ptrfile, sep="")
        ptrlist = read.csv(ptrinfo, header=T)
        ptrloc = which(ptrlist$ptrnm == ptrnm)


        if(syear_mme <= eyear_sim){

          mme = MWRObs.Extract.Predictor.Grid(rnldir, ptrfile, cmonth, LTime, syear_mme, eyear_sim)
          mme = mme[ , minrow:maxrow, , ptrloc]

          # Read bounday file
          BndDFile = paste(indir, "/MWRObs-", varnm, "-", cmonnm, "-", ltime, "-", ptrnm, "-", dirnm, ".csv", sep="")
          bbnd = as.data.frame(read.csv(BndDFile, header=T)/100)

          yearmon = as.data.frame(paste(seq(syear_mme, eyear_sim), "-", sprintf("%02d", cmonth), sep=""))
          colnames(yearmon) = "yearmon"
          nyr = length(mme[1,1,])
          ptrdata = array(NA,dim=list(nyr))
          for(tt in 1:nyr){
            bmme = mme[,,tt] * bbnd
            bmme[bmme == 0] = NA
            ptrdata[tt] = mean(unlist(bmme), na.rm=T)
          }
          ptrdata = as.data.frame(ptrdata)
          colnames(ptrdata) = c("APCC")
          tsadd = cbind(yearmon, ptrdata )
          tsadd$Variable = varnm
          #tsadd$model = mdlnm
          tsadd$LTime = LTime
          tsadd$predictor = ptrnm
          tsadd$month = cmonth

          tsout = rbind(tsdata, tsadd)
          write.csv(tsout, TsDFile, row.names=F)
        }

      }

    }

    cat(sprintf("     MWR: Time-series of best predictor of MWR have been updated.\n"))

  }

}

#' @export
MWRObs.Run.Best.Regression <- function(EnvList, tyear) {

  #prjdir, rnldir, ptrfile, vardir, varfile, syear_mme, eyear_obs, eyear_sim, tyear
  prjdir = EnvList$prjdir
  rnldir = EnvList$rnldir
  ptrfile = EnvList$ptrfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  eyear_sim = as.numeric(EnvList$eyear_sim)
  #mmetype = "6MON"

  options(stringsAsFactors = FALSE)

  #indir = paste(prjdir, "/MWRObs/", mmetype, "/", tyear, "/02_BPredictor_TSeries", sep="")
  #outdir = paste(prjdir, "/MWRObs/", mmetype, "/", tyear, "/03_Run_Regression", sep="")
  indir = paste(prjdir, "/MWRObs/", tyear, "/02_BPredictor_TSeries", sep="")
  outdir = paste(prjdir, "/MWRObs/", tyear, "/03_Run_Regression", sep="")
  SetWorkingDir(outdir)

  srchstr = "MWRObs-*TSeries.csv"
  flist = list.files(indir, pattern = glob2rx(srchstr), full.names = F)

  fcnt = length(flist)
  if(fcnt > 0){
    for(i in 1:fcnt){

      fname = flist[i]
      TsDFile = paste(indir, "/", fname, sep="")

      ptr = read.csv(TsDFile, header=T)
      varnm = unique(ptr$Variable)
      #mdlnm = unique(ptr$model)
      ptrnm = unique(ptr$predictor)
      cmonth = as.numeric(unique(ptr$month))
      ptrdata = ptr[c("yearmon", "APCC")]

      VarDFile = paste(vardir, "/", varfile, sep="")
      var = read.csv(VarDFile, header=T)
      vardata = var[c("yearmon", varnm)]

      data = merge(vardata, ptrdata, by="yearmon", all=T)

      regtbl = data[which(as.numeric(substr(data$yearmon,1,4)) >= syear_mme & as.numeric(substr(data$yearmon,1,4)) <= eyear_obs & as.numeric(substr(data$yearmon,6,7)) == cmonth), c("yearmon", varnm, "APCC")]

      cnt=1

      ############ Observed period  ###############################
      for(j in syear_mme:eyear_obs){

        curx = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) == j), c("APCC")]

        cvaldata = regtbl[which(!as.numeric(substr(regtbl$yearmon,1,4)) == j), c("yearmon", varnm, "APCC")]
        #Exclude target year
        cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]

        cvaldata = cvaldata[,c(varnm, "APCC")]
        colnames(cvaldata) = c("y", "x")
        fit = lm(y ~ x, data=cvaldata)

        summary = summary(fit)
        coeff = coefficients(fit)
        R2 = summary$r.squared
        sigma = summary$sigma
        pVal = anova(fit)$'Pr(>F)'[1]
        if(any(is.na(curx))){
          cury = NA
        } else {
          cury = coeff[1] + coeff[2] * curx
        }

        if(cnt == 1){
          valout = cbind(varnm, j, cmonth, coeff[1], coeff[2], R2, sigma, pVal, curx, cury)
          colnames(valout) = c("name", "RmYear", "month", "intercept", "slope", "R2", "Sigma", "p-value", "X", "estimation")
          cnt = cnt + 1
        } else {
          valtmp = cbind(varnm, j, cmonth, coeff[1], coeff[2], R2, sigma, pVal, curx, cury)
          colnames(valtmp) = c("name", "RmYear", "month", "intercept", "slope", "R2", "Sigma", "p-value", "X", "estimation")
          valout = rbind(valout, valtmp)
          cnt = cnt + 1
        }

      }

      ############ split simulation ###############################

      regtbl = data[which(as.numeric(substr(data$yearmon,1,4)) >= syear_mme & as.numeric(substr(data$yearmon,1,4)) <= eyear_sim & as.numeric(substr(data$yearmon,6,7)) == cmonth), c("yearmon", varnm, "APCC")]
      # Exclude target year
      regtbl = regtbl[which(!(as.numeric(substr(regtbl$yearmon,1,4)) == tyear)), ]

      for(k in (eyear_obs+1):eyear_sim){

        #curx = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) == k), c("APCC")]
        curx = data[which(as.numeric(substr(data$yearmon,1,4)) == k & as.numeric(substr(data$yearmon,6,7)) == cmonth), c("APCC")]

        cvaldata = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) >= syear_mme & as.numeric(substr(regtbl$yearmon,1,4)) <= eyear_obs), c(varnm, "APCC")]
        colnames(cvaldata) = c("y", "x")
        fit = lm(y ~ x, data=cvaldata)

        summary = summary(fit)
        coeff = coefficients(fit)
        R2 = summary$r.squared
        sigma = summary$sigma
        pVal = anova(fit)$'Pr(>F)'[1]
        if(any(is.na(curx))){
          cury = NA
        } else {
          cury = coeff[1] + coeff[2] * curx
        }

        if(cnt == 1){
          valout = cbind(varnm, k, cmonth, coeff[1], coeff[2], R2, sigma, pVal, curx, cury)
          colnames(valout) = c("name", "RmYear", "month", "intercept", "slope", "R2", "Sigma", "p-value", "X", "estimation")
          cnt = cnt + 1
        } else {
          valtmp = cbind(varnm, k, cmonth, coeff[1], coeff[2], R2, sigma, pVal, curx, cury)
          colnames(valtmp) = c("name", "RmYear", "month", "intercept", "slope", "R2", "Sigma", "p-value", "X", "estimation")
          valout = rbind(valout, valtmp)
          cnt = cnt + 1
        }

      }

      varnm = unlist(strsplit(fname, "-"))[2]
      curmon = unlist(strsplit(fname, "-"))[3]
      #mdlnm = unlist(strsplit(fname, "-"))[4]
      ltime = unlist(strsplit(fname, "-"))[4]

      outfname = paste("MWRObs-", varnm, "-", curmon, "-", ltime, ".csv", sep="")
      OutDFile = paste(outdir, "/", outfname, sep="")
      write.csv(valout, OutDFile, row.names=F)

    }

  }

}

#' @export
MWRObs.Draw.Best.Predictor.Map <- function(EnvList, tyear){

  #prjdir, mmetype, tyear, minlat, maxlat
  prjdir = EnvList$prjdir
  minlat = as.numeric(EnvList$minlat)
  maxlat = as.numeric(EnvList$maxlat)
  #mmetype = "6MON"


  monnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  #indir = paste(prjdir, "/MWRObs/", mmetype, "/", tyear, "/01_Regression-all", sep="")
  #bestdir = paste(prjdir, "/MWRObs/", mmetype, "/", tyear, "/02_BPredictor_TSeries", sep="")
  indir = paste(prjdir, "/MWRObs/", tyear, "/01_Regression-all", sep="")
  bestdir = paste(prjdir, "/MWRObs/", tyear, "/02_BPredictor_TSeries", sep="")

  srchstr = paste("*-TSeries.csv", sep="")
  flist = list.files(bestdir, pattern = glob2rx(srchstr), full.names = F)
  fcnt = length(flist)

  if(fcnt > 0){
    for(i in 1:fcnt){

      fname = flist[i]

      varnm = unlist(strsplit(fname, "-"))[2]
      monnm = monnms[as.numeric(substr(unlist(strsplit(fname, "-"))[3], 2, 3))]
      #mdlnm = unlist(strsplit(fname, "-"))[4]
      ltime = as.numeric(substr(unlist(strsplit(fname, "-"))[4], 3, 4))
      ptrnm = unlist(strsplit(fname, "-"))[5]
      signnm = unlist(strsplit(fname, "-"))[6]

      InDFile = paste(indir, "/",  substr(fname, 1, nchar(fname)-12), ".csv", sep="")
      PngDFile = paste(bestdir, "/",  substr(fname, 1, nchar(fname)-12), ".png", sep="")

      png(PngDFile,width=800,height=(400*(maxlat-minlat)/180+60),bg="white")

      data = read.csv(InDFile, header=T)
      rdata = t(data)
      rdata = rdata[nrow(rdata):1,]
      rdata2 = rdata[,c(137:144,1:136)]

      r = raster(rdata2)
      image(r, zlim=c(100,100))
      prj = "+proj=longlat +datum=WGS84 +no_defs"
      projection(r) = prj
      xmin(r) = -20 -1.25
      xmax(r) = 340 + 1.25
      ymin(r) = minlat -1.25
      ymax(r) = maxlat + 1.25

      par(mar=c(1, 1, 2, 0))
      par(oma=c(1,1,1,1))
      image(r, zlim=c(100,100))

      world.map("world", col="transparent",bg="gray",fill=TRUE, ylim=c(max(minlat,-60),maxlat), mar=c(0,0,0,0), add=T)
      title(main = paste("Best predictor (", ptrnm, ":", signnm, ") for Var=", varnm, " Month=", monnm, " LagTime=", ltime, " months", sep=""), font.main = 4, line=1)

      dev.off()


    }

  }

}

#' @export
MWRObs.Create.Summary.Table <- function(EnvList){

  prjdir = EnvList$prjdir
  rnldir = EnvList$rnldir
  ptrfile = EnvList$ptrfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  minlat = as.numeric(EnvList$minlat)
  maxlat = as.numeric(EnvList$maxlat)
  AcuMonths = as.numeric(EnvList$AcuMonths)
  precopt = as.logical(EnvList$precopt)
  CRAdj = as.double(EnvList$CRAdj)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)

  monthnms = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

  options(warn=-1)
  options(stringsAsFactors = FALSE)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  #ltcnt = as.numeric(substr(mmetype, 1, 1))
  ltcnt = 6

  # Empty existing summary files
  outdir = paste(prjdir, "/0_Analysis/MWRObs/", sep="")
  SetWorkingDir(outdir)
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  monoutdir = paste(prjdir, "/0_Analysis/MWRObs/Monthly-TSeries", sep="")
  SetWorkingDir(monoutdir)
  flist = list.files(monoutdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  cmpoutdir = paste(prjdir, "/0_Analysis/MWRObs/Composite-Grid", sep="")
  SetWorkingDir(cmpoutdir)
  flist = list.files(cmpoutdir, pattern = glob2rx("MWRObs*.*"), full.names = T)
  file.remove(flist)


  if(AcuMonths > ltcnt){
    AcuMonths = ltcnt
    cat(sprintf("     MWR: Accumulation month is greater than max. lead time!"))
  }

  for(i in 1:varcnt){

    #for(k in 1:ltcnt){
    varnm = varnms[i]
    #lagtime = k

    ### Combine forecasted output into table and convert unit
    smry = MWRObs.Combine.Forecast.Output (EnvList, varnm, lagtime=1)
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
        fname = paste("MWRObs-", varnm, sprintf("-AM%02d",acumon), "-Summary.csv", sep="")
        OutDFile = paste(outdir, "/", fname, sep="")
        write.csv(acusmry, OutDFile, row.names=F)

        for(j in smonth:emonth){

          ctdata = acusmry[which(as.numeric(substr(acusmry$yearmon, 6,7)) == j), ]

          #colcnt = ncol(ctdata) - 2 # include MME
          colcnt = ncol(ctdata) - 3
          outdata = as.data.frame(ctdata$yearmon); colnames(outdata) = "yearmon"
          for(jj in 2:colcnt){
            colnm = names(ctdata)[jj]
            #mdlnm = substr(colnm, 1, (nchar(colnm)-2))
            monnm = sprintf("M%02d", j)
            ltnm = paste("LT", substr(colnm, (nchar(colnm)-1), nchar(colnm)), sep="")
            errdata = ctdata[,c(colnm,"OBS")]

            ###### Decision rule
            rsltthold = nrow(errdata[!is.na(errdata$OBS),]) * 0.8
            varthold = nrow(errdata[!is.na(errdata$OBS),]) * 0.5
            errdata = na.omit(errdata)
            colnames(errdata) = c("SIM", "OBS")

            cor = format(cor(errdata$SIM, errdata$OBS, method="pearson"), digits=2)

            CR = critical.r(nrow(errdata))

            if(nrow(errdata) >= rsltthold & cor >= CR * CRAdj){

              srchstr2 = paste("MWRObs-", varnm, "-", monnm, "-", ltnm, "*.png", sep="")
              #srchdir2 = paste(prjdir, "/MWRObs/", mmetype, sep="")
              srchdir2 = paste(prjdir, "/MWRObs/", sep="")
              flist = list.files(srchdir2, pattern = glob2rx(srchstr2), recursive = T)
              fnames = matrix(unlist(strsplit(flist[], "/")), nrow=length(flist), byrow=T)[, 3]
              vnames = matrix(unlist(strsplit(fnames[], "-")), nrow=length(fnames), byrow=T)[, 5]
              vnmtbl = as.data.frame(table(vnames))
              vnmtbl = vnmtbl[order(vnmtbl$Freq, decreasing=T), ]
              mainvarcnt = vnmtbl[1,c("Freq")]
              mainvarnm = vnmtbl[1,c("vnames")]
              vnmcnt = nrow(vnmtbl)
              # correlation direction
              #cdirnm = substr(matrix(unlist(strsplit(fnames[], "-")), nrow=length(fnames), byrow=T)[1, 6],1,8)
              imsi = matrix(unlist(strsplit(fnames[], "-")), nrow=length(fnames), byrow=T)
              cdirnm = substr(imsi[which(imsi[, 5] == mainvarnm), 6],1,8)[1]

              if(mainvarcnt >= varthold & vnmcnt <= 3){
                outtmp = ctdata[, c("yearmon", colnm)]
                outdata = merge(outdata, outtmp, by="yearmon", all=T)

                ######## Copy Composite Grid Files
                srchstr2 = paste("MWRObs-", varnm, "-", monnm, "-", ltnm, "-", mainvarnm, "-", cdirnm, ".csv", sep="")
                #srchdir2 = paste(prjdir, "/MWRObs/", mmetype, sep="")
                srchdir2 = paste(prjdir, "/MWRObs/", sep="")

                ## When output file does not exists
                flist = list.files(srchdir2, pattern = glob2rx(srchstr2), recursive = T)
                for(kk in 1:length(flist)){
                  GrdDFile  = paste(srchdir2, "/", flist[kk], sep="")
                  data = read.csv(GrdDFile, header=T)
                  if(kk == 1){
                    grdout = data/100
                  } else {
                    grdout = grdout + data/100
                  }
                }

                OutDFile = paste(cmpoutdir, "/", srchstr2, sep="")
                write.csv(grdout, OutDFile, row.names=F)

                # draw map
                fname = srchstr2

                varnm = unlist(strsplit(fname, "-"))[2]
                monnm = monthnms[as.numeric(substr(unlist(strsplit(fname, "-"))[3], 2, 3))]
                #mdlnm = unlist(strsplit(fname, "-"))[4]
                ltime = as.numeric(substr(unlist(strsplit(fname, "-"))[4], 3, 4))
                ptrnm = unlist(strsplit(fname, "-"))[5]
                signnm = unlist(strsplit(fname, "-"))[6]

                PngDFile = paste(cmpoutdir, "/",  substr(fname, 1, nchar(fname)-4), ".png", sep="")

                png(PngDFile,width=800,height=(400*(maxlat-minlat)/180+60),bg="white")

                data = grdout
                rdata = t(data)
                rdata = rdata[nrow(rdata):1,]
                rdata2 = rdata[,c(137:144,1:136)]

                r = raster(rdata2)
                image(r, zlim=c(100,100))
                prj = "+proj=longlat +datum=WGS84 +no_defs"
                projection(r) = prj
                xmin(r) = -20 -1.25
                xmax(r) = 340 + 1.25
                ymin(r) = minlat -1.25
                ymax(r) = maxlat + 1.25

                par(mar=c(1, 1, 2, 0))
                par(oma=c(1,1,1,1))
                image(r, zlim=c(1,40))

                world.map("world", col="transparent",bg="gray",fill=TRUE, ylim=c(max(minlat,-60),maxlat), mar=c(0,0,0,0), add=T)
                title(main = paste("Best predictor (", ptrnm, ":", signnm, ") for Var=", varnm, " Month=", monnm, " LeadTime=", ltime, " months", sep=""), font.main = 4, line=1)

                dev.off()

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

            outfile = paste("MWRObs-", varnm, sprintf("-M%02d", j), sprintf("-AM%02d",acumon), ".csv", sep="")
            MonDFile = paste(monoutdir, "/", outfile, sep="")
            write.csv(outdata, MonDFile, row.names=F)

          } # End IF

        } # Month Loop

      } # Accumulation Loop


    } # If there is no output

  } # Variable Loop

}

#' @export
MWRObs.Combine.Forecast.Output <- function(EnvList, varnm, lagtime){

  prjdir = EnvList$prjdir
  rnldir = EnvList$rnldir
  ptrfile = EnvList$ptrfile
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)

  lagcnt = 8

  syear = syear_mme
  eyear = as.numeric(substr(Sys.Date(),1,4))+1

  cnt = 1
  out = NA
  for(i in lagtime:lagcnt){

    mcnt = 1
    for(ii in syear:eyear){

      bcadir = paste(prjdir, "/MWRObs/", ii, "/03_Run_Regression", sep="")
      flist = list.files(bcadir, pattern = glob2rx("*.csv"))

      if(length(flist) > 0){

        for(k in smonth:emonth){

          srchstr = paste("MWRObs-", varnm, sprintf("-M%02d",k), "-", sprintf("LT%02d",i), ".csv", sep="")
          OutDFile = list.files(bcadir, pattern = glob2rx(srchstr), full.names = T)

          if(length(OutDFile) == 1){

            data = read.csv(OutDFile, header=T)
            data$yearmon = sprintf("%s-%02d", data$RmYear, data$month)
            data = data[which(data$RmYear == ii), c("yearmon", "estimation")]

            if(mcnt == 1){
              mout = data
              mcnt = mcnt + 1
            } else {
              mtmp = data
              mout = rbind(mout, mtmp)
              mcnt = mcnt + 1
            }
          } # if(length(OutDFile) = 1)

        } # Month Loop

      } # IF file.exists
    } # Year Loop

    if(cnt == 1 & mcnt > 1){
      colnames(mout) = c("yearmon", sprintf("LTime%02d",i))
      out = mout
      cnt = cnt + 1
    } else if(cnt > 1 & mcnt > 1) {
      colnames(mout) = c("yearmon", sprintf("LTime%02d",i))
      tmp = mout
      out = merge(out, tmp, by="yearmon", all=T)
      cnt = cnt + 1
    }

  } # LeadTime Loop


  # Overall MME by considering different models and lead times
  if(!is.na(out)){
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
    obs = obs[which(as.numeric(substr(obs$yearmon,1,4)) >= syear_mme), ]
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



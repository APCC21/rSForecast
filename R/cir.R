#' @export
CIReg.Copy.Hindcast.Results <- function(EnvList, tyear){

  prjdir = EnvList$prjdir
  eyear_obs = as.numeric(EnvList$eyear_obs)

  if(tyear > eyear_obs){

    srcdir = paste(prjdir, "/CIReg/", eyear_obs, "/01_Regression-all", sep="")
    dstdir = paste(prjdir, "/CIReg/", tyear, sep="")
    SetWorkingDir(dstdir)

    file.copy(srcdir, dstdir, recursive = T)

  }

}

#' @export
CIReg.Create.Summary.Table <- function(EnvList){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  AcuMonths = as.numeric(EnvList$AcuMonths)
  precopt = as.logical(EnvList$precopt)
  CRAdj = as.double(EnvList$CRAdj)


  options(warn=-1)
  options(stringsAsFactors = FALSE)

  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  ltcnt = 6

  # remove all existing summary files
  outdir = paste(prjdir, "/0_Analysis/CIReg", sep="")
  SetWorkingDir(outdir)
  flist = list.files(outdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  monoutdir = paste(prjdir, "/0_Analysis/CIReg/Monthly-TSeries", sep="")
  SetWorkingDir(monoutdir)
  flist = list.files(monoutdir, pattern = glob2rx("*.csv"), full.names = T)
  file.remove(flist)

  if(AcuMonths > ltcnt){
    AcuMonths = ltcnt
    cat(sprintf("     CIR: Accumulation month is greater than max. lead time!"))
  }

  for(i in 1:varcnt){

    #for(k in 1:ltcnt){
    varnm = varnms[i]
    #lagtime = k

    ### Combine forecasted output into table and convert unit
    smry = CIReg.Combine.Forecast.Output(EnvList, varnm, lagtime=1, ltcnt)

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
        fname = paste("CIReg-", varnm, sprintf("-AM%02d",acumon), "-Summary.csv", sep="")
        OutDFile = paste(outdir, "/", fname, sep="")
        write.csv(acusmry, OutDFile, row.names=F)

        for(j in 1:12){

          #ctdata = smry[which(as.numeric(substr(smry$yearmon, 6,7)) == j),c("yearmon", "MME", "OBS")]
          ctdata = acusmry[which(as.numeric(substr(acusmry$yearmon, 6,7)) == j), ]

          colcnt = ncol(ctdata) - 3
          outdata = as.data.frame(ctdata$yearmon); colnames(outdata) = "yearmon"
          for(jj in 2:colcnt){
            colnm = names(ctdata)[jj]
            monnm = sprintf("M%02d", j)
            ltnm = paste("MinLag", substr(colnm, (nchar(colnm)-1), nchar(colnm)), sep="")
            errdata = ctdata[,c(colnm,"OBS")]
            idxthold = nrow(errdata[!is.na(errdata$OBS),]) * 0.5
            errdata = na.omit(errdata)
            colnames(errdata) = c("SIM", "OBS")

            CR = critical.r(nrow(errdata))
            cor = format(cor(errdata$SIM, errdata$OBS, method="pearson"), digits=2)

            if(cor >= CR * CRAdj){

              bmdir = paste(prjdir, "/CIReg", sep="")
              dlist = list.dirs(bmdir, recursive=F, full.names=T)
              dcnt = length(dlist)

              for(k in 1:dcnt){
                dirnm = dlist[k]
                BmdlDFile = paste(dirnm, "/02_Best_multi-model/02_BestModel/", varnm, "-", ltnm, "-Observed.csv", sep="")
                bmdata = read.csv(BmdlDFile, header=T)
                bmdata$sel = paste(bmdata$x1, bmdata$l1, sep="")
                idxsel = bmdata[j, c("sel")]
                if(k == 1){
                  idxout = idxsel
                } else {
                  idxtmp = idxsel
                  idxout = rbind(idxout, idxtmp)
                }
              } # Year Loop

              idxtbl = as.data.frame(table(idxout))
              idxtbl = idxtbl[order(idxtbl$Freq, decreasing=T), ]
              mainidxcnt = idxtbl[1,c("Freq")]
              idxcnt = nrow(idxtbl)

              if(mainidxcnt >= idxthold){
                outtmp = ctdata[, c("yearmon", colnm)]
                outdata = merge(outdata, outtmp, by="yearmon", all=T)
              }

            } # TCC is greater Critical Value

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


            outfile = paste("CIReg-", varnm, sprintf("-M%02d", j), sprintf("-AM%02d",acumon), ".csv", sep="")

            MonDFile = paste(monoutdir, "/", outfile, sep="")
            write.csv(ctdata, MonDFile, row.names=F)

          } # End if

        } # Month Loop

      } # Accumulation Loop

    }

  } # Variable Loop

}

#' @export
CIReg.Monthly.Cross.Split.Validation <- function(VMode, currow, BestModels, regtbl, varnm, blagnms, syear_mme, eyear_obs, eyear_sim, tyear, precopt){

  ptrcnt = ncol(BestModels)/3
  bptrcnt = length(blagnms)

  lagnms = as.vector(t(BestModels[currow,(1:ptrcnt)]))
  lagmons = as.vector(t(BestModels[currow,(1+ptrcnt):(ptrcnt*2)]))
  bestptrs = as.vector(t(BestModels[currow,(1+ptrcnt*2):(ptrcnt*3)]))
  ptrcoeff = bestptrs

  # Check the minimum value of the variable
  minval = min(regtbl[, c(varnm)], na.rm=T)

  if(VMode == "Cross"){
    regtbl = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) >= syear_mme & as.numeric(substr(regtbl$yearmon,1,4)) <= eyear_obs), c("yearmon", varnm, blagnms)]
  } else {
    regtbl = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) >= syear_mme & as.numeric(substr(regtbl$yearmon,1,4)) <= eyear_sim), c("yearmon", varnm, blagnms)]
  }
  colnames(regtbl) = c("yearmon", "var", blagnms)

  curmon = as.numeric(unique(substr(regtbl$yearmon,6,7)))

  #### Cross / Split Validation
  cnt = 1
  if(VMode == "Cross"){
    startyear = syear_mme
    endyear = eyear_obs
  } else {
    startyear = (eyear_obs + 1)
    endyear = eyear_sim
  }
  for(m in startyear:endyear){

    curxs = regtbl[which(as.numeric(substr(regtbl$yearmon,1,4)) == m), c(blagnms)]

    cvaldata = regtbl[which(!as.numeric(substr(regtbl$yearmon,1,4)) == m), c("yearmon", "var", blagnms)]
    #Exclude target year
    cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]
    cvaldata = cvaldata[,c("var", blagnms)]

    fit = lm(var ~ ., data=cvaldata)

    summary = summary(fit)
    coeff = coefficients(fit)
    R2 = summary$r.squared
    sigma = summary$sigma
    pVal = anova(fit)$'Pr(>F)'[1]
    if(any(is.na(curxs)) | length(curxs) == 0){
      cury = NA
    } else {
      cury = coeff[1] + as.vector(coeff[2:(bptrcnt+1)]) %*% as.vector(t(curxs))
    }

    # save estimated regression model parameters
    for(kk in 1:ptrcnt){
      bcoln = which(blagnms[] == lagnms[kk])
      if(length(bcoln) != 0){ ptrcoeff[kk] = coeff[bcoln+1] }
    }

    # change negative values to zero
    if((varnm == "prec" | precopt == T) & (!is.na(cury)) & (cury <= 0)) { cury = 0}

    if(cnt == 1){
      valout = cbind(varnm, m, curmon, t(lagnms), t(lagmons), coeff[1], t(ptrcoeff), R2, sigma, pVal, cury)
      colnames(valout) = c("name", "RmYear", "month", rep(paste("idx",1:ptrcnt, sep="")), rep(paste("lag",1:ptrcnt, sep="")), "intercept", rep(paste("x",1:ptrcnt, sep="")), "R2", "Sigma", "p-value", "estimation")
      cnt = cnt + 1
    } else {
      valtmp = cbind(varnm, m, curmon, t(lagnms), t(lagmons), coeff[1], t(ptrcoeff), R2, sigma, pVal, cury)
      colnames(valtmp) = c("name", "RmYear", "month", rep(paste("idx",1:ptrcnt, sep="")), rep(paste("lag",1:ptrcnt, sep="")), "intercept", rep(paste("x",1:ptrcnt, sep="")), "R2", "Sigma", "p-value", "estimation")
      valout = rbind(valout, valtmp)
      cnt = cnt + 1
    }

  }

  return(valout)

}

#' @export
CIReg.Extract.Best.Model.Parameters <- function(BestModels, curmon){

  BestModels[BestModels == "NA"] = NA

  ptrcnt = ncol(BestModels)/3

  lagnms = as.vector(t(BestModels[curmon,(1:ptrcnt)]))
  lagmons = as.vector(t(BestModels[curmon,(1+ptrcnt):(ptrcnt*2)]))
  bestptrs = as.vector(t(BestModels[curmon,(1+ptrcnt*2):(ptrcnt*3)]))

  # remove predictors which is not selected
  ptrtbl = rbind(lagnms, lagmons, bestptrs)
  ptrtbl = ptrtbl[,colSums(is.na(ptrtbl)) == 0]

  # if one perdictor is selected
  if(is.null(ncol(ptrtbl))){
    bptrcnt = 1
    blagnms = ptrtbl[1]
    blagmons = ptrtbl[2]
  }else{
    bptrcnt = ncol(ptrtbl)
    blagnms = ptrtbl[1,]
    blagmons = ptrtbl[2,]
  }

  outList = list("blagnms"=blagnms, "blagmons"=blagmons, "bestptrs"=bestptrs)

  return(outList)
}

#' @export
CIReg.Select.BestModel.and.Run.Best.Regression <- function(EnvList, tyear){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  idxdir = EnvList$idxdir
  idxfile = EnvList$idxfile
  NBest = as.numeric(EnvList$NBest)
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_mme = as.numeric(EnvList$eyear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  eyear_sim = as.numeric(EnvList$eyear_sim)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)
  precopt = as.logical(EnvList$precopt)

  MaxLTime =6

  eyear_mme

  regdir = paste(prjdir, "/CIReg/", tyear, "/01_Regression-all", sep="")
  bmmdir = paste(prjdir, "/CIReg/", tyear, "/02_Best_multi-model", sep="")
  SetWorkingDir(bmmdir)
  tbldir = paste(bmmdir, "/01_PtrTable", sep="")
  SetWorkingDir(tbldir)
  bestdir = paste(bmmdir, "/02_BestModel", sep="")
  SetWorkingDir(bestdir)

  # read a sample file and get the number of predictors
  ptrcnt = NBest

  # Get variable(predictand) names
  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms = colnames(var)[2:ncol(var)]
  varcnt = length(varnms)

  # Get index names (predictors)
  IdxDFile = paste(idxdir, "/", idxfile, sep="")
  idx = read.csv(IdxDFile, header=T)
  idxnms = colnames(idx)
  idxnms = idxnms[2:length(idxnms)]

  #outfilecheck = paste(outdir, "/", varnms[varcnt], "-MinLag", sprintf("%02d", MaxLTime), "-Observed.csv", sep="")
  #if(!file.exists(outfilecheck)){
  ########## Combine Predictor Tables ############################################
  for(i in 1:varcnt){
    varnm = varnms[i]

    #for(j in 0:MaxLTime){
    for(j in 1:MaxLTime){
      MinLagTime = j

      #### Write Predictor Table
      infile = paste(regdir, "/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")
      outfile = paste(bmmdir, "/01_PtrTable/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")

      data = read.csv(infile, header=T)
      PtrTable = array(NA, dim=c(12, ptrcnt*2))  #Empty array
      colnames(PtrTable) = c(rep(paste("x",1:ptrcnt, sep="")), rep(paste("l",1:ptrcnt, sep="")))
      for(ii in 1:12){
        tem = data[which(data$month == ii), c("index", "lag")]
        for(jj in 1:ptrcnt){
          PtrTable[ii,jj] = as.character(tem$index[jj])
          PtrTable[ii,jj+ptrcnt] = tem$lag[jj]
        }
      }
      write.csv(PtrTable, outfile, row.names=F)

    }
  }
  cat(sprintf("     CIR: Predictors have in combined into table\n"))

  ######## Select Best Model######################################################
  varcnt = length(varnms)
  for(i in 1:varcnt){
    varnm = varnms[i]

    for(j in 1:MaxLTime){
      #for(j in 0:MaxLTime){
      MinLagTime = j

      #### Write Predictor Table
      ptrfile = paste(bmmdir, "/01_PtrTable/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")
      bestfile = paste(bmmdir, "/02_BestModel/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")

      ########
      BestModels = array(NA, dim=c(12, ptrcnt*3))  #Empty array
      colnames(BestModels) = c(rep(paste("x",1:ptrcnt, sep="")), rep(paste("l",1:ptrcnt, sep="")), rep(paste("b",1:ptrcnt, sep="")))

      PtrTable = read.csv(ptrfile, header=T)
      for(ii in smonth:emonth){
        curmon = ii
        lagnms = as.vector(t(PtrTable[ii,(1:ptrcnt)]))
        lagmons = as.vector(t(PtrTable[ii,(1+ptrcnt):(ptrcnt*2)]))

        lagnms2 = lagnms[!is.na(lagnms)]
        lagmons2 = lagmons[!is.na(lagmons)]

        if(length(lagnms2) == 0){

          BestModels[ii,1:ptrcnt] = NA
          BestModels[ii,(ptrcnt+1):(ptrcnt*2)] = NA
          BestModels[ii,(ptrcnt*2+1):(ptrcnt*3)] = ""

        } else {
          #### Combine predictors based on lags
          data = CIReg.Combine.Multiple.Observed.CIndex(VarDFile, IdxDFile, curmon, varnms, lagnms2, lagmons2)

          # Exclude target year
          data = data[which(!(as.numeric(substr(data$yearmon,1,4)) == tyear)), ]


          #data = data[which(as.numeric(substr(data$yearmon,1,4)) >= syear_mme & as.numeric(substr(data$yearmon,1,4)) <= eyear_mme), c(varnm, lagnms2)]
          data = data[which(as.numeric(substr(data$yearmon,1,4)) >= syear_mme & as.numeric(substr(data$yearmon,1,4)) <= eyear_obs), c(varnm, lagnms2)]


          colnames(data) = c("var", lagnms2)

          # leaps : http://www.inside-r.org/packages/cran/leaps/docs/regsubsets
          # http://rstudio-pubs-static.s3.amazonaws.com/2897_9220b21cfc0c43a396ff9abf122bb351.html
          if(length(lagnms2) > 1){
            leaps = regsubsets(var~., data, nbest=1)
            #plot(leaps, scale="adjr2")

            sumry = summary(leaps)
            brow = which.max(sumry$adjr2)
            if(length(brow) == 0){
              bmodel = array("", dim=c(NBest))
            } else {
              bmodel = as.vector(t(as.data.frame(sumry$outmat)[brow,]))
              imsi = array("", dim=c(NBest))
              imsi[match(lagnms2, lagnms)] = bmodel
              bmodel = imsi
            }

            BestModels[ii,1:ptrcnt] = lagnms
            BestModels[ii,(ptrcnt+1):(ptrcnt*2)] = lagmons[1:ptrcnt]
            BestModels[ii,(ptrcnt*2+1):(ptrcnt*3)] = bmodel
          } else {
            BestModels[ii,1:ptrcnt] = lagnms
            BestModels[ii,(ptrcnt+1):(ptrcnt*2)] = lagmons[1:ptrcnt]
            BestModels[ii,(ptrcnt*2+1)] = "*"
          }

        } # If lagnms2 = NA

      } # Month Loop

      write.csv(BestModels, bestfile, row.names=F)
    }
  }
  cat(sprintf("     CIR: Best Multivariate Regression Model has been selected.\n"))

  ######## Run Regression Models #################################################
  varcnt = length(varnms)
  for(i in 1:varcnt){
    varnm = varnms[i]

    for(j in 1:MaxLTime){
      #for(j in 0:MaxLTime){
      MinLagTime = j

      #### Write Predictor Table
      bestfile = paste(bmmdir, "/02_BestModel/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")
      regfile = paste(bmmdir, "/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")

      BestModels = read.csv(bestfile, header=T, na.strings=c(" "))

      cnt=0
      for(ii in smonth:emonth){

        curmon = ii

        # extract best model parameters from Best Model Table
        out = CIReg.Extract.Best.Model.Parameters(BestModels, curmon)
        blagnms = out$blagnms; blagmons = out$blagmons; bestptrs = out$bestptrs

        blagnms[blagnms == "NA"] = NA
        blagmons[blagmons == "NA"] = NA

        blagnms2 = blagnms[!is.na(blagnms)]
        blagmons2 = blagmons[!is.na(blagmons)]

        if(length(blagnms2) > 0) {

          cnt = cnt + 1

          #### Combine predictors and predictands based on lags
          regtbl = CIReg.Combine.Multiple.Observed.CIndex(VarDFile, IdxDFile, curmon, varnms, blagnms2, blagmons2)

          #### Do cross validation
          cvalout = CIReg.Monthly.Cross.Split.Validation(VMode="Cross", curmon, BestModels, regtbl, varnm, blagnms2, syear_mme, eyear_obs, eyear_sim, tyear, precopt=F)

          ### Do split validation
          if(eyear_sim > eyear_mme){
            svalout = CIReg.Monthly.Cross.Split.Validation(VMode="Split", curmon, BestModels, regtbl, varnm, blagnms2, syear_mme, eyear_obs, eyear_sim, tyear, precopt=F)
          }


          if(cnt == 1){
            if(eyear_sim > eyear_mme){
              sumout = rbind(cvalout, svalout)
            } else {
              sumout = cvalout
            }

          } else {
            if(eyear_sim > eyear_mme){
              imsi = rbind(cvalout, svalout)
              sumout = rbind(sumout, imsi)
            } else {
              imsi = cvalout
              sumout = rbind(sumout, imsi)
            }

          }

        } # if(length(blagnms2) > 0)

      } # Loop for month

      cat(sprintf("     CIR: Cross and Split Validation has been finished: Var=%s Lag=%s\n", varnm, MinLagTime))
      write.csv(sumout, regfile, row.names=F)

    }
  }

  #}# IF file exists

}

#' @export
CIReg.Combine.Multiple.Observed.CIndex <- function(VarDFile, IdxDFile, curmon, varnms, lagnms, lagmons) {

  var = read.csv(VarDFile, header=T)
  var = var[c("yearmon", varnms)]

  idx = read.csv(IdxDFile, header=T)
  idx = idx[c("yearmon", lagnms)]

  idxcnt = length(lagnms)
  for(i in 1:idxcnt){

    lagval = as.numeric(lagmons[i])
    colnum = i + 1

    if(lagval < 0){
      srowt = 1- lagval
      erowt = nrow(idx)
      srows = 1
      erows = nrow(idx) + lagval

      idx[1:srowt-1, colnum] =  NA
      idx[srowt:erowt, colnum] =  idx[srows:erows, colnum]

    } else {
      srowt = 1
      erowt = nrow(idx) - lagval
      srows = 1 + lagval
      erows = nrow(idx)

      idx[erowt+1:erows, colnum] =  NA
      idx[srowt:erowt, colnum] =  idx[srows:erows, colnum]
    }

  }

  out = merge(var, idx, by="yearmon", all=T)

  if(curmon <= 12){
    out = out[which(as.numeric(substr(out$yearmon,6,7)) == curmon),]
  }


  return(out)

}

#' @export
CIReg.Calculate.All.Regression.Split.Validation <- function(EnvList, tyear, precopt) {

  #prjdir, vardir, varfile, idxdir, idxfile, NBest, mmetype, syear_mme, eyear_mme, eyear_obs, tyear, smonth, emonth, precopt
  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  idxdir = EnvList$idxdir
  idxfile = EnvList$idxfile
  NBest = as.numeric(EnvList$NBest)
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_mme = as.numeric(EnvList$eyear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)
  smonth = as.numeric(EnvList$smonth)
  emonth = as.numeric(EnvList$emonth)
  precopt = as.logical(EnvList$precopt)


  #MaxLTime = 6; syear_trn=1983; eyear_trn=2007; precopt=F
  #LTime = as.numeric(substr(mmetype, 1, 1))
  #if(LTime == 3) {MaxLag = -6; MaxLTime =3}
  #if(LTime == 6) {MaxLag = -8; MaxLTime =6}
  LTime = 6
  MaxLag = -8; MaxLTime =6

  options(stringsAsFactors = FALSE)

  outdir = paste(prjdir, "/CIReg/", tyear, "/01_Regression-all", sep="")
  SetWorkingDir(outdir)

  # Get variable(predictand) names
  VarDFile = paste(vardir, "/", varfile, sep="")
  var = read.csv(VarDFile, header=T)
  varnms_cir = colnames(var)[2:ncol(var)]

  # Get index names (predictors)
  IdxDFile = paste(idxdir, "/", idxfile, sep="")
  idx = read.csv(IdxDFile, header=T)
  idxnms = colnames(idx)
  idxnms = idxnms[2:length(idxnms)]

  varcnt = length(varnms_cir)
  idxcnt = length(idxnms)
  #moncnt = (emonth - smonth + 1)

  for(i in 1:varcnt){
    varnm = varnms_cir[i]

    for(j in smonth:emonth){

      curmon = j

      outfile = paste(outdir, "/Reg_", varnm, "_", sprintf("%02d", j), "-SplitTest.csv", sep="")

      #if(!file.exists(outfile)){
      CR = critical.r(round((eyear_mme - syear_mme + 1)*2/3, digits = 0))
      #CR = critical.r(eyear_mme - syear_mme + 1)
      MaxCR = CR
      repeat{

        cnt = 1
        #curmon = j

        for(k in 1:idxcnt){
          idxnm = idxnms[k]

          # Need 2 month lag for forecasting next month (data is avilalble at previous month)
          for(ii in MaxLag:-2){
            #for(ii in MaxLag:0){
            lagmon = ii

            #data = CIReg.combine.Obs_Cindex(VarDFile, IdxDFile, curmon, varnm, idxnm, lagmon)
            data = CIReg.Combine.Multiple.Observed.CIndex(VarDFile, IdxDFile, curmon, varnm, idxnm, lagmon)
            data = data[which(as.numeric(substr(data$yearmon,1,4)) >= syear_mme & as.numeric(substr(data$yearmon,1,4)) <= eyear_mme), ]
            # Exclude target_year for cross validation
            data = data[which(!(as.numeric(substr(data$yearmon,1,4)) == tyear)), ]
            data = na.omit(data)
            colnames(data) = c("yearmon", "y", "x")

            if(nrow(data) >= 15){

              SdTCC = 0; MeanTCC = 0

              # 2/3 of data
              cend = round(nrow(data)*2/3, digits=0)
              vend = nrow(data)
              vstart = vend - cend + 1

              cdata = data[1:cend,]
              CTcc = cor(cdata$x, cdata$y, method="pearson")

              vdata = data[vstart:vend,]
              VTcc = cor(vdata$x, vdata$y, method="pearson")

              nrsmpls = (vend-cend) * 3
              for(tt in 1:nrsmpls){
                rdatanm = paste("rdata", tt, sep="")
                RTccnm = paste("RTcc", tt, sep="")
                assign(rdatanm, data[sample(1:vend, cend, replace=F),])
                assign(RTccnm, cor(get(rdatanm)[c("x")], get(rdatanm)[c("y")], method="pearson"))
              }

              # Check the correlation sign is on the sameway and greater than critical TCC.
              ntsmpls = nrsmpls + 2
              tcc = array(0, dim = c(ntsmpls))
              mcnt = 0; pcnt = 0
              if(CTcc > 0 & CTcc >= MaxCR){
                pcnt = pcnt + 1
              } else if(CTcc < 0 & abs(CTcc) >= MaxCR) {
                mcnt = mcnt + 1
              }
              tcc[1] = CTcc

              if(VTcc > 0 & CTcc >= MaxCR){
                pcnt = pcnt + 1
              } else if(CTcc < 0 & abs(CTcc) >= MaxCR)  {
                mcnt = mcnt + 1
              }
              tcc[2] = VTcc

              for(tt in 1:nrsmpls){
                RTccnm = paste("RTcc", tt, sep="")
                if(get(RTccnm) > 0 & get(RTccnm) >= MaxCR) {
                  pcnt = pcnt + 1
                } else if(get(RTccnm) < 0 & abs(get(RTccnm)) >= MaxCR)  {
                  mcnt = mcnt + 1
                }
                tcc[tt+2] = get(RTccnm)
              }


              if(pcnt == ntsmpls | mcnt == ntsmpls){
                SdTCC = sd(tcc)
                MeanTCC = mean(tcc)
                if(cnt == 1){
                  splitout = cbind(varnm, curmon, idxnm, lagmon, SdTCC, MeanTCC)
                  cnt = cnt + 1
                } else {
                  splittmp = cbind(varnm, curmon, idxnm, lagmon, SdTCC, MeanTCC)
                  splitout = rbind(splitout, splittmp)
                  cnt = cnt + 1
                }
              }
            }
          } # End of Loop for lag-time for index
        } # End of Index loop

        if(cnt > 1){

          splitout = as.data.frame(splitout)
          colnames(splitout) = c("name", "month", "index", "lag", "SdTCC", "MeanTCC")
          splitout = splitout[order(abs(as.numeric(splitout$MeanTCC)), decreasing = T),][1:nrow(splitout),]
          #splitout = splitout[order(abs(as.numeric(splitout$MeanTCC)), decreasing = T),][1:MaxLTime,]

          bcnt = 1
          selcnt = nrow(splitout)
          for(k in 1:selcnt){

            idxnm = splitout[k, c("index")]
            lagmon = as.numeric(splitout[k, c("lag")])

            regdata = CIReg.Combine.Multiple.Observed.CIndex(VarDFile, IdxDFile, curmon, varnm, idxnm, lagmon)

            data = regdata[which(as.numeric(substr(regdata$yearmon,1,4)) >= syear_mme & as.numeric(substr(regdata$yearmon,1,4)) <= eyear_mme), ]

            # if there is at least one of NA in predictor, skip the regression
            #if(!any(is.na(data[c(idxnm)]))){
            cnt = 1
            for(jj in syear_mme:eyear_mme){

              curx = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c(idxnm)]

              if(!is.na(curx)){

                cvaldata = data[which(!as.numeric(substr(data$yearmon,1,4)) == jj), c("yearmon", varnm, idxnm)]
                cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]

                colnames(cvaldata) = c("yearmon", "y", "x")
                fit = lm(y ~ x, data=cvaldata)

                summary = summary(fit)
                coeff = coefficients(fit)
                R2 = summary$r.squared
                sigma = summary$sigma
                pVal = anova(fit)$'Pr(>F)'[1]
                cury = coeff[1] + coeff[2]*curx
                if((varnm == "prec" | precopt == T) & (!is.na(cury)) & (cury <= 0)) { cury = 0}
              } else {
                coeff[1:2] = NA; R2 = NA; sigma = NA; pVal = NA; cury = NA
              }

              if(cnt == 1){
                sumout = cbind(varnm, curmon, idxnm, lagmon, jj, t(coeff), R2, sigma, pVal, cury)
                cnt = cnt + 1
              } else {
                sumtmp = cbind(varnm, curmon, idxnm, lagmon, jj, t(coeff), R2, sigma, pVal, cury)
                sumout = rbind(sumout, sumtmp)
                cnt = cnt + 1
              }

            } # End of Year Loop for calibration

            cvaldata = data
            cvaldata = cvaldata[which(!(as.numeric(substr(cvaldata$yearmon,1,4)) == tyear)), ]
            colnames(cvaldata) = c("yearmon", "y", "x")
            fit = lm(y ~ x, data=cvaldata)
            summary = summary(fit)
            coeff = coefficients(fit)
            R2 = summary$r.squared
            sigma = summary$sigma
            pVal = anova(fit)$'Pr(>F)'[1]
            # Remaining period for verification
            data = regdata[which(as.numeric(substr(regdata$yearmon,1,4)) >= (eyear_mme + 1) & as.numeric(substr(regdata$yearmon,1,4)) <= eyear_obs), ]
            for(jj in (eyear_mme + 1):eyear_obs){


              curx = data[which(as.numeric(substr(data$yearmon,1,4)) == jj), c(idxnm)]

              if(!is.na(curx)){
                coeff = coefficients(fit)
                R2 = summary$r.squared
                sigma = summary$sigma
                pVal = anova(fit)$'Pr(>F)'[1]
                cury = coeff[1] + coeff[2]*curx
                if((varnm == "prec" | precopt == T) & (!is.na(cury)) & (cury <= 0)) { cury = 0}
              } else {
                coeff[1:2] = NA; R2 = NA; sigma = NA; pVal = NA; cury = NA
              }

              if(cnt == 1){
                sumout = cbind(varnm, curmon, idxnm, lagmon, jj, t(coeff), R2, sigma, pVal, cury)
                cnt = cnt + 1
              } else {
                sumtmp = cbind(varnm, curmon, idxnm, lagmon, jj, t(coeff), R2, sigma, pVal, cury)
                sumout = rbind(sumout, sumtmp)
                cnt = cnt + 1
              }

            } # End of Year Loop for Cross Validation

            # Check TCC for training and verification period
            sumout = as.data.frame(sumout)
            colnames(sumout) = c("name", "month", "index", "lag", "RmYear", "Intercept", "X", "R2", "Sigma", "p-value", "estimation")
            sumout$yearmon = paste(sumout$RmYear, "-", sprintf("%02d", as.numeric(sumout$month)), sep="")

            sim = sumout[c("yearmon", "estimation")]
            obs = var[c("yearmon", varnm)]
            colnames(obs) = c("yearmon", "observed")

            outdata = merge(sim, obs, by="yearmon")
            trndata = outdata[which(as.numeric(substr(outdata$yearmon,1,4)) >= syear_mme & as.numeric(substr(outdata$yearmon,1,4)) <= eyear_mme), ]
            veridata = outdata[which(as.numeric(substr(outdata$yearmon,1,4)) >= (eyear_mme + 1) & as.numeric(substr(outdata$yearmon,1,4)) <= eyear_obs), ]

            TrnTcc = cor(as.numeric(trndata$estimation), trndata$observed, method="pearson")
            VeriTcc = cor(as.numeric(veridata$estimation), veridata$observed, method="pearson")

            if(TrnTcc * VeriTcc > 0){
              if(bcnt == 1){
                bptr = splitout[k,]
                bcnt = bcnt + 1
              } else {
                btmp = splitout[k,]
                bptr = rbind(bptr, btmp)
                bcnt = bcnt + 1
              }
            }

          } # End of SelCnt

        } else {
          bcnt = 1
        } # IF (cnt > 1)

        if(bcnt >1 | MaxCR < 0){
          #if(bcnt > MaxLTime | MaxCR < CR){

          if(MaxCR < CR){
            cat(sprintf("     CIR: Correlation Coefficient is smaller than Critical CC: Critical=%3.2f, Selected=%3.2f, Year=%d, Month=%s, Variable=%s\n", CR, MaxCR, tyear, curmon, varnm))
          }

          break
        }
        MaxCR = MaxCR - 0.05
      } # End of Repeat

      if(bcnt > 1){
        cat(sprintf("     CIR: Analysis has been finished Current Month=%s, Variable=%s\n", curmon, varnm))
        outfile = paste(outdir, "/Reg_", varnm, "_", sprintf("%02d", j), "-SplitTest.csv", sep="")
        #outfile = paste(outdir, "/Reg_", varnm, "_", sprintf("%02d", j), ".csv", sep="")
        write.csv(bptr, outfile, row.names=F)
      } else {
        cat(sprintf("     CIR: Predictors are not available : Year=%d, Month=%s, Variable=%s\n", tyear, curmon, varnm))
      }

    } # End of Month loop
  } # End of Station Loop

  #### Select the best N models for each lag time #################
  for(ii in 1:varcnt){
    varnm = varnms_cir[ii]

    for(jj in 1:MaxLTime){
      #for(jj in 0:MaxLTime){

      MinLagTime = jj

      #### calculate mean of R2 for each combination
      cnt = 0
      for(j in smonth:emonth){
        #regfile = paste(outdir, "/Reg_", varnm, "_", sprintf("%02d", j), ".csv", sep="")
        regfile = paste(outdir, "/Reg_", varnm, "_", sprintf("%02d", j), "-SplitTest.csv", sep="")
        if(file.exists(regfile)){
          regdata = read.csv(regfile, header=T)
          cnt = cnt + 1

          if(cnt == 1){
            rdata = regdata
          } else {
            temp = regdata
            rdata = rbind(rdata, temp)
          }
        }

      } # Month Loop

      cnt = 0
      for(j in smonth:emonth){
        #mondata = data[which(data$month == j & data$lag <= -MinLagTime), ]

        if(exists("rdata")){

          mondata = rdata[which(rdata$month == j & rdata$lag <= (-MinLagTime+1)), ]

          srtdata = mondata[order(mondata$MeanTCC, decreasing=T),]
          #if same index is slected, exclude it
          srtdata = subset(srtdata, !duplicated(srtdata[,c("index")]))

          if(nrow(srtdata)>0){
            cnt = cnt + 1
            if(cnt == 1){
              bestdata = srtdata[1:NBest,]
              bestdata$rank = seq(1:NBest)
              cnt = cnt + 1
            } else {
              tempdata = srtdata[1:NBest,]
              tempdata$rank = seq(1:NBest)
              bestdata = rbind(bestdata, tempdata)
              cnt = cnt + 1
            }
          }
        } # File exists
      } # Month Loop

      if(exists("bestdata")){
        bestdata$FcstMode = "Obs"

        outfile = paste(outdir, "/", varnm, "-MinLag", sprintf("%02d", MinLagTime), "-Observed.csv", sep="")
        write.csv(bestdata, outfile, row.names=F)
        cat(sprintf("     CIR: Best Index been selected for Variable=%s, MinLagtime=%s\n", varnm, jj))
      }


    } # Lag Loop
  } # Var Loop

}

#' @export
CIReg.Combine.Forecast.Output <- function(EnvList, varnm, lagtime, MaxLTime){

  prjdir = EnvList$prjdir
  vardir = EnvList$vardir
  varfile = EnvList$varfile
  syear_mme = as.numeric(EnvList$syear_mme)
  eyear_obs = as.numeric(EnvList$eyear_obs)

  syear = syear_mme
  eyear = as.numeric(substr(Sys.Date(),1,4))+1
  #eyear = eyear_obs

  # Max LeadTime
  cnt = 1
  out = NA
  for(i in lagtime:MaxLTime){

    mcnt = 1
    for(k in syear:eyear){
      #bcadir = paste(prjdir, "/CIReg/", mmetype, "/", k, "/02_Best_multi-model", sep="")
      bcadir = paste(prjdir, "/CIReg/", k, "/02_Best_multi-model", sep="")
      if(file.exists(bcadir)){
        srchstr = paste(varnm, sprintf("-MinLag%02d",i), "-Observed.csv", sep="")
        OutDFile = list.files(bcadir, pattern = glob2rx(srchstr), full.names = T)

        data = read.csv(OutDFile, header=T)
        data$yearmon = sprintf("%s-%02d", data$RmYear, data$month)
        data = data[which(data$RmYear == k), c("yearmon", "estimation")]

        if(mcnt == 1){
          mout = data
          mcnt = mcnt + 1
        } else {
          mtmp = data
          mout = rbind(mout, mtmp)
        }

      } # IF file.exists
    } # Year Loop

    colnames(mout) = c("yearmon", sprintf("Lag%02d",i))

    if(cnt == 1){
      out = mout
      cnt = cnt + 1
    } else {
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

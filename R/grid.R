#' @export
Aphrodite.Area.Averages.Monthly <- function(EnvList, pointopt, weightopt){

  aphrodir = EnvList$aphrodir
  bnddir = EnvList$bnddir
  bndfile = EnvList$bndfile
  vardir = EnvList$vardir

  #library(ncdf4)
  #library(maptools)
  #library(raster)

  outfile = paste("aphro_", bndfile, "_month.csv", sep="")

  outdir = vardir
  SetWorkingDir(outdir)
  OutDFile = paste(outdir, "/", outfile, sep="")

  if(file.exists(OutDFile)){

    cat(sprintf("File=%s already exists!\n\n", outfile))

  } else {

    # read boundary shape file (gis file should have "ID" column)
    GisDFile = paste(bnddir, "/", bndfile, sep="")
    if(pointopt == 'T') {
      bond = readShapePoints(GisDFile, IDVar="ID")
    } else {
      bond = readShapePoly(GisDFile, IDvar="ID")
    }

    srchstr = paste("APHRO_*monthly*.nc", sep="")
    flist = list.files(aphrodir, pattern = glob2rx(srchstr), full.names = F)
    years = as.numeric(substr(flist, (nchar(flist)-6),(nchar(flist)-3)))
    syear = min(years); eyear = max(years)

    cnt = 1
    for (j in syear:eyear){

      nstn = length(bond$ID)
      df = data.frame(matrix(nrow=0, ncol=nstn))

      srchstr = paste("*monthly*", j, ".nc", sep="")
      ncfile = list.files(aphrodir, pattern = glob2rx(srchstr), full.names = F)

      setwd(aphrodir)
      fin = ncdf4::nc_open(ncfile)
      x = ncvar_get(fin, "longitude")
      y = ncvar_get(fin, "latitude")
      var = ncvar_get(fin,"precip")

      sdate = as.Date(paste(j, "-1-1", sep=""))
      edate = as.Date(paste(j, "-12-31", sep=""))
      yearmon = as.data.frame(substr(seq(sdate, edate, by="mon"),1,7))
      colnames(yearmon) = "yearmon"

      for (k in 1:12) {
        varg = var[,,k]

        r = GIS.ncVar2raster(x, y, varg)

        v <- extract(r, bond, fun=mean, weights=weightopt, small=T)

        df =rbind(df, t(v))

      }

      names(df) = bond$ID
      out = cbind(yearmon, df)

      cat(sprintf("Completed File=%s\n", ncfile))

      if(cnt == 1){
        outdata = out
        cnt = cnt + 1
      } else {
        temp = out
        outdata =rbind(outdata, temp)
        cnt = cnt + 1
      }

    }

    write.csv(outdata, OutDFile, row.names=F)
    cat(sprintf("File=%s has been successfully created!\n\n", outfile))

  }

}


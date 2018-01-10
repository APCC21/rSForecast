# Setting working environment
EnvList = Set.Working.Environment (basedir="G:/SCService-RM", prjname = "Korea_SForecast")

#ghcn.daily.update (EnvList, cntry = "KR")

# Download required climate information
#Update.Climate.Information (EnvList, updatemode = F)

# Preparing observed data
#Cal.Monthly.Station.Mean (EnvList, varnms = c("prec", "t2m"))

# Model construction run for hindcast period
Hindcast.Forecast.Model.Construction (EnvList)

# Combine individual forecast model output
Integrate.Individual.Forecast.Model (EnvList)

library("sforecast")

# Setting working environment
EnvList = Set.Working.Environment (basedir="G:/SCService-RM", prjname = "Korea_SForecast")

# Update required climate information
Update.Climate.Information (EnvList, updatemode = T)

# Update KMA ASOS observed data https://data.kma.go.kr/data/grnd/selectAsosList.do?pgmNo=34
# 기상변수 선택시 "기사" 항목 제외 (최근에 추가된 항목)
KMA.ASOS.Observation.Update (EnvList, fillingmode = T)

# Real-time forecast
Run.RealTime.Forecast.Model (EnvList, fiyearmon = "1990-04")

# Temporal downscaling
Daily.Hourly.Sampling (EnvList, fiyearmon = "1990-04", smplmoncnt=6)

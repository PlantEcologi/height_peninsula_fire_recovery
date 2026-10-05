##################################################
##data prep and model fitting for the joint MODIS NDVI + LiDAR height model
##(focal plots only) - based on fit_model2026.R and fit_model_Landsat2026.R
##################################################

####################
##setup
####################

### Load libraries
libs=c(
  "doParallel",
  "sf",
  "knitr",
  "rmarkdown",
  "ggplot2",
  "tidyr",
  "tidyverse",
  "lubridate",
  "rjags",
  "dclone",
  "terra")
lapply(libs, require, character.only=T)

#file locations and names
mdatwd <- "data/"
mname <- "peninsulaHeightMODIS2026" #model name for file naming

# Calculate the number of cores
no_cores <- detectCores() - 1

# # Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cl)

###########################################################
###Load data
###########################################################

#set extent and projection
cp <- ext(18.37, 18.48, -34.362, -34.19)
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Get points of interest, including the LiDAR height data (one column per survey year)
height_years <- c("2015", "2017", "2021", "2023")
pts <- read.csv("data/focal_plots_n50_coords_height.csv", check.names = FALSE) |>
  mutate(geometry = str_remove_all(geometry, pattern = "[c()]")) |>
  separate_wider_delim(cols = geometry, names = c("Longitude", "Latitude"), delim = ",") |>
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) |>
  rename(Site = name)

#modis ndvi from MODIS/061/MOD13Q1
NDVI <- rast("data/NDVI_stack_2001_2026.tif") #new NDVI data
crs(NDVI) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
NAflag(NDVI) <- -9999
NDVI <- trim(NDVI)
NDVI <- project(NDVI, y = crs(projection)) #Make sure it's the same projection
NDVI <- crop(NDVI, cp)
NDVIdates <- read.csv("data/NDVI_dates_2001_2026.csv") # Load dates
time(NDVI) <- as.Date(NDVIdates$date) # Assign time
NDVI <- NDVI[[time(NDVI) < as.Date("2022-06-01")]] # Select only NDVI within the fire observation period

#spatial covariates
vegtype <- vect("data/veg/Vegetation_Indigenous_Remnants.shp") |>
  project(crs(projection)) |>
  crop(cp)
vegtype$National_ <- str_replace_all(vegtype$National_,
                                    c("Beach" = "beach",
                                      "Cape Flats Dune Strandveld - False Bay" = "dune",
                                      "Cape Lowland Freshwater Wetlands" = "wetland",
                                      "Hangklip Sand Fynbos" = "sand",
                                      "Peninsula Granite Fynbos - South" = "granite",
                                      "Peninsula Sandstone Fynbos" = "sandstone"))
vegtype <- rasterize(vegtype, NDVI, field="National_", background=NA, na.rm = T)
names(vegtype) <- "vegtype"

sta <- rast("data/STATIC_stack.tif")
NAflag(sta) <- -9999
sta <- trim(sta)
crs(sta) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
sta <- project(sta, NDVI) #Make sure it's the same projection
sta$northness <- cos(sta$aspect*pi/180)
sta$eastness <- sin(sta$aspect*pi/180)
sta$vegtype <- vegtype

#vegetation age (days since fire; one layer per NDVI date, layer names are dates)
rfi <- rast("data/veldage2022.tif") #new fire age data
rfi <- project(rfi, NDVI, method = "near")

###########################################################
###Extract raster data at the focal points
###########################################################

# NDVI time-series
ndat <- extract(NDVI, pts, method = "simple", na.rm = F) |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!c("ID")) |>
  pivot_longer(cols = -c("UIJ"), names_to = "Date", values_to = "NDVI") |>
  mutate(Date = as.Date(substr(Date,1,10), format = "%Y_%m_%d"))

# Fire age time-series
fdat <- extract(rfi, pts, method = "simple", na.rm = F) |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!c("ID")) |>
  pivot_longer(cols = -c("UIJ"), names_to = "Date", values_to = "Age") |>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

cdat <- inner_join(na.omit(ndat), na.omit(fdat), by=c("UIJ","Date"))

# Covariates
cov <- extract(sta, pts, method = "simple", na.rm = F) |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!c("ID")) |>
  mutate(UI = paste(pts$Longitude, pts$Latitude, sep = "_")) |>
  mutate(vegtype = as.character(vegtype))

###########################################################
###Height data: long format with fire age at the time of each LiDAR survey
###########################################################

# Age (days) at the survey: the age on the last available date on or before
# 1 July of the survey year plus the days elapsed since then (valid as long as
# there is no further fire, i.e. beyond the end of the fire data)
age_dates <- as.Date(names(rfi))
age_mat <- as.data.frame(extract(rfi, pts, method = "simple", na.rm = F))[,-1]

hdat <- pts |>
  st_drop_geometry() |>
  select(UIJ = plot, all_of(height_years)) |>
  pivot_longer(cols = all_of(height_years), names_to = "Year", values_to = "Height") |>
  mutate(Year = as.integer(Year),
         SurveyDate = as.Date(paste0(Year, "-07-01"))) |>
  rowwise() |>
  mutate(Age = {
    p <- match(UIJ, pts$plot)
    k <- max(which(age_dates <= SurveyDate))
    age_mat[p, k] + as.numeric(SurveyDate - age_dates[k])
  }) |>
  ungroup() |>
  mutate(DA = Age/365.25)

###########################################################
###Filter by a few NB criteria and trim covariates to match
###########################################################

#drop NDVI values < 0 or points with <50 NDVI observations
cdat <- cdat |> filter(NDVI > 0) |>
  group_by(UIJ) |>
  filter(n() >= 50) |>
  ungroup()

#trim covariates and temporal data to match
cov <- cov |>
  filter(complete.cases(cov)) |>
  filter(vegtype %in% c("dune", "granite", "sandstone", "sand")) |> #drop beach and wetland veg types
  filter(UIJ %in% cdat$UIJ) |>
  droplevels()

hdat <- hdat |> filter(!is.na(Height), !is.na(DA), DA > 0)

#only keep sites that have both NDVI and height data
keep <- intersect(intersect(cov$UIJ, cdat$UIJ), hdat$UIJ)
cov <- filter(cov, UIJ %in% keep)
cdat <- filter(cdat, UIJ %in% keep)
hdat <- filter(hdat, UIJ %in% keep)

rm(list = c("NDVI", "NDVIdates", "sta", "rfi", "vegtype", "libs", "ndat", "fdat", "age_mat"))
gc()

###########################################################
###Add columns for month of fire and age in years
###########################################################

cdat <- cdat |> mutate(firemonth = month(Date - Age))
cdat$DA <- cdat$Age/365.25

###########################################################
###Quick look at raw data
###########################################################

cdat |>
  ggplot(aes(x = DA, y = NDVI)) +
  geom_point() +
  geom_point(data = hdat, aes(y = Height), colour = "red") +
  facet_wrap(~UIJ)

###########################################################
### Create dummy variables for veg type and select and scale environmental data
###########################################################

tveg <- as.numeric(as.factor(cov$vegtype)) - 1
dummies <- model.matrix(~as.factor(tveg))
dummies <- dummies[,-1, drop = FALSE]

#select vars and scale data
envars <- c("elevation", "slope", "tpi", "northness", "eastness")
scaled <- scale(as.matrix(cov[,envars]))
env <- as.data.frame(cbind(intercept=1, scaled, tveg))
env$UI <- cov$UI
env$UIJ <- cov$UIJ
env <- cbind(env,dummies)

#save the scaling parameters to convert fitted coefficients back to metric units later
beta.mu=c(intercept=0,attr(scaled,"scaled:center"))
beta.sd=c(intercept=1,attr(scaled,"scaled:scale"))
rm(scaled)

###########################################################
###Format data for JAGS
###########################################################

tdat <- cdat
tdat$NDIN <- tdat$NDVI

#create new id that goes from 1 to nGrid (to order env and the temporal data in the same way)
env$jag_id <- as.integer(as.factor(env$UIJ))
jtab <- data.frame(UIJ=env$UIJ, jag_id=env$jag_id, stringsAsFactors=F)
tdat <- left_join(tdat, jtab, by='UIJ')
hdat <- left_join(hdat, jtab, by='UIJ')

#arrange temporal and env data into same order
drop.cols <- c('jag_id', 'UI','UIJ','tveg')
env <- env[order(env$jag_id),]
save(env, file=paste(mdatwd,mname,"_envdata.Rdata", Sys.Date(),sep=""))
env <- env %>% dplyr::select(-one_of(drop.cols))
env <- as.matrix(env)
tdat <- tdat[order(tdat$jag_id),]
hdat <- hdat[order(hdat$jag_id),]

#final check
if(length(unique(tdat$jag_id)) != nrow(env) | length(unique(hdat$jag_id)) != nrow(env))  print("sites not matching between spatial and temporal data!")

#save data for later analysis
save(tdat, hdat, file=paste(mdatwd,mname, Sys.Date(),"_inputdata_small.Rdata",sep=""))

###########################################################
###Prep JAGS inputs
###########################################################

nGrid=nrow(env)       ;nGrid
nBeta=ncol(env)       ;nBeta

#ages (years) at which to predict height
pred.age <- seq(0, 40, by = 1)

data=list(
  #NDVI
  age.ndvi=tdat$DA,
  ndvi=tdat$NDIN,
  id.ndvi=tdat$jag_id,
  firemonth.ndvi=tdat$firemonth,
  nObsNDVI=nrow(tdat),
  #height
  age.height=hdat$DA,
  height=hdat$Height,
  id.height=hdat$jag_id,
  nObsHeight=nrow(hdat),
  #shared
  env=env,
  nGrid=nGrid,
  nBeta=nBeta,
  pred.age=pred.age,
  nPredAge=length(pred.age)
)

#function to generate initial values
gen.inits=function(nGrid,nBeta) { list(
  alpha=runif(nGrid,0.1,0.5),
  gamma=runif(nGrid,0.1,.9),
  A=runif(nGrid,0.1,.9),
  lambda.ndvi=runif(nGrid,0.2,1),
  H0=runif(nGrid,0.1,1),
  Hmax=runif(nGrid,1,3),
  alpha.mu=runif(1,0.1,0.2),
  gamma.beta=runif(nBeta,0,1),
  gamma.tau=runif(1,1,5),
  alpha.tau=runif(1,1,5),
  lambda.ndvi.beta=runif(nBeta,0,2),
  lambda.ndvi.tau=runif(1,0.5,2),
  A.beta=runif(nBeta,0,1),
  A.tau=runif(1,1,5),
  Hmax.beta=runif(nBeta,0,1),
  Hmax.tau=runif(1,1,5),
  H0.tau=runif(1,1,5),
  tau.ndvi=runif(1,1,5),
  tau.height=runif(1,1,5)
)
}

#list of parameters to monitor (save)
params=c("phi","gamma.beta","gamma.sigma","A.beta","A.sigma","alpha","gamma","A",
         "alpha.mu","alpha.sigma","lambda.ndvi","lambda.height","lambda.ndvi.beta","lambda.ndvi.sigma",
         "H0","Hmax","H0.mu","H0.sigma","Hmax.beta","Hmax.sigma",
         "delta.mu","delta.sigma","sigma.ndvi","sigma.height","height.pred")

rm(list = ls()[-which(ls() %in% c("mname", "data", "params", "cl", "mdatwd", "gen.inits"))])
gc()

###########################################################
###Run JAGS
###########################################################

foutput=paste0(mdatwd, mname, Sys.Date(), "_modeloutput.Rdata")

m <- jags.parfit(cl = cl, #runs chains in parallel with library(dclone)
                 data = data,
                 params = params,
                 model = "Model_height.R",
                 inits = gen.inits(data$nGrid,data$nBeta),
                 n.chains = 3,
                 n.adapt=10000,n.update=10000,
                 thin = 5, n.iter = 10000
)

save(m,file=foutput)

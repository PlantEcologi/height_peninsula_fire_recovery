##################################################
## Data prep and model fitting for the joint MODIS NDVI + LiDAR height model
## (Model_height.R) at the focal plot locations only.
## Based on fit_model2026.R and fit_model_Landsat2026.R
##################################################

### Load libraries
libs=c("doParallel", "sf", "ggplot2", "tidyr", "tidyverse", "lubridate",
       "rjags", "dclone", "terra")
lapply(libs, require, character.only=T)

#file locations and names
mdatwd <- "data/"
mname <- "peninsulaHeight2026" #model name for file naming

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cl)

###########################################################
### Load point locations and height data
###########################################################

projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Points of interest, with LiDAR heights (columns named by survey year)
pts <- read.csv("data/focal_plots_n50_coords_height.csv", check.names = FALSE) |>
  mutate(geometry = str_remove_all(geometry, pattern = "[c()]")) |>
  separate_wider_delim(cols = geometry, names = c("Longitude", "Latitude"), delim = ",") |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  rename(Site = name)

height_years <- c("2015", "2017", "2021", "2023")

hdat <- pts |>
  st_drop_geometry() |>
  select(UIJ = plot, all_of(height_years)) |>
  pivot_longer(cols = all_of(height_years), names_to = "Year", values_to = "Height") |>
  mutate(Year = as.integer(Year)) |>
  filter(!is.na(Height), Height > 0)

###########################################################
### Load raster data
###########################################################

#modis ndvi from MODIS/061/MOD13Q1
NDVI <- rast("data/NDVI_stack_2001_2026.tif")
crs(NDVI) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
NAflag(NDVI) <- -9999
NDVI <- trim(NDVI)
NDVI <- project(NDVI, y = crs(projection))
NDVIdates <- read.csv("data/NDVI_dates_2001_2026.csv")
time(NDVI) <- as.Date(NDVIdates$date)
NDVI <- NDVI[[time(NDVI) <= as.Date("2023-12-31")]] # observation period (extends beyond the fire record, which ends in 2022)
names(NDVI) <- format(time(NDVI), "%Y-%m-%d")

#spatial covariates
sta <- rast("data/STATIC_stack.tif")
NAflag(sta) <- -9999
sta <- trim(sta)
crs(sta) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
sta <- project(sta, NDVI)
sta$northness <- cos(sta$aspect*pi/180)
sta$eastness <- sin(sta$aspect*pi/180)

#vegetation age
rfi <- rast("data/veldage2022.tif")
rfi <- project(rfi, NDVI)

###########################################################
### Extract raster data at the points
###########################################################

# NDVI time series
ndat <- extract(NDVI, pts, method = "simple") |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!ID) |>
  pivot_longer(cols = -UIJ, names_to = "Date", values_to = "NDVI") |>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

# Vegetation age time series
fdat <- extract(rfi, pts, method = "simple") |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!ID) |>
  pivot_longer(cols = -UIJ, names_to = "Date", values_to = "Age") |>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))

fdat <- na.omit(fdat)

# The fire age data end in 2022, so for later NDVI dates age is calculated from
# the date of the most recent fire (assumes no fires after the fire record ends)
last_fire <- fdat |>
  group_by(UIJ) |>
  slice_max(Date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(FireDate = Date - Age) |>
  select(UIJ, FireDate, LastDate = Date)

adat <- ndat |>
  select(UIJ, Date) |>
  left_join(fdat, by = c("UIJ", "Date")) |>
  left_join(last_fire, by = "UIJ") |>
  mutate(Age = if_else(is.na(Age) & Date > LastDate,
                       as.numeric(Date - FireDate), as.numeric(Age))) |>
  select(UIJ, Date, Age)

cdat <- inner_join(na.omit(ndat), na.omit(adat), by = c("UIJ", "Date"))

# Covariates
cov <- extract(sta, pts, method = "simple") |>
  as.data.frame() |>
  mutate(UIJ = pts$plot) |>
  select(!ID)

###########################################################
### Filter and trim to matching sites
###########################################################

#drop NDVI values < 0 or points with <50 NDVI observations
cdat <- cdat |>
  filter(NDVI > 0) |>
  group_by(UIJ) |>
  filter(n() >= 50) |>
  ungroup()

envars <- c("elevation", "slope", "tpi", "northness", "eastness")

cov <- cov |>
  filter(if_all(all_of(envars), ~ !is.na(.x))) |>
  filter(UIJ %in% cdat$UIJ)

# Assign the height data from each year to the NDVI date closest to mid-January
# of that year, using the vegetation age on that date
ndvi_dates <- sort(unique(ndat$Date))
height_dates <- tibble(
  Year = as.integer(height_years),
  Date = ndvi_dates[sapply(as.Date(paste0(height_years, "-01-15")),
                           function(d) which.min(abs(as.numeric(ndvi_dates - d))))]
)

hdat <- hdat |>
  inner_join(height_dates, by = "Year") |>
  inner_join(na.omit(adat), by = c("UIJ", "Date")) |>
  filter(Age > 0)

keep <- intersect(unique(cov$UIJ), unique(hdat$UIJ))
cov <- filter(cov, UIJ %in% keep)
cdat <- filter(cdat, UIJ %in% keep)
hdat <- filter(hdat, UIJ %in% keep)

rm(list = c("NDVI", "NDVIdates", "sta", "rfi", "ndat", "fdat", "adat", "last_fire")); gc()

###########################################################
### Add columns for month of fire and age in years
###########################################################

cdat <- cdat |> mutate(firemonth = month(Date - Age))
cdat$DA <- cdat$Age / 365.25

hdat <- hdat |> mutate(DA = Age / 365.25, firemonth = month(Date - Age))

###########################################################
### Select and scale environmental data
###########################################################

scaled <- scale(as.matrix(cov[, envars]))
env <- as.data.frame(cbind(intercept = 1, scaled))
env$UIJ <- cov$UIJ

beta.mu <- c(intercept = 0, attr(scaled, "scaled:center"))
beta.sd <- c(intercept = 1, attr(scaled, "scaled:scale"))
rm(scaled)

###########################################################
### Format data for JAGS
###########################################################

#create id that goes from 1 to nGrid (to order env and data in the same way)
env$jag_id <- as.integer(as.factor(env$UIJ))
jtab <- data.frame(UIJ = env$UIJ, jag_id = env$jag_id, stringsAsFactors = F)
tdat <- left_join(cdat, jtab, by = "UIJ")
hdat <- left_join(hdat, jtab, by = "UIJ")

env <- env[order(env$jag_id), ]
save(env, file = paste0(mdatwd, mname, "_envdata.Rdata", Sys.Date()))
env <- as.matrix(env[, !(names(env) %in% c("jag_id", "UIJ"))])
tdat <- tdat[order(tdat$jag_id), ]
hdat <- hdat[order(hdat$jag_id), ]

if(length(unique(tdat$jag_id)) != nrow(env))  print("sites not matching between spatial and NDVI data!")
if(length(unique(hdat$jag_id)) != nrow(env))  print("sites not matching between spatial and height data!")

save(tdat, hdat, file = paste0(mdatwd, mname, Sys.Date(), "_inputdata_small.Rdata"))

###########################################################
### Prep JAGS inputs
###########################################################

nGrid <- nrow(env)       ;nGrid
nBeta <- ncol(env)       ;nBeta

#ages (years) at which to predict height
pred.age <- seq(0, 25, by = 0.5)

data <- list(
  # NDVI
  age.ndvi = tdat$DA,
  ndvi = tdat$NDVI,
  id.ndvi = tdat$jag_id,
  firemonth.ndvi = tdat$firemonth,
  nObsNDVI = nrow(tdat),
  # Height
  age.height = hdat$DA,
  height = hdat$Height,
  id.height = hdat$jag_id,
  nObsHeight = nrow(hdat),
  # Environment
  env = env,
  nGrid = nGrid,
  nBeta = nBeta,
  # Prediction
  pred.age = pred.age,
  nPredAge = length(pred.age)
)

#function to generate initial values
gen.inits <- function(nGrid, nBeta) { list(
  alpha = runif(nGrid, 0.1, 0.5),
  gamma = runif(nGrid, 0.1, 0.9),
  A = runif(nGrid, 0.1, 0.9),
  lambda.ndvi = runif(nGrid, 0.2, 1),
  H0 = runif(nGrid, 0.1, 0.5),
  Hmax = runif(nGrid, 1, 3),
  alpha.mu = runif(1, 0.1, 0.2),
  gamma.beta = runif(nBeta, 0, 1),
  lambda.ndvi.beta = runif(nBeta, 0, 2),
  A.beta = runif(nBeta, 0, 1),
  Hmax.beta = runif(nBeta, 0, 1),
  gamma.tau = runif(1, 1, 5),
  alpha.tau = runif(1, 1, 5),
  lambda.ndvi.tau = runif(1, 1, 5),
  A.tau = runif(1, 1, 5),
  H0.tau = runif(1, 1, 5),
  Hmax.tau = runif(1, 1, 5),
  tau.ndvi = runif(1, 0, 2),
  tau.height = runif(1, 0, 2)
)
}

#parameters to monitor
params <- c("phi", "gamma.beta", "gamma.sigma", "A.beta", "A.sigma",
            "lambda.ndvi.beta", "lambda.ndvi.sigma", "Hmax.beta", "Hmax.sigma",
            "alpha", "gamma", "lambda.ndvi", "lambda.height", "A", "H0", "Hmax",
            "alpha.mu", "alpha.sigma", "H0.mu", "H0.sigma",
            "delta.mu", "delta.sigma", "sigma.ndvi", "sigma.height", "height.pred")

rm(list = ls()[-which(ls() %in% c("mname", "data", "params", "cl", "mdatwd", "gen.inits", "beta.mu", "beta.sd"))])
gc()

###########################################################
### Run JAGS
###########################################################

foutput <- paste0(mdatwd, mname, Sys.Date(), "_modeloutput.Rdata")

m <- jags.parfit(cl = cl,
                 data = data,
                 params = params,
                 model = "Model_height.R",
                 inits = gen.inits(data$nGrid, data$nBeta),
                 n.chains = 3,
                 n.adapt = 10000, n.update = 10000,
                 thin = 5, n.iter = 10000
)

save(m, file = foutput)

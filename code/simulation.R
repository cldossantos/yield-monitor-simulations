library(pacu)
library(nlme)
library(exactextractr)
library(nlraa)
library(gstat)
library(ggplot2)
library(sf)
library(terra)

## setting the number of cores
cores <- 6

## reading the data
raw.data <- st_read('./raw-data/raw-yield.shp')
trial.design <- st_read('./raw-data/trial-design.shp')
boundary <- st_read('./raw-data/boundary.shp')

## converting to utm
raw.data <- pa_2utm(raw.data)
trial.design <- pa_2utm(trial.design)
boundary <- pa_2utm(boundary)



## assigning experimental units to headland and trial
trial.design$type <- 'trial'
trial.design$type[c(21, 26, 17, 12, 7, 8, 3,
                    145, 146, 76:80,
                    147:149,
                    seq(54, 74, 5),
                    seq(45, 75, 5), 150)] <- 'headland'

## defining the nitrogen rates used in the simulation
nrates <- rep(seq(0, 300, length.out = 6), 20)

## defining a simulation grid
sim.grid <- st_make_grid(boundary, cellsize = c(10, 10))
sim.grid <- st_as_sf(sim.grid)
st_crs(sim.grid) <- st_crs(boundary)

## defining the number of simulations
nsim <- 600

## simulating a Gaussian random field
geospatial.process <- gstat(formula = ~ 1,
                            dummy = TRUE,
                            beta =  0,
                            model = vgm(0.015, 'Sph', 600))

geospatial.preds <- predict(geospatial.process,
                            newdata = st_centroid(sim.grid),
                            nsim = nsim,
                            debug.level = -1)



## building the harvest polygons
harvest.polygons <- pa_make_vehicle_polygons(points = st_geometry(raw.data),
                                             swath = raw.data$Swth_Wdth_ * 0.3048, ## converting to meter
                                             distance = raw.data$Distance_f * 0.3048, ## converting to meter
                                             angle = raw.data$Track_deg_,
                                             cores = cores,
                                             verbose = TRUE)

## adjusting harvest polygons based on overlap
tassellated.polygons <- pa_adjust_obs_effective_area(polygons = harvest.polygons,
                                                     obs.vector = rep(1, length(harvest.polygons)), ## dummy variable
                                                     var.label = 'mass',
                                                     overlap.threshold = -1,
                                                     cores = cores,
                                                     verbose = TRUE)

harvest.polygons <- st_as_sf(harvest.polygons)
harvest.polygons$id <- 1:nrow(harvest.polygons)
harvest.polygons$area <- as.numeric(st_area(harvest.polygons))

## removing invalid geometries from the tassellated polygons
## this is done to prevent errors when extracting values from
## raster later in the workflow
tassellated.polygons <- subset(tassellated.polygons, area.ratio > 0)

## checking which observations were removed
which(!(harvest.polygons$id %in% tassellated.polygons$id))


## beginning iterative simulations

i = 1
while(i <= nsim){

  cat('Iteration #', i, '\n')
  ## randomizing nitrogen rates
  new.design <- trial.design
  new.design$nrate <- 180
  new.design$nrate[new.design$type == 'trial'] <- sample(nrates, 120)

  ## defining the quadratic-plateau relationship
  sim <- st_geometry(sim.grid)
  sim <- st_as_sf(sim)
  sim$a <- 7
  sim$b <- 0.048
  sim$c <- -0.0001067
  sim$epsilon <- geospatial.preds[[i]]

  ## rasterizing the model parameters
  sim.rast <- rast(ext(sim.grid),
                   crs = st_crs(new.design)$wkt,
                   resolution = 0.25)
  a.rast <- rasterize(x = sim, y = sim.rast, field = 'a')
  b.rast <- rasterize(x = sim, y = sim.rast, field = 'b')
  c.rast <- rasterize(x = sim, y = sim.rast, field = 'c')
  nitrogen.rast <- rasterize(x = new.design, y = sim.rast, field = 'nrate')

  ## joining quadratic-plateau parameters in one raster
  qp.rast <- c(nitrogen.rast, a.rast, b.rast, c.rast)

  ## rasterizing the error
  epsilon.rast <- rasterize(x = sim, y = sim.rast, field = 'epsilon')

  ## simulating the yield response to nitrogen
  yield.rast <- app(qp.rast, function(x){quadp3(x[['nrate']], x[['a']], x[['b']], x[['c']])})

  ## adding the effect of the geospatial process to the simulated yield
  ymu <- mean(values(yield.rast), na.rm = TRUE)
  epsilon.rast <- epsilon.rast * ymu
  yi <- yield.rast + epsilon.rast

  ## simulate harvesting
  harvested.yield <- exact_extract(x = yi,
                                   y = tassellated.polygons,
                                   function(values, coverage_fraction)
                                     weighted.mean(values, coverage_fraction, na.rm = TRUE))

  harvested <- st_geometry(tassellated.polygons)
  harvested <- st_as_sf(harvested)
  harvested$id <- tassellated.polygons$id
  harvested$yield <- harvested.yield

  ## converting yield to g/m2
  harvested$yield_gm2 <- harvested$yield * 100

  ## converting yield to mass harvested within each polygon
  harvested$mass <- harvested$yield_gm2 * as.numeric(st_area(harvested))


  ## creating the synthetic data set and entering the
  ## information from the simulation
  synthetic.data <- raw.data
  synthetic.data$id <- 1:nrow(synthetic.data)
  harvested.df <- as.data.frame(harvested)
  synthetic.data <- merge(synthetic.data,
                          harvested.df[c('id', 'mass')],
                          by = 'id',
                          all.x = TRUE)
  synthetic.data$mass[is.na(synthetic.data$mass)] <- 0

  ## converting mass to mass to lbs
  synthetic.data$mass <- synthetic.data$mass / 453.592

  ## converting mass to crop flow (lbs/s)
  synthetic.data$Crop_Flw_M <- synthetic.data$mass / synthetic.data$Duration_s

  ## adding random error to the mass flow data
  random.noise <- median(synthetic.data$Crop_Flw_M) * rnorm(nrow(synthetic.data), 0, 0.05)
  synthetic.data$Crop_Flw_M <-  synthetic.data$Crop_Flw_M + random.noise
  synthetic.data$Crop_Flw_M[synthetic.data$Crop_Flw_M < 0] <- 0

  ## adjusting moisture from 15.5% to the measured in the field
  synthetic.data$Crop_Flw_M <- pacu:::.pa_moisture(synthetic.data$Crop_Flw_M,
                                                   15.5,
                                                   synthetic.data$Moisture__)

  ## calculating yield in bu/ac
  synthetic.data$Yield <- (4046 * synthetic.data$mass) / (harvest.polygons$area * 56)

  ## removing unnecessary columns from the synthetic data
  synthetic.data <- synthetic.data[-9]

  ## removing headland plots to speed up the simulations
  new.design <- subset(new.design, type == 'trial')

  ## processing the yield data using the simple algorithm
  yield.simple <- pa_yield(input = synthetic.data,
                           grid = new.design,
                           algorithm = 'simple',
                           unit.system = 'metric',
                           var.label = 'yield.simple',
                           clean = TRUE,
                           smooth.method = 'none',
                           moisture.adj = 15.5,
                           lbs.per.bushel = 56,
                           boundary = boundary)

  new.design$experimental.unit <- paste0('EU', new.design$FID)
  new.design$experimental.unit <- as.factor(new.design$experimental.unit)

  yield.ritas <- pa_yield(input = synthetic.data,
                          grid = new.design,
                          algorithm = 'ritas',
                          unit.system = 'metric',
                          data.units = c('flow' = 'lb/s'),
                          var.label = 'yield.ritas',
                          steps = TRUE,
                          smooth.method = 'krige',
                          formula = z ~ experimental.unit,
                          moisture.adj = 15.5,
                          cores = cores,
                          remove.crossed.polygons = TRUE,
                          boundary = boundary,
                          verbose = 2)



  ## joining the two processed dasta frames into one
  cmp <- st_join(yield.ritas$yield, yield.simple$yield, join = st_equals)

  ## computing the correctly specified model

  csm.yield <- exact_extract(x = yi,
                             y = cmp,
                             function(values, coverage_fraction)
                               weighted.mean(values, coverage_fraction, na.rm = TRUE))

  cmp$yield.csm <- csm.yield

  ## including the geographical coordinates in the data set so that they can be
  ## modeled
  cmp <- cbind(cmp, st_coordinates(st_centroid(cmp)))

  ## removing the headland plots
  cmp <- st_join(x = cmp, y = new.design[c('type', 'nrate')], join = st_equals)
  cmp <- subset(cmp, type == 'trial')

  ## fitting models


  fit.csm <- try(gnls(yield.csm ~ SSquadp3(x = nrate, a, b, c),
                      correlation = corSpher(form = ~ X + Y),
                      data = cmp),
                 silent = TRUE)

  fit.ritas <- try(gnls(yield.ritas ~ SSquadp3(x = nrate, a, b, c),
                        correlation = corSpher(form = ~ X + Y),
                        data = cmp),
                   silent = TRUE)

  fit.simple <- try(gnls(yield.simple ~ SSquadp3(x = nrate, a, b, c),
                         correlation = corSpher(form = ~ X + Y),
                         data = cmp),
                    silent = TRUE)


  ## checking if any models failed
  model.fail <- any(sapply(list(fit.csm, fit.ritas, fit.simple),
                           function(x){inherits(x, 'try-error')}))


  if (model.fail){
    cat('Failed to fit models\n')
    next
  }

  plot(cmp$nrate, cmp$yield.csm, pch = 16)
  lines(0:300, predict(fit.csm, newdata = list(nrate = 0:300)))
  points(cmp$nrate, cmp$yield.ritas, col = 'blue')
  lines(0:300, predict(fit.ritas, newdata = list(nrate = 0:300)), col = 'blue')
  points(cmp$nrate, cmp$yield.simple, col = 'red')
  lines(0:300, predict(fit.simple, newdata = list(nrate = 0:300)), col = 'red')
  legend('bottomright', pch = c(16, 1, 1), col = c('black', 'red', 'blue'),
         legend = c('csm', 'simple', 'ritas'))

  ## exporting the simulated data
  models <- list(fit.csm, fit.ritas, fit.simple)
  names(models) <- c('csm', 'ritas', 'simple' )
  out.name <- paste0('simulation-',
                     strftime(Sys.time(), '%Y%m%dT%H%M'),
                     '.rds')

  out.data <- list(data = list(cmp),
                   models = models)
  saveRDS(out.data, file.path('./processed-data/', out.name))


  i <- i + 1
}


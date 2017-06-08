# Run Skynet Proof of Concept -------
# Author: Dan Ovando
# Project: Skynet
# Summary: This script runs a serious of alternative models
# for estimating density from FishData using fishing behavior from
# Global Fishing Watch
rm(list = ls())
# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(randomForest)
library(lme4)
library(bigrquery)
library(stringr)
library(modelr)
library(broom)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(sf)
library(ggmap)
library(VAST)
library(TMB)
library(caret)
library(stats4)
library(extrafont)

# run options -------------------------------------------------------------

run_name <- 'testing'

run_description <- 'development of scripts'

run_dir <- file.path('results', run_name, '')

if (dir.exists(run_dir) == F) {
  dir.create(run_dir)
}

write(run_description, file = paste0(run_dir, 'description.txt'))


# run settings ------------------------------------------------------------

sci_to_use <- 'Gadus_chalcogrammus'

# figure options ----------------------------------------------------------

fig_theme <- theme_ipsum(base_size = 14)

theme_set(fig_theme)


# prepare data ------------------------------------------------------------

# download eastern bering sea trawl survey data and set VAST options
ebs_trawl <-
  FishData::download_catch_rates(survey = "EBSBTS", species_set = 50) %>%
  set_names(tolower(colnames(.)))

Version = "VAST_v2_4_0"

Method = c("Grid", "Mesh")[2]

grid_size_km = 25

n_x = c(100) #c(100, 250, 500, 1000, 2000)[1] # Number of stations

Kmeans_Config = list("randomseed" = 1,
                     "nstart" = 200,
                     "iter.max" = 1e3)

FieldConfig = c(
  "Omega1" = 1,
  "Epsilon1" = 1,
  "Omega2" = 1,
  "Epsilon2" = 1
)

RhoConfig = c(
  "Beta1" = 0,
  "Beta2" = 0,
  "Epsilon1" = 0,
  "Epsilon2" = 0
)

OverdispersionConfig = c("Vessel" = 0, "VesselYear" = 0)

ObsModel = c(2, 0)

Options =  c(
  "SD_site_density" = 0,
  "SD_site_logdensity" = 0,
  "Calculate_Range" = 1,
  "Calculate_evenness" = 0,
  "Calculate_effective_area" = 1,
  "Calculate_Cov_SE" = 0,
  'Calculate_Synchrony' = 0,
  'Calculate_Coherence' = 0
)

Region <- 'Eastern_Bering_Sea'

DateFile = paste0(run_dir, 'VAST_output/')

dir.create(DateFile)


ebs_species <- ebs_trawl %>%
  filter(sci %in% sci_to_use) %>%
  dplyr::rename(
    Year = year,
    Lat = lat,
    Lon = long,
    Catch_KG = wt,
    Sci = sci
  ) %>%
  mutate(AreaSwept_km2 = 1, Vessel = 'missing') %>%
  select(Year, Lat, Lon, Vessel, AreaSwept_km2, Catch_KG)

qmplot(
  x = Lon,
  y = Lat,
  data = ebs_species,
  color = log10(Catch_KG)
) +
  scale_color_viridis() +
  facet_wrap(~ Year) +
  theme_classic()

Record = ThorsonUtilities::bundlelist(
  c(
    "ebs_species",
    "Version",
    "Method",
    "grid_size_km",
    "n_x",
    "FieldConfig",
    "RhoConfig",
    "OverdispersionConfig",
    "ObsModel",
    "Kmeans_Config"
  )
)
save(Record, file = paste0(DateFile, "Record.RData"))
capture.output(Record, file = paste0(DateFile, "Record.txt"))
strata.limits <- data.frame('STRATA' = "All_areas")

Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region = Region, strata.limits = strata.limits)

ebs_species <- ebs_species %>%
  na.omit() %>%
  mutate(Vessel = as.factor(Vessel))

Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(
  grid_size_km = grid_size_km,
  n_x = n_x,
  Method = Method,
  Lon = ebs_species[, 'Lon'],
  Lat = ebs_species[, 'Lat'],
  Extrapolation_List = Extrapolation_List,
  randomseed = Kmeans_Config[["randomseed"]],
  nstart = Kmeans_Config[["nstart"]],
  iter.max = Kmeans_Config[["iter.max"]],
  DirPath = DateFile,
  Save_Results = FALSE
)
# Add knots to Data_Geostat
ebs_species = cbind(ebs_species, "knot_i" = Spatial_List$knot_i)

TmbData = Data_Fn(
  "Version" = Version,
  "FieldConfig" = FieldConfig,
  "OverdispersionConfig" = OverdispersionConfig,
  "RhoConfig" = RhoConfig,
  "ObsModel" = ObsModel,
  "c_i" = rep(0, nrow(ebs_species)),
  "b_i" = ebs_species[, 'Catch_KG'],
  "a_i" = ebs_species[, 'AreaSwept_km2'],
  "v_i" = as.numeric(ebs_species[, 'Vessel']) - 1,
  "s_i" = ebs_species[, 'knot_i'] - 1,
  "t_iz" = ebs_species[, 'Year'],
  "a_xl" = Spatial_List$a_xl,
  "MeshList" = Spatial_List$MeshList,
  "GridList" = Spatial_List$GridList,
  "Method" = Spatial_List$Method,
  "Options" = Options
)

TmbList = Build_TMB_Fn(
  "TmbData" = TmbData,
  "RunDir" = DateFile,
  "Version" = Version,
  "RhoConfig" = RhoConfig,
  "loc_x" = Spatial_List$loc_x,
  "Method" = Method
)
Obj = TmbList[["Obj"]]

Opt = TMBhelper::Optimize(
  obj = Obj,
  lower = TmbList[["Lower"]],
  upper = TmbList[["Upper"]],
  getsd = TRUE,
  savedir = DateFile,
  bias.correct = FALSE
)

Report = Obj$report()
Save = list(
  "Opt" = Opt,
  "Report" = Report,
  "ParHat" = Obj$env$parList(Opt$par),
  "TmbData" = TmbData
)
save(Save, file = paste0(DateFile, "Save.RData"))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn(
  "Region" = Region,
  "NN_Extrap" = Spatial_List$PolygonList$NN_Extrap,
  "Extrapolation_List" = Extrapolation_List
)
# Decide which years to plot
Year_Set = seq(min(ebs_species[, 'Year']), max(ebs_species[, 'Year']))
Years2Include = which(Year_Set %in% sort(unique(ebs_species[, 'Year'])))

Dens_xt <-
  SpatialDeltaGLMM::PlotResultsOnMap_Fn(
    plot_set = c(3),
    MappingDetails = MapDetails_List[["MappingDetails"]],
    Report = Report,
    Sdreport = Opt$SD,
    PlotDF = MapDetails_List[["PlotDF"]],
    MapSizeRatio = MapDetails_List[["MapSizeRatio"]],
    Xlim = MapDetails_List[["Xlim"]],
    Ylim = MapDetails_List[["Ylim"]],
    FileName = ,
    ,
    ,
    Year_Set = Year_Set,
    Years2Include = Years2Include,
    Rotate = MapDetails_List[["Rotate"]],
    Cex = MapDetails_List[["Cex"]],
    Legend = MapDetails_List[["Legend"]],
    zone = MapDetails_List[["Zone"]],
    mar = c(0, 0, 2, 0),
    oma = c(3.5, 3.5, 0, 0),
    cex = 1.8,
    plot_legend_fig = FALSE
  )

# note this only works for one species... need to work on this.
#
ebs_species_densities = cbind(
  "density" = as.vector(Dens_xt),
  "year" = Year_Set[col(Dens_xt)],
  "e_km" = Spatial_List$MeshList$loc_x[row(Dens_xt), 'E_km'],
  "n_km" = Spatial_List$MeshList$loc_x[row(Dens_xt), 'N_km']
) %>%
  as_data_frame() %>%
  mutate(knot = as.numeric(factor(paste(e_km, n_km)))) %>%
  arrange(knot, year)

project <- "ucsb-gfw"

fishing_connection <-
  src_bigquery(project, "skynet") # This function initiliazes a connection with a BQ dataset


ebs_raw <- fishing_connection %>%
  tbl("ebs_w_vessels") %>%
  collect(n = Inf)

ebs <- ebs_raw %>%
  select(-b_mmsi,-b_year) %>%
  set_names(str_replace_all(colnames(.), '(a_)|(b_)', '')) %>%
  filter(is.na(on_fishing_list_nn) == F) %>%
  arrange(year, mmsi)

qmplot(
  x = rounded_lon,
  y = rounded_lat,
  data = ebs,
  color = log10(total_hours)
) +
  scale_color_viridis() +
  facet_wrap(~ year) +
  theme_classic()

utm_coords <-
  SpatialDeltaGLMM::Convert_LL_to_UTM_Fn(
    Lon = ebs$rounded_lon,
    Lat = ebs$rounded_lat,
    zone = Extrapolation_List$zone,
    flip_around_dateline = Extrapolation_List$flip_around_dateline
  ) %>%
  dplyr::rename(e_km = X, n_km = Y) %>%
  select(e_km, n_km)

ebs <- ebs %>%
  bind_cols(utm_coords)

ebs %>%
  ggplot(aes(best_label)) +
  geom_bar()

knots <- ebs_species_densities %>%
  select(knot, e_km, n_km) %>%
  unique()

nearest_knot <-
  RANN::nn2(knots %>% select(-knot), utm_coords, k = 1)

ebs$knot <- knots$knot[nearest_knot$nn.idx]

ebs$distance <- nearest_knot$nn.dists


ebs %>%
  filter(distance <= quantile(ebs$distance, 0.75)) %>%
  select(e_km, n_km, knot) %>%
  unique() %>%
  ggplot(aes(e_km, n_km, label = knot)) +
  geom_text() +
  geom_text(data = knots, aes(e_km, n_km, label = knot, color = 'red'))


trawl_fishing_by_knot <- ebs %>%
  filter(best_label == 'trawlers',
         distance <= quantile(ebs$distance, 0.75)) %>%
  dplyr::rename(th = total_hours) %>%
  group_by(year, knot) %>%
  dplyr::summarise(
    total_hours = sum(th, na.rm = T),
    total_engine_hours = sum(th * inferred_engine_power, na.rm = T),
    num_vessels = length(unique(mmsi)),
    dist_from_shore = mean(mean_distance_from_shore, na.rm = T),
    dist_from_port = mean(mean_distance_from_port, na.rm = T),
    mean_vessel_length = mean(inferred_length, na.rm = T)
  )

skynet_data <- trawl_fishing_by_knot %>%
  left_join(ebs_species_densities, by = c('year', 'knot')) %>%
  ungroup() %>%
  mutate(vessel_hours = total_engine_hours * num_vessels) %>%
  group_by(knot) %>%
  arrange(knot, year) %>%
  mutate(
    total_engine_hours_lag1 = lag(total_engine_hours, 1),
    total_engine_hours_lag2 = lag(total_engine_hours, 2),
    total_engine_hours_lag3 = lag(total_engine_hours, 3)
  ) %>%
  ungroup() %>%
  filter(total_engine_hours > 0,
         density > 0) #fix this later, need to include zeros


# prepare training data ---------------------------------------------------
# expand this more later, for now let's focus on getting all the models up and running well
skynet_sandbox <- modelr::crossv_kfold(skynet_data, 5)


# prep maching learning ---------------------------------------------------

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

gbm_model <- train(
  density ~ total_engine_hours + dist_from_port + dist_from_shore + num_vessels + mean_vessel_length,
  data = skynet_sandbox$train[[1]] %>% as_data_frame(),
  method = "gbm",
  trControl = fitControl
)



# fit simple model --------------------------------------------------------


lm_model <-
  lm(
    density ~ total_engine_hours + dist_from_port + dist_from_shore + num_vessels + mean_vessel_length,
    data = skynet_sandbox$train[[1]] %>% as_data_frame()
  )



# fit structural model ----------------------------------------------------

fit_sffs <- function(beta_distance,
                     beta_size,
                     beta_interaction,
                     q,
                     sigma,
                     dat,
                     price,
                     marginal_profits = 0,
                     use = 'mle') {
  old_sigma <- sigma
  sigma <-  exp(sigma)

  old_q <- q

  q <- exp(q)

  beta_distance <-  exp(beta_distance)

  beta_size <- exp(beta_size)

  beta_interaction <- exp(beta_interaction)

  # beta_intercept <- exp(beta_intercept)

  cost <-
    beta_distance * dat$dist_from_port + beta_size * dat$mean_vessel_length + beta_interaction * dat$dist_from_port * dat$mean_vessel_length

  estimated_abundance <-
    ((cost + marginal_profits)) / (price  * exp(-q * dat$total_engine_hours))
  if (use == 1) {
    out <-
      -sum(stats::dlnorm(dat$density, log(estimated_abundance), sigma, log = T))
  } else {
    out <- estimated_abundance
  }
  print(out)
  return(out)

}



mle_model <- stats4::mle(
  fit_sffs,
  start = list(
    beta_distance = log(10),
    beta_size = 10 %>% log,
    beta_interaction = 10 %>% log,
    sigma = log(sd(skynet_data$density) / 10),
    q = log(1e-7)
  ),
  fixed = list(
    price = 765,
    dat = skynet_data %>% filter(total_engine_hours > 0),
    marginal_profits = 0,
    use = 1
  ),
  upper = list(sigma = log(100),
               q = log(1.5e-5))
)


mle_coefs <- coef(mle_model)


#  assess performance -----------------------------------------------------

performance <- skynet_sandbox$test[[1]] %>% as_data_frame()

performance <- performance %>%
  mutate(
    gbm = predict(gbm_model, newdata = .),
    lm = predict(lm_model, newdata = .),
    structural = fit_sffs(
      beta_distance = mle_coefs['beta_distance'],
      beta_size = mle_coefs['beta_size'],
      beta_interaction = mle_coefs['beta_interaction'],
      q = mle_coefs['q'],
      sigma = mle_coefs['sigma'],
      dat = .,
      price = 765,
      marginal_profits = 10,
      use = 0
    )
  ) %>%
  gather('model', 'density_hat', gbm:structural)


performance <- performance %>%
  mutate(se = (density - density_hat) ^ 2)

performance %>%
  group_by(model) %>%
  dplyr::summarise(rmse = sqrt(mean(se)))

performance_plot <- performance %>%
  ggplot(aes(density, density_hat)) +
  geom_abline(aes(intercept = 0, slope = 1),
              color = 'red',
              linetype = 2) +
  geom_point(alpha = 0.75) +
  facet_wrap( ~ model) +
  labs(x = 'EBS Trawl Survey Alaska Pollock Density',
       y = 'Model Predicted Density',
       caption = 'Models trained on separate data')

ggsave(filename = paste0(run_dir,'performance_plot.pdf'), performance_plot)
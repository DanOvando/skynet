#' Apply VAST to survey data
#'
#' \code{vasterize_index} takes survey data and applies the VAST model
#' to provide standardized spatio-temporal abundance indicies
#'
#' @param raw_data the survey data in question
#' @param region the general region, see VAST documentation
#' @param run_dir the directory results will be stored in
#' @param version the vast version to use
#' @param method grid or mesh method
#' @param strata
#' @param n_x
#' @param randomseed
#' @param nstart
#' @param iter.max
#' @param omega1
#' @param epsilon1
#' @param omega2
#' @param epsilon2
#' @param beta1
#' @param beta2
#' @param vessel
#' @param vessel_year
#' @param obs_model a vector with the statistical model for the "delta" model
#' @param sd_site_density
#' @param sd_site_logdensity
#' @param calculate_range
#' @param calculate_evenness
#' @param calculate_effective_area
#' @param calculate_cov_se
#' @param calculate_synchrony
#' @param calculate_coherence
#'
#' @return a list of cleaned up outputs from VAST
#' @export
#'
vasterize_index <- function(raw_data,
                            region = 'California_current',
                            run_dir = '../results/',
                            version = "VAST_v2_4_0",
                            method = 'Mesh',
                            strata = "All_areas",
                            n_x = 1000,
                            spp,
                            grid_size_km = 25,
                            randomseed = 1,
                            nstart = 100,
                            iter.max = 1e3,
                            omega1 = 1,
                            epsilon1 = 1,
                            omega2 = 1,
                            epsilon2 = 1,
                            beta1 = 0,
                            beta2 = 0,
                            vessel = 0,
                            vessel_year = 0,
                            obs_model = c(2, 0),
                            sd_site_density = 0,
                            sd_site_logdensity = 0,
                            calculate_range = 1,
                            calculate_evenness = 0,
                            calculate_effective_area = 1,
                            calculate_cov_se = 0,
                            calculate_synchrony = 0,
                            calculate_coherence = 0) {
  raw_data <- as.data.frame(raw_data) %>%
    mutate(vessel = as.factor(vessel),
           spp = as.factor(spp))

  years <- min(raw_data$year):max(raw_data$year)

  Kmeans_Config = list("randomseed" = randomseed,
                       "nstart" = nstart,
                       #the number of times that the k-means algorithm is run while searching for the best solution (default=100)
                       "iter.max" = iter.max)

  # need to get a better understanding of what these do
  FieldConfig = c(
    "Omega1" = omega1,
    "Epsilon1" = epsilon1,
    "Omega2" = omega2,
    "Epsilon2" = epsilon2
  )

  RhoConfig = c(
    "Beta1" = 0,
    "Beta2" = 0,
    "Epsilon1" = 0,
    "Epsilon2" = 0
  )
  # overdispersion_config = c("vessel" = vessel, "vessel_year" = vessel_year)

  overdispersion_config = c("Eta1"=0, "Eta2"=0)

  options =  c(
    "SD_site_density" = sd_site_density,
    "SD_site_logdensity" = sd_site_logdensity,
    "Calculate_Range" = calculate_range,
    "Calculate_Evenness" = calculate_evenness,
    "Calculate_effective_area" = calculate_effective_area,
    "Calculate_Cov_SE" = calculate_cov_se,
    'Calculate_Synchrony' = calculate_synchrony,
    'Calculate_Coherence' = calculate_coherence
  )

  vast_file <-  paste0(getwd(), '/', run_dir, 'VAST_output')

  if (dir.exists(vast_file) == F) {
    dir.create(vast_file)
  }

  strata_limits <- data.frame('STRATA' = strata)

  extrapolation_list <-
    SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region = region, strata.limits = strata_limits)

  spatial_list <-  SpatialDeltaGLMM::Spatial_Information_Fn(
    n_x = n_x,
    grid_size_km = grid_size_km,
    Method = method,
    Lon = raw_data[, 'lon'],
    Lat = raw_data[, 'lat'],
    Extrapolation_List = extrapolation_list,
    randomseed = Kmeans_Config[["randomseed"]],
    nstart = Kmeans_Config[["nstart"]],
    iter.max = Kmeans_Config[["iter.max"]],
    DirPath = vast_file,
    Save_Results = FALSE
  )


  raw_data <-  cbind(raw_data, "knot_i" = spatial_list$knot_i)

  tmb_data <-  VAST::Data_Fn(
    "Version" = version,
    "FieldConfig" = FieldConfig,
    "OverdispersionConfig" = overdispersion_config,
    "RhoConfig" = RhoConfig,
    "ObsModel_ez" = obs_model,
    "c_iz" = as.numeric(raw_data[, 'spp']) - 1,
    "b_i" = raw_data[, 'catch_kg'],
    "a_i" = raw_data[, 'areaswept_km2'],
    "v_i" = as.numeric(raw_data[, 'vessel']) - 1,
    "s_i" = raw_data[, 'knot_i'] - 1,
    "t_iz" = raw_data[, 'year'],
    "a_xl" = spatial_list$a_xl,
    "MeshList" = spatial_list$MeshList,
    "GridList" = spatial_list$GridList,
    "Method" = spatial_list$Method,
    "Options" = options
  )
  tmb_list <-  VAST::Build_TMB_Fn(
    "TmbData" = tmb_data,
    "RunDir" = vast_file,
    "Version" = version,
    "RhoConfig" = RhoConfig,
    "loc_x" = spatial_list$loc_x,
    "Method" = method
  )

  obj <-  tmb_list[["Obj"]]

  opt <-  TMBhelper::Optimize(
    obj = obj,
    lower = tmb_list[["Lower"]],
    upper = tmb_list[["Upper"]],
    getsd = TRUE,
    savedir = vast_file,
    bias.correct = FALSE
  )

  # dial in TMB fit
  opt <-  TMBhelper::Optimize(
    obj = obj,
    start = obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par), 0, 0.2),
    lower = tmb_list[["Lower"]],
    upper = tmb_list[["Upper"]],
    getsd = TRUE,
    savedir = vast_file,
    bias.correct = FALSE,
    newtonsteps = 3
  )

  if(any(abs(opt$diagnostics$final_gradient) > 0.0001)){
    warning("VAST probably didn't converge after tuning")
  }

  report <-  obj$report()

  map_details_list <-  SpatialDeltaGLMM::MapDetails_Fn(
    "Region" = region,
    "NN_Extrap" = spatial_list$PolygonList$NN_Extrap,
    "Extrapolation_List" = extrapolation_list
  )
  # Decide which years to plot

  spatial_densities <- report$D_xcy %>%
    reshape2::melt() %>%
    as_data_frame() %>%
    set_names(c('knot', 'species', 'year', 'density')) %>%
    mutate(
      species = factor(species, labels = unique(raw_data$spp) %>% as.character()) %>% as.character(),
      year = factor(year, labels = years) %>% as.character() %>% as.numeric()
    ) %>%
    left_join(spatial_list$loc_x_lat_long, by = 'knot')


  spatial_biomass <-  report$Index_xcyl %>%
    drop() %>%
    reshape2::melt() %>%
    as_data_frame() %>%
    set_names(c('knot', 'year', 'biomass')) %>%
    mutate(
      species = factor(spp, labels = spp %>% as.character()) %>% as.character(),
      year = factor(year, labels = years) %>% as.character() %>% as.numeric()
    ) %>%
    left_join(spatial_list$loc_x_lat_long, by = 'knot')


  spatial_densities <- spatial_densities %>%
    left_join(spatial_biomass %>% select(knot, species, year, biomass),
              by = c("knot","species","year")) %>%
    mutate(area = biomass / density)

  year_set <- years

  years_2_include <-
    which(year_set %in% years)

  diagnostics <-
    opt$diagnostics[, c('Param', 'Lower', 'MLE', 'Upper', 'final_gradient')]

  # get overall index

  index <-
    (
      SpatialDeltaGLMM::PlotIndex_Fn(
        DirName = vast_file,
        TmbData = tmb_data,
        Sdreport = opt[["SD"]],
        Year_Set = year_set,
        Years2Include = years_2_include,
        use_biascorr = TRUE
      )
    )


  vast_index <- index$Table %>%
    select(Year, Unit, Estimate_metric_tons) %>%
    rename(abundance = Estimate_metric_tons) %>%
    mutate(source = 'vast')


  outlist <-
    list(
      extrapolation_list = extrapolation_list,
      spatial_list = spatial_list,
      spatial_densities = spatial_densities,
      time_index = vast_index,
      report = report,
      obj = obj,
      opt = opt,
      diagnostics = diagnostics
    )

  return(outlist)


}
vasterize_index <- function(raw_data,
                            region = 'California_current',
                            run_dir = '../results/',
                            version = "VAST_v2_4_0",
                            method = 'Mesh',
                            strata = "All_areas",
                            n_x = 1000,
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
  Kmeans_Config = list("randomseed" = randomseed,
                       "nstart" = nstart,
                       #the number of times that the k-means algorithm is run while searching for the best solution (default=100)
                       "iter.max" = iter.max)

  # Really no idea what almost any of these actually do
  FieldConfig = c(
    "Omega1" = omega1,
    "Epsilon1" = epsilon1,
    "Omega2" = omega2,
    "Epsilon2" = epsilon2
  )

  RhoConfig = c(
    "Beta1" = beta1,
    "Beta2" = beta2,
    "Epsilon1" = epsilon1,
    "Epsilon2" = epsilon2
  )

  overdispersion_config = c("vessel" = vessel, "vessel_year" = vessel_year)

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


  vast_file <-  paste0(run_dir, 'VAST_output/')

  if (dir.exists(vast_file) == F) {
    dir.create(vast_file)
  }

  strata_limits <- data.frame('STRATA' = strata)

  extrapolation_list <-
    SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region = region, strata.limits = strata_limits)


  spatial_list <-  SpatialDeltaGLMM::Spatial_Information_Fn(
    n_x = n_x,
    Method = Method,
    Lon = raw_data[, 'Lon'],
    Lat = raw_data[, 'Lat'],
    Extrapolation_List = extrapolation_list,
    randomseed = Kmeans_Config[["randomseed"]],
    nstart = Kmeans_Config[["nstart"]],
    iter.max = Kmeans_Config[["iter.max"]],
    DirPath = vast_file,
    Save_Results = FALSE
  )


  raw_data <-  cbind(raw_data, "knot_i" = Spatial_List$knot_i)


  tmb_data <-  Data_Fn(
    "Version" = Version,
    "FieldConfig" = FieldConfig,
    "OverdispersionConfig" = overdispersion_config,
    "RhoConfig" = RhoConfig,
    "ObsModel" = obs_model,
    "c_i" = as.numeric(raw_data[, 'spp']) - 1,
    "b_i" = raw_data[, 'Catch_KG'],
    "a_i" = raw_data[, 'AreaSwept_km2'],
    "v_i" = as.numeric(raw_data[, 'Vessel']) - 1,
    "s_i" = raw_data[, 'knot_i'] - 1,
    "t_iz" = raw_data[, 'Year'],
    "a_xl" = spatial_list$a_xl,
    "MeshList" = spatial_list$MeshList,
    "GridList" = spatial_list$GridList,
    "Method" = spatial_list$Method,
    "Options" = options
  )
  tmb_list <-  Build_TMB_Fn(
    "TmbData" = tmb_data,
    "RunDir" = vast_file,
    "Version" = Version,
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

  report <-  obj$report()

  map_details_list <-  SpatialDeltaGLMM::MapDetails_Fn(
    "Region" = region,
    "NN_Extrap" = spatial_list$PolygonList$NN_Extrap,
    "Extrapolation_List" = extrapolation_list
  )
  # Decide which years to plot
  year_set <-
    seq(min(raw_data[, 'Year']), max(raw_data[, 'Year']))

  years_2_include <-
    which(year_set %in% sort(unique(raw_data[, 'Year'])))

  #god damn it Dens_xt is in log space
  dens_xt <-
    SpatialDeltaGLMM::PlotResultsOnMap_Fn(
      plot_set = c(3),
      MappingDetails = map_details_list[["MappingDetails"]],
      Report = report,
      Sdreport = opt$SD,
      PlotDF = map_details_list[["PlotDF"]],
      MapSizeRatio = map_details_list[["MapSizeRatio"]],
      Xlim = map_details_list[["Xlim"]],
      Ylim = map_details_list[["Ylim"]],
      FileName = ,
      ,
      ,
      Year_Set = year_set,
      Years2Include = years_2_include,
      Rotate = map_details_list[["Rotate"]],
      Cex = map_details_list[["Cex"]],
      Legend = map_details_list[["Legend"]],
      zone = map_details_list[["Zone"]],
      mar = c(0, 0, 2, 0),
      oma = c(3.5, 3.5, 0, 0),
      cex = 1.8,
      plot_legend_fig = FALSE
    )

  diagnostics <-
    opt$diagnostics[, c('Param', 'Lower', 'MLE', 'Upper', 'final_gradient')]

  # get overall index
  #
  index <-
    SpatialDeltaGLMM::PlotIndex_Fn(
      DirName = vast_file,
      TmbData = tmb_data,
      Sdreport = opt[["SD"]],
      Year_Set = year_set,
      Years2Include = years_2_include,
      use_biascorr = TRUE
    )


  vast_index <- index$Table %>%
    select(Year, Estimate_metric_tons) %>%
    rename(abundance = Estimate_metric_tons) %>%
    mutate(source = 'vast')

  outlist <-
    list(
      spatial_index = dens_xt,
      time_index = vast_index,
      report = report,
      obj = obj,
      opt = opt
    )

  return(outlist)


}
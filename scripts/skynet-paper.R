library(modelr)
library(broom)
library(hrbrthemes)
library(viridis)
library(sf)
library(ggmap)
library(caret)
library(extrafont)
library(patchwork)
library(wesanderson)
library(ggridges)
library(taxize)
library(tidyposterior)
library(recipes)
library(rstan)
library(brms)
library(tidyverse)
extrafont::loadfonts()

# Sys.setenv(R_MAX_VSIZE = 9e9)

functions <- list.files(here::here("functions"))

walk(functions, ~ here::here("functions", .x) %>% source()) # load local functions

# filters -----------------------------------------------------------------


bad_surveys = "wcghl"

test_sets_of_interest <- c("random",
                           "historic",
                           "spatial",
                           "wa-or")

# run options -------------------------------------------------------------

run_name <- 'd1.1'

run_dir <- here::here('results',run_name)

fig_theme <- theme_ipsum(base_size = 14, axis_title_size = 18)

theme_set(fig_theme)

load(here::here("results", run_name,'skynet_data.Rdata'))

load(here::here("results", run_name,'skynet_models.Rdata'))

load(here::here("data","vast_fish.Rdata"))

load(file = here::here("data","global_map.Rdata"))

load(file = here::here('data','fish_prices.Rdata'))

load(file = here::here('data','gfw_data.Rdata'))

map_size <- 1


skynet_data <- candidate_data %>%
  filter(fished_only == T, survey_months_only ==T) %>%
  select(skynet_data) %>%
  unnest()

# prep things

pacific_map <- global_map %>%
  as("Spatial") %>%
  maptools::nowrapRecenter() %>%
  sf::st_as_sf()

map_theme <-   theme_get() +theme(legend.key.height = unit(1.5,"cm"),
                                  axis.text.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  plot.margin = unit(c(0,0,0,0),"cm"),
                                  legend.box.margin =  ggplot2::margin(0,0,0,0))

vast_fish <- vast_fish[!map_lgl(vast_fish$vasterized_data, is.null), ]


skynet_models <- skynet_models %>%
  mutate(index = 1:nrow(.)) %>%
  mutate(error = map(fitted_model, "error")) %>%
  mutate(no_error = map_lgl(error, is.null)) %>%
  filter(no_error) #%>%
  # filter(test_sets %in% c("random",
  #                         "spatial",
  #                         "california"))

diagnostic_plot_foo <- function(data,
                                r2,
                                dep_var,
                                test_region,
                                train_region,
                                data_set) {
  data %>%
    filter(surveyed_year == T) %>%
    ggplot(aes_(as.name(dep_var), ~ pred)) +
    geom_abline(aes(intercept = 0, slope = 1),
                color = "red",
                linetype = 2) +
    geom_point(alpha = 0.75) +
    labs(
      title = paste0(
        "tested on ",
        test_region,
        "- trained on ",
        train_region,
        "; R2 = ",
        r2
      ),
      subtitle = paste0("data subset is: ", data_set)
    )
} # close diagnostic functions

trim_skynet <- function(x, object,vars = c("survey","surveyed_year","year", "rounded_lat","rounded_lon","knot")){

  y <- x$result %>% pluck(object)

  out <- y %>%
    dplyr::select(-contains("interaction"),-"random_var",-contains("log"))
  # dplyr::select(vars, contains("density"),"pred")

}


skynet_models <- skynet_models %>%
  mutate(
    test_data = map(fitted_model, trim_skynet, object = "test_predictions"),
    training_data = map(fitted_model,trim_skynet, "training_predictions")
  ) %>%
  select(-fitted_model) %>%
  mutate(
    var_y = map2_dbl(test_data, dep_var, ~ var(.x[, .y])),
    r2 = map2_dbl(test_data, dep_var, ~ yardstick::rsq(.x %>% filter(surveyed_year == T), .y, "pred")),
    rmse = map2_dbl(test_data, dep_var, ~ yardstick::rmse(.x %>% filter(surveyed_year == T), .y, "pred")),
    ccc = map2_dbl(test_data, dep_var, ~ yardstick::ccc(.x %>% filter(surveyed_year == T), .y, "pred")),
    r2_training = map2_dbl(training_data, dep_var, ~ yardstick::rsq(.x %>% filter(surveyed_year == T), .y, "pred")),
    rmse_training = map2_dbl(training_data, dep_var, ~ yardstick::rmse(.x %>% filter(surveyed_year == T), .y, "pred")),
    ccc_training = map2_dbl(training_data, dep_var, ~ yardstick::ccc(.x %>% filter(surveyed_year == T), .y, "pred"))
  ) #%>%
  # mutate(
  #   test_plot = pmap(
  #     list(
  #       data = test_data,
  #       r2 = r2,
  #       test_region = test_sets,
  #       train_region = train_set,
  #       data_set = data_subset,
  #       dep_var = dep_var
  #     ),
  #     diagnostic_plot_foo
  #   ),
  #   training_plot = pmap(
  #     list(
  #       data = training_data,
  #       r2 = r2_training,
  #       test_region = paste0(test_sets),
  #       train_region = paste0(train_set, "- training plot"),
  #       data_set = data_subset,
  #       dep_var = dep_var
  #     ),
  #     diagnostic_plot_foo
  #   )
  # )

skynet_models <- skynet_models %>%
    mutate(unfished_only = str_detect(data_subset, "unfished"))


# data summaries ----------------------------------------------------------

fish_summary <- vast_fish %>%
  filter(!survey %in% bad_surveys) %>%
  select(survey, data,spp) %>%
  unnest() %>%
  group_by(survey, spp) %>%
  summarise(n_seen = sum(catch_kg > 0)) %>%
  group_by(spp) %>%
  mutate(tc = sum(n_seen)) %>%
  ungroup() %>%
  mutate(spp = forcats::fct_reorder(spp, tc)) %>%
  group_by(survey) %>%
  top_n(10, tc) %>%
  ungroup()

fish_summary %>%
  ggplot(aes(spp, tc, fill = survey)) +
  geom_col() +
  coord_flip() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank()) +
  labs(y = "# of times observed") +
  scale_fill_viridis_d()

# how does effort relate to abundance -------------------------------------


downscaled_performance <-
  cross_df(list(
    resolution =  seq(25, 200, by = 25),
    training_data = list(skynet_data %>% filter(no_take_mpa == FALSE, restricted_use_mpa == FALSE))
  ))

downscaled_performance <-  downscaled_performance %>%
  mutate(
    new_grid  = map2(
      training_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  )

downscaled_performance <-  downscaled_performance %>%
  mutate(
    new_data = map2(
      training_data,
      new_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat
    )
  )


e_v_d_plot <- downscaled_performance %>%
  select(resolution, new_data) %>%
  unnest() %>%
  filter(agg_total_hours > 0) %>%
  ggplot(aes((agg_total_hours), (agg_mean_density), fill = approx_survey)) +
  geom_point(alpha = 0.5, shape = 21,size = 2) +
  labs(x = "Fishing Hours", y = bquote("Fish Density"~(ton/km^2)),
       caption = "Note log-10 scale of axes") +
  scale_fill_viridis_d(name = "Survey Region") +
  scale_color_viridis_d(name = "Survey Region") +
  geom_abline(aes(intercept = 0, slope = 1),linetype = 2, alpha = 0.75) +
  geom_smooth(method = "lm", aes(color = approx_survey), se = FALSE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  facet_wrap(~resolution)


total_e_v_d <- downscaled_performance %>%
  select(resolution, new_data) %>%
  unnest() %>%
  filter(agg_total_hours > 0,
         resolution == 200) %>%
  group_by(approx_survey) %>%
  summarise(available_rev = sum(agg_mean_economic_density * resolution),
            total_hours = sum(agg_total_hours)) %>%
  mutate(approx_survey = fct_reorder(approx_survey, total_hours)) %>%
  ggplot(aes(approx_survey,available_rev)) +
  geom_point()


# hack in here to deal with wrong price units in analysis
e_v_ed_plot <- downscaled_performance %>%
  select(resolution, new_data) %>%
  unnest() %>%
  filter(agg_total_hours > 0,
         resolution %in% c(100)) %>%
  ggplot(aes((agg_total_hours), (total_revenue * 1000), fill = approx_survey)) +
  geom_point(alpha = 0.5, shape = 21,size = 2) +
  labs(x = "Fishing Hours", y = "Available Revenues (USD)",
       caption = "Note log-10 scale of axes") +
  scale_fill_viridis_d(name = "Survey Region") +
  scale_color_viridis_d(name = "Survey Region") +
  geom_abline(aes(intercept = 0, slope = 1),linetype = 2, alpha = 0.75) +
  geom_smooth(method = "lm", aes(color = approx_survey), se = FALSE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::dollar) +
  facet_wrap(~resolution)

e_v_ed_plot


# catch per unit effort ---------------------------------------------------

comma_foo <- function(x){

  if (str_detect(x, ", ")){

    split <- str_split(x, ", ", simplify = T)

    x <- glue::glue("{split[2]} {split[1]}")

  }

  return(x)
}

nmfs_landings <- read_csv("~/Box Sync/Databases/NMFS/nmfs_pacific_state_landings.csv") %>%
  mutate(common_name = map_chr(common_name, comma_foo))

lookup_nmfs_names <-  TRUE
if (lookup_nmfs_names == TRUE){
nmfs_sci_names <- data_frame(common_name = unique(nmfs_landings$common_name)) %>%
  mutate(scientific_name = map(common_name,~taxize::comm2sci(.x, db = "worms")))

nmfs_sci_names <- nmfs_sci_names %>%
  mutate(scientific_name = map_chr(scientific_name, ~ifelse(length(.x[[1]])>0,.x[[1]], NA)))

saveRDS(nmfs_sci_names, here::here("processed_data","nmfs_names.RDS"))
} else{

  nmfs_sci_names <- readRDS(here::here("processed_data","nmfs_names.RDS"))

}

nmfs_landings <- nmfs_landings %>%
  left_join(nmfs_sci_names, by = "common_name")

load(here::here("data","DBdata.RData")) # as of copy, RAM v4.40 with model fits

fao_capture <- read_csv("~/Box Sync/Databases/FAO/fao_2018_capture_statistics.csv")

included_species <- vast_fish %>%
  filter(survey %in% c("ebsbts"))

included_species <- unique(included_species$spp) %>%
  str_replace_all("_"," ")

ram_biomass <- tb.data %>%
  mutate(year = as.numeric(rownames(.))) %>%
  select(year, everything()) %>%
  gather(stockid, biomass, -year) %>%
  left_join(stock, by = "stockid") %>%
  left_join(assessment, by = "stockid")


ram_catch <- tcbest.data %>%
  mutate(year = as.numeric(rownames(.))) %>%
  select(year, everything()) %>%
  gather(stockid, catch, -year) %>%
  left_join(stock, by = "stockid") %>%
  left_join(assessment, by = "stockid")


nmfs_catch <- nmfs_landings %>%
  filter(state == "Alaska", scientific_name %in% included_species) %>%
  group_by(year) %>%
  summarise(total_nmfs_catch = sum(metric_tons, na.rm = T))

fao_catch <- fao_capture %>%
  filter(fao_area_code == 67, scientific_name %in% included_species) %>%
  group_by(year) %>%
  summarise(total_fao_catch = sum(capture, na.rm = T))

ebsbts_catch <- ram_catch %>%
  filter(region == "US Alaska") %>%
  filter(!str_detect(stocklong.x, "Gulf of Alaska")) %>%
  filter(scientificname %in% included_species) %>%
  group_by(year) %>%
  summarise(total_catch = sum(catch, na.rm = T))


ebsbts_biomass<- ram_biomass %>%
  filter(region == "US Alaska") %>%
  filter(!str_detect(stocklong.x, "Gulf of Alaska")) %>%
  filter(scientificname %in% included_species) %>%
  group_by(year) %>%
  summarise(total_biomass = sum(biomass, na.rm = T))

ebsbts_effort <- gfw_data %>%
  filter(survey %in% c("ebsbts")) %>%
  unnest() %>%
  filter(
    (
      inferred_label_allyears == "trawlers" |
        inferred_label_allyears == "pots_and_traps" |
        str_detect(inferred_label_allyears, "set_") |
        inferred_label_allyears == "trawlers"
    ) &
      !inferred_label_allyears %in% c("cargo_or_tanker", "passenger")
  ) %>%
  group_by(year) %>%
  summarise(total_effort = sum(total_hours, na.rm = T))

survey_biomass <- candidate_data %>%
  filter(fished_only == TRUE, unfished_only == FALSE,
         survey_months_only == TRUE) %>%
  select(total_fish_data) %>%
  unnest() %>%
  filter(survey %in% "ebsbts") %>%
  unnest() %>%
  group_by(year) %>%
  summarise(total_survey_biomass = sum(biomass, na.rm = T))

ebsbts_cpue <- ebsbts_catch %>%
  left_join(ebsbts_effort, by = "year") %>%
  left_join(ebsbts_biomass, by = "year") %>%
  left_join(survey_biomass, by = "year") %>%
  left_join(fao_catch, by = "year") %>%
  left_join(nmfs_catch, by = "year") %>%
  mutate(ram_cpue = total_catch / total_effort,
         fao_cpue = total_fao_catch / total_effort,
         nmfs_cpue = total_nmfs_catch / total_effort)

ebsbts_cpue_plot <- ebsbts_cpue %>%
  filter(!is.na(ram_cpue)) %>%
  gather(metric, value, -year) %>%
  filter(value > 0, str_detect(metric, "cpue|biomass")) %>%
  mutate(Index = ifelse(str_detect(metric, "cpue"), "CPUE", "Abudance")) %>%
  group_by(metric) %>%
  mutate(scaled_value = value / max(value, na.rm = T)) %>%
  ungroup() %>%
  ggplot(aes(year, scaled_value, color = Index, linetype = metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis_d(option = "E") +
  labs(title = "Eastern Bering Sea", x = "", y = "")

# west coast

wc_included_species <- vast_fish %>%
  filter(survey %in% c("wcgbts"))

wc_included_species <- unique(wc_included_species$spp) %>%
  str_replace_all("_"," ")


wc_nmfs_catch <- nmfs_landings %>%
  filter(
    state %in% c("California", "Washington", "Oregon"),
    scientific_name %in% wc_included_species
  ) %>%
  group_by(year) %>%
  summarise(total_nmfs_catch = sum(metric_tons, na.rm = T))


wc_fao_catch <- fao_capture %>%
  filter(fao_area_code == 77, scientific_name %in% wc_included_species) %>%
  group_by(year) %>%
  summarise(total_fao_catch = sum(capture, na.rm = T))

wc_ram_catch <- ram_catch %>%
  filter(region == "US West Coast") %>%
  filter(!str_detect(stocklong.x, "Gulf of Alaska")) %>%
  filter(scientificname %in% wc_included_species) %>%
  group_by(year) %>%
  summarise(total_catch = sum(catch, na.rm = T))

wc_biomass<- ram_biomass %>%
  filter(region == "US West Coast") %>%
  filter(!str_detect(stocklong.x, "Gulf of Alaska")) %>%
  filter(scientificname %in% included_species) %>%
  group_by(year) %>%
  summarise(total_biomass = sum(biomass, na.rm = T))

wc_effort <- gfw_data %>%
  filter(survey %in% c("wcgbts")) %>%
  unnest() %>%
  filter(
    (
      inferred_label_allyears == "trawlers" |
        inferred_label_allyears == "pots_and_traps" |
        str_detect(inferred_label_allyears, "set_") |
        inferred_label_allyears == "trawlers"
    ) &
      !inferred_label_allyears %in% c("cargo_or_tanker", "passenger")
  ) %>%
  group_by(year) %>%
  summarise(total_effort = sum(total_hours, na.rm = T))

wc_survey_biomass <- vast_fish %>%
  filter(survey %in% "wcgbts",
         spp %in% str_replace_all(included_species," ","_")) %>%
  mutate(abundance = map(vasterized_data,"time_index")) %>%
  select(survey, spp, abundance) %>%
  unnest() %>%
  group_by(Year) %>%
  summarise(total_survey_biomass = sum(abundance)) %>%
  rename(year = Year)

wc_cpue <- wc_ram_catch %>%
  left_join(wc_effort, by = "year") %>%
  left_join(wc_biomass, by = "year") %>%
  left_join(wc_survey_biomass, by = "year") %>%
  mutate(cpue = total_catch / total_effort)

wc_cpue <- wc_ram_catch %>%
  left_join(wc_effort, by = "year") %>%
  left_join(wc_biomass, by = "year") %>%
  left_join(survey_biomass, by = "year") %>%
  left_join(wc_fao_catch, by = "year") %>%
  left_join(wc_nmfs_catch, by = "year") %>%
  mutate(ram_cpue = total_catch / total_effort,
         fao_cpue = total_fao_catch / total_effort,
         nmfs_cpue = total_nmfs_catch / total_effort)

wc_cpue_plot <- wc_cpue %>%
  filter(!is.na(ram_cpue)) %>%
  gather(metric, value, -year) %>%
  filter(value > 0, str_detect(metric, "cpue|biomass")) %>%
  mutate(Index = ifelse(str_detect(metric, "cpue"), "CPUE", "Abudance")) %>%
  group_by(metric) %>%
  mutate(scaled_value = value / max(value, na.rm = T)) %>%
  ungroup() %>%
  ggplot(aes(year, scaled_value, color = Index, linetype = metric)) +
  geom_line(size = 1.5, show.legend = F) +
  scale_color_viridis_d(option = "E") +
  labs(title = "US West Coast", y = "Scaled Value", x = "Year")



# cost coefficients -------------------------------------------------------

structural <- skynet_models %>%
  filter(model == "structural",
         dep_var == "biomass",
         data_subset == "skynet_100km",
         variables == "gfw_only") %>%
  select(train_set, training_data, r2_training)

strct_ovp_plot <- structural %>%
  unnest() %>%
  filter(train_set %in% c("random","spatial")) %>%
  mutate(train_set = fct_reorder(train_set, rev(r2_training))) %>%
  group_by(train_set) %>%
  mutate(biomass = biomass / max(biomass),
         pred = pred / max(pred)) %>%
  ggplot(aes(biomass, pred)) +
  geom_point() +
  geom_text(aes(0.9, .75, label = glue::glue("R2 = {round(r2_training,2)}"))) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~train_set) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Observed", y = "Predicted")

strct_r2 <- structural %>%
  ggplot(aes(r2_training)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  labs(x = bquote(R^2), y = "Count", title = "B)")


# structural fits ---------------------------------------------------------


# latent variables

skynet_100km <- rescale_data(skynet_data, resolution = 50) %>%
  na.omit()

latent_data <- skynet_100km %>%
  mutate(location = paste0(survey, knot,sep = '_')) %>%
  filter(year == max(year))

cs_latent_data <-
  recipe(
    total_hours ~ dist_from_port + restricted_use_mpa + m_below_sea_level + dist_from_shore + location,
    data = latent_data %>% filter(total_hours > 0)
  ) %>%
  step_log(all_outcomes()) %>%
  step_center(all_numeric(), -all_outcomes()) %>%
  step_scale(all_numeric(), -all_outcomes()) %>%
  prep(data = latent_data, retain = TRUE) %>%
  juice() %>%
  na.omit()

latent_model <- brm(total_hours ~ 1 + dist_from_port + restricted_use_mpa + m_below_sea_level + dist_from_shore + (1|location),
data = cs_latent_data)

location_effects <- broom::tidy(latent_model) %>%
  filter(str_detect(term, "location")) %>%
  mutate(location = str_extract(term, pattern = "(?<=\\[).*(?=_)")) %>%
  mutate(knot = str_replace_all(location,"\\D","") %>% as.numeric(),
         survey = str_replace_all(location,"\\d","")) %>%
  mutate(location_effect = exp(estimate + std.error^2/2)) %>%
  select(survey, knot, location_effect)


a <- latent_data %>%
  left_join(location_effects, by = c("survey", "knot"))

linear_latent_plot <- a %>%
  mutate(density = density / max(density,na.rm = T),
         biomass = biomass / max(biomass, na.rm = T),
         location_effect = location_effect / max(location_effect, na.rm = T)) %>%
  ggplot(aes(biomass, location_effect)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2, color = "red") +
  geom_smooth(method = "lm", color = "blue", se = F) +
  labs(x = "Observed biomass", y = "Space Effects")


skynet_50km <- rescale_data(skynet_data %>% filter(survey == "wcgbts"), resolution = 50) %>%
  na.omit()

latent_data <- skynet_50km %>%
  mutate(location = paste(survey, knot,sep = 'x')) %>%
  filter(year == max(year)) %>%
  na.omit() %>%
  rename(price = aggregate_price,
         PortDist = dist_from_port)

cs_latent_data <-
  recipe(
    total_hours ~ PortDist + restricted_use_mpa + m_below_sea_level + dist_from_shore + location +
      price + mean_vessel_length,
    data = latent_data %>% filter(total_hours > 0)
  ) %>%
  # step_log(all_outcomes()) %>%
  step_center(all_numeric(), -all_outcomes(),-price) %>%
  step_scale(all_numeric(), -all_outcomes(),-price) %>%
  prep(data = latent_data, retain = TRUE) %>%
  juice() %>%
  na.omit()



oa_model <- bf(total_hours ~ 1 / exp(q) * log((price * exp(q)*exp(cpue)) / ( exp(cost))),
               q ~ 1,
               cost ~ 1 + PortDist + m_below_sea_level + dist_from_shore + mean_vessel_length,
               cpue ~ (1|location), nl = TRUE)

oa_prior <- c(
  prior(normal(0,1), nlpar = "q"),
  prior(normal(0,1), nlpar = "cpue"),
  prior(normal(0,1), nlpar = "cost")
)

oa_fit <-
  brm(
    formula = oa_model,
    data = cs_latent_data,
    family = Gamma(),
    prior = oa_prior,
    control = list(adapt_delta = 0.95),
    chains = 1,
    iter = 10000
  )



tidy_oa <- broom::tidy(oa_fit)

oa_location_effects <- broom::tidy(oa_fit)%>%
  filter(str_detect(term, "location")) %>%
  mutate(location = str_extract(term, pattern = "(?<=\\[).*(?=\\,)")) %>%
  mutate(location = str_replace_all(location,"x","")) %>%
  mutate(knot = str_replace_all(location,"\\D","") %>% as.numeric(),
         survey = str_replace_all(location,"\\d","")) %>%
  mutate(location_effect = exp(estimate + std.error^2/2)) %>%
  select(survey, knot, location_effect)


oa_preds <- latent_data %>%
  left_join(oa_location_effects, by = c("survey", "knot"))

strct_latent_plot <- oa_preds %>%
  mutate(density = density / max(density,na.rm = T),
         biomass = biomass / max(biomass),
         location_effect = location_effect / max(location_effect, na.rm = T)) %>%
  ggplot(aes(biomass, location_effect)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 2, color = "red") +
  geom_smooth(method = "lm", color = "blue", se = F) +
  labs(x = "Observed biomass", y = "Space Effects")


# predicting effort -------------------------------------------------------



# predicting density ------------------------------------------------------

### Which ML Model is Best?

kfold_preds <- readRDS(file = here::here("results",run_name, "kfold_preds.RDS"))

kfold_preds_summary <-  kfold_preds %>%
  filter(variables == "gfw_only") %>%
  select(-fitted_model) %>%
  unnest() %>%
  mutate(model = str_replace_all(model, '-gfw_only',"")) %>%
  select(Resample, model, variables, obs, pred,) %>%
  # mutate(model = glue::glue("{model}-{variables}")) %>%
  select(-variables) %>%
  group_by(Resample, model) %>%
  summarise(model_rmse = sqrt(mean((obs - pred)^2))) %>%
  rename(id = Resample) %>%
  spread(model, model_rmse)

rmse_mod <- perf_mod(kfold_preds_summary, seed = 4344, iter = 5000)

rmse_post <- tidy(rmse_mod)

stacked_summary <- kfold_preds_summary %>%
  gather(model, model_rmse, -id)

best_ml_model_plot <- rmse_post %>%
  ggplot(rmse_post) +
  geom_point(
    data = stacked_summary,
    aes(x = model, y = model_rmse),
    alpha = .5, col = "blue"
  ) +
  coord_flip() +
  labs(y = "Posterior Probability Distribution of RMSE")


best_ml_model <- rmse_mod$stan %>%
  tidy() %>%
  filter(str_detect(term, "model")) %>%
  mutate(model = str_replace_all(term, "model","")) %>%
  arrange((estimate))


best_ml_model <- best_ml_model$model[[1]]


# best data resolution

ml_resolution_plot <- skynet_models %>%
  filter(model == "ranger",
         dep_var == "biomass",
         variables == "gfw_only",
         test_sets %in% c("random","spatial"),
         data_subset %in% c("skynet", "skynet_25km", "skynet_100km")) %>%
  group_by(data_subset) %>%
  mutate(mean_r2_training = mean(r2_training)) %>%
  ungroup() %>%
  mutate(data_subset = fct_reorder(data_subset, (mean_r2_training))) %>%
  ggplot(aes(data_subset, r2_training, fill = test_sets)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  ggsci::scale_fill_simpsons()



# best model performance --------------------------------------------------

ranger <- skynet_models %>%
  filter(model == "ranger",
         dep_var == "biomass",
         data_subset == "skynet_100km",
         variables == "gfw_only") %>%
  select(train_set, training_data, r2_training)


ranger_ovp_plot <- ranger %>%
  unnest() %>%
  filter(train_set %in% c("random","spatial")) %>%
  mutate(train_set = fct_reorder(train_set, rev(r2_training))) %>%
  group_by(train_set) %>%
  mutate(biomass = biomass / max(biomass),
         pred = pred / max(pred)) %>%
  ggplot(aes(biomass, pred)) +
  geom_point() +
  geom_text(aes(0.9, .75, label = glue::glue("R2 = {round(r2_training,2)}"))) +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~train_set) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Observed", y = "Predicted")

ranger_r2 <- ranger %>%
  ggplot(aes(r2_training)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  labs(x = bquote(R^2), y = "Count", title = "B)")


# best_ml_model <- "gbm"

best_fits <- skynet_models %>%
  filter(model == best_ml_model,
         dep_var == "biomass",
         variables == "gfw_only") %>%
  mutate(data_subset = if_else(data_subset == "lag_1_skynet_data", "skynet", data_subset))

# comparison so structural ------------------------------------------------


sub_skynet_models <- skynet_models %>%
  mutate(
    run_name = glue::glue(
      "{.$dep_var}--{.$model}--{.$weight_surveys}--{.$test_sets}--{.$data_subset}"
    )
  ) %>%
  filter(model %in% c(best_ml_model, "structural"),
         dep_var == "biomass",
         data_subset == "skynet",
         test_set == "random",
         variables == "gfw_only") %>%
  mutate(training_data = map(training_data, ~.x %>% filter(surveyed_year == T)))


downscaled_performance <-
  cross_df(list(
    run_name = sub_skynet_models$run_name,
    resolution =  seq(25, 200, by = 50)
  )) %>%
  left_join(
    sub_skynet_models %>% select(
      run_name,
      training_data,
      dep_var,
      model,
      weight_surveys,
      test_sets,
      data_subset
    ),
    by = "run_name"
  ) %>%
  arrange(run_name)

downscaled_performance <-  downscaled_performance %>%
  mutate(
    new_grid  = map2(
      training_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  )

downscaled_performance <-  downscaled_performance %>%
  mutate(
    new_data = map2(
      training_data,
      new_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat
    )
  )

r2foo <- function(x) {
  r2 <-
    yardstick::rsq(x, truth = total_biomass, estimate = agg_pred)

}


downscaled_performance <- downscaled_performance %>%
  mutate(oob_r2 = map_dbl(new_data, r2foo))


downscale_plot <- downscaled_performance %>%
  group_by(model) %>%
  mutate(meanr2 = mean(oob_r2, na.rm = T)) %>%
  ungroup() %>%
  mutate(model = fct_reorder(model, -meanr2)) %>%
  ggplot(aes(resolution, oob_r2, color = model)) +
  geom_line(show.legend = F) +
  geom_point(size = 2, show.legend = F) +
  labs(y = bquote(R^2), x = bquote("Resolution"~(km^2)), title = "A)") +
  theme(axis.title.y = element_text(angle = 0))


upscale_plot <- skynet_models %>%
  ungroup() %>%
  filter(test_sets == "random", data_subset %in% c("skynet","skynet_25km","skynet_100km"), dep_var == "biomass",!str_detect(data_subset,"delta"),
         model %in% c(best_ml_model,"structural"),
         variables == "gfw_only") %>%
  group_by(model) %>%
  mutate(meanr2 = mean(r2_training, na.rm = T)) %>%
  ungroup() %>%
  mutate(model = fct_reorder(model, -meanr2)) %>%
  mutate(data_subset = fct_relevel(data_subset,c("skynet","skynet_25km","skynet_100km"))) %>%
  ggplot(aes(data_subset,r2_training, fill = model)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  scale_x_discrete(labels = c("Raw", bquote("25"~km^2), bquote("100"~km^2))) +
  labs(title = "B)", y = bquote(R^2)) +
  theme(axis.title.y = element_blank()) +
  coord_flip()

upscale_plot <-  upscale_plot + coord_flip()

strct_v_ml_plot <- downscale_plot + upscale_plot


# value of information ----------------------------------------------------

#first with training

voi_plot <- skynet_models %>%
  filter(dep_var == "biomass", model == "ranger",
         test_sets %in% c("random"),
         data_subset %in% c("skynet","skynet_25km","skynet_100km")) %>%
  select(variables, data_subset, test_sets, r2_training) %>%
  spread(variables, r2_training) %>%
  mutate(delta_g_and_e = gfw_and_enviro - enviro_only,
         delta_g = gfw_only - enviro_only) %>%
  select(-enviro_only, -gfw_and_enviro, -gfw_only) %>%
  gather(comparison, delta, delta_g_and_e, delta_g) %>%
  ggplot(aes(comparison, delta, fill = data_subset)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  scale_x_discrete(labels = c("Effort vs Environment", "Effort and Evironment vs. Environment")) +
  coord_flip() +
  labs(y = bquote("Difference in"~R^2), x = "") +
  ggsci::scale_fill_startrek(name = "Data Resolution")


  ggplot(aes(delta, fill = comparison)) +
  geom_vline(aes(xintercept = 0), linetype = 2, color = "red") +
  geom_density(alpha = 0.5) +
  ggsci::scale_fill_npg(labels = c("GFW Only", "GFW and Environmental Data"))

#then with testing
#
# suggests that in general just enviro would be your best bet, how does that stand up
# to testing data

test_training_plot <- skynet_models %>%
  filter(dep_var == "biomass",
         data_subset == "skynet_100km",
         model %in% c(best_ml_model, 'structural'),
         test_sets %in% c("random","spatial","california")) %>%
  filter(!(model == "structural" & variables != "gfw_only")) %>%
  mutate(variables = fct_reorder(variables,r2 )) %>%
  select(variables, data_subset, test_sets,model, r2_training, r2) %>%
  gather(`Split`,r2, r2_training, r2) %>%
  ggplot(aes(variables, r2)) +
  geom_line() +
  geom_point(aes(fill = `Split`),shape = 21, size = 4) +
  labs(y = bquote(R^2),x = "Variables") +
  ggsci::scale_fill_rickandmorty(labels = c("Testing", "Training")) +
  coord_flip() +
  facet_grid(test_sets~model)


# methods and data --------------------------------------------------------------------

knots <- vast_fish %>%
  filter(!survey %in% bad_surveys) %>%
  select(survey,spp,vasterized_data) %>%
  mutate(knots = map(vasterized_data, c('spatial_list', 'loc_x_lat_long'))) %>%
  select(-vasterized_data,-spp) %>%
  unnest() %>%
  unique() %>%
  mutate(recenter_lon = ifelse(approx_long < 0, 180 + (180 - abs(approx_long)), approx_long))

map_knots <-  knots %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, approx_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(map_knots)

knot_plot <- map_knots %>%
  ggplot() +
  geom_sf(data = pacific_map, fill = 'grey60') +
  geom_sf(aes(color = survey), size = 1, alpha = 0.5) +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  scale_color_viridis_d() +
  theme(legend.key.height = unit(1.5,"cm"))


fish_summary <- vast_fish %>%
  filter(!survey %in% bad_surveys) %>%
  select(survey, data,spp) %>%
  unnest() %>%
  group_by(survey, spp) %>%
  summarise(n_seen = sum(catch_kg > 0)) %>%
  group_by(spp) %>%
  mutate(tc = sum(n_seen)) %>%
  ungroup() %>%
  mutate(spp = forcats::fct_reorder(spp, tc)) %>%
  group_by(survey) %>%
  top_n(10, tc) %>%
  ungroup()

fish_plot <- fish_summary %>%
  ggplot(aes(spp, tc, fill = survey)) +
  geom_col() +
  coord_flip() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank()) +
  labs(y = "# of times observed") +
  scale_fill_viridis_d()


fish_abundance <- skynet_data %>%
  filter(!survey %in% bad_surveys) %>%
  group_by(survey,knot, rounded_lat, rounded_lon) %>%
  summarise(total_density = mean(density),
            total_revenue = mean(economic_density)) %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon))

fish_abundance <-  fish_abundance %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()


bbox <- sf::st_bbox(map_knots)

fish_map_plot <- fish_abundance %>%
  ggplot() +
  geom_sf(aes(color = total_density), size = 1, alpha = 0.5) +
  geom_sf(data = pacific_map, fill = 'grey60') +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  # scale_color_viridis(name =bquote("Density"~(ton/km^2)))+
  scale_color_viridis(trans = "log10", name =bquote("Density"~(ton/km^2)))+
  theme(legend.key.height = unit(1.5,"cm"))

revenue_map_plot <- fish_abundance %>%
  ggplot() +
  geom_sf(aes(color = total_revenue), size = 1, alpha = 0.5) +
  geom_sf(data = pacific_map, fill = 'grey60') +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  # scale_color_viridis(name =bquote("Density"~(ton/km^2)))+
  scale_color_viridis(trans = "log10", name =bquote("Revenue Density"~(ton/km^2)), option = "A")+
  theme(legend.key.height = unit(1.5,"cm"))


price_plot <- species_prices %>%
  filter(is.na(mean_exvessel_price) == F) %>%
  mutate(species = fct_reorder(species, mean_exvessel_price))  %>%
  ggplot(aes(species, mean_exvessel_price)) +
  geom_col() +
  scale_y_continuous(labels = scales::dollar) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank()) +
  labs(y = "Mean ex-vessel price (USD)")


gfw_effort <- gfw_data %>%
  unnest() %>%
  filter(inferred_label_allyears == "trawlers",
         survey != "wcgbts" & survey != "wcghl") %>%
  group_by(year, rounded_lat, rounded_lon) %>%
  summarise(total_effort = pmin(10,sum(total_hours, na.rm = T))) %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) %>%
  filter(total_effort > 0, year == max(year))

gfw_map <-  gfw_effort %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()


bbox <- sf::st_bbox(gfw_map)

gfw_plot <- gfw_map %>%
  ggplot() +
  geom_sf(aes(color = total_effort), size = 0.1, alpha = 0.25) +
  geom_sf(data = pacific_map, fill = 'grey60') +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin'], bbox['ymax'])) +
  # scale_color_viridis(name = "Fishing Hours", labels = scales::comma) +
  scale_color_viridis(name = "Fishing Hours", labels = scales::comma) +
  theme(legend.key.height = unit(1.5,"cm"))


# save --------------------------------------------------------------------



plots <- ls()[str_detect(ls(),"_plot")]


save(file = paste0(run_dir,"/skynet-plots.Rdata"), list = plots)



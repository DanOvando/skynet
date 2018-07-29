library(gbm)
library(randomForest)
library(stringr)
library(modelr)
library(broom)
library(hrbrthemes)
library(viridis)
library(sf)
library(ggmap)
library(caret)
library(extrafont)
library(patchwork)
library(scico)
library(wesanderson)
library(ggridges)
library(taxize)
library(tidyposterior)
library(recipes)
library(rstan)
library(tidyverse)


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
  ) %>%
  mutate(
    test_plot = pmap(
      list(
        data = test_data,
        r2 = r2,
        test_region = test_sets,
        train_region = train_set,
        data_set = data_subset,
        dep_var = dep_var
      ),
      diagnostic_plot_foo
    ),
    training_plot = pmap(
      list(
        data = training_data,
        r2 = r2_training,
        test_region = paste0(test_sets),
        train_region = paste0(train_set, "- training plot"),
        data_set = data_subset,
        dep_var = dep_var
      ),
      diagnostic_plot_foo
    )
  )

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

strct_ovp <- structural %>%
  unnest() %>%
  mutate(train_set = fct_reorder(train_set, rev(r2_training))) %>%
  group_by(train_set) %>%
  mutate(biomass = biomass / max(biomass),
         pred = pred / max(pred)) %>%
  ggplot(aes(biomass, pred)) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "red") +
  facet_wrap(~train_set) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Observed", y = "Predicted", title = "A)")

strct_r2 <- structural %>%
  ggplot(aes(r2_training)) +
  geom_histogram(fill = "lightgrey", color = "black") +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  labs(x = bquote(R^2), y = "Count", title = "B)")


# structural fits ---------------------------------------------------------




# predicting effort -------------------------------------------------------

data <- skynet_data %>%
  select(any_fishing, dist_from_port,dist_from_shore,
         no_take_mpa,restricted_use_mpa, density, mean_analysed_sst,
         mean_chla )

fish_recipe <- recipe(any_fishing ~ . , data = data) %>%
  step_center(all_predictors(), -no_take_mpa,-restricted_use_mpa) %>%
  step_scale(all_predictors(), -no_take_mpa,-restricted_use_mpa)

prepped_data <- prep(fish_recipe, data = data, retain = T)

prepped_data <- prepped_data %>% juice()

test <- glm(any_fishing ~., data = prepped_data, family = binomial)

summary(test)


# latent spatial variables


a <- candidate_data %>%
  filter(fished_only == TRUE,
         unfished_only == FALSE,
         survey_months_only == TRUE)

data <- a$skynet_data[[1]] %>%
  filter(survey == "wcgbts") %>% mutate(state = ifelse(
    survey == "wcgbts",
    if_else(
      rounded_lat > 46.25,
      "washington",
      if_else(rounded_lat > 42, "oregon", "california")
    ),
    "Alaska"
  ))

oregon <- data %>%
  filter(state == "oregon", total_hours > 0) %>%
  filter(year == max(year))

oregon_cost <- oregon %>%
  select(dist_from_port, dist_from_shore, m_below_sea_level)

cost_recipe <- recipe(oregon_cost) %>%
  step_log(all_numeric())
# step_center(all_numeric()) %>%
# step_scale(all_numeric())
#
oregon_cost <-
  prep(cost_recipe, oregon_cost, retain = T) %>% juice() %>%
  mutate(intercept = 1)


location_year <- oregon %>%
  select(knot, year) %>%
  mutate(knot_year = paste(knot, year,sep = '_')) %>%
  select(knot_year) %>%
  mutate(marker = 1) %>%
  mutate(index = 1:nrow(.)) %>%
  spread(knot_year, marker, fill = 0) %>%
  select(-index)

warmups <- 5000

total_iterations <- 10000

max_treedepth <-  12

adapt_delta <-  0.8

chains <- 1

struct_data <- list(
  cost_data = oregon_cost,
  location_data = location_year,
  log_effort = as.numeric(oregon$total_hours %>% log()),
  price = oregon$aggregate_price,
  mp = 0,
  n = nrow(oregon),
  max_q = .8/max(oregon$total_hours)
)

struct_data$n_betas <- ncol(struct_data$cost_data)

struct_data$n_cpue_betas <- ncol(struct_data$location_data)

# inits <- list(list(cpue_betas = rep(25, struct_data$n_cpue_betas)))


# stan_fit <-
#   stan(
#     file = here::here("src", "fit_effort_model.stan"),
#     data = struct_data,
#     chains = chains,
#     warmup = warmups,
#     iter = total_iterations,
#     cores = 1,
#     refresh = 250,
#     control = list(max_treedepth = max_treedepth,
#                    adapt_delta = adapt_delta)
#   )
#
# fits <- stan_fit %>%
#   tidy()
#
# knots <- fits %>%
#   filter(str_detect(term, "knot"))
#
# effort <- fits %>%
#   filter(str_detect(term, "effort") & !str_detect(term, "log")) %>%
#   mutate(effort = oregon$total_hours)
#
# effort %>%
#   ggplot(aes(effort, estimate)) +
#   geom_point()


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


best_ml_model <- str_extract(best_ml_model$model[[1]]
                             ,".*(?=-)")


# best model performance --------------------------------------------------

ranger <- skynet_models %>%
  filter(model == "ranger",
         dep_var == "biomass",
         data_subset == "skynet_100km",
         variables == "gfw_only") %>%
  select(train_set, training_data, r2_training)

ranger_ovp <- ranger %>%
  unnest() %>%
  mutate(train_set = fct_reorder(train_set, rev(r2_training))) %>%
  group_by(train_set) %>%
  mutate(biomass = biomass / max(biomass),
         pred = pred / max(pred)) %>%
  ggplot(aes(biomass, pred)) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2, color = "red") +
  facet_wrap(~train_set) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Observed", y = "Predicted", title = "A)")

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
         test_sets %in% test_sets_of_interest) %>%
  select(variables, data_subset, test_sets, r2_training) %>%
  spread(variables, r2_training) %>%
  mutate(delta_g_and_e = gfw_and_enviro - enviro_only,
         delta_g = gfw_only - enviro_only) %>%
  select(-enviro_only, -gfw_and_enviro, -gfw_only) %>%
  gather(comparison, delta, delta_g_and_e, delta_g) %>%
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
         model == 'ranger',
         test_sets %in% test_sets_of_interest) %>%
  mutate(variables = fct_reorder(variables,r2 )) %>%
  select(variables, data_subset, test_sets, r2_training, r2) %>%
  gather(`Split`,r2, r2_training, r2) %>%
  ggplot(aes(variables, r2)) +
  geom_line() +
  geom_point(aes(fill = `Split`),shape = 21, size = 4) +
  labs(y = bquote(R^2),x = "Variables") +
  ggsci::scale_fill_rickandmorty(labels = c("Training", "Testing")) +
  coord_flip() +
  facet_wrap(~test_sets)


# which data subset

best_data_subset <- best_fits %>%
  filter(dep_var == "density", unfished_only == FALSE) %>%
  group_by(data_subset) %>%
  mutate(median_r2 = median(r2_training, na.rm = T)) %>%
  ungroup() %>%
  mutate(data_subset = fct_reorder(data_subset, median_r2)) %>%
  arrange(desc(median_r2))

best_data_subset %>%
  ggplot() +
  geom_point(position = "dodge",aes(data_subset,r2_training), alpha = 0.5) +
  geom_point(aes(data_subset, median_r2), color = "red", size = 2) +
  coord_flip()

best_data_subset <- "skynet"


best_fits <- skynet_models %>%
  filter(model == best_ml_model | model == "structural",
         dep_var == "density",
         variables == "gfw_only") %>%
  mutate(data_subset = if_else(data_subset == "lag_1_skynet_data", "skynet", data_subset))

best_fits <- best_fits %>%
  filter(data_subset == best_data_subset)

best_fits %>%
  group_by(test_set) %>%
  mutate(mean_r2 = mean(r2, na.rm = T)) %>%
  ungroup() %>%
  mutate(test_set = fct_reorder(test_set, mean_r2)) %>%
  gather(r2_source, r2_value, r2, r2_training) %>%
  ggplot(aes(test_set, r2_value)) +
  geom_line() +
  geom_point(aes(color = r2_source), size = 4) +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  labs(y = bquote(R ^ 2)) +
  scale_color_manual(
    name = "Data Split",
    labels = c("Testing", "Training"),
    values = wes_palette("Zissou1")
  ) +
  facet_wrap(~model)


# spatial performance -----------------------------------------------------

spatial_performance <- best_fits %>%
  filter(test_set == "spatial")  %>% {
    .$test_data[[1]]
  } %>%
  filter(surveyed_year == T) %>%
  group_by(survey) %>%
  mutate(sd_density = sd(density, na.rm = T)) %>%
  group_by(rounded_lat, rounded_lon, survey) %>%
  summarise(resid = mean(pred - density),
            scaled_resid = pmax(-2,pmin(2,resid / mean(sd_density))))

spatial_performance_plot <- spatial_performance %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(spatial_performance_plot)

alaska_bbox <- sf::st_bbox(spatial_performance_plot %>% filter(survey != "wcgbts"))

wc_bbox <- sf::st_bbox(spatial_performance_plot %>% filter(survey == "wcgbts"))


resid_hist <- spatial_performance %>%
  ggplot(aes(scaled_resid)) +
  geom_density(color = "black", fill = "grey", alpha = 0.5) +
  labs(x = "Scaled Residuals", caption = "Residuals divided by standard deviation of observed densities") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text.x = element_text(size = 8))



alaska_spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(survey != "wcgbts"), aes(color = scaled_resid),size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map, shape = 21) +
  coord_sf(xlim = c(alaska_bbox['xmin'], alaska_bbox['xmax']),
           ylim = c(alaska_bbox['ymin']*.95, alaska_bbox['ymax'] * 1.1), expand = F) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "grey", name = "Scaled Residuals", guide = guide_colorbar(frame.colour = "black")) +
  map_theme +
  labs(caption = "Alaska") +
  theme(legend.key.height = unit(1,"cm"))


wc_spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(survey == "wcgbts"), aes(color = scaled_resid), size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map) +
  coord_sf(xlim = c(wc_bbox['xmin']*.98, wc_bbox['xmax']),
           ylim = c(wc_bbox['ymin'], wc_bbox['ymax'] * 1.05), expand = F) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "grey", name = "Residuals", guide = "none") +
  map_theme +
  labs(caption = "West Coast")


spatial_residual_plot <- (
  wc_spatial_residual_plot + alaska_spatial_residual_plot
) / resid_hist  + plot_layout(nrow = 2, ncol = 1, widths = c(2,1),heights = c(2.5,1)) & theme(plot.margin = unit(rep(0.01,4), "points"))


# more spatial ------------------------------------------------------------

train_performance <- best_fits %>%
  filter(test_set == "california", data_subset == "skynet")  %>%{
    .$training_data[[1]]
  } %>%
  mutate(resid = pred - density) %>%
  mutate(split = "Training")

test_performance <- best_fits %>%
  filter(test_set == "california", data_subset == "skynet")  %>% {
    .$test_data[[1]]
  } %>%
  mutate(resid = pred - density) %>%
  mutate(split = "Testing")

spatial_performance <- train_performance %>%
  bind_rows(test_performance) %>%
  filter(surveyed_year == T) %>%
  group_by(survey) %>%
  mutate(scaled_resid = resid / sd(density))

spatial_performance_plot <- spatial_performance %>%
  mutate(recenter_lon = ifelse(rounded_lon < 0, 180 + (180 - abs(rounded_lon)), rounded_lon)) %>%
  dplyr::mutate(geometry = purrr::map2(recenter_lon, rounded_lat, ~ sf::st_point(x = c(.x, .y), dim = 'XY'))) %>%
  ungroup() %>%
  mutate(geometry = sf::st_sfc(geometry, crs =
                                 "+proj=longlat +datum=WGS84 +no_defs")) %>%
  sf::st_sf()

bbox <- sf::st_bbox(spatial_performance_plot)


resid_hist <- spatial_performance %>%
  ggplot(aes(scaled_resid,fill = split)) +
  geom_density(alpha = 0.75) +
  labs(x = "Scaled Residuals", title = "A)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  theme(legend.position = "right") +
  scale_fill_discrete(name = '')


spatial_residual_plot <- ggplot() +
  geom_sf(data = spatial_performance_plot %>% filter(split == "Testing"), aes(color = scaled_resid), size = 0.5, alpha = 0.75) +
  geom_sf(data = pacific_map) +
  coord_sf(xlim = c(bbox['xmin'], bbox['xmax']),
           ylim = c(bbox['ymin']*.75, bbox['ymax'] * 1.1), expand = F) +
  scale_shape(guide = FALSE) +
  scale_color_gradient2(low = "tomato", high = "steelblue",midpoint = 0, mid = "lightgrey", name = "Residuals") +
  map_theme +
  labs(title = "B)")


resid_hist + spatial_residual_plot + plot_layout(ncol = 2, nrow = 1, widths = c(1,2))



# no idea -----------------------------------------------------------------

resolution_km2 <- 100

sub_skynet_models <- skynet_models %>%
  mutate(
    run_name = glue::glue(
      "{.$dep_var}--{.$model}--{.$weight_surveys}--{.$test_sets}--{.$data_subset}"
    )
  ) %>%
  filter(model %in% c(best_ml_model),
         dep_var == "density",
         data_subset == "skynet",
         test_set %in% c("random","spatial","alaska", "random_west_coast", "west_coast","year_greq_than_2014"),
         variables == "gfw_only")  %>%
  mutate(test_data = map(test_data, ~.x %>% filter(surveyed_year == T)),
         training_data = map(training_data, ~.x %>% filter(surveyed_year == T)))


downscaled_results <-
  cross_df(list(
    run_name = sub_skynet_models$run_name,
    resolution = resolution_km2
  )) %>%
  left_join(
    sub_skynet_models %>% select(
      run_name,
      test_data,
      training_data,
      dep_var,
      model,
      weight_surveys,
      train_set,
      test_set,
      data_subset
    ),
    by = "run_name"
  ) %>%
  arrange(run_name)

downscaled_results <-  downscaled_results %>%
  mutate(
    training_grid  = map2(
      training_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  ) %>%
  mutate(
    training_data = map2(
      training_data,
      training_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat,
      dep_var = "density"
    )
  )

downscaled_results <-  downscaled_results %>%
  mutate(
    testing_grid  = map2(
      test_data,
      resolution,
      create_grid,
      lon_name = rounded_lon,
      lat_name = rounded_lat
    )
  ) %>%
  mutate(
    test_data = map2(
      test_data,
      testing_grid,
      snap_to_grid,
      old_lon_name = rounded_lon,
      old_lat_name = rounded_lat,
      new_lon_name = lon,
      new_lat_name = lat,
      dep_var = "density"
    )
  )

train_performance <- downscaled_results %>%
  select(test_set,train_set, training_data) %>%
  unnest() %>%
  group_by(train_set) %>%
  mutate(rel_pred = agg_pred / max(agg_pred),
         rel_density = agg_mean_density / max(agg_mean_density)) %>%
  mutate(split = "Training")

test_performance <- downscaled_results %>%
  select(test_set,train_set, test_data) %>%
  unnest() %>%
  group_by(test_set) %>%
  mutate(rel_pred = agg_pred / max(agg_pred),
         rel_density = agg_mean_density / max(agg_mean_density)) %>%
  mutate(split = "Testing")

spatial_performance <- train_performance %>%
  bind_rows(test_performance) %>%
  ungroup() %>%
  filter(!str_detect(test_set, "year_")) %>%
  mutate(set_order = as.numeric(factor(test_set))) %>%
  mutate(test_set = fct_reorder(test_set, set_order),
         train_set = fct_reorder(train_set, set_order))

trained_plot <- spatial_performance %>%
  filter(split == "Training") %>%
  ggplot(aes(rel_density, rel_pred)) +
  geom_point() +
  facet_wrap(~ train_set, strip.position = "left", nrow = n_distinct(spatial_performance$test_set), ncol = 1 ) +
  labs(title = "Trained on...", y = "Predicted")  +
  theme(axis.title.x = element_blank())


tested_plot <- spatial_performance %>%
  filter(split == "Testing") %>%
  ggplot(aes(rel_density, rel_pred)) +
  geom_point() +
  facet_wrap(~ test_set, strip.position = "right", nrow = n_distinct(spatial_performance$test_set), ncol = 1 ) +
  labs(title = "Tested on...", x = "Observed") +
  theme(axis.title.y = element_blank())

trained_plot + tested_plot &
  theme(
    panel.spacing = unit(10, "points"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 7),
    plot.margin = unit(rep(.1, 4), "lines")
  )


## Value of information

skynet_models %>%
  filter(data_subset == "skynet", model == best_ml_model,
         test_set == "spatial", dep_var == "density") %>%
  select(variables, r2, r2_training) %>%
  gather(r2_source, r2_value, -variables) %>%
  ggplot(aes(variables, r2_value, fill = r2_source)) +
  geom_col(position = "dodge")

skynet_models %>%
  filter(data_subset == "skynet", model == best_ml_model,
         test_set == "spatial", dep_var == "density") %>%
  select(variables, rmse, rmse_training) %>%
  gather(rmse_source, rmse_value, -variables) %>%
  ggplot(aes(variables, rmse_value, fill = rmse_source)) +
  geom_col(position = "dodge")



## Supplementary predictions

generate_test_training <-
  function(dat,
           test_set,
           prop = 0.75,
           cut_year = 2014,
           keep_every = 25) {

     if (test_set == "spatial_alaska") {

       temp_dat <- dat %>%
         filter(!(dat$survey %in% c('wcgbts', 'wcghl'))) %>%
         arrange(rounded_lon, rounded_lat) %>%
         mutate(marker = 1:nrow(.)) %>%
         mutate(keepers = (marker %% keep_every) == 0) %>%
         filter(keepers)

       train_dex <- temp_dat$index

       test_dex <- dat$index[!(dat$index %in% train_dex) & (!dat$survey %in% c('wcgbts', 'wcghl')) &
                               dat$surveyed_year == TRUE]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'spatial_alaska') %>%
        mutate(test_set = 'spatial_alaska')
    }

    if (test_set == "spatial_west_coast") {


      temp_dat <- dat %>%
        filter((dat$survey %in% c('wcgbts', 'wcghl'))) %>%
        arrange(rounded_lon, rounded_lat) %>%
        mutate(marker = 1:nrow(.)) %>%
        mutate(keepers = (marker %% keep_every) == 0) %>%
        filter(keepers)

      train_dex <- temp_dat$index

      test_dex <- dat$index[!(dat$index %in% train_dex) & (dat$survey %in% c('wcgbts', 'wcghl')) &
                              dat$surveyed_year == TRUE]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'spatial_west_coast') %>%
        mutate(test_set = 'spatial_west_coast')
    }

    if (test_set == "spatial") {

      temp_dat <- dat %>%
        arrange(rounded_lat, rounded_lon) %>%
        mutate(marker = 1:nrow(.)) %>%
        mutate(keepers = (marker %% keep_every) == 0) %>%
        filter(keepers)
# browser()
# coordinates(temp_dat) = ~rounded_lon + rounded_lat
#
#       test <- variogram(density ~ rounded_lon + rounded_lat, temp_dat)
#
#
#       plot(test)
      train_dex <- temp_dat$index

      test_dex <- dat$index[!(dat$index %in% train_dex) &
                              dat$surveyed_year == TRUE]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'spatial') %>%
        mutate(test_set = 'spatial')
    }

    if (test_set == 'random') {

      splits <-
        rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <- dat$index[splits$in_id]

      surveyed_index <- dat$index[dat$surveyed_year == T]

      train_dex <- train_dex[train_dex %in% surveyed_index]

      test_dex <-
        dat$index[!dat$index %in% train_dex]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'random') %>%
        mutate(test_set = 'random')

    }

    if (test_set == 'random_alaska') {
      dat <- dat %>%
        filter(!dat$survey %in% c('wcgbts', 'wcghl'))

      splits <-
        rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <-
        dat$index[splits$in_id]

      surveyed_index <- dat$index[dat$surveyed_year == T]

      train_dex <- train_dex[train_dex %in% surveyed_index]

      test_dex <- dat$index[!dat$index %in% train_dex]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'random_alaska') %>%
        mutate(test_set = 'random_alaska')

    }
    if (test_set == 'random_west_coast') {
      dat <- dat %>%
        filter(dat$survey %in% c('wcgbts', 'wcghl'))

      splits <-
        rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <-
        dat$index[splits$in_id]

      surveyed_index <- dat$index[dat$surveyed_year == T]

      train_dex <- train_dex[train_dex %in% surveyed_index]

      test_dex <- dat$index[!dat$index %in% train_dex]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'random_west_coast') %>%
        mutate(test_set = 'random_west_coast')

    }

    if (test_set == 'california') {

      train <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'), surveyed_year == TRUE, rounded_lat >= 42)

      test <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'), rounded_lat < 42)

      train_dex <- train$index

      test_dex <- test$index

      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'california'
        ) %>%
        mutate(train_set = 'wa-or')
    }

    if (test_set == 'wa-or') {

      train <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'), surveyed_year == TRUE, rounded_lat < 42)

      test <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'), rounded_lat >= 42)

      train_dex <- train$index

      test_dex <- test$index

      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'wa-or'
        ) %>%
        mutate(train_set = 'california')
    }


    if (test_set == 'west_coast') {

      train <- dat %>%
        filter(!survey %in% c('wcgbts', 'wcghl'), surveyed_year == TRUE)

      test <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'))

      train_dex <- train$index

      test_dex <- test$index

      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'west_coast'
        ) %>%
        mutate(train_set = 'not_west_coast')
    }

    if (test_set == 'alaska') {

      train <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'), surveyed_year == TRUE)

      test <- dat %>%
        filter(!survey %in% c('wcgbts', 'wcghl'))

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'alaska'
        ) %>%
        mutate(train_set = 'west_coast')


    }
    if (test_set == 'historic') {
      train <- dat %>%
        filter(year < cut_year, surveyed_year == TRUE)

      test <- dat %>%
        filter(year >= cut_year)

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = paste0('year_greq_than_', cut_year)
        ) %>%
        mutate(train_set =  paste0('year_l_than_', cut_year))


    }

    if (test_set == 'ebs-ai') {
      train <- dat %>%
        filter(survey %in% c('goabts'), surveyed_year == TRUE)

      test <- dat %>%
        filter(survey %in% c('ebsbts', 'aibts'))

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'ebs-ai'
        ) %>%
        mutate(train_set = 'goabts')

    }

    if (test_set == 'goa-ai') {
      train <- dat %>%
        filter(survey %in% c('ebsbts'), surveyed_year == TRUE)

      test <- dat %>%
        filter(survey %in% c('goabts', 'aibts'))

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = 'goa-ai'
        ) %>%
        mutate(train_set = 'ebsbts')

    }
    if (test_set == 'historic_alaska') {
      train <- dat %>%
        filter(year < cut_year,
               !survey %in% c('wcgbts', 'wcghl'),
               surveyed_year == TRUE)

      test <- dat %>%
        filter(year >= cut_year,!survey %in% c('wcgbts', 'wcghl'))

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = paste0('alaska_year_greq_than_', cut_year)
        ) %>%
        mutate(train_set =  paste0('alaska_year_l_than_', cut_year))


    }
    if (test_set == 'historic_west_coast') {
      train <- dat %>%
        filter(year < cut_year,
               survey %in% c('wcgbts', 'wcghl'),
               surveyed_year == TRUE)

      test <- dat %>%
        filter(year >= cut_year, survey %in% c('wcgbts', 'wcghl'))

      train_dex <- train$index

      test_dex <- test$index


      out <-
        data_frame(
          train = list(train_dex),
          test = list(test_dex),
          test_set = paste0('west_coast_year_greq_than_', cut_year)
        ) %>%
        mutate(train_set =  paste0('west_coast_year_l_than_', cut_year))


    }


    return(out)
  }
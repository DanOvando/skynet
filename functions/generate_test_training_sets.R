generate_test_training <-
  function(dat,
           test_set,
           prop = 0.75,
           cut_year = 2014) {

    if (test_set == "spatial_alaska"){
      dat <- dat %>%
        mutate(lat_dex = plyr::round_any(rounded_lat,0.25) %% 0.5 == 0,
               lon_dex = plyr::round_any(rounded_lon,0.25) %% 0.5 == 0)

      train_dex <-
        dat$index[dat$lat_dex == TRUE &
                    dat$lon_dex == TRUE & (!dat$survey %in% c('wcgbts', 'wcghl'))]

      test_dex <-
        dat$index[!(dat$index %in% train_dex) &
                    !(dat$survey %in% c('wcgbts', 'wcghl'))]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'spatial_alaska') %>%
        mutate(test_set = 'spatial_alaska')
    }

    if (test_set == "spatial_west_coast"){

      dat <- dat %>%
        mutate(lat_dex = plyr::round_any(rounded_lat,0.25) %% 0.5 == 0,
               lon_dex = plyr::round_any(rounded_lon,0.25) %% 0.5 == 0)
      train_dex <-
        dat$index[dat$lat_dex == TRUE &
                    dat$lon_dex == TRUE & (dat$survey %in% c('wcgbts', 'wcghl'))]

      test_dex <-
        dat$index[!(dat$index %in% train_dex) &
                    (dat$survey %in% c('wcgbts', 'wcghl'))]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'spatial_west_coast') %>%
        mutate(test_set = 'spatial_west_coast')
    }


     if (test_set == 'random') {

      splits <- rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <- dat$index[splits$in_id]

      test_dex <- dat$index[!dat$index %in% train_dex]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'random') %>%
        mutate(test_set = 'random')

     }

    if (test_set == 'random_alaska') {

      dat <- dat %>%
        filter(!dat$survey %in% c('wcgbts', 'wcghl'))

      splits <- rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <- dat$index[splits$in_id]

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

      splits <- rsample::initial_split(dat, prop = prop, strata = "survey")

      train_dex <- dat$index[splits$in_id]

      test_dex <- dat$index[!dat$index %in% train_dex]

      out <-
        data_frame(train = list(train_dex),
                   test = list(test_dex)) %>%
        mutate(train_set = 'random_west_coast') %>%
        mutate(test_set = 'random_west_coast')

    }


    if (test_set == 'west_coast') {

      train <- dat %>%
        filter(!survey %in% c('wcgbts', 'wcghl'))


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
        filter(survey %in% c('wcgbts', 'wcghl'))

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
        filter(year < cut_year)

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
        filter(survey %in% c('goabts'))

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
        filter(survey %in% c('ebsbts'))

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
        filter(year < cut_year, !survey %in% c('wcgbts', 'wcghl'))

      test <- dat %>%
        filter(year >= cut_year, !survey %in% c('wcgbts', 'wcghl'))

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
        filter(year < cut_year, survey %in% c('wcgbts', 'wcghl'))

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
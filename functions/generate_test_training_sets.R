generate_test_training <-
  function(dat,
           test_set,
           kfolds = 2,
           cut_year = 2014) {
    if (test_set == 'random') {
    browser()
      out <- modelr::crossv_kfold(dat, k = kfolds, id = 'test_set') %>%
        mutate(train_set = 'random') %>%
        mutate(test_set = 'random') %>%
        slice(1)

      out <- rsample::initial_split(data, prop = 0.75) %>%
        rsample::training()

      # out <- modelr::crossv_kfold(dat, k = kfolds, id = 'test_set') %>%
      #   mutate(train_set = 'random') %>%
      #   mutate(test_set = 'random') %>%
      #   slice(1)

    }
    if (test_set == 'west_coast') {
      train <- dat %>%
        filter(!survey %in% c('wcgbts', 'wcghl'))

      test <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'))

      out <-
        data_frame(
          train = list(resample(train, 1:nrow(train))),
          test = list(resample(test, 1:nrow(test))),
          test_set = 'west_coast'
        ) %>%
        mutate(train_set = 'not_west_coast')


    }

    if (test_set == 'alaska') {
      train <- dat %>%
        filter(survey %in% c('wcgbts', 'wcghl'))

      test <- dat %>%
        filter(!survey %in% c('wcgbts', 'wcghl'))

      out <-
        data_frame(
          train = list(resample(train, 1:nrow(train))),
          test = list(resample(test, 1:nrow(test))),
          test_set = 'alaska'
        ) %>%
        mutate(train_set = 'west_coast')


    }
    if (test_set == 'historic') {
      train <- dat %>%
        filter(year < cut_year)

      test <- dat %>%
        filter(year >= cut_year)

      out <-
        data_frame(
          train = list(resample(train, 1:nrow(train))),
          test = list(resample(test, 1:nrow(test))),
          test_set = paste0('year_greq_than_', cut_year)
        ) %>%
        mutate(train_set =  paste0('year_l_than_', cut_year))


    }

    if (test_set == 'ebs-ai') {
      train <- dat %>%
        filter(survey %in% c('goabts'))

      test <- dat %>%
        filter(survey %in% c('ebsbts', 'aibts'))

      out <-
        data_frame(
          train = list(resample(train, 1:nrow(train))),
          test = list(resample(test, 1:nrow(test))),
          test_set = 'ebs-ai'
        ) %>%
        mutate(train_set = 'goabts')

    }

    if (test_set == 'goa-ai') {
      train <- dat %>%
        filter(survey %in% c('ebsbts'))

      test <- dat %>%
        filter(survey %in% c('goabts', 'aibts'))

      out <-
        data_frame(
          train = list(resample(train, 1:nrow(train))),
          test = list(resample(test, 1:nrow(test))),
          test_set = 'goa-ai'
        ) %>%
        mutate(train_set = 'ebsbts')



    }

    return(out)
  }
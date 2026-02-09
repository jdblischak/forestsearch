
########################
# SECTION on De-bugging
########################

run_debug <- FALSE


if(run_debug){

  fs.est <- fs
  # Capture what one fold actually produces
  fs_args <- fs.est$args_call_all
  fs_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  fs_args$details <- TRUE
  fs_args$plot.sg <- FALSE

  # Replicate the data prep from Section 3
  confounders.name <- fs_args$confounders.name
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  get_names <- c(confounders.name, outcome.name, event.name, id.name, treat.name)

  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  # Replicate one fold
  set.seed(8316951 + 1000)
  id_sample <- sample(seq_len(nrow(dfnew)), replace = FALSE)
  df_scrambled <- dfnew[id_sample, ]
  folds <- cut(seq_len(nrow(df_scrambled)), breaks = 5, labels = FALSE)

  x.train <- df_scrambled[folds != 1, ]
  fs_args$df.analysis <- x.train

  result <- try(do.call(forestsearch, fs_args), silent = FALSE)
  print(result)


  # Continue from where we left off — simulate a full 5-fold run
  resCV_list <- vector("list", 5)

  for (cv_index in 1:5) {
    x.test <- df_scrambled[folds == cv_index, ]
    x.train <- df_scrambled[folds != cv_index, ]

    fs_args$df.analysis <- x.train
    fs.train <- suppressWarnings(try(do.call(forestsearch, fs_args), TRUE))

    if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {
      cat("Fold", cv_index, ": subgroup FOUND\n")
    } else {
      cat("Fold", cv_index, ": no subgroup\n")
      x.test$treat.recommend <- 1.0
      x.test$sg1 <- NA
      x.test$sg2 <- NA
      x.test$cvindex <- cv_index
    }

    resCV_list[[cv_index]] <- x.test
  }

  # Now test the part that actually fails:
  resCV <- do.call(rbind, resCV_list)
  res <- list(
    resCV = as.data.frame(resCV),
    Kfolds = 5,
    cv_args = fs_args,
    sg_analysis = fs.est$sg.harm,
    sg0.name = "Not recommend",
    sg1.name = "Recommend"
  )

  out <- try(forestsearch_KfoldOut(res = res, outall = FALSE, details = TRUE), silent = FALSE)
  print(out)



  # In your console, run ONE simulation manually without error swallowing:
  fs_args <- fs.est$args_call_all
  fs_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  fs_args$details <- FALSE
  fs_args$plot.sg <- FALSE

  confounders.name <- fs_args$confounders.name
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  get_names <- c(confounders.name, outcome.name, event.name, id.name, treat.name)

  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  Kfolds <- 5

  # Simulate one iteration
  df_scrambled <- data.table::copy(dfnew)
  set.seed(8316951 + 1000)
  id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
  df_scrambled <- df_scrambled[id_sample, ]
  folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)

  resCV_list <- vector("list", Kfolds)

  for (cv_index in seq_len(Kfolds)) {
    testIndexes <- which(folds == cv_index)
    x.test <- df_scrambled[testIndexes, ]
    x.train <- df_scrambled[-testIndexes, ]

    fs_args$df.analysis <- x.train
    fs.train <- suppressWarnings(try(do.call(forestsearch, fs_args), TRUE))

    if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {

      print(fs.train$sg.harm)

      cat("Fold", cv_index, ": sg.harm =", fs.train$sg.harm, "\n")
      cat("Fold", cv_index, ": x.test cols =", head(names(x.test), 20), "\n")

      sg1 <- fs.train$sg.harm[1]
      sg2 <- if (length(fs.train$sg.harm) > 1) fs.train$sg.harm[2] else NA
      df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
      df.test$cvindex <- cv_index
      df.test$sg1 <- sg1
      df.test$sg2 <- sg2
    } else {
      df.test <- x.test
      df.test$cvindex <- cv_index
      df.test$sg1 <- NA
      df.test$sg2 <- NA
      df.test$treat.recommend <- 1.0
    }

    # This is the line from the actual function — does it error?
    resCV_list[[cv_index]] <- df.test[, c(id.name, outcome.name, event.name, treat.name,
                                          "treat.recommend", "treat.recommend.original",
                                          "cvindex", "sg1", "sg2")]

    cat("Fold", cv_index, ": OK -", ncol(resCV_list[[cv_index]]), "cols,",
        nrow(resCV_list[[cv_index]]), "rows\n")
  }

  resCV <- do.call(rbind, resCV_list)
  cat("rbind OK:", nrow(resCV), "rows\n")
}

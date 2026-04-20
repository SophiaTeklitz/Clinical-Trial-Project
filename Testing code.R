simRandomisation <- function(n)
{
  block  <- rep(seq(1:n), each = 4, length.out = n)
  trt    <- rep(0:1, length.out = n)
  random <- runif(n)
  
  data         <- data.frame(block, trt, random)
  data         <- data[order(data$block, data$random),]
  data$obs_no  <- 1:n
  data         <- data[, c('obs_no', 'trt')]
  
  return(data)
}

simAccrual <- function(n, accrual_period)
{
  # Recruitment times uniform over accrual period (in years)
  accrual_time <- sort(runif(n, min = 0, max = accrual_period))
  
  return(accrual_time)
}


simTrialData <- function(n, accrual_period, followup_period,
                         shape, scale_control, scale_treatment)
{
  max_time <- accrual_period + followup_period
  
  data             <- simRandomisation(n)
  data$accrual_time <- simAccrual(n, accrual_period)
  
  # Assign scale parameter based on treatment arm
  scale_vec <- ifelse(data$trt == 0, scale_control, scale_treatment)
  
  # Generate Weibull event times
  # rweibull uses shape and scale where S(t) = exp(-(t/scale)^shape)
  data$event_time <- rweibull(n, shape = shape, scale = scale_vec)
  
  # Calendar time of event = accrual time + event time
  data$event_calendar <- data$accrual_time + data$event_time
  
  # Administrative censoring: each participant followed until max_time
  data$censor_calendar <- max_time
  
  # Observed time (from recruitment) and event indicator
  # A participant is censored if their event occurs after study closes
  data$observed_time <- pmin(data$event_time,
                             data$censor_calendar - data$accrual_time)
  
  data$event <- as.integer(data$event_calendar <= data$censor_calendar)
  
  return(data)
}

simInterimData <- function(data, events_at_interims, accrual_period, followup_period)
{
  max_time <- accrual_period + followup_period
  
  # Order by calendar time of event/censoring to find when events accumulate
  data <- data[order(data$event_calendar), ]
  
  # Cumulative events in calendar time order
  data$cum_events <- cumsum(data$event)
  
  for(i in seq_along(events_at_interims)){
    
    interim_ind_name  <- paste0("interim_ind_", i)
    event_name        <- paste0("event_interim_", i)
    obs_time_name     <- paste0("obs_time_interim_", i)
    interim_time_name <- paste0("interim_cal_time_", i)
    
    # Find the calendar time when the i-th interim event count is reached
    rows_with_event <- which(data$event == 1 &
                               data$cum_events == events_at_interims[i])
    
    if(length(rows_with_event) == 0){
      warning(paste("Interim", i, "event count not reached in this simulation"))
      next
    }
    
    interim_cal_time <- data$event_calendar[min(rows_with_event)]
    
    # Participants enrolled by the interim calendar time
    enrolled <- data$accrual_time <= interim_cal_time
    
    # For enrolled participants, observed time is capped at interim calendar time
    obs_time_interim <- ifelse(
      enrolled,
      pmin(data$event_time, interim_cal_time - data$accrual_time),
      NA
    )
    
    # Event status at interim: event occurred before interim calendar time
    event_interim <- ifelse(
      enrolled,
      as.integer(data$event_calendar <= interim_cal_time),
      NA
    )
    
    data[[interim_time_name]] <- interim_cal_time
    data[[interim_ind_name]]  <- sum(enrolled)
    data[[obs_time_name]]     <- obs_time_interim
    data[[event_name]]        <- event_interim
  }
  
  return(data)
}

analyseData <- function(data, alpha_interims, alpha_final, events_at_interims)
{
  library(survival)
  
  n_interims      <- length(alpha_interims)
  interim_results <- list()
  stopped_early   <- FALSE
  
  for(i in 1:n_interims){
    
    event_col    <- paste0("event_interim_", i)
    obs_time_col <- paste0("obs_time_interim_", i)
    ind_col      <- paste0("interim_ind_", i)
    cal_time_col <- paste0("interim_cal_time_", i)
    
    if(!event_col %in% names(data)) next
    
    interim_data <- data[!is.na(data[[event_col]]), ]
    
    # Fit Cox model at interim
    surv_obj  <- Surv(time  = interim_data[[obs_time_col]],
                      event = interim_data[[event_col]])
    log_rank <- survdiff(surv_obj ~ trt, data = interim_data)
    
    interim_p    <- log_rank[["pvalue"]]
    interim_stop <- as.integer(interim_p < alpha_interims[i])
    
    interim_results[[i]] <- data.frame(
      interim          = i,
      interim_cal_time = unique(data[[cal_time_col]]),
      interim_n        = unique(data[[ind_col]]),
      interim_events   = events_at_interims[i],
      interim_p        = interim_p,
      interim_stop     = interim_stop
    )
    
    if(interim_stop == 1){
      stopped_early <- TRUE
      sample_size   <- unique(data[[ind_col]])
      break
    }
  }
  
  interim_summary <- do.call(rbind, interim_results)
  
  # Final analysis

    sample_size <- nrow(data)
    
    surv_final  <- Surv(time = data$observed_time, event = data$event)
    log_rank_final   <- survdiff(surv_final ~ trt, data = data)
    
    final_p     <- log_rank_final[["pvalue"]]
    final_stop  <- as.integer(final_p < alpha_final)

  stop     <- ifelse(stopped_early, 1, final_stop)
  flipflop <- ifelse(stopped_early & final_stop == 0, 1, 0)
  
  results <- list(
    interim_summary = interim_summary,
    final_p         = final_p,
    final_stop      = final_stop,
    sample_size     = sample_size,
    stop            = stop,
    flipflop        = flipflop
  )
  
  return(results)
}

runTrial <- function(n, accrual_period, followup_period, shape,
                     scale_control, scale_treatment,
                     events_at_interims, alpha_interims, alpha_final)
{
  data    <- simTrialData(n, accrual_period, followup_period,
                          shape, scale_control, scale_treatment)
  
  data    <- simInterimData(data, events_at_interims,
                            accrual_period, followup_period)
  
  results <- analyseData(data, alpha_interims, alpha_final, events_at_interims)
  
  return(list(data = data, results = results))
}


runMultipleTrials <- function(simno, seed, n, accrual_period, followup_period,
                              shape, scale_control, scale_treatment,
                              events_at_interims, alpha_interims, alpha_final)
{
  seeds <- seed + seq(1:simno)
  
  multiple_trials <- lapply(1:simno, function(x) {
    set.seed(seeds[x])
    y <- runTrial(n, accrual_period, followup_period,
                  shape, scale_control, scale_treatment,
                  events_at_interims, alpha_interims, alpha_final)
    return(y)
  })
  
  # Extract scalar results for summary
  results_all <- do.call(rbind, lapply(multiple_trials, function(x){
    data.frame(
      interim1_p = x[["results"]][["interim_summary"]][["interim_p"]][1],
      interim1_stop = x[["results"]][["interim_summary"]][["interim_stop"]][1],
      interim2_p = x[["results"]][["interim_summary"]][["interim_p"]][2],
      interim2_stop = x[["results"]][["interim_summary"]][["interim_stop"]][2],
      interim3_p = x[["results"]][["interim_summary"]][["interim_p"]][3],
      interim3_stop = x[["results"]][["interim_summary"]][["interim_stop"]][3],
      final_p     = x$results$final_p,
      final_stop  = x$results$final_stop,
      sample_size = x$results$sample_size,
      stop        = x$results$stop,
      flipflop    = x$results$flipflop
    )
  }))
  
  results_summary <- apply(results_all, 2, summary)
  
  return(list(
    results_all     = results_all,
    results_summary = results_summary,
    seeds           = seeds
  ))
}

seed               <- 48376491
n                  <- 2334
events_at_interims <- c(78, 116, 155)   # two interims
alpha_interims     <- c(0.004, 0.009, 0.016) ## the same cutoffs used by the study
alpha_final        <- 0.0414

# Weibull shape parameter
shape <- 0.6

# 5-year disease free survival in control group
S0_5yr <- 0.901

# Hazard ratio (treatment vs control)
hr <- 0.655

# Accrual period (years)
accrual_period <- 5

# Follow-up period after accrual (years)
followup_period <- 3

# Maximum study time
max_time <- accrual_period + followup_period

# Derive Weibull scale parameter for control group from S(t) = exp(-(t/scale)^shape)
# S(5) = 0.901 => scale = 5 / (-log(0.901))^(1/shape)
scale_control <- 5 / (-log(S0_5yr))^(1/shape)

# Scale parameter for treatment group
# Under proportional hazards with Weibull, scaling the hazard by HR
# means scale_treatment = scale_control * HR^(-1/shape)
scale_treatment <- scale_control * hr^(-1/shape)

cat("Scale parameter (control):  ", scale_control, "\n")
cat("Scale parameter (treatment):", scale_treatment, "\n")


# Single trial
set.seed(seed)
results_single <- runTrial(n, accrual_period, followup_period,
                           shape, scale_control, scale_treatment,
                           events_at_interims, alpha_interims, alpha_final)

# Multiple trials
results_multi <- runMultipleTrials(
  simno = 5000, seed = seed, n = n,
  accrual_period = accrual_period, followup_period = followup_period,
  shape = shape, scale_control = scale_control,
  scale_treatment = scale_treatment,
  events_at_interims = events_at_interims,
  alpha_interims = alpha_interims,
  alpha_final = alpha_final
)

results_multi$results_summary

## Using O'Brien Flemming cut offs

alpha_interims     <- c(0.00062, 0.00419, 0.01121) 
alpha_final        <- 0.02058
set.seed(seed)
results_single <- runTrial(n, accrual_period, followup_period,
                           shape, scale_control, scale_treatment,
                           events_at_interims, alpha_interims, alpha_final)

# Multiple trials
results_multi <- runMultipleTrials(
  simno = 5000, seed = seed, n = n,
  accrual_period = accrual_period, followup_period = followup_period,
  shape = shape, scale_control = scale_control,
  scale_treatment = scale_treatment,
  events_at_interims = events_at_interims,
  alpha_interims = alpha_interims,
  alpha_final = alpha_final
)

results_multi$results_summary

## Pocock cutoffs

alpha_interims     <- c(0.0097, 0.0097, 0.0097) ## the same cutoffs used by the study
alpha_final        <- 0.0097
set.seed(seed)
results_single <- runTrial(n, accrual_period, followup_period,
                           shape, scale_control, scale_treatment,
                           events_at_interims, alpha_interims, alpha_final)

# Multiple trials
results_multi <- runMultipleTrials(
  simno = 5000, seed = seed, n = n,
  accrual_period = accrual_period, followup_period = followup_period,
  shape = shape, scale_control = scale_control,
  scale_treatment = scale_treatment,
  events_at_interims = events_at_interims,
  alpha_interims = alpha_interims,
  alpha_final = alpha_final
)

results_multi$results_summary



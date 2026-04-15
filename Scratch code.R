
library(survival)

simRandomisation <- function(n)
{
  # A sequence indicating block.
  block <- rep(seq(1:n), each = 4, length.out = n)
  
  # A sequence indicating treatment.
  trt <- rep(0:1, length.out = n)
  
  # A random number on unit interval.
  random <- runif(n)
  
  # Create a data frame
  data <- data.frame(block, trt, random)
  
  # Order by block and then by random to create block randomised treatment
  data <- data[order(data$block, data$random),]
  
  data$obs_no <- 1:n
  
  data <- data[,c('obs_no', 'trt')]
  
  return(data)
}

simAccrual <- function(n, recruit_period)
{
  # Generate recruitment times: Simulate trial-time that patient enters the trial.
  # Adding 0.5 ensures the recruitment times are greater than day 1 when rounded.
  accrual_time <- round(runif(n) * recruit_period + 0.5)
  
  accrual_time <- sort(accrual_time)
  
  return(accrual_time)
}

simTrialData <- function(n, recruit_period, p)
{
  # Simulate random allocation.
  data <- simRandomisation(n)
  
  # Simulate recruitment times
  data$accrual_time <- simAccrual(n, recruit_period)
  
  # Simulate events from binomial distribution with respective probabilities
  # of events for control and treatment arms
  data$event <- rweibull(n, shape = lambda[data$trt+1], scale = 0.6)
  
  # Return the simulated trial data.
  return(data)
}

simInterimData <- function(data, events_at_interim)
{
  # Obtain the cumulative number of events
  data$cum_events <- cumsum(data$event)
  
  # Loop over each interim timepoint
  for(i in seq_along(events_at_interims)){
    
    interim_name <- paste0("interim_ind_", i)
    event_name   <- paste0("event_interim_", i)
    
    # observation number at which this interim occurs
    interim_ind <- data$obs_no[min(which(data$cum_events == events_at_interims[i]))]
    
    data[[interim_name]] <- interim_ind
    
    # Events available at this interim (NA for participants not yet recruited)
    data[[event_name]] <- ifelse(data$obs_no <= interim_ind, data$event, NA)
  }
  
  return(data)
}

analyseData <- function(data, alpha_interim, alpha_final)
{
  nevents0 <- sum(data$event[data$trt == 0])
  nevents1 <- sum(data$event[data$trt == 1])
  pevents0 <- nevents0/sum(data$trt == 0)
  pevents1 <- nevents1/sum(data$trt == 1)
  
  ###############################################
  # Interim analysis: logistic regression
  data$event_interim <- factor(data$event_interim)
  modellogit_int <- glm(event_interim ~ trt, data = data, family = "binomial")
  conf_int <- confint(modellogit_int)
  
  # Results at interim from the summary object
  interim_or  <- exp(coef(modellogit_int)["trt"])
  interim_lci <- exp(conf_int["trt", "2.5 %"])
  interim_uci <- exp(conf_int["trt", "97.5 %"])
  interim_p   <- coef(summary(modellogit_int))["trt", "Pr(>|z|)"]
  
  # A stop for trial success at interim?
  interim_stop <- ifelse(interim_p < alpha_interim, 1, 0)
  
  # The proportion of events if the trial stopped at the interim
  if(interim_stop == 1){
    nevents0 <- sum(data$trt == 0 & !is.na(data$event_interim))
    nevents1 <- sum(data$trt == 1 & !is.na(data$event_interim))
    pevents0 <- nevents0/sum(data$trt == 0 & !is.na(data$event_interim))
    pevents1 <- nevents1/sum(data$trt == 1 & !is.na(data$event_interim))
  }
  
  ###############################################
  # Final analysis: logistic regression
  data$event <- factor(data$event)
  modellogit <- glm(event ~ trt, data = data, family = "binomial")
  conf       <- confint(modellogit)
  
  # Results at final analysis from the summary object
  final_or  <- exp(coef(modellogit)["trt"])
  final_lci <- exp(conf["trt", "2.5 %"])
  final_uci <- exp(conf["trt", "97.5 %"])
  final_p   <- coef(summary(modellogit))["trt", "Pr(>|z|)"]
  
  final_stop <- ifelse(final_p < alpha_final, 1, 0)
  
  # Whether the trial is conclusive to determine type I error or Power depending on the scenario
  stop <- ifelse(interim_stop == 1, interim_stop, final_stop)
  
  # Sample size used to determine average sample size
  if(interim_stop == 1){
    sample_size <- unique(data$interim_ind)
  } else {
    sample_size <- nrow(data)
  }
  
  # To determine trial flip-flop probability
  flipflop <- ifelse(interim_stop == 1 & final_stop == 0, 1, 0)
  
  # results
  results <- data.frame(nevents0, nevents1,
                        pevents0, pevents1,
                        sample_size,
                        interim_time = unique(data$interim_ind),
                        interim_or, interim_lci, interim_uci, interim_p,
                        interim_stop,
                        final_or, final_lci, final_uci, final_p,
                        final_stop, stop, flipflop)
  return(results)
}

analyseData <- function(data, alpha_interims, alpha_final)
{
  # alpha_interims is now a vector matching length of events_at_interims
  # e.g. c(0.001, 0.005) for two interims
  
  n_interims <- length(alpha_interims)
  interim_results <- list()
  
  for(i in 1:n_interims){
    
    event_col <- paste0("event_interim_", i)
    ind_col   <- paste0("interim_ind_", i)
    
    interim_data <- data[!is.na(data[[event_col]]), ]
    interim_data$event_interim <- interim_data[[event_col]]
    interim_data$event_interim <- factor(interim_data$event_interim)
    
    model     <- glm(event_interim ~ trt, data = interim_data, family = "binomial")
    conf      <- confint(model)
    interim_p <- coef(summary(model))["trt", "Pr(>|z|)"]
    interim_or  <- exp(coef(model)["trt"])
    interim_lci <- exp(conf["trt", "2.5 %"])
    interim_uci <- exp(conf["trt", "97.5 %"])
    interim_stop <- ifelse(interim_p < alpha_interims[i], 1, 0)
    
    interim_results[[i]] <- data.frame(
      interim        = i,
      interim_time   = unique(data[[ind_col]]),
      interim_or     = interim_or,
      interim_lci    = interim_lci,
      interim_uci    = interim_uci,
      interim_p      = interim_p,
      interim_stop   = interim_stop
    )
  }
  
  # Combine interim results
  interim_summary <- do.call(rbind, interim_results)
  
  # Final analysis (only runs if not stopped early)
    sample_size   <- nrow(data)
    data$event    <- factor(data$event)
    model_final   <- glm(event ~ trt, data = data, family = "binomial")
    conf_final    <- confint(model_final)
    final_p       <- coef(summary(model_final))["trt", "Pr(>|z|)"]
    final_or      <- exp(coef(model_final)["trt"])
    final_lci     <- exp(conf_final["trt", "2.5 %"])
    final_uci     <- exp(conf_final["trt", "97.5 %"])
    final_stop    <- ifelse(final_p < alpha_final, 1, 0)
    # Use interim model results as final if stopped early
    last_interim  <- interim_summary[nrow(interim_summary), ]
    final_or      <- last_interim$interim_or
    final_lci     <- last_interim$interim_lci
    final_uci     <- last_interim$interim_uci
    final_p       <- last_interim$interim_p
    final_stop    <- last_interim$interim_stop
  
  stop     <- ifelse(stopped_early, 1, final_stop)
  flipflop <- ifelse(stopped_early & final_stop == 0, 1, 0)
  
  results <- list(
    interim_summary = interim_summary,
    final_or        = final_or,
    final_lci       = final_lci,
    final_uci       = final_uci,
    final_p         = final_p,
    final_stop      = final_stop,
    sample_size     = sample_size,
    stop            = stop,
    flipflop        = flipflop
  )
  
  return(results)
}


## Parameters from the example

# The random seed to ensure reproducibility
seed          <- 48376491

# Recruitment period = Days in 2.5 years. 690 patients in 3 years (365.25*3days)
recruit_period <- 365.25*3*584/690

# Maximum trial sample size.
n             <- 584

# The number of events at the interim: half recruitment (584/2 = 292; 292*(.1+.04)/2=20 events)
events_at_interim <- ## number of events at the 20%, 40%, 60%, and 80% of expected events

# event probabilities
# Event probability in Na77 at 72 hours
lambda0 <- -log(0.901)/(5^0.6)

# The event probability for Na140 arm
lambda1    <- -log(0.934)/(5^0.6)

# vector of event probabilities
lambda     <- c(p0, p1)

## scratch work solving for survival time

hist(((-log(runif(1000))/(lambda0*0.655)))^(1/0.6))

hist(((-log(runif(1000))/(lambda0)))^(1/0.6))

# Decision thresholds/boundaries (alpha)
# At final analysis
alpha_final    <- 0.045

# At the interim
alpha_interim  <- 0.005

runTrial <- function(n, recruit_period, p, events_at_interim, alpha_interim, alpha_final)
{
  # Step 1: simulate the trial data
  data    <- simTrialData(n, recruit_period, p)
  
  # Step 2: create interim data
  data    <- simInterimData(data, events_at_interim)
  
  # Step 3: analyse the data
  results <- analyseData(data, alpha_interim, alpha_final)
  
  return(list(data = data, results = results))
}

set.seed(seed)
results_single_trial <- runTrial(n, recruit_period, p,
                                 events_at_interim,
                                 alpha_interim,
                                 alpha_final)
results_single_trial$results

runMultipleTrials <- function(simno, seed, n, recruit_period, p,
                              events_at_interim, alpha_interim, alpha_final)
{
  # Random seeds; length should be equal to the simno
  seeds <- seed + seq(1:simno)
  
  # Simulate datasets
  multiple_trials <- lapply(1:simno, function(x) {
    set.seed(seeds[x])
    y <- runTrial(n, recruit_period, p, events_at_interim, alpha_interim,
                  alpha_final)
    return(y)
  })
  
  # Summarise results
  results     <- lapply(multiple_trials, function(x) return(x$results))
  results_all <- do.call(rbind, results)
  saveRDS(multiple_trials, file = 'results_multitrials.rds')
  
  # Remove any simulations with non-estimable CIs and calculate the summary
  x <- which(apply(results_all, 1, function(x) any(is.na(x))))
  
  if(length(x) > 0){
    results_summary <- apply(results_all[-x,], 2, summary)
  } else {
    results_summary <- apply(results_all, 2, summary)
  }
  
  return(list(results_all = results_all,
              results_summary = results_summary,
              seeds = seeds))
}

simno <- 5000

results_multi_trials <- runMultipleTrials(simno, seed, n, recruit_period,
                                          p, events_at_interim,
                                          alpha_interim, alpha_final)
results_multi_trials$results_summary














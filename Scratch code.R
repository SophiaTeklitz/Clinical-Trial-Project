
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
  
  # observation number at which interim occurs
  data$interim_ind <- data$obs_no[min(which(data$cum_events == events_at_interim))]
  
  # Events at interim
  data$event_interim <- with(data, ifelse(obs_no <= interim_ind, event, NA))
  
  # you can also extract the interim data set and output it separately as below
  # data_interim <- subset(data, !is.na(data$event_interim))
  # return(data_interim)
  
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


## Parameters from the example

# The random seed to ensure reproducibility
seed          <- 48376491

# Recruitment period = Days in 2.5 years. 690 patients in 3 years (365.25*3days)
recruit_period <- 365.25*3*584/690

# Maximum trial sample size.
n             <- 584

# The number of events at the interim: half recruitment (584/2 = 292; 292*(.1+.04)/2=20 events)
events_at_interim <- 20

# event probabilities
# Event probability in Na77 at 72 hours
lambda0    <- 5/(-log(0.901)^(1/0.6))

# The event probability for Na140 arm
lambda1    <- 5/(-log(0.934)^(1/0.6))

# vector of event probabilities
lambda     <- c(p0, p1)

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














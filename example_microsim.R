# =====================================================
# example estimation of undiagnosed cases for EQUITY-MS
# -----------------------------------------------------

# inputs

days <- 1095 # 1095 = 3 years
mean_cases_per_day <- 0.5
mean_delay <- c(300, 200, 100) # delay in days for three groups
group_sizes <- c(0.5, 0.3, 0.2)

# generate data (in the real analysis this will be taken from observed patient data)

cases <- rpois(days, mean_cases_per_day) # cases on each day
n <- sum(cases) # total number of cases
group <- sample(1:3, n, replace = T, prob = group_sizes)
duration_targets <- mean_delay[group]
durations <- rnbinom(n, size = 0.5, mu = duration_targets)

# describe duration by group

mean_duration_by_group <- tapply(durations, INDEX = group, FUN = mean)
tapply(durations, INDEX = group, FUN = quantile, probs = c(0.25, 0.5, 0.75))
par(mfrow = c(1, 3))
lapply(split(durations, group), hist, xlim = c(0, 2000), xlab = "Days", breaks = 20)

# typical number in the community (simple analysis using averages)
# number in community = rates of cases * duration of wait
# e.g if you diagnose 1/day, and the wait is 1 day, there is always 1 in the community
# e.g.if you diagnose 5/day, and the wait is 5 days, there are always 25 in the community

(n/days) * mean(durations) # total in the community
rate_by_group <- table(group) / days
rate_by_group * mean_duration_by_group # total by group

# microsim: ungrouped
# incidence adds to the pool, diagnosis removes them

sim_days <- 20*365
sim_cases <- rpois(sim_days, mean_cases_per_day)
sim_n <- sum(sim_cases)
incident_day <- rep(seq_len(sim_days), sim_cases)
diagnosis_day <- incident_day + sample(durations, sim_n, replace = T)
day_index <- min(incident_day):sim_days
waiting_counts <- sapply(day_index, function(t) {
  sum(incident_day <= t & diagnosis_day > t)
})
results <- data.frame(day = day_index, waiting = waiting_counts)
dev.off()
plot(results$waiting, type = 'l')
mean(tail(results$waiting, 365))

# key additions for main analysis:
# 1. do simulation by group to show equity
# 2. relax assumption that reported delays represent delays in the community
# (eg. some cases never get diagnosed)
# 3. consider modelling diagnosis as a stochastic risk conditional on equity-based characteristics
# 4. repeat simulation multiple times to understand uncertainty

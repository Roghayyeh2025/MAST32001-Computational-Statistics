
### TASK 2 ###

############
## PART A ##
############
# LATENT-VARIABLE REPRESENTATION AND GIBBS SAMPLER


# 1. STATISTICAL MODEL
# We consider the two-component Gaussian mixture
#   Y_i | Z_i = 1, mu_1, sigma^2 ~ N(mu_1, sigma^2)
#   Y_i | Z_i = 2, mu_2, sigma^2 ~ N(mu_2, sigma^2)
#
#   P(Z_i = 1) = pi
#   P(Z_i = 2) = 1 - pi
#
# Inverse-Gamma parameterisation used in this script:
#
#   p(sigma^2) proportional to
#     (sigma^2)^(-(alpha0 + 1)) * exp(-beta0 / sigma^2)
#
# Therefore, if X ~ InvGamma(alpha, beta), generate it as:
#
#   X = 1 / rgamma(1, shape = alpha, rate = beta)
#
# The latent allocations Z_i make Gibbs sampling possible because the
# full conditional distributions have standard forms.

# 2. PRIOR CHOICES

# The data in the example below are standardised, so these priors are
# specified on the standardised scale.
#
#   mu_k ~ N(0, 5^2)
#   pi ~ Beta(1, 1)
#   sigma^2 ~ InvGamma(2.5, 1.5)
#
# Beta(1, 1) is uniform on [0, 1].
# Under the stated Inverse-Gamma parameterisation, the prior mean of
# sigma^2 is beta0 / (alpha0 - 1) = 1.5 / 1.5 = 1.

m0 <- 0
s0_sq <- 25
a <- 1
b <- 1
alpha0 <- 2.5
beta0 <- 1.5


# 3. HELPER FUNCTION: STANDARDISE DATA


standardise_data <- function(y) {
  y_mean <- mean(y)
  y_sd <- sd(y)
  
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("The data must have a positive, finite standard deviation.")
  }
  
  list(
    y_std = (y - y_mean) / y_sd,
    mean = y_mean,
    sd = y_sd
  )
}


# 4. FULL CONDITIONAL DISTRIBUTIONS AND THEIR DERIVATION


# The full conditionals follow from:
# posterior kernel proportional to likelihood * prior.
# Terms that do not depend on the parameter being updated are omitted.
#
# 4.1 Latent allocation Z_i
#
# Define:
#   p1_i = P(Z_i = 1 | y_i, parameters)
#
# Then:
#
#   p1_i =
#     pi * phi(y_i; mu_1, sigma^2) /
#     [pi * phi(y_i; mu_1, sigma^2)
#       + (1 - pi) * phi(y_i; mu_2, sigma^2)]
#
#   Z_i | rest ~ Bernoulli(p1_i), represented by labels 1 and 2.
#
# Derivation for Z_i:
# Conditional on the parameters, Bayes' rule gives the probability of
# each allocation as its mixture weight times its Gaussian density,
# divided by the sum of those two weighted densities.
#
# 4.2 Component means
#
# For component k, combine the likelihood kernel
#   exp[-sum_{i:Z_i=k}(y_i - mu_k)^2 / (2*sigma^2)]
# with the Normal prior kernel
#   exp[-(mu_k - m0)^2 / (2*s0_sq)].
# Completing the square in mu_k gives the Normal full conditional below.
#
# Let n_k be the number of observations allocated to component k, and
# S_k be the sum of those observations. Then:
#
#   mu_k | rest ~ N(m_k_star, v_k_star)
#
#   v_k_star = 1 / (n_k / sigma^2 + 1 / s0_sq)
#
#   m_k_star = v_k_star * (S_k / sigma^2 + m0 / s0_sq)
#
# 4.3 Mixture probability
#
# Given the allocations, the likelihood for pi is proportional to
# pi^n1 * (1 - pi)^n2. Multiplying by the Beta(a, b) prior gives the
# Beta(a + n1, b + n2) full conditional.
#
# Let n_1 = number of observations with Z_i = 1 and
#     n_2 = number of observations with Z_i = 2.
#
#   pi | rest ~ Beta(a + n_1, b + n_2)
#
# 4.4 Common variance
#
# Conditional on the allocations and component means, the Gaussian
# likelihood kernel in sigma^2 is proportional to
# (sigma^2)^(-n/2) * exp[-sum_i (y_i - mu_Zi)^2 / (2*sigma^2)].
# Multiplying this by the specified Inverse-Gamma prior kernel adds
# n/2 to the shape parameter and half the residual sum of squares to
# the second parameter.
#
# Let mu_Zi equal mu_1 when Z_i = 1 and mu_2 when Z_i = 2.
# Under the prior specified above:
#
#   sigma^2 | rest ~ InvGamma(
#     alpha0 + n / 2,
#     beta0 + 0.5 * sum_i (y_i - mu_Zi)^2
#   )

# Now, I have drived the full conditional distributions.
# The following section implements the Gibbs sampler,
# stores posterior draws, and assesses Markov chain behaviour.

# 5. GIBBS SAMPLER IMPLEMENTED FROM SCRATCH


gibbs_mix <- function(
    y,
    n_iter = 12000,
    burn = 4000,
    thin = 4,
    init = NULL,
    m0 = 0,
    s0_sq = 25,
    a = 1,
    b = 1,
    alpha0 = 2.5,
    beta0 = 1.5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.numeric(y) || length(y) < 2 || any(!is.finite(y))) {
    stop("y must be a finite numeric vector with at least two observations.")
  }
  if (n_iter <= burn) stop("n_iter must be larger than burn.")
  if (thin < 1) stop("thin must be at least 1.")
  if (s0_sq <= 0 || a <= 0 || b <= 0 ||
      alpha0 <= 0 || beta0 <= 0) {
    stop("Prior variance and distribution parameters must be positive.")
  }
  
  n <- length(y)
  
# Initial values
if (is.null(init)) {
    z <- sample(1:2, n, replace = TRUE)
    
# Ensure both components have at least one observation initially.
    if (length(unique(z)) < 2) {
      z[1] <- 1L
      z[2] <- 2L
    }
    
    mu1 <- mean(y[z == 1])
    mu2 <- mean(y[z == 2])
    
    if (mu1 == mu2) {
      mu1 <- mu1 - 0.5
      mu2 <- mu2 + 0.5
    }
    
    pi <- min(max(mean(z == 1), 0.05), 0.95)
    sigma2 <- var(y)
  } else {
    required <- c("mu1", "mu2", "pi", "sigma2")
    if (!all(required %in% names(init))) {
      stop("init must contain mu1, mu2, pi, and sigma2.")
    }
    
    mu1 <- init$mu1
    mu2 <- init$mu2
    pi <- init$pi
    sigma2 <- init$sigma2
    
    if (!is.null(init$z)) {
      z <- init$z
      if (length(z) != n || any(!z %in% c(1, 2))) {
        stop("init$z must have length length(y) and contain only 1 or 2.")
      }
    } else {
      z <- sample(1:2, n, replace = TRUE)
    }
  }
  
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("The initial sigma2 must be positive and finite.")
  }
  if (!is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("The initial pi must be strictly between 0 and 1.")
  }
  
  keep <- seq(from = burn + 1, to = n_iter, by = thin)
  n_keep <- length(keep)
  
  # Store parameter draws and the latent allocations for posterior
  # inference and examination of Markov-chain behaviour.
  out <- list(
    mu1 = numeric(n_keep),
    mu2 = numeric(n_keep),
    pi = numeric(n_keep),
    sigma2 = numeric(n_keep),
    z = matrix(NA_integer_, nrow = n_keep, ncol = n)
  )
  
  keep_counter <- 0L
  
  for (iter in seq_len(n_iter)) {
    
    # STEP 1: Update each latent allocation Z_i.
    # Work on the log scale to reduce numerical underflow.
    log_p1 <- log(pi) +
      dnorm(y, mean = mu1, sd = sqrt(sigma2), log = TRUE)
    
    log_p2 <- log1p(-pi) +
      dnorm(y, mean = mu2, sd = sqrt(sigma2), log = TRUE)
    
    p1 <- plogis(log_p1 - log_p2)
    z <- ifelse(runif(n) < p1, 1L, 2L)
    
    # STEP 2: Update mu_1.
    n1 <- sum(z == 1)
    v1 <- 1 / (n1 / sigma2 + 1 / s0_sq)
    m1 <- v1 * (sum(y[z == 1]) / sigma2 + m0 / s0_sq)
    mu1 <- rnorm(1, mean = m1, sd = sqrt(v1))
    
    # STEP 3: Update mu_2.
    n2 <- sum(z == 2)
    v2 <- 1 / (n2 / sigma2 + 1 / s0_sq)
    m2 <- v2 * (sum(y[z == 2]) / sigma2 + m0 / s0_sq)
    mu2 <- rnorm(1, mean = m2, sd = sqrt(v2))
    
    # STEP 4: Update pi.
    pi <- rbeta(1, shape1 = a + n1, shape2 = b + n2)
    
    # STEP 5: Update sigma^2.
    residuals <- ifelse(z == 1, y - mu1, y - mu2)
    alpha_post <- alpha0 + n / 2
    beta_post <- beta0 + 0.5 * sum(residuals^2)
    
    sigma2 <- 1 / rgamma(
      1,
      shape = alpha_post,
      rate = beta_post
    )
    
    # Save draws after burn-in, applying the requested thinning.
    if (iter %in% keep) {
      keep_counter <- keep_counter + 1L
      out$mu1[keep_counter] <- mu1
      out$mu2[keep_counter] <- mu2
      out$pi[keep_counter] <- pi
      out$sigma2[keep_counter] <- sigma2
      out$z[keep_counter, ] <- z
    }
  }
  
  out$n_iter <- n_iter
  out$burn <- burn
  out$thin <- thin
  out
}



# 6. BASIC CHECK OF THE GIBBS SAMPLER


# A simple simulated mixture is used only to check that the sampler
# runs and stores the expected quantities.
# n_iter = 6000 and burn = 2000 are initial pilot settings, not a claim
# that these values guarantee convergence. We must increase the run length if
# chains have not mixed adequately or effective sample sizes are too low.
# Burn-in is chosen by inspecting the early transient behaviour, thinning
# is used here only to reduce the number of stored draws, not to create
# convergence. 

set.seed(1001)
y_check <- c(
  rnorm(100, mean = -2, sd = 1),
  rnorm(100, mean = 2, sd = 1)
)

y_check_info <- standardise_data(y_check)
y_check_std <- y_check_info$y_std

check_init <- list(
  mu1 = -1,
  mu2 = 1,
  pi = 0.5,
  sigma2 = 1,
  z = sample(1:2, length(y_check_std), replace = TRUE)
)

check_draws <- gibbs_mix(
  y = y_check_std,
  n_iter = 6000,
  burn = 2000,
  thin = 4,
  init = check_init,
  m0 = m0,
  s0_sq = s0_sq,
  a = a,
  b = b,
  alpha0 = alpha0,
  beta0 = beta0,
  seed = 1002
)

# Basic posterior summary for the stored parameter draws.
posterior_summary <- function(x) {
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  c(
    Mean = mean(x),
    SD = sd(x),
    Q2.5 = unname(qs[1]),
    Median = unname(qs[2]),
    Q97.5 = unname(qs[3])
  )
}

summary_table <- rbind(
  mu1 = posterior_summary(check_draws$mu1),
  mu2 = posterior_summary(check_draws$mu2),
  pi = posterior_summary(check_draws$pi),
  sigma2 = posterior_summary(check_draws$sigma2)
)


# Part A is complete: I specified initialization and preliminary iteration
# settings.
# Part B evaluates chain behaviour to assess their adequacy.
# I will assess their adequacy using convergence diagnostics
# in Part B.

############
## PART B ##
############

#  SIMULATION SETTINGS

set.seed(20261012)
n_sim <- 300
true_pi <- 0.5
true_sigma <- 1

# Clearly separated: the mean difference is four standard deviations.

scenario_separated <- list(
  name = "Clearly separated",
  mu1 = -2,
  mu2 =  2,
  pi = true_pi,
  sigma = true_sigma
)

# Substantially overlapping: the mean difference is 1.5 standard deviations.

scenario_overlapping <- list(
  name = "Substantially overlapping",
  mu1 = -0.75,
  mu2 =  0.75,
  pi = true_pi,
  sigma = true_sigma
)

# Generate observations from the known mixture and retain the generating
# component indicators for the later allocation-recovery assessment.

simulate_mixture <- function(n, mu1, mu2, pi, sigma) {
  z_true <- ifelse(runif(n) < pi, 1L, 2L)
  y <- numeric(n)
  y[z_true == 1L] <- rnorm(sum(z_true == 1L), mean = mu1, sd = sigma)
  y[z_true == 2L] <- rnorm(sum(z_true == 2L), mean = mu2, sd = sigma)
  
  list(y = y, z_true = z_true)
}

sim_sep <- simulate_mixture(
  n_sim, scenario_separated$mu1, scenario_separated$mu2,
  scenario_separated$pi, scenario_separated$sigma
)

sim_overlap <- simulate_mixture(
  n_sim, scenario_overlapping$mu1, scenario_overlapping$mu2,
  scenario_overlapping$pi, scenario_overlapping$sigma
)

# Visualise the generated data and the true component densities.

plot_simulated_mixture <- function(sim, pars, main) {
  hist(sim$y, breaks = 25, probability = TRUE,
       col = "grey85", border = "white",
       main = main, xlab = "y", ylim = c(0, max(
         density(sim$y)$y,
         dnorm(seq(min(sim$y), max(sim$y), length.out = 300),
               pars$mu1, pars$sigma),
         dnorm(seq(min(sim$y), max(sim$y), length.out = 300),
               pars$mu2, pars$sigma)
       ) * 1.15))
  xx <- seq(min(sim$y) - 0.5, max(sim$y) + 0.5, length.out = 500)
  lines(xx, pars$pi * dnorm(xx, pars$mu1, pars$sigma), lwd = 2)
  lines(xx, (1 - pars$pi) * dnorm(xx, pars$mu2, pars$sigma), lwd = 2, lty = 2)
  lines(xx, pars$pi * dnorm(xx, pars$mu1, pars$sigma) +
          (1 - pars$pi) * dnorm(xx, pars$mu2, pars$sigma), lwd = 2, lty = 3)
  legend("topright",
         legend = c("Component 1", "Component 2", "Mixture"),
         lty = c(1, 2, 3), lwd = 2, bty = "n")
}

par(mfrow = c(1, 2))
plot_simulated_mixture(sim_sep, scenario_separated,
                       "Part B: clearly separated mixture")
plot_simulated_mixture(sim_overlap, scenario_overlapping,
                       "Part B: overlapping mixture")
par(mfrow = c(1, 1))


#  MULTIPLE CHAINS FROM DIFFERENT INITIAL VALUES
#
# The initial means are deliberately dispersed and include reversed
# component labels. This helps reveal sensitivity to initial values and
# label switching. The sampler itself is not constrained to keep mu1 < mu2.
#
# These are pilot settings. Diagnostics below should be inspected before
# treating the selected run length as adequate.

n_iter_B <- 12000
burn_B <- 4000
thin_B <- 2
n_chains_B <- 3

make_initial_values <- function(y, chain_id) {
# Different starting values for each chain. The allocation vector is
# independently randomised, and the initial variance is kept positive.
  starts <- list(
    c(-1.5,  1.5),
    c( 1.5, -1.5),
    c(-0.25, 0.25)
  )
  list(
    mu1 = starts[[chain_id]][1],
    mu2 = starts[[chain_id]][2],
    pi = c(0.35, 0.65, 0.50)[chain_id],
    sigma2 = c(1.5, 0.7, 2.0)[chain_id],
    z = sample(1:2, length(y), replace = TRUE)
  )
}

run_scenario_chains <- function(sim, scenario_name) {
  chains <- vector("list", n_chains_B)
  
  for (ch in seq_len(n_chains_B)) {
    chains[[ch]] <- gibbs_mix(
      y = sim$y,
      n_iter = n_iter_B,
      burn = burn_B,
      thin = thin_B,
      init = make_initial_values(sim$y, ch),
      m0 = m0,
      s0_sq = s0_sq,
      a = a,
      b = b,
      alpha0 = alpha0,
      beta0 = beta0,
      seed = 3100 + ch
    )
  }
  
  names(chains) <- paste0("chain", seq_len(n_chains_B))
  attr(chains, "scenario") <- scenario_name
  chains
}

chains_separated <- run_scenario_chains(sim_sep, "Clearly separated")
chains_overlapping <- run_scenario_chains(sim_overlap, "Substantially overlapping")


#  LABEL-INVARIANT SUMMARIES
#
# The likelihood is unchanged when the two component labels are exchanged.
# For comparing component means across chains, order the means within each
# saved draw. The corresponding allocation labels and mixture weight are
# relabelled in the same draw. This is a post-processing device for the
# simulation comparison, not a constraint imposed on the Gibbs sampler.

align_chain_draws <- function(chain) {
  m <- length(chain$mu1)
  mu_low <- pmin(chain$mu1, chain$mu2)
  mu_high <- pmax(chain$mu1, chain$mu2)
  component_1_is_low <- chain$mu1 <= chain$mu2
  
  pi_low <- ifelse(component_1_is_low, chain$pi, 1 - chain$pi)
  z_low <- chain$z
  for (i in seq_len(m)) {
    if (!component_1_is_low[i]) z_low[i, ] <- 3L - z_low[i, ]
  }
  
  list(
    mu_low = mu_low,
    mu_high = mu_high,
    pi_low = pi_low,
    sigma2 = chain$sigma2,
    z_low = z_low
  )
}

aligned_separated <- lapply(chains_separated, align_chain_draws)
aligned_overlapping <- lapply(chains_overlapping, align_chain_draws)


#  POSTERIOR SUMMARIES AND RECOVERY OF GENERATING VALUES

posterior_summary_B <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  c(
    Mean = mean(x),
    SD = sd(x),
    Q2.5 = unname(q[1]),
    Median = unname(q[2]),
    Q97.5 = unname(q[3])
  )
}

combine_parameter_draws <- function(aligned_chains) {
  do.call(rbind, lapply(aligned_chains, function(ch) {
    data.frame(
      mu_low = ch$mu_low,
      mu_high = ch$mu_high,
      pi_low = ch$pi_low,
      sigma2 = ch$sigma2
    )
  }))
}

draws_sep <- combine_parameter_draws(aligned_separated)
draws_overlap <- combine_parameter_draws(aligned_overlapping)

summary_separated <- rbind(
  lower_mean = posterior_summary_B(draws_sep$mu_low),
  upper_mean = posterior_summary_B(draws_sep$mu_high),
  mixing_weight_for_lower_mean = posterior_summary_B(draws_sep$pi_low),
  common_variance = posterior_summary_B(draws_sep$sigma2)
)

summary_overlapping <- rbind(
  lower_mean = posterior_summary_B(draws_overlap$mu_low),
  upper_mean = posterior_summary_B(draws_overlap$mu_high),
  mixing_weight_for_lower_mean = posterior_summary_B(draws_overlap$pi_low),
  common_variance = posterior_summary_B(draws_overlap$sigma2)
)

cat("\nPART B: posterior summaries, clearly separated mixture\n")
print(round(summary_separated, 3))
cat("\nTrue values (ordered means, weight of lower-mean component, variance):\n")
print(c(mu_low = scenario_separated$mu1,
        mu_high = scenario_separated$mu2,
        pi_low = scenario_separated$pi,
        sigma2 = scenario_separated$sigma^2))

cat("\nPART B: posterior summaries, overlapping mixture\n")
print(round(summary_overlapping, 3))
cat("\nTrue values (ordered means, weight of lower-mean component, variance):\n")
print(c(mu_low = scenario_overlapping$mu1,
        mu_high = scenario_overlapping$mu2,
        pi_low = scenario_overlapping$pi,
        sigma2 = scenario_overlapping$sigma^2))

#  TRACE PLOTS
#
# Plot ordered means and the variance. Ordered means are used here to avoid
# mistaking label exchange for a change in the underlying mixture. The
# un-ordered mu1 and mu2 traces can also be inspected to see label switching.

plot_chain_traces <- function(aligned_chains, main_prefix) {
  par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
  vars <- c("mu_low", "mu_high", "sigma2")
  for (v in vars) {
    all_y <- unlist(lapply(aligned_chains, function(ch) ch[[v]]))
    plot(aligned_chains[[1]][[v]], type = "l", col = 1,
         xlab = "Stored post-burn-in draw", ylab = v,
         main = paste(main_prefix, "-", v),
         ylim = range(all_y, finite = TRUE))
    if (length(aligned_chains) > 1) {
      for (ch in 2:length(aligned_chains)) {
        lines(aligned_chains[[ch]][[v]], col = ch)
      }
    }
  }
  legend("topright", legend = names(aligned_chains),
         col = seq_along(aligned_chains), lty = 1, bty = "n")
  par(mfrow = c(1, 1))
}

plot_chain_traces(aligned_separated, "Separated mixture")
plot_chain_traces(aligned_overlapping, "Overlapping mixture")


#  AUTOCORRELATION AND EFFECTIVE SAMPLE SIZE (ESS)

# ESS is estimated separately for each chain from its autocorrelation sequence
# using a positive-pair truncation rule, then summed across chains. The result
# is an approximate diagnostic, not a substitute for a validated MCMC package.
# MCSE for the posterior mean is approximated by posterior SD / sqrt(ESS).
# The trace plots and between-chain comparisons should be considered too.

ess_one_chain <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3 || !is.finite(sd(x)) || sd(x) == 0) return(NA_real_)
  if (is.null(max_lag)) max_lag <- min(n - 1, 1000)
  
  ac <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE,
                       demean = TRUE)$acf)[-1]
  if (!length(ac)) return(n)
  
# Sum autocorrelations in adjacent pairs and stop at the first
# non-positive pair, a simple approximation to the initial-positive-sequence
# approach used in MCMC diagnostics.
  pair_sums <- numeric(0)
  j <- 1
  while (j <= length(ac)) {
    pair_sum <- ac[j]
    if (j + 1 <= length(ac)) pair_sum <- pair_sum + ac[j + 1]
    if (!is.finite(pair_sum) || pair_sum <= 0) break
    pair_sums <- c(pair_sums, pair_sum)
    j <- j + 2
  }
  
  tau <- max(1, 1 + 2 * sum(pair_sums))
  min(n, n / tau)
}

ess_multi_chain <- function(aligned_chains, variable) {
  sum(vapply(aligned_chains, function(ch) ess_one_chain(ch[[variable]]),
             numeric(1)), na.rm = TRUE)
}

diagnostic_ess_table <- function(aligned_chains) {
  variables <- c("mu_low", "mu_high", "pi_low", "sigma2")
  all_draws <- combine_parameter_draws(aligned_chains)
  ess_values <- vapply(variables, function(v) {
    ess_multi_chain(aligned_chains, v)
  }, numeric(1))
  posterior_sd <- vapply(variables, function(v) sd(all_draws[[v]]), numeric(1))
  
  data.frame(
    parameter = variables,
    stored_draws = nrow(all_draws),
    ESS_total = unname(ess_values),
    posterior_SD = unname(posterior_sd),
    MCSE_mean_approx = unname(posterior_sd / sqrt(ess_values))
  )
}

ess_separated <- diagnostic_ess_table(aligned_separated)
ess_overlapping <- diagnostic_ess_table(aligned_overlapping)

cat("\nApproximate ESS: clearly separated mixture\n")
print(ess_separated, row.names = FALSE)
cat("\nApproximate ESS: overlapping mixture\n")
print(ess_overlapping, row.names = FALSE)

# ACF plots for representative label-invariant parameters from the first
# chain. Slow decay indicates stronger serial dependence.

par(mfrow = c(2, 2))
for (v in c("mu_low", "mu_high", "pi_low", "sigma2")) {
  acf(aligned_separated[[1]][[v]], main = paste("Separated:", v))
}
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
for (v in c("mu_low", "mu_high", "pi_low", "sigma2")) {
  acf(aligned_overlapping[[1]][[v]], main = paste("Overlapping:", v))
}
par(mfrow = c(1, 1))

# BASIC MULTI-CHAIN R-HAT
# This is the classical Gelman-Rubin potential scale reduction estimate
# calculated from equal-length chains. The reported estimate is bounded below
# by 1. Values close to 1 are reassuring,
# but no single threshold proves convergence. Interpret it together with
# trace plots, ESS, autocorrelation, and stability under longer runs.

rhat_basic <- function(chains) {
  n <- min(vapply(chains, length, integer(1)))
  m <- length(chains)
  if (m < 2 || n < 2) return(NA_real_)
  
  mat <- sapply(chains, function(x) as.numeric(x[seq_len(n)]))
  chain_means <- colMeans(mat)
  chain_vars <- apply(mat, 2, var)
  W <- mean(chain_vars)
  B <- n * var(chain_means)
  
  if (!is.finite(W) || W <= 0) return(NA_real_)
  var_hat <- ((n - 1) / n) * W + B / n
# The conventional reported potential scale reduction factor is
# bounded below by 1; small-sample estimates can otherwise fall below 1.
  max(1, sqrt(var_hat / W))
}

diagnostic_rhat_table <- function(aligned_chains) {
  variables <- c("mu_low", "mu_high", "pi_low", "sigma2")
  data.frame(
    parameter = variables,
    Rhat_classical = vapply(variables, function(v) {
      rhat_basic(lapply(aligned_chains, function(ch) ch[[v]]))
    }, numeric(1))
  )
}
rhat_separated <- diagnostic_rhat_table(aligned_separated)
rhat_overlapping <- diagnostic_rhat_table(aligned_overlapping)

cat("\nR-hat: clearly separated mixture\n")
print(rhat_separated, row.names = FALSE)

cat("\nR-hat: overlapping mixture\n")
print(rhat_overlapping, row.names = FALSE)


# RECOVERY OF LATENT ALLOCATIONS
#
# Allocation recovery is evaluated after aligning labels within each draw.
# For each observation, estimate the posterior probability of belonging to
# the lower-mean component, then classify at 0.5. Accuracy is descriptive:
# it is expected to be lower in the overlapping case because the components
# generate observations in the same region.

allocation_recovery <- function(aligned_chains, z_true) {
  z_draws <- do.call(rbind, lapply(aligned_chains, function(ch) ch$z_low))
  prob_lower <- colMeans(z_draws == 1L)
  z_hat <- ifelse(prob_lower >= 0.5, 1L, 2L)
  data.frame(
    posterior_mean_probability_lower = mean(prob_lower),
    classification_accuracy = mean(z_hat == z_true),
    mean_allocation_uncertainty = mean(2 * pmin(prob_lower, 1 - prob_lower)),
    n_observations = length(z_true)
  )
}

allocation_sep <- allocation_recovery(aligned_separated, sim_sep$z_true)
allocation_overlap <- allocation_recovery(aligned_overlapping, sim_overlap$z_true)

cat("\nLatent-allocation recovery: clearly separated mixture\n")
print(allocation_sep, row.names = FALSE)
cat("\nLatent-allocation recovery: overlapping mixture\n")
print(allocation_overlap, row.names = FALSE)


# POSTERIOR DENSITIES AGAINST GENERATING VALUES
#
# These plots compare posterior uncertainty with the known data-generating
# values. The means are ordered, and the mixture weight refers to the
# lower-mean component after relabelling.

plot_posterior_recovery <- function(draws, pars, main_prefix) {
  true_values <- c(
    mu_low = min(pars$mu1, pars$mu2),
    mu_high = max(pars$mu1, pars$mu2),
    pi_low = if (pars$mu1 <= pars$mu2) pars$pi else 1 - pars$pi,
    sigma2 = pars$sigma^2
  )
  par(mfrow = c(2, 2))
  for (v in names(true_values)) {
    plot(density(draws[[v]]), main = paste(main_prefix, "-", v),
         xlab = v, lwd = 2)
    abline(v = true_values[[v]], lty = 2, lwd = 2)
  }
  par(mfrow = c(1, 1))
}

plot_posterior_recovery(draws_sep, scenario_separated,
                        "Separated posterior")
plot_posterior_recovery(draws_overlap, scenario_overlapping,
                        "Overlapping posterior")


# EXTENSION CHECK FOR ITERATION CHOICE

# A practical check is to compare posterior summaries from the first and
# second halves of the retained draws. Large differences suggest that the
# chain may not yet be stable. This is a descriptive check and does not
# replace multi-chain diagnostics. If diagnostics are unsatisfactory,
# increase n_iter and repeat the analysis rather than relying on thinning.


# HALF-RUN STABILITY CHECK WITHIN EACH CHAIN

# Compare the first and second halves of each chain separately.
# Do not combine chains before splitting the retained draws.

half_summary_difference <- function(aligned_chains) {
  
  vars <- c("mu_low", "mu_high", "pi_low", "sigma2")
  
  results <- lapply(seq_along(aligned_chains), function(ch_id) {
    
    ch <- aligned_chains[[ch_id]]
    n <- length(ch$mu_low)
    
    if (n < 2) {
      stop("Each chain must contain at least two retained draws.")
    }
    
    half_n <- floor(n / 2)
    
    first_idx <- seq_len(half_n)
    second_idx <- seq.int(half_n + 1, n)
    
    first <- lapply(vars, function(v) {
      mean(ch[[v]][first_idx])
    })
    
    second <- lapply(vars, function(v) {
      mean(ch[[v]][second_idx])
    })
    
    first <- unlist(first)
    second <- unlist(second)
    
    data.frame(
      chain = names(aligned_chains)[ch_id],
      parameter = vars,
      first_half_mean = unname(first),
      second_half_mean = unname(second),
      absolute_difference = abs(first - second)
    )
  })
  
  do.call(rbind, results)
}

# Apply the check separately to each simulation scenario.
stability_separated <- half_summary_difference(aligned_separated)
stability_overlapping <- half_summary_difference(aligned_overlapping)

cat("\nHalf-run stability: clearly separated mixture\n")
print(stability_separated, row.names = FALSE)

cat("\nHalf-run stability: overlapping mixture\n")
print(stability_overlapping, row.names = FALSE)








# I assessed iteration adequacy using ESS, classical R-hat, trace plots,
# and autocorrelation diagnostics.
# For the separated mixture, I found R-hat values close to one, ESS above 10,000,
# and stable trace plots, supporting the adequacy of 12,000 iterations per chain
# with 4,000 burn-in iterations. For the overlapping mixture, I observed lower 
# ESS, pronounced excursions, and irregular posterior distributions, indicating
# less efficient sampling. I therefore concluded that R-hat alone was 
# insufficient to confirm the adequacy of the same settings for the overlapping 
# mixture.
# Part B is complete.

############
## PART C ##
############

# The mixture likelihood is invariant under
# (mu1, mu2, pi) <-> (mu2, mu1, 1 - pi).
# The sampler is not constrained to keep mu1 < mu2. I first inspect
# the raw draws, then apply ordering as a post-processing step.

# 1. Raw label-order diagnostics.
# A change in the sign of mu1 - mu2 indicates that the component labels
# exchanged their ordering between consecutive stored draws. The
# proportion mu1 < mu2 describes the share of draws in each ordering.
label_order_diagnostics <- function(chains, scenario_name) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    lower_first <- ch$mu1 < ch$mu2
    data.frame(
      scenario = scenario_name,
      chain = names(chains)[i],
      stored_draws = length(ch$mu1),
      proportion_mu1_less_mu2 = mean(lower_first),
      ordering_switches = sum(lower_first[-1] != lower_first[-length(lower_first)]),
      stringsAsFactors = FALSE
    )
  }))
}

label_order_separated <- label_order_diagnostics(
  chains_separated, "Clearly separated mixture"
)
label_order_overlapping <- label_order_diagnostics(
  chains_overlapping, "Substantially overlapping mixture"
)
label_order_results <- rbind(label_order_separated, label_order_overlapping)
cat("\nPART C: raw label-order diagnostics\n")
print(label_order_results, row.names = FALSE)

# 2. Plot the raw, unconstrained parameter traces for the overlapping case.
# These plots are produced before any relabelling so that the label issue
# is not concealed by post-processing.
par(mfrow = c(3, 1), mar = c(3.5, 4, 2.5, 1))
for (v in c("mu1", "mu2", "pi")) {
  all_values <- unlist(lapply(chains_overlapping, function(ch) ch[[v]]))
  plot(
    chains_overlapping[[1]][[v]], type = "l",
    ylim = range(all_values, finite = TRUE),
    xlab = "Stored post-burn-in draw", ylab = v,
    main = paste("Overlapping mixture: raw", v)
  )
  if (length(chains_overlapping) > 1) {
    for (j in 2:length(chains_overlapping)) {
      lines(chains_overlapping[[j]][[v]], col = j)
    }
    legend("topright", legend = names(chains_overlapping),
           col = seq_along(chains_overlapping), lty = 1, bty = "n")
  }
}
par(mfrow = c(1, 1))

# 3. Raw posterior summaries. These are shown before relabelling to make
# clear that mu1, mu2, and pi do not necessarily have stable component-
# specific interpretations when labels exchange.
raw_parameter_summary <- function(chains, scenario_name) {
  draws <- do.call(rbind, lapply(chains, function(ch) {
    data.frame(mu1 = ch$mu1, mu2 = ch$mu2, pi = ch$pi, sigma2 = ch$sigma2)
  }))
  summary <- t(vapply(names(draws), function(v) {
    posterior_summary(draws[[v]])
  }, numeric(5)))
  data.frame(
    scenario = scenario_name,
    parameter = rownames(summary),
    summary,
    row.names = NULL,
    check.names = FALSE
  )
}

raw_summary_separated_C <- raw_parameter_summary(
  chains_separated, "Clearly separated mixture"
)
raw_summary_overlapping_C <- raw_parameter_summary(
  chains_overlapping, "Substantially overlapping mixture"
)
raw_summary_C <- rbind(raw_summary_separated_C, raw_summary_overlapping_C)
cat("\nPART C: raw posterior summaries before relabelling\n")
print(raw_summary_C, row.names = FALSE)

# 4. Post-process each draw by ordering the component means. This is an
# identifiability convention for interpretation, not a restriction placed
# on the Gibbs sampler itself. When the means are swapped, pi and the
# allocation labels are swapped consistently as well.
relabel_chain_C <- function(ch) {
  mu1_is_lower <- ch$mu1 <= ch$mu2
  mu_low <- pmin(ch$mu1, ch$mu2)
  mu_high <- pmax(ch$mu1, ch$mu2)
  pi_low <- ifelse(mu1_is_lower, ch$pi, 1 - ch$pi)
  z_ordered <- ch$z
  switched_rows <- which(!mu1_is_lower)
  if (length(switched_rows) > 0) {
    z_ordered[switched_rows, ] <- 3L - z_ordered[switched_rows, , drop = FALSE]
  }
  list(mu_low = mu_low, mu_high = mu_high, pi_low = pi_low,
       sigma2 = ch$sigma2, z_ordered = z_ordered)
}

summarise_relabelled_C <- function(chains, scenario_name) {
  ordered_chains <- lapply(chains, relabel_chain_C)
  combined <- do.call(rbind, lapply(ordered_chains, function(ch) {
    data.frame(mu_low = ch$mu_low, mu_high = ch$mu_high,
               pi_low = ch$pi_low, sigma2 = ch$sigma2)
  }))
  summary <- t(vapply(names(combined), function(v) {
    posterior_summary(combined[[v]])
  }, numeric(5)))
  data.frame(
    scenario = scenario_name,
    parameter = rownames(summary),
    summary,
    row.names = NULL,
    check.names = FALSE
  )
}

relabelled_summary_separated_C <- summarise_relabelled_C(
  chains_separated, "Clearly separated mixture"
)
relabelled_summary_overlapping_C <- summarise_relabelled_C(
  chains_overlapping, "Substantially overlapping mixture"
)
relabelled_summary_C <- rbind(
  relabelled_summary_separated_C,
  relabelled_summary_overlapping_C
)
cat("\nPART C: posterior summaries after ordering mu_low < mu_high\n")
print(relabelled_summary_C, row.names = FALSE)

# 5. Label-invariant posterior quantities. The absolute component
# separation and the mixture mean are unchanged by exchanging labels.
label_invariant_summary_C <- function(chains, scenario_name) {
  draws <- do.call(rbind, lapply(chains, function(ch) {
    data.frame(
      separation = abs(ch$mu1 - ch$mu2),
      mixture_mean = ch$pi * ch$mu1 + (1 - ch$pi) * ch$mu2
    )
  }))
  summary <- t(vapply(names(draws), function(v) {
    posterior_summary(draws[[v]])
  }, numeric(5)))
  data.frame(
    scenario = scenario_name,
    quantity = rownames(summary),
    summary,
    row.names = NULL,
    check.names = FALSE
  )
}

invariant_separated_C <- label_invariant_summary_C(
  chains_separated, "Clearly separated mixture"
)
invariant_overlapping_C <- label_invariant_summary_C(
  chains_overlapping, "Substantially overlapping mixture"
)
label_invariant_C <- rbind(invariant_separated_C, invariant_overlapping_C)
cat("\nPART C: label-invariant posterior quantities\n")
print(label_invariant_C, row.names = FALSE)



############
## PART D ##
############

data(faithful)
y <- faithful$eruptions

#Visualise the observed distribution before fitting the mixture model.

hist(
  y,
  breaks = 20,
  probability = TRUE,
  main = "Old Faithful eruption durations",
  xlab = "Eruption duration (minutes)"
)

# Standardise the eruption durations using the helper from Part A.
# This can make prior specification and MCMC computation easier.
y_info_D <- standardise_data(y)
y_std_D <- y_info_D$y_std

# Set reproducible initial values for the Gibbs sampler.
set.seed(4101)

init_D <- list(
  mu1 = -0.8,
  mu2 = 0.8,
  pi = 0.5,
  sigma2 = 1,
  z = sample(1:2, length(y_std_D), replace = TRUE)
)
# Fit the two-component Gaussian mixture to the standardised data.
set.seed(4102)
draws_D <- gibbs_mix(
  y = y_std_D,
  n_iter = 12000,
  burn = 4000,
  thin = 2,
  init = init_D,
  m0 = m0,
  s0_sq = s0_sq,
  a = a,
  b = b,
  alpha0 = alpha0,
  beta0 = beta0,
  seed = 4102
)

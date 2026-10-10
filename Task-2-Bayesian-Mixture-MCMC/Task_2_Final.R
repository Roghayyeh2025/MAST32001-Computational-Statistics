
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

# Next, I will work on implementing the Gibbs sampler,
# storing posterior draws, assessing Markov chain behaviour,
# and justifying the choice of the number of iterations.

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

# Remaining consideration:
# The chosen number of iterations and burn-in are preliminary.
#I think we should justify them using convergence diagnostics,  
#such as trace plots, autocorrelation, ESS, and comparisons 
#between multiple chains, and check whether the posterior  
#summaries remain stable. We can work on this as 
#the next step to complete Part A.



############
## PART B ##
############




############
## PART C ##
############




############
## PART D ##
############

data(faithful)
y <- faithful$eruptions




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



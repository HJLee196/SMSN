rm(list = ls())

source("./SMSN_sample base code.R")

# Generate Data
mu1 = c(0, 0)
Sigma1 = matrix(c(1, 0,
                  0, 1), nrow = 2)


mu2 = c(1, 1)
Sigma2 = matrix(c(0.1, 0,
                  0, 0.1), nrow = 2)
Lambda2 = matrix(c(3, 0,
                   0, 3), nrow = 2)

pi_true_set = c(0.5, 0.5)

n_sample = 200

set.seed(2020)

X = NULL
label_true = NULL

X = 
  rbind(
    LaplacesDemon::rmvl(n = n_sample*pi_true_set[1], mu = mu1, Sigma = Sigma1),
    rmsn_sahu(n = n_sample*pi_true_set[2], xi = mu2, Sigma = Sigma2, Lambda = Lambda2))

label_true = c(rep(1, times = n_sample*pi_true_set[1]), rep(2, times = n_sample*pi_true_set[2]))

X_df = cbind(X, label_true) %>% as.data.frame()
colnames(X_df) = c("x1", "x2", "label_true")
X_df$label_true = 
  as.character(X_df$label_true)

X_df %>%
  ggplot(aes(x = x1, y = x2, color = label_true)) + 
  geom_point() +
  theme_bw()

# Generate Initial Values
Q_support_initial_list = 
  Q_support_initial_generator(x = X, label_est = label_true)

# Fit the Model
result_CNM = 
  CNM_SMSN_mixture(X = X, 
                   Q_support_initial_list = Q_support_initial_list)

# Plot a Result
x1_set = seq(from = min(X[,1]), to = max(X[,1]), length.out = 300)
x2_set = seq(from = min(X[,2]), to = max(X[,2]), length.out = 300)
x_set = data.frame(x1 = rep(x1_set, each = length(x2_set)),
                   x2 = rep(x2_set, times = length(x1_set)))

FSMSN_density = x_set
FSMSN_density$density = 
  density_mixture_Q_support_list(x = x_set, Q_support_list = result_CNM$Q_support_list)

FSMSN_density %>%
  ggplot(aes(x = x2, y = x1)) + 
  geom_contour(aes(z = density), alpha = 0.7, bins = 20) +
  geom_point(data = X_df, aes(x = x1, y = x2, color = label_true), alpha = 0.5) +
  theme_bw()

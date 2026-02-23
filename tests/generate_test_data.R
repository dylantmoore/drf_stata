#!/usr/bin/env Rscript
# Generate reference test data for drf_stata validation
# Pin: install drf from CRAN or github at a fixed version

library(drf)

set.seed(42)

# --- Test 1: Basic conditional mean ---
n <- 1000
p <- 4
X <- matrix(rnorm(n * p), nrow = n)
colnames(X) <- paste0("x", 1:p)

# Known signal: y = 3*x1 + 2*x2 + x3 + 0.3*x1^2 + noise
beta <- c(3, 2, 1, 0.5)
y <- X %*% beta + 0.3 * X[,1]^2 + rnorm(n, 0, 0.5)

# Split into train/test
n_train <- 800
n_test <- 200
X_train <- X[1:n_train, ]
Y_train <- matrix(y[1:n_train], ncol = 1)
X_test <- X[(n_train+1):n, ]
y_test <- y[(n_train+1):n]

# Fit DRF
forest <- drf(X_train, Y_train, num.trees = 500, seed = 42)

# Predict conditional mean
pred_mean <- predict(forest, X_test, functional = "mean")
ref_mean <- pred_mean$mean

# Predict conditional quantiles (0.1, 0.5, 0.9)
pred_q10 <- predict(forest, X_test, functional = "quantile", quantiles = 0.1)
pred_q50 <- predict(forest, X_test, functional = "quantile", quantiles = 0.5)
pred_q90 <- predict(forest, X_test, functional = "quantile", quantiles = 0.9)

ref_q10 <- pred_q10$quantile
ref_q50 <- pred_q50$quantile
ref_q90 <- pred_q90$quantile

# Predict weights (for a few test points)
pred_weights <- predict(forest, X_test[1:5, , drop=FALSE])
# pred_weights$weights is a sparse 5 x 800 matrix

# Build output CSV
df_test <- data.frame(
  X_test,
  y_true = y_test,
  ref_mean = as.numeric(ref_mean),
  ref_q10 = as.numeric(ref_q10),
  ref_q50 = as.numeric(ref_q50),
  ref_q90 = as.numeric(ref_q90)
)
write.csv(df_test, "test_condmean_test.csv", row.names = FALSE)

df_train <- data.frame(
  X_train,
  y = as.numeric(Y_train)
)
write.csv(df_train, "test_condmean_train.csv", row.names = FALSE)

# Save weights for first 5 test points (dense, for validation)
w_dense <- as.matrix(pred_weights$weights)
df_weights <- data.frame(w_dense[1:5, ])
colnames(df_weights) <- paste0("w_train_", 1:n_train)
write.csv(df_weights, "test_weights_5pts.csv", row.names = FALSE)

cat("Test 1 (conditional mean, univariate):\n")
cat("  R mean prediction correlation with truth:", cor(ref_mean, y_test), "\n")
cat("  R q50 prediction correlation with truth:", cor(ref_q50, y_test), "\n")

# --- Test 2: Multivariate response (conditional covariance) ---
set.seed(123)
n2 <- 500
p2 <- 3
X2 <- matrix(rnorm(n2 * p2), nrow = n2)
colnames(X2) <- paste0("x", 1:p2)

# Bivariate response with x-dependent correlation
y1 <- 2 * X2[,1] + rnorm(n2, 0, 0.5)
rho <- 0.3 * X2[,1]  # correlation depends on x1
y2 <- rho * y1 + sqrt(1 - rho^2) * rnorm(n2, 0, 0.5) + X2[,2]
Y2 <- cbind(y1, y2)

n_train2 <- 400
n_test2 <- 100
X2_train <- X2[1:n_train2, ]
Y2_train <- Y2[1:n_train2, ]
X2_test <- X2[(n_train2+1):n2, ]
Y2_test <- Y2[(n_train2+1):n2, ]

forest2 <- drf(X2_train, Y2_train, num.trees = 500, seed = 123)

# Conditional mean for multivariate
pred_mean2 <- predict(forest2, X2_test, functional = "mean")
ref_mean2_y1 <- pred_mean2$mean[, 1]
ref_mean2_y2 <- pred_mean2$mean[, 2]

# Conditional correlation
pred_cor <- predict(forest2, X2_test, functional = "cor")
ref_cor12 <- sapply(1:n_test2, function(i) pred_cor$cor[[i]][1, 2])

df_test2 <- data.frame(
  X2_test,
  y1_true = Y2_test[,1],
  y2_true = Y2_test[,2],
  ref_mean_y1 = as.numeric(ref_mean2_y1),
  ref_mean_y2 = as.numeric(ref_mean2_y2),
  ref_cor12 = ref_cor12
)
write.csv(df_test2, "test_multivar_test.csv", row.names = FALSE)

df_train2 <- data.frame(
  X2_train,
  y1 = Y2_train[,1],
  y2 = Y2_train[,2]
)
write.csv(df_train2, "test_multivar_train.csv", row.names = FALSE)

cat("\nTest 2 (multivariate, conditional correlation):\n")
cat("  R mean y1 corr with truth:", cor(ref_mean2_y1, Y2_test[,1]), "\n")
cat("  R mean y2 corr with truth:", cor(ref_mean2_y2, Y2_test[,2]), "\n")

# --- Test 3: Small stress test ---
set.seed(999)
n3 <- 200
p3 <- 20
X3 <- matrix(rnorm(n3 * p3), nrow = n3)
colnames(X3) <- paste0("x", 1:p3)
beta3 <- c(3, 2, 1, rep(0, p3 - 3))
y3 <- X3 %*% beta3 + rnorm(n3, 0, 1)

n_train3 <- 150
X3_train <- X3[1:n_train3, ]
Y3_train <- matrix(y3[1:n_train3], ncol = 1)
X3_test <- X3[(n_train3+1):n3, ]
y3_test <- y3[(n_train3+1):n3]

forest3 <- drf(X3_train, Y3_train, num.trees = 500, seed = 999)
pred_mean3 <- predict(forest3, X3_test, functional = "mean")

df_test3 <- data.frame(X3_test, y_true = y3_test, ref_mean = as.numeric(pred_mean3$mean))
write.csv(df_test3, "test_highp_test.csv", row.names = FALSE)

df_train3 <- data.frame(X3_train, y = as.numeric(Y3_train))
write.csv(df_train3, "test_highp_train.csv", row.names = FALSE)

cat("\nTest 3 (high-dimensional p=20):\n")
cat("  R mean corr with truth:", cor(pred_mean3$mean, y3_test), "\n")

cat("\nReference data generation complete.\n")
cat("R package version:", as.character(packageVersion("drf")), "\n")

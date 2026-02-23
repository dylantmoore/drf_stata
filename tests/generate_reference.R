# generate_reference.R -- Generate R reference data for validating Stata DRF plugin
# Runs the R drf package and saves inputs + outputs for comparison

library(drf)

cat("=== Generating DRF reference data ===\n")

# ---- Scenario 1: Linear signal (conditional mean) ----
set.seed(42)
n <- 500
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
y <- 3 * x1 + 2 * x2 + rnorm(n, 0, 0.5)

X <- cbind(x1, x2, x3)

# Fit DRF
fit1 <- drf(X = X, Y = matrix(y), num.trees = 200, seed = 42,
            min.node.size = 15, num.features = 10,
            sample.fraction = 0.5, honesty = TRUE,
            honesty.fraction = 0.5)

# Predict (in-sample with OOB weights)
w1 <- predict(fit1, newdata = X)$weights
pred_mean1 <- as.numeric(w1 %*% y)

# Save
df1 <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r_pred_mean = pred_mean1)
write.csv(df1, "tests/ref_linear.csv", row.names = FALSE)

cat("Scenario 1 (linear mean):\n")
cat("  Correlation(pred, y) =", cor(pred_mean1, y), "\n")
cat("  Mean pred =", mean(pred_mean1), "\n")
cat("  SD pred =", sd(pred_mean1), "\n")

# ---- Scenario 2: Nonlinear signal ----
set.seed(123)
n2 <- 1000
x1 <- rnorm(n2)
x2 <- rnorm(n2)
x3 <- rnorm(n2)
y2 <- x1^2 + sin(x2) + rnorm(n2, 0, 0.3)

X2 <- cbind(x1, x2, x3)

fit2 <- drf(X = X2, Y = matrix(y2), num.trees = 300, seed = 123,
            min.node.size = 15, num.features = 10,
            sample.fraction = 0.5, honesty = TRUE,
            honesty.fraction = 0.5)

w2 <- predict(fit2, newdata = X2)$weights
pred_mean2 <- as.numeric(w2 %*% y2)

df2 <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y2, r_pred_mean = pred_mean2)
write.csv(df2, "tests/ref_nonlinear.csv", row.names = FALSE)

cat("\nScenario 2 (nonlinear mean):\n")
cat("  Correlation(pred, y) =", cor(pred_mean2, y2), "\n")
cat("  Mean pred =", mean(pred_mean2), "\n")
cat("  SD pred =", sd(pred_mean2), "\n")

# ---- Scenario 3: Quantile (median) on linear signal ----
# Re-use scenario 1 data
w3 <- predict(fit1, newdata = X)$weights
# Compute weighted median for each observation
pred_q50 <- numeric(n)
for (i in 1:n) {
    wi <- w3[i, ]
    # Sort by y, compute weighted CDF, find median
    ord <- order(y)
    y_sorted <- y[ord]
    w_sorted <- wi[ord]
    cum_w <- cumsum(w_sorted)
    cum_w <- cum_w / cum_w[length(cum_w)]
    # Linear interpolation at 0.5
    idx <- which(cum_w >= 0.5)[1]
    if (idx == 1) {
        pred_q50[i] <- y_sorted[1]
    } else {
        # Linear interpolation
        f0 <- cum_w[idx - 1]
        f1 <- cum_w[idx]
        pred_q50[i] <- y_sorted[idx - 1] + (0.5 - f0) / (f1 - f0) * (y_sorted[idx] - y_sorted[idx - 1])
    }
}

df3 <- data.frame(x1 = df1$x1, x2 = df1$x2, x3 = df1$x3, y = df1$y,
                  r_pred_q50 = pred_q50)
write.csv(df3, "tests/ref_quantile.csv", row.names = FALSE)

cat("\nScenario 3 (median quantile):\n")
cat("  Correlation(pred_q50, y) =", cor(pred_q50, y), "\n")
cat("  Mean pred_q50 =", mean(pred_q50), "\n")
cat("  SD pred_q50 =", sd(pred_q50), "\n")

cat("\n=== Reference data written to tests/ref_*.csv ===\n")

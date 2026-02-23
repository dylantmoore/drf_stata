# generate_reference_large.R
# Larger, higher-dimensional test scenarios for DRF validation

library(drf)
cat("drf version:", as.character(packageVersion("drf")), "\n\n")

outdir <- "tests/ref"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

run_scenario <- function(name, X, y, ntrees, seed, functional = "mean",
                         quantile_prob = 0.5, ...) {
    cat("=== Scenario:", name, "===\n")
    cat("  n =", nrow(X), " p =", ncol(X), " trees =", ntrees, "\n")

    t0 <- proc.time()
    fit <- drf(X = X, Y = matrix(y), num.trees = ntrees, seed = seed,
               ci.group.size = 1, ...)
    t1 <- proc.time()

    w <- predict(fit, newdata = X)$weights
    t2 <- proc.time()

    if (functional == "mean") {
        pred <- as.numeric(w %*% y)
    } else if (functional == "quantile") {
        pred <- numeric(nrow(X))
        for (i in 1:nrow(X)) {
            wi <- w[i, ]
            ord <- order(y)
            y_sorted <- y[ord]
            w_sorted <- wi[ord]
            cum_w <- cumsum(w_sorted)
            cum_w <- cum_w / cum_w[length(cum_w)]
            idx <- which(cum_w >= quantile_prob)[1]
            if (idx == 1) {
                pred[i] <- y_sorted[1]
            } else {
                f0 <- cum_w[idx - 1]
                f1 <- cum_w[idx]
                pred[i] <- y_sorted[idx - 1] +
                    (quantile_prob - f0) / (f1 - f0) *
                    (y_sorted[idx] - y_sorted[idx - 1])
            }
        }
    }
    t3 <- proc.time()

    r_corr <- cor(pred, y)
    cat("  r(pred, y) =", round(r_corr, 4), "\n")
    cat("  mean(pred) =", round(mean(pred), 4),
        " sd(pred) =", round(sd(pred), 4), "\n")
    cat("  Train time:", round((t1 - t0)[3], 1), "s",
        " Predict time:", round((t2 - t1)[3], 1), "s\n\n")

    df <- as.data.frame(X)
    colnames(df) <- paste0("x", 1:ncol(X))
    df$y <- y
    df$r_pred <- pred

    write.csv(df, file.path(outdir, paste0(name, ".csv")), row.names = FALSE)
}

# ---- Scenario A: n=10000, p=20, sparse linear ----
# Only 5 of 20 predictors matter
set.seed(100)
n <- 10000
p <- 20
X <- matrix(rnorm(n * p), n, p)
beta <- rep(0, p)
beta[c(1, 3, 7, 12, 18)] <- c(3, -2, 1.5, -1, 0.8)
y <- X %*% beta + rnorm(n, 0, 1)
run_scenario("large_sparse_linear", X, as.numeric(y),
             ntrees = 300, seed = 100)

# ---- Scenario B: n=10000, p=30, nonlinear + interactions ----
set.seed(200)
n <- 10000
p <- 30
X <- matrix(rnorm(n * p), n, p)
y <- 2 * X[,1]^2 + sin(3 * X[,2]) + X[,3] * X[,4] -
     abs(X[,5]) + 0.5 * X[,6] * (X[,7] > 0) + rnorm(n, 0, 0.5)
run_scenario("large_nonlinear_interactions", X, as.numeric(y),
             ntrees = 300, seed = 200)

# ---- Scenario C: n=20000, p=15, heteroskedastic (quantile test) ----
# Variance depends on x1: good test for distributional estimation
set.seed(300)
n <- 20000
p <- 15
X <- matrix(rnorm(n * p), n, p)
y <- 2 * X[,1] + X[,2] + (1 + abs(X[,1])) * rnorm(n)
# q10 and q90 should diverge more for large |x1|
run_scenario("large_heteroskedastic_mean", X, as.numeric(y),
             ntrees = 300, seed = 300)
run_scenario("large_heteroskedastic_q10", X, as.numeric(y),
             ntrees = 300, seed = 300,
             functional = "quantile", quantile_prob = 0.1)
run_scenario("large_heteroskedastic_q90", X, as.numeric(y),
             ntrees = 300, seed = 300,
             functional = "quantile", quantile_prob = 0.9)

cat("=== All large reference data written ===\n")

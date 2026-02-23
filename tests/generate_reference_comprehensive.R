# generate_reference_comprehensive.R
# Run matching scenarios through R's drf and save inputs + predictions
# for comparison against Stata

library(drf)
cat("drf version:", as.character(packageVersion("drf")), "\n\n")

outdir <- "tests/ref"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

run_scenario <- function(name, X, y, ntrees, seed, functional = "mean",
                         quantile_prob = 0.5, honesty = TRUE,
                         min.node.size = 15, sample.fraction = 0.5,
                         num.features = 10, honesty.fraction = 0.5,
                         bandwidth = NULL) {
    cat("=== Scenario:", name, "===\n")

    fit <- drf(X = X, Y = matrix(y), num.trees = ntrees, seed = seed,
               min.node.size = min.node.size, num.features = num.features,
               sample.fraction = sample.fraction, honesty = honesty,
               honesty.fraction = honesty.fraction, bandwidth = bandwidth,
               ci.group.size = 1)

    # Get OOB weights
    w <- predict(fit, newdata = X)$weights

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

    r_corr <- cor(pred, y)
    cat("  n =", nrow(X), " p =", ncol(X), " trees =", ntrees, "\n")
    cat("  r(pred, y) =", round(r_corr, 4), "\n")
    cat("  mean(pred) =", round(mean(pred), 4), " sd(pred) =", round(sd(pred), 4), "\n\n")

    df <- as.data.frame(X)
    colnames(df) <- paste0("x", 1:ncol(X))
    df$y <- y
    df$r_pred <- pred

    write.csv(df, file.path(outdir, paste0(name, ".csv")), row.names = FALSE)
}

# ---- Scenario 1: Single predictor, conditional mean ----
set.seed(1)
n <- 300; x <- rnorm(n); y <- 2*x + rnorm(n, 0, 0.5)
run_scenario("single_pred_mean", matrix(x), y, ntrees = 100, seed = 1)

# ---- Scenario 2: Many predictors (p=10), conditional mean ----
set.seed(2)
n <- 500
X <- matrix(rnorm(n * 10), n, 10)
y <- 3*X[,1] + 2*X[,2] - X[,3] + rnorm(n, 0, 1)
run_scenario("many_pred_mean", X, y, ntrees = 200, seed = 2)

# ---- Scenario 3: nohonesty ----
set.seed(3)
n <- 300; X <- matrix(rnorm(n * 2), n, 2)
y <- 2*X[,1] + X[,2] + rnorm(n, 0, 0.5)
run_scenario("nohonesty", X, y, ntrees = 100, seed = 3, honesty = FALSE)

# ---- Scenario 4: Quantile 0.1 ----
set.seed(4)
n <- 500; X <- matrix(rnorm(n * 2), n, 2)
y <- 3*X[,1] + rnorm(n, 0, 1)
run_scenario("quantile_10", X, y, ntrees = 200, seed = 4,
             functional = "quantile", quantile_prob = 0.1)

# ---- Scenario 5: Quantile 0.9 ----
run_scenario("quantile_90", X, y, ntrees = 200, seed = 4,
             functional = "quantile", quantile_prob = 0.9)

# ---- Scenario 6: Quantile 0.5 (median) ----
run_scenario("quantile_50", X, y, ntrees = 200, seed = 4,
             functional = "quantile", quantile_prob = 0.5)

# ---- Scenario 7: Custom bandwidth ----
set.seed(7)
n <- 300; X <- matrix(rnorm(n * 2), n, 2)
y <- 2*X[,1] + rnorm(n, 0, 0.5)
run_scenario("custom_bw", X, y, ntrees = 100, seed = 7, bandwidth = 1.5)

# ---- Scenario 8: Custom minnodesize ----
run_scenario("custom_minnodesize", X, y, ntrees = 100, seed = 7,
             min.node.size = 5)

# ---- Scenario 9: Custom samplefrac ----
run_scenario("custom_samplefrac", X, y, ntrees = 100, seed = 7,
             sample.fraction = 0.7)

# ---- Scenario 10: Custom numfeatures ----
run_scenario("custom_numfeatures", X, y, ntrees = 100, seed = 7,
             num.features = 20)

# ---- Scenario 11: Custom honestyfrac ----
run_scenario("custom_honestyfrac", X, y, ntrees = 100, seed = 7,
             honesty.fraction = 0.3)

# ---- Scenario 12: Small dataset (n=30) ----
set.seed(14)
n <- 30; X <- matrix(rnorm(n * 2), n, 2)
y <- 3*X[,1] + rnorm(n, 0, 0.5)
run_scenario("small_n30", X, y, ntrees = 50, seed = 14, min.node.size = 3)

# ---- Scenario 13: Larger dataset (n=5000) ----
set.seed(15)
n <- 5000; X <- matrix(rnorm(n * 3), n, 3)
y <- 2*X[,1] - X[,2] + 0.5*X[,3] + rnorm(n, 0, 1)
run_scenario("large_n5000", X, y, ntrees = 200, seed = 15)

# ---- Scenario 14: Nonlinear signal ----
set.seed(123)
n <- 1000; X <- matrix(rnorm(n * 3), n, 3)
y <- X[,1]^2 + sin(X[,2]) + rnorm(n, 0, 0.3)
run_scenario("nonlinear", X, y, ntrees = 300, seed = 123)

# ---- Scenario 15: Auto-like real data pattern ----
set.seed(22)
n <- 200; X <- matrix(rnorm(n * 4), n, 4)
y <- 5000 + 1000*X[,1] - 500*X[,2] + 200*X[,3]*X[,4] + rnorm(n, 0, 500)
run_scenario("auto_like", X, y, ntrees = 200, seed = 22)

cat("=== All reference data written to", outdir, "===\n")

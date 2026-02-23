/*
 * drf_plugin_cpp.cpp -- Stata plugin wrapping the drf C++ library
 *
 * This wraps the actual drf/grf C++ code (Cevid, Michel, Meinshausen,
 * Buhlmann 2022) to produce OOB distributional random forest predictions
 * (conditional mean or conditional quantiles) directly in Stata.
 *
 * Compile (macOS arm64):
 *   See Makefile
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>

/* Stata plugin interface -- must be included with C linkage */
extern "C" {
#include "stplugin.h"
}

/* drf C++ library headers */
#include "commons/DefaultData.h"
#include "commons/globals.h"
#include "forest/Forest.h"
#include "forest/ForestOptions.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "prediction/collector/SampleWeightComputer.h"
#include "prediction/collector/TreeTraverser.h"

/* ================================================================
 * Median heuristic for bandwidth (matches R's medianHeuristic)
 * ================================================================
 * Computes sqrt(median(sqrt(dist(Y)/2))) where dist is pairwise
 * Euclidean distance. Uses subsample of 5000 if n > 5000.
 */
static double median_heuristic(const double* Y, int n, int d) {
    /* Subsample if too large */
    int use_n = n;
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    if (n > 5000) {
        /* Simple deterministic subsample: take evenly spaced */
        use_n = 5000;
        std::vector<int> sub(use_n);
        for (int i = 0; i < use_n; i++) {
            sub[i] = (int)((long long)i * n / use_n);
        }
        indices = sub;
    }

    /* Compute all pairwise sqrt(dist/2) values */
    size_t npairs = (size_t)use_n * (use_n - 1) / 2;
    std::vector<double> dists;
    dists.reserve(npairs);

    for (int i = 0; i < use_n; i++) {
        for (int j = i + 1; j < use_n; j++) {
            double dist2 = 0.0;
            for (int k = 0; k < d; k++) {
                double diff = Y[indices[i] * d + k] - Y[indices[j] * d + k];
                dist2 += diff * diff;
            }
            dists.push_back(std::sqrt(dist2 / 2.0));
        }
    }

    /* Compute median */
    if (dists.empty()) return 1.0;
    size_t mid = dists.size() / 2;
    std::nth_element(dists.begin(), dists.begin() + mid, dists.end());
    double median = dists[mid];

    return std::sqrt(median);
}

/* ================================================================
 * Weighted quantile (matches R's weighted.quantile)
 * ================================================================ */
static double weighted_quantile(const std::vector<std::pair<double, double>>& sorted_vw,
                                 double prob) {
    /* sorted_vw: pairs of (value, weight), sorted by value, weights sum to 1 */
    if (sorted_vw.empty()) return 0.0;
    if (sorted_vw.size() == 1) return sorted_vw[0].first;

    /* Cumulative weights */
    double cum = 0.0;
    for (size_t i = 0; i < sorted_vw.size(); i++) {
        cum += sorted_vw[i].second;
        if (cum >= prob) {
            return sorted_vw[i].first;
        }
    }
    return sorted_vw.back().first;
}

/* ================================================================
 * stata_call -- Main Entry Point
 * ================================================================
 *
 * argv[] layout (passed from drf.ado):
 *   [0]  n_trees          (int)
 *   [1]  seed             (int)
 *   [2]  mtry             (int, 0 = auto)
 *   [3]  min_node_size    (int)
 *   [4]  bandwidth        (double, 0 = auto median heuristic)
 *   [5]  sample_fraction  (double)
 *   [6]  alpha            (double)
 *   [7]  num_features     (int, random Fourier features)
 *   [8]  functional_type  (int: 0=mean, 1=quantile)
 *   [9]  quantile_prob    (double)
 *   [10] honesty          (int: 0=no, 1=yes)
 *   [11] honesty_fraction (double)
 *
 * Variable layout:
 *   Var 1:              Y (depvar)
 *   Var 2 .. nvar-1:    X (indepvars)
 *   Var nvar:           output (predictions written here)
 */
extern "C" STDLL stata_call(int argc, char *argv[])
{
    char msg[512];

    /* ----------------------------------------------------------
     * Step 0: Validate arguments
     * ---------------------------------------------------------- */
    if (argc < 12) {
        snprintf(msg, sizeof(msg),
                 "DRF error: expected 12 arguments, got %d\n", argc);
        SF_error(msg);
        return 198;
    }

    int n_trees        = atoi(argv[0]);
    int seed           = atoi(argv[1]);
    int mtry           = atoi(argv[2]);
    int min_node_size  = atoi(argv[3]);
    double bandwidth   = atof(argv[4]);
    double sample_frac = atof(argv[5]);
    double alpha       = atof(argv[6]);
    int num_features   = atoi(argv[7]);
    int func_type      = atoi(argv[8]);
    double quant_prob  = atof(argv[9]);
    int honesty        = atoi(argv[10]);
    double honesty_frac= atof(argv[11]);

    /* Defaults */
    if (n_trees <= 0)       n_trees = 500;
    if (min_node_size <= 0) min_node_size = 15;
    if (sample_frac <= 0.0 || sample_frac > 1.0) sample_frac = 0.5;
    if (alpha < 0.0 || alpha > 0.5) alpha = 0.05;
    if (num_features <= 0)  num_features = 10;
    if (func_type < 0 || func_type > 1) func_type = 0;
    if (quant_prob <= 0.0 || quant_prob >= 1.0) quant_prob = 0.5;
    if (honesty_frac <= 0.0 || honesty_frac >= 1.0) honesty_frac = 0.5;
    if (seed <= 0) seed = 42;

    /* ----------------------------------------------------------
     * Step 1: Read dimensions
     * ---------------------------------------------------------- */
    int nvar = SF_nvars();
    int p = nvar - 2;  /* number of X variables (Y + Xs + output) */
    if (p < 1) {
        SF_error("DRF error: need at least 1 predictor variable.\n");
        return 198;
    }

    /* Count usable observations */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    int n = 0;
    for (ST_int i = obs1; i <= obs2; i++) {
        if (SF_ifobs(i)) n++;
    }

    if (n < 2) {
        SF_error("DRF error: need at least 2 non-missing observations.\n");
        return 2000;
    }

    /* Auto mtry */
    if (mtry <= 0) {
        mtry = (int)std::ceil(std::sqrt((double)p) + 20.0);
        if (mtry > p) mtry = p;
    }
    if (mtry > p) mtry = p;

    snprintf(msg, sizeof(msg),
             "DRF (C++ grf engine): n=%d, p=%d, trees=%d, mtry=%d, "
             "min_node=%d, features=%d, honesty=%d\n",
             n, p, n_trees, mtry, min_node_size, num_features, honesty);
    SF_display(msg);

    /* ----------------------------------------------------------
     * Step 2: Read data from Stata into arrays
     * ---------------------------------------------------------- */
    /* Y values (original, for applying the functional) */
    std::vector<double> Y_orig(n);
    /* Y values (scaled, for the forest training) */
    std::vector<double> Y_scaled(n);
    /* X values in row-major */
    std::vector<double> X_row(n * p);
    /* Observation mapping */
    std::vector<int> obs_map(n);

    int idx = 0;
    for (ST_int i = obs1; i <= obs2; i++) {
        if (!SF_ifobs(i)) continue;

        double val;
        ST_retcode rc = SF_vdata(1, i, &val);
        if (rc || SF_is_missing(val)) continue;

        Y_orig[idx] = val;
        obs_map[idx] = (int)i;

        for (int j = 0; j < p; j++) {
            rc = SF_vdata(j + 2, i, &val);
            if (rc || SF_is_missing(val)) {
                val = 0.0;
            }
            X_row[idx * p + j] = val;
        }
        idx++;
    }
    n = idx;
    if (n < 2) {
        SF_error("DRF error: fewer than 2 complete observations.\n");
        return 2000;
    }

    /* Scale Y (z-score standardization, matching R's response.scaling=TRUE) */
    double y_mean = 0.0, y_sd = 0.0;
    for (int i = 0; i < n; i++) y_mean += Y_orig[i];
    y_mean /= n;
    for (int i = 0; i < n; i++) {
        double d = Y_orig[i] - y_mean;
        y_sd += d * d;
    }
    y_sd = std::sqrt(y_sd / (n - 1));
    if (y_sd < 1e-15) y_sd = 1.0;

    for (int i = 0; i < n; i++) {
        Y_scaled[i] = (Y_orig[i] - y_mean) / y_sd;
    }

    /* Auto bandwidth using median heuristic on scaled Y */
    if (bandwidth <= 0.0) {
        /* For univariate Y, compute median heuristic */
        bandwidth = median_heuristic(Y_scaled.data(), n, 1);
        if (bandwidth <= 0.0) bandwidth = 1.0;
        snprintf(msg, sizeof(msg), "  Auto bandwidth (median heuristic): %.6f\n", bandwidth);
        SF_display(msg);
    }

    /* ----------------------------------------------------------
     * Step 3: Create drf::DefaultData
     * ---------------------------------------------------------- */
    /* drf's DefaultData stores column-major: data[col * num_rows + row]
     * Columns: Y_scaled, X1, X2, ..., Xp
     * Total columns: 1 + p
     */
    size_t num_cols = (size_t)(1 + p);
    size_t num_rows = (size_t)n;

    std::vector<double> data_vec(num_cols * num_rows);
    /* Column 0: Y_scaled */
    for (size_t i = 0; i < num_rows; i++) {
        data_vec[0 * num_rows + i] = Y_scaled[i];
    }
    /* Columns 1..p: X */
    for (int j = 0; j < p; j++) {
        for (size_t i = 0; i < num_rows; i++) {
            data_vec[(j + 1) * num_rows + i] = X_row[i * p + j];
        }
    }

    std::unique_ptr<drf::Data> data(new drf::DefaultData(data_vec, num_rows, num_cols));

    /* Set outcome index (column 0 = Y) */
    std::vector<size_t> outcome_index = {0};
    data->set_outcome_index(outcome_index);

    /* Sort data (required for splitting rule efficiency) */
    data->sort();

    /* ----------------------------------------------------------
     * Step 4: Configure and train forest
     * ---------------------------------------------------------- */
    drf::ForestTrainer trainer = drf::fourier_trainer(outcome_index.size());

    /* ForestOptions parameters */
    size_t ci_group_size = 1;          /* No CI groups needed for OOB predictions */
    bool honesty_bool = (honesty != 0);
    bool honesty_prune_leaves = true;  /* Match R default */
    double imbalance_penalty = 0.0;    /* Match R default */
    drf::uint num_threads = 0;         /* 0 = auto (hardware concurrency) */
    std::vector<size_t> clusters;
    drf::uint samples_per_cluster = 0;
    unsigned int node_scaling = 0;     /* Match R default */

    drf::ForestOptions options(
        (drf::uint)n_trees,
        ci_group_size,
        sample_frac,
        (drf::uint)mtry,
        (drf::uint)min_node_size,
        honesty_bool,
        honesty_frac,
        honesty_prune_leaves,
        alpha,
        imbalance_penalty,
        num_threads,
        (drf::uint)seed,
        clusters,
        samples_per_cluster,
        (size_t)num_features,
        bandwidth,
        node_scaling
    );

    SF_display("  Training forest...\n");
    drf::Forest forest = trainer.train(*data, options);
    SF_display("  Forest trained.\n");

    /* ----------------------------------------------------------
     * Step 5: Compute OOB weights and predictions
     * ---------------------------------------------------------- */
    SF_display("  Computing OOB predictions...\n");

    /* Get leaf nodes and valid trees for OOB */
    drf::uint actual_threads = drf::ForestOptions::validate_num_threads(num_threads);
    drf::TreeTraverser tree_traverser(actual_threads);
    drf::SampleWeightComputer weight_computer;

    bool oob = true;
    std::vector<std::vector<size_t>> leaf_nodes_by_tree =
        tree_traverser.get_leaf_nodes(forest, *data, oob);
    std::vector<std::vector<bool>> trees_by_sample =
        tree_traverser.get_valid_trees_by_sample(forest, *data, oob);

    /* For each sample, compute OOB weights and apply the functional */
    int n_pred = 0;
    for (size_t sample = 0; sample < num_rows; sample++) {
        /* Compute weights for this sample */
        std::unordered_map<size_t, double> weights =
            weight_computer.compute_weights(
                sample, forest, leaf_nodes_by_tree, trees_by_sample,
                0, forest.get_trees().size());

        if (weights.empty()) {
            /* No OOB trees for this sample -- leave as missing */
            continue;
        }

        double prediction;

        if (func_type == 0) {
            /* Conditional mean: weighted average of Y_orig */
            prediction = 0.0;
            for (auto& entry : weights) {
                size_t neighbor = entry.first;
                double weight = entry.second;
                prediction += weight * Y_orig[neighbor];
            }
        } else {
            /* Conditional quantile */
            /* Collect (value, weight) pairs and sort by value */
            std::vector<std::pair<double, double>> vw;
            vw.reserve(weights.size());
            for (auto& entry : weights) {
                size_t neighbor = entry.first;
                double weight = entry.second;
                vw.push_back({Y_orig[neighbor], weight});
            }
            std::sort(vw.begin(), vw.end());
            prediction = weighted_quantile(vw, quant_prob);
        }

        /* Write prediction back to Stata */
        ST_retcode rc = SF_vstore(nvar, obs_map[sample], prediction);
        if (rc == 0) n_pred++;
    }

    snprintf(msg, sizeof(msg), "  Wrote %d OOB predictions.\n", n_pred);
    SF_display(msg);

    return 0;
}

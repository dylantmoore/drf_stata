/*
 * drf_plugin.c -- Distributional Random Forests for Stata
 *
 * Implements the DRF algorithm of Cevid, Michel, Meinshausen, and
 * Buhlmann (2022) as a Stata C plugin.  A single stata_call() entry
 * point reads data, trains a forest with FourierMMD splits, computes
 * OOB predictions (weighted mean or weighted quantile), and writes
 * them back into Stata.
 *
 * Compile (macOS arm64):
 *   gcc -O3 -fPIC -DSYSTEM=APPLEMAC -arch arm64 -bundle \
 *       -o drf_plugin.darwin-arm64.plugin drf_plugin.c stplugin.c -lm
 *
 * Compile (Linux x86_64):
 *   gcc -O3 -fPIC -DSYSTEM=OPUNIX -shared \
 *       -o drf_plugin.linux-x86_64.plugin drf_plugin.c stplugin.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <float.h>
#include "stplugin.h"

/* ================================================================
 * Section 1: XorShift128+ RNG
 * ================================================================
 * Fast, high-quality PRNG used throughout.  Each tree gets its own
 * state seeded deterministically from (base_seed + tree_index).
 */

typedef struct {
    uint64_t s[2];
} xorshift128p_state;

static uint64_t splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void xorshift128p_seed(xorshift128p_state *state, uint64_t seed)
{
    state->s[0] = splitmix64(&seed);
    state->s[1] = splitmix64(&seed);
    if (state->s[0] == 0 && state->s[1] == 0) {
        state->s[0] = 1;
    }
}

static uint64_t xorshift128p_next(xorshift128p_state *state)
{
    uint64_t s1 = state->s[0];
    uint64_t s0 = state->s[1];
    uint64_t result = s0 + s1;
    state->s[0] = s0;
    s1 ^= s1 << 23;
    state->s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return result;
}

static double xorshift128p_double(xorshift128p_state *state)
{
    /* Uniform in [0, 1) with 53-bit precision */
    return (double)(xorshift128p_next(state) >> 11) * 0x1.0p-53;
}

/* ================================================================
 * Section 2: Box-Muller Normal Generator
 * ================================================================ */

static double box_muller_normal(xorshift128p_state *rng, int *has_spare,
                                double *spare)
{
    if (*has_spare) {
        *has_spare = 0;
        return *spare;
    }

    double u, v, s;
    do {
        u = 2.0 * xorshift128p_double(rng) - 1.0;
        v = 2.0 * xorshift128p_double(rng) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    *spare = v * s;
    *has_spare = 1;
    return u * s;
}

/* ================================================================
 * Section 3: DRF Tree Data Structure
 * ================================================================ */

typedef struct {
    int *left_child;    /* 0 = leaf                              */
    int *right_child;   /* 0 = leaf                              */
    int *split_var;     /* variable index (0-based into X)       */
    double *split_value;

    int *leaf_samples;  /* flat array of sample indices in leaves */
    int *leaf_offset;   /* leaf_offset[node] = start index       */
    int *leaf_count;    /* leaf_count[node] = count              */

    char *in_bag;       /* in_bag[i] == 1 if obs i was used      */
    int n_nodes;
    int capacity;
    int n_total;        /* total number of observations           */
} DRFTree;

static int drf_tree_init(DRFTree *tree, int capacity, int n_total)
{
    tree->capacity = capacity;
    tree->n_nodes = 0;
    tree->n_total = n_total;

    tree->left_child  = (int *)calloc((size_t)capacity, sizeof(int));
    tree->right_child = (int *)calloc((size_t)capacity, sizeof(int));
    tree->split_var   = (int *)calloc((size_t)capacity, sizeof(int));
    tree->split_value = (double *)calloc((size_t)capacity, sizeof(double));
    tree->leaf_offset = (int *)calloc((size_t)capacity, sizeof(int));
    tree->leaf_count  = (int *)calloc((size_t)capacity, sizeof(int));
    tree->in_bag      = (char *)calloc((size_t)n_total, sizeof(char));

    /* leaf_samples allocated later after tree is built */
    tree->leaf_samples = NULL;

    if (!tree->left_child || !tree->right_child || !tree->split_var ||
        !tree->split_value || !tree->leaf_offset || !tree->leaf_count ||
        !tree->in_bag) {
        return -1;
    }
    return 0;
}

static void drf_tree_free(DRFTree *tree)
{
    free(tree->left_child);
    free(tree->right_child);
    free(tree->split_var);
    free(tree->split_value);
    free(tree->leaf_samples);
    free(tree->leaf_offset);
    free(tree->leaf_count);
    free(tree->in_bag);
    memset(tree, 0, sizeof(DRFTree));
}

static int drf_tree_add_node(DRFTree *tree)
{
    if (tree->n_nodes >= tree->capacity) return -1;
    int id = tree->n_nodes++;
    tree->left_child[id]  = 0;
    tree->right_child[id] = 0;
    tree->split_var[id]   = -1;
    tree->split_value[id] = 0.0;
    tree->leaf_offset[id] = 0;
    tree->leaf_count[id]  = 0;
    return id;
}

/* ================================================================
 * Section 4: Utility -- Sort helpers
 * ================================================================ */

/* Argsort: sort indices by corresponding values */
typedef struct {
    double value;
    int index;
} ValueIndex;

static int cmp_value_index(const void *a, const void *b)
{
    double va = ((const ValueIndex *)a)->value;
    double vb = ((const ValueIndex *)b)->value;
    if (va < vb) return -1;
    if (va > vb) return  1;
    return 0;
}

/* Sort two parallel arrays (y, w) by y in place */
static void sort_by_y(double *y, double *w, int n)
{
    /* Simple insertion sort for small n; otherwise use stdlib qsort
       via ValueIndex approach.  We sort indices, then permute. */
    if (n <= 1) return;

    ValueIndex *vi = (ValueIndex *)malloc((size_t)n * sizeof(ValueIndex));
    if (!vi) return;  /* degrade gracefully */

    for (int i = 0; i < n; i++) {
        vi[i].value = y[i];
        vi[i].index = i;
    }
    qsort(vi, (size_t)n, sizeof(ValueIndex), cmp_value_index);

    double *tmp_y = (double *)malloc((size_t)n * sizeof(double));
    double *tmp_w = (double *)malloc((size_t)n * sizeof(double));
    if (tmp_y && tmp_w) {
        for (int i = 0; i < n; i++) {
            tmp_y[i] = y[vi[i].index];
            tmp_w[i] = w[vi[i].index];
        }
        memcpy(y, tmp_y, (size_t)n * sizeof(double));
        memcpy(w, tmp_w, (size_t)n * sizeof(double));
    }
    free(tmp_y);
    free(tmp_w);
    free(vi);
}

/* ================================================================
 * Section 5: Quickselect for Median
 * ================================================================ */

static double quickselect_median(double *arr, int n, xorshift128p_state *rng)
{
    if (n <= 0) return 0.0;
    if (n == 1) return arr[0];

    /* Work on a copy so we don't destroy the input */
    double *buf = (double *)malloc((size_t)n * sizeof(double));
    if (!buf) {
        /* Fallback: sort in place */
        qsort(arr, (size_t)n, sizeof(double), cmp_value_index);
        return arr[n / 2];
    }
    memcpy(buf, arr, (size_t)n * sizeof(double));

    int lo = 0, hi = n - 1;
    int target = n / 2;
    while (lo < hi) {
        /* Random pivot */
        int pi = lo + (int)(xorshift128p_double(rng) * (hi - lo + 1));
        if (pi > hi) pi = hi;
        double pivot = buf[pi];

        /* Move pivot to end */
        double tmp = buf[pi]; buf[pi] = buf[hi]; buf[hi] = tmp;

        int store = lo;
        for (int i = lo; i < hi; i++) {
            if (buf[i] < pivot) {
                tmp = buf[i]; buf[i] = buf[store]; buf[store] = tmp;
                store++;
            }
        }
        tmp = buf[store]; buf[store] = buf[hi]; buf[hi] = tmp;

        if (store == target) break;
        else if (store < target) lo = store + 1;
        else hi = store - 1;
    }

    double result = buf[target];
    free(buf);
    return result;
}

/* ================================================================
 * Section 6: Weighted Quantile
 * ================================================================ */

static double weighted_quantile(double *y, double *w, int n, double prob)
{
    if (n <= 0) return 0.0;
    if (n == 1) return y[0];

    /* Make copies to sort */
    double *yc = (double *)malloc((size_t)n * sizeof(double));
    double *wc = (double *)malloc((size_t)n * sizeof(double));
    if (!yc || !wc) {
        free(yc); free(wc);
        return y[0]; /* fallback */
    }
    memcpy(yc, y, (size_t)n * sizeof(double));
    memcpy(wc, w, (size_t)n * sizeof(double));

    sort_by_y(yc, wc, n);

    /* Compute cumulative weight CDF */
    double total_w = 0.0;
    for (int i = 0; i < n; i++) total_w += wc[i];
    if (total_w <= 0.0) { free(yc); free(wc); return yc[0]; }

    double cum = 0.0;
    int left = 0;
    for (int i = 0; i < n; i++) {
        cum += wc[i];
        double F_i = cum / total_w;
        if (F_i <= prob) {
            left = i;
        } else {
            break;
        }
    }

    /* Recompute F[left] and F[left+1] for interpolation */
    double cum_left = 0.0;
    for (int i = 0; i <= left; i++) cum_left += wc[i];
    double F_left = cum_left / total_w;

    double result;
    if (F_left < prob && left + 1 < n) {
        double cum_right = cum_left + wc[left + 1];
        double F_right = cum_right / total_w;
        double frac = (prob - F_left) / (F_right - F_left);
        result = yc[left] + (yc[left + 1] - yc[left]) * frac;
    } else {
        result = yc[left];
    }

    free(yc);
    free(wc);
    return result;
}

/* ================================================================
 * Section 7: Median Heuristic Bandwidth
 * ================================================================
 * For univariate Y (d=1):
 *   1. Optionally subsample to 5000 observations.
 *   2. Compute all pairwise |Y_i - Y_j|.
 *   3. Apply sqrt(dist/2) to each pair distance.
 *   4. Take median via quickselect.
 *   5. bandwidth = sqrt(median).
 */

static double median_heuristic_bandwidth(const double *Y_scaled, int n,
                                         xorshift128p_state *rng)
{
    int nsub = n;
    int *sub_idx = NULL;

    if (n > 5000) {
        nsub = 5000;
        sub_idx = (int *)malloc((size_t)n * sizeof(int));
        if (!sub_idx) { nsub = n; }
        else {
            for (int i = 0; i < n; i++) sub_idx[i] = i;
            /* Fisher-Yates to select 5000 */
            for (int i = n - 1; i > 0 && i >= n - nsub; i--) {
                int j = (int)(xorshift128p_double(rng) * (i + 1));
                if (j > i) j = i;
                int tmp = sub_idx[i]; sub_idx[i] = sub_idx[j]; sub_idx[j] = tmp;
            }
        }
    }

    /* Number of pairs */
    long long npairs = (long long)nsub * (nsub - 1) / 2;
    if (npairs <= 0) {
        free(sub_idx);
        return 1.0;
    }

    /* Cap pairwise computation to avoid excessive memory use */
    int cap_nsub = nsub;
    if (npairs > 50000000LL) {
        /* Reduce subsample to keep pairs manageable */
        cap_nsub = 10000;
        if (cap_nsub > nsub) cap_nsub = nsub;
        npairs = (long long)cap_nsub * (cap_nsub - 1) / 2;
    }

    double *dists = (double *)malloc((size_t)npairs * sizeof(double));
    if (!dists) {
        free(sub_idx);
        return 1.0; /* fallback */
    }

    long long k = 0;
    for (int i = 0; i < cap_nsub; i++) {
        int ii = sub_idx ? sub_idx[n - 1 - i] : i;
        double yi = Y_scaled[ii];
        for (int j = i + 1; j < cap_nsub; j++) {
            int jj = sub_idx ? sub_idx[n - 1 - j] : j;
            double d = fabs(yi - Y_scaled[jj]);
            dists[k++] = sqrt(d / 2.0);
        }
    }

    double med = quickselect_median(dists, (int)k, rng);
    double bandwidth = sqrt(med);
    if (bandwidth < 1e-10) bandwidth = 1e-10;

    free(dists);
    free(sub_idx);
    return bandwidth;
}

/* ================================================================
 * Section 8: Fisher-Yates Shuffle (for subsampling)
 * ================================================================ */

static void fisher_yates_shuffle(int *arr, int n, xorshift128p_state *rng)
{
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(xorshift128p_double(rng) * (i + 1));
        if (j > i) j = i;
        int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }
}

/* ================================================================
 * Section 9: FourierMMD Splitting Rule
 * ================================================================
 *
 * For a node with k samples, this finds the best (variable, threshold)
 * split by maximizing the MMD criterion computed via random Fourier
 * features.
 *
 * Steps:
 *   1. Generate num_features omega scalars ~ N(0, 1/bandwidth^2).
 *   2. For each sample, compute cos/sin Fourier features.
 *   3. For each candidate variable (mtry sampled from p predictors):
 *      a. Sort node samples by X[variable].
 *      b. Sweep left-to-right accumulating cos/sin sums.
 *      c. Compute weighted MMD at each split point.
 *      d. Track the best split point.
 *   4. Return the variable and threshold with the highest MMD.
 */

static int fourier_mmd_split(
    const int *node_samples, int node_n,
    const double *Y_scaled,
    const double *X,         /* X[i * p + j] = predictor j of obs i */
    int p,                   /* number of predictors */
    int mtry,
    int num_features,
    double bandwidth,
    int min_node_size,
    double alpha,
    xorshift128p_state *rng,
    int *best_var_out,
    double *best_threshold_out)
{
    *best_var_out = -1;
    *best_threshold_out = 0.0;

    if (node_n < 2) return 0;

    /* Alpha-based minimum child size (NOT min_node_size — that's the
       BFS stopping criterion, handled in build_tree). */
    int min_child = (int)ceil(alpha * node_n);
    if (min_child < 1) min_child = 1;

    int has_spare = 0;
    double spare = 0.0;

    /* Step 1: Generate omega scalars */
    double *omega = (double *)malloc((size_t)num_features * sizeof(double));
    if (!omega) return 0;

    double sigma = 1.0 / (bandwidth * bandwidth);
    for (int l = 0; l < num_features; l++) {
        omega[l] = box_muller_normal(rng, &has_spare, &spare) * sqrt(sigma);
    }

    /* Step 2: Compute Fourier features for node samples */
    double *cos_feat = (double *)malloc((size_t)num_features * node_n * sizeof(double));
    double *sin_feat = (double *)malloc((size_t)num_features * node_n * sizeof(double));
    if (!cos_feat || !sin_feat) {
        free(omega); free(cos_feat); free(sin_feat);
        return 0;
    }

    for (int l = 0; l < num_features; l++) {
        for (int j = 0; j < node_n; j++) {
            double theta = omega[l] * Y_scaled[node_samples[j]];
            cos_feat[l * node_n + j] = cos(theta);
            sin_feat[l * node_n + j] = sin(theta);
        }
    }

    /* Pre-compute total cos/sin sums for each feature */
    double *cos_total = (double *)calloc((size_t)num_features, sizeof(double));
    double *sin_total = (double *)calloc((size_t)num_features, sizeof(double));
    if (!cos_total || !sin_total) {
        free(omega); free(cos_feat); free(sin_feat);
        free(cos_total); free(sin_total);
        return 0;
    }

    for (int l = 0; l < num_features; l++) {
        for (int j = 0; j < node_n; j++) {
            cos_total[l] += cos_feat[l * node_n + j];
            sin_total[l] += sin_feat[l * node_n + j];
        }
    }

    /* Step 3: Sample mtry candidate variables */
    int *cand_vars = (int *)malloc((size_t)p * sizeof(int));
    if (!cand_vars) {
        free(omega); free(cos_feat); free(sin_feat);
        free(cos_total); free(sin_total);
        return 0;
    }
    for (int j = 0; j < p; j++) cand_vars[j] = j;
    fisher_yates_shuffle(cand_vars, p, rng);
    int n_cand = mtry < p ? mtry : p;

    /* Working buffer for sorted indices */
    ValueIndex *sorted = (ValueIndex *)malloc((size_t)node_n * sizeof(ValueIndex));
    if (!sorted) {
        free(omega); free(cos_feat); free(sin_feat);
        free(cos_total); free(sin_total); free(cand_vars);
        return 0;
    }

    double best_mmd = -1.0;

    for (int c = 0; c < n_cand; c++) {
        int var = cand_vars[c];

        /* Sort node samples by X[var] */
        for (int j = 0; j < node_n; j++) {
            sorted[j].value = X[node_samples[j] * p + var];
            sorted[j].index = j; /* index into node_samples */
        }
        qsort(sorted, (size_t)node_n, sizeof(ValueIndex), cmp_value_index);

        /* Sweep left-to-right accumulating Fourier sums */
        double *cos_left = (double *)calloc((size_t)num_features, sizeof(double));
        double *sin_left = (double *)calloc((size_t)num_features, sizeof(double));
        if (!cos_left || !sin_left) {
            free(cos_left); free(sin_left);
            continue;
        }

        for (int j = 0; j < node_n - 1; j++) {
            int idx = sorted[j].index;

            for (int l = 0; l < num_features; l++) {
                cos_left[l] += cos_feat[l * node_n + idx];
                sin_left[l] += sin_feat[l * node_n + idx];
            }

            int n_left = j + 1;
            int n_right = node_n - n_left;

            /* Enforce minimum child size */
            if (n_left < min_child || n_right < min_child) continue;

            /* Skip tied X values */
            if (sorted[j].value == sorted[j + 1].value) continue;

            /* Compute MMD across Fourier features */
            double mmd = 0.0;
            for (int l = 0; l < num_features; l++) {
                double cr = cos_total[l] - cos_left[l]; /* cos_right */
                double sr = sin_total[l] - sin_left[l]; /* sin_right */

                double diff_real = cos_left[l] / n_left - cr / n_right;
                double diff_imag = sin_left[l] / n_left - sr / n_right;

                mmd += ((double)n_left * n_right / (double)node_n) *
                       (diff_real * diff_real + diff_imag * diff_imag) /
                       num_features;
            }

            if (mmd > best_mmd) {
                best_mmd = mmd;
                *best_var_out = var;
                *best_threshold_out = 0.5 * (sorted[j].value +
                                             sorted[j + 1].value);
            }
        }

        free(cos_left);
        free(sin_left);
    }

    free(omega);
    free(cos_feat);
    free(sin_feat);
    free(cos_total);
    free(sin_total);
    free(cand_vars);
    free(sorted);

    return (*best_var_out >= 0) ? 1 : 0;
}

/* ================================================================
 * Section 10: Tree Building (BFS)
 * ================================================================
 *
 * Uses a BFS queue to split nodes.  Each node holds a contiguous
 * slice of the working sample-index buffer.  On a split, samples
 * are partitioned in-place into left and right halves.
 */

typedef struct {
    int node_id;
    int start;   /* index into sample buffer */
    int count;   /* number of samples in this node */
} BFSItem;

static int build_tree(
    DRFTree *tree,
    int *samples,        /* working buffer; indices into global data */
    int n_samples,
    const double *Y_scaled,
    const double *X,
    int p,
    int mtry,
    int num_features,
    double bandwidth,
    int min_node_size,
    double alpha,
    xorshift128p_state *rng)
{
    /* BFS queue */
    int queue_cap = tree->capacity;
    BFSItem *queue = (BFSItem *)malloc((size_t)queue_cap * sizeof(BFSItem));
    if (!queue) return -1;

    int root = drf_tree_add_node(tree);
    if (root < 0) { free(queue); return -1; }

    int q_head = 0, q_tail = 0;
    queue[q_tail].node_id = root;
    queue[q_tail].start   = 0;
    queue[q_tail].count   = n_samples;
    q_tail++;

    /* Temporary per-node storage for leaf sample indices.
       We accumulate them into a dynamic buffer and finalize later. */
    int leaf_buf_cap = n_samples * 2;
    int *leaf_buf = (int *)malloc((size_t)leaf_buf_cap * sizeof(int));
    int leaf_buf_len = 0;
    if (!leaf_buf) { free(queue); return -1; }

    while (q_head < q_tail) {
        BFSItem item = queue[q_head++];
        int node = item.node_id;
        int start = item.start;
        int count = item.count;

        /* Leaf if too few samples */
        if (count <= min_node_size) {
            tree->leaf_offset[node] = leaf_buf_len;
            tree->leaf_count[node] = count;
            /* Ensure space */
            if (leaf_buf_len + count > leaf_buf_cap) {
                leaf_buf_cap = (leaf_buf_len + count) * 2;
                int *tmp = (int *)realloc(leaf_buf,
                                          (size_t)leaf_buf_cap * sizeof(int));
                if (!tmp) { free(queue); free(leaf_buf); return -1; }
                leaf_buf = tmp;
            }
            for (int i = 0; i < count; i++)
                leaf_buf[leaf_buf_len++] = samples[start + i];
            continue;
        }

        /* Attempt split */
        int split_var = -1;
        double split_val = 0.0;
        int found = fourier_mmd_split(
            samples + start, count,
            Y_scaled, X, p, mtry, num_features, bandwidth,
            min_node_size, alpha, rng,
            &split_var, &split_val);

        if (!found) {
            /* Make this a leaf */
            tree->leaf_offset[node] = leaf_buf_len;
            tree->leaf_count[node] = count;
            if (leaf_buf_len + count > leaf_buf_cap) {
                leaf_buf_cap = (leaf_buf_len + count) * 2;
                int *tmp = (int *)realloc(leaf_buf,
                                          (size_t)leaf_buf_cap * sizeof(int));
                if (!tmp) { free(queue); free(leaf_buf); return -1; }
                leaf_buf = tmp;
            }
            for (int i = 0; i < count; i++)
                leaf_buf[leaf_buf_len++] = samples[start + i];
            continue;
        }

        /* Partition samples: left goes <= split_val, right goes > */
        tree->split_var[node] = split_var;
        tree->split_value[node] = split_val;

        int lo = start;
        int hi = start + count - 1;
        while (lo <= hi) {
            if (X[samples[lo] * p + split_var] <= split_val) {
                lo++;
            } else {
                int tmp = samples[lo];
                samples[lo] = samples[hi];
                samples[hi] = tmp;
                hi--;
            }
        }
        int n_left = lo - start;
        int n_right = count - n_left;

        /* Safety: if partition is degenerate, make leaf */
        if (n_left == 0 || n_right == 0) {
            tree->split_var[node] = -1;
            tree->leaf_offset[node] = leaf_buf_len;
            tree->leaf_count[node] = count;
            if (leaf_buf_len + count > leaf_buf_cap) {
                leaf_buf_cap = (leaf_buf_len + count) * 2;
                int *tmp = (int *)realloc(leaf_buf,
                                          (size_t)leaf_buf_cap * sizeof(int));
                if (!tmp) { free(queue); free(leaf_buf); return -1; }
                leaf_buf = tmp;
            }
            for (int i = 0; i < count; i++)
                leaf_buf[leaf_buf_len++] = samples[start + i];
            continue;
        }

        /* Create child nodes */
        int left_id = drf_tree_add_node(tree);
        int right_id = drf_tree_add_node(tree);
        if (left_id < 0 || right_id < 0) {
            free(queue); free(leaf_buf); return -1;
        }

        tree->left_child[node] = left_id;
        tree->right_child[node] = right_id;

        /* Enqueue children */
        if (q_tail + 2 > queue_cap) {
            queue_cap *= 2;
            BFSItem *tmp = (BFSItem *)realloc(queue,
                                              (size_t)queue_cap * sizeof(BFSItem));
            if (!tmp) { free(queue); free(leaf_buf); return -1; }
            queue = tmp;
        }
        queue[q_tail].node_id = left_id;
        queue[q_tail].start   = start;
        queue[q_tail].count   = n_left;
        q_tail++;

        queue[q_tail].node_id = right_id;
        queue[q_tail].start   = start + n_left;
        queue[q_tail].count   = n_right;
        q_tail++;
    }

    /* Finalize leaf_samples */
    tree->leaf_samples = (int *)malloc((size_t)leaf_buf_len * sizeof(int));
    if (!tree->leaf_samples && leaf_buf_len > 0) {
        free(queue); free(leaf_buf); return -1;
    }
    if (leaf_buf_len > 0)
        memcpy(tree->leaf_samples, leaf_buf, (size_t)leaf_buf_len * sizeof(int));

    free(queue);
    free(leaf_buf);
    return 0;
}

/* ================================================================
 * Section 11: Honesty -- Re-populate Leaves
 * ================================================================
 *
 * After building the tree on the "growing" subsample, push the
 * "honest" subsample through and replace leaf contents with honest
 * sample indices.
 */

static int traverse_tree(const DRFTree *tree, const double *X, int p,
                         int obs_idx)
{
    int node = 0; /* root */
    while (tree->left_child[node] != 0) {
        int var = tree->split_var[node];
        double val = tree->split_value[node];
        if (X[obs_idx * p + var] <= val) {
            node = tree->left_child[node];
        } else {
            node = tree->right_child[node];
        }
    }
    return node;
}

static int honesty_repopulate(DRFTree *tree, const int *honest_samples,
                              int n_honest, const double *X, int p)
{
    /* For each leaf node, collect honest samples that land there. */
    /* First pass: count samples per leaf */
    int *leaf_new_count = (int *)calloc((size_t)tree->n_nodes, sizeof(int));
    int *obs_leaf = (int *)malloc((size_t)n_honest * sizeof(int));
    if (!leaf_new_count || !obs_leaf) {
        free(leaf_new_count); free(obs_leaf);
        return -1;
    }

    for (int i = 0; i < n_honest; i++) {
        int leaf = traverse_tree(tree, X, p, honest_samples[i]);
        obs_leaf[i] = leaf;
        leaf_new_count[leaf]++;
    }

    /* Compute new offsets */
    int total = 0;
    for (int nd = 0; nd < tree->n_nodes; nd++) {
        if (tree->left_child[nd] == 0) { /* leaf */
            tree->leaf_offset[nd] = total;
            tree->leaf_count[nd] = leaf_new_count[nd];
            total += leaf_new_count[nd];
        } else {
            tree->leaf_offset[nd] = 0;
            tree->leaf_count[nd] = 0;
        }
    }

    /* Allocate new leaf_samples */
    free(tree->leaf_samples);
    tree->leaf_samples = NULL;
    if (total > 0) {
        tree->leaf_samples = (int *)malloc((size_t)total * sizeof(int));
        if (!tree->leaf_samples) {
            free(leaf_new_count); free(obs_leaf);
            return -1;
        }
    }

    /* Fill using a write cursor per leaf */
    int *cursor = (int *)calloc((size_t)tree->n_nodes, sizeof(int));
    if (!cursor) { free(leaf_new_count); free(obs_leaf); return -1; }

    for (int i = 0; i < n_honest; i++) {
        int leaf = obs_leaf[i];
        int pos = tree->leaf_offset[leaf] + cursor[leaf];
        tree->leaf_samples[pos] = honest_samples[i];
        cursor[leaf]++;
    }

    free(leaf_new_count);
    free(obs_leaf);
    free(cursor);
    return 0;
}

/* ================================================================
 * Section 12: Forest Training
 * ================================================================ */

static int train_forest(
    DRFTree *forest,       /* pre-allocated array of n_trees trees */
    int n_trees,
    const double *Y_scaled,
    const double *Y_orig,
    const double *X,
    int n, int p,
    int mtry,
    int min_node_size,
    double bandwidth,
    double sample_fraction,
    double alpha,
    int num_features,
    int honesty,
    double honesty_fraction,
    uint64_t base_seed)
{
    int n_subsample = (int)(n * sample_fraction);
    if (n_subsample < 2) n_subsample = 2;
    if (n_subsample > n) n_subsample = n;

    int tree_capacity = 2 * n_subsample + 2;

    char msg[256];

    for (int t = 0; t < n_trees; t++) {

        /* Progress display */
        if (t % 100 == 0) {
            snprintf(msg, sizeof(msg),
                     "DRF: training tree %d / %d\n", t + 1, n_trees);
            SF_display(msg);
        }

        /* Seed RNG for this tree */
        xorshift128p_state rng;
        xorshift128p_seed(&rng, base_seed + (uint64_t)t);

        /* Initialize tree */
        if (drf_tree_init(&forest[t], tree_capacity, n) != 0) {
            snprintf(msg, sizeof(msg),
                     "DRF error: memory allocation failed for tree %d\n", t);
            SF_error(msg);
            return -1;
        }

        /* Create index array and shuffle for subsampling */
        int *all_idx = (int *)malloc((size_t)n * sizeof(int));
        if (!all_idx) {
            SF_error("DRF error: memory allocation failed for index array\n");
            return -1;
        }
        for (int i = 0; i < n; i++) all_idx[i] = i;
        fisher_yates_shuffle(all_idx, n, &rng);

        /* Mark in-bag observations */
        for (int i = 0; i < n_subsample; i++) {
            forest[t].in_bag[all_idx[i]] = 1;
        }

        if (honesty) {
            /* Split subsample into growing and honest halves */
            int n_grow = (int)(n_subsample * honesty_fraction);
            if (n_grow < 2) n_grow = 2;
            int n_honest = n_subsample - n_grow;
            if (n_honest < 1) n_honest = 1;
            n_grow = n_subsample - n_honest;

            int *grow_samples = (int *)malloc((size_t)n_grow * sizeof(int));
            int *honest_samples = (int *)malloc((size_t)n_honest * sizeof(int));
            if (!grow_samples || !honest_samples) {
                free(all_idx); free(grow_samples); free(honest_samples);
                SF_error("DRF error: memory allocation failed for honesty split\n");
                return -1;
            }
            memcpy(grow_samples, all_idx, (size_t)n_grow * sizeof(int));
            memcpy(honest_samples, all_idx + n_grow,
                   (size_t)n_honest * sizeof(int));

            /* Build tree on growing half */
            int rc = build_tree(&forest[t], grow_samples, n_grow,
                                Y_scaled, X, p, mtry, num_features,
                                bandwidth, min_node_size, alpha, &rng);
            if (rc != 0) {
                free(all_idx); free(grow_samples); free(honest_samples);
                snprintf(msg, sizeof(msg),
                         "DRF error: tree build failed for tree %d\n", t);
                SF_error(msg);
                return -1;
            }

            /* Re-populate leaves with honest samples */
            rc = honesty_repopulate(&forest[t], honest_samples, n_honest,
                                    X, p);
            if (rc != 0) {
                free(all_idx); free(grow_samples); free(honest_samples);
                snprintf(msg, sizeof(msg),
                         "DRF error: honesty repopulate failed for tree %d\n", t);
                SF_error(msg);
                return -1;
            }

            free(grow_samples);
            free(honest_samples);
        } else {
            /* No honesty: build tree on full subsample */
            /* Make a working copy since build_tree partitions in-place */
            int *work = (int *)malloc((size_t)n_subsample * sizeof(int));
            if (!work) {
                free(all_idx);
                SF_error("DRF error: memory allocation failed\n");
                return -1;
            }
            memcpy(work, all_idx, (size_t)n_subsample * sizeof(int));

            int rc = build_tree(&forest[t], work, n_subsample,
                                Y_scaled, X, p, mtry, num_features,
                                bandwidth, min_node_size, alpha, &rng);
            free(work);
            if (rc != 0) {
                free(all_idx);
                snprintf(msg, sizeof(msg),
                         "DRF error: tree build failed for tree %d\n", t);
                SF_error(msg);
                return -1;
            }
        }

        free(all_idx);
    }

    /* Final progress message */
    snprintf(msg, sizeof(msg), "DRF: all %d trees trained.\n", n_trees);
    SF_display(msg);

    return 0;
}

/* ================================================================
 * Section 13: OOB Prediction
 * ================================================================
 *
 * For each observation i:
 *   1. Find all trees where in_bag[i] == 0.
 *   2. Traverse each such tree to find i's leaf.
 *   3. For each training sample j in that leaf, accumulate
 *      weight[j] += 1.0 / leaf_count.
 *   4. Normalize weights.
 *   5. Apply functional: mean or weighted quantile.
 */

static int oob_predict(
    const DRFTree *forest,
    int n_trees,
    const double *Y_orig,
    const double *X,
    int n, int p,
    int functional_type,   /* 0 = mean, 1 = quantile */
    double quantile_prob,
    double *predictions)   /* output: n predictions */
{
    /* Weight buffer: reused per observation */
    double *weights = (double *)calloc((size_t)n, sizeof(double));
    /* For quantile: need arrays of (y, w) pairs with nonzero weight */
    double *q_y = NULL;
    double *q_w = NULL;
    if (functional_type == 1) {
        q_y = (double *)malloc((size_t)n * sizeof(double));
        q_w = (double *)malloc((size_t)n * sizeof(double));
        if (!q_y || !q_w) {
            free(weights); free(q_y); free(q_w);
            return -1;
        }
    }

    if (!weights) { free(q_y); free(q_w); return -1; }

    char msg[256];
    int n_no_oob = 0;

    for (int i = 0; i < n; i++) {

        if (i % 5000 == 0 && i > 0) {
            snprintf(msg, sizeof(msg),
                     "DRF: OOB prediction for observation %d / %d\n",
                     i, n);
            SF_display(msg);
        }

        /* Reset weights */
        memset(weights, 0, (size_t)n * sizeof(double));

        double total_weight = 0.0;
        int n_oob_trees = 0;

        for (int t = 0; t < n_trees; t++) {
            if (forest[t].in_bag[i]) continue; /* skip in-bag trees */
            n_oob_trees++;

            /* Traverse tree to find leaf */
            int leaf = traverse_tree(&forest[t], X, p, i);
            int lc = forest[t].leaf_count[leaf];
            if (lc == 0) continue;

            double w = 1.0 / (double)lc;
            int off = forest[t].leaf_offset[leaf];
            for (int k = 0; k < lc; k++) {
                int j = forest[t].leaf_samples[off + k];
                weights[j] += w;
                total_weight += w;
            }
        }

        if (n_oob_trees == 0 || total_weight <= 0.0) {
            /* No OOB trees for this observation: set missing */
            predictions[i] = SV_missval;
            n_no_oob++;
            continue;
        }

        /* Normalize weights */
        for (int j = 0; j < n; j++) {
            weights[j] /= total_weight;
        }

        if (functional_type == 0) {
            /* Weighted mean */
            double pred = 0.0;
            for (int j = 0; j < n; j++) {
                if (weights[j] > 0.0) {
                    pred += weights[j] * Y_orig[j];
                }
            }
            predictions[i] = pred;
        } else {
            /* Weighted quantile */
            int nz = 0;
            for (int j = 0; j < n; j++) {
                if (weights[j] > 0.0) {
                    q_y[nz] = Y_orig[j];
                    q_w[nz] = weights[j];
                    nz++;
                }
            }
            predictions[i] = weighted_quantile(q_y, q_w, nz, quantile_prob);
        }
    }

    if (n_no_oob > 0) {
        snprintf(msg, sizeof(msg),
                 "DRF warning: %d observations had no OOB trees "
                 "(set to missing).\n", n_no_oob);
        SF_display(msg);
    }

    free(weights);
    free(q_y);
    free(q_w);
    return 0;
}

/* ================================================================
 * Section 14: stata_call() -- Main Entry Point
 * ================================================================
 *
 * argv[] layout:
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

STDLL stata_call(int argc, char *argv[])
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

    uint64_t base_seed = (uint64_t)seed;

    /* ----------------------------------------------------------
     * Step 1: Read dimensions
     * ---------------------------------------------------------- */
    int nvar = SF_nvars();  /* variables in the plugin varlist, not the whole dataset */
    int p = nvar - 2;   /* number of X variables */
    if (p < 1) {
        SF_error("DRF error: need at least 1 predictor variable.\n");
        return 198;
    }

    /* Count usable observations (non-missing, in if/in range) */
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

    /* Auto mtry: min(ceil(sqrt(p) + 20), p) */
    if (mtry <= 0) {
        mtry = (int)ceil(sqrt((double)p) + 20.0);
        if (mtry > p) mtry = p;
    }
    if (mtry > p) mtry = p;

    snprintf(msg, sizeof(msg),
             "DRF: n=%d, p=%d, trees=%d, mtry=%d, min_node=%d, "
             "features=%d, honesty=%d\n",
             n, p, n_trees, mtry, min_node_size, num_features, honesty);
    SF_display(msg);

    /* ----------------------------------------------------------
     * Step 2: Read data from Stata into C arrays
     * ---------------------------------------------------------- */
    double *Y_orig   = (double *)malloc((size_t)n * sizeof(double));
    double *Y_scaled = (double *)malloc((size_t)n * sizeof(double));
    double *X        = (double *)malloc((size_t)n * p * sizeof(double));
    int    *obs_map  = (int *)malloc((size_t)n * sizeof(int));

    if (!Y_orig || !Y_scaled || !X || !obs_map) {
        SF_error("DRF error: memory allocation failed for data arrays.\n");
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        return 2000;
    }

    int idx = 0;
    int any_missing = 0;
    for (ST_int i = obs1; i <= obs2; i++) {
        if (!SF_ifobs(i)) continue;

        /* Read Y (variable 1, Stata is 1-indexed) */
        double val;
        ST_retcode rc = SF_vdata(1, i, &val);
        if (rc || SF_is_missing(val)) {
            any_missing = 1;
            continue;
        }
        Y_orig[idx] = val;

        /* Read X (variables 2 .. nvar-1) */
        int skip = 0;
        for (int j = 0; j < p; j++) {
            rc = SF_vdata(j + 2, i, &val);
            if (rc || SF_is_missing(val)) { skip = 1; break; }
            X[idx * p + j] = val;
        }
        if (skip) { any_missing = 1; continue; }

        obs_map[idx] = i;  /* map back to Stata obs index */
        idx++;
    }
    n = idx;  /* update to actual count of complete cases */

    if (any_missing) {
        snprintf(msg, sizeof(msg),
                 "DRF: dropped observations with missing values, "
                 "using n=%d complete cases.\n", n);
        SF_display(msg);
    }

    if (n < 2) {
        SF_error("DRF error: fewer than 2 complete observations.\n");
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        return 2000;
    }

    /* ----------------------------------------------------------
     * Step 3: Scale Y (center and standardize)
     * ---------------------------------------------------------- */
    double y_mean = 0.0;
    for (int i = 0; i < n; i++) y_mean += Y_orig[i];
    y_mean /= n;

    double y_var = 0.0;
    for (int i = 0; i < n; i++) {
        double d = Y_orig[i] - y_mean;
        y_var += d * d;
    }
    y_var /= (n - 1);
    double y_sd = sqrt(y_var);
    if (y_sd < 1e-15) y_sd = 1.0;  /* constant Y: avoid division by zero */

    for (int i = 0; i < n; i++) {
        Y_scaled[i] = (Y_orig[i] - y_mean) / y_sd;
    }

    snprintf(msg, sizeof(msg),
             "DRF: Y mean=%.6f, sd=%.6f\n", y_mean, y_sd);
    SF_display(msg);

    /* ----------------------------------------------------------
     * Step 4: Compute bandwidth (median heuristic if auto)
     * ---------------------------------------------------------- */
    xorshift128p_state bw_rng;
    xorshift128p_seed(&bw_rng, base_seed);

    if (bandwidth <= 0.0) {
        bandwidth = median_heuristic_bandwidth(Y_scaled, n, &bw_rng);
        snprintf(msg, sizeof(msg),
                 "DRF: auto bandwidth (median heuristic) = %.6f\n",
                 bandwidth);
        SF_display(msg);
    } else {
        snprintf(msg, sizeof(msg),
                 "DRF: user bandwidth = %.6f\n", bandwidth);
        SF_display(msg);
    }

    /* ----------------------------------------------------------
     * Step 5: Train the forest
     * ---------------------------------------------------------- */
    DRFTree *forest = (DRFTree *)calloc((size_t)n_trees, sizeof(DRFTree));
    if (!forest) {
        SF_error("DRF error: memory allocation failed for forest.\n");
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        return 2000;
    }

    int rc = train_forest(forest, n_trees,
                          Y_scaled, Y_orig, X, n, p,
                          mtry, min_node_size, bandwidth,
                          sample_frac, alpha, num_features,
                          honesty, honesty_frac, base_seed);

    if (rc != 0) {
        /* Cleanup partial forest */
        for (int t = 0; t < n_trees; t++) drf_tree_free(&forest[t]);
        free(forest);
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        return 2000;
    }

    /* ----------------------------------------------------------
     * Step 6: OOB predictions
     * ---------------------------------------------------------- */
    double *predictions = (double *)malloc((size_t)n * sizeof(double));
    if (!predictions) {
        for (int t = 0; t < n_trees; t++) drf_tree_free(&forest[t]);
        free(forest);
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        SF_error("DRF error: memory allocation for predictions failed.\n");
        return 2000;
    }

    SF_display("DRF: computing OOB predictions...\n");

    rc = oob_predict(forest, n_trees, Y_orig, X, n, p,
                     func_type, quant_prob, predictions);

    if (rc != 0) {
        free(predictions);
        for (int t = 0; t < n_trees; t++) drf_tree_free(&forest[t]);
        free(forest);
        free(Y_orig); free(Y_scaled); free(X); free(obs_map);
        SF_error("DRF error: OOB prediction failed.\n");
        return 2000;
    }

    /* ----------------------------------------------------------
     * Step 7: Write predictions back to Stata
     * ---------------------------------------------------------- */
    int out_var = nvar;  /* last variable is the output */
    int n_written = 0;

    for (int i = 0; i < n; i++) {
        ST_retcode src = SF_vstore(out_var, obs_map[i], predictions[i]);
        if (src) {
            snprintf(msg, sizeof(msg),
                     "DRF warning: failed to store prediction "
                     "for obs %d (rc=%d)\n", obs_map[i], src);
            SF_display(msg);
        } else {
            n_written++;
        }
    }

    snprintf(msg, sizeof(msg),
             "DRF: wrote %d OOB predictions to output variable.\n",
             n_written);
    SF_display(msg);

    /* ----------------------------------------------------------
     * Step 8: Save summary scalars for Stata
     * ---------------------------------------------------------- */
    SF_scal_save("e(N)",        (double)n);
    SF_scal_save("e(n_trees)",  (double)n_trees);
    SF_scal_save("e(mtry)",     (double)mtry);
    SF_scal_save("e(bandwidth)",(double)bandwidth);
    SF_scal_save("e(p)",        (double)p);

    /* ----------------------------------------------------------
     * Cleanup
     * ---------------------------------------------------------- */
    free(predictions);
    for (int t = 0; t < n_trees; t++) drf_tree_free(&forest[t]);
    free(forest);
    free(Y_orig);
    free(Y_scaled);
    free(X);
    free(obs_map);

    SF_display("DRF: done.\n");
    return 0;
}

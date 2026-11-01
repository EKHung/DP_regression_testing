# Run from this folder: Rscript run_sensitivity.R
# Output: sensitivity.Rdata, containing sensitivity$bandwidth, $blocks,
# $constants and $coordinate. Each has a summary and error matrices.
source('test_functions.R')
cfg <- default_settings()
script_dir <- getwd()

# Settings: change these directly. For a quick execution check, use n_rep <- 2.
n <- 4000
d <- 500
epsilon <- 1
n_rep <- 250                    # 250 null AND 250 alternative repetitions
seed <- 1
set.seed(seed)                 # Set once; all random draws follow sequentially.
x_bandwidths <- c(0.5, 1, 2)     # Multipliers of sqrt(d)
y_bandwidths <- c(0.5, 1, 2, 4)
block_counts <- seq_len(choose_blocks(n,d))
do_bandwidth <- TRUE
do_blocks <- TRUE
do_constants <- TRUE
do_coordinate <- TRUE

# Transformed-mean constants sensitivity (fixed sample size).
cfg$constant_n <- 4000L
cfg$constant_d <- c(Gaussian=500L,Exponential=50L)
cfg$constant_epsilon <- c(Gaussian=1,Exponential=4)
cfg$constant_C_sep <- c(3,5,7)
cfg$constant_C_clip <- c(0.25,0.5,0.75)
cfg$constant_C_tol <- c(0.025,0.05,0.075)

# Gaussian coordinate-test sensitivity, beta=e1; no sparsity-dependent tuning.
cfg$coordinate_sensitivity_n <- 4000L
cfg$coordinate_sensitivity_d <- 500L
cfg$coordinate_sensitivity_epsilon <- 1
cfg$coordinate_sensitivity_c_B <- c(0.1,0.25,1,5,10)
cfg$coordinate_sensitivity_c_lambda <- c(0.1,0.25,1,5,10)
output_file <- file.path(script_dir,"sensitivity.Rdata")

# These settings are used by the functions loaded from test_functions.R.
cfg$significance_level <- 0.05
cfg$beta_norm <- 1
cfg$permutations <- 99
cfg$subset_alpha <- 0.25
cfg$pb_clip_quantile <- 0.99
cfg$pb_calibration_draws <- 10000

# Only plotting tables and matrices go into this list.
# Each new run starts afresh; this elementary version does not resume old runs.
sensitivity <- list(format_version='regression-v4',settings=cfg,
  bandwidth=NULL,blocks=NULL,constants=NULL,coordinate=NULL)

# Summarise individual decisions, using the same MCSE as the main script.
plot_summary <- function(decisions, columns) {
  s <- summarise_by(decisions, columns)
  s[, !names(s) %in% c("null_rejections", "alt_rejections"), drop = FALSE]
}

# 1. dpHSIC bandwidth sensitivity.
# Within each trial, bandwidth pairs share data but draw fresh permutations and noise.
if (do_bandwidth) {
  bandwidth_decisions <- NULL

  for (r in seq_len(n_rep)) {
    for (h in 0:1) {             # h=0: beta=0; h=1: beta=e1
      X <- matrix(rnorm(n * d), n, d)
      data <- list(X = X, Y = h * cfg$beta_norm * X[, 1] + rnorm(n))
      rm(X)
      G <- tcrossprod(data$X)

      decisions <- hsic_bandwidth_grid(
        G, data$Y, d, epsilon, x_bandwidths, y_bandwidths
      )
      decisions$n <- n
      decisions$d <- d
      decisions$epsilon <- epsilon
      decisions$rep <- r
      decisions$h <- h
      bandwidth_decisions <- rbind(bandwidth_decisions, decisions)
      rm(data, G)
    }

    # Save after each completed null/alternative pair.
    s <- plot_summary(bandwidth_decisions,
      c("n", "d", "epsilon", "x_multiplier", "sigma_X", "sigma_Y"))
    sensitivity$bandwidth <- list(
      summary = s,
      total_error_matrix = xtabs(total_error ~ x_multiplier + sigma_Y, data = s),
      power_matrix = xtabs(power ~ x_multiplier + sigma_Y, data = s),
      type_I_matrix = xtabs(type_I ~ x_multiplier + sigma_Y, data = s),
      type_II_matrix = xtabs(type_II ~ x_multiplier + sigma_Y, data = s)
    )
    save(sensitivity, file = output_file)
  }
  message(sprintf("Bandwidth sensitivity: n=%d complete",n))
}

# 2. ToT-F and PB-LRT block-count sensitivity.
# Fresh datasets are generated for the block-count experiment.
# Each method and block count draws a separate random partition and privacy noise.
if (do_blocks) {
  if (any(floor(n / block_counts) <= d)) {
    stop("Every block must have more than d observations.")
  }

  # Calibrate each PB configuration once, before generating the test datasets.
  # U and the private null critical value depend on the actual block sizes.
  pb_calibration <- vector("list", length(block_counts))
  tot_cutoff <- numeric(length(block_counts))
  for (j in seq_along(block_counts)) {
    m <- block_counts[j]
    sizes <- block_sizes(n, m)
    U <- pb_clip_bound(sizes, d)
    reference <- pb_null_reference(sizes, d, U, epsilon, cfg$pb_calibration_draws)
    pb_calibration[[j]] <- list(diagnostics = list(
      U = U, critical = unname(quantile(reference, 1 - cfg$significance_level, type = 7))
    ))
    tot_cutoff[j] <- tot_critical(m, epsilon)
  }

  block_decisions <- NULL
  for (r in seq_len(n_rep)) {
    for (h in 0:1) {
      X <- matrix(rnorm(n * d), n, d)
      data <- list(X = X, Y = h * cfg$beta_norm * X[, 1] + rnorm(n))
      rm(X)

      for (j in seq_along(block_counts)) {
        m <- block_counts[j]

        fits_tot <- subset_stats(data$X, data$Y, m)

        # ToT-F: count significant subset F-tests, then add Tulap noise.
        votes <- sum(fits_tot[, 1] <= cfg$subset_alpha)

        b <- exp(-epsilon)
        noise <- rgeom(1, 1-b) - rgeom(1, 1-b) + runif(1, -0.5, 0.5)
        reject_tot <- as.integer(votes + noise > tot_cutoff[j])

        # PB-LRT: clip subset log likelihood ratios, average and add Laplace noise.
        cal <- pb_calibration[[j]]$diagnostics

        fits_pb <- subset_stats(data$X, data$Y, m)
        statistic <- mean(pmin(fits_pb[, 2], cal$U)) +
        rlaplace(1, cal$U / (m * epsilon))
        reject_pb <- as.integer(statistic > cal$critical)

        decisions <- data.frame(
          n = n, d = d, epsilon = epsilon, rep = r, h = h, blocks = m,
          method = c("ToT-F", "PB-LRT"), reject = c(reject_tot, reject_pb)
        )
        block_decisions <- rbind(block_decisions, decisions)
      }
      rm(data, fits_tot, fits_pb)
    }

    s <- plot_summary(block_decisions, c("n", "d", "epsilon", "blocks", "method"))
    sensitivity$blocks <- list(
      summary = s,
      total_error_matrix = xtabs(total_error ~ blocks + method, data = s),
      power_matrix = xtabs(power ~ blocks + method, data = s),
      type_I_matrix = xtabs(type_I ~ blocks + method, data = s),
      type_II_matrix = xtabs(type_II ~ blocks + method, data = s)
    )
    save(sensitivity, file = output_file)
  }
  message(sprintf("Block sensitivity: n=%d complete",n))
}

# 3. Transformed-mean constants sensitivity.
if (do_constants) {
  raw <- NULL
  grid <- expand.grid(C_sep=cfg$constant_C_sep,C_clip=cfg$constant_C_clip,C_tol=cfg$constant_C_tol,
    KEEP.OUT.ATTRS=FALSE)
  n <- cfg$constant_n
  for(model in c('Gaussian','Exponential')) {
    d <- unname(cfg$constant_d[model])
    epsilon <- unname(cfg$constant_epsilon[model])
    for (rep in seq_len(n_rep)) {
      for (h in 0:1) {
        if(model=='Gaussian') {
          data <- gaussian_data(n,d,h)
          X <- data$X
          Y <- data$Y
          rm(data)
        }else {

          X <- matrix(sample(c(-1,1),n*d,replace=TRUE),n,d)
          linear_predictor <- cfg$exponential_intercept+h*cfg$beta_norm*rowSums(X)/sqrt(d)
          Y <- rexp(n,rate=exp(-linear_predictor))
        }
        G <- tcrossprod(X)
        rows <- lapply(seq_len(nrow(grid)),function(i) {
            g <- grid[i,]
            p <- parameters(n,d,epsilon,g$C_sep,g$C_clip,g$C_tol)

            decision <- transformed_mean(G,Y,d,epsilon,C_sep=g$C_sep,C_clip=g$C_clip,C_tol=g$C_tol)
            data.frame(model=model,n=n,d=d,epsilon=epsilon,rep=rep,h=h,
            C_sep=g$C_sep,C_clip=g$C_clip,C_tol=g$C_tol,eta=p$eta,alpha=p$alpha,J=p$J,r_tilde=p$r_tilde,reject=decision)
        })
        raw <- rbind(raw,do.call(rbind,rows))

        rm(X,Y,G)
        invisible(gc(FALSE))
      }
    }
    summary <- plot_summary(raw,c('model','n','d','epsilon','C_sep','C_clip','C_tol','eta','alpha','J','r_tilde'))
    # One named three-dimensional array per model: [C_sep, C_clip, C_tol].
    # For a table at C_tol=0.05, select array[, , "0.05"].
    # Keeping C_tol as a dimension prevents summing errors across its values.
    total_error <- lapply(split(summary,summary$model),function(g)
      xtabs(total_error~C_sep+C_clip+C_tol,data=g))
    type_I <- lapply(split(summary,summary$model),function(g)
      xtabs(type_I~C_sep+C_clip+C_tol,data=g))
    type_II <- lapply(split(summary,summary$model),function(g)
      xtabs(type_II~C_sep+C_clip+C_tol,data=g))
    power <- lapply(split(summary,summary$model),function(g)
      xtabs(power~C_sep+C_clip+C_tol,data=g))
    sensitivity$constants <- list(summary=summary, settings=cfg,
      algorithm_version="paper-alpha-log-names-v2",
      total_error=total_error,type_I=type_I,type_II=type_II,power=power)
    save(sensitivity,file=output_file)
    message(sprintf('Constants: %s, n=%d complete',model,n))
  }
}

# 4. Coordinate-test sensitivity: same data, fresh Laplace noise per grid cell.
# Vary working C_X=c_B*sqrt(8/3) and lambda_bar=c_lambda, keeping N(0,I_d) fixed.
# Changing c_B recomputes B, b_B, q_B and privacy noise; c_lambda changes only q_B.
# Values below valid bounds assess misspecification; every fixed setting is DP.
if(do_coordinate) {
  coordinate_raw <- NULL
  coordinate_n <- cfg$coordinate_sensitivity_n
  coordinate_d <- cfg$coordinate_sensitivity_d
  coordinate_epsilon <- cfg$coordinate_sensitivity_epsilon
  for(r in seq_len(n_rep)) {
    for(h in 0:1) {
      data <- gaussian_data(coordinate_n,coordinate_d,h)
      decisions <- coordinate_grid(data$X,data$Y,coordinate_epsilon,
        c_B=cfg$coordinate_sensitivity_c_B,
        c_lambda=cfg$coordinate_sensitivity_c_lambda)
      decisions$n <- coordinate_n
      decisions$d <- coordinate_d
      decisions$epsilon <- coordinate_epsilon
      decisions$sparsity <- 1L
      decisions$rep <- r
      decisions$h <- h
      coordinate_raw <- rbind(coordinate_raw,decisions)
      rm(data)
    }
    summary <- plot_summary(coordinate_raw,c('n','d','epsilon','sparsity',
      'c_B','c_lambda','C_X_working','lambda_bar_working',
      'B','b_B','q_B','noise_scale','threshold'))
    summary$coordinate_threshold <- 'Bernstein'
    sensitivity$coordinate <- list(summary=summary,settings=cfg,
      total_error_matrix=xtabs(total_error~c_lambda+c_B,data=summary),
      type_I_matrix=xtabs(type_I~c_lambda+c_B,data=summary),
      type_II_matrix=xtabs(type_II~c_lambda+c_B,data=summary),
      power_matrix=xtabs(power~c_lambda+c_B,data=summary))
    save(sensitivity,file=output_file)
  }
  message(sprintf('Coordinate sensitivity: n=%d complete',coordinate_n))
}

# Reload later with load("sensitivity.Rdata"). For example:
# sensitivity$bandwidth$total_error_matrix
# sensitivity$blocks$summary
# sensitivity$constants$summary  # All 27 combinations per model when complete.
# sensitivity$constants$total_error$Gaussian[, , "0.05"]
# sensitivity$coordinate$summary              # 5-by-5 grid, s*=1, n=4000
# sensitivity$coordinate$type_I_matrix        # rows c_lambda, columns c_B
# sensitivity$coordinate$type_II_matrix
# sensitivity$coordinate$total_error_matrix
# For a separate table with C_tol columns and (C_sep,C_clip) rows:
# xtabs(total_error ~ C_sep + C_clip + C_tol,
#       data=subset(sensitivity$constants$summary, model=="Gaussian"))
# No individual decisions or calibration records are saved in sensitivity.Rdata.
cat("Plotting sensitivity saved in", output_file, "\n")

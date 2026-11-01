# Run from this folder: Rscript run_gaussian.R
# Output: gauss.RData, containing gauss$summary (five methods).
source('test_functions.R')
cfg <- default_settings()
cfg$repetitions <- 250L
set.seed(cfg$seed)

run_comparison <- function() {
  raw <- NULL
  gauss <- list(format_version="regression-v4",settings=cfg,summary=NULL)
  d <- cfg$comparison_d
  epsilon <- cfg$comparison_epsilon
  for(n in cfg$comparison_n) {
    blocks <- choose_blocks(n,d) # Same block count; each method draws its own partition
    cutoff_tot <- tot_critical(blocks,epsilon)
    cal <- calibrate_pb(n,d,blocks,epsilon)
    U <- cal$diagnostics$U
    cutoff_pb <- cal$diagnostics$critical
    for (rep in seq_len(cfg$repetitions)) {
      for (h in 0:1) {
        data <- gaussian_data(n,d,h)
        X <- data$X
        Y <- data$Y
        G <- tcrossprod(X)

        reject_tm <- transformed_mean(G,Y,d,epsilon)

        default_hsic <- hsic_bandwidth_grid(G,Y,d,epsilon,1,1)

        # Each aggregation method draws its own random block partition.

        fits_tot <- subset_stats(X,Y,blocks)
        votes <- sum(fits_tot[,1L]<=cfg$subset_alpha)

        b <- exp(-epsilon)
        tulap_noise <- rgeom(1,1-b)-rgeom(1,1-b)+runif(1,-.5,.5)
        reject_tot <- as.integer(votes+tulap_noise>cutoff_tot)

        fits_pb <- subset_stats(X,Y,blocks)
        pb_statistic <- mean(pmin(fits_pb[,2L],U))+rlaplace(1,U/(blocks*epsilon))
        reject_pb <- as.integer(pb_statistic>cutoff_pb)
        rows <- data.frame(n=n,d=d,rep=rep,h=h,epsilon=epsilon,
          method=c('Transformed mean','dpHSIC','ToT-F','PB-LRT','Coordinate test'),
          blocks=c(NA,NA,blocks,blocks,NA),
          reject=c(reject_tm,default_hsic$reject,reject_tot,reject_pb,
          coordinate_test(X,Y,epsilon)))
        raw <- rbind(raw,rows)

        rm(data,X,Y,G,fits_tot,fits_pb)
        invisible(gc(FALSE))
      }
    }
    summary <- summarise_by(raw,c('n','d','epsilon','method'))
    own <- summary$method=='Transformed mean'
    for(key in c('C_sep','C_clip','C_tol','eta'))summary[[key]] <- ifelse(own,cfg[[key]],NA_real_)
    for(key in c('J','r_tilde','alpha'))summary[[key]] <- vapply(seq_len(nrow(summary)),function(i)
      if(own[i])parameters(summary$n[i],d,epsilon,cfg$C_sep,cfg$C_clip,cfg$C_tol)[[key]] else NA_real_,numeric(1))
    coordinate <- summary$method=='Coordinate test'
    summary$eta[coordinate] <- cfg$eta
    for(key in c('coordinate_c_B','coordinate_c_lambda','coordinate_C_X','coordinate_lambda_bar'))
    summary[[key]] <- ifelse(coordinate,cfg[[key]],NA_real_)
    summary$coordinate_threshold <- ifelse(coordinate,'Bernstein',NA_character_)
    summary$blocks <- ifelse(summary$method %in% c('ToT-F','PB-LRT'),floor(summary$n/(d+1)),NA)
    gauss$summary <- summary
    save(gauss,file=file.path(cfg$output,'gauss.RData'))
    message(sprintf('Gaussian comparison: n=%d complete',n))
  }
  message("Gaussian comparison saved in gauss.RData")
  invisible(gauss)
}


gauss <- run_comparison()

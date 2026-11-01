# Run from this folder: Rscript run_illustrations.R
# Output: illustrations.RData, object illustrations with three components:
# distributions (Cauchy and Exponential), response_dimension, coordinate_sparsity.
source('test_functions.R')
cfg <- default_settings()
cfg$repetitions <- 250L
cfg$response_n <- 2000L
cfg$response_dimensions <- c(50L,100L,500L,1000L,2000L,5000L,10000L)
cfg$response_epsilons <- c(0.1,0.15,0.2,0.3)
run_distributions <- TRUE
run_response_dimension <- TRUE
run_sparsity <- TRUE
set.seed(cfg$seed)
illustrations <- list(format_version='regression-v4',settings=cfg,
                      distributions=NULL,response_dimension=NULL,coordinate_sparsity=NULL)
save_illustration <- function(name,summary) {
  illustrations[[name]] <<- list(summary=summary)
  save(illustrations,file=file.path(cfg$output,'illustrations.RData'))
}

run_illustrations <- function() {
  raw <- NULL
  for(model in c('Cauchy','Exponential')) {
    d <- unname(cfg$illustration_d[model])
    gamma <- if(model=='Exponential')cfg$exponential_intercept else 0
    for (n in cfg$illustration_n) {
      for (rep in seq_len(cfg$repetitions)) {
        for (h in 0:1) {

          X <- matrix(sample(c(-1,1),n*d,replace=TRUE),n,d)
          linear_predictor <- gamma+h*cfg$beta_norm*rowSums(X)/sqrt(d)
          # Exponential rate exp(-linear_predictor) gives mean exp(linear_predictor).
          Y <- if(model=='Cauchy')linear_predictor+rcauchy(n)else rexp(n,rate=exp(-linear_predictor))
          G <- tcrossprod(X)
          rows <- lapply(cfg$illustration_epsilon,function(ep) {

              data.frame(model=model,n=n,d=d,rep=rep,h=h,epsilon=ep,
              intercept=gamma,beta_norm=h*cfg$beta_norm,reject=transformed_mean(G,Y,d,ep))
          })
          raw <- rbind(raw,do.call(rbind,rows))

          rm(X,Y,G)
          invisible(gc(FALSE))
          if(rep==cfg$repetitions && h==1L)
          message(sprintf('%s illustration: n=%d complete',model,n))

        }
      }
    }
    summary <- summarise_by(raw,c('model','n','d','epsilon','intercept'))
    for(key in c('C_sep','C_clip','C_tol','eta'))summary[[key]] <- cfg[[key]]
    for(key in c('J','r_tilde','alpha'))summary[[key]] <- vapply(seq_len(nrow(summary)),function(i)
      parameters(summary$n[i],summary$d[i],summary$epsilon[i],cfg$C_sep,cfg$C_clip,cfg$C_tol)[[key]],numeric(1))
    save_illustration('distributions',summary)
  }
  message("Response-model illustrations complete; plotting results saved.")
}

run_coordinate_sparsity <- function() {
  raw <- NULL
  d <- cfg$sparsity_d
  epsilon <- cfg$sparsity_epsilon
  levels <- cfg$sparsity_levels
  stopifnot(length(levels)>0,all(levels>=1),all(levels<=d),!anyDuplicated(levels))
  for(n in cfg$sparsity_n) {
    for (rep in seq_len(cfg$repetitions)) {
      for (h in 0:1) {
        X <- matrix(rnorm(n*d),n,d)
        error <- rnorm(n)
        if(h==0L) {
          reject <- vapply(levels,function(s)coordinate_test(X,error,epsilon),integer(1))
        }else {
          reject <- vapply(levels,function(s) {
              Y <- cfg$beta_norm*rowSums(X[,seq_len(s),drop=FALSE])/sqrt(s)+error
              coordinate_test(X,Y,epsilon)
          },integer(1))
        }
        raw <- rbind(raw,data.frame(n=n,d=d,epsilon=epsilon,rep=rep,h=h,
          sparsity=levels,beta_norm=h*cfg$beta_norm,reject=reject))
        rm(X,error)
      }
    }
    summary <- summarise_by(raw,c('n','d','epsilon','sparsity'))
    summary$coordinate_c_B <- cfg$coordinate_c_B
    summary$coordinate_c_lambda <- cfg$coordinate_c_lambda
    summary$coordinate_C_X <- cfg$coordinate_C_X
    summary$coordinate_lambda_bar <- cfg$coordinate_lambda_bar
    summary$coordinate_threshold <- 'Bernstein'
    summary$eta <- cfg$eta
    summary$beta_norm <- cfg$beta_norm
    save_illustration('coordinate_sparsity',summary)
    message(sprintf('Coordinate sparsity curves: n=%d complete',n))
  }
}

# Same exponential rate exp(-(gamma+X^T beta)) as the transformed-mean example.
# Exact binomial projection avoids allocating an n-by-d Rademacher matrix.
run_dimension_study <- function() {
  raw <- NULL
  n <- cfg$response_n
  for(d in cfg$response_dimensions) {
    for(rep in seq_len(cfg$repetitions)) {
      for(h in 0:1) {
        predictor <- if(h==0L)rep.int(0,n) else
          cfg$beta_norm*(2*rbinom(n,d,0.5)-d)/sqrt(d)
        log_y <- log(rexp(n))+cfg$exponential_intercept+predictor
        statistic <- response_statistic(log_y)
        for(epsilon in cfg$response_epsilons) {
          reject <- response_test(log_y,epsilon,eta=cfg$eta,statistic=statistic)
          raw <- rbind(raw,data.frame(n=n,d=d,epsilon=epsilon,rep=rep,h=h,reject=reject))
        }
      }
    }
    summary <- summarise_by(raw,c('n','d','epsilon'))
    save_illustration('response_dimension',summary)
    message(sprintf('Response-only dimension study: d=%d complete',d))
  }
}

if(run_distributions)run_illustrations()
if(run_response_dimension)run_dimension_study()
if(run_sparsity)run_coordinate_sparsity()

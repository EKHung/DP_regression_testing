# Shared test functions, summaries and default settings. No simulations run here.
# Rcpp and a C++ compiler are required for the mean and competitor kernels.
# Runners define cfg <- default_settings() before calling the tests.

default_settings <- function() {
  list(
  seed = 1L,
  output = getwd(),
  repetitions = 250L,             # per hypothesis, per sample size
  significance_level = 0.05,     # competitor nominal level
  eta = 0.05,                    # own-test failure parameter, not the empirical level
  beta_norm = 1,

  # Gaussian competitor comparison: X~N(0,I), error~N(0,1), beta=e1.
  comparison_n = seq(2000L,5000L,300L),
  comparison_d = 500L,
  comparison_epsilon = 1,

  # BOTH ToT-F and PB-LRT use the largest feasible integer block count.
  # Every actual block must contain at least d+1 observations.
  # The count is determined by n,d, without inspecting simulated power.
  subset_alpha = 0.25,

  # PB clipping is our specified rule, not a universal default in the paper.
  pb_clip_quantile = 0.99,
  pb_calibration_draws = 10000L,

  # Fixed main dpHSIC kernels: sigma_X=sqrt(d), sigma_Y=1.
  permutations = 99L,

  # Transformed-mean constants; internal epsilon/(6J), coefficients 4,57.
  C_sep = 5,
  C_clip = 0.5,
  C_tol = 0.05,

  # Coordinate test: Bernstein threshold, with no empirical threshold multiplier.
  # c_B and c_lambda scale working distributional bounds in sensitivity checks.
  coordinate_c_B = 1,
  coordinate_c_lambda = 1,
  coordinate_C_X = sqrt(8/3),    # Gaussian psi_2 norm: E exp(X_j^2/C_X^2)=2
  coordinate_lambda_bar = 1,     # Gaussian covariance bound

  # Four Gaussian sparsity curves; ||beta||_2 remains cfg$beta_norm.
  sparsity_n = seq(300L,3000L,300L),
  sparsity_levels = c(1L,5L,10L,20L),
  sparsity_d = 500L,
  sparsity_epsilon = 1,

  # Own-test-only illustrations: dense beta, independent Rademacher X.
  illustration_n = seq(2000L,5000L,300L),
  illustration_epsilon = c(1,2,3,4),
  illustration_d = c(Cauchy=500L,Exponential=50L),
  exponential_intercept = 0.2

)
}

if(!requireNamespace('Rcpp',quietly=TRUE))
stop('Install Rcpp first: install.packages("Rcpp"). A C++ compiler is also required.')

# 2. Compiled pairwise arithmetic (exact, without kernel approximation).
# Exact pairwise calculations in C++; no random-feature approximation.
Rcpp::sourceCpp(code='
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector clipped_rows_cpp(NumericMatrix G, NumericVector t, double R) {
 int n=G.nrow(); NumericVector rows(n);
 for(int j=1;j<n;++j) for(int i=0;i<j;++i) {
  double v=std::max(-R,std::min(R,G(i,j)*t[i]*t[j]));
  rows[i]+=v; rows[j]+=v;
 }
 for(int i=0;i<n;++i) rows[i]/=(n-1.0);
 return rows;
}
// [[Rcpp::export]]
double weighted_cpp(NumericMatrix G, NumericVector t, NumericVector w,
                    double R, double U) {
 int n=G.nrow(); double z=0;
 for(int j=1;j<n;++j) for(int i=0;i<j;++i) {
  double v=std::max(-R,std::min(R,G(i,j)*t[i]*t[j]));
  z+=std::min(w[i],w[j])*(v-U);
 }
 return U+2*z/(double(n)*(n-1.0));
}
// [[Rcpp::export]]
NumericVector perm_sums_cpp(NumericMatrix K, NumericMatrix L, IntegerMatrix P) {
 int n=K.nrow(), B=P.ncol(); NumericVector out(B);
 for(int b=0;b<B;++b) {
  double z=0;
  for(int j=0;j<n;++j) {
   int pj=P(j,b); z+=K(j,j)*L(pj,pj);
   for(int i=0;i<j;++i) z+=2*K(i,j)*L(P(i,b),pj);
  }
  out[b]=z; Rcpp::checkUserInterrupt();
 }
 return out;
}
')

# 3. Noise generation and plotting output.
rlaplace <- function(n,scale=1) {
  u <- runif(n,-.5,.5)
  -scale*sign(u)*log1p(-2*abs(u))
}
# Independent null and alternative samples; plug-in Monte Carlo standard error.
error_rates <- function(g) {
  z <- g$reject[g$h==0L]
  o <- g$reject[g$h==1L]
  if(!length(z)||!length(o))stop('Summary requires both hypotheses.')
  type_I <- mean(z)
  power <- mean(o)
  data.frame(null_reps=length(z),alt_reps=length(o),type_I=type_I,
    power=power,type_II=1-power,total_error=type_I+1-power,
    mcse=sqrt(type_I*(1-type_I)/length(z)+power*(1-power)/length(o)))
}
# keys identify a configuration (for example n, epsilon and method).
# Each group contains its null and alternative decisions; no groups are pooled.
summarise_by <- function(raw,keys) {
  labels <- do.call(interaction,c(raw[keys],list(drop=TRUE,lex.order=TRUE)))
  groups <- split(raw,labels)
  ans <- do.call(rbind,lapply(groups,function(g)cbind(g[1,keys,drop=FALSE],error_rates(g))))
  rownames(ans) <- NULL
  ans
}
gaussian_data <- function(n,d,h) {

  X <- matrix(rnorm(n*d),n,d)
  list(X=X,Y=h*cfg$beta_norm*X[,1L]+rnorm(n))
}

# 4. Transformed-mean test: each scale clips the pairwise kernel, computes the
# weighted correction, adds smooth-sensitivity noise, then compares with s^2/2.
parameters <- function(n,d,epsilon,C_sep,C_clip,C_tol,eta=cfg$eta) {
  stopifnot(n>=2, d>=1, epsilon>0, C_sep>0, C_clip>0, C_tol>0, eta>0, eta<0.25)
  r_tilde <- C_sep*d^.25*(n^(-.5)+n^(-.75)/sqrt(epsilon)+1/(n*epsilon)+
    n^(-1.5)/epsilon^1.5)
  if (r_tilde>1) stop('Starting radius exceeds one; no implemented scale is available.')
  J <- max(1L,ceiling(log2(1/r_tilde)))
  s <- r_tilde*2^(0:(J-1L))
  alpha <- eta/(3*J)
  # Paper: l_s = 1 + log((d+1)*(n+1)*(1/s+1)/alpha).
  l_s <- 1+log(d+1)+log(n+1)+log1p(1/s)-log(alpha)
  list(s=s,R=C_clip*(sqrt(d*l_s)+l_s),xi=C_tol*l_s*(s+sqrt(d/n)),
    J=J,r_tilde=r_tilde,l_s=l_s,eta=eta,alpha=alpha)
}
smooth_bound <- function(L,n,ch,xi,e) {
  a0 <- 4*xi/n+57*ch/(n*n*e)
  a1 <- 4*xi*e/n+57*ch/n^2
  a2 <- 57*ch/n^2*(e+1/n)
  roots <- polyroot(c(a1-e*a0,2*a2-e*a1,-e*a2))
  real <- Re(roots[abs(Im(roots))<1e-8])
  u <- unique(c(0,ceiling(2/e),floor(real-L),ceiling(real-L)))
  u <- u[u>=0 & u<=ceiling(2/e)]
  v <- L+u
  max(exp(-e*u)*(a0+a1*v+a2*v*v))
}
transformed_mean <- function(G,Y,d,epsilon,C_sep=cfg$C_sep,C_clip=cfg$C_clip,C_tol=cfg$C_tol,eta=cfg$eta) {
  n <- length(Y)
  t <- tanh(Y)
  p <- parameters(n,d,epsilon,C_sep,C_clip,C_tol,eta)
  e <- epsilon/(6*p$J) # uncapped experimental implementation, coefficients 4,57
  Q <- sample(c(-1,1),p$J,replace=TRUE)*
  (rgamma(p$J,.25)/rgamma(p$J,.75))^.25
  for (j in seq_len(p$J)) {
    R <- p$R[j]
    xi <- p$xi[j]
    ch <- 2*R
    rows <- clipped_rows_cpp(G,t,R)
    U <- mean(rows)
    dev <- abs(rows-U)
    ordered <- sort(dev,decreasing=TRUE)
    L <- which(c(ordered[-1L],-Inf)<=xi+12*ch*seq_len(n)/n)[1L]
    w <- pmax(0,1-n*e/(12*ch)*pmax(0,dev-xi-12*ch*L/n))
    A <- if (all(w==1)) U else weighted_cpp(G,t,w,R,U)
    if (A+smooth_bound(L,n,ch,xi,e)/e*Q[j]>p$s[j]^2/2) return(1L)
  }
  0L
}

# 4b. Coordinate test: clip X, average X*tanh(Y), take the maximum, add Laplace.
# c_B scales the working C_X in BOTH the clipping level and the bias allowance.
# c_lambda scales the working covariance bound; it changes only q_B.
coordinate_parameters <- function(n,d,epsilon,c_B=cfg$coordinate_c_B,
  c_lambda=cfg$coordinate_c_lambda,eta=cfg$eta,
  C_X=cfg$coordinate_C_X,lambda_bar=cfg$coordinate_lambda_bar) {
  stopifnot(n>=2,d>=1,epsilon>0,c_B>0,c_lambda>0,C_X>0,
    lambda_bar>0,eta>0,eta<0.25)
  C_X_working <- c_B*C_X
  lambda_bar_working <- c_lambda*lambda_bar
  B <- C_X_working*sqrt(log(2*n))
  b_B <- C_X_working^2/B*exp(-B^2/C_X_working^2)
  q_B <- sqrt(2*lambda_bar_working*log(4*d/eta)/n)+4*B*log(4*d/eta)/(3*n)
  noise_scale <- 2*B/(n*epsilon)
  list(C_X_working=C_X_working,lambda_bar_working=lambda_bar_working,
    B=B,b_B=b_B,q_B=q_B,noise_scale=noise_scale,
    threshold=b_B+q_B+noise_scale*log(1/eta))
}
coordinate_statistic <- function(X,Y,B) {
  stopifnot(is.matrix(X),nrow(X)==length(Y),ncol(X)>=1,B>0)
  # Put the matrix first in pmax/pmin to preserve its dimensions in base R.
  clipped_X <- pmin(pmax(X,-B),B)
  max(abs(drop(crossprod(clipped_X,tanh(Y)))/nrow(X)))
}
coordinate_test <- function(X,Y,epsilon,c_B=cfg$coordinate_c_B,
  c_lambda=cfg$coordinate_c_lambda,eta=cfg$eta,
  C_X=cfg$coordinate_C_X,lambda_bar=cfg$coordinate_lambda_bar,noise=NULL) {
  p <- coordinate_parameters(nrow(X),ncol(X),epsilon,c_B=c_B,
    c_lambda=c_lambda,eta=eta,C_X=C_X,lambda_bar=lambda_bar)
  if(is.null(noise))noise <- rlaplace(1)
  as.integer(coordinate_statistic(X,Y,p$B)+p$noise_scale*noise>p$threshold)
}
# Same dataset, with fresh Laplace noise for every sensitivity-grid cell.
# Projection/coordinate clipping is computed only once for each c_B.
coordinate_grid <- function(X,Y,epsilon,c_B,c_lambda,eta=cfg$eta,
  C_X=cfg$coordinate_C_X,lambda_bar=cfg$coordinate_lambda_bar,noise=NULL) {
  stopifnot(length(c_B)>0,length(c_lambda)>0,all(is.finite(c_B)),
    all(is.finite(c_lambda)),all(c_B>0),all(c_lambda>0),
    !anyDuplicated(c_B),!anyDuplicated(c_lambda))
  rows <- lapply(c_B,function(cb) {
      p <- coordinate_parameters(nrow(X),ncol(X),epsilon,c_B=cb,
        c_lambda=1,eta=eta,C_X=C_X,lambda_bar=lambda_bar)
      statistic <- coordinate_statistic(X,Y,p$B)
      cells <- lapply(c_lambda,function(cl) {
        p <- coordinate_parameters(nrow(X),ncol(X),epsilon,c_B=cb,
          c_lambda=cl,eta=eta,C_X=C_X,lambda_bar=lambda_bar)
        cell_noise <- if(is.null(noise))rlaplace(1)else noise
        data.frame(c_B=cb,c_lambda=cl,C_X_working=p$C_X_working,
          lambda_bar_working=p$lambda_bar_working,B=p$B,b_B=p$b_B,q_B=p$q_B,
          noise_scale=p$noise_scale,threshold=p$threshold,
          reject=as.integer(statistic+p$noise_scale*cell_noise>p$threshold))
      })
      do.call(rbind,cells)
  })
  do.call(rbind,rows)
}
# Optional small deterministic check; does not alter the simulation RNG state.
coordinate_self_test <- function() {
  X <- matrix(seq(-3,3,length.out=60),20,3)
  Y <- seq(-1,1,length.out=20)
  p <- coordinate_parameters(20,3,1)
  truth <- max(abs(colSums(sweep(pmin(pmax(X,-p$B),p$B),1,tanh(Y),'*'))/20))
  statistic <- coordinate_statistic(X,Y,p$B)
  X2 <- X
  Y2 <- Y
  X2[1,] <- c(1000,-1000,1000)
  Y2[1] <- -1000
  grid <- coordinate_grid(X,Y,1,c_B=c(.5,1),c_lambda=c(.5,1,2),noise=0)
  stopifnot(abs(statistic-truth)<1e-12,
    abs(statistic-coordinate_statistic(X2,Y2,p$B))<=2*p$B/20+1e-12,
    all(vapply(seq_len(nrow(grid)),function(i)
    coordinate_test(X,Y,1,c_B=grid$c_B[i],c_lambda=grid$c_lambda[i],
      noise=0)==grid$reject[i],logical(1))))
  invisible(TRUE)
}

# 5. ToT-F and subset regression fits.
# Tulap(epsilon,delta=0): difference of geometric variables + Uniform(-.5,.5).
tulap_cdf <- function(x,epsilon) {
  b <- exp(-epsilon)
  k <- floor(x+.5)
  u <- x-k+.5
  ifelse(k<=0,b^(-k)/(1+b)*(b+(1-b)*u),
    1-b^k/(1+b)*(1-(1-b)*u))
}
tot_critical <- function(S,epsilon) {
  pmf <- dbinom(0:S,S,cfg$subset_alpha)
  uniroot(function(z) sum(pmf*tulap_cdf((0:S)-z,epsilon))-cfg$significance_level,
    c(-100,S+100),tol=1e-10)$root
}
block_sizes <- function(n,S) rep(n%/%S,S)+as.integer(seq_len(S)<=n%%S)
# No intercept: the Gaussian null is Y=error and the alternative is Y=X beta+error.
subset_stats <- function(X,Y,S) {
  n <- nrow(X)
  d <- ncol(X)
  sizes <- block_sizes(n,S)
  groups <- split(sample.int(n),rep(seq_len(S),sizes))
  out <- matrix(NA_real_,S,2L)
  for (j in seq_len(S)) {
    ids <- groups[[j]]
    m <- length(ids)
    xb <- X[ids,,drop=FALSE]
    yb <- Y[ids]
    fit <- qr(xb)
    if (fit$rank!=d) { out[j,] <- c(1,0)
      next }
    s0 <- sum(yb^2)
    s1 <- max(sum(qr.resid(fit,yb)^2),.Machine$double.xmin)
    proj <- max(0,s0-s1)
    out[j,1L] <- pf((proj/d)/(s1/(m-d)),d,m-d,lower.tail=FALSE)
    out[j,2L] <- max(0,m*log(s0/s1)/2)
  }
  out
}

# Check the compiled pairwise routines against direct R calculations.
self_test <- function() {
  X <- matrix(rnorm(60),20,3)
  G <- tcrossprod(X)
  t <- tanh(rnorm(20))
  R <- .7
  H <- pmax(pmin(G*outer(t,t),R),-R)
  diag(H) <- 0
  rows <- rowSums(H)/19
  U <- mean(rows)
  w <- runif(20)
  W <- outer(w,w,pmin)
  diag(W) <- 0
  stopifnot(max(abs(rows-clipped_rows_cpp(G,t,R)))<1e-12,
    abs(weighted_cpp(G,t,w,R,U)-(U+sum(W*(H-U))/(20*19)))<1e-12)
  K <- centre(G)
  LL <- centre(outer(t,t))
  P <- cbind(0:19,sample.int(20)-1L)
  truth <- apply(P,2,function(p)sum(K*LL[p+1L,p+1L]))
  stopifnot(max(abs(perm_sums_cpp(K,LL,P)-truth))<1e-9,
    abs(tulap_cdf(0,cfg$comparison_epsilon)-.5)<1e-12)
}

# Same deterministic count rule for both aggregation tests.
choose_blocks <- function(n,d) {
  m <- floor(n/(d+1L))
  if(m<1L)stop('No feasible block count for this n,d.')
  as.integer(m)
}

# 6. PB-LRT: clip subset log likelihood ratios, average, add Laplace noise.
# Conditional on a full-rank design and Gaussian errors under the null,
# Uncentred R^2 in a block of size b is Beta(d/2,(b-d)/2). No covariate-law estimate.
pb_clip_bound <- function(sizes,d) {
  max(-sizes/2*log(qbeta(1-cfg$pb_clip_quantile,(sizes-d)/2,d/2)))
}
pb_null_reference <- function(sizes,d,U,epsilon,draws) {
  values <- numeric(draws)
  M <- length(sizes)
  for(b in sizes) {
    R2 <- rbeta(draws,d/2,(b-d)/2)
    log_LR <- -b*log1p(-R2)/2
    values <- values+pmin(log_LR,U)/M
  }
  values+rlaplace(draws,U/(M*epsilon))
}
calibrate_pb <- function(n,d,M,epsilon) {
  sizes <- block_sizes(n,M)
  U <- pb_clip_bound(sizes,d)

  reference <- pb_null_reference(sizes,d,U,epsilon,cfg$pb_calibration_draws)
  critical <- unname(quantile(reference,1-cfg$significance_level,type=7))

  diagnostics <- data.frame(n=n,d=d,epsilon=epsilon,blocks=M,
    smallest_block=min(sizes),largest_block=max(sizes),U=U,critical=critical,
    calibration_draws=length(reference))
  list(diagnostics=diagnostics,reference=reference)
}

# 7. dpHSIC: fixed Gaussian kernels, square-root biased HSIC, private permutations.
# Kernel bandwidths are sigma_X=x_multiplier*sqrt(d) and sigma_Y=y_bandwidth.
# Every grid cell draws fresh permutations and Laplace noise.
centre <- function(M) {
  sweep(sweep(M,1L,rowMeans(M),'-'),2L,colMeans(M),'-')+mean(M)
}
hsic_bandwidth_grid <- function(G,Y,d,epsilon,x_multipliers,y_bandwidths) {
  n <- length(Y)
  norms <- diag(G)
  distance_X <- pmax(outer(norms,norms,'+')-2*G,0)
  distance_Y <- outer(Y,Y,'-')^2
  rows <- list()
  i <- 0L
  for(mx in x_multipliers) {
    K <- centre(exp(-distance_X/(2*d*mx^2)))
    for(sy in y_bandwidths) {
      P <- cbind(0:(n-1L),replicate(cfg$permutations,sample.int(n)-1L))
      noise <- rlaplace(cfg$permutations+1L,8*(n-1)/(n*n*epsilon))
      L <- centre(exp(-distance_Y/(2*sy^2)))
      statistic <- sqrt(pmax(0,perm_sums_cpp(K,L,P)))/n
      private_statistic <- statistic+noise
      p <- (1+sum(private_statistic[-1L]>=private_statistic[1L]))/length(statistic)
      i <- i+1L
      rows[[i]] <- data.frame(x_multiplier=mx,sigma_X=mx*sqrt(d),sigma_Y=sy,
        p_value=p,reject=as.integer(p<=cfg$significance_level),statistic=statistic[1L])
    }
  }
  do.call(rbind,rows)
}

# Response-only exponential test. Input log_y avoids overflow for large responses.
# Full U-statistic; each unordered pair is counted exactly once.
response_statistic <- function(log_y) {
  n <- length(log_y)
  stopifnot(n>=2,all(is.finite(log_y)))
  z <- sort(log_y)
  pair_count <- sum(findInterval(z-log(2),z,left.open=TRUE))
  pair_count/choose(n,2)-2/3
}
response_test <- function(log_y,epsilon,eta=0.05,noise=NULL,statistic=NULL) {
  n <- length(log_y)
  stopifnot(n>=2,epsilon>0,eta>0,eta<1)
  if(is.null(noise))noise <- rlaplace(1)
  noise_scale <- 2/(n*epsilon)
  threshold <- sqrt(2*log(2/eta)/n)+noise_scale*log(1/eta)
  if(is.null(statistic))statistic <- response_statistic(log_y)
  as.integer(statistic+noise_scale*noise>=threshold)
}

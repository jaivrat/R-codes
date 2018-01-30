library(nloptr)
wbef <- c(0.1, 0.9)

wrec <- c(0.8, 0.2)
threshold <- 10000



#z will be coming as vector of z+(1..N), z-(1..N) 
#f(z+,z-) = 
N <- length(wbef)

eval_f <- function(z) {
  v1 <- z[1:N]
  v2 <- z[(N+1):(2*N)]
  v <- v1 - v2 + wbef - wrec
  J <- sum(v * v)
  #gradients
  Jac <- 2 * v * c(rep(1,N), rep(-1,N))
  return( list( "objective" = J,
                "gradient" = Jac))
}




# constraint functions
# inequalities
eval_g_ineq <- function( z ) {
  v1 <- z[1:N]
  v2 <- z[(N+1):(2*N)]
  
  res_thresh <- sum(v1 + v2) - threshold
  res_wts_gt_0 <- matrix(-1 * v1 + 1 * v2 - wbef, ncol=1)
  res_wts_lt_0 <- matrix(v1 - v2 + wbef - 1, ncol=1)
    
  constr <- rbind(res_thresh,
                  res_wts_gt_0,
                  res_wts_lt_0)
  
  g_res_thresh <- c(rep(1,N), rep(1,N))
  g_res_wts_gt_0 <- matrix(c(-1,1), nrow=1) %x% diag(N)
  g_res_wts_lt_0 <- matrix(c(1,-1), nrow=1) %x% diag(N)
  #gradients
  grad   <- rbind( g_res_thresh,
                   g_res_wts_gt_0,
                   g_res_wts_lt_0)
  return( list( "constraints"=constr, "jacobian"=grad ) )
}


# equalities
eval_g_eq <- function( z ) {
  v1 <- z[1:N]
  v2 <- z[(N+1):(2*N)]
  
  res_sum_wts_is_1 <- sum(v1 -  v2 + wbef) -1 
  
  g_res_sum_wts_is_1 <- matrix(c(rep(1,N), rep(-1,N)), nrow=1)
  
  return( list( "constraints"=res_sum_wts_is_1, "jacobian"=g_res_sum_wts_is_1 ) )
}

lb <- rep(0, 2*N)
ub <- rep(Inf, 2*N)


# initial values
z0 <- runif(2*N)
#z0 <- rep(100, 2*N)

local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel"  = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel"  = 1.0e-7,
              "maxeval"   = 1000,
              "local_opts" = local_opts )


res <- nloptr( x0=z0,
               eval_f=eval_f,
               lb=lb,
               ub=ub,
               eval_g_ineq=eval_g_ineq,
               eval_g_eq=eval_g_eq,
               opts=opts)

res
#Before wts : wbef
wbef
#Recommended: wrec
wrec
#Solution: res$
res$solution
zPlus <- res$solution[1:N]
zMinus <- res$solution[(N+1):(2*N)]
final <- wbef + zPlus - zMinus
final
#Turnover
sum(abs(final-wbef))



# Bat speed distribution functions

my_dskt <- function(x, df, gamma, mu, s) {
  return(skewt::dskt((x - mu) / s, df, gamma) / s)
}

my_pskt <- function(x, df, gamma, mu, s) {
  return(skewt::pskt((x - mu) / s, df, gamma))
}

my_qskt <- function(p, df, gamma, mu, s) {
  return(mu + s * skewt::qskt(p, df, gamma))
}

my_dmix <- function(x, pi, df1, gamma1, mu1, s1, df2, gamma2, mu2, sadd) {
  lik_skt1 <- my_dskt(x, df1, gamma1, mu1, s1)
  lik_skt2 <- my_dskt(x, df2, gamma2, mu2, s1 + sadd)
  return(pi * lik_skt1 + (1 - pi) * lik_skt2)
}

est_bs_dist <- function(bs) {
  nllik <- function(par) {
    dmix_args <- c(list(bs), as.list(par))
    return(-sum(log(do.call(my_dmix, dmix_args))))
  }
  # Parameters: pi; df, gamma, mu, s (x2)
  initpar <- c(0.9, 5, 1, mean(bs), sd(bs), 5, 1, 30, 5)
  opt_res <- optim(initpar, nllik, method = "L-BFGS-B",
                   lower = c(0.01,   1, 0.1, 60,  1,   1, 0.1,  0, 0),
                   upper = c(0.99, 500,   1, 90, 20, 500,  10, 60, 30))
  opt_res$par[9] <- opt_res$par[5] + opt_res$par[9]
  lik_mat <- tibble(
    bat_speed = bs,
    lik_full = opt_res$par[1] * do.call(my_dskt, c(list(bs), as.list(opt_res$par[2:5]))),
    lik_weak = (1 - opt_res$par[1]) * do.call(my_dskt, c(list(bs), as.list(opt_res$par[c(6:9)])))
  ) |>
    mutate(LR = lik_full / lik_weak) |>
    mutate(p_full = LR / (1 + LR))
  return(list(
    pi = opt_res$par[1],
    df1 = opt_res$par[2],
    gamma1 = opt_res$par[3],
    mu1 = opt_res$par[4],
    s1 = opt_res$par[5],
    df2 = opt_res$par[6],
    gamma2 = opt_res$par[7],
    mu2 = opt_res$par[8],
    sadd = opt_res$par[9],
    lik_mat = lik_mat
  ))
}

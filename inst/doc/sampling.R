## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4)

## -----------------------------------------------------------------------------
library(coxdev)
library(survival)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 500
p <- 3
X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, c("x1", "x2", "x3")))
beta <- c(0.5, -0.5, 0.25)
time <- rexp(n, exp(drop(X %*% beta)))
cens <- rexp(n, 0.5)
stop <- ceiling(pmin(time, cens) * 4) / 4
status <- as.integer(time <= cens)
start <- floor(runif(n, 0, 0.5) * stop * 8) / 8
y <- Surv(start, stop, status)

cox <- CoxDeviance$new(y, tie_breaking = "efron")
cox

## -----------------------------------------------------------------------------
b <- coef(coxph(y ~ X, ties = "efron"))
b

## ----draws--------------------------------------------------------------------
eta_hat <- drop(X %*% b)
methods <- c("multinomial", "fisher", "bernoulli", "wild")
draws <- lapply(methods, function(m) {
  crossprod(X, cox$sample_score(eta_hat, method = m, n_draws = 4000, seed = 1))
})
names(draws) <- methods

## -----------------------------------------------------------------------------
breslow <- CoxDeviance$new(y, tie_breaking = "breslow")
info_breslow <- crossprod(X, breslow$information(eta_hat) %*% X)
summary_table <- t(sapply(draws, function(d) {
  c(rowMeans(d) / sqrt(diag(info_breslow)), diag(cov(t(d))) / diag(info_breslow))
}))
colnames(summary_table) <- c(paste("mean", colnames(X)), paste("variance", colnames(X)))
knitr::kable(summary_table, digits = 2, caption = paste(
  "Means of the draws of the score, in units of sqrt(Breslow information),",
  "and variances relative to the Breslow information"))

## -----------------------------------------------------------------------------
U <- cox$score_residuals(eta_hat, X)
sandwich <- crossprod(U)
fit_robust <- coxph(y ~ X, ties = "efron", id = seq_len(n), robust = TRUE)
I_hat <- solve(fit_robust$naive.var)
rbind(wild = diag(cov(t(draws$wild))), sandwich = diag(sandwich),
      coxph = diag(I_hat %*% fit_robust$var %*% I_hat))

## ----fig.alt = "Densities of 4000 draws of the score of x1, by method, with the normal density with the Breslow information as variance."----
cols <- c(multinomial = "black", fisher = "firebrick", bernoulli = "darkorange", wild = "steelblue")
dens <- lapply(draws, function(d) density(d[1, ]))
plot(NULL, xlim = c(-1, 1) * 4 * sqrt(info_breslow[1, 1]),
     ylim = c(0, 1.35 * max(sapply(dens, function(d) max(d$y)))),
     xlab = "score of x1", ylab = "density")
for (m in methods) lines(dens[[m]], col = cols[m], lwd = 2)
curve(dnorm(x, sd = sqrt(info_breslow[1, 1])), add = TRUE, lty = 2)
legend("topleft", c(methods, "N(0, Breslow info.)"), cex = 0.8,
       col = c(cols, "black"), lty = c(rep(1, 4), 2), lwd = c(rep(2, 4), 1), bty = "n")

## -----------------------------------------------------------------------------
all.equal(cox$wild_score(eta_hat, rep(1, n)), -cox$compute(eta_hat)$gradient / 2)
rademacher <- replicate(4000, {
  drop(crossprod(X, cox$wild_score(eta_hat, sample(c(-1, 1), n, replace = TRUE))))
})
rbind(rademacher = diag(cov(t(rademacher))), sandwich = diag(sandwich))

## -----------------------------------------------------------------------------
identical(cox$sample_score(eta_hat, method = "fisher", seed = 7),
          cox$sample_score(eta_hat, method = "fisher", seed = 7))


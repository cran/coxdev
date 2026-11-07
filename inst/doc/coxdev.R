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
b <- rep(0, p)
for (iter in 1:10) {
  eta <- drop(X %*% b)
  res <- cox$compute(eta)
  score <- -drop(crossprod(X, res$gradient)) / 2
  info <- crossprod(X, cox$information(eta) %*% X)
  step <- solve(info, score)
  b <- b + step
  if (max(abs(step)) < 1e-10) break
}
iter
b

## -----------------------------------------------------------------------------
fit <- coxph(y ~ X, ties = "efron")
rbind(coxdev = b, coxph = coef(fit))
rbind(coxdev = sqrt(diag(solve(info))), coxph = sqrt(diag(vcov(fit))))

## -----------------------------------------------------------------------------
eta_hat <- drop(X %*% b)
U <- cox$score_residuals(eta_hat, X)
info_hat <- crossprod(X, cox$information(eta_hat) %*% X)
robust <- solve(info_hat, t(solve(info_hat, crossprod(U))))
fit_robust <- coxph(y ~ X, ties = "efron", id = seq_len(n), robust = TRUE)
rbind(coxdev = sqrt(diag(robust)), coxph = sqrt(diag(vcov(fit_robust))))

## -----------------------------------------------------------------------------
H <- cox$baseline_hazard(eta_hat)
head(H)

## ----include = FALSE----------------------------------------------------------
# survival::basehaz() fails on some R-devel builds (October 2026, in
# stats::delete.response(), as in survival's own CRAN check): the comparison
# below is then left out. TODO: remove once survival's R-devel checks pass.
basehaz_works <- !inherits(try(basehaz(fit, centered = FALSE), silent = TRUE), "try-error")

## ----eval = basehaz_works-----------------------------------------------------
bh <- basehaz(fit, centered = FALSE)
all.equal(H$cumulative_hazard, bh$hazard[match(H$time, bh$time)])

## ----fig.alt = "Estimated survival curves for x1 = -1, 0 and 1, with x2 = x3 = 0."----
plot(NULL, xlim = c(0, max(H$time)), ylim = c(0, 1), xlab = "time", ylab = "survival")
for (x1 in c(-1, 0, 1)) {
  surv <- exp(-H$cumulative_hazard * exp(b[1] * x1))
  lines(c(0, H$time), c(1, surv), type = "s", lty = x1 + 2)
}
legend("topright", paste("x1 =", c(-1, 0, 1)), lty = 1:3, bty = "n")

## ----timing-------------------------------------------------------------------
set.seed(2)
N <- 1e5
XN <- matrix(rnorm(N * p), N, p)
etaN <- drop(XN %*% beta)
stopN <- ceiling(rexp(N, exp(etaN)) * 20) / 20
statusN <- rbinom(N, 1, 0.7)
startN <- floor(runif(N, 0, 0.5) * stopN * 40) / 40
yN <- Surv(startN, stopN, statusN)
etas <- lapply(1:50, function(i) etaN + rnorm(N, sd = 0.1))

per_call <- function(expr) {
  t <- system.time(expr)[["elapsed"]]
  sprintf("%.1f ms", 1000 * t / 50)
}
setup <- system.time(coxN <- CoxDeviance$new(yN))[["elapsed"]]
timings <- c(
  "one object, reused" = per_call(for (e in etas) coxN$compute(e)),
  "a new object for each call" = per_call(for (e in etas) CoxDeviance$new(yN)$compute(e)),
  "one object, same eta again (cached)" = per_call(for (e in etas) coxN$compute(etas[[1]])),
  "one information product" = per_call({
    I <- coxN$information(etas[[1]])
    for (e in etas) I %*% XN[, 1]
  })
)
knitr::kable(data.frame(time = timings), caption = sprintf(
  "Time per call, n = %d. Creating the object (the sort) took %.0f ms.",
  N, 1000 * setup))


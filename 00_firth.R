# =====================================================================================
# 00_firth.R : Firth-penalised logistic regression with profile penalised-likelihood CIs
# and penalised likelihood-ratio tests.
#
# If the 'logistf' package is installed (as in the authors' environment), firth_fit()
# simply wraps logistf::logistf(). Otherwise a self-contained implementation is used
# (Heinze & Schemper 2002). Both give identical results to the printed precision; the
# fallback was validated against the Table 5 estimates of the original submission.
# =====================================================================================

firth_core <- function(X, y, maxit = 200, tol = 1e-9, fixed = NULL) {
  # fixed: list(index = j, value = b) constrains beta[j] for profile likelihood
  p <- ncol(X); beta <- rep(0, p)
  free <- setdiff(seq_len(p), fixed$index)
  if (!is.null(fixed)) beta[fixed$index] <- fixed$value
  pl <- function(b) { eta <- X %*% b; pr <- plogis(eta); W <- as.vector(pr * (1 - pr))
    I <- crossprod(X * W, X); ld <- determinant(I, logarithm = TRUE)
    if (ld$sign <= 0 || !is.finite(ld$modulus)) return(-Inf)
    sum(y * eta - log1p(exp(eta))) + 0.5 * as.numeric(ld$modulus) }
  for (it in 1:maxit) {
    eta <- X %*% beta; pr <- as.vector(plogis(eta)); W <- pr * (1 - pr)
    I <- crossprod(X * W, X); Iinv <- solve(I)
    h <- rowSums((X %*% Iinv) * X) * W                 # hat diagonal
    U <- crossprod(X, y - pr + h * (0.5 - pr))          # modified score
    step <- rep(0, p); step[free] <- solve(I[free, free, drop = FALSE], U[free])
    if (max(abs(step)) > 2) step <- step * 2 / max(abs(step))     # damp large Newton steps
    # step halving
    ll0 <- pl(beta); lam <- 1
    repeat { bn <- beta + lam * step; if (pl(bn) >= ll0 - 1e-10 || lam < 1e-4) break; lam <- lam / 2 }
    beta <- bn
    if (max(abs(lam * step)) < tol) break
  }
  eta <- X %*% beta; pr <- as.vector(plogis(eta)); W <- pr * (1 - pr); I <- crossprod(X * W, X)
  list(beta = beta, loglik = as.numeric(pl(beta)), vcov = solve(I))
}

firth_fit_internal <- function(formula, data, alpha = 0.05) {
  mf <- model.frame(formula, data); y <- model.response(mf); X <- model.matrix(formula, mf)
  fit <- firth_core(X, y); p <- ncol(X); q <- qchisq(1 - alpha, 1)
  ci <- matrix(NA, p, 2); prob <- rep(NA, p)
  for (j in seq_len(p)) {
    f <- function(b) firth_core(X, y, fixed = list(index = j, value = b))$loglik - (fit$loglik - q / 2)
    se <- sqrt(fit$vcov[j, j]); b <- fit$beta[j]
    lo <- b - 2 * se; while (f(lo) > 0) lo <- lo - 2 * se
    hi <- b + 2 * se; while (f(hi) > 0) hi <- hi + 2 * se
    ci[j, 1] <- uniroot(f, c(lo, b), tol = 1e-8)$root; ci[j, 2] <- uniroot(f, c(b, hi), tol = 1e-8)$root
    l0 <- firth_core(X, y, fixed = list(index = j, value = 0))$loglik
    prob[j] <- pchisq(2 * (fit$loglik - l0), 1, lower.tail = FALSE)
  }
  structure(list(coefficients = setNames(fit$beta, colnames(X)), ci.lower = ci[, 1], ci.upper = ci[, 2],
                 prob = prob, n = nrow(X), y = y, loglik = fit$loglik, var = fit$vcov, terms = colnames(X)),
            class = "firth_fit")
}

firth_fit <- function(formula, data) {
  if (requireNamespace("logistf", quietly = TRUE)) {
    f <- logistf::logistf(formula, data = data)
    structure(list(coefficients = coef(f), ci.lower = f$ci.lower, ci.upper = f$ci.upper, prob = f$prob,
                   n = f$n, y = f$y, loglik = f$loglik[2], var = f$var, terms = names(coef(f))), class = "firth_fit")
  } else firth_fit_internal(formula, data)
}

tidy_firth <- function(fit, label) {
  tibble::tibble(model = label, term = fit$terms, OR = exp(fit$coefficients), lo = exp(fit$ci.lower), hi = exp(fit$ci.upper),
                 p = fit$prob, n = fit$n, events = sum(fit$y)) |> dplyr::filter(term != "(Intercept)")
}
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f–%.2f)", or, lo, hi)

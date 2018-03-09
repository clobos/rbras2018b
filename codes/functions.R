# =====================================================================
# Fitting flexible models for count data
#                                                        Eduardo Junior
#                                                    edujrrib@gmail.com
#                                                            2018-03-05
# =====================================================================

#-------------------------------------------
# Negative of log-likelihood function for COM-Poisson (RibeiroJr2018)
llcmp <- function(params, beta = NULL,
                   X, Z, y, sumto = 500L) {
    # Obtem os parametros
    beta <- params[-1]
    phi <- params[ 1]
    nu <- exp(phi)
    xb <- X %*% beta
    rlambda <- exp(xb) + (nu - 1) / (2 * nu)
    lrlambda <- log(rlambda)
    # Obtem as constantes
    j <- 1:sumto
    zi <- vapply(lrlambda, FUN = function(x) {
        1 + sum(exp(nu * j * x - nu * lfactorial(j)))
    }, FUN.VALUE = numeric(1))
    # Calcula o negativo do log da função de verossimilhança
    ll <- sum(y * nu * lrlambda) -
        sum(nu * lfactorial(y)) -
        sum(log(zi))
    return(-ll)
}

# Probability mass function for COM-Poisson
dcmp <- function(y, mu, phi, sumto = 500L) {
    vapply(y, function(yi) {
        exp(-llcmp(params = c(phi, log(mu)),
                   X = cbind(1L),
                   y = yi))
    }, numeric(1))
}

# Moments for COM-Poisson
moments_cmp <- function(mu, phi, sumto = 500L, tol = 1e-6) {
    nu <- exp(phi)
    ap_stde <- sqrt(mu / nu)
    ymax <- ceiling(mu + 5 * ap_stde)
    pmax <- dcmp(ymax, mu, phi, sumto = sumto)
    # Verifica se prob(ymax) é pequena o suficiente
    while (pmax > tol) {
        ymax <- ymax + 1L
        pmax <- dcmp(ymax, mu, phi, sumto = sumto)
    }
    yrange <- 1:ymax
    expec_y2 <- sum(yrange^2 * dcmp(yrange, mu, phi, sumto = sumto))
    variance <- expec_y2 - mu^2
    return(list("mean" = mu, "variance" = variance))
}

#-------------------------------------------
# Negative of log-likelihood function for Gamma-Count (Zeviani2014)
llgct <- function (params, X, y) {
    # Obtem os parametros
    beta <- params[-1]
    gama <- params[ 1]
    alpha <- exp(gama)
    xb <- X %*% beta
    eta <- exp(xb)
    alpha_eta <- alpha * eta
    alpha_y <- alpha * y
    ll <- sum(log(
        pgamma(1L, shape = alpha_y, rate = alpha_eta) -
        pgamma(1L, shape = alpha_y + alpha, rate = alpha_eta)
    ))
    return(-ll)
}

# Probability mass function for Gamma-Count
dgct <- function(y, kapa, gama) {
    vapply(y, function(yi) {
        exp(-llgct(params = c(gama, log(kapa)),
                   X = cbind(1L),
                   y = yi))
    }, numeric(1))
}

# Moments for Gamma-Count
moments_gct <- function(kapa, gama, tol = 1e-6) {
    alpha <- exp(gama)
    ap_mean <- kapa / alpha
    ap_stde <- sqrt(kapa) / alpha
    ymax <- ceiling(ap_mean + 5 * ap_stde)
    pmax <- dgct(ymax, kapa, gama)
    # Verifica se prob(ymax) é pequena o suficiente
    while (pmax > tol) {
        ymax <- ymax + 1L
        pmax <- dgct(ymax, kapa, gama)
    }
    yrange <- 1:ymax
    expec_y2 <- sum(yrange^2 * dgct(yrange, kapa, gama))
    expectat <- sum(yrange   * dgct(yrange, kapa, gama))
    variance <- expec_y2 - expectat^2
    return(list("mean" = expectat, "variance" = variance))
}

#-------------------------------------------
# Negative of log-likelihood function for Gamma-Count (Zamani2012)
llgpo <- function (params, X, y) {
    beta <- params[-1]
    sigma <- params[ 1]
    xb <- X %*% beta
    mu <- exp(X %*% params[-1])
    sigma_mu <- 1 + sigma * mu
    sigma_y <- 1 + sigma * y
    ll <- sum(y * (xb - log(sigma_mu)) +
              (y - 1) * log(sigma_y) -
              mu * (sigma_y / sigma_mu) -
              lfactorial(y))
    return(-ll)
}

# Probability mass function for Gamma-Count
dgpo <- function(y, mu, sigma) {
    vapply(y, function(yi) {
        exp(-llgpo(params = c(sigma, log(mu)),
                   X = cbind(1L),
                   y = yi))
    }, numeric(1))
}

# Moments for Poisson generalizada
moments_gpo <- function(mu, sigma) {
    psi <- exp(sigma)
    variance <- mu * (1 - mu * psi)^2
    return(list("mean" = mu, "variance" = variance))
}

#-------------------------------------------
# Functions for Poisson-Tweedie distribution
moments_ptw <- function(mu, omega, p) {
    variance <- mu * (1 + omega * mu^(p-1))
    return(list("mean" = mu, "variance" = variance))
}

#-------------------------------------------
# Find the parameters such lead to fixed DI and Expectation
system_equation <- function(param,
                            di,
                            expect,
                            moments_function,
                            ...,
                            trace = FALSE) {
    par1 <- param[1]
    par2 <- param[2]
    if (trace) print(param)
    moments <- moments_function(par1, par2, ...)
    eq1 <- moments[["mean"]] -  expect
    eq2 <- moments[["variance"]] / moments[["mean"]] - di
    return(c(eq1, eq2))
}

get_parameters <- function(expect, di, ...,
                           model = c("CMP", "GCT", "GPo", "PTw")) {
    if (model %in% c("CMP", "GCT")) {
        fun <- if (model == "CMP") moments_cmp else moments_gct
        safely_fun <- function(x) list(root = c(NA, NA))
        out <- tryCatch(
            rootSolve::multiroot(
                system_equation,
                di = di,
                expect = expect,
                start = c(expect, -log(di)),
                moments_function = fun),
            error   = safely_fun,
            warning = safely_fun
        )
        out <- out$root
    }
    if (model == "GPo") {
        if (sqrt(di) < 1L) return(c(NA, NA))
        out <- c(expect, log((sqrt(di) - 1)) - log(expect))
    }
    if (model == "PTw") {
        p <- list(...)$p
        if (di < 1L) return(c(NA, NA))
        out <- c(expect, (di - 1)/expect^(p - 1))
    }
    return(out)
}

#-----------------------------------------------------------------------
# Framework for fitting flexible count models
fit_cm <- function(formula,
                   data,
                   model = c("CMP", "GCT", "GPo"),
                   start = NULL,
                   sumto = 500L) {
    model <- match.arg(model)
    settings <- switch(
        model,
        "CMP" = list(pname = "phi",   ll = llcmp),
        "GCT" = list(pname = "gamma", ll = llgct),
        "GPo" = list(pname = "sigma", ll = llgpo)
    )
    #-------------------------------------------
    # Obtendo as matrizes
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    X <- model.matrix(terms, frame)
    y <- model.response(frame)
    #-------------------------------------------
    # Parametros iniciais
    if (is.null(start)) {
        m0 <- glm.fit(x = X, y = y, family = poisson())
        start <- c(0, m0$coefficients)
        names(start)[1L] <- settings[["pname"]]
    }
    #-------------------------------------------
    # Ajuste do modelo
    compute_loglik <- settings[["ll"]]
    bbmle::parnames(compute_loglik) <- names(start)
    fixed <- list(X = X, y = y, formula = formula)
    if (model == "CMP") fixed[["sumto"]] <- sumto
    fit <- suppressWarnings(
        bbmle::mle2(compute_loglik,
                    start = start,
                    data = fixed,
                    method = "BFGS",
                    vecpar = TRUE))
    return(fit)
}

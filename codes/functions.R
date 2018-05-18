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
dcmp <- function(y, mu, phi, sumto = 500L, log = FALSE) {
    nu <- exp(phi)
    lrlambda <- log(mu + (nu - 1)/(2*nu))
    j <- 1:sumto
    Z <- 1 + sum(exp(j * nu * lrlambda - nu * lfactorial(j)))
    out <- y * nu * lrlambda - nu * lfactorial(y) - log(Z)
    if (!log) out <- exp(out)
    return(out)
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
dgct <- function(y, kapa, gama, log = FALSE) {
    alpha <- exp(gama)
    out <- pgamma(q = 1L,
                  shape = alpha * y,
                  rate = alpha * kapa) -
        pgamma(q = 1L,
               shape = alpha * (y + 1),
               rate = alpha * kapa)
    if (log) out <- log(out)
    return(out)
}

# Moments for Gamma-Count
moments_gct <- function(kapa, gama, tol = 1e-6) {
    alpha <- exp(gama)
    ap_stde <- sqrt(kapa / alpha)
    ymax <- ceiling(kapa + 5 * ap_stde)
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
dgpo <- function(y, mu, sigma, log = FALSE) {
    sigma_mu <- 1 + sigma * mu
    sigma_y <- 1 + sigma * y
    out <- suppressWarnings(
        y * (log(mu) - log(sigma_mu)) +
        (y - 1) * log(sigma_y) -
        mu * (sigma_y / sigma_mu) -
        lfactorial(y))
    out[is.nan(out)] <- -Inf
    if (!log) out <- exp(out)
    return(out)
}

# Moments for Poisson generalizada
moments_gpo <- function(mu, sigma) {
    variance <- mu * (1 + mu * sigma)^2
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

get_parameters <- function(expect, di,
                           model = c("CMP", "GCT", "GPo", "PTw"),
                           ...) {
    model <- match.arg(model)
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
        # if (sqrt(di) < 1L) return(c(NA, NA))
        out <- c(expect, (sqrt(di) - 1)/expect)
    }
    if (model == "PTw") {
        p <- list(...)$p
        if (di < 1L) return(c(NA, NA))
        out <- c(expect, (di - 1)/expect^(p - 1))
    }
    return(out)
}

#-----------------------------------------------------------------------
# Check if pmf is a genuine probability function
check_pmf <- function(params, y = 0:500, tol = 1e-2,
                      model = c("CMP", "GCT", "GPo"),
                      ...) {
    model <- match.arg(model)
    pmf <- switch(model,
                  "CMP" = dcmp,
                  "GCT" = dgct,
                  "GPo" = dgpo)
    spr <- sum(pmf(y, params[1], params[2], ...))
    out <- ifelse(abs(1 - spr) > tol, 0L, 1L)
    return(out)
}

#-----------------------------------------------------------------------
# Auxiliary function to build mean and variance data
make_mvdata <- function(settings, model) {
    # Parameters
    lista <- settings[[model]]
    p1 <- lista[[1]]
    p2 <- lista[[2]]
    aux <- expand.grid(
        par1 = seq(p1[1], p1[2], length.out = 70),
        par2 = seq(p2[1], p2[2], length.out = 20),
        KEEP.OUT.ATTRS = FALSE)
    wrap_moments <- function(x, ...) {
        names(x) <- names(lista)[1:2]
        args <- c(x, lista[-c(1:2, length(lista))])
        out <- do.call(lista$fun, args)
        unlist(out)
    }
    moments <- t(apply(aux, 1L, wrap_moments))
    colnames(moments) <- c("mean", "variance")
    out <- cbind(aux, moments)
    out <- transform(out, di = variance / mean)
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
    fixed <- list(y = y, X = X, formula = formula, model = model)
    if (model == "CMP") fixed[["sumto"]] <- sumto
    fit <- suppressWarnings(
        bbmle::mle2(compute_loglik,
                    start = start,
                    data = fixed,
                    method = "BFGS",
                    vecpar = TRUE))
    return(fit)
}

#-----------------------------------------------------------------------
# Fitted values and confidence intervals
eta2mean_gct <- function(eta, gama, tol = 1e-5) {
    vapply(eta, function(etai) {
        alpha <- exp(gama)
        kapa <- exp(etai)
        ap_stde <- sqrt(kapa / alpha)
        ymax <- ceiling(kapa + 5 * ap_stde)
        pmax <- dgct(ymax, kapa, gama)
        # Verifica se prob(ymax) é pequena o suficiente
        while (pmax > tol) {
            ymax <- ymax + 1L
            pmax <- dgct(ymax, kapa, gama)
        }
        yrange <- 1:ymax
        expectat <- sum(yrange   * dgct(yrange, kapa, gama))
        return(expectat)
    }, FUN.VALUE = double(1))
}

predict_ptw <- function(object,
                        newdata,
                        type = c("response", "link"),
                        interval = c("confidence", "none"),
                        level = 0.95,
                        augment_data = TRUE) {
    type <- match.arg(type)
    interval <- match.arg(interval)
    #-------------------------------------------
    Vcov <- vcov(object)
    indb <- grep("beta", rownames(Vcov))
    Vmar <- as.matrix(Vcov[indb, indb])
    #-------------------------------------------
    formula <- object$linear_pred[[1]]
    formula[[2]] <- NULL
    beta <- coef(object)[indb, "Estimates"]
    #-------------------------------------------
    X <- model.matrix(formula, newdata)
    est <- X %*% beta
    if (interval == "none") {
        out <- data.frame("fit" = est)
    }
    if (interval == "confidence") {
        qn <- -qnorm((1 - level[1])/2)
        std <- sqrt(diag(tcrossprod(X %*% Vmar, X)))
        out <- data.frame("fit" = est,
                          "lwr" = est - qn * std,
                          "upr" = est + qn * std)
    }
    if (type == "response") {
        out <- data.frame(apply(out, 2L, exp))
    }
    if (augment_data) out <- cbind(newdata, out)
    return(out)
}

predict_cm <- function(object,
                       newdata,
                       type = c("response", "link"),
                       interval = c("confidence", "none"),
                       level = 0.95,
                       augment_data = TRUE) {
    #-------------------------------------------
    if (class(object) == "mcglm") {
        out <- predict_ptw(object, newdata, type,
                           interval, level, augment_data)
        return(out)
    }
    #-------------------------------------------
    type <- match.arg(type)
    interval <- match.arg(interval)
    model <- object@data$model
    #-------------------------------------------
    Vcov <- vcov(object)
    Vbeta <- Vcov[-1, -1, drop = FALSE]
    Vdisp <- Vcov[ 1,  1, drop = FALSE]
    Vbedi <- Vcov[ 1, -1]
    #-------------------------------------------
    Vcond <- Vbeta - tcrossprod(Vbedi %*% solve(Vdisp), Vbedi)
    formula <- object@data$formula
    formula[[2]] <- NULL
    beta <- object@coef[-1]
    #-------------------------------------------
    X <- model.matrix(formula, newdata)
    est <- X %*% beta
    if (interval == "none") {
        out <- data.frame("fit" = est)
    }
    if (interval == "confidence") {
        qn <- -qnorm((1 - level[1])/2)
        std <- sqrt(diag(tcrossprod(X %*% Vcond, X)))
        out <- data.frame("fit" = est,
                          "lwr" = est - qn * std,
                          "upr" = est + qn * std)
    }
    if (type == "response") {
        if (model == "GCT")
            out <- data.frame(apply(out, 2L, eta2mean_gct,
                                    gama = object@coef[1]))
        else
            out <- data.frame(apply(out, 2L, exp))
    }
    if (augment_data) out <- cbind(newdata, out)
    return(out)
}


#-----------------------------------------------------------------------
# Compute log-Liklihood for Poisson-Tweedie models (functions by
# Bonat2018)
library(tweedie)
library(statmod)

integrand <- function(x, y, mu, phi, power) {
    int = dpois(y, lambda = x)*dtweedie(x, mu = mu,
                                        phi = phi, power = power)
    return(int)
}

# Numerical integration using Gauss-Laguerre method
gauss_laguerre <- function(integrand, n_pts, y, mu, phi, power) {
    pts <- gauss.quad(n_pts, kind="laguerre")
    integral <- sum(pts$weights*integrand(pts$nodes, y = y, mu = mu,
                                          phi = phi, power = power)/
                        exp(-pts$nodes))
    return(integral)
}
gauss_laguerre_vec <- Vectorize(FUN = gauss_laguerre,
                                vectorize.args = "y")

# Numerical integration using Monte Carlo method
monte_carlo <- function(integrand, n_pts, y, mu, phi, power) {
    pts <- rtweedie(n_pts, mu = mu, phi = phi, power = power)
    norma <- dtweedie(pts, mu = mu, phi = phi, power = power)
    integral <- mean(integrand(pts, y = y, mu = mu, phi = phi,
                               power = power)/norma)
    return(integral)
}

# Probability mass function Poisson-Tweedie
dptweedie_aux <- function(y, mu, phi, power, n_pts, method) {
    if(method == "laguerre" | y > 0) {
        pmf <- gauss_laguerre(integrand = integrand, n_pts = n_pts,
                              y = y, mu = mu, phi = phi, power = power)
    }
    if(method == "laguerre" & y == 0) {
        v.y <- round(10*sqrt(mu + phi*mu^power),0)
        if(v.y > 1000) {v.y <- 1000}
        #print(v.y)
        y.grid <- 0:v.y
        temp <- gauss_laguerre_vec(integrand = integrand, n_pts = n_pts,
                                   y = y.grid, mu = mu, phi = phi,
                                   power = power)
        pmf <- 1-sum(temp)+temp[1]
    }
    if(method == "MC") {
        pmf <- monte_carlo(integrand = integrand, n_pts = n_pts,
                           y = y, mu = mu, phi = phi, power = power)
    }
    return(pmf)
}

# Vectorize version
dptweedie <- Vectorize(FUN = dptweedie_aux, vectorize.args = c("y","mu"))

# LogLik function
logLik.mcglm <- function(object, n_pts = 50L, method = "laguerre"){
    y <- as.vector(object$observed)
    mu <- object$fitted
    phi <- coef(object, type = "tau")[["Estimates"]]
    power <- coef(object, type = "power")[["Estimates"]]
    # Conditions
    if (object$variance != "poisson_tweedie" | phi < 0) {
        stop("Method is not valid")
    }
    ll <- sum(log(dptweedie(y = y,
                            mu = mu,
                            phi = phi,
                            power = power,
                            n_pts = n_pts,
                            method = method)))
    attr(ll, "df") <- nrow(coef(object))
    attr(ll, "nobs") <- length(mu)
    class(ll) <- "logLik"
    return(ll)
}

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

# Negative of log-likelihood function for Gamma-Count (Zamani2012)
llgpo <- function (params, X, y) {
    beta <- params[-1]
    sigma <- params[ 1]
    psi <- exp(sigma)
    xb <- X %*% beta
    mu <- exp(X %*% params[-1])
    psi_mu <- 1 + psi * mu
    psi_y <- 1 + psi * y
    ll <- sum(y * (xb - log(psi_mu)) +
              (y - 1) * log(psi_y) -
              mu * (psi_y / psi_mu) -
              lfactorial(y))
    return(-ll)
}

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

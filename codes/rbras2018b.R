#-------------------------------------------
# Setting session

library(xtable)
options(xtable.caption.placement = "top",
        xtable.booktabs = TRUE,
        xtable.sanitize.text.function = identity)

library(dplyr)
library(purrr)
library(lattice)
library(latticeExtra)
library(gridExtra)

source("functions.R")
source("lattice-panels.R")

#-----------------------------------------------------------------------
# Study distributions flexibility

## ---- compute_indexes

# Find parameters to fixed E(Y) = 40 and DI(Y) = 4 and 0.25
par_cmp <- vapply(c(4, 0.25), get_parameters, expect = 50,
                  model = "CMP", FUN.VALUE = double(2))
par_gct <- vapply(c(4, 0.25), get_parameters, expect = 50,
                  model = "GCT", FUN.VALUE = double(2))
par_gpo <- vapply(c(4, 0.25), get_parameters, expect = 50,
                  model = "GPo", FUN.VALUE = double(2))
par_ptw1 <- vapply(c(4, 0.25), get_parameters, expect = 50, p = 1.1,
                   model = "PTw", FUN.VALUE = double(2))
par_ptw2 <- vapply(c(4, 0.25), get_parameters, expect = 50, p = 2,
                   model = "PTw", FUN.VALUE = double(2))
par_ptw3 <- vapply(c(4, 0.25), get_parameters, expect = 50, p = 3,
                   model = "PTw", FUN.VALUE = double(2))

# Stack the parameters in a setting list
settings <- list(
    "CMP" = list("mu" = c(2, 50), "phi" = par_cmp[2, ],
                 "fun" = moments_cmp),
    "GCT" = list("kapa" = c(2, 50), "gama" = par_gct[2, ],
                 "fun" = moments_gct),
    "GPo" = list("mu" = c(2, 50), "sigma" = par_gpo[2, ],
                 "fun" = moments_gpo),
    "PTw1" = list("mu" = c(2, 50), "omega" = c(0, par_ptw1[2, 1]),
                  "p" = 1.1, "fun" = moments_ptw),
    "PTw2" = list("mu" = c(2, 50), "omega" = c(0, par_ptw2[2, 1]),
                  "p" = 2, "fun" = moments_ptw),
    "PTw3" = list("mu" = c(2, 50), "omega" = c(0, par_ptw3[2, 1]),
                  "p" = 3, "fun" = moments_ptw)
)

# Compute mean, variance and dispersion indexes
models <- names(settings); names(models) <- models
damv <- map_dfr(models, make_mvdata, settings = settings,
                .id = "model")

# Colors for the graphics
pallete <- brewer.pal(9, "RdBu")
pallete[4:6] <- c("gray80", "gray90", "gray80")
myreg <- colorRampPalette(pallete)
fvltex <- parse(text = c("'CMP'*(mu*','~phi)",
                         "'GCT'*(kappa*','~gamma)",
                         "'GPo'*(mu*','~sigma)",
                         "'PTw'[1*','*1]*(mu*','~omega)",
                         "'PTw'[2]*(mu*','~omega)",
                         "'PTw'[3]*(mu*','~omega)"))
names(fvltex) <- models

#-------------------------------------------
# Mean and variance

## ---- plot-mean-variance
xys <- map(models, function(m) {
    da <- filter(damv, model == m)
    colp <- myreg(length(unique(da$par2)))
    args <- list(key = list(
                     space = "top",
                     col = colp,
                     at = unique(da$par2),
                     lab = list(rot = 0, cex = 0.9),
                     draw = FALSE))
    xy <- xyplot(variance ~ mean | model,
                 groups = par2,
                 type = "l",
                 xlab = expression("E"*(Y)),
                 ylab = expression("Var"*(Y)),
                 scales = list(
                     y = list(rot = 90)),
                 legend = list(
                     top = list(
                         fun = draw.colorkey,
                         args = args)),
                 strip = strip.custom(
                     factor.levels = fvltex[m]),
                 par.settings = list(
                     superpose.line = list(
                         col = colp),
                     layout.widths = list(
                         right.padding = 0)),
                 data = da)
})

layout <- matrix(seq_len(6), nrow = 2, ncol = 3, byrow = TRUE)
marrangeGrob(xys, top = "", layout_matrix = layout)

#-------------------------------------------
# Dispersion index

## ---- plot-dispersion-index
xys <- map(models, function(m) {
    da <- filter(damv, model == m)
    colp <- myreg(length(unique(da$par2)))
    args <- list(key = list(
                     space = "top",
                     col = colp,
                     at = unique(da$par2),
                     draw = FALSE))
    xy <- xyplot(di ~ mean | model,
                 groups = par2,
                 type = "l",
                 xlab = expression("E"*(Y)),
                 ylab = "DI",
                 ylim = extendrange(c(0, 4.5)),
                 scales = list(
                     y = list(rot = 90)),
                 legend = list(
                     top = list(
                   fun = draw.colorkey,
                   args = args)),
                 strip = strip.custom(
                     factor.levels = fvltex[m]),
                 par.settings = list(
                     superpose.line = list(
                         col = colp),
                     layout.widths = list(
                         right.padding = 0)),
                 panel = function(x, y, ...) {
                     panel.xyplot(x, y, ...)
                     panel.curve(1 + 0*x, min(x), max(x), lty = 2)
                 },
                 data = da)
})

layout <- matrix(seq_len(6), nrow = 2, ncol = 3, byrow = TRUE)
marrangeGrob(xys, top = "", layout_matrix = layout)

#-----------------------------------------------------------------------
# Case study

# ---- fit-models
sitophilus <- data.frame(
    "extract" = rep(c("Folha", "Ramo", "Semente", "Controle"),
                    each = 10),
    "ninsect" = c(19, 20, 36, 32, 18, 47, 38, 31, 32, 40, 20, 34, 41,
                  29, 31, 15, 31, 33, 45, 20, 4, 0, 1, 1, 1, 0, 2, 0, 2,
                  0, 35, 26, 41, 34, 23, 29, 39, 34, 16, 38)
)

library(bbmle)
library(mcglm)
m1cmp <- fit_cm(ninsect ~ extract, data = sitophilus,
                model = "CMP", sumto = 50)
m1gct <- fit_cm(ninsect ~ extract, data = sitophilus,
                model = "GCT")
m1gpo <- fit_cm(ninsect ~ extract, data = sitophilus,
                model = "GPo")
m1ptw <- mcglm(linear_pred = c(ninsect ~ extract),
               matrix_pred = list(mc_id(sitophilus)),
               variance = "poisson_tweedie",
               power_fixed = FALSE,
               link = "log",
               data = sitophilus,
               control_algorithm =
                   list(max_iter = 100, correct = FALSE))

#-----------------------------------------------------------------------
# Fitted values

## ---- fitted-values

# Fitted mean of the number of insects
mnames <- c("CMP", "GCT", "GPo", "PTw")
pred <- unique(sitophilus[, "extract", drop = FALSE])
models <- list(m1cmp, m1gct, m1gpo, m1ptw)
names(models) <- mnames

pred <- map_dfr(models,
                predict_cm,
                newdata = pred,
                .id = "model") %>%
    mutate(model = as.factor(model)) %>%
    arrange(extract, model)

# Fitted values and fitted mean-variance relationship
museq <- seq(min(pred$fit), max(pred$fit), length = 50)

variance_cmp <- function(...) moments_cmp(...)$variance
vacmp <- map_dbl(museq, variance_cmp, phi = coef(m1cmp)[1])
vaseq <- c(vacmp)

variance_gct <- function(...) moments_gct(...)$variance
vagct <- map_dbl(museq, variance_gct, gama = coef(m1gct)[1])
vaseq <- c(vaseq, vagct)

variance_gpo <- function(...) moments_gpo(...)$variance
vagpo <- map_dbl(museq, variance_gpo, sigma = coef(m1gpo)[1])
vaseq <- c(vaseq, vagpo)

variance_ptw <- function(...) moments_ptw(...)$variance
vaptw <- map_dbl(museq, variance_ptw, omega = coef(m1ptw)[6, 1],
                 p = coef(m1ptw)[5, 1])
vaseq <- c(vaseq, vaptw)

mvfit <- data.frame(model = rep(mnames, each = length(museq)),
                    mean = museq, variance = vaseq)
mvobs <- sitophilus %>%
    group_by(extract) %>%
    summarise(mu = mean(ninsect), va = var(ninsect))

# Make graphs

## ---- fitted-plot
cols <- trellis.par.get("superpose.line")$col[1:4]
pchs <- c(8, 15, 17, 19)
keyt <- parse(text = c("'COM-Poisson'",
                       "italic('Gamma-Count')",
                       "'Poisson generalizada'",
                       "'Poisson-Tweedie'"))
key <- list(columns = 1,
            type = "o",
            divide = 1,
            lines = list(pch = pchs),
            text = list(keyt))

xy1 <- segplot(extract ~ lwr + upr,
               as.table = TRUE,
               centers = fit,
               groups = model,
               horizontal = FALSE,
               draw = FALSE,
               ylab = "Número de insetos emergentes",
               sub = "(a)",
               lwd = 1.5,
               gap = 0.2,
               cex = 1,
               pch = pchs,
               key = key,
               ylim = c(-3, 50),
               panel = panel.groups.segplot,
               data = pred) +
    xyplot(ninsect ~ extract,
           spread = 0.08,
           panel = panel.beeswarm,
           data = sitophilus,
           alpha = 0.3)

xy2 <- xyplot(variance ~ mean,
              groups = model,
              type = "l",
              xlab = "Esperança",
              ylab = "Variância",
              sub = "(b)",
              auto.key = list(
                  columns = 1,
                  points = FALSE,
                  lines = TRUE,
                  text = keyt
              ),
              par.settings = list(
                  superpose.line = list(
                      lty = 1:7,
                      col = 1,
                      lwd = 1.5
                  )
              ),
              data = mvfit) +
    as.layer(xyplot(va ~ mu, data = mvobs, col = "gray50"))

print(xy1, split = c(1, 1, 2, 1), more = TRUE)
print(xy2, split = c(2, 1, 2, 1), more = FALSE)


## ---- noinclude

############
### BART ###
############

library(BART)

BART_autoparams <- function(d.train) {
    autoparams(start = list(theta = 0, # zero means random (learned from the data)
                            omega = 0, # zero means random
                            a = .5,
                            b = 1,
                            rho = 150,
                            k = 2.5,
                            ntree = 200L),
               range = list(theta = NULL,
                            omega = NULL,
                            a = NULL, # c(.5, 1),
                            b = NULL, # c(0, 10),
                            rho = c(5, ncol(d.train) - 3),
                            k = c(0, 5),
                            ntree = NULL),
               max_iter = 2,
               loss = function(aggr) {
                   r <- aggr$pred_perf_pval
                   if (is.na(r) || is.null(r))
                       return(1)
                   else
                       return(r)
               })
}

BART_model <- function(d.train, d.test, params, strata = NULL, return.fit = FALSE, verbose = FALSE) {
    if (!is.null(strata))
        stop("Strata are not supported.")

    # Grow the trees
    bart.fit <- wbart(
        x.train = as.matrix(d.train[, !c("GENE", "IID", "EXPRESSION")]),
        y.train = d.train$EXPRESSION,
        x.test = if (is.null(d.test)) NULL else as.matrix(d.test[, !c("GENE", "IID", "EXPRESSION")]),
        sparse = FALSE,
        augment = FALSE,
        usequants = TRUE,
        cont = TRUE,
        # rm.const = TRUE,
        theta = params$theta,
        omega = params$omega,
        a = params$a,
        b = params$b,
        rho = params$rho,
        k = params$k,
        ntree = params$ntree,
        # nskip = params$nskip,
        # ndpost = params$ndpost,
        # keepevery = params$keepevery,
        printevery = 500
    )

    if (return.fit)
        return(bart.fit)

    ypred <- bart.fit$yhat.test.mean
    pearson <- compute_pearson(d.test$EXPRESSION, ypred, by = d.test$GENE)
    pearson <- cbind(
        theta = params$theta,
        omega = params$omega,
        a = params$a,
        b = params$b,
        rho = params$rho,
        k = params$k,
        ntree = params$ntree,
        data.name = pearson$data.name,
        mse = mean(bart.fit$sigma), # technically not the mse
        rsq = cvrsq(d.test$EXPRESSION, ypred),
        pearson[, c("parameter.df", "estimate.cor", "p.value")]
    )

    pearson
}


###########
### PCR ###
###########

library(pls)

PCR_autoparams <- function(d.train) {
    autoparams(start = list(ncomp = "NULL"), # NULL here means to use the max possible ncomp
               range = list(ncomp = NULL),
               max_iter = 0,
               loss = function(aggr) {
                   r <- aggr$pred_perf_pval
                   if (is.na(r) || is.null(r))
                       return(1)
                   else
                       return(r)
               })
}

PCR_model <- function(d.train, d.test, params.grid = NULL, params.auto = NULL, n.folds = 10,
                      folds.strata = NULL, model.strata = NULL, return.fit = FALSE, threads = 1, verbose = FALSE) {
    # Validate the arguments
    if (is.null(params.grid) && is.null(params.auto))
        stop("At least one of params.grid and params.auto must be specified.")
    if (!is.null(params.grid) && !is.null(params.auto))
        stop("At most one of params.grid and params.auto can be specified.")
    if (!is.null(folds.strata))
        stop("Folds strata are not supported.")
    if (!is.null(model.strata))
        stop("Model strata are not supported.")

    if (params.auto@best$ncomp == "NULL")
        ncomp <- min(nrow(d.train) - 1, ncol(d.train) - 3)
    else
        ncomp <- params.auto@best$ncomp

    try({
    if (!is.null(params.auto)) {
        pcr.fit <- pls::pcr(
            EXPRESSION ~ .,
            data = d.train[, !c("GENE", "IID")],
            ncomp = ncomp,
            segments = n.folds,
            scale = FALSE,
            validation = "CV")
    } else {
        stop("params.grid is not supported for PCR_model.")
    }

    if (return.fit)
        return(pcr.fit)

    cv.pearsons <- lapply(1:pcr.fit$ncomp, function(p) {
        pearson <- compute_pearson(d.train$EXPRESSION, pcr.fit$validation$pred[, 1, p], by = d.train$GENE)
        cbind(
            ncomp = p,
            data.name = pearson$data.name,
            mse = pls::MSEP(pcr.fit, estimate = "CV")$val[1, 1, p + 1],
            rsq = pls::R2(pcr.fit, estimate = "CV")$val[1, 1, p + 1],
            pearson[, c("parameter.df", "estimate.cor", "p.value")]
        )
    })
    cv.ncomps <- which.min(sapply(cv.pearsons, function(p) mean(p$p.value)))
    cv.pearsons <- do.call("rbind", cv.pearsons)

    # NOTE: The MSEP() and R2() functions are VERY slow... Consider replacing them
    ypred <- predict(pcr.fit, newdata = d.test[, !c("GENE", "IID", "EXPRESSION")], comps = cv.ncomps)
    pearson <- compute_pearson(d.test$EXPRESSION, ypred[, 1], by = d.test$GENE)
    pearson <- cbind(
        ncomp = cv.ncomps,
        data.name = pearson$data.name,
        mse = pls::MSEP(pcr.fit, estimate = "test", newdata = d.test[, !c("GENE", "IID")], comps = cv.ncomps)$val[1, 1, 1],
        rsq = pls::R2(pcr.fit, estimate = "test", newdata = d.test[, !c("GENE", "IID")], comps = cv.ncomps)$val[1, 1, 1],
        pearson[, c("parameter.df", "estimate.cor", "p.value")]
    )

    list(inner.cv = cv.pearsons, pearson = pearson)
    })
}


##############
### RANGER ###
##############

library(ranger)

RANGER_autoparams <- function(d.train) {
    autoparams(start = list(mtry = max(floor((ncol(d.train) - 3) / 3), 1),
                            nodesize = max(floor(nrow(d.train) * 4/5 * 0.66), 5),
                            ntree = 500,
                            maxnodes = "NULL"), # NULL here means the default
               range = list(mtry = c(20, (ncol(d.train) - 3)),
                            nodesize = c(2, nrow(d.train)),
                            ntree = c(1, 1000),
                            maxnodes = NULL),
               max_iter = c(5, 4, 3, 0),
               loss = function(aggr) {
                   r <- aggr$pred_perf_pval
                   if (is.na(r) || is.null(r))
                       return(1)
                   else
                       return(r)
               })
}

# Call randomForest on the specified data and return the pearson's
# correlation data.table() with all the genes
RANGER_model <- function(d.train, d.test, params, strata = NULL, return.fit = is.null(d.test), verbose = FALSE) {
    if (!is.null(strata))
        stop("Strata are not supported.")

    if (params$nodesize == "NULL")
        params$nodesize <- NULL
    if (params$maxnodes == "NULL")
        params$maxnodes <- NULL
    params <- lapply(params, as.numeric)

    # Grow the trees
    rf.fit <- ranger(
        EXPRESSION ~ .,
        data = d.train[, !c("GENE", "IID")],
        num.trees = params$ntree,
        mtry = params$mtry,
        min.node.size = params$nodesize,
        max.depth = params$maxnodes,
        splitrule = "maxstat",
        quantreg = TRUE,
        # seed = 2019,
        num.threads = 1,
        verbose = verbose
    )

    if (return.fit)
        return(rf.fit)

    # Analyse the correlation
    pred <- predict(rf.fit, d.test[, !c("GENE", "IID", "EXPRESSION")])$predictions
    pearson <- compute_pearson(d.test$EXPRESSION, pred, by = d.test$GENE)

    # Add info about the parameters and the mse/rsq of the model
    if (is.null(params$nodesize)) params$nodesize <- "NULL"
    if (is.null(params$maxnodes)) params$maxnodes <- "NULL"
    cbind(
        ntree = params$ntree,
        maxnodes = params$maxnodes,
        mtry = params$mtry,
        nodesize = params$nodesize,
        data.name = pearson$data.name,
        mse = cvmse(d.test$EXPRESSION, pred),
        rsq = cvrsq(d.test$EXPRESSION, pred),
        pearson[, c("parameter.df", "estimate.cor", "p.value")]
    )
}


#############
### RIDGE ###
#############

library(glmnet)

RIDGE_autoparams <- function(d.train) {
    autoparams(start = list(alpha = 0,
                            lambda = 2^seq(-7, 25, length = 100)),
               range = list(alpha = NULL,
                            lambda = NULL),
               max_iter = 0,
               loss = function(aggr) {
                   r <- aggr$pred_perf_pval
                   if (is.na(r) || is.null(r))
                       return(1)
                   else
                       return(r)
               })
}

RIDGE_model <- function(d.train, d.test, params.grid = NULL, params.auto = NULL, n.folds = 10,
                      folds.strata = NULL, model.strata = NULL, return.fit = FALSE, threads = 1, verbose = FALSE) {
    # Validate the arguments
    if (!is.null(params.grid))
        stop("params.grid is not supported yet.")
    if (!is.null(folds.strata))
        stop("Folds strata are not supported.")
    if (!is.null(model.strata))
        stop("Model strata are not supported.")
    if (is.null(params.grid) && is.null(params.auto))
        stop("At least one of params.grid and params.auto must be specified.")
    if (!is.null(params.grid) && !is.null(params.auto))
        stop("At most one of params.grid and params.auto can be specified.")

    params <- list(lambda = params.auto@best$lambda,
                   alpha = params.auto@best$alpha)

    glm.fit <- cv.glmnet(
        x = as.matrix(d.train[, !c("GENE", "IID", "EXPRESSION")]),
        y = d.train$EXPRESSION,
        nfolds = n.folds,
        alpha = params$alpha,
        lambda = params$lambda,
        standardize = FALSE, # the variables are already pretty much the same thing...
        keep = TRUE
    )

    if (return.fit)
        return(glm.fit)

    cv.pearsons <- lapply(seq_along(glm.fit$lambda), function(p) {
        pearson <- compute_pearson(d.train$EXPRESSION, glm.fit$fit.preval[, p], by = d.train$GENE)
        cbind(
            alpha = params$alpha,
            lambda = glm.fit$lambda[p],
            data.name = pearson$data.name,
            mse = cvmse(d.train$EXPRESSION, glm.fit$fit.preval[, p]),
            rsq = cvrsq(d.train$EXPRESSION, glm.fit$fit.preval[, p]),
            pearson[, c("parameter.df", "estimate.cor", "p.value")]
        )
    })
    cv.lambda <- glm.fit$lambda[which.min(sapply(cv.pearsons, function(p) mean(p$p.value)))]
    cv.pearsons <- do.call("rbind", cv.pearsons)

    ypred <- predict(glm.fit, newx = as.matrix(d.test[, !c("GENE", "IID", "EXPRESSION")]), s = cv.lambda)
    pearson <- compute_pearson(d.test$EXPRESSION, ypred[, 1], by = d.test$GENE)
    pearson <- cbind(
        alpha = params$alpha,
        lambda = cv.lambda,
        data.name = pearson$data.name,
        mse = cvmse(d.test$EXPRESSION, ypred),
        rsq = cvrsq(d.test$EXPRESSION, ypred),
        pearson[, c("parameter.df", "estimate.cor", "p.value")]
    )

    list(inner.cv = cv.pearsons, pearson = pearson)
}

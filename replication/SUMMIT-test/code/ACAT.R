ACAT <- function(Pvals, Weights = NULL) {
    if (sum(is.na(Pvals)) > 0) {
        stop("Cannot have NAs in the p-values!")
    }
    if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0){
        stop("P-values must be between 0 and 1!")
    }
    is.zero <- (sum(Pvals == 0) >= 1)
    is.one  <- (sum(Pvals == 1) >= 1)
    if (is.zero && is.one) {
        return(-1)
    }
    if (is.zero) {
        return(0)
    }
    if (is.one) {
        return(1)
    }

    if (is.null(Weights)) {
        Weights <- rep(1 / length(Pvals), length(Pvals))
    } else if (length(Weights) != length(Pvals)) {
        stop("The length of weights should be the same as that of the p-values")
    } else if (sum(Weights < 0) > 0){
        stop("All the weights must be positive!")
    } else {
        Weights <- Weights / sum(Weights)
    }

    is.small <- (Pvals < 1e-16)
    if (sum(is.small) == 0){
        cct.stat <- sum(Weights * tan((0.5 - Pvals) * pi))
    } else {
        cct.stat <- sum((Weights[is.small] / Pvals[is.small]) / pi)
        cct.stat <- cct.stat + sum(Weights[!is.small] * tan((0.5 - Pvals[!is.small]) * pi))
    }

    if (cct.stat > 1e15){
        pval <- (1 / cct.stat) / pi
    } else {
        pval <- pcauchy(cct.stat, lower.tail = F)
    }

    return(pval)
}
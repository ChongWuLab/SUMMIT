# Package
suppressMessages(library(BEDMatrix))
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(genio)
library(lassosum)
suppressMessages(library(Rcpp))

library(optparse)

setwd("/gpfs/research/chongwu/zichenzhang/SUMMIT-test/")

###########
# Options #
###########

# Options
option_list <- list(
    make_option("--h2_e", type = "numeric", default = NA, action = "store", help = "Heritability of expression"
    ),
    make_option("--h2_p", type = "numeric", default = NA, action = "store", help = "Heritability of phenotype"
    ),
    make_option("--p_causal", type = "numeric", default = NA, action = "store", help = "Percentage of causal SNPs"
    ),
    make_option("--n", type = "numeric", default = NA, action = "store", help = "Sample size"
    ),
    make_option("--sumstats", type = "logical", default = FALSE, action = "store", help = "Use summary-level data or not"
    ),
    make_option("--gene_ENSG", type = "character", default = FALSE, action = "store", help = "Gene's Ensembl ID"
    ),
    make_option("--UKB", type = "logical", default = FALSE, action = "store", help = "Using UKBioBank genotype data or not"
    ),
    make_option("--folder_output", type = "character", default = "", action = "store", help = "Folder to store the results"
    ),
    make_option("--t1e", type = "logical", default = FALSE, action = "store", help = "Explore Type-I error rate or not"
    ),
    make_option("--seed", type = "numeric", default = 0, action = "store", help = "Seed"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Pass options to variables
h2.e      <- opt$h2_e
h2.p      <- opt$h2_p
causal    <- opt$p_causal
n         <- opt$n
do.ss     <- opt$sumstats
gene.ENSG <- opt$gene_ENSG
UKB       <- opt$UKB
dir.out   <- opt$folder_output
T1E       <- opt$t1e

seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) + (opt$seed - 1) * 1000

#############
# Parameter #
#############

list.method <- c("MCP", "LASSO", "ElNet", "MNet", "SCAD")
s.array     <- 0.1 * (1:9)

# Exterior executables
path.Plink <- "plink"

# Get basic information
load(paste0("summary-statistics/eQTLGen/", gene.ENSG, ".RData"))

chr      <- ss$SNPChr[1]
list.SNP <- ss$SNP

##########################
# User-defined functions #
##########################

# ACAT
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

# Cpp functions
suppressMessages(sourceCpp("code/DetectAndSolve.cpp")) 

# FindOptimalResult
FindOptimalResult <- function(input, response.true) {
    index.optimal <- 0
    max           <- -1
    out           <- list()

    for (t in 1:length(input)) {
        if (norm(input[[t]][[1]], type = "2") <= 1) {
            # Compute the predicted response
            response.pred <- genotype.tune %*% input[[t]][[1]]

            # Do the regression
            reg <- summary(lm(response.true ~ response.pred))

            # Keep the optimal result
            if (reg$adj.r.sq > max) {
                max            <- reg$adj.r.sq
                result.optimal <- list(input[[t]][[1]], input[[t]][[2]])
                t.optimal      <- t
            }
        }
    }

    out[[1]] <- result.optimal
    out[[2]] <- t.optimal

    return(out)
}

# FindOptimalResult.Lassosum
FindOptimalResult.Lassosum <- function(response.true) {
    index.optimal <- 0
    max           <- -1
    out           <- list()

    for (t in 1:ncol(temp.weight)) {
        if (norm(temp.weight[, t], type = "2") <= 1) {
            # Compute the predicted response
            response.pred <- genotype.tune %*% temp.weight[, t]

            # Do the regression
            reg <- summary(lm(response.true ~ response.pred))

            # Keep the optimal result
            if (reg$adj.r.sq > max) {
                max            <- reg$adj.r.sq
                result.optimal <- temp.weight[, t]
                t.optimal      <- t
            }
        }
    }

    out[[1]] <- result.optimal
    out[[2]] <- t.optimal

    return(out)
}

# PatchUp
source("code/PatchUp.R")

# Standardize
Standardize <- function(M) {
    # Centralize
    M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)

    # Standardize
    M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")

    return(M)
}

# Marginal Z-score
MarginalZ <- function(genotype) {
    test <- summary(lm(phenotype ~ genotype))
    if (nrow(test$coef) == 1) {
        output <- 0
    } else {
        output <- test$coef[2, "Estimate"] / test$coef[2, "Std. Error"]
    }

    return(output)
}

####################################
# Extract the genotype information #
####################################

if (UKB) {
    temp.outdir <- "/gpfs/home/zz17/current-project/temp/"
    temp.snpfile <- paste0(temp.outdir, gene.ENSG, ".txt")

    if (!file.exists(paste0("/gpfs/home/zz17/current-project/temp/", gene.ENSG, ".bed"))) {
        # Export a SNP list for Plink
        write.table(list.SNP, file = temp.snpfile, quote = FALSE, col.names = FALSE, row.names = FALSE)

        # Extract genotype from UKBioBank datasets per the SNP list
        command <- path.Plink
        command <- paste0(command, " --bed /gpfs/research/chongwu/shared/UKBiobank/genetic/processed_white_british_", chr, ".bed")
        command <- paste0(command, " --bim /gpfs/home/zz17/current-project/temp/processed_white_british_", chr, ".bim")
        command <- paste0(command, " --fam /gpfs/research/chongwu/shared/UKBiobank/genetic/processed_white_british_", chr, ".fam")
        command <- paste0(command, " --extract ", temp.snpfile)
        command <- paste0(command, " --memory 1000")
        command <- paste0(command, " --make-bed")
        command <- paste0(command, " --out ", temp.outdir, gene.ENSG)

        quiet(
            system(command, intern = TRUE)
        )
    }

    # Extract genotype information
    quiet(
        seq.ref <- BEDMatrix(paste0(temp.outdir, gene.ENSG), simple_names = TRUE)
    )

    # Update the SNP list
    list.SNP <- colnames(seq.ref)
    ss       <- data.frame(list.SNP)
    size     <- length(list.SNP)

    names(ss) <- "SNP"

    genotype <- seq.ref[1:(n + 369 + 10000 + 10000), ] %>% PatchUp() %>% Standardize()
} else {
    set.seed(123)
    temp <- runif((n + 369 + 10000 + 10000) * size)
    a <- which(temp <= 0.33)
    b <- which(temp <= 0.66 & temp > 0.33)
    c <- which(temp > 0.66)
    temp[a] <- 0
    temp[b] <- 1
    temp[c] <- 2
    genotype <- matrix(temp, nrow = (n + 369 + 10000 + 10000), ncol = size)
}

# Split the data into four parts
genotype.train <- genotype[1:n, ]
genotype.tune  <- genotype[(n + 1):(n + 369), ]
genotype.test  <- genotype[(n + 369 + 1):(n + 369 + 10000), ]
genotype.ref   <- genotype[(n + 369 + 10000 + 1):(n + 369 + 10000 + 10000), ]

##############################
# Main iteration starts here #
##############################

#################
# Simulate data #
#################

# Set seed
set.seed(seed)

# Generate w
w <- numeric(size)
w[sample(1:size, floor(size * causal))] <- rnorm(floor(size * causal), 0, 1)

# Rescale w
w <- w * sqrt(h2.e / as.numeric(var(genotype %*% w)))

# Generate noise.e
noise.e <- rnorm(nrow(genotype), 0, sqrt(1 - h2.e))

# Generate expression levels
expression <- genotype %*% w + noise.e

##########################
# The original LD matrix #
##########################

if (do.ss) {
    X <- genotype.ref
} else {
    X <- genotype.train
}

X <- Standardize(X)
matrix.LD <- t(X) %*% X / nrow(X)

############
# r vector #
############

# Compute r vector
r <- numeric(size)
p <- numeric(size)
z <- numeric(size)

for (asd in 1:size) {
    model.temp <- summary(lm(expression[1:n] ~ genotype.train[, asd]))
    if (nrow(model.temp$coef) != 1) {
        z.temp <- model.temp$coef[2, "Estimate"] / model.temp$coef[2, "Std. Error"]
        r[asd] <- z.temp / sqrt(n - 1 + z.temp ^ 2)
        p[asd] <- model.temp$coef[2, "Pr(>|t|)"]
        z[asd] <- z.temp
    } else {
        r[asd] <- 0
        p[asd] <- 1
        z[asd] <- 0
    }
}

# Handle p-valuels that are too small
p[p == 0] <- 2e-16

##################
# Main iteration #
##################

result <- list()

for (sdf in 1:5) {
    method <- list.method[sdf]

    for (k in 1:length(s.array)) {
        s <- s.array[k]

        #########################################
        # Constructing the initial lambda array #
        #########################################

        # Get the maximum for lambda
        z.temp <- numeric(size)
        for (m in 1:size) {
            z.temp[m] <- abs(r[m])
        }

        lambda.max   <- max(z.temp)
        lambda.min   <- lambda.max * 1E-3
        lambda.array <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 10000)
        rm(z.temp)

        ###########################################
        # Explore the possible minimum for lambda #
        ###########################################

        problem <- numeric(10000)
        start   <- 1
        end     <- 10000

        max.iteration   <- 500

        threshold <- 100
        alpha     <- 0.5

        if (method == "SCAD") {
            gamma <- 3.7
        } else {
            gamma <- 3
        }

        dummy.detect <- TRUE

        while (dummy.detect) {
            # Stop if all lambdas are OK
            if (sum(problem == rep(1, 10000)) == 10000) {
                print("All lambdas check out!")
                break
            }

            start <- min(which(problem == 0))
            end   <- max(which(problem == 0))

            try    <- ceiling((start + end) /2)
            lambda <- lambda.array[try]

            if (method == "MCP") {
                problem[try] <- MCPDetect(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold, try)
            }
            if (method == "LASSO") {
                problem[try] <- ElNetDetect(r, matrix.LD, lambda.array, s, 1, max.iteration, threshold, try)
            }
            if (method == "ElNet") {
                problem[try] <- ElNetDetect(r, matrix.LD, lambda.array, s, 0.5, max.iteration, threshold, try)
            }
            if (method == "MNet") {
                problem[try] <- MNetDetect(r, matrix.LD, lambda.array, s, alpha, gamma, max.iteration, threshold, try)
            }
            if (method == "SCAD") {
                problem[try] <- SCADDetect(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold, try)
            }

            if (problem[try] == 1) {
                problem[1:try] <- 1
            } else {
                problem[try:10000] <- 2
            }

            if (problem[try] == 1 & problem[min(10000, try + 1)] == 2) {
                print(paste0("Problem detected for lambda=", lambda.array[try + 1], ". Finer interval has been saved!"))
                lambda.min   <- lambda
                dummy.detect <- FALSE
            }

            if (problem[try] == 2 & problem[try - 1] == 1) {
                print(paste0("Problem detected for lambda=", lambda.array[try], ". Finer interval has been saved!"))
                lambda.min   <- lambda.array[try - 1]
                dummy.detect <- FALSE
            }
        }

        # Update the lambda array
        lambda.array <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 100)

        #######################################################################
        # Optimizing with the chosen penalty (using the updated lambda array) #
        #######################################################################

        #####################
        # Iterate by lambda #
        #####################

        beta <- t(matrix(0, nrow = 1, ncol = size))

        max.iteration <- 500

        threshold <- 1e-6
        alpha     <- 0.5

        if (method == "SCAD") {
            gamma <- 3.7
        } else {
            gamma <- 3
        }

        if (method == "MCP") {
            res <- MCP(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold)
        }
        if (method == "LASSO") {
            res <- ElNet(r, matrix.LD, lambda.array, s, 1, max.iteration, threshold)
        }
        if (method == "ElNet") {
            res <- ElNet(r, matrix.LD, lambda.array, s, 0.5, max.iteration, threshold)
        }
        if (method == "MNet") {
            res <- MNet(r, matrix.LD, lambda.array, s, alpha, gamma, max.iteration, threshold)
        }
        if (method == "SCAD") {
            res <- SCAD(r, matrix.LD, lambda.array, s, gamma, max.iteration, threshold)
        }

        for (dfg in 1:100) {
            result[[length(result) + 1]] <- list(res[, dfg], lambda.array[dfg])
        }
    }
}

weight <- matrix(0, nrow = size, ncol = 8)
r2     <- numeric(8)

################################
# SUMMIT: select optimal model #
################################

result.optimal <- list()

# Get the true response vector
response.true <- expression[(n + 1):(n + 369)]

# Find the best model for each method
for (a in 1:5) {
    temp <- FindOptimalResult(input = result[(900 * a - 899):(900 * a)], response.true)

    weight[, a] <- temp[[1]][[1]]
}

##################################
# SUMMIT: test the optimal model #
##################################

# Get the response vector
response.true <- expression[(n + 369 + 1):(n + 369 + 10000)]

# Test
for (a in 1:5) {
    response.pred <- genotype.test %*% weight[, a]
    reg           <- summary(lm(response.true ~ response.pred))
    r2[a]         <- reg$adj.r.squared
}

############
# Lassosum #
############

# Correlation
cor <- p2cor(
    p = p,
    n = n,
    sign = r
)

# Assemble information for Lassosum
if (UKB) {
    if (do.ss) {
        ref.Lassosum <- seq.ref[(n + 369 + 10000 + 1):(n + 369 + 10000 + 10000), ]
    } else {
        ref.Lassosum <- seq.ref[1:n, ]
    }
} else {
    set.seed(123)
    temp <- runif((n + 369 + 10000 + 10000) * size)
    a <- which(temp <= 0.33)
    b <- which(temp <= 0.66 & temp > 0.33)
    c <- which(temp > 0.66)
    temp[a] <- 0
    temp[b] <- 1
    temp[c] <- 2
    temp <- matrix(temp, nrow = (n + 369 + 10000 + 10000), ncol = size)
    set.seed(seed)

    if (do.ss) {
        ref.Lassosum <- temp[(n + 369 + 10000 + 1):(n + 369 + 10000 + 10000), ]
    } else {
        ref.Lassosum <- temp[1:n, ]
    }
}

write_plink(
    file = paste0("data/", gene.ENSG, seed, "-Lassosum"),
    t(ref.Lassosum),
    bim = NULL,
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
)

bim <- paste0("data/", gene.ENSG, seed, "-Lassosum.bim") %>% fread() %>% as.data.frame()
fam <- paste0("data/", gene.ENSG, seed, "-Lassosum.fam") %>% fread() %>% as.data.frame()

# Run Lassosum
out <- lassosum.pipeline(
    cor = cor,
    chr = bim$V1,
    pos = bim$V4,
    A1 = bim$V5,
    A2 = bim$V6,
    ref.bfile = paste0("data/", gene.ENSG, seed, "-Lassosum"),
    LDblocks = "EUR.hg19",
    max.ref.bfile.n = 40000
)

temp.weight <- cbind(out$beta$`0.2`, out$beta$`0.5`) %>% cbind(., out$beta$`0.9`) %>% cbind(., out$beta$`1`)

temp.ss <- out$sumstats

match <- match(bim$V4, temp.ss$pos)

temp.ss <- temp.ss[match, ]

temp.weight <- temp.weight[match, ]

# Tune: Get the true response vector
response.true <- expression[(n + 1):(n + 369)]

# Tune: Find the best model
result.optimal <- FindOptimalResult.Lassosum(response.true)

# Tune: Little tweak to make beta vector an one-column matrix
result.optimal[[1]] <- as.matrix(result.optimal[[1]])

# Test: Get the response vectors
response.true <- expression[(n + 369 + 1):(n + 369 + 10000)]
response.pred <- genotype.test %*% result.optimal[[1]]

# Test: Regression
reg <- summary(lm(response.true ~ response.pred))

weight[, 6] <- result.optimal[[1]]
r2[6]       <- reg$adj.r.squared

###############
# TWAS-fusion #
###############

identifier <- seed

write_plink(
    file = paste0("data/", gene.ENSG, seed, "-fusion"),
    t(seq.ref[1:670, ]),
    bim = NULL,
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
)

temp <- paste0("data/", gene.ENSG, seed, "-fusion.fam") %>% fread() %>% as.data.frame()
temp[, 6] <- expression[1:670]

write.table(
    temp,
    file = paste0("data/", gene.ENSG, seed, "temp-fusion.fam"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)

file.copy(paste0("data/", gene.ENSG, seed, "-fusion.bed"), paste0("data/", gene.ENSG, seed, "temp-fusion.bed"))
file.copy(paste0("data/", gene.ENSG, seed, "-fusion.bim"), paste0("data/", gene.ENSG, seed, "temp-fusion.bim"))

command.twas <- "Rscript software/TWAS-fusion/fusion_twas-master/FUSION.compute_weights.R"
command.twas <- paste0(command.twas, " --bfile data/", gene.ENSG, seed, "temp-fusion")
command.twas <- paste0(command.twas, " --tmp data/TEMP-", identifier)
command.twas <- paste0(command.twas, " --out data/twas-", identifier)
command.twas <- paste0(command.twas, " --models top1,lasso,enet")
command.twas <- paste0(command.twas, " --PATH_plink software/plink")
command.twas <- paste0(command.twas, " --PATH_gcta software/gcta_nr_robust")
command.twas <- paste0(command.twas, " --hsq_set 0.1")
command.twas <- paste0(command.twas, " --crossval 5")

system(command.twas)

load(paste0("data/twas-", identifier, ".wgt.RDat"))
wgt.matrix[is.na(wgt.matrix)] <- 0

r2.twas <- numeric(3)
for (i in 1:3) {
    # Test
    reg        <- summary(lm(expression[(n + 369 + 1):(n + 369 + 10000)] ~ genotype.test %*% wgt.matrix[, i]))
    r2.twas[i] <- reg$adj.r.squared
}

weight[, 7] <- wgt.matrix[, which.max(r2.twas)]
weight[, 8] <- wgt.matrix[, 3]

r2[7] <- max(r2.twas, na.rm = TRUE)
r2[8] <- r2.twas[3]

colnames(weight) <- c(list.method, "Lassosum", "TWAS-fusion", "PrediXcan")
names(r2)        <- c(list.method, "Lassosum", "TWAS-fusion", "PrediXcan")

#######################
# Association studies #
#######################

runs <- 20

# Out
out <- matrix(0, nrow = runs, ncol = 9)

for (rep in 1:runs) {
    if (T1E) {
        # Generate beta
        beta <- 0

        # Generate phenotype levels
        phenotype <- rnorm(n, 0, 1)

        # Create ss
        Z <- apply(genotype.train, 2, MarginalZ)
    } else {
        # Generate beta
        beta <- rnorm(1, 0, 1)

        # Rescale beta
        beta <- beta * sqrt(h2.p / as.numeric(var(expression[1:n] * beta)))

        # Generate noise.p
        noise.p <- rnorm(n, 0, sqrt(1 - h2.p))

        # Generate phenotype levels
        phenotype <- expression[1:n] * beta + noise.p

        # Create ss
        Z <- apply(genotype.train, 2, MarginalZ)
    }

    # Iteration by method
    p <- numeric(9)

    for (w in 1:ncol(weight)) {
        # Settings
        weights <- weight[, w]

        # Skip if weight is a zero vector
        if (sum(weights) == 0) {
            p[w] <- NA
            next
        }

        # Keep the non-zero components of weights vector
        keep <- (weights != 0)
        weights <- weights[keep]

        # Compute TWAS z-score, r2, and p-value
        z.twas  <- as.numeric(weights %*% Z[keep])
        r2.twas <- as.numeric(weights %*% matrix.LD[keep, keep] %*% weights)

        p[w] <- 2 * (pnorm(abs(z.twas / sqrt(r2.twas)), lower.tail = F))
    }

    index.ACAT <- !is.na(r2[1:5]) & (r2[1:5] > 0)
    if (sum(index.ACAT) != 0) {
        p[9] <- ACAT(p[1:5][index.ACAT], r2[1:5][index.ACAT] / sum(r2[1:5][index.ACAT]))
    } else {
        p[9] <- NA
    }
    out[rep, ] <- p
}

out <- as.data.frame(out)

names(out) <- c(colnames(weight), "ACAT")

########
# Save #
########

# Create folders
dir.create(dir.out, showWarnings = FALSE)

# Save
save(
    weight,
    r2,
    h2.e,
    h2.p,
    causal,
    n,
    do.ss,
    seed,
    out,
    T1E,
    file = paste0(dir.out, "/", seed, ".RData")
)
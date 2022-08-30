#!/usr/bin/env Rscript
array  <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id.job <- as.numeric(array)

# Packages
library(BEDMatrix)
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(optparse)

# Options
option_list <- list(
    make_option("--models",   type = "character", default = FALSE, action = "store", help = "Path of the model folder"
    ),
    make_option("--path.ref", type = "character", default = FALSE, action = "store", help = "Path of the reference panel plus the prefixes"
    ),
    make_option("--path.ss",  type = "character", default = FALSE, action = "store", help = "Path of the summary statistics"
    ),
    make_option("--trait"   , type = "character", default = FALSE, action = "store", help = "Name of the summary statistics"
    ),
    make_option("--path.out", type = "character", default = FALSE, action = "store", help = "Path of the output file plus the prefixes"
    ),
    make_option("--parallel", type = "numeric",   default = 100,   action = "store", help = "The number of parallel instances"
    )    
)

opt <- parse_args(OptionParser(option_list = option_list))

path.model <- opt$models
path.ref   <- opt$path.ref
path.ss    <- opt$path.ss
trait      <- opt$trait
path.out   <- opt$path.out
pieces     <- opt$parallel

# User-defined functions
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

# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

# Object to store the output
out <- data.frame(
    gene_symbol = character(),
    gene_id     = character(),
    chromosome  = numeric(),
    model_best  = character(),
    r2_test     = numeric(),
    p_SCAD      = numeric(),
    p_ElNet     = numeric(),
    p_LASSO     = numeric(),
    p_MCP       = numeric(),
    p_MNet      = numeric(),
    z_SCAD      = numeric(),
    z_ElNet     = numeric(),
    z_LASSO     = numeric(),
    z_MCP       = numeric(),
    z_MNet      = numeric(),
    p_ACAT      = numeric(),
    gene_pos    = numeric(),
    runtime     = numeric(),
    stringsAsFactors = FALSE
)

# Total number of jobs should be 16884
files <- dir(path.model)

size.pieces <- ceiling(length(files) / pieces)

# Main iteration
for (gene.index in (1 + (id.job - 1) * size.pieces):(min(length(files), id.job * size.pieces))) {
    print(gene.index)

    # Start tracking runtime
    time.start <- proc.time()[3]

    # The vector to store all the updates in this iteration
    update <- rep(NA, 18)

    # Load weight
    load(paste0(path.model, files[gene.index]))
    update[5] <- max(out.r2[1, ], na.rm = TRUE)
    if (update[5] <= 0.005) {
        next
    }

    # Get the chromosome
    chr <- snps$SNPChr[1]

    # Read the summary statistics file
    if (file.exists(paste0(path.ss, trait, "-", chr, ".sumstats"))) {
        ss <- paste0(path.ss, trait, "-", chr, ".sumstats") %>% fread() %>% as.data.frame()
    } else {
        next
    }
    names(ss) <- colnames(ss) %>% tolower()
    ss["ID"] <- paste0(ss$chr, "_", ss$pos)

    # Load the reference panel
    quiet(
        bim.ref <- as.data.frame(fread(paste0(path.ref, chr, ".bim")))
    )
    quiet(
        genotype.ref <- BEDMatrix(paste0(path.ref, chr), simple_names = TRUE)
    )

    names(bim.ref)[c(2, 5, 6)] <- c("SNP" ,"a1", "a2")

    # Find the best model
    model.best <- which.max(out.r2[1, ])

    # Create new identifier for reference panel
    bim.ref["ID"] <- paste0(bim.ref$V1, "_", bim.ref$V4)

    # Find the common snps in all three data sets
    list.common <- intersect(bim.ref$ID, snps$ID) %>% intersect(., ss$ID)

    # Skip if no common snps found
    if (length(list.common) == 0) {
        next
    }

    # Trim genotype.ref
    genotype.temp <- genotype.ref[, bim.ref$ID %in% list.common]
    bim.temp      <- bim.ref[bim.ref$ID %in% list.common, ]

    # Fix the NAs in reference panel
    if (sum(is.na(genotype.temp)) != 0) {
        genotype.temp <- PatchUp(genotype.temp)
    }

    # Trim ss
    ss.temp <- ss[ss$ID %in% list.common, ]

    # Trim snps and wgt.matrix
    index.temp <- snps$ID %in% list.common

    snps       <- snps[index.temp, ]
    out.weight <- out.weight[index.temp, ]
    if ("numeric" %in% class(out.weight)) {
        out.weight <- t(as.matrix(out.weight))
    }

    rm(index.temp)

    # New
    # Re-arrange every data sets
    m.1 <- match(ss.temp$ID, bim.temp$ID)
    m.2 <- match(ss.temp$ID, snps$ID)

    bim.temp      <- bim.temp[m.1, ]
    genotype.temp <- genotype.temp[, m.1]

    snps       <- snps[m.2, ]
    out.weight <- out.weight[m.2, ]

    # Align the mismatched alleles
    problem.1 <- ss.temp$a1 != bim.temp$a1
    genotype.temp[, problem.1] <- 2 - genotype.temp[, problem.1]

    problem.2 <- ss.temp$a1 != snps$a1
    # flipped <- snps$SNP[problem.2]
    out.weight[problem.2, ] <- -1 * out.weight[problem.2, ]

    # # Old
    # # Re-arrange every data sets
    # m.1 <- match(bim.temp$ID, ss.temp$ID)
    # m.2 <- match(bim.temp$ID, snps$ID)

    # ss.temp    <- ss.temp[m.1, ]
    # snps       <- snps[m.2, ]
    # out.weight <- out.weight[m.2, ]

    # # Align the mismatched alleles
    # problem <- ss.temp$a1 != snps$a1
    # ss.temp$a1 <- snps$a1
    # ss.temp$a2 <- snps$a2
    # ss.temp$beta[problem] <- -1 * ss.temp$beta[problem]
    # ss.temp$z[problem] <- -1 * ss.temp$z[problem]

    # Compute LD matrix
    genotype.temp <- scale(genotype.temp)
    matrix.LD  <- t(genotype.temp) %*% genotype.temp / (nrow(genotype.temp) - 1)

    # Catch: When there is only one row in wgt.matrix
    if ("numeric" %in% class(out.weight)) {
        out.weight <- out.weight %>% as.matrix() %>% t() %>% as.data.frame()
    }

    # Catch: Over-fitting
    for (qwe in 1:5) {
        if (max(abs(out.weight[, qwe])) >= 10) {
            out.weight[, qwe] <- 0
        }
    }

    # Iteration by method
    for (w in 1:ncol(out.weight)) {
        # Settings
        weights <- out.weight[, w]

        # Skip if weight is a zero vector
        if (sum(weights) == 0) {
            update[5 + w] <- NA
            next
        }

        # Keep the non-zero components of weights vector
        keep <- (weights != 0)
        print(sum(keep))
        weights <- weights[keep]
        #print(sum(names(weights) %in% flipped))

        # Compute TWAS z-score, r2, and p-value
        z.twas  <- as.numeric(weights %*% ss.temp$z[keep])
        r2.twas <- as.numeric(weights %*% matrix.LD[keep, keep] %*% weights)

        update[5 + w]  <- 2 * (pnorm(abs(z.twas / sqrt(r2.twas)), lower.tail = F))
        update[10 + w] <- z.twas
    }

    # ACAT
    check.na   <- !is.na(update[6:10])
    check.sign <- out.r2[1, ] > 0

    check.final <- check.na & check.sign

    if (sum(check.final) == 0) {
        update[16] <- NA
    } else {
        update[16] <- ACAT(update[6:10][check.final], out.r2[1, check.final] / sum(out.r2[1, check.final]))
    }

    # Stop tracking runtime
    time.end <- proc.time()[3]

    #################
    # Output format #
    #################
    # 1.     Gene symbol
    # 2.     ENSG id
    # 3.     Chromosome
    # 4.     Best model
    # 5.     Best model's R^2 on testing data (Skipped)
    # 6-10.  P-value of TWAS from weights constructed by (SCAD, ElNet, LASSO, MCP, and MNet)
    # 11-15. Z-score of TWAS from weights constructed by (SCAD, ElNet, LASSO, MCP, and MNet)
    # 16.    P-value from ACAT
    # 17.    Gene position
    # 18.    Runtime

    # Update
    update[1]  <- snps$GeneSymbol[1]
    update[2]  <- snps$Gene[1]
    update[3]  <- chr
    update[4]  <- names(model.best)
    update[17] <- snps$GenePos[1]
    update[18] <- time.end - time.start

    out[nrow(out) + 1, ] <- update
}

# Write the result
dir.create(paste0(path.out, trait))
write.table(out, file = paste0(path.out, trait, "/", trait, "-", id.job), row.names = FALSE, quote = FALSE)

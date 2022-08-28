# Packages
library(BEDMatrix)
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(optparse)

# Options
option_list <- list(
    make_option("--path.ref", type = "character", default = FALSE, action = "store", help = "Path of the reference panel plus the prefixes"
    ),
    make_option("--trait"   , type = "character", default = FALSE, action = "store", help = "Name of the trait"
    ),
    make_option("--path.out", type = "character", default = FALSE, action = "store", help = "Path of the output file plus the prefixes"
    ),
    make_option("--parallel", type = "numeric", default = FALSE, action = "store", help = "Path of the output file plus the prefixes"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

path.ref <- opt$path.ref
trait    <- opt$trait
path.out <- opt$path.out
parallel <- opt$parallel

# Modify this line based on your environment
w <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# ACAT
source("code/ACAT.R")

# PatchUp
source("code/PatchUp.R")

# Options for test runs
test <- FALSE
if (test) {
    path.ref <- "data/1000G.EUR.ALLSNP.QC.CHR"
    trait    <- "COVID-B2-V5"
    path.out <- "COVID-B2-V5"
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

# Main iteration
files <- dir("model")

index <- ((w - 1) * ceiling(length(files) / parallel) + 1):min(w * ceiling(length(files) / parallel), length(files))

for (gene.index in index) {
    # Start tracking runtime
    time.start <- proc.time()[3]

    # The vector to store all the updates in this iteration
    update <- rep(NA, 18)

    # Load weight
    load(paste0("model/", files[gene.index]))
    update[5] <- max(out.r2[1, ], na.rm = TRUE)
    if (update[5] <= 0.005) {
        next
    }

    # Get the chromosome
    chr <- snps$SNPChr[1]

    # Read the summary statistics file
    if (file.exists(paste0("summary-statistics/", trait, "-", chr, ".sumstats"))) {
        ss <- paste0("summary-statistics/", trait, "-", chr, ".sumstats") %>% fread() %>% as.data.frame()
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
    out.weight[problem.2, ] <- -1 * out.weight[problem.2, ]

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
        weights <- weights[keep]

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
dir.create(path.out)
write.table(out, file = paste0(path.out, "/", trait, "-", w), row.names = FALSE, quote = FALSE)

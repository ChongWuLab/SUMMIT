#########################
# Package and parameter #
#########################

# Package
library(BEDMatrix)
library(compiler)
library(data.table)
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
library(optparse)
library(Rcpp)

# Parameter
do.adjust <- TRUE
s.array   <- 0.1 * (1:9)

###########
# Options #
###########

# Options
option_list <- list(
    make_option("--name_batch", type = "character", default = FALSE, action = "store", help = "Name of this batch"
    ),
    make_option("--method"    , type = "character", default = FALSE, action = "store", help = "Method to use"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Pass the options to variables
batch  <- opt$name_batch
method <- opt$method

##########################
# User-defined functions #
##########################

# Cpp functions
suppressMessages(sourceCpp("code/DetectAndSolve.cpp"))

# FindOptimalResult
FindOptimalResult <- function(response.true) {
    # Find the optimal out of 900 results
    index.optimal <- 0
    max           <- -1

    for (t in 1:900) {
        if (norm(result[[t]][[1]], type = "2") <= 100) {
            # Compute the predicted response
            response.pred <- genotype.train %*% result[[t]][[1]]

            # Do the regression
            reg <- summary(lm(response.true ~ response.pred))

            # Keep only the optimal result
            if (reg$adj.r.sq > max) {
                max            <- reg$adj.r.sq
                result.optimal <- list(result[[t]][[1]], result[[t]][[2]], reg$adj.r.sq, reg$coef)
                t.optimal      <- t
            }
        }
    }

    # Compute the P-value of R square
    quiet(
        test <- cor.test(as.vector(response.true), as.vector(response.pred), method = "pearson", alternative = "greater")
    )
    result.optimal[[5]] <- test$p.value

    return(result.optimal)
}

FindOptimalResult <- cmpfun(FindOptimalResult)

# PatchUp
source("code/PatchUp.R")

PatchUp <- cmpfun(PatchUp)

# Standardize
Standardize <- function(M) {
    # Centralize
    M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)

    # Standardize
    M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")

    return(M)
}

Standardize <- cmpfun(Standardize)

# Translate
list.translation           <- readRDS("data/gencode.v26.hg19.genes.rds")
list.translation$gene_id   <- substr(list.translation$gene_id, 1, 15)
list.translation$gene_name <- as.character(list.translation$gene_name)

list.supplementary <- readRDS("data/list.supplementary.rds")

Translate <- function(gene.ENSG) {
    temp.1 <- list.translation[which(list.translation$gene_id == gene.ENSG), ]

    if (nrow(temp.1) == 0) {
        temp.2 <- list.supplementary[which(list.supplementary$ensembl_gene_id == gene.ENSG), ]

        if (nrow(temp.2) == 0) {
            gene.proper <- gene.ENSG
        } else {
            gene.proper <- temp.2[1, 3]
        }
    } else {
        gene.proper <- temp.1[1, 12]
    }

    if (gene.proper == "") {
        gene.proper <- gene.ENSG
    }

    return(gene.proper)
}

Translate <- cmpfun(Translate)

##############################
# Read GTEx-7's subject list #
##############################

load("data/list.7.RData")

##########################################################################
# Read GTEx-8's expression data (response.8.RData is already translated) #
##########################################################################

load("data/response.8.RData")

####################
# List of subjects #
####################

# Split subjects in GTEx-8 into two groups
list.8 <- rownames(response.8)

list.train <- list.8[which(list.8 %in% list.7)]
list.valid <- list.8[which(!(list.8 %in% list.7))]

####################################
# There are 19,227 genes in total! #
####################################

filenames <- dir("summary-statistics/eQTLGen")

for (i in 1:length(filenames)) {
    # Start keeping track of runtime
    time.start <- proc.time()[3]

    #################
    # Preprocess ss #
    #################

    load(paste0("summary-statistics/eQTLGen/", filenames[i]))

    # Get the name and chromosome of to-be-processed gene
    gene.ENSG   <- ss$Gene[1]
    gene.proper <- Translate(gene.ENSG)
    chr         <- ss$SNPChr[1]

    ##########################
    # Clumping and filtering #
    ##########################

    do.clumping <- FALSE

    # Obsolete
    if (do.clumping) {
        ###########
        # Extract #
        ###########

        # Create a temporary SNP list for Plink
        temp.outdir  <- "/gpfs/home/zz17/current-project/temp/"
        temp.snpfile <- paste0(temp.outdir, gene.ENSG, ".txt")
        write.table(ss[, 1], file = temp.snpfile, quote = FALSE, col.names = FALSE, row.names = FALSE)

        # Put together the command and execute
        command <- path.Plink
        command <- paste0(command)
        command <- paste0(command, " --bfile /gpfs/home/zz17/resource/1000-genome/sequence/1000G.EUR.ALLSNP.QC.CHR", chr)
        command <- paste0(command, " --chr ", chr)
        command <- paste0(command, " --extract ", temp.snpfile)
        command <- paste0(command, " --make-bed")
        command <- paste0(command, " --out ", temp.outdir, gene.ENSG)

        quiet(
            system(command, intern = TRUE)
        )

        ###########
        # Process #
        ###########

        # Create a temporary ss file and add a p-value column
        ss$pval     <- 2 * pnorm(abs(ss$Zscore), lower.tail = FALSE)
        temp.ssfile <- paste0(temp.outdir, gene.ENSG, "_ss.txt")
        write.table(ss, file = temp.ssfile, quote = FALSE, col.names = TRUE, row.names = FALSE)

        # Put together the command and execute
        command <- path.Plink
        command <- paste0(command, " --bfile ", temp.outdir, gene.ENSG)
        command <- paste0(command, " --clump ", temp.ssfile)
        command <- paste0(command, " --clump-field pval")
        command <- paste0(command, " --clump-kb 250")
        command <- paste0(command, " --clump-p1 1")
        command <- paste0(command, " --clump-r2 0.95")
        command <- paste0(command, " --clump-snp-field SNP")
        command <- paste0(command, " --out ", temp.outdir, gene.ENSG)

        quiet(
            system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
        )

        clumped.snp <- fread(paste0(temp.outdir, gene.ENSG, ".clumped"))
        clumped.snp <- as.data.frame(clumped.snp)
        clumped.snp <- clumped.snp[, "SNP"]

        # Clean up and update
        system(paste0("rm ", temp.outdir, gene.ENSG, "*"))

        ss <- ss[ss$SNP %in% clumped.snp, ]
    }

    ############
    # Quick QC #
    ############

    # Skip this round of iteration if ss is empty to begin with
    if (nrow(ss) == 0) {
        next
    }

    # Read the bim files of GTEx-8 and reference panel
    load("data/bim.RData")

    # Only keep SNPs that are simultaneously IN (reference panel, GTEx-8)
    ss <- ss[ss$SNP %in% bim.ref$V2 & ss$SNP %in% bim.8$V2, ]

    # Sieve out problematic SNPs
    list.1  <- bim.8$V2[duplicated(bim.8$V2)]
    list.2  <- bim.8$V2[nchar(bim.8$V5) > 1 | nchar(bim.8$V6) > 1]
    list.3  <- bim.ref$V2[duplicated(bim.ref$V2)]
    list.4  <- bim.ref$V2[nchar(bim.ref$V5) > 1 | nchar(bim.ref$V6) > 1]
    list.5  <- ss$SNP[duplicated(ss$SNP)]
    problem <- ss$SNP %in% list.1 | ss$SNP %in% list.2 | ss$SNP %in% list.3 | ss$SNP %in% list.4 | ss$SNP %in% list.5

    ss <- ss[!problem, ]

    rm(list.1)
    rm(list.2)
    rm(list.3)
    rm(list.4)
    rm(list.5)
    rm(problem)

    # Skip this round of iteration if ss has few row left
    if (nrow(ss) <= 1) {
        next
    }

    ########################################################
    # Make sure that the gene of the hour is in response.8 #
    ########################################################

    if (!(gene.proper %in% colnames(response.8))) {
        next
    }

    #################################################
    # Preprocess genotype matrix of reference panel #
    #################################################

    # Extract genotype information
    quiet(
        seq.ref <- BEDMatrix(paste0("data/1000G.EUR.ALLSNP.QC.CHR", chr), simple_names = TRUE)
    )

    genotype.ref <- seq.ref[, ss$SNP]
    rm(seq.ref)

    # Skip this round of iteration if genotype is empty
    if (ncol(genotype.ref) == 0) {
        next
    }

    # Compare to see if all (A1 and A2) pairs are well-aligned
    # Always use reference panel's A1/A2 as the standard order
    names(bim.ref)[2] <- "SNP"
    bim.temp <- subset(bim.ref, SNP %in% ss$SNP)
    ss.temp  <- left_join(ss, bim.temp, by = "SNP")
    problem  <- !(ss.temp$a1 == ss.temp$V5)

    if (sum(problem) != 0) {
        # Flip the problematic pairs
        ss$Zscore[problem] <- -1 * ss$Zscore[problem]

        ss$a1 <- ss.temp$V5
        ss$a2 <- ss.temp$V6
    }

    rm(bim.temp)
    rm(ss.temp)
    rm(problem)

    # Patch up the NAs and centralize genotype matrix
    genotype.ref <- PatchUp(genotype.ref)
    genotype.ref <- Standardize(genotype.ref)

    #######################################
    # Preprocess genotype matrix of GTEx8 #
    #######################################

    # The big genotype matrix
    load("data/genotype.8.RData")

    # Skip this round of iteration if genotype is empty
    if (ncol(genotype.8) == 0) {
        next
    }

    # Compare to see if all (a1 and a2) are well-aligned
    names(bim.ref)[2] <- "SNP"
    names(bim.8)[2]   <- "SNP"
    bim.temp1 <- subset(bim.ref, SNP %in% ss$SNP)
    bim.temp2 <- subset(bim.8, SNP %in% ss$SNP)
    bim.temp3 <- left_join(bim.temp1, bim.temp2, by = "SNP")
    problem   <- !(bim.temp3$V5.x == bim.temp3$V5.y)

    if (sum(problem) != 0) {
        genotype.8[, problem] <- 2 - genotype.8[, problem]
    }

    rm(bim.temp1)
    rm(bim.temp2)
    rm(bim.temp3)
    rm(problem)

    #######################################################
    # EXTREMELY ANNOYING ISSUE FOUND!!!!!!!!!!!!!!!!!!!!! #
    #######################################################

    # Subject number is NOT the same for all chromosomes!
    genotype.train <- genotype.8[which(rownames(genotype.8) %in% list.train), ]
    genotype.train <- scale(PatchUp(genotype.train))

    genotype.valid <- genotype.8[which(rownames(genotype.8) %in% list.valid), ]
    genotype.valid <- scale(PatchUp(genotype.valid))

    genotype.train[is.na(genotype.train)] <- 0
    genotype.valid[is.na(genotype.valid)] <- 0

    ##########################
    # The original LD matrix #
    ##########################

    matrix.LD <- t(genotype.ref) %*% genotype.ref / nrow(genotype.ref)

    #########################
    # Additional adjustment #
    #########################

    if (do.adjust) {
        # Some given parameters
        cutoff        <- 0.001
        matrix.adjust <- matrix(0, nrow(ss), nrow(ss))
        n.effective   <- 11400
        n.size        <- 183

        temp.cof <- (-2 * n.effective) / n.size

        # Read distance information and create an identifer
        distance <- as.data.frame(fread(paste0("data/chr", chr, ".OMNI.interpolated_genetic_map")))
        names(distance) <- c("SNP", "position", "distance")

        # Merge and process
        ss.temp <- left_join(ss, distance, by = "SNP")

        for (col in 1:(nrow(ss.temp) - 1)) {
            for (row in (col + 1):nrow(ss.temp)) {
                entry.temp <- exp(temp.cof * abs(ss.temp[col, "distance"] - ss.temp[row, "distance"]))
                entry.temp <- entry.temp * (entry.temp >= cutoff)
                matrix.adjust[row, col] <- entry.temp
            }
        }

        matrix.adjust <- t(matrix.adjust) + matrix.adjust + diag(nrow(ss.temp))

        rm(distance)
        rm(ss.temp)
        rm(temp.cof)

        ###################
        # Final LD matrix #
        ###################

        # "*" stands for elementwise multiplication
        matrix.LD <- matrix.LD * matrix.adjust
    }

    ############
    # r vector #
    ############

    # Compute r vector
    r <- ss$Zscore / sqrt(ss$NrSamples - 1 + ss$Zscore ^2)

    # Get size
    size <- nrow(ss)

    ##########################
    # Iterate By parameter s #
    ##########################

    result <- list()

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

        ######################################
        # What do 0, 1, and 2 stand for?     #
        # 0: Not tested                      #
        # 1: tested and no problem detected  #
        # 2: tested and problem is detected  #
        ######################################

        problem <- numeric(10000)
        start   <- 1
        end     <- 10000

        max.iteration <- 300

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

        max.iteration <- 300

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

        for (b in 1:100) {
            result[[length(result) + 1]] <- list(res[, b], lambda.array[b])
        }
    }

    ############################################################################################################################
    #         -----.                   osssso            ---------..`                                                          #
    #        `ssssso                   osssso           `sssssssssssss+:`                                                      #
    #        `ssssso       `...`       osssso   ``      `ssssssssssssssss/        `...`                ```          `...`      #
    #        `ssssso   `/osssssss+:`   osssso:osssso/`  `sssss/  `.:osssss+   `:+sssssss+:`   +sssso.+sssso/`    -+sssssss+-   #
    #        `ssssso  /sssssssssssso-  ossssssssssssss- `sssss/     `ssssss` :sssssssssssss-  +ssssssssssssso  .ossso/::ossso` #
    #        `ssssso -sssss:``./sssss` osssso-``./sssss``sssss/      osssss..sssss:``./sssss. +sssss-``osssss` ossss/---:ssss+ #
    #        .ssssso /ssss+    `sssss. ossss/    `sssss.`sssss/    `/ssssso :sssso     sssss- +sssss   /sssss``sssssssssssssso #
    # .:/oso/osssss/ .sssss+--:osssso  osssss/--:osssso `ssssso+++ossssss+` `sssss+:-:+sssso` +sssss   /sssss` +ssss/`  `-`    #
    # .ossssssssss+`  .osssssssssss+`  ossssosssssssso. `ssssssssssssss+-    .osssssssssss+.  +sssss   /sssss` `/sssssssssso/` #
    #   -/+osso+/.      ./+oosoo+:`    /++++.-/+oo+/-   `++++++++//:-.`        .:+oosoo+:.    :+++++   :+++++`   `:/oosoo+/-`  #
    ############################################################################################################################

    #########################
    # Select optimal lambda #
    #########################

    #############################
    # Calculate goodness of fit #
    #############################

    ###############################################
    # Structure of the list result/result.optimal #
    ###############################################
    # index:                                      #
    # [[1]]  beta vector                          #
    # [[2]]  lambda                               #
    # [[3]]  r square adjusted (to be added)      #
    # [[4]]  regression coefficient               #
    # [[5]]  p-value of r square                  #
    # [[6]]  r-square projected by GTEx8 - GTEx7  #
    # [[7]]  runtime for the gene                 #
    # [[8]]  average sample size                  #
    # [[9]]  iteration index                      #
    # [[10]] correlation                          #
    ###############################################

    # Get the true response vector
    response.true <- response.8[list.train, which(colnames(response.8) == gene.proper)]

    # Figure out if response.true is just a vector or multiple vectors and find the optimal
    if (class(response.true) == "numeric") {
        result.optimal <- FindOptimalResult(response.true)
    } else {
        max <- -1
        for (asd in 1:ncol(response.true)) {
            temp <- FindOptimalResult(response.true[, asd])[[3]]
            if (temp > max) {
                max <- temp
                n.col <- asd
            }
        }

        response.true  <- response.true[, n.col]
        result.optimal <- FindOptimalResult(response.true)
    }

    # Compute the P-value of R square
    response.pred <- genotype.train %*% result.optimal[[1]]
    quiet(
        test <- cor.test(as.vector(response.true), as.vector(response.pred), method = "pearson", alternative = "greater")
    )
    result.optimal[[5]] <- test$p.value

    # Little tweak to make beta vector an one-column matrix
    result.optimal[[1]] <- as.matrix(result.optimal[[1]])

    ##############################
    # Validate model with GTex-8 #
    ##############################

    # Get the response vectors
    response.true <- response.8[list.valid, which(colnames(response.8) == gene.proper)]

    if (exists("n.col")) {
        response.true <- response.true[, n.col]
        rm(n.col)
    }

    response.pred <- genotype.valid %*% result.optimal[[1]]

    # Cross validation with subject that are NOT INCLUDED in GTEx7
    reg                 <- summary(lm(response.true ~ response.pred))
    result.optimal[[6]] <- reg$adj.r.squared

    # Keep track of runtime
    time.end <- proc.time()[3]
    result.optimal[[7]] <- time.end - time.start

    # Keep track of average sample size and miscellaneous items
    result.optimal[[8]]  <- mean(ss$NrSamples)
    result.optimal[[9]]  <- i
    result.optimal[[10]] <- cor(response.true, response.pred)

    #############################
    # Save (TWAS-fusion format) #
    #############################

    # cv.performance
    cv.performance           <- matrix(nrow = 2, ncol = 1)
    colnames(cv.performance) <- method
    rownames(cv.performance) <- c("rsq", "pval")
    cv.performance[1, 1]     <- reg$adj.r.squared
    quiet(
        test <- cor.test(as.vector(response.true), as.vector(response.pred), method = "pearson", alternative = "greater")
    )
    cv.performance[2, 1]     <- test$p.value

    # snps
    snps <- ss

    # wgt.matrix
    wgt.matrix           <- result.optimal[[1]]
    colnames(wgt.matrix) <- method
    rownames(wgt.matrix) <- ss$SNP

    dir.create(paste0("output_twas_", batch), showWarnings = FALSE)
    save.twas.chr <- paste0("output_twas_", batch, "/Whole_Blood.", gene.ENSG, ".wgt.RData")
    save(cv.performance, snps, wgt.matrix, file = save.twas.chr)
}

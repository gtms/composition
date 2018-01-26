## Given a concatenated expression matrix (with genes as rows and samples as
## columns) of at least two datasets; a factor describing each dataset; and a
## factor of clinical subtypes as derived by an external classifier, returns an
## integrated matrix of expression, the result of stacking ComBat outputs for
## the integration of all datasets split by molecular clinical subtypes.

## * doSubtypeCombat
doSubtypeCombat <- function(mtx, # concatenated expression matrix of at least
                                        # two datasets
                            d.fac, # factor of length ncol(mtx) describing
                                        # datasets, with at least two levels
                            subtypePreds.fac, # factor of subtype predictions,
                                        # of length ncol(mtx), with at least two
                                        # levels
                            nMin = 3) { # minimum number of samples required per
                                        # subtype-specific batch
    sampleInfo.dtb <- data.table(sn = colnames(mtx),
                                 b = d.fac,
                                 st = subtypePreds.fac)
    ## Some requirements have to be met for ComBat to run:
    ## 1. Each subtype-specific batch must have at least 2 samples (although we
    ## require nMin = 3 to be more stringent) to draw summary statistics from:
    ## remove any sample that contibutes to violate this condition.
    tableSampleInfo.dtb <- sampleInfo.dtb[, .N, by = .(b, st)]
    setkey(tableSampleInfo.dtb, b, st)
    setkey(sampleInfo.dtb, b, st)
    sn1 <- sampleInfo.dtb[tableSampleInfo.dtb[N >= nMin]][, sn]
    ## 2. Each subtype must encompass at least two batches, i.e., no given
    ## subtype can only integrate samples from one batch: remove all samples
    ## that contribute to violate this condition.
    ## tableSubtypeCrossBatch.dtb <- tableSampleInfo.dtb[, .N, by = st]
    setkey(sampleInfo.dtb, sn)
    tableSubtypeCrossBatch.dtb <- sampleInfo.dtb[sn1, .N, by = .(b, st)][, .N, by = st]
    setkey(tableSubtypeCrossBatch.dtb, st)
    setkey(sampleInfo.dtb, st)
    sn2 <- sampleInfo.dtb[tableSubtypeCrossBatch.dtb[N > 1]][, sn]
    ## intersect both sets of samples
    validSn <- intersect(sn1, sn2)
    ## define index of samples that meet required conditions
    idx2kp <- colnames(mtx) %in% validSn
    ## redefine list of subtype-specific batches only with samples that meet
    ## required conditions
    rdx.lst <- tapply(1:ncol(mtx[, idx2kp]),
                      subtypePreds.fac[idx2kp], function(x) {
                          list(d = mtx[, idx2kp][, x],
                               b = droplevels(d.fac[idx2kp][x]))
                      })
    ## remove potential dropped out subtype-specific batches
    rdx.lst <- rdx.lst[!sapply(rdx.lst, is.null)]
    ## perform subtype combat
    cbt.mtx <- do.call(cbind, lapply(rdx.lst, function(sbt) {
        ComBat(dat = sbt$d,
               batch = sbt$b)
    }))
    ## return integrated matrix, with removed samples as NA
    out.mtx <- matrix(NA,
                      ncol = ncol(mtx),
                      nrow = nrow(mtx),
                      dimnames = list(rownames(mtx),
                                      colnames(mtx)))
    out.mtx[, colnames(cbt.mtx)] <- cbt.mtx
    out.mtx
}

## test code
## 1. generate four datasets of 60 samples * 1000 features
## d.lst <- lapply(1:4, function(x) {
##     matrix(rnorm(60 * 1e3),
##            ncol = 60,
##            dimnames = list(sprintf("g%d", 1:1e3),
##                            sprintf("%ss%d", sprintf("d%d", x), 1:60)))
## })
## names(d.lst) <- sprintf("d%d", 1:4)

## 2. add corresponding batch effect to each of the four experiments
## batchEffect <- c(0, 2, -2, 5)
## batchEffect <- c(0, 2, -2, 5) * 1e2
## d.lst <- lapply(1:4, function(x) {
##     d.lst[[x]] <- d.lst[[x]] + matrix(rnorm(60 * 1e3,
##                                             mean = batchEffect[x],
##                                             sd = 5),
##                                       ncol = 60)
## })
## boxplot(do.call(cbind, d.lst))

## 3. randomly assign one of five subtypes to each dataset
## sbt.lst <- list()
## subtypes <- sprintf("sbt%d", 1:5)
## 3.1 d1 has all five subtypes in roughly same proportions
## sbt.lst[["d1"]] <- sample(subtypes, 60, replace = TRUE, prob = c(rep(1/5, 5)))
## 3.2 d2 lacks s5 and has an excess of s1
## sbt.lst[["d2"]] <- sample(subtypes, 60, replace = TRUE, prob = c(.4, .2, .2, .2, 0))
## 3.3 d3 lacks s1 and has an excess of s5
## sbt.lst[["d3"]] <- sample(subtypes, 60, replace = TRUE, prob = c(0, .2, .2, .2, .4))
## 3.4 d4 lacks s3 and only has one sample of s5
## sbt.lst[["d4"]] <- sample(subtypes, 60, replace = TRUE, prob = c(.3, .3, 0, .4, 0))
## sbt.lst[["d4"]][1] <- "sbt5"

## 4. concatenate
## c.mtx <- do.call(cbind, d.lst)
## sbt.fac <- ordered(do.call(c, sbt.lst), levels = subtypes)
## dsets.fac <- ordered(rep(sprintf("d%d", 1:4), each = 60),
##                      levels = sprintf("d%d", 1:4))

## 5. add subtype noise
## subtypeNoise <- c(0, 6, -3, 7, -5)
## subtypeNoise <- c(0, 6, -3, 7, -5) * 1e2
## splitBySbt.lst <- tapply(colnames(c.mtx), sbt.fac, function(clnms) {
##     c.mtx[, clnms]
## })
## splitBySbt.lst <- lapply(1:5, function(x) {
##     nSmpls <- ncol(splitBySbt.lst[[x]])
##     splitBySbt.lst[[x]] <- splitBySbt.lst[[x]] + matrix(rnorm(nSmpls * 1e3,
##                                                               mean = subtypeNoise[x],
##                                                               sd = 2),
##                                                         ncol = nSmpls)
## })
## mtx <- do.call(cbind, splitBySbt.lst)
## mtx <- mtx[, colnames(c.mtx)]
## boxplot(mtx, col = sbt.fac)

## 6. integrate with ComBat
## library(data.table)
## library(sva)
## cb.mtx <- ComBat(mtx,
##                  batch = dsets.fac)
## boxplot(cb.mtx, col = sbt.fac)
## tapply(colnames(cb.mtx), sbt.fac, function(smplNms) {
##     sd(as.numeric(cb.mtx[, smplNms]))
## })

## 7. integrate with SABC
## sabc.mtx <- doSubtypeCombat(mtx,
##                             dsets.fac,
##                             sbt.fac)
## boxplot(sabc.mtx, col = sbt.fac)
## tapply(colnames(sabc.mtx), sbt.fac, function(smplNms) {
##     sd(as.numeric(sabc.mtx[, smplNms]))
## })

## 8. integrate with ComBat with covariate
## cbCov.mtx <- ComBat(mtx,
##                     batch = dsets.fac,
##                     mod = model.matrix(~ as.numeric(sbt.fac)))
## boxplot(cbCov.mtx, col = sbt.fac)

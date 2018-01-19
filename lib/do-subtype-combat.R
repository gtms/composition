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

################################################################################
# Given a set of GBRS (Genotyping by RNA-Seq) haplotype reconstructions, compare
# them to the DOQTL MUGA based reconstructions.
#
# Daniel Gatti
# dan.gatti@jax.org
# Mar. 6, 2017
################################################################################
# NOTE: The reconstructions must be on the same marker grid. Use DOQTL function
# interpolate.markers() to change the grid.
# NOTE: The sample names must match exactly between the GBRS and MUGA data.
# Arguments: gbrs: 3-dimensional, numeric array containing the GBRS (8-state)
#                  haplotype probabilities. Dimensions must have names.
#                  rows: samples (sample IDs)
#                  columns: 8 founders (A-H)
#                  dim 3: markers (marker IDs)
#            muga: 3-dimensional, numeric array containing the GBRS (8-state)
#                  haplotype probabilities. Dimensions must have names.
#                  rows: samples (sample IDs)
#                  columns: 8 founders (A-H)
#                  dim 3: markers (marker IDs)
#
# Returns: data.frame with nrow(muga) rows and 5 columns:
#          muga.sample: The MUGA sample ID.
#          self.cor: the correlation of the MUGA and GBRS samples with 
#                    the same sample ID. NA if the sample does not occur in 
#                    both data sets.
#          best.gbrs.sample: The GBRS sample with the best correlation to 
#                    the MUGA sample.
#          cross.cor: The best correlation between the MUGA sample and all
#                     GBRS samples.
#          mismatch: TRUE if the MUGA and GBRS samples don't match. NA if the
#                    sample ID does not occur in both data sets.
# 
sample_mismatch = function(gbrs, muga) {

  if(missing(gbrs)) {
    stop("Please supply non-null gbrs data.")
  } else {
    if(length(dim(gbrs)) < 3) {
      stop("The gbrs data must be a 3 dimensional, numeric array.")
    } # if(length(dim(gbrs)) < 3)
  } # else

  if(missing(muga)) {
    stop("Please supply non-null muga data.")
  } else {
    if(length(dim(muga)) < 3) {
      stop("The muga data must be a 3 dimensional, numeric array.")
    } # if(length(dim(muga)) < 3)
  } # else

  # Make sure the dimensions match.
  if(ncol(muga) != ncol(gbrs)) {
    stop(paste0("The number of columns does not match between MUGA (",
         ncol(muga), ") and GBRS (", ncol(gbrs), ")."))
  } # if(ncol(muga) != ncol(gbrs))

  if(dim(muga)[3] != dim(gbrs)[3]) {
    stop(paste0("The number of markers does not match between MUGA (",
         dim(muga)[3], ") and GBRS (", dim(gbrs)[3], ")."))
  } # if(dim(muga)[3] != dim(gbrs)[3])

  # Make sure that the two data sets contain at least some sample in common.
  samples = intersect(rownames(muga), rownames(gbrs))
  if(length(samples) == 0) {
    stop(paste("The two data sets contain no sample names in common. Please make",
         "sure that there are some sample IDs that match between the two data sets."))
  } # if(length(samples) == 0)

  # Convert each dataset into a matrix.
  message("Converting matrices to vectors...")
  muga = matrix(as.vector(muga), nrow = nrow(muga), 
         ncol = ncol(muga) * dim(muga)[3], dimnames = list(rownames(muga), NULL))

  gbrs = matrix(as.vector(gbrs), nrow = nrow(gbrs), 
         ncol = ncol(gbrs) * dim(gbrs)[3], dimnames = list(rownames(gbrs), NULL))

  muga = t(muga)
  gbrs = t(gbrs)

  # Use the Pearson correlation to compare samples.
  # MUGA samples will be in rows. GBRS samples in columns.
  message("Calculating corrlations...")
  geno.cor = cor(muga, gbrs)

  # Get the GBRS sample with the maximum correlation for each MUGA sample.
  muga.max.cor = apply(geno.cor, 1, max)
  muga.max.cor.col = apply(geno.cor, 1, which.max)

  # Collate the results.
  message("Collating results...")
  result = data.frame(muga.sample = colnames(muga),
                      self.cor    = rep(NA, ncol(muga)),
                      best.gbrs.sample = colnames(geno.cor)[muga.max.cor.col],
                      cross.cor    = geno.cor[matrix(c(1:nrow(geno.cor), muga.max.cor.col),
                                    ncol = 2)],
                      mismatch = rep(NA, ncol(muga)),
                      stringsAsFactors = F)

  # If the MUGA sample was in the GBRS data, then add the correlation of the
  # MUGA sample to the GBRS sample of the same name.
  m = match(rownames(geno.cor), colnames(geno.cor))
  muga.row = which(!is.na(m))
  result$self.cor[muga.row] = geno.cor[matrix(c(muga.row, m[muga.row]), ncol = 2)]

  wh = which(!is.na(result$self.cor))
  result$mismatch[wh] = result$muga.sample[wh] != result$best.gbrs.sample[wh]

  return(result)

} # sample_mismatch()



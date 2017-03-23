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

  # Use the Pearson correlation of the haplotype probabilities as a vector
  # to compare samples.
  message("Calculating corrlations...")
  result = data.frame(muga.sample = colnames(muga),
                      self.cor    = rep(NA, ncol(muga)),
                      best.gbrs.sample = rep(NA, ncol(muga)),
                      best.gbrs.cor    = rep(NA, ncol(muga)),
                      mismatch = rep(NA, ncol(muga)),
                      stringsAsFactors = F)
  m = intersect(colnames(muga), colnames(gbrs))
  muga.m = match(m, colnames(muga))
  gbrs.m = match(m, colnames(gbrs))
  result.m = match(m, result$muga.sample)
  for(i in 1:length(m)) {

    result$self.cor[result.m[i]] = cor(muga[,muga.m[i]], gbrs[,gbrs.m[i]])

  } # for(i)

  ok = which(result$self.cor >= 0.6)
  result$best.gbrs.sample[ok] = result$muga.sample[ok]
  result$best.gbrs.cor[ok]    = result$self.cor[ok]

  # For any sample with self.cor < 0.6, look for a sample with a better
  # match.
  mismatch = which(result$self.cor < 0.6)
  for(i in 1:length(mismatch)) {

    tmp.cor  = cor(muga[,mismatch[i]], gbrs)
    best.cor = tmp.cor[1,which.max(tmp.cor[1,])]
    result$best.gbrs.sample[mismatch[i]] = names(best.cor)
    result$best.gbrs.cor[mismatch[i]]    = best.cor

  } # for(i)

  # If the best GBRS correlation is still low, cross-check against all
  # MUGA samples.
  mismatch = which(result$best.gbrs.cor < 0.6)
  for(i in 1:length(mismatch)) {

    gbrs.column = which(colnames(gbrs) == result$muga.sample[mismatch[i]])
    tmp.cor  = cor(muga, gbrs[,gbrs.column])
    best.cor = tmp.cor[which.max(tmp.cor[,1]),1]
    result$best.gbrs.sample[mismatch[i]] = names(best.cor)
    result$best.gbrs.cor[mismatch[i]]    = best.cor

  } # for(i)  

  # Set mismatch column to TRUE if the self-correlation column is 
  # less than 0.65.
  result$mismatch = result$self.cor < 0.65

  return(result)

} # sample_mismatch()



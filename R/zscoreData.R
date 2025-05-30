#' Z-transform data
#'
#' The z-transform is meant to work with normalized,
#' imputed data
#' @param data - Normalized, imputed data (not logged). Data matrix includes
#'               features as rows, samples as columns.
#' @param ref - Normalized, imputed reference sample data (not logged). Data
#'              includes features as rows, samples as columns.
#' @return zscored.data - Z-transformed data.
#' @importFrom stats quantile qqnorm lm
#' @examples
#' dis_data = matrix(rexp(500), ncol=100)
#' rownames(dis_data)=sprintf("Feature%d",seq_len(nrow(dis_data)))
#' colnames(dis_data)=sprintf("Sample%d",seq_len(ncol(dis_data)))
#' ref_data = matrix(rexp(500), ncol=100)
#' rownames(ref_data)=sprintf("Feature%d",seq_len(nrow(ref_data)))
#' colnames(ref_data)=sprintf("Sample%d",seq_len(ncol(ref_data)))
#' zscored.data=data.zscoreData(dis_data,ref_data)
#' @export zscoreData
zscoreData = function(data, ref) {
  print("zscoreData() called.")

  # Only metabolites that also occur in the reference population can be z-scored
  data = data[which(rownames(data) %in% rownames(ref)),]
  ref = ref[which(rownames(ref) %in% rownames(data)),]
  data = as.matrix(data[sort(rownames(data)),])
  ref = as.matrix(ref[sort(rownames(ref)),])

  # Log transform data
  # data = log(data)
  # ref = log(data.matrix(ref))

  zscore.data = matrix(NA, nrow=nrow(data), ncol=ncol(data))
  rownames(zscore.data) = rownames(data)
  colnames(zscore.data) = colnames(data)
  for (met in 1:nrow(data)) {
    met_data = as.numeric(ref[met,])
    rmSamples = unique(c(which(is.na(met_data)), which(is.infinite(met_data))))
    if (length(rmSamples)>0) {
      x = met_data[-rmSamples]
    } else {
      x = met_data
    }
    if (all(is.na(x))) {

    } else {
      #x = x[intersect(which(x>quantile(x, 0.025)), which(x<quantile(x, .975)))]
      d = qqnorm(x, plot.it = FALSE);
      x = as.numeric(d$y)
      z = as.numeric(d$x)
      df = data.frame(x=x,z=z)
      t = lm(x~z, data=df)
      mn.est = as.numeric(t$coefficients[1])
      sd.est = as.numeric(t$coefficients[2])
      rm(d,x,z,df,t)
      zscore.data[met,] = (data[met, ]-mn.est)/sd.est
    }
  }

  return(zscore.data)
}

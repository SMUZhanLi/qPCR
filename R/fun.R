#' Single threshold of real-time PCR
#' @param x a dataframe or matrix with two or more columns.
#' The name of first columns has to be "Cycle",and the names from the
#' second column to the last column are "Well".
#' @examples
#' ## threshold(x)
#' @export

threshold <- function(x) {
  llp_start <- function(rfu) {
    bl_end <- function(rfu) {
      len <- length(rfu)
      for (i in seq(3, len - 3)) {
        end <- i
        bl <- rfu[2:end]
        t <- 10 * sd(bl)
        if (rfu[i + 2] > t & all(max(bl) < rfu[(i + 2):len])) break
      }
      return(end)
    }

    end <- bl_end(rfu)
    len <- length(rfu)
    t <- NA
    for (i in seq(end, len - 3)) {
      x <- i:(i + 3)
      y <- rfu[i:(i + 3)]
      r <- cor(x, y)
      if (r^2 > 0.99) {
        t <- rfu[i]
        break
      }
    }
    return(t)
  }

  multi_t <- apply(x[-1], 2, llp_start)
  single_t <- quantile(multi_t, 0.8, na.rm = T)
  return(single_t)
}

#' The Ct value of real-time PCR
#' @param x a dataframe or matrix with two or more columns.
#' The name of first columns has to be "Cycle",and the names from the
#' second column to the last column are "Well".
#' @param s_t The single threshold of qPCR.
#' @examples
#' ## cycle_threshold(x, s_t)
#' @export

cycle_threshold <- function(x, s_t) {

  get_ct <- function(rfu, t) {
    len <- length(rfu)
    ct <- NA
    for (i in seq(1, len - 1)) {
      r1 <- rfu[i]
      r2 <- rfu[i + 1]
      if (r1 < t & r2 > t) {
        ct <- round(i + (t - r1) / (r2 - r1), 3)
        break
      }
    }
    return(ct)
  }

  r <- apply(x[-1], 2, get_ct, t = s_t)
  df <- data.frame(Well = names(r), Ct = r)
  rownames(df) <- NULL
  return(df)
}

#' Comparative Ct method
#' @param x a dataframe or matrix.It must have the following columns:
#' Target, Sample, Ct
#' @param goi gene of interest
#' @param ref reference gene
#' @param ctrl control sample of relative quantification in qPCR
#' @examples
#' ## ddct(x, goi, ref, ctrl)
#' @export

ddct <- function(x, goi, ref, ctrl) {
  cct_goi <- x[x$Sample == ctrl & x$Target == goi, "Ct"]
  cct_ref <- x[x$Sample == ctrl & x$Target == ref, "Ct"]
  dct_ctrl <- cct_goi - cct_ref
  for (i in unique(x$Sample)) {
    ct_goi <- x[x$Sample == i & x$Target == goi, "Ct"]
    ct_ref <- x[x$Sample == i & x$Target == ref, "Ct"]
    dct <- ct_goi - ct_ref
    x[x$Sample == i & x$Target == goi, "Rq"] <- round(2^(dct_ctrl - dct), 3)
  }
  return(x)
}

#' Positive calling
#' @param x a dataframe or matrix with two or more columns.
#' The name of first columns has to be "Cycle",and the names from the
#' second column to the last column are "Well".
#' @param n End cycles of every well
#' @param negative one or more wells which are negative control sample
#' @param method The calculation mode of tolerance.
#' @param p Percent of range or absolute RFU
#' @details
#' If method is POR(the default), the range of p is from 0.1(the default) to 0.9.
#' If method is RFU, p is a absolute RFU.
#' @examples
#' ## neg <- c("w1", "w2")
#' ## positive(x, 5, neg)
#' ## method is RFU
#' ## positive(x, 5, neg, "RFU", 2500)
#' @export

positive <- function(x, n, negative, method = "POR", p = 0.1) {
  AER <- apply(x, 2, function(rfu) {mean(tail(rfu,n))})
  AER_Neg <- mean(AER[negative])
  if (method == "POR") {
    tolerance <- (max(AER) - AER_Neg) * p
  }
  else if (method == "RFU") {
    tolerance <- p
  }
  cutoff <- AER_Neg + tolerance
  call <- ifelse(AER > cutoff, "Positive", NA)
  r <- data.frame(Well = names(AER), Call = call)
  rownames(r) <- NULL
  return(r)
}

#' Amplification efficiency
#' @param x a dataframe or matrix.It must have the following columns:
#' Target, Sample, Ct, Log.Starting.Quantity
#' @import dplyr
#' @export

efficiency <- function(x){
  r <- dplyr::summarise(
        dplyr::group_by(x, Target, Sample),
        slope = lm(Ct~Log.Starting.Quantity)$coefficients[[2]],
        `R^2` = cor(Log.Starting.Quantity, Ct)^2,
        `eff(%)` = round((10^(-1/slope) - 1)*100, 2)
      )
  return(r)
}




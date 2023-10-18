
#' Title original simu_zi function
#' @description  original simu_zi function, generating simualtion data from ZINB or ZIP distribution.
#' @param family ZINB or ZIP.
#' @param subject.n  sample size.
#' @param zi.p zero-inflated proportion.
#' @param mu if ZINB, it is the mu parameter, if ZIP, it is the lambda paramter.
#' @param size if ZINB, it is the size parameter, if ZIP, it is the unused argument.
#' @return  a vector. The simulation data.
#' @export
simu_zi=function (family, subject.n, zi.p = 0.5, mu = 0.5, size = 0.25) {
  Y = rep(NA, subject.n)
  ind.mix <- rbinom(length(Y), 1, zi.p)
  if (family == "ZIP") {
    Y[which(ind.mix != 0)] <- 0
    Y[which(ind.mix == 0)] = rpois(n = sum(ind.mix == 0),
                                   lambda = mu)
  }
  if (family == "ZINB") {
    Y[which(ind.mix != 0)] <- 0
    Y[which(ind.mix == 0)] = rnbinom(n = sum(ind.mix == 0),
                                     size = size, mu = mu)
  }
  return(Y)
}




#' Title  对数正态的的norm-----
#' @description  2023_1004 new simu function, generating simulation data from normal Normal distribution.
#' 2^(x) follow logNormal. Skewness，with cut(the threshold).
#' @param subject.n unused argument
#' @param zi.p unused argument
#' @param mu mu of Normal.
#' @param sd sd of Normal.
#' @param cut the threshold.Ensuring that the expression of SVG has obviously high values.
#' @return a vector. The simulation data.
#' @export
simu_part_norm=function (subject.n, zi.p = 0.5, mu = 0.5, sd = 0.25,cut=5) {
  Y = rep(NA, subject.n)
  ind.mix <- rbinom(length(Y), 1, zi.p)
  Y[which(ind.mix != 0)] <- 0

  random_numbers <- numeric(0)  # 创建一个空的向量来存储随机数
  while (length(random_numbers) <sum(ind.mix ==  0)) {   ##还不够时
    # 生成一个正态分布的随机数
    random_value <- rnorm(sum(ind.mix ==  0), mean =mu, sd = sd)
    random_value=random_value[which(random_value>=cut)]
    random_numbers <- c(random_numbers, random_value)
  }
  # Y[which(ind.mix == 0)] = rnorm(n = sum(ind.mix ==  0), mean = mu, sd = sd)
  Y[which(ind.mix == 0)] =  round(random_numbers[sample(length(random_numbers),sum(ind.mix==0))])
  return(Y)
}

#' Title zero-inflated Normal distribution, without cut.
#' @description 2023_1004 new simu function, generating simulation data from zero-inflated Normal distribution,
#' without the cut (threshold).
#' @param subject.n unused argument
#' @param zi.p unused argument
#' @param mu mu of Normal.
#' @param sd sd of Normal.
#' @return a vector. The simulation data.
#' @export
simu_norm2=function (subject.n, zi.p = 0.5, mu = 0.5, sd = 0.25)
{
  Y = rep(NA, subject.n)
  ind.mix <- rbinom(length(Y), 1, zi.p)
  Y[which(ind.mix != 0)] <- 0
  Y[which(ind.mix == 0)] = rnorm(n = sum(ind.mix ==  0), mean = mu, sd = sd)
  return(Y)
}









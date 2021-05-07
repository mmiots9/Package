#' @name par.2aov
#' @title Parametric/non parametric test evaluation for multivariate data
#' @description Get info about normality and variances of the given vectors in ordere to decide wether use 2-way ANOVA or Kruskal-Wallis
#' @usage mrt.par.2aov(x, y, z, type.of.int = "+")
#' @param x Numeric vector (or dataframe column)
#' @param y A vector of factors
#' @param z A vector of factors
#' @param tyoe.of.int "+" to compute the analysis on the two categorial variables only or "*" to compute also the analysis on the interaction
#' @returns
#' \item{Shapiro categorie}{A dataframe containing W and p-value coming from shapiro.test of x in the two different categoric conditions}
#' \item{Bartlett}{A dataframe containing k-squared value, df and p-value coming from bartlett.test of x in the two different categoric conditions and, if type.of.int = "*", also from the interaction of the two categorical}
#' \item{Shapiro interazione}{optional, only if type.of.int = "*". Return the W and p-value coming from the shapiro test of all the possible interaction of y and z over x}
#'
#' @examples
#' par.2aov(rnorm(100), sample(rep(c(0,1), 50)), sample(rep(c("a", "b", "c", "d"), 25)), type.of.int = "*")
#' @author Matteo Miotto
#'

#' @export




par.2aov <- function(x, y, z, type.of.int = NA){

  # 1. controlli
  # x, y, z numeriche
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }

  # valid type of int
  if (is.na(type.of.int)){
    stop("type.of.int not specified with no default")}
  if (type.of.int != "+" && type.of.int != "*") {
    stop("type.of.int must be '+' or '*'" )
  }

  # 2. shapiro categorie
  shapy <- tapply(x, y, shapiro.test)
  shapz <- tapply(x, z, shapiro.test)

  # 3. Inizializzo vettori utili
  categy <- c()
  fatty  <- c()
  wy     <- c()
  pvaly  <- c()
  categz <- c()
  fattz  <- c()
  wz     <- c()
  pvalz  <- c()

  # 4. Creo nomi categorie ingresso
  caty <- deparse1(substitute(y))
  catz <- deparse1(substitute(z))

  # 5. creo loop per estrarre cose importanti da shapy e shapz
  # shapy
  for (i in 1:length(shapy)){
    categy[i] <- caty
    fatty[i]  <- names(shapy[i])
    wy[i]     <- shapy[[i]]$statistic
    pvaly[i]  <- shapy[[i]]$p.value
  }

  # shapz
  for (i in 1:length(shapz)){
    categz[i] <- catz
    fattz[i]  <- names(shapz[i])
    wz[i]     <- shapz[[i]]$statistic
    pvalz[i]  <- shapz[[i]]$p.value
  }

  # 6. creo dataframe
  shcat <- data.frame(categoric = c(categy, categz), fattore = c(fatty, fattz), W = c(wy, wz), pval = c(pvaly, pvalz))

  #7. Bartlett
  # categorie
  bary <- bartlett.test(x, y)
  barz <- bartlett.test(x, z)

  # estraggo valori
  barcat    <- caty
  ksq       <- bary$statistic
  bardf     <- bary$parameter
  barpv     <- bary$p.value
  barcat[2] <- catz
  ksq[2]    <- barz$statistic
  bardf[2]  <- barz$parameter
  barpv[2]  <- barz$p.value

  # se interazione
  if(type.of.int == "*"){
    bartint <- bartlett.test(x, interaction(y,z))
    barcat[3] <- "interaction"
    ksq[3]    <- bartint$statistic
    bardf[3]  <- bartint$parameter
    barpv[3]  <- bartint$p.value
  } else {
    barcat[3] <- NA
    ksq[3]    <- NA
    bardf[3]  <- NA
    barpv[3]  <- NA
  }

  # creo dataframe
  bartot <- data.frame(barcat, ksq, bardf, barpv)
  names(bartot) <- c("categoric", "k-squared", "df", "pval")
  if(type.of.int != "*") {
  bartot <- bartot[-3,]
  }

  # 8. Shapiro se *
  if(type.of.int == "*"){

    # creo lista con tutti i risultati
    mat <- tapply(x, list(y, z), shapiro.test)
    nm1 <- outer(rownames(mat), colnames(mat), FUN=paste)
    nm2 <- setNames(c(mat), nm1)
    nom <- strsplit(names(nm2), split = " ") # creo lista con fattori categorie

    # creo vettori interesse
    prima   <- c()
    seconda <- c()
    wint    <- c()
    pvalint <- c()

    # loop per riempire vettori
    for(i in 1:length(nm2)){
      prima[i]   <- nom[[i]][1]
      seconda[i] <- nom[[i]][3]
      wint[i]    <- nm2[[i]]$statistic
      pvalint[i] <- nm2[[i]]$p.value
    }

    # creo dataframe
    shint <- data.frame(prima, seconda, wint, pvalint)
    names(shint) <- c(caty, catz, "W", "pval")

    # return
    res <- list(shcat, shint, bartot)
    names(res) <- c("Shapiro categorie", "Shapiro interazione", "Bartlett")
  }
  else {
    res <- c(list(shcat, bartot))
    names(res) <- c("Shapiro categorie", "Bartlett")
  }
  return(res)
}

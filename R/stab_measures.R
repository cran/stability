#' @name    stab_measures
#' @aliases stab_measures
#' @title  Stability Measures for Genotypes by Environment Interaction (GEI)
#' @description Stability Measures for Genotypes by Environment Interaction (GEI)
#' @param data data.frame
#' @param G Genotypes Factor
#' @param E Environment Factor
#' @param Y Response Variable
#' @return Stability Measures
#' @author
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#'
#' @importFrom matrixStats rowVars rowSds
#' @importFrom reshape2 acast
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.StabMeasures <- stab_measures(data = ge_data, G = Gen, E = Env, Y = Yield)
#' Yield.StabMeasures
#'
stab_measures <-
  function(data, G, E, Y){
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))

    GE.Means <-
      reshape2::acast(
        data          = data
        , formula       = as.formula(paste(G, E, sep = "~"))
        , fun.aggregate = mean
        , margins       = FALSE
        , value.var     = Y
      )

    GGE.Effects <-
      sweep(
        x = GE.Means
        , MARGIN = 2
        , STATS = colMeans(GE.Means)
      )

    GE.Effects <-
      sweep(
        x = GGE.Effects
        , MARGIN = 1
        , STATS = colMeans(GGE.Effects)
      )

    Gen.Means   <- rowMeans(GE.Means)
    Gen.SS      <- (ncol(GE.Means)-1)*matrixStats::rowVars(GE.Means)
    Gen.Var     <- matrixStats::rowVars(GE.Means)
    Gen.CV      <- matrixStats::rowSds(GE.Means)/rowMeans(GE.Means)*100
    Ecovalence  <- matrixStats::rowVars(GE.Effects)*(ncol(GE.Effects)-1)
    GE.SS       <- sum(Ecovalence)
    GE.MSE      <- sum(Ecovalence)/((nrow(GE.Means)-1)*(ncol(GE.Means)-1))
    Shukla.Var  <- nrow(GE.Means)/((nrow(GE.Means)-2)*(ncol(GE.Means)-1)) * Ecovalence - GE.MSE/(nrow(GE.Means)-2)


    StabMeasures <-
      data.frame(
          Mean        = Gen.Means
        , GenSS       = Gen.SS
        , Var         = Gen.Var
        , CV          = Gen.CV
        , Ecov        = Ecovalence
        , "ShuklaVar" = Shukla.Var
        )
    return(StabMeasures)
  }

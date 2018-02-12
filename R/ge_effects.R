#' @name    ge_effects
#' @aliases ge_effects
#' @title Genotype by Environment Interaction Effects
#' @description Provides Genotype by Environment Interaction Effects
#' @param data data.frame
#' @param G Genotypes Factor
#' @param E Environment Factor
#' @param Y Response Variable
#' @return Effects
#'
#' @author
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.Effects <- ge_effects(data = ge_data, G = Gen, E = Env, Y = Yield)
#' names(Yield.Effects)
#'
#' Yield.Effects$ge_means
#' Yield.Effects$ge_effects
#' Yield.Effects$gge_effects


ge_effects <-
  function(data, G, E, Y){
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))

    ge_means <-  tapply(data$Y, list(data$G, data$E), mean)

    gge_effects <-
      sweep(
          x      = ge_means
        , MARGIN = 2
        , STATS  = colMeans(ge_means)
      )

    ge_effects <-
      sweep(
          x      = gge_effects
        , MARGIN = 1
        , STATS  = rowMeans(gge_effects)
      )

    return(list(
        "ge_means"    = ge_means
      , "gge_effects" = gge_effects
      , "ge_effects"  = ge_effects
    ))
  }

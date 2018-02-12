#' @name    ge_means
#' @aliases ge_means
#' @title Genotype by Environment Interaction Means and Ranks
#' @description Provides Genotype by Environment Interaction Means along with their Ranks
#' @param data data.frame
#' @param G Genotypes Factor
#' @param E Environment Factor
#' @param Y Response Variable
#' @return Means and Ranks
#'
#' @author
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom reshape2 acast
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#'
#' Yield.ge_means <- ge_means(data = ge_data, G = Gen, E = Env, Y = Yield)
#'
#' Yield.ge_means$ge_means1
#' Yield.ge_means$ge_means
#' Yield.ge_means$ge_ranks
#'
ge_means <-
  function(data, G, E, Y){
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))

    ge_means1 <-
      reshape2::acast(
          data          = data
        , formula       = as.formula(paste(G, E, sep = "~"))
        , fun.aggregate = mean
        , margins       = TRUE
        , value.var     = Y
      )

    ge_means2 <-
      reshape2::acast(
        data          = data
        , formula       = as.formula(paste(G, E, sep = "~"))
        , fun.aggregate = mean
        , margins       = FALSE
        , value.var     = Y
      )

    g_means <-
      data %>%
      dplyr::group_by(data$G) %>%
      dplyr::summarise(mean(data$Y, na.rm = TRUE))
    names(g_means) <- c("G", "Y")

    e_means <-
      data %>%
      dplyr::group_by(data$E) %>%
      dplyr::summarise(mean(data$Y, na.rm = TRUE))
    names(e_means) <- c("E", "Y")

    Eg_means  <- t(ge_means2)
    ge_ranks  <- matrix(data = NA, nrow = nrow(Eg_means), ncol = ncol(Eg_means))
    for(i in 1:nrow(Eg_means)){
      ge_ranks[i, ] <- names(sort(Eg_means[i, ], decreasing = TRUE))
      dimnames(ge_ranks) <- list(rownames(Eg_means), c(1:ncol(Eg_means)))
    }

    return(
      list(
          "ge_means1" = ge_means1
        , "ge_means"  = ge_means2
        , "g_means"   = g_means
        , "e_means"   = e_means
        , "ge_ranks"  = ge_ranks
      ))
  }

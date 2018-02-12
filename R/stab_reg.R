#' @name    stab_reg
#' @aliases stab_reg
#' @title Individual Regression for each Genotype
#' @description Individual Regression for each Genotype in Genotypes by Environment Interaction (GEI)
#' @param data data.frame
#' @param Rep Replication Factor
#' @param G Genotypes Factor
#' @param E Environment Factor
#' @param Y Response Variable
#' @return Additive ANOVA
#' @author
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#' @importFrom lme4 lmList
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' Yield.StabReg <- stab_reg(data = ge_data, G = Gen, E = Env, Y = Yield)
#' Yield.StabReg
#'
stab_reg <-
  function(data, Rep, G, E, Y){
    Rep <- deparse(substitute(Rep))
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))

    EnvMean  <- ave(data$Y, data$E)
    GEMean   <- ave(data$Y, data$G:data$E)
    DataNew  <- data.frame(GEMean, EnvMean, Gen = data$Gen)
    g <- length(levels(data$G))
    e <- length(levels(data$E))
    r <- length(levels(data$E))

    IndvReg <- lme4::lmList(GEMean ~ EnvMean|Gen, data = DataNew)
    IndvRegFit <- summary(IndvReg)

    StabIndvReg <-
      data.frame(
        "Slope"  = coef(IndvRegFit)[ , , 2][ ,1]
        , "LCI"   = confint(IndvReg)[ , ,2][ ,1]
        , "UCI"   = confint(IndvReg)[ , ,2][ ,2]
        , "R.Sqr" = IndvRegFit$r.squared
        , "RMSE"  = IndvRegFit$sigma
        , "SSE"   = IndvRegFit$sigma^2*IndvRegFit$df[,2]
        , "Delta" = IndvRegFit$sigma^2*IndvRegFit$df[,2]/r
      )

    return(StabIndvReg)
  }


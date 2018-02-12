#' @name    add_anova
#' @aliases add_anova
#' @title Additive ANOVA for Genotypes by Environment Interaction (GEI) model
#' @description Additive ANOVA for Genotypes by Environment Interaction (GEI) model
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
#' @importFrom stats anova as.formula ave coef confint lm pf terms
#'
#' @export
#'
#' @examples
#' data(ge_data)
#' YieldANOVA <-
#'      add_anova(
#'         data = ge_data
#'       , Rep  = Rep
#'       , G    = Gen
#'       , E    = Env
#'       , Y    = Yield
#'       )[[1]]
#' YieldANOVA
#'
add_anova <-
  function(data, Rep, G, E, Y){
    Rep <- deparse(substitute(Rep))
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))
    fm1 <- lm(formula=terms(data$Y ~ data$E + data$Rep:data$E + data$G + data$G:data$E, keep.order = TRUE), data = data)
    fm1ANOVA <- anova(fm1)
    rownames(fm1ANOVA) <- c("Env", "Rep(Env)", "Gen", "Gen:Env", "Residuals")
    fm1ANOVA[1, 4] <- fm1ANOVA[1, 3]/fm1ANOVA[2, 3]
    fm1ANOVA[2, 4] <- NA
    fm1ANOVA[1, 5] <- 1-pf(as.numeric(fm1ANOVA[1, 4]), fm1ANOVA[1, 1], fm1ANOVA[2, 1])
    fm1ANOVA[2, 5] <- 1-pf(as.numeric(fm1ANOVA[2, 4]), fm1ANOVA[2, 1], fm1ANOVA[5, 1])
    class(fm1ANOVA) <- c("anova", "data.frame")

    return(
      list(
        "fm1ANOVA" = fm1ANOVA
      ))
  }

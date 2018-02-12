#' @name    Eberhart_Russel_anova
#' @aliases Eberhart_Russel_anova
#' @title  ANOVA of Eberhart & Russel’s Model
#' @description ANOVA of Eberhart & Russel’s Model
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
#' @export
#'
#' @examples
#' data(ge_data)
#' Yield.ERANOVA <-
#'      Eberhart_Russel_anova(
#'           data = ge_data
#'        , G     = Gen
#'        , E     = Env
#'        , Y     = Yield
#'        )[[1]]
#' Yield.ERANOVA
#'
Eberhart_Russel_anova <-
  function(data, Rep, G, E, Y){
    Rep <- deparse(substitute(Rep))
    G   <- deparse(substitute(G))
    E   <- deparse(substitute(E))
    Y   <- deparse(substitute(Y))

    g <- length(levels(data$Gen))
    e <- length(levels(data$Env))
    r <- length(levels(data$Rep))

    EnvMean  <- ave(data$Y, data$E)
    GEMean   <- ave(data$Y, data$G:data$E)
    DataNew  <- data.frame(GEMean, EnvMean, Gen = data$Gen)

    IndvReg <- lme4::lmList(GEMean ~ EnvMean|Gen, data = DataNew)
    IndvRegFit <- summary(IndvReg)

    StabIndReg <-
      data.frame(
        "Slope"  = coef(IndvRegFit)[ , , 2][ ,1]
        , "LCI"   = confint(IndvReg)[ , ,2][ ,1]
        , "UCI"   = confint(IndvReg)[ , ,2][ ,2]
        , "R.Sqr" = IndvRegFit$r.squared
        , "RMSE"  = IndvRegFit$sigma
        , "SSE"   = IndvRegFit$sigma^2*IndvRegFit$df[,2]
        , "Delta" = IndvRegFit$sigma^2*IndvRegFit$df[,2]/r
      )


    fm1 <- lm(formula=terms(data$Y ~ data$E + data$Rep:data$E + data$G + data$G:data$E, keep.order = TRUE), data = data)
    fm1ANOVA <- anova(fm1)
    rownames(fm1ANOVA) <- c("Env", "Rep(Env)", "Gen", "Gen:Env", "Residuals")
    fm1ANOVA[1, 4] <- fm1ANOVA[1, 3]/fm1ANOVA[2, 3]
    fm1ANOVA[2, 4] <- fm1ANOVA[2, 3]/fm1ANOVA[5, 3]
    fm1ANOVA[1, 5] <- 1-pf(as.numeric(fm1ANOVA[1, 4]), fm1ANOVA[1, 1], fm1ANOVA[2, 1])
    fm1ANOVA[2, 5] <- 1-pf(as.numeric(fm1ANOVA[2, 4]), fm1ANOVA[2, 1], fm1ANOVA[5, 1])
    class(fm1ANOVA) <- c("anova", "data.frame")

    Df <-
      c(
        g*e-1
        , g-1
        , g*(e-1)
        , 1
        , g-1
        , g*(e-2)
        , rep(e-2, g)
        , e*g*(r- 1)
      )

    PooledError <- fm1ANOVA["Residuals", "Sum Sq"]/r
    TotalSS <- (fm1ANOVA["Env", "Sum Sq"] + fm1ANOVA["Rep(Env)", "Sum Sq"] + fm1ANOVA["Gen:Env", "Sum Sq"])/r
    GenSS   <- fm1ANOVA["Rep(Env)", "Sum Sq"]/r
    EnvGESS <- (fm1ANOVA["Env", "Sum Sq"] + fm1ANOVA["Gen:Env", "Sum Sq"])/r
    EnvL    <- fm1ANOVA[1, "Sum Sq"]/r
    GELSS   <- EnvGESS - EnvL - sum(StabIndReg$Delta)
    PooledDev  <- sum(StabIndReg$Delta)
    GensSS  <- StabIndReg$Delta


    SS <-
      c(
        TotalSS
        , GenSS
        , EnvGESS
        , EnvL
        , GELSS
        , PooledDev
        , GensSS
        , PooledError
      )

    MS <- SS/Df

    F <-
      c(
        NA
        , MS[2]/MS[6]
        , NA
        , NA
        , MS[5]/MS[6]
        , NA
        , MS[7:(length(MS) - 1)]/MS[length(MS)]
        , NA
      )

    PLines <- 1 - pf(F[7:(length(MS) - 1)], Df[7], Df[length(Df)])
    pval <- c(NA, 1 - pf(F[2], Df[2], Df[6]), NA, NA, 1 - pf(F[5], Df[5], Df[6]), NA, PLines, NA)

    ANOVA <- data.frame(Df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F, `Pr(>F)` = pval, check.names = FALSE)

    rownames(ANOVA) <-
      c(
        "Total"
        , "Gen"
        , "Env + (Gen x Env)"
        , "  Env (linear)"
        , "  Gen x Env(linear)"
        , "  Pooled deviation"
        , paste0("    ", levels(data$Gen))
        , "Pooled error"
      )

    class(ANOVA) <- c("anova", "data.frame")

    return(list(
      "ANOVA" = ANOVA
    )
    )
  }


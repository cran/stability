#' @name    stab_par
#' @aliases stab_par
#' @title  Stability Parameters for Genotypes by Environment Interaction (GEI)
#' @description Stability Parameters for Genotypes by Environment Interaction (GEI)
#' @param x PARAM_DESCRIPTION
#' @param rep No of Replicates
#' @param MSE Mean Square Error
#' @param alpha Level of Significance, default is 0.1
#' @param Env.Cov Environmental Covariate, default is NULL
#' @return Stability Parameters
#' @author
#' Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'
#' @references
#'  Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'
#'
#' @importFrom stats anova as.formula ave coef confint lm pf terms qf qt
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#'
#' data(ge_data)
#' YieldANOVA <- add_anova(ge_data, Rep, Gen, Env, Yield)[[1]]
#' Yield.ge_means <- ge_means(data = ge_data, G = Gen, E = Env, Y = Yield)
#' Yield.ge_means$ge_means
#'
#' Yield.StabPar <-
#' sta_par(
#'   x     = as.matrix(Yield.ge_means$ge_means)
#' , rep   = length(levels(ge_data$Rep))
#' , MSE   = YieldANOVA["Residuals", "Mean Sq"]
#' , alpha = 0.1
#' )
#'
#' Yield.StabPar
#'
sta_par <-
  function(x, rep, MSE, alpha = 0.1, Env.Cov = NULL) {
    x <- as.matrix(x)
    if(is.null(Env.Cov)==TRUE)  Env.Cov <- colMeans(x) - mean(x) else Env.Cov <- Env.Cov

    x1 <- as.data.frame(x)
    DataNew <-
      data.frame(
        "Gen" = rownames(x1)
        , reshape2::melt(data=x1, measure.vars=dput(names(x1)), variable.name = "Env", value.name = "Y"))

    fm1ANOVA <- anova(lm(formula=Y~Gen*Env, data = DataNew))

    g <- nrow(x)
    e <- ncol(x)
    G.Effects <-
      sweep(
        x = x
        , MARGIN = 1
        , STATS = rowMeans(x)
      )

    GE.Effects <-
      sweep(
        x = G.Effects
        , MARGIN = 2
        , STATS = colMeans(G.Effects)
      )


    Gen.Mean <- rowMeans(x)
    GenSS    <- rowVars(x)*(e - 1)
    Gen.Var  <- rowVars(x)
    Gen.CV   <- rowSds(x)/rowMeans(x)*100
    Ecov     <- rowVars(GE.Effects)*(e - 1)
    GE.SS    <- sum(Ecov)
    GE.MSE   <- sum(Ecov)/((g-1)*(e-1))

    Sigma   <- g/((g-2)*(e-1)) * Ecov - GE.MSE/(g-2)
    Beta    <- GE.Effects %*% Env.Cov/sum(Env.Cov^2)
    New.GE.Effects <- GE.Effects - matrix(data = Env.Cov, nrow=g, ncol=e, byrow=TRUE) * matrix(data = Beta, nrow = g, ncol = e, byrow=FALSE)
    SP      <- (g/((g-2)*(e-2)))*(rowSums(New.GE.Effects^2)- (sum(rowSums(New.GE.Effects^2))/(g*(g-1))))
    SSRes  <- ((g - 1)*(e - 2)*sum(SP))/g
    SSHetro   <- fm1ANOVA["Gen:Env", "Sum Sq"] - SSRes


    Df <- c(sum(fm1ANOVA[-4, "Df"]), fm1ANOVA[-4, "Df"], fm1ANOVA["Gen", "Df"], fm1ANOVA["Gen:Env", "Df"] - fm1ANOVA["Gen", "Df"], fm1ANOVA["Gen", "Df"]^2*(fm1ANOVA["Env", "Df"]+1))
    SumS  <- c(sum(fm1ANOVA[-4, "Sum Sq"]), fm1ANOVA[-4, "Sum Sq"], SSHetro, SSRes, NA)*rep
    MeanSS <- c(NA, SumS[2:6]/Df[2:6], MSE)

    F <-
      c(
        NA
        , MeanSS[2]/MeanSS[6]
        , MeanSS[3]/MeanSS[7]
        , MeanSS[4]/MeanSS[7]
        , MeanSS[5]/MeanSS[6]
        , MeanSS[6]/MeanSS[7]
        , NA
      )

    pval <-
      c(
        NA
        , 1 - pf(F[2], Df[2], Df[6])
        , 1 - pf(F[3], Df[3], Df[7])
        , 1 - pf(F[4], Df[4], Df[7])
        , 1 - pf(F[5], Df[5], Df[6])
        , 1 - pf(F[6], Df[6], Df[7])
        , NA
      )

    ANOVA <- data.frame(Df, `Sum Sq` = SumS, `Mean Sq` = MeanSS,
                        `F value` = F, `Pr(>F)` = pval, check.names = FALSE)

    rownames(ANOVA) <-
      c(
        "Total"
        , "Gen"
        , "Env"
        , "Gen x Env"
        , "  Heterogeneity"
        , "  Residual"
        , "Pooled error"
      )

    class(ANOVA) <- c("anova", "data.frame")



    dfShukla <- e*(g - 1)*(rep - 1) # df for Shukla Test
    FM0 <- qf(1 - alpha, e-1, dfShukla)
    F05 <- qf(0.95, e-1, dfShukla)
    F01 <- qf(0.99, e-1, dfShukla)
    MS    <- SSRes/(fm1ANOVA["Gen", "Df"]^2*(fm1ANOVA["Env", "Df"]+1))
    DMV05 <- qt(0.95, dfShukla) * sqrt(2 * MSE/(rep * e))
    MES   <-  mean(x)
    FS    <- Sigma*rep/MSE
    FSS   <- SP*rep/MSE

    NN <-
      ifelse(
        test = FS < F05
        , yes  = "ns"
        , no   =  ifelse(
          test = FS < F01 & FS >= F05
          , yes  = "*"
          , no = "**"
        )
      )


    MMM <-
      ifelse(
        test = FSS < F05
        , yes  = "ns"
        , no   =  ifelse(
          test = FSS < F01 & FSS >= F05
          , yes  = "*"
          , no = "**"
        )
      )

    Stb.Results <-
      data.frame(
        Mean  = Gen.Mean
        , GenSS = GenSS/rep
        , Var   = Gen.Var/rep
        , CV    = Gen.CV/rep
        , Ecov  = Ecov/rep
        , GE.SS = GE.SS/rep
        , GE.MSE = GE.MSE/rep
        , "Sigma" = Sigma/rep
        , "."     = NN
        , "SP"    = SP/rep
        , ".."     = MMM
      )
    rownames(Stb.Results) <- rownames(x)

    # Simultaneous Selection
    F1 <- as.numeric(as.character(
      cut(
        x = FS
        , breaks = c(-Inf, FM0, F05, F01,  Inf)
        , labels = c(0, 2, -4, -8)
        , right=TRUE
      )
    ))


    MV1 <- as.numeric(as.character(
      cut(
        x = Gen.Mean
        , breaks = c(-Inf, MES-2*DMV05, MES-DMV05, MES, MES+DMV05, MES+2*DMV05, Inf)
        , labels = c(-3, -2, -1, 1, 2, 3)
        , right=FALSE
      )
    ))

    R    <- rank(x = Gen.Mean, na.last = TRUE, ties.method = "min")
    GY   <- R + MV1
    GYS  <- GY + F1
    SGYS <- sum(GYS)
    MGYS <- SGYS/g
    GYY  <- ifelse(test = GYS > MGYS, yes = "+", no =  "-")

    SimulSel <- data.frame(Gen.Mean, R, MV1, GY, Sigma/rep, F1, GYS, GYY)
    rownames(SimulSel) <- rownames(DataNew$Y)
    names(SimulSel) <- c("Yield", "Rank", "Adjustment", "Adj.Rank", "Sigma", "Stab.Rating", "YSi", "Select")

    return(list(
      ANOVA     = ANOVA
      , StabPar   = Stb.Results
      , SimultSel = SimulSel
    )
    )

  }

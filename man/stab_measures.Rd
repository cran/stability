% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stab_measures.R
\name{stab_measures}
\alias{stab_measures}
\alias{stab_measures.default}
\title{Stability Measures for Genotypes by Environment Interaction (GEI)}
\usage{
stab_measures(.data, .y, .gen, .env)

\method{stab_measures}{default}(.data, .y, .gen, .env)
}
\arguments{
\item{.data}{data.frame}

\item{.y}{Response Variable}

\item{.gen}{Genotypes Factor}

\item{.env}{Environment Factor}
}
\value{
Stability Measures
}
\description{
Stability Measures for Genotypes by Environment Interaction (GEI)
}
\examples{

data(ge_data)
Yield.StabMeasures <- stab_measures(
                .data  = ge_data
               , .y    = Yield
               , .gen  = Gen
               , .env  = Env
               )
Yield.StabMeasures



}
\references{
Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
             New Delhi: Kalyani.
}
\author{
\enumerate{
         \item Muhammad Yaseen (\email{myaseen208@gmail.com})
         \item Kent M. Edkridge (\email{keskridge1@unl.edu})
         }
}

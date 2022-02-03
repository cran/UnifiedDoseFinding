\name{select_mtd_gBOIN_continuous}
\alias{select_mtd_gBOIN_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD) for single agent trials using gBOIN design
}
\description{
Select the maximum tolerated dose (MTD) when the trial is completed using gBOIN design (Mu et al. 2017)
}
\usage{
select_mtd_gBOIN_continuous(target, npts, ntox)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the continuous target score
}
  \item{npts}{
the number of patients enrolled at each dose level
}
  \item{ntox}{
the toxicity score at each dose level
}
}

\value{
\code{select_mtd_gBOIN_continuous()} returns the selected dose
}
\references{
Rongji Mu, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, Jun Yin: gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points. Royal Statistical Society 2019
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
target <- 1.47
n <- c(3, 3, 3, 9, 0, 0)
y <- c(0.1951265, 1.5434317, 2.1967343, 13.9266838, 0, 0)
select_mtd_gBOIN_continuous(target = target, npts = n, ntox = y)
}

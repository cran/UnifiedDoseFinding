\name{select_mtd_gBOIN_TB}
\alias{select_mtd_gBOIN_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD) defined by Toxicity Burden (TB) Score for single agent trials using gBOIN design
}
\description{
Select the maximum tolerated dose (MTD) defined by the toxicity burden (BT) score proposed by Bekele et al. (2004) when the trial is completed using the generalized Bayesian optimal interval (gBOIN) design (Mu et al. 2017). The algorithm of this function is exactly same to the Select_mtd_gBOIN.Continuous() just the input parameter is used by the TB score
}
\usage{
select_mtd_gBOIN_TB(target, npts, ntox)
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
\code{select_mtd_gBOIN_TB()} returns the selected dose
}
\references{
B. Nebiyou Bekele & Peter F Thall (2004) Dose-Finding Based on Multiple Toxicities in a Soft Tissue Sarcoma Trial, Journal of the American Statistical Association

Mu, Rongji, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, and Jun Yin. "gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points." Journal of the Royal Statistical Society. Series C: Applied Statistics 68, no. 2 (2019): 289-308.

}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
target <- 3.344
n <- c(3, 9, 6, 0, 0, 0, 0, 0, 0, 0)
y <- c(5.5, 26.95, 25.3, 0, 0, 0, 0, 0, 0, 0)
select_mtd_gBOIN_TB(target = target, npts = n, ntox = y)
}

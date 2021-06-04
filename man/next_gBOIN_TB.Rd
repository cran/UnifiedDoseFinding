\name{next_gBOIN_TB}
\alias{next_gBOIN_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a maximum tolerated dose (MTD) defined by Toxicity Burden (TB) Score  using gBOIN design
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find the MTD defined by the toxicity burden (BT) score proposed by Bekele et al. (2004) using the generalized Bayesian optimal interval (gBOIN) design (Mu et al. 2017) . The algorithm of this function is exactly same to the next_mtd_gBOIN_Continuous() just the input parameter is used by the TB score
}
\usage{
next_gBOIN_TB(target, n, y, d, mu_1 = 0.6 * target, mu_2 = 1.4 * target)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target TB score
  }
  \item{n}{
the number of patients enrolled at each dose level
  }
  \item{y}{
the toxicity score at each dose level
  }
  \item{d}{
the current dose level
  }
  \item{mu_1}{
the lower bound. The default value is 0.6 * target
  }
  \item{mu_2}{
the upper bound. The default value is 1.4 * target
  }
}

\value{
\code{next_gBOIN_TB()} returns recommended dose level for the next cohort as a numeric value under ordinal measure
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
d <- 2
next_gBOIN_TB(target = target, n = n, y = y, d = d)
}

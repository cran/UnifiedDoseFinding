\name{next_gBOIN_continuous}
\alias{next_gBOIN_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a maximum tolerated dose (MTD) using gBOIN design
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MTD under continuous measure using gBOIN design (Mu et al., 2017)
}
\usage{
next_gBOIN_continuous(target, n, y, d, mu_1 = 0.6 * target, mu_2 = 1.4 * target)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the continuous target score
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
\code{next_gBOIN_continuous()} returns recommended dose level for the next cohort as a numeric value under continuous measure
}
\references{
Mu, Rongji, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, and Jun Yin. "gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points." Journal of the Royal Statistical Society. Series C: Applied Statistics 68, no. 2 (2019): 289-308.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
target <- 1.47
n <- c(3, 3, 3, 9, 0, 0)
y <- c(0.1951265, 1.5434317, 2.1967343, 13.9266838, 0, 0)
d <- 4
next_gBOIN_continuous(target = target, n = n, y = y, d = d)
}

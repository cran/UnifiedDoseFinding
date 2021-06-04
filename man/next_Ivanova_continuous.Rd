\name{next_Ivanova_continuous}
\alias{next_Ivanova_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients using Inanova design
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MTD
}
\usage{
next_Ivanova_continuous(target, eps, c_resp, n, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity score
  }
  \item{eps}{
the decision criterion
  }
  \item{c_resp}{
the list object. Each element contains continuous value for each measurement at the certain dose level
  }
  \item{n}{
the number of patients enrolled at each dose level
  }
  \item{d}{
the current dose level
  }
}

\value{
\code{next_Ivanova_continuous()} returns recommended dose level for the next cohort as a numeric value
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}

\examples{
target <- 1.47
eps <- 1
c_resp <- list(c(0, 0.05475884, 0.12446843, 0.10131912),
               c(0, 0.4716962, 0.2792428, 0.3296575),
               c(0, 0.3931168, 1.6116607, 0.1642561),
               c(0, 0.9410027, 1.6021326, 1.6115235,
                 1.1735981, 2.5575655, 1.6513679, 1.4269044,
                 0.8983843, 2.2209587),
               0,
               0)
n <- c(3, 3, 3, 9, 0, 0)
d <- 4
next_Ivanova_continuous(target = target, eps = eps, c_resp = c_resp,
                        n = n, d = d)
}


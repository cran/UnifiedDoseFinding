\name{next_Ivanova_binary}
\alias{next_Ivanova_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients of binary endpoint using design by Ivanova et al (2009)
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials of binary endpoint that aim to find a MTD using design by Ivanova et al (2009)
}
\usage{
next_Ivanova_binary(target, eps, y, n, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity rate
  }
  \item{eps}{
the decision criterion
  }
  \item{y}{
the number of toxicity patients at each dose level
  }
  \item{n}{
the number of patients enrolled at each dose level

  }
  \item{d}{
the current dose level
  }
}

\value{
\code{next_Ivanova_binary()} returns recommended dose level for the next cohort as a numeric value
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}


\examples{
target <- 0.3
eps <- 1
y <- c(0, 4, 0, 0, 0, 0)
n <- c(3, 15, 0, 0, 0, 0)
d <- 2
next_Ivanova_binary(target = target, eps = eps, y = y, n = n, d = d)

}


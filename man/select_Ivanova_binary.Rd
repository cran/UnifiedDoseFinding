\name{select_Ivanova_binary}
\alias{select_Ivanova_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD) of binary endpoint for single agent trials using design by Ivanova et al (2009)
}
\description{
Select the maximum tolerated dose (MTD) when the trial is completed for binary endpoint using design by Ivanova et al (2009)
}
\usage{
select_Ivanova_binary(target, y, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity rate
  }
  \item{y}{
the number of toxicity patients at each dose level
  }
  \item{n}{
the number of patients enrolled at each dose level
  }
}

\value{
\code{select_Ivanova_binary()} returns a list object including: (1) dose selected (2) patients treated at each dose level
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}

\examples{
target<- 0.3
y <- c(0, 9, 0, 0, 0, 0)
n <- c(3, 22, 0, 0, 0, 0)
select_Ivanova_binary(target = target, y = y, n = n)

}


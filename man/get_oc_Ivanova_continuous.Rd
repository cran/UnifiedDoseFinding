\name{get_oc_Ivanova_continuous}
\alias{get_oc_Ivanova_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD) of continuous endpoint using design by Ivanova et al (2009)
}
\description{
Obtain the operating characteristics of the dose-finding design of continuous endpoint by Ivanova et al (2009)
}
\usage{
get_oc_Ivanova_continuous(target, eps = 1, ptox, ncohort,
                          cohortsize, n.earlystop = 100,
                          ntrial, startdose = 1, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the continuous target score
}
  \item{eps}{
the decision criterion. The default value is eps = 1
  }
  \item{ptox}{
the true mean value of the continuous measure
}
  \item{ncohort}{
the number of cohorts
}
  \item{cohortsize}{
the cohort size
}
  \item{n.earlystop}{
the early stopping parameter. The default value is n.earlystop = 100
}
  \item{ntrial}{
the number of simulated trials
}
  \item{startdose}{
the starting dose level. The default value is startdose = 1
}
  \item{seed}{
the seed. The default value is seed = 100
  }
}

\value{
\code{get_oc_Ivanova_continuous()} returns the operating characteristics of Ivanova design as a list object, including: (1) selection percentage at each dose level (2) patients treated at each dose level
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Fang Wang, Haitao Pan, Rongji Mu
}



\examples{
target <- 1.47
ptox <- c(0.11, 0.25, 0.94, 1.47, 2.38, 2.40)
ncohort <- 10
cohortsize <- 3
ntrial <- 4000
get_oc_Ivanova_continuous(target = target, ptox = ptox, ncohort = ncohort,
                          cohortsize = cohortsize, ntrial = ntrial)

}


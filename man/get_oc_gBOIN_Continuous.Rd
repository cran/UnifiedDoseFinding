\name{get_oc_gBOIN_Continuous}
\alias{get_oc_gBOIN_Continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD) using gBOIN design
}
\description{
Obtain the operating characteristics of the general Bayesian optimal interval (gBOIN) design (Mu et al. 2017) for maximum tolerated dose (MTD)-based dosing-finding trials under the continuous measure
}
\usage{
get_oc_gBOIN_Continuous(target, c_true, ncohort, cohortsize, ntrial, mu_1 = 0.6 * target,
                  mu_2 = 1.4 * target, startdose = 1, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the continuous target score
  }
  \item{c_true}{
the true mean value of the continuous measure
  }
  \item{ncohort}{
the number of cohorts
  }
  \item{cohortsize}{
the cohort size
  }
  \item{ntrial}{
the number of simulated trials
  }
  \item{mu_1}{
the lower bound. The default value is mu_1 = 0.6 * target
  }
  \item{mu_2}{
the upper bound. The default value is mu_2 = 1.4 * target
  }
  \item{startdose}{
the starting dose level. The default value is startdose = 1
  }
  \item{seed}{
the seed. The default value is seed = 100
  }

}
\value{
\code{get_oc_gBOIN_Continuous()} returns the operating characteristics of generalized Bayesian optimal interval design (gBOIN) as a list object, including: (1) selection percentage of each dose, (2) the average number of patients treated at each dose

}
\references{
Mu, Rongji, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, and Jun Yin. "gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points." Journal of the Royal Statistical Society. Series C: Applied Statistics 68, no. 2 (2019): 289-308.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}

\examples{
target <- 1.47
c_true <- c(0.11, 0.25, 0.94, 1.47, 2.38, 2.40)
ncohort <- 15
cohortsize <- 1
ntrial <- 4000
get_oc_gBOIN_Continuous(target = target, c_true = c_true, ncohort = ncohort,
                  cohortsize = cohortsize, ntrial = ntrial)
}

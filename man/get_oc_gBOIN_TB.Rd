\name{get_oc_gBOIN_TB}
\alias{get_oc_gBOIN_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD) defined by Toxicity Burden (TB) Score  using gBOIN design
}
\description{
Obtain the operating characteristics of the generalized Bayesian optimal interval (gBOIN) design (Mu et al. 2017) for maximum tolerated dose (MTD) (defined by the toxicity burden (BT) score proposed by Bekele et al. (2004))-based dosing-finding trials using. The algorithm of this function is exactly same to the get_oc_gBOIN_Continuous() just the input parameter is used by the TB score
}
\usage{
get_oc_gBOIN_TB(target, pmat, weight, ncohort, cohortsize, ntrial, mu_1 = 0.6 * target,
               mu_2 = 1.4 * target, startdose = 1, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target TB score
  }
  \item{pmat}{
pmat is a list. Each element is a matrix, representing the probability of different toxicity type and scale under
different dose levels.
  }
  \item{weight}{
the severity weight
  }
  \item{ncohort}{
the number of cohort
  }
  \item{cohortsize}{
the cohort size
  }
  \item{ntrial}{
the number of simulated trial
  }
  \item{mu_1}{
the lower bound. The default value is p.saf = 0.6 * target
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
\code{get_oc_gBOIN_TB()} returns the operating characteristics of  generalized Bayesian optimal interval design as a list object, including: (1) selection percentage of each dose, (2) the average number of patients treated at each dose
}
\references{
Bekele, B. Nebiyou, and Peter F. Thall. "Dose-finding based on multiple toxicities in a soft tissue sarcoma trial." Journal of the American Statistical Association 99, no. 465 (2004): 26-35.

Rongji Mu, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, Jun Yin: gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points. Royal Statistical Society 2019

}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
target <- 3.344
ncohort <- 9
cohortsize <- 4
ntrial <- 1000
rate <- 1.1
weight <- rate * rbind(c(0,1,1.5,5,6), c(0,2.5,6,rep(0,2)), c(0,2,3,6,0),
                       c(0,1.5,2,0,0), c(0,0.5,1,0,0))
pmat <- list()
pmat[[1]] <- rbind(c(0.5,0.5,rep(0,3)),
                   c(1,rep(0,4)),
                   c(1,rep(0,4)),
                   c(1,rep(0,4)),
                   c(0.5,0,0.5,0,0))
pmat[[2]] <- rbind(c(0.5,0,0.5,0,0),
                   c(1,rep(0,4)),
                   c(0.5,0.5,0,0,0),
                   c(0.5,0.5,rep(0,3)),
                   c(0.46,0,0.54,rep(0,2)))
pmat[[3]] <- rbind(c(0.5,0,0.5,0,0),
                   c(0.4,0.6,0,0,0),
                   c(0.25,0.75,0,0,0),
                   c(0.5,0.5,0,0,0),
                   c(1,0,0,0,0))
pmat[[4]] <- rbind(c(0.5,0,0.5,0,0),
                   c(0.4,0.6,0,0,0),
                   c(0.25,0.75,0,0,0),
                   c(0.5,0.5,0,0,0),
                   c(0.5,0,0.5,0,0))
pmat[[5]] <- rbind(c(0.5,0,0.5,0,0),
                   c(0,1,0,0,0),
                   c(0.25,0.75,0,0,0),
                   c(0.5,0.5,0,0,0),
                   c(0.5,0,0.5,0,0))
pmat[[6]] <- rbind(c(0,0.5,0.5,0,0),
                   c(0,1,0,0,0),
                   c(0,1,0,0,0),
                   c(0.5,0.5,0,0,0),
                   c(0.5,0,0.5,0,0))
pmat[[7]] <- rbind(c(0,0.5,0.5,0,0),
                   c(0,1,0,0,0),
                   c(0,1,0,0,0),
                   c(0,0.5,0.5,0,0),
                   c(0.5,0,0.5,0,0))
pmat[[8]] <- rbind(c(0,0.5,0.5,0,0),
                   c(0,1,0,0,0),
                   c(0,0,1,0,0),
                   c(0,0.5,0.5,0,0),
                   c(0.5,0,0.5,0,0))
pmat[[9]] <- rbind(c(0,0,1,0,0),
                   c(0,1,0,0,0),
                   c(0,0,1,0,0),
                   c(0,0,1,0,0),
                   c(0,0,1,0,0))
pmat[[10]] <- rbind(c(0,0,1,0,0),
                    c(0,1,0,0,0),
                    c(1/3,0,0,2/3,0),
                    c(0,0,1,0,0),
                    c(0,0,1,0,0))
get_oc_gBOIN_TB(target = target, pmat = pmat, weight = weight, ncohort = ncohort,
               cohortsize = cohortsize, ntrial = ntrial)
}


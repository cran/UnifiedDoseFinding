\name{get_oc_Quasi_CRM}
\alias{get_oc_Quasi_CRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD) defined by Equivalent Score (ET) using Quasi-CRM design
}
\description{
Obtain the operating characteristics of Quasi-CRM design (Yuan et al. 2007) and Robust-Quasi-CRM design (Pan et al. 2014) for finding the maximum tolerated dose (MTD) using Equivalent Score (ET) derived from toxicity grade information
}
\usage{
get_oc_Quasi_CRM(ptox, skeletons, target, score, cohortsize, ncohort,
                 start.dose = 1, mselection = 1, cutoff.eli = 0.90,
                 ntrial = 10, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptox}{
true toxicity probability at each dose level
  }
  \item{skeletons}{
a matrix to provide multiple skeletons with each row presenting a skeleton. If just one row, the function implements the Quasi-CRM design; if >=2 rows, the function implements the Robust-Quasi-CRM designn
  }
  \item{target}{
the target toxicity score
  }
  \item{score}{
the vector weight for ordinal toxicity levels
  }
  \item{cohortsize}{
the cohort size
  }
  \item{ncohort}{
the number of cohort
  }
  \item{start.dose}{
the starting dose level. The default value is start.dose = 1
  }
  \item{mselection}{
mselection = 1 (or 0) indicate to use Bayesian model selection (or mode averaging) to make inference across multiple skeletons. The default value is mselection = 1. It only applies to the Robust-Quasi-CRM design
  }
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. The default value is cutoff.eli = 0.90
  }
  \item{ntrial}{
the number of simulated trials. The default value is ntrial = 10
  }
  \item{seed}{
the seed. The default value is seed = 100
  }
}

\value{
\code{get_oc_Quasi_CRM()} returns the operating characteristics of (Robust)-Quasi-CRM design as a list object, including: (1) selection percentage at each dose level (2) patients treated at each dose level
}
\references{
Yuan, Z., R. Chappell, and H. Bailey. "The continual reassessment method for multiple toxicity grades: a Bayesian quasi-likelihood approach." Biometrics 63, no. 1 (2007): 173-179.

Pan, Haitao, Cailin Zhu, Feng Zhang, Ying Yuan, Shemin Zhang, Wenhong Zhang, Chanjuan Li, Ling Wang, and Jielai Xia. "The continual reassessment method for multiple toxicity grades: a Bayesian model selection approach." PloS one 9, no. 5 (2014): e98147.

}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
### Scenario 1 in Yuan et al. (2007) and Pan et al. (2014)
target <- 0.47
score <- c(0, 0.5, 1, 1.5)
cohortsize <- 1
ncohort <- 20

ptox <- matrix(nrow = 4, ncol = 6)
ptox[1,] <- c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19)
ptox[2,] <- c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11)
ptox[3,] <- c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11)
ptox[4,] <- c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)


### specify one skeleton (Quasi-CRM design)
p1 <- c(0.11, 0.25, 0.40, 0.55, 0.75, 0.85)

get_oc_Quasi_CRM(ptox = ptox, skeletons = p1, target = target,
                 score = score, cohortsize = cohortsize,
                 ncohort = ncohort)




###########################################

### specify three skeletons (Quasi-CRM design)
p1 <- c(0.11, 0.25, 0.40, 0.55, 0.75, 0.85)
p2 <- c(0.05, 0.10, 0.15, 0.25, 0.40, 0.65)
p3 <- c(0.20, 0.40, 0.60, 0.75, 0.85, 0.95)
skeletons <- rbind(p1, p2, p3)


get_oc_Quasi_CRM(ptox = ptox, skeletons = skeletons, target = target,
                 score = score, cohortsize = cohortsize,
                 ncohort = ncohort)


}

\name{next_mtd_Quasi_CRM}
\alias{next_mtd_Quasi_CRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients using Quasi-CRM design
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MTD defined by the Equivalent Score (ET) using Quasi-CRM design (Yuan et al. 2007) and Robust-Quasi-CRM design (Pan et al. 2014)
}
\usage{
next_mtd_Quasi_CRM(target, n, y, dose.curr, score, skeleton,
                          cutoff.eli = 0.90, mselection = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity score
  }
  \item{n}{
the number of patients treated at each dose level
  }
  \item{y}{
the toxicity score at each dose level
  }
  \item{dose.curr}{
the current dose level
  }
  \item{score}{
the vector weight for ordinal toxicity levels
  }
  \item{skeleton}{
a matrix to provide multiple skeletons with each row presenting a skeleton
  }
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. The default value is cutoff.eli = 0.90
  }
  \item{mselection}{
mselection = 1 (or 0) indicate to use Bayesian model selection (or mode averaging) to make inference across multiple skeletons. The default value is mselection = 1
  }
}

\value{
\code{next_mtd_Quasi_CRM()} returns recommended dose level for the next cohort as a numeric value
}
\references{
Yuan, Z., R. Chappell, and H. Bailey. "The continual reassessment method for multiple toxicity grades: a Bayesian quasi-likelihood approach." Biometrics 63, no. 1 (2007): 173-179.

Pan, Haitao, Cailin Zhu, Feng Zhang, Ying Yuan, Shemin Zhang, Wenhong Zhang, Chanjuan Li, Ling Wang, and Jielai Xia. "The continual reassessment method for multiple toxicity grades: a Bayesian model selection approach." PloS one 9, no. 5 (2014): e98147.

}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}
\examples{
### Implement Robust-Quasi-CRM design (Pan et al. 2014) with pre-specifying 3 skeletons
target <- 0.47
score <- c(0, 0.5, 1, 1.5)
p1 <- c(0.11, 0.25, 0.40, 0.55, 0.75, 0.85)
p2 <- c(0.05, 0.10, 0.15, 0.25, 0.40, 0.65)
p3 <- c(0.20, 0.40, 0.60, 0.75, 0.85, 0.95)
skeletons <- rbind(p1, p2, p3)
n <- c(1, 1, 1, 1, 15, 1)
y <- c(0.000000, 0.000000, 0.000000, 0.000000, 4.666667, 1.000000)
next_mtd_Quasi_CRM (target = target, n = n, y = y, dose.curr = 5,
                          score = score, skeleton = skeletons)

### Implement Quasi-CRM design (Yuan et al. 2007) with pre-specifying/using 1 skeletons
next_mtd_Quasi_CRM (target = target, n = n, y = y, dose.curr = 5,
                          score = score, skeleton = p1)

}

\name{next_mtd_QuasiBOIN}
\alias{next_mtd_QuasiBOIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients based on equivalent score (ET)-defined target using gBOIN design
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MTD defined by the Equivalent Score (ET) in Quasi-CRM design (Yuan et al. 2007) and Robust-Quasi-CRM design (Pan et al. 2014) using the gBOIN design (Mu et al. 2017)
}
\usage{
next_mtd_QuasiBOIN(target, n, y, d, p.saf = 0.6 * target, p.tox = 1.4 * target,
                   cutoff.eli = 0.95, extrasafe = FALSE, n.earlystop = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target DLT rate
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
  \item{p.saf}{
the lower bound. The default value is p.saf = 0.6 * target
  }
  \item{p.tox}{
the upper bound. The default value is p.tox = 1.4 * target
  }
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. The default value is cutoff.eli = 0.95
  }
  \item{extrasafe}{
extrasafe	set extrasafe = TRUE to impose a more stringent stopping . The default value is extrasafe = FALSE
  }
  \item{n.earlystop}{
the early stopping parameter. The default value is n.earlystop = 100
  }
}
\value{
\code{next_mtd_QuasiBOIN()} returns recommended dose level for the next cohort as a numeric value under quasi-binary measure
}
\references{
Yuan, Z., R. Chappell, and H. Bailey. "The continual reassessment method for multiple toxicity grades: a Bayesian quasi-likelihood approach." Biometrics 63, no. 1 (2007): 173-179.

Pan, Haitao, Cailin Zhu, Feng Zhang, Ying Yuan, Shemin Zhang, Wenhong Zhang, Chanjuan Li, Ling Wang, and Jielai Xia. "The continual reassessment method for multiple toxicity grades: a Bayesian model selection approach." PloS one 9, no. 5 (2014): e98147.

Mu, Rongji, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, and Jun Yin. "gBOIN: a unified model-assisted phase I trial design accounting for toxicity grades, and binary or continuous end points." Journal of the Royal Statistical Society. Series C: Applied Statistics 68, no. 2 (2019): 289-308.

}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}



\examples{
target <- 0.47 / 1.5
n <- c(3, 3, 6, 9, 9, 0)
y <- c(0.000000, 0.000000, 1.333333, 2.333333, 3.666667, 0.0000000)
d <- 4
next_mtd_QuasiBOIN(target = target, n = n, y = y, d = d)

}

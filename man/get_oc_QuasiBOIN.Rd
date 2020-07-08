\name{get_oc_QuasiBOIN}
\alias{get_oc_QuasiBOIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD) defined by Equivalent Score (ET) using Quasi-CRM design using gBOIN
}
\description{
Obtain the operating characteristics of Quasi-CRM design (Yuan et al. 2007) and Robust-Quasi-CRM design (Pan et al. 2014) for finding the maximum tolerated dose (MTD) using Equivalent Score (ET) derived from toxicity grade information using the gBOIN design (Mu et al. 2017)
}
\usage{
get_oc_QuasiBOIN(target, p.true, score, ncohort, cohortsize, n.earlystop = 100,
                 startdose = 1, p.saf = 0.6 * target, p.tox = 1.4 * target,
                 cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05,
                 ntrial = 1000, seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target DLT rate
  }
  \item{p.true}{
the true toxicity probability at each dose level
  }
  \item{score}{
the default value is score = seq(0, 1.5, by = 0.5) / 1.5
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
  \item{startdose}{
the starting dose level. The default value is startdose = 1
  }
  \item{p.saf}{
lower bound. The default value is p.saf = 0.6 * target
  }
  \item{p.tox}{
upper bound. The default value is p.tox = 1.4 * target
  }
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. The default value is cutoff.eli = 0.95
  }
  \item{extrasafe}{
extrasafe	set extrasafe = TRUE to impose a more stringent stopping rule. The default value is extrasafe = FALSE
  }
  \item{offset}{
when extrasafe = TRUE will have effect. The default value is offset = 0.05
  }
  \item{ntrial}{
the number of simulated trials
  }
  \item{seed}{
the seed. The default value is seed = 100
  }
}

\value{
\code{get_oc_QuasiBOIN()} returns the operating characteristics of Bayesian optimal interval design as a list object,
       including:
       (1) the target DLT rate,
       (2) the true DLT rate at different scale for each dose level,
       (3) number of cohort,
       (4) cohortsize,
       (5) starting dose level,
       (6) lower bound,
       (7) upper bound,
       (8) selection percentage of each dose level,
       (9) the average number of patients treated at each dose,
       (10) the average number of patients responded to toxicity at each dose level
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
p.true <- matrix(c(0.83, 0.12, 0.04, 0.01,
                   0.75, 0.15, 0.07, 0.03,
                   0.62, 0.18, 0.11, 0.09,
                   0.51, 0.19, 0.14, 0.16,
                   0.34, 0.16, 0.15, 0.35,
                   0.19, 0.11, 0.11, 0.59), ncol = 4, byrow = TRUE)
score <- seq(0, 1.5,by = 0.5) / 1.5
ncohort <- 10
cohortsize <- 3
ntrial <- 4000
get_oc_QuasiBOIN(target = target, p.true = p.true, score = score, ncohort = ncohort,
                 cohortsize = cohortsize, ntrial = ntrial)

}


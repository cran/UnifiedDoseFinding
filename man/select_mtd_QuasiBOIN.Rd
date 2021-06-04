\name{select_mtd_QuasiBOIN}
\alias{select_mtd_QuasiBOIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD)-defined by equivalent score (ET) using gBOIN design
}
\description{
Select the maximum tolerated dose (MTD) defined by the Equivalent Score (ET) in Quasi-CRM design (Yuan et al. 2007) and Robust-Quasi-CRM design (Pan et al. 2014) when the trial is completed using the gBOIN design (Mu et al. 2017)
}
\usage{
select_mtd_QuasiBOIN(target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE,
                     offset = 0.05, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target DLT rate
}
  \item{npts}{
the number of patients enrolled at each dose level
}
  \item{ntox}{
the toxicity score at each dose level
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
  \item{print}{
print the additional result or not. The default value is print = FALSE
}
}

\value{
\code{select_mtd_QuasiBOIN()} returns the selected dose
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
y <- c(0, 0, 1.333333, 2.333333, 3.666667, 0)
select_mtd_QuasiBOIN(target = target, npts = n, ntox = y)
}

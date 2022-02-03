\name{select_mtd_Ivanova_continuous}
\alias{select_mtd_Ivanova_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD) for single agent trials of continuous endpoint using design by Ivanova et al (2009)
}
\description{
Select the maximum tolerated dose (MTD) when the trial is completed for continuous endpoint using design by Ivanova et al (2009)
}
\usage{
select_mtd_Ivanova_continuous(target, c_resp, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
the target toxicity score
  }
  \item{c_resp}{
list object. Each element contains continuous value for each measurement
  }
  \item{n}{
the number of patients enrolled at each dose level
  }
}

\value{
\code{select_mtd_Ivanova_continuous()} returns a list object including: (1) dose selected (2) patients treated at each dose level
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}

\examples{
target <- 1.47
c_resp <- list(c(0, 0.05475884, 0.12446843, 0.10131912),
               c(0, 0.4716962, 0.2792428, 0.3296575),
               c(0, 0.3931168, 1.6116607, 0.1642561),
               c(0, 0.9410027, 1.6021326, 1.6115235,
                 1.1735981, 2.5575655, 1.6513679, 1.4269044,
                 0.8983843, 2.2209587),
               0,
               0)
n <- c(3, 3, 3, 9, 0, 0)
select_mtd_Ivanova_continuous(target = target, c_resp = c_resp, n = n)
}


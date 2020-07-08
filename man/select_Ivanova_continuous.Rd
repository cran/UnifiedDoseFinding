\name{select_Ivanova_continuous}
\alias{select_Ivanova_continuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select the maximum tolerated dose (MTD) for single agent trials of continuous endpoint using design by Ivanova et al (2009)
}
\description{
Select the maximum tolerated dose (MTD) when the trial is completed for continuous endpoint using design by Ivanova et al (2009)
}
\usage{
select_Ivanova_continuous(target, c_resp, n)
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
\code{select_Ivanova_continuous()} returns a list object including: (1) dose selected (2) patients treated at each dose level
}
\references{
Ivanova, Anastasia, and Se Hee Kim. "Dose finding for continuous and ordinal outcomes with a monotone objective function: a unified approach." Biometrics 65, no. 1 (2009): 307-315.
}
\author{
Chia-Wei Hsu, Haitao Pan, Rongji Mu
}

\examples{
target <- 0.3
c_resp <- list(c(0.14934229, 0.33945935, 0.27632487, 0.12546279,
                 0.51435981, 0.05242217, 0.23950981, 0.52195215,
                 0.33701385, 0.02585574, 0.99308905, 0.16857301,
                 0.52921818, 0.37858839, 0.53202138, 0.05568626,
                 0.16846483, 0.08393353, 0.36928336, -0.04731884, 0.37412280),
               c(0.7547139, 0.4467885, 0.5274520, 0.2560551, 0.4359545,
                 0.4385098, 0.3882733, 0.2444583, 0.6043425),
               0,
               0,
               0,
               0)
n <- c(21, 9, 0, 0, 0, 0)
select_Ivanova_continuous(target = target, c_resp = c_resp, n = n)
}


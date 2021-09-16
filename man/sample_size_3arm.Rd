\name{sample_size_3arm}
\alias{sample_size_3arm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the sample size for each arm
}
\description{
Calculate the sample size for each arm in a three-arm trial
}
\usage{
sample_size_3arm(kappa, m0, mA, mB, delta, ta, tf, P, diff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
the shape parameter of weibull distribution
  }
  \item{m0}{
the median survival time of historical control with respect to both arm A, B and C
  }
  \item{mA}{
the median survival time of arm A
  }
  \item{mB}{
the median survival time of arm B
  }
  \item{delta}{
the hazard ratio of arm C to arm A
  }
  \item{ta}{
the accrual duration
  }
  \item{tf}{
the follow-up duration
  }
  \item{P}{
the chance of correctly selecting the superior arm
  }
  \item{diff}{
the equivalence margin
  }
}

\value{
\code{sample_size_3arm()} returns required sample size for each arm
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Two-Stage Screened Selection Designs for Randomized
Phase II Trials with Time-to-Event Endpoints." Submitted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}

\examples{
sample_size_3arm(kappa = 1, m0 = 1, mA = 2, mB = 2, delta = 0.6,
                 ta = 24, tf = 3, P = 0.90, diff = 0)

sample_size_3arm(kappa = 1, m0 = 1, mA = 2, mB = 2, delta = 0.6,
                 ta = 24, tf = 3, P = 0.90, diff = 0.25)
}

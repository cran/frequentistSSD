\name{sample_size}
\alias{sample_size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the sample size for each arm
}
\description{
Calculate the sample size for each arm in a two-arm trial
}
\usage{
sample_size(kappa, mA, hr, ta, tf, diff, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
the shape parameter of weibull distribution
  }
  \item{mA}{
the median survival time of arm A
  }
  \item{hr}{
the hazard ratio of arm B to arm A
  }
  \item{ta}{
the accrual duration
  }
  \item{tf}{
the follow-up duration
  }
  \item{diff}{
the equivalence margin
  }
  \item{P}{
the chance of correctly selecting the superior arm
  }
}

\value{
\code{sample_size()} returns required sample size for each arm
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Two-Stage Screened Selection Designs for Randomized
Phase II Trials with Time-to-Event Endpoints." Submitted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}


\examples{
sample_size(kappa = 1, mA = 2.4, hr = 0.60, ta = 24,
            tf = 3, diff = 0.25, P = 0.9)
}

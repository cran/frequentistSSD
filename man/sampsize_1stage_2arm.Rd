\name{sampsize_1stage_2arm}
\alias{sampsize_1stage_2arm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the sample size for each arm
}
\description{
Provide the required sample size (per arm) based on the one-stage selection design given the type I error rate and power. Details see the reference paper
}
\usage{
sampsize_1stage_2arm(kappa, alpha, beta, m0, m1, ta, tf, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
the shape parameter of weibull distribution
  }
  \item{alpha}{
type I error rate
  }
  \item{beta}{
type II error rate
  }
  \item{m0}{
the median survival time of historical control with respect to both arm A and B
  }
  \item{m1}{
the median survival time of arm A
  }
  \item{ta}{
the accrual duration
  }
  \item{tf}{
the follow-up duration
  }
  \item{delta}{
the hazard ratio of arm B to arm A
  }

}

\value{
\code{sampsize_1stage_2arm()} returns required sample size for each arm
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Two-Stage Screened Selection Designs for Randomized
Phase II Trials with Time-to-Event Endpoints." Submitted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}


\examples{
sampsize_1stage_2arm(kappa = 1, alpha = 0.05, beta = 0.1, m0 = 1,
                     m1 = 2.5, delta = 0.65, ta = 24, tf = 9)
}

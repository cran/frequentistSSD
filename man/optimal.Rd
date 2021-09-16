\name{optimal}
\alias{optimal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find optimal parameters
}
\description{
Find the optimal parameters used in the \code{get_oc()} function
}
\usage{
optimal(shape, m0, ms, tf, ta,
        tot_size, dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
the shape parameter of weibull distribution
  }
  \item{m0}{
the median survival time of historical data
  }
  \item{ms}{
the minimal clinical meaningful median survival time
  }
  \item{tf}{
the follow-up duration
  }
    \item{ta}{
the accrual duration
  }
  \item{tot_size}{
the required sample size for each arm
  }
  \item{dist}{
the distribution
  }

}

\value{
\code{optimal()} returns a list containing two blocks:

\code{$param}:
(1) S0: historical survival probability at the landmark time point x0 (2) hr: hazard ratio of the hypothetical arm with “minimal clinical meaningful median survival time” to the historical arm (3) rate: accrual rate (calculated as tot_size/ta) \cr \cr
\code{$Two_stage}:
(1) alpha: type I error (2) beta: type II error (3) n1: the required sample size for the first stage (4) c1: the cutoff point at the first stage (5) n: the required sample size for each arm (is equal to tot_size) (6) c: the cutoff point at the final stage (7) t1: the expected time of interim analysis (first stage) (8) MTSL: the maximum total study length (ta + tf) (9) ESS: the expected sample size for each arm (10) PS: the probability of early stopping
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Two-Stage Screened Selection Designs for Randomized
Phase II Trials with Time-to-Event Endpoints." Submitted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}
\note{
(1) For \code{optimal()} function, if the ms is set to be equal to the less effective arm (e.g., say, the arm A is less effective than the arm B), with the rationale, by this way, the computed boundaries, c1, n1, c, would control the pre-specified type I (since both arms A and B are all under the same null case) while the power for the better arm B, would have at least pre-specified power 1-beta

(2) \code{optimal()} funcition will be applied to either two or three-arm cases to find the optimal cutoff bounaries, c1, n1, c
}


\examples{
## For 2 arms
\dontrun{
shape <- 1
m0 <- 1
ms <- 2.4
tf <- 3
ta <- 24
tot_size <- 19
dist <- "WB"

optimal(shape = shape, m0 = m0, ms = ms, tf = tf,
        ta = ta, tot_size = tot_size, dist = dist)
}

## For 3 arms
\dontrun{
shape <- 1
m0 <- 1
ms <- 2
tf <- 3
ta <- 24
tot_size <- 21
dist <- "WB"

optimal(shape = shape, m0 = m0, ms = ms, tf = tf,
        ta = ta, tot_size = tot_size, dist = dist)
}
}

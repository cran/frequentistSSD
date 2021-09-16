\name{get_oc}
\alias{get_oc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for Two-Stage Screened Selection Design for Randomized Phase II Trials with Time-to-Event Endpoints
}
\description{
Obtain the operating characteristics of Two-Stage Screened Selection Design for Randomized Phase II Trials with Time-to-Event Endpoints. The arguments for this function are from outputs of the functions of \code{optimal()} and \code{sample_size()}
}
\usage{
get_oc(shape, m0, mA, hr, frac, ta, tf, c1, c,
       diff, n, nsim, seed = 2483)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
the shape parameter of weibull distribution
  }
  \item{m0}{
the median survival time of historical data
  }
  \item{mA}{
the median survival time of arm A
  }
  \item{hr}{
the hazard ratio of arm B to arm A
  }
  \item{frac}{
the ratio of sample size in stage 1 compared to total sample size (e.g., n1/n)
  }
  \item{ta}{
the accrual duration
  }
  \item{tf}{
the follow-up duration
  }
  \item{c1}{
the critical value at stage 1
  }
  \item{c}{
the critical value at stage 2
  }
  \item{diff}{
the equivalence margin
  }
  \item{n}{
the total sample size for each arm
  }
  \item{nsim}{
the number of simulated trials
  }
  \item{seed}{
the seed. The default value is seed = 2483
  }
}

\value{
\code{get_oc()} returns: (1) n: total sample size for each arm (2) SSD.Arm.A: selection probability of Arm A (3) SSD.Arm.B: selection probability of Arm B (4) SSD.No.Arm: the probability of no arms selected (5) diff: the equivalence margin (6) Mean.N.Arm.A: the average number of patients allocated to Arm A (7) Mean.N.Arm.B: the average number of patients allocated to Arm B
}
\references{
Jianrong Wu, Haitao Pan, Chia-Wei Hsu (2021). "Two-Stage Screened Selection Designs for Randomized
Phase II Trials with Time-to-Event Endpoints." Submitted
}
\author{
Chia-Wei Hsu, Haitao Pan, Jianrong Wu
}

\examples{
## c1, n, c they all come from optimal function. Details are shown as follows:
\dontrun{optimal(shape = 1, m0 = 1, ms = 2.4, tf = 3,
                 ta = 24, tot_size = 19, dist = "WB")
| alpha | beta | n1 |   c1   |  n |    c   |   t1    | MTSL |   ESS   |   PS
|  0.05 | 0.08 | 15 | 0.0508 | 19 | 1.6434 | 17.9751 |  27  | 16.5183 | 0.5203
}

get_oc(shape = 1, m0 = 1, mA = 2.4, hr = 0.6, frac = 15/19,
       ta = 24, tf = 3, c1 = 0.0508, c = 1.6434,
       diff = 0.25, n = 19, nsim = 10000)
}

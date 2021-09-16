\name{get_oc_3arm}
\alias{get_oc_3arm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for Two-Stage Screened
Selection Design for Randomized Phase II Trials
with Time-to-Event Endpoints for 3 arms
}
\description{
Obtain the operating characteristics of Two-Stage Screened Selection Design for Randomized Phase II Trials with Time-to-Event Endpoints for 3 arms. The arguments for this function are from outputs of the functions of \code{optimal()} and \code{sample_size_3arm()}
}
\usage{
get_oc_3arm(shape, m0, mA, hr2, hr3, frac, ta, tf, c1, c,
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
  \item{hr2}{
the hazard ratio of arm B to arm A
  }
  \item{hr3}{
the hazard ratio of arm C to arm A
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
\code{get_oc_3arm()} returns: (1) n: total sample size for each arm (2) SSD.Arm.A: selection probability of Arm A (3) SSD.Arm.B: selection probability of Arm B (4) SSD.Arm.C: selection probability of Arm C (5) SSD.No.Arm: the probability of no arms selected (6) diff: the equivalence margin (7) Mean.N.Arm.A: the average number of patients allocated to Arm A (8) Mean.N.Arm.B: the average number of patients allocated to Arm B (9) Mean.N.Arm.C: the average number of patients allocated to Arm C
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
\dontrun{optimal(shape = 1, m0 = 1, ms = 2, tf = 3,
                 ta = 24, tot_size = 21, dist = "WB")
| alpha | beta | n1 |    c1   |  n |   c   |   t1    | MTSL |   ESS   |   PS
|  0.05 | 0.14 | 16 | -0.0463 | 21 | 1.644 | 17.8248 |  27  | 18.3983 | 0.4815
}

get_oc_3arm(shape = 1, m0 = 1, mA = 2, hr2 = 1, hr3 = 0.6,
            frac = 16/21, ta = 24, tf = 3, c1 = -0.0463,
            c = 1.644, diff = 0, n = 21, nsim = 10000,
            seed = 2483)
}

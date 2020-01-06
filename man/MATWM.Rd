\name{MATWM}
\alias{MATWM}
\title{
PERICLIMv1.0--MATWM
}
\description{
Derivation of palaeo-air temperature characteristics from thaw depth with the magnitude of annual air temperature oscillations defined by mean air temperature of the warmest month (MATWM)
}

\usage{
MATWM(z, kt, vmc, nt, MATWM, P, showInputs = TRUE)
}

\arguments{
  \item{z}{Thaw depth [m]}
  \item{kt}{Thawed ground thermal conductivity [W/m/K]}
  \item{vmc}{Volumetric ground moisture content [--]}
  \item{nt}{Thawing n-factor [--]}
  \item{MATWM}{Mean air temperature of the warmest month [degC]}
  \item{P}{Period of air temperature oscillations [d]}
  \item{showInputs = TRUE}{Should the inputs be included in the output (TRUE by default)?}
}

\details{
All the numerical inputs must be greater than 0. The moisture content has to be given as a decimal number with a maximum value of 1.
}

\value{
  \item{MAAT}{Mean annual air temperature [degC]}
  \item{MATWM}{Mean air temperature of the warmest month [degC]}
  \item{MATCM}{Mean air temperature of the coldest month [degC]}
  \item{MATTS}{Mean air temperature of the thawing season [degC]}
  \item{MATFS}{Mean air temperature of the freezing season [degC]}
  \item{Ita}{Air thawing index [degC.d]}
  \item{Ifa}{Air freezing index [degC.d]}
  \item{Lt}{Length of the thawing season [d]}
  \item{Lf}{Length of the freezing season [d]}
  \item{Its}{Ground-surface thawing index [degC.d]}
}

\section{Warning}{
This solution should be utilized cautiously because it is highly sensitive to uncertainties in the air thawing index and the mean air temperature of the warmest month.
}

\examples{
# Example code run with sample data:

PERICLIMv1.0::MATWM(z = 0.5, kt = 0.6, vmc= 0.3, nt = 1.2, MATWM = 3, P = 365, showInputs = TRUE)

# If executed correctly, the output and input variables should be:

# MAAT  = -5.301689
# MATWM = 3.000000
# MATCM = -13.603379
# MATTS = 1.973573
# MATFS = -8.123979
# Ita   = 201.340664
# Ifa   = -2136.457259
# Lt    = 102.018369
# Lf    = 262.981631
# Its   = 241.608796
# z     = 0.5
# kt    = 0.6
# wmc   = 0.3
# nt    = 1.2
# Aa    = 16.603379
# P     = 365
}

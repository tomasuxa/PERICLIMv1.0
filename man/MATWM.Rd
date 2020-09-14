\name{MATWM}

\alias{MATWM}

\title{PERICLIMv1.0::MATWM}

\description{
Derivation of palaeo-air temperature characteristics from active-layer thickness with the magnitude of annual air temperature oscillations defined by the mean air temperature of the warmest month (MATWM), or rather its difference from the mean annual air temperature.
}

\usage{MATWM(z, vmc, dbd, q, fc=c('fine','coarse'), nt, MATWM, showInputs = TRUE)}

\arguments{
  \item{z}{Active-layer thickness [m]}
  \item{vmc}{Volumetric ground moisture content [--]}
  \item{dbd}{Dry ground bulk density [kg/m3]}
  \item{q}{Ground quartz content [--]}
  \item{fc}{Ground grain-size class ['fine' or 'coarse']}
  \item{nt}{Ground-surface thawing n-factor [--]}
  \item{MATWM}{Mean air temperature of the warmest month [degC]}
  \item{showInputs = TRUE}{Should the inputs be included among the outputs (TRUE by default)?}
}

\details{
All the numerical inputs must be greater than 0. The moisture and quartz contents must be given as a decimal number with a maximum value of 1. A maximum allowed value for the bulk density is set at 2700 kg/m3.
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
This solution should be utilized cautiously because it is highly sensitive to uncertainties in the air thawing index as well as the mean air temperature of the warmest month.
}

\examples{
# Example code run with sample data:

PERICLIMv1.0::MATWM(z = 0.5, vmc = 0.3, dbd = 1500, q = 0.4, fc = 'fine', nt = 1.2, MATWM = 1.9, showInputs = TRUE)

# If executed correctly, the output and input variables should be:

# MAAT  = -8.06965559495226
# MATWM = 1.9
# MATCM = -18.0393111899045
# MATTS = 1.25823408572454
# MATFS = -10.3984330754712
# Ita   = 91.7503875413873
# Ifa   = -3037.17467969896
# Lt    = 72.9199666281127
# Lf    = 292.080033371887
# Its   = 110.100465049665
# z     = 0.5
# vmc   = 0.3
# dbd   = 1500
# q     = 0.4
# fc    = 'fine'
# kt    = 1.31666362819073
# nt    = 1.2
# Aa    = 19.9393111899045
}

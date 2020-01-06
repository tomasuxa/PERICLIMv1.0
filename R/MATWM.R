################################################################################################################
#
# Derivation of palaeo-air temperature characteristics from thaw depth
# with the magnitude of annual air temperature oscillations defined by mean air temperature of the warmest month
#
################################################################################################################

MATWM <- function(z, kt, vmc, nt, MATWM, P, showInputs = TRUE) {

  # Checking for physically feasible values of the input parameters
  if(z <= 0 |               # Thaw depth [m]
     kt <= 0 |              # Thawed ground thermal conductivity [W/m/K]
     vmc <= 0 | vmc > 1 |   # Volumetric ground moisture content [-]
     nt <= 0 |              # Thawing n-factor [-]
     MATWM <= 0 |           # Mean air temperature of the warmest month [degC]
     P <= 0) {              # Period of air temperature oscillations [d]

    return(NA)   # If any of the inputs is outside the physically feasible values then the solution returns NA

  } else {

    # Calculation of ground-surface (Its) and air (Ita) thawing index required to reach the specific thaw depth through an inverse solution of the Stefan equation
    It <- function(z, kt, vmc, nt) {
      Its <- (z^2 * 334000 * vmc * 1000)/(2 * kt * 86400)   # Ground-surface thawing index [degC.d]
      Ita <- Its/nt                                         # Air thawing index [degC.d]
      It <- list(Its, Ita)
      return(It)
    }

    # Calculation of mean annual air temperature (MAAT) based on the air thawing index (Ita) and the mean air temperature of the warmest month (MATWM)
    MAAT <- function(Ita, MATWM) {
      f4 <- function(T) {
        f3 <- function(T) {
          f2 <- function(T) {
            f1 <- function(t) {T + (MATWM - T) * sin(2 * pi * t/P)}
          }
          integrate(f2(T), lower = asin(-T/(MATWM - T)) * P/(2 * pi), upper = (pi - asin(-T/(MATWM - T))) * P/(2 * pi))$value
        }
        Ita - f3(T)
      }
      # Since the annual air temperature amplitude is unknown at this stage, MAAT is searched in the interval [-100,0],
      # which roughly meets the conditions for the presence of permafrost (the lower limit of -100 substitutes
      # negative infinity, which speeds up the solution); if the search fails then the solution returns NA
      MAAT <- tryCatch(uniroot(f4, lower = -100, upper = 0, tol = 0.001)$root, error=function(MAAT) NA)
      return(MAAT)
    }

    # Calculation of air temperature characteristics
    It <- as.numeric(It(z, kt, vmc, nt))                       # Thawing indices [degC.d]
    Its <- It[1]                                               # Ground-surface thawing index [degC.d]
    Ita <- It[2]                                               # Air thawing index [degC.d]
    MAAT <- MAAT(Ita, MATWM)                                   # Mean annual air temperature [degC]
    Ifa <- MAAT * P - Ita                                      # Air freezing index [degC.d]
    MATCM <- MAAT - (MATWM - MAAT)                             # Mean air temperature of the coldest month [degC]
    Aa <- MATWM - MATCM                                        # Annual air temperature amplitude [degC]
    Lt <- (pi - 2 * asin(-MAAT/(MATWM - MAAT))) * P/(2 * pi)   # Length of the thawing season [d]
    Lf <- P - Lt                                               # Length of the freezing season [d]
    MATTS <- Ita/Lt                                            # Mean air temperature of the thawing season [degC]
    MATFS <- Ifa/Lf                                            # Mean air temperature of the freezing season [degC]

    # List of model outputs and inputs
    if(showInputs==TRUE) {   # Output includes the input parameters

      OutputsInputs <- c(MAAT, MATWM, MATCM, MATTS, MATFS, Ita, Ifa, Lt, Lf, Its, z, kt, vmc, nt, Aa, P)

    } else {   # Output excludes the input parameters

      OutputsInputs <- c(MAAT, MATWM, MATCM, MATTS, MATFS, Ita, Ifa, Lt, Lf, Its)
    }

    if(any(is.na(OutputsInputs))) {

      return(NA)   # If any of the outputs or inputs is NA then the solution returns NA

    } else {

      return(OutputsInputs)
    }
  }
}

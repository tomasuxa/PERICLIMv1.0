############################################################################################
#
# Derivation of palaeo-air temperature characteristics from thaw depth
# with the magnitude of annual air temperature oscillations defined by temperature amplitude
#
############################################################################################

Aa <- function(z, kt, vmc, nt, Aa, P, showInputs = TRUE) {

  # Checking for physically feasible values of the input parameters
  if(z <= 0 |               # Thaw depth [m]
     kt <= 0 |              # Thawed ground thermal conductivity [W/m/K]
     vmc <= 0 | vmc > 1 |   # Volumetric ground moisture content [-]
     nt <= 0 |              # Thawing n-factor [-]
     Aa <= 0 |              # Annual air temperature amplitude [degC]
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

    # Calculation of mean annual air temperature (MAAT) based on the air thawing index (Ita) and the annual air temperature amplitude (Aa)
    MAAT <- function(Ita, Aa) {
      f4 <- function(T) {
        f3 <- function(T) {
          f2 <- function(T) {
            f1 <- function(t) {T + Aa/2 * sin(2 * pi * t/P)}
          }
          integrate(f2(T), lower = asin(-T/(Aa/2)) * P/(2 * pi), upper = (pi - asin(-T/(Aa/2))) *P/(2 * pi))$value
        }
        Ita - f3(T)
      }
      # MAAT is searched in the right-closed interval (-Aa/2,0], which roughly meets the conditions for the presence of permafrost,
      # but at the same time it ensures that both negative and positive air temperatures have occurred during the annual period,
      # which is an essential prerequisite for the active layer to form; if the search fails then the solution returns NA
      MAAT <- tryCatch(uniroot(f4, lower = -Aa/2 + 0.001, upper = 0, tol = 0.001)$root, error=function(MAAT) NA)
      return(MAAT)
    }

    # Calculation of air temperature characteristics
    It <- as.numeric(It(z, kt, vmc, nt))               # Thawing indices [degC.d]
    Its <- It[1]                                       # Ground-surface thawing index [degC.d]
    Ita <- It[2]                                       # Air thawing index [degC.d]
    MAAT <- MAAT(Ita, Aa)                              # Mean annual air temperature [degC]
    Ifa <- MAAT * P - Ita                              # Air freezing index [degC.d]
    MATWM <- MAAT + Aa/2                               # Mean air temperature of the warmest month [degC]
    MATCM <- MAAT - Aa/2                               # Mean air temperature of the coldest month [degC]
    Lt <- (pi - 2 * asin(-MAAT/(Aa/2))) * P/(2 * pi)   # Length of the thawing season [d]
    Lf <- P - Lt                                       # Length of the freezing season [d]
    MATTS <- Ita/Lt                                    # Mean air temperature of the thawing season [degC]
    MATFS <- Ifa/Lf                                    # Mean air temperature of the freezing season [degC]

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

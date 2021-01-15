#####################################################################################################################
#
# Derivation of palaeo-air temperature characteristics from palaeo-active-layer thickness
# with the magnitude of annual air temperature oscillations defined by the mean air temperature of the warmest month,
# or rather its difference from the mean annual air temperature
#
#####################################################################################################################

MATWM <- function(xi, phi, rho, q, fc = c('fine','coarse'), nt, MATWM, showInputs = TRUE) {

  # Checking for physically feasible values of the input variables
  if(xi <= 0 | is.na(xi) == TRUE |                  # Palaeo-active-layer thickness [m]
     phi <= 0 | phi > 1 | is.na(phi) == TRUE |      # Volumetric ground moisture content [-]
     rho <= 0 | rho > 2700 | is.na(rho) == TRUE |   # Dry ground bulk density [kg/m3]
     q < 0 | q > 1 | is.na(q) == TRUE |             # Ground quartz content [-]
     fc != 'fine' & fc != 'coarse' |                # Ground grain-size class [-]
     nt <= 0 | is.na(nt) |                          # Ground-surface thawing n-factor [-]
     MATWM <= 0 | is.na(MATWM)) {                   # Mean air temperature of the warmest month [degC]

    return(NA)   # If any of the inputs is outside the range of physically feasible values or the grain-size class is not set as 'fine' or 'coarse' then the solution returns NA

  } else {

    # Calculation of thermal conductivity of thawed ground based on the Johansen's (1977) thermal-conductivity model [W/m/K]
    kt <- function(phi, rho, q, fc = c('fine','coarse')) {

      n <- 1 - rho / 2700   # Ground porosity as a function of dry bulk density and typical density of solids [-]
      S <- phi / n          # Degree of saturation as a proportion of volumetric ground moisture content and ground porosity [-]

      # Calculations specific for thermal conductivity of thawed fine-grained ground
      if(any(fc == 'fine')) {
        if(S <= 0.1 | S > 1) {
          return(NA)           # If the degree of saturation is outside the given range then the solution returns NA
        } else {
          Ke <- log10(S) + 1   # Kersten number for thawed fine-grained ground [-]
        }

      } else {

        # Calculations specific for thermal conductivity of thawed coarse-grained ground
        if(S <= 0.05 | S > 1) {
          return(NA)                 # If the degree of saturation is outside the given range then the solution returns NA
        } else {
          Ke <- 0.7 * log10(S) + 1   # Kersten number for thawed coarse-grained ground [-]
        }
      }

      kw <- 0.57                                            # Thermal conductivity of water [W/m/K]
      kq <- 7.7                                             # Thermal conductivity of quartz [W/m/K]
      ifelse(q < 0.2 & fc == 'coarse', ko <- 3, ko <- 2)    # Thermal conductivity of non-quartz minerals [W/m/K]
      ks <- kq^q * ko^(1 - q)                               # Thermal conductivity of solids as a function of thermal conductivity of quartz and non-quartz minerals [W/m/K]
      kdry <- (0.135 * rho + 64.7) / (2700 - 0.947 * rho)   # Semi-empirical equation for dry thermal conductivity of ground [W/m/K]
      ksat <- ks^(1 - n) * kw^n                             # Saturated thermal conductivity of thawed ground as a function of thermal conductivity of solids and water [W/m/K]

      kt <- kdry + (ksat - kdry) * Ke                       # Final calculation of thawed ground thermal conductivity [W/m/K]
      return(kt)
    }

    # Calculation of ground-surface (Its) and air (Ita) thawing index required to reach the given palaeo-active-layer thickness through an inverse solution of the Stefan equation
    It <- function(xi, kt, phi, nt) {
      Its <- (xi^2 * 334000 * phi * 1000) / (2 * kt * 86400)   # Ground-surface thawing index [degC.d]
      Ita <- Its/nt                                            # Air thawing index [degC.d]
      It <- list(Its, Ita)
      return(It)
    }

    # Calculation of mean annual air temperature (MAAT) based on the air thawing index (Ita) and the mean air temperature of the warmest month (MATWM)
    MAAT <- function(Ita, MATWM) {
      f4 <- function(T) {
        f3 <- function(T) {
          f2 <- function(T) {
            f1 <- function(t) {T + (MATWM - T) * sin(2 * pi * t / 365)}
          }
          integrate(f2(T), lower = asin(-T / (MATWM - T)) * 365 / (2 * pi), upper = (pi - asin(-T / (MATWM - T))) * 365 / (2 * pi))$value
        }
        Ita - f3(T)
      }
      # MAAT is searched in the interval [-100,0] (the lower limit of -100 substitutes negative infinity, which speeds up the solution); if the search fails then the solution returns NA
      MAAT <- tryCatch(uniroot(f4, lower = -100, upper = 0, tol = 0.001)$root, error=function(MAAT) NA)
      return(MAAT)
    }

    # Calculation of temperature characteristics
    kt <- kt(phi, rho, q, fc)                                        # Thermal conductivity of thawed ground [W/m/K]
    It <- as.numeric(It(xi, kt, phi, nt))                            # Thawing indices [degC.d]
    Its <- It[1]                                                     # Ground-surface thawing index [degC.d]
    Ita <- It[2]                                                     # Air thawing index [degC.d]
    MAAT <- MAAT(Ita, MATWM)                                         # Mean annual air temperature [degC]
    Ifa <- MAAT * 365 - Ita                                          # Air freezing index [degC.d]
    MATCM <- MAAT - (MATWM - MAAT)                                   # Mean air temperature of the coldest month [degC]
    Aa <- MATWM - MATCM                                              # Annual air temperature range [degC]
    Lt <- (pi - 2 * asin(-MAAT / (MATWM - MAAT))) * 365 / (2 * pi)   # Length of the thawing season [d]
    Lf <- 365 - Lt                                                   # Length of the freezing season [d]
    MATTS <- Ita / Lt                                                # Mean air temperature of the thawing season [degC]
    MATFS <- Ifa / Lf                                                # Mean air temperature of the freezing season [degC]

    # List of model outputs and inputs
    if(showInputs == TRUE) {   # Output includes the input variables (including calculated thermal conductivity of thawed ground)

      OutputsInputs <- c(MAAT, MATWM, MATCM, MATTS, MATFS, Ita, Ifa, Lt, Lf, Its, xi, phi, rho, q, fc, kt, nt, Aa)

    } else {   # Output excludes the input variables

      OutputsInputs <- c(MAAT, MATWM, MATCM, MATTS, MATFS, Ita, Ifa, Lt, Lf, Its)
    }

    if(any(is.na(OutputsInputs))) {

      return(NA)   # If any of the outputs or inputs is NA then the solution returns NA

    } else {

      return(OutputsInputs)
    }
  }
}

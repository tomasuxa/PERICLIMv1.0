## PERICLIMv1.0: A model deriving palaeo-air temperatures from thaw depth in past permafrost regions
PERICLIMv1.0 is a simple modelling scheme that allows to derive palaeo-air temperature characteristics associated with relict periglacial features indicative of the palaeo-active-layer thickness, and thus it is supposed to be used by periglacial geomorphologists interested in reconstructions of Quaternary palaeo-environments. It principally builds on an inverse solution of the Stefan equation, which is among the most commonly used analytical tools to estimate the thickness of the active layer. This way it deduces thawing-season temperature conditions, which it further converts into annual as well as freezing-season air temperature attributes based on the assumed annual air temperature range.

The package is a supplementary material to the below publication, which details the model, evaluates its performance against modern data, tests its viability to derive palaeo-air temperature characteristics, and discusses its uncertainties and applicability with respect to other palaeo-proxy records and/or model products, and thus the publication is advised to be read carefully before running the code:

Uxa, T., Křížek, M., and Hrbáček, F.: PERICLIMv1.0: a model deriving palaeo-air temperatures from thaw depth in past permafrost regions, Geosci. Model Dev., 14, 1865–1884, https://doi.org/10.5194/gmd-14-1865-2021, 2021.

Please cite the paper when using the package.

### Installation instructions
install.packages("devtools")

library(devtools)

install_github("tomasuxa/PERICLIMv1.0")

### PERICLIMv1.0 schemes
The package consists of two alternate schemes deriving past air temperature conditions that differ in how they define the magnitude of annual air temperature oscillations, but both provide identical outcomes if driven with compatible data.

`PERICLIMv1.0::Aa(xi, phi, rho, q, fc = c('fine','coarse'), nt, Aa, showInputs = TRUE)` utilizes the air temperature range (Aa), that is, the difference between the mean air temperature of the warmest and coldest month, to define the magnitude of annual air temperature oscillations and its input variables are as follows:

Variable | Symbol | Unit/Class
:-------- | :------ | :----
Palaeo-active-layer thickness | xi | m
Volumetric ground moisture content | phi | –
Dry ground bulk density | rho | kg/m3
Ground quartz content | q | –
Ground grain-size class | fc | –
Ground-surface thawing n-factor | nt | –
Annual air temperature range | Aa | degC
Should the inputs be included among the outputs (TRUE by default)? | showInputs | –

`PERICLIMv1.0::MATWM(xi, phi, rho, q, fc = c('fine','coarse'), nt, MATWM, showInputs = TRUE)` utilizes the mean air temperature of the warmest month (MATWM), or rather its difference from the mean annual air temperature, to define the magnitude of annual air temperature oscillations and its input variables are as follows:

Variable | Symbol | Unit/Class
:-------- | :------ | :----
Palaeo-active-layer thickness | xi | m
Volumetric ground moisture content | phi | –
Dry ground bulk density | rho | kg/m3
Ground quartz content | q | –
Ground grain-size class | fc | –
Ground-surface thawing n-factor | nt | –
Mean air temperature of the warmest month | MATWM | degC
Should the inputs be included among the outputs (TRUE by default)? | showInputs | –

Because both schemes calculate ten below-listed air temperature characteristics that differ only in that `PERICLIMv1.0::Aa(xi, phi, rho, q, fc = c('fine','coarse'), nt, Aa, showInputs = TRUE)` computes (among other attributes) the mean air temperature of the warmest month, whereas `PERICLIMv1.0::MATWM(xi, phi, rho, q, fc = c('fine','coarse'), nt, MATWM, showInputs = TRUE)` computes (among other attributes) the annual air temperature range, the outputs are, for convenience, provided in the same order for both schemes together with their inputs (if showInputs = TRUE) as follows:

Variable | Symbol | Unit
:-------- | :------ | :----
Mean annual air temperature | MAAT | degC
Mean air temperature of the warmest month | MATWM | degC
Mean air temperature of the coldest month | MATCM | degC
Mean air temperature of the thawing season | MATTS | degC
Mean air temperature of the freezing season | MATFS | degC
Air thawing index | Ita | degC.d
Air freezing index | Ifa | degC.d
Length of the thawing season | Lt | d 
Length of the freezing season | Lf | d
Ground-surface thawing index | Its | degC.d
Palaeo-active-layer thickness | xi | m
Volumetric ground moisture content | phi | –
Dry ground bulk density | rho | kg/m3
Ground quartz content | q | –
Ground grain-size class | fc | –
Thermal conductivity of thawed ground | kt | W/m/K
Ground-surface thawing n-factor | nt | –
Annual air temperature range | Aa | degC

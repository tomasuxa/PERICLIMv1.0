## PERICLIMv1.0: A model deriving palaeo-air temperatures from thaw depth in past permafrost regions
PERICLIMv1.0 is a simple modelling scheme that is designed to derive palaeo-air temperature characteristics associated with relict periglacial features indicative of the base of the palaeo-active layer, and thus it supposed to be used by periglacial geomorphologists interested in reconstructions of Quaternary palaeo-environments. It principally builds on an inverse solution of the Stefan equation, which is among the most commonly used analytical tools to estimate the thickness of the active layer. This way it deduces summer temperature conditions, which it further converts into annual and winter air temperature attributes through an annual air temperature amplitude or a mean air temperature of the warmest month assuming a sine air temperature curve.

The package is a supplementary material to the below publication, which details the model, evaluates its performance against modern data, and discusses its uncertainties and applicability, and thus it is advised to be read carefully before running the code:

Uxa, T., Křížek, M., and Hrbáček, F.: PERICLIMv1.0: A model deriving palaeo-air temperatures from thaw depth in past permafrost regions, Geosci. Model Dev. Discuss., Manuscript ID: gmd-2020-278, in review, 2020.

Please cite the paper when using the package.

### Installation instructions
install.packages("devtools")

library(devtools)

install_github("tomasuxa/PERICLIMv1.0")

### PERICLIMv1.0 schemes
The package consists of two alternate schemes to derive past air temperature conditions that differ in how they define the magnitude of annual air temperature oscillations, but both provide identical outcomes if driven with compatible data.

`PERICLIMv1.0::Aa(z, kt, vmc, nt, Aa, P, showInputs = TRUE)` utilizes the temperature amplitude (Aa) to define the magnitude of annual air temperature oscillations and its input parameters are as follows:

Variable | Symbol | Unit
:-------- | :------ | :----
Thaw depth | z | m
Thawed ground thermal conductivity | kt | W/m/K
Volumetric ground moisture content | vmc | –
Thawing n-factor | nt | –
Annual air temperature amplitude | Aa | °C
Period of air temperature oscillations | P | d
Should the inputs be included in the output (TRUE by default)? | showInputs | –

`PERICLIMv1.0::MATWM(z, kt, vmc, nt, MATWM, P, showInputs = TRUE)` utilizes the mean air temperature of the warmest month (MATWM) to define the magnitude of annual air temperature oscillations and its input parameters are as follows:

Variable | Symbol | Unit
:-------- | :------ | :----
Thaw depth | z | m
Thawed ground thermal conductivity | kt | W/m/K
Volumetric ground moisture content | vmc | –
Thawing n-factor | nt | –
Mean air temperature of the warmest month | MATWM | °C
Period of air temperature oscillations | P | d
Should the inputs be included in the output (TRUE by default)? | showInputs | –

Because both schemes calculate ten air temperature characteristics that differ only in that `PERICLIMv1.0::Aa(z, kt, vmc, nt, Aa, P, showInputs = TRUE)` computes the mean air temperature of the warmest month, whereas `PERICLIMv1.0::MATWM(z, kt, vmc, nt, MATWM, P, showInputs = TRUE)` computes the annual air temperature amplitude, the outputs are, for convenience, provided in the same order for both schemes together with their inputs (if showInputs = TRUE) as follows:

Variable | Symbol | Unit
:-------- | :------ | :----
Mean annual air temperature | MAAT | °C
Mean air temperature of the warmest month | MATWM | °C
Mean air temperature of the coldest month | MATCM | °C
Mean air temperature of the thawing season | MATTS | °C
Mean air temperature of the freezing season | MATFS | °C
Air thawing index | Ita | °C.d
Air freezing index | Ifa | °C.d
Length of the thawing season | Lt | d 
Length of the freezing season | Lf | d
Ground-surface thawing index | Its | °C.d
Thaw depth | z | m
Thawed ground thermal conductivity | kt | W/m/K
Volumetric ground moisture content | vmc | –
Thawing n-factor | nt | –
Annual air temperature amplitude | Aa | °C
Period of air temperature oscillations | P | d

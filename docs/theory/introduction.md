#  pywellfast calculation methods

## introduction



pyfastwell is dependent on two python libraries developed by TNO: pythermonomics, and pywellgeo, which do the techno-economic calculations, and defined the well architectures respectively.

pyfastwell is an extension of  these packages. It includes:
- Analytical Element Method (AEM) (Fokker et al., 2005, Egberts et al., 2013) to calculate the well flow performance of advanced well architectures for homogeneous reservoir properties
- coldfront tracking in the reservoir, based on finite difference solution
- effects of brine properties (heat capacity, viscosity, density) on well flow performance
- wellbore heat losses, frictional losses and thermosiphon effects 
- geothermal doublet techno-economic performance, including well flow results, coldfront evolution, and thermal performance
- evaluation of Skin effects to be adopted in doubletCalc1D to obtain similar results


pywelllfast has been  benchmarked against semi-analytical (and numerical) well models. It includes many tests and examples,
including examples for Monte Carlo, sensitivity simulations to analyse uncertainties in reservoir properties, and well architectures


##  Approach

The geothermal power which can be produced takes into account the production flow rate Q, 
and the cooling of the produced brine in the heat conversion process:

E = 𝜂 Q C ΔT

where: 

- E is the converted power [W]
- 𝜂 is the conversion efficiency [-], which depends on the heat conversion system, band here is assumed to be 1
- Q is the flow rate of the produced brine [m³/s]
- C is the volumetric heat capacity of the brine or circulation fluid [J m⁻³ K⁻¹]
- ∆T is the temperature difference between the producer and injector at the topside at the heat exchanger [K]


##  Doublet well architectures

The well architecture can be defined in a flexible way, using the pywellresult package.


##  Doublet flow performance

Achievable flow rates  are calculated with the [AEM model](aemmodel\aem_model.md) in pyfastwell. The AEM model calculates the well inflow performance based on the well architecture, reservoir properties, and fluid properties.
The coldfront evolution in the reservoir is calculated based on a finite difference solution, taking into account the flow rate, reservoir properties, and thermal properties of the fluid.
The [AEM model](aemmodel\aem_model.md)  section provides more details on the AEM implementation in pyfastwell.
and has been extensively benchmarked against semi-analytical and numerical models as described in the [AEM benchmark](aemmodel/aem_benchmark.md) section.

##  Doublet  coldfront evolution

The coldfront evolution in the reservoir is calculated based on a finite difference solution, taking into account the flow rate, reservoir properties, and thermal properties of the fluid.
the coldfront evolution is described in [coldfront evolution](coldfront.md).

## Geothermal Resources

The geothermal resource is considered as a permeable subsurface layer, 
which can be corresponding a permeable sedimentary rock or fractured basement or magmatic rock. The geothermal resource parameters 
include reservoir depth, thickness, temperature, brine compositions, and reservoir flow properties such as permeability and porosity. 


## Economic Model

- The economic model performs a NPV and LCOE calculation, incorporating a discounted cash flow approach:
  - CAPEX (Capital Expenditure) and OPEX (Operational Expenditure) are calculated based on the well design, reservoir properties
  - The LCOE is calculated as the ratio of the total discounted cost to the total discounted energy produced over the lifetime of the geothermal system.
  - The economic model allows for sensitivity analysis on key parameters such as reservoir properties, asscoiated flow rate, 
    and conversion efficiency.

## Key Performance indicators
- [key performance indicators](kpiperformance.md): The key performance indicators (KPIs) are calculated based on the technical performance and economic model, 
  providing insights into the feasibility and profitability of geothermal projects. These include:
  - net power  (power in MWth)
  - Leveleized Cost of Energy (LCOE in  €ct/kWh)
  - Net present value (NPV in million €)
    

## References
1. Egberts, P., Shatyrbayeva, I., Fokker, P.A., 2013. *Well inflow modelling for wells not aligned to a numerical grid. SPE 165986.*
2. Fokker, P. A., Verga, F., & Egberts, P. J., 2005 *New semianalytic technique to determine horizontal well productivity index in fractured reservoirs. SPE Reservoir Evaluation & Engineering, 8(02), 123-131.*
3. Van Wees et al., 2012. *Geothermal aquifer performance assessment for direct heat production–Methodology and application to Rotliegend aquifers. Netherlands Journal of Geosciences, 91(4), 651-665. https://doi.org/10.1017/S0016774600000433*
4. Mijnlieff, 2020. *Introduction to the geothermal play and reservoir geology of the Netherlands. Netherlands Journal of Geosciences.https://doi.org/10.1017/njg.2020.2*

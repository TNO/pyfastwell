# Analytical Element Method (AEM) model for geothermal reservoir simulation


The well productivity and injectivity is modelled by a semi-analytical model approach 
based on the Analytical Element Model (AEM) developed by Fokker et al., (2005) and Egberts et al. (2013). 
The model has been modified for geothermal wells, assuming a mass balance between produced and injected brines. 
The mass balance and pressure difference between the injector and producer results in a set of equations 
which is one row and column larger than the number of well segments involved in the 
Analytical element solution proposed. In the AEM models described here we assume a 
laterally infinite aquifer of constant thickness and homogeneous (but possibly anisotropic) permeability. 
In order to explain the AEM approach and in order compare/benchmark the 3 dimensional results with analytical solutions (i.e. Van Wees et al., 2012, DoubletCalc1D), we first introduce the AEM for the two-dimensional case and later extend to three dimensions.

## Two-Dimensional Solution

The steady state pressure conditions at distance r from the well radius is given by Thiem’s solution for an infinite aquifer. In the well, a correction for the skin is applied:

$$
p(r) = 
\left\\{ 
\begin{matrix}
\frac{\mu q}{2 \pi k H} \ln\frac{r_e}{r}                     & \text{(in reservoir)}\\
\frac{\mu q}{2 \pi k H} \left( \ln\frac{r_e}{r_w} + S\right) & \text{(in well)}
\end{matrix}
\right.
$$

The pressure differences can be superposed for N wells. When choosing a very large reference radius, we can write the superposing effects of pressure for a well i with additional surrounding wells j as:

$$
P_i = \frac{\mu}{2 \pi k H} \left( q_i\left[\ln\frac{r_e}{r_{wi}} + S_i\right]  + \sum_{\substack{j=1 \\ j \ne i}}^{N} q_j\ln\frac{r_e}{L_j}  \right)
$$

$r_{wi}$ is the radius of the well completion, and  $L_j$ is distance of well j to the well i. As the sum of flow rates is zero, this can be simplified to:



$$
P_i = \frac{\mu}{2 \pi k H} \left( \sum_{j=1}^{N} q_j\ln r_e  + q_i\left[\ln\frac{1}{r_{wi}} + S_i\right]  + \sum_{\substack{j=1 \\ j \ne i}}^{N} q_j\ln\frac{1}{L_j}  \right)=
\frac{\mu}{2 \pi k H} \left( q_i\left[\ln\frac{1}{r_{wi}} + S_i\right] + \sum_{\substack{j=1 \\ j \ne i}}^{N} q_j\ln\frac{1}{L_j}  \right)
$$


When the injectors have a pressure $P_i=P_{inj}$ and the producers are known to have the same pressure $P_P$, we have (with $P_i=0$ for the producers) the following system of equations to solve:

$$
\begin{bmatrix}
I_{11} & I_{12} & \cdots & I_{1N} & D_1    \\
I_{22} & I_{22} & \cdots & I_{2N} & D_2    \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
I_{N1} & I_{N2} & \cdots & I_{NN} & D_N    \\
1      & 1      & \cdots & 1      & 0
\end{bmatrix}
\begin{bmatrix}
q_1 \\
q_2 \\
\vdots \\
q_N \\
P_P
\end{bmatrix}=
\begin{bmatrix}
P_1 \\
P_2 \\
\vdots \\
P_N \\
0
\end{bmatrix}
$$

(eq. 4)

with


$$
I_{ij} = \left\\{
\begin{matrix}
\frac{\mu_i}{2 \pi k H} \left[ \ln\frac{1}{r_{wi}} + S_i \right] & (i=j) \\
\frac{\mu_i}{2 \pi k H}  \ln\frac{1}{r_{ij}} & (i\ne j)
\end{matrix}
\right.
$$

where $\mu_i$ is viscosity at well $i$; $r_{wi}$ is the well radius and $r_{ij}$ the distance between the wells i and j

$D_i$: 0 when the well is injector, -1 when producer

$P_i$: reference injection pressure (i.e. 1 bar) when the well is injector, 0 when producer.

Results from solving the system of equations become:

$q_i$: flow rate for the wells, which can be used to determine injectivity index  [q_i/P_i] and productivity index in [-q_i/P_i], 
dependent of the sign of .

$P_p$: resulting production pressure to maintain the steady state mass balance.

For a doublet configuration with equal viscosities, equal wellbore radius and equal skin, the system of equations reduces to the well-known steady state solution (e.g. Van Wees et al., 2012):

$$\Delta p = \frac{\mu q}{ \pi k H} \left[ \ln\frac{L}{r_w} + S\right]$$

Please note that the viscosity is constant for each row of the equations(cf. eq. 3). In order to include the effects of the viscosity contrasts, the viscosity for each row is set to the viscosity around the well under consideration.

## Three-Dimensional Solution

In three dimensions the system of equations is setup for N well line segments instead of N wells. 
The structure of the system of equations remains the same, 
and the system of equations and solution approach (excluding the last row and column) 
is explained in detail in the annex of Egberts et al. (2013). 
The M submatrix is a summation over line segment contribution of injectors and producers (k = 0) 
and multiple images below the base and above the top of the reservoir, marked by the indices k = -1 and k = 1. 
These image contributions enforce no flow boundary conditions at the top and base of the reservoir. 
The i and j denote pairs of control point location for each of the line segments (chosen at the middle) and is at least $r_w$ away from the well bore segment i.

$$M_{ij} = \sum_{k=-1}^{1} \frac{\mu_i}{2 \pi k H} (\ln\left( \frac{1}{r_{ijk}} +S_i) \right)$$

Such that the full system of equations becomes:


$$
\begin{bmatrix}
M_{ij} & \cdots & D_1 \\\\
\vdots & \ddots & \vdots \\\\
M_{ij} & \cdots & D_N \\\\
1 & \cdots & 0
\end{bmatrix}
\begin{bmatrix}
q_1 \\\\
q_2 \\\\
\vdots \\\\
q_N \\\\
P_p
\end{bmatrix}=
\begin{bmatrix}
P_1 \\\\
P_2 \\\\
\vdots \\\\
P_N \\\\
0
\end{bmatrix}
$$


The method of Egberts et al. (2013) allows to incorporate an anisotropic permeability field.
In addition additional images are included to enforce spatial constrains in flow (i.e. no flow barriers). 
The AEM also allows to be extended to incorporate frictional losses in the right hand side based on the solution for .
This requires an iterative solution and needs to model the exact well segment pressure. Please note that the solution is in 3D is based a single reference  well pressure  and  relative to the ambient pressure, and neglects potential vertical variation in pressure in the well bore compared to the ambient reservoir pressure. This is considered a valid  assumption for geothermal reservoirs.

## Benchmark

The outcomes of the AEM model have been benchmarked for deviated and horizontal wells against 
analytical models including Doubletcalc1D in [AEM benchmark](aem_benchmark.md)

## Using and comparing results in Doubletcalc1D

The AEM implementation in pyfastwell adopts the geothermal gradeint and  brine viscosity formulation as adopted 
in  [DoubletCalc1D](https://www.nlog.nl/en/tools), adopting Batzle and Wang viscosity

pyfastwell contains a function to calculate the equivalent Skins for the injector and 
producer wells in [DoubletCalc1D](https://www.nlog.nl/en/tools). These have been validated
in the [AEM benchmark](aem_benchmark.md)

## References

1. Fokker, P. A., Verga, F., & Egberts, P. J., 2005 *New semianalytic technique to determine horizontal well productivity index in fractured reservoirs. SPE Reservoir Evaluation & Engineering, 8(02), 123-131.*
2. Egberts, P., Shatyrbayeva, I., Fokker, P.A., 2013. *Well inflow modelling for wells not aligned to a numerical grid. SPE 165986.*
3. Van Wees, J. D., Kronimus, A., Van Putten, M., Pluymaekers, M., Mijnlieff, H., Van Hooff, P., Obdam, A., & Kramers, L., 2012. *Geothermal aquifer performance assessment 
for direct heat production–Methodology and application to Rotliegend aquifers. Netherlands Journal of Geosciences 91 (04): 651-665.*

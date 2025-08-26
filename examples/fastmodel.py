"""
Wellcostmodel computes multilateral well costs 
based on well trajectory and well construction
"""

import pyfastwell as pwf
from pythermonomics.geothermal_economics import GeothermalEconomics

import sys
sys.setrecursionlimit(3000)

settingfile = 'input/npv_thermogis2024_138.yml'

#trajectoryfile = 'input/inputsMultilateral3legsDetailed.yml'
trajectoryfile = 'input/inputsStandardDetailed.yml'
#trajectoryfile = 'input/inputsDetailedTNOhor.yml'
#trajectoryfile ='input/inputXYZspline.yml'
#trajectoryfile ='input/inputXYZ3legs.yml'
#trajectoryfile ='input/inputXYZvertical.yml'

economics = GeothermalEconomics.from_trajectory (settingfile, trajectoryfile)
npv, lcoe, cashflow, simdata, wells, well_results = economics.compute_economics()
print(f"Well cost {cashflow['wellcosts'].iloc[0]/1e6} MEUR")


fastmodel = pwf.WellFastModel(economics, k=100, khkv=3, segmentlength=50, verbose=False)
fastmodel.set_dP_BHP(60)
npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
print(f"NPV {npv}, LCOE {lcoe_val}, COP {cop}, Power {power} MW, Treservoir {fastmodel.production_temperature} C, Tsurface {Tsurface} C, DT_eff {DT_eff} C, DP_eff {DP_eff} bar")
print(f"Well cost {cashflow['wellcosts'].iloc[0]/1e6} MEUR")

for i,w in enumerate(fastmodel.wells):
    print(f"Well cumulative AHD {w.cumulative_ahd()}")

L= 1300
dia_inch = 8
rw = dia_inch * 0.0254  * 0.5
skin_inj, skin_prd, ratio = fastmodel.getSkinFactors_dc1d(L, rw)
print(f"for use in DC1D with L ={L} [m]  and diameter={dia_inch} [inch] Skin factors: inj {skin_inj}, prd {skin_prd}, ratio {ratio}")

skin = fastmodel.skin_from_deviation(fastmodel.H,  rw, fastmodel.khkv, deviation=45.0, type='cinco')
print(f"for H {fastmodel.H}, khkv {fastmodel.khkv}, diameter={dia_inch} [inch] Skin factors: inj {skin}")

fastmodel.run_coldfront(simyears=30, rockdens=2700.0, rockcap=1000.0, porosity=0.21)
fastmodel.plot_coldfront()


fastmodel.plot_trajectories()



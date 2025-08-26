import unittest

from pywellgeo.well_data.dc1dwell import *
from pyfastwell.wellflow.well_fastmodel import WellFastModel
from pythermonomics.geothermal_economics import GeothermalEconomics
import sys
sys.setrecursionlimit(3000)

class MyTestCase(unittest.TestCase):


    def test_trajectory_types(self):
        settingfile = './input/npv_thermogis2024_138.yml'
        trajectory_files = [
                    './input/inputsMultilateral3legsDetailed.yml',
                 './input/inputsDetailedTNOhor.yml',
                  './input/inputXYZspline.yml',
                  './input/inputXYZ3legs.yml'
        ]
        outdir = 'output/fastmodel/types'

        for trajectoryfile in trajectory_files:
            print(f"Testing trajectory file: {trajectoryfile}")
            economics = GeothermalEconomics.from_trajectory(settingfile, trajectoryfile)
            npv, lcoe_val, cashflow, *_ = economics.compute_economics()
            fastmodel = WellFastModel(economics, k=100, khkv=1, segmentlength=50)
            fastmodel.set_flowrate(400)
            print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {economics.cop}, Power: {economics.power}")
            npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
            fastmodel.plot_trajectories(outdir=outdir)
            fastmodel.run_coldfront(simyears=30, rockdens=2700.0, rockcap=1000.0, porosity=0.21,  ncell=200)
            fastmodel.plot_coldfront(outdir=outdir)


    def test_trajectory_flowrates(self):
        trajectoryfile = './input/inputsDetailedTNOhor.yml'
        settingfile = './input/npv_thermogis2024_138.yml'
        outdir = 'output/fastmodel/flowrates'
        economics = GeothermalEconomics.from_trajectory (settingfile, trajectoryfile)
        fastmodel = WellFastModel(economics,k=100, khkv=1, segmentlength=50)

        print()

        fastmodel.set_flowrate(400)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")

        fastmodel.set_flowrate(100)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff  = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")

        fastmodel.set_flowrate(200)
        fastmodel.set_k_khkv(50, 1.0, tryscale=True)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff  = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")
        fastmodel.plot_trajectories(outdir=outdir)

        fastmodel.run_coldfront(simyears=30, rockdens=2700.0, rockcap=1000.0, porosity=0.21)
        fastmodel.plot_coldfront(outdir=outdir)

    def test_inputXYZvertical(self):
        trajectoryfile = './input/inputXYZvertical.yml'
        settingfile = './input/npv_thermogis2024_138.yml'
        outdir = 'output/fastmodel/inputXYZvertical'
        economics = GeothermalEconomics.from_trajectory(settingfile, trajectoryfile)
        npv, lcoe_val, cashflow, *_ = economics.compute_economics()
        fastmodel = WellFastModel(economics, k=100, khkv=1, segmentlength=50)
        fastmodel.set_flowrate(106)
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {economics.cop}, Power: {economics.power}")
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
        fastmodel.plot_trajectories(outdir=outdir)
        fastmodel.run_coldfront(simyears=30, rockdens=2700.0, rockcap=1000.0, porosity=0.21, ncell=200, maxlength=5)
        fastmodel.plot_coldfront(outdir=outdir) #, title =f'Cold Front Vertical: k=100, khkv=1,  Q={fastmodel.flowrate} m3/h, DP= {round(DP_eff)} bar')

    def test_DC1D_trajectory(self):
        dc1dsettings = './input/dc1dwell_skin0.yml'
        settingfile = './input/npv_thermogis2024_138.yml'
        outdir =  'output/fastmodel/dc1d'
        dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        dc1dwell.qvol =   -125/3600
        dc1dwell.dp = 30
        dc1dwell.calculateDP_qvol()
        economics = GeothermalEconomics.from_dc1d(settingfile, dc1dsettings, dc1dwell)
        # reference economics call (not using pyfastwell)
        npv, lcoe_val, cashflow, *_ = economics.compute_economics()
        cop = economics.cop
        power = economics.power
        flowrate = economics.simresults.res_params.flowrate
        print(f"flowrate: {flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")

        fastmodel = WellFastModel(economics, k=500, khkv=1, segmentlength=50)
        # check if pyfastwell solution gives same result with flowrate set to same as in economics
        fastmodel.set_flowrate(136)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff  = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")


        # check if energy losses are appied correctly by testing for setting pressure BHP difference
        # corrected for pressure drop in well trajector (dpsum)
        dpsum = fastmodel.economics.welltrajectory.trajectoryinput.dpsum
        fastmodel.set_dP_BHP(30-dpsum[0])
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff  = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")


        # check if pyfastwell solution gives same result with khkv changed (should not change result)
        fastmodel.set_k_khkv(500, 10.0, tryscale=True)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
        print(f"flowrate: {fastmodel.flowrate}, NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power}")
        fastmodel.plot_trajectories(outdir=outdir)

        fastmodel.run_coldfront(simyears=100, rockdens=2700.0, rockcap=1000.0, porosity=0.21, maxlength=5, ncell=100)
        fastmodel.plot_coldfront(outdir=outdir)






if __name__ == '__main__':
    unittest.main()
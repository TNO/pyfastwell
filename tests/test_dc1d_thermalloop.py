import os
import unittest

import numpy as np

from pywellgeo.well_data.dc1dwell import Dc1dwell
from pythermonomics.geothermal_economics import GeothermalEconomics
import pyfastwell.plot.plotwells as plotwells



class MyTestCase(unittest.TestCase):
    settingfile = './input/npv_thermogis2024.yml'

    def test_DC1D(self):
        """
        this tests the DC1D well trajectory and economics calculation
        it includes the calculation of the pressure drop and the volume flow rate based
        on an imposed pressure drop

        it generates the temperature profile along the well, based on the thermal gradient and surface
        temperature.

        it also includes the calculation of the heat loss along the well
        and the calculation of the production temperature at surface

        the heat loss (ca 1.15 C) is lower than the DC1D values (close to 1.4C). This is
        because the DC1D values are based on AHD linearly correlated with the TVD, whereas
        the heat loss in this test is based on the actual temperature values related to TVD of the wellpath

        """

        dc1dsettings = './input/dc1dwell.yml'
        dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        dc1dwell.qvol = -1
        dc1dwell.dp = 30
        g = dc1dwell.get_params()
        g["ahd"] = g["tvd"] + (np.array(g["ahd"]) - np.array(g["tvd"]))
        dc1dwell.update_params(**g)
        dc1dwell.calculateDP_qvol()

        gt_economics = GeothermalEconomics.from_dc1d(self.settingfile, dc1dsettings, dc1dwell)
        npv, lcoe, cashflow, *_ = gt_economics.compute_economics()

        inj, prod =  gt_economics.welltrajectory.tw['INJ1']['welltree'], gt_economics.welltrajectory.tw['PRD1']['welltree']
        ahd = inj.cumulative_ahd()
        fname = os.path.basename(dc1dsettings)

        outputnamebase = os.path.join('output/benchmark_aem/thermalloop', str(fname).split('.')[0])
        plotwells.plot_wells( [inj, prod],outputnamebase )
        #economics.plot_wells_panels(fig=fig, ax=ax,doplot=(i == len(wells) - 1),tofile=outputbasename + '3d.png')
        fnameout = outputnamebase + "_cashflow.csv"
        with open(fnameout, 'w') as f:
            cashflow.to_csv(f, lineterminator='\n', index=True)

        self.assertAlmostEqual(gt_economics.power, 3.68, delta=0.01)
        self.assertAlmostEqual(gt_economics.cop, 21.3, delta=0.1)
        self.assertAlmostEqual(dc1dwell.qvol, 125 / 3600, delta=1e-3)
        self.assertAlmostEqual(dc1dwell.dp, 30, delta=1e-1)
        self.assertAlmostEqual(ahd, 2140, delta=1.0)
        self.assertAlmostEqual(lcoe, 9.96, delta=0.1)
        self.assertAlmostEqual(npv, -6.09e6, delta=1e5)




    def test_DC1D_fixedheatloss(self):
        """
        this tests the DC1D well trajectory and economics calculation
        it includes the calculation of the pressure drop and the volume flow rate based
        on an imposed pressure drop

        the heat loss is assumed to be zero, so the temperature is constant along the well
        and the production temperature at surface (and reservoir) is equal to temperature
        specified in the input file.

        """
        dc1dsettings = './input/dc1dwell_fixedheatloss.yml'
        dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        dc1dwell.qvol = -1 # negative means use dp instead
        dc1dwell.dp = 30
        g = dc1dwell.get_params()
        ahd_min_tvd = np.array(g['ahd']) - np.array(g['tvd'])

        #g['tvd'] = [ 3000,3000]
        g['ahd'] = g['tvd'] + ahd_min_tvd
        dc1dwell.update_params(**g)

        dc1dwell.calculateDP_qvol()

        economics= GeothermalEconomics.from_dc1d(self.settingfile, dc1dsettings, dc1dwell)
        npv, lcoe, cashflow, *_ = economics.compute_economics()
        #gt_lcoe = Geothermal_Lcoe.from_dc1d(self.settingfile, dc1dsettings, dc1dwell)
        #npv, lcoe, cashflow, simdata, wells = gt_lcoe.compute_npv()
        #ahd = gt_lcoe.welltrajectory.tw['INJ1']['welltree'].cumulative_ahd()

        inj, prod =  economics.welltrajectory.tw['INJ1']['welltree'], economics.welltrajectory.tw['PRD1']['welltree']
        ahd = inj.cumulative_ahd()
        fname = os.path.basename(dc1dsettings)
        outputnamebase = os.path.join('output/benchmark_aem/thermalloop', str(fname).split('.')[0])

        plotwells.plot_wells( [inj, prod],outputnamebase )
        #economics.plot_wells_panels(fig=fig, ax=ax,doplot=(i == len(wells) - 1),tofile=outputbasename + '3d.png')
        fnameout = outputnamebase + "_cashflow.csv"
        with open(fnameout, 'w') as f:
            cashflow.to_csv(f, lineterminator='\n', index=True)
        cop = economics.cop
        power = economics.power


        self.assertAlmostEqual(dc1dwell.qvol, 125/3600, delta=1e-3)
        self.assertAlmostEqual(dc1dwell.dp, 30, delta=1e-3)
        self.assertAlmostEqual(dc1dwell.tprod, 56.14, delta=1e-3)
        self.assertAlmostEqual(ahd, 2140, delta=1.0)
        self.assertAlmostEqual(power, 3.61, delta=1e-2)
        self.assertAlmostEqual(cop, 20.8, delta=0.25)
        self.assertAlmostEqual(lcoe, 9.96, delta=0.1)
        self.assertAlmostEqual(npv, -6.09e6, delta=6e4)


if __name__ == '__main__':
    unittest.main()

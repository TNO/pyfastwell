"""
benchmark test for the AEM well model.
It checks the AEM well model against the Thiem solution and the DC1D solution

for two verticals wells in a homogeneous reservoir seperated at 1300 m distance.
The parameters are chosen such that they agree with the example2300.xml of the DC1D model with a skin factor of 0.0 (for vertical wells)
and n applied pressure difference of 60 bar for the loop resulting in a flowrate of ca 106 m3/h.

Please note that the DC1D input uses a well radius of 8 inch (0.1018m), corresponding to slightly higher radius than the 0.1 m
used here.

In adddition the tests plots the well trajectories and the pressure distribution along a line between the two wells
for the Thiem solution, the DC1D solution and the AEM solution.

"""

import matplotlib.pyplot as plt
import pywellgeo.well_data.water_properties as watprop
from pywellgeo.well_data.dc1dwell import *
from pywellgeo.well_data.names_constants import Constants
from pywellgeo.well_tree.well_tree_tno import WellTreeTNO
from pyfastwell.wellflow.wi_tno import WiTNO
from pyfastwell.plot.plotwells import plot_wells
from pyfastwell.wellflow.well_fastmodel import WellFastModel
from pythermonomics.geothermal_economics import GeothermalEconomics
import numpy as np
import sys
sys.setrecursionlimit(3000)
import unittest
import os

def thiems(x, xw, rw, qw, mu=1e-3, h=100, k=100):
    r = abs(x-xw)
    rmin = rw + r*0
    r = np.maximum(r,rmin)
    kH = k*h* Constants.DARCY * 1e-3  # convert to mD*m
    dp = qw* mu * np.log(r/rw)/ (2 * np.pi * kH)
    return dp/1e5




class MyTestCase(unittest.TestCase):

    def test_benchmark_aem(self):

        well1 = WellTreeTNO.from_vertical(-650, 0, -2450, 0.1)
        well2 = WellTreeTNO.from_vertical(650,0,-2450, 0.1)
        wells = [well1, well2]
        production_temperature = 2350*0.031 + 10
        injection_temperature  = 50
        salinity = 140000
        flowrate = 106
        topres = [2300, 2300]
        botres = [2400, 2400]
        print('AHD INJ1 ', wells[0].cumulative_ahd())
        print('AHD PROD1 ', wells[1].cumulative_ahd())

        muprod = watprop.viscosity( production_temperature, salinity * 1e-6)
        muinj = watprop.viscosity(injection_temperature, salinity * 1e-6)

        k = 100
        H = 100
        kxx, kyy, kzz = k, k, k*0.1
        wi = WiTNO(wells[0], wells[1], -topres[0], -botres[0], kxx, kyy,kzz, muinj=muinj, muprod=muprod,  segmentlength=10)
        res,ii, pi = wi.setupmatrix(n=0)

        print ('q/m segments res ', res)
        print ('II m3/s /bar ', ii)
        print ('PI m3/s /bar ', pi)



        outdir = 'output/benchmark_aem'
        plotname = outdir + '/verticalwellWI'
        wells = [well1, well2]
        plot_wells(wells, plotname)

        xw = np.asarray([ -650, 650])
        q = flowrate/3600
        qw = np.asarray([q, -q])
        xmu = np.asarray([ muinj, muprod])
        x = np.arange(-1000,1000,1)
        #x= np.asarray(xw)
        dp = x*0
        rw = xw * 0
        rw = rw +0.1


        for i,xwell in enumerate(xw):
                dp = dp + thiems(x, xw[i], rw[i], qw[i], mu=np.flip(xmu)[i], h=100, k=100)

        dpw_dc1d = xw*0.0
        L = xw[1]-xw[0]
        for i, xwell in enumerate(xw):
            dpw_dc1d[i] = qw[i] / WellFastModel.pi_dc1d(L, rw[i], k*H, xmu[i], 0.0)

        ii_pi = [ ii, pi]
        dpw_aem = xw* 0.0
        for i, xwell in enumerate(xw):
            dpw_aem[i] = qw[i]/ ii_pi[i]

        self.assertAlmostEqual(dpw_aem[0], 35.9, delta=1)
        self.assertAlmostEqual(-dpw_aem[1], 24.7, delta=1)
        self.assertAlmostEqual(dpw_dc1d[0], 35.9, delta=1)
        self.assertAlmostEqual(-dpw_dc1d[1], 24.7, delta=1)

        fig, axes = plt.subplots(1, 1, figsize=(10,10))
        plt.xlabel('distance [m]')
        plt.ylabel('dp [bar]')
        plt.title (f"kH: 10 Dm, L: 1300, rw: 0.1, flowrate: {flowrate} m3/h, DP: {dpw_aem[0]-dpw_aem[1]} bar")
        plt.plot(x, -dp,  color='red', marker='x', label='Thiem solution')
        plt.scatter(xw, dpw_dc1d, color='green', marker='o', facecolors='none', label='DC1D solution')
        plt.scatter(xw, dpw_aem, color='magenta',  marker='o', facecolors='none', label='AEM solution')
        plt.legend()

        plotname = outdir + '/benchmark_DC1D_AEM_THIEMS'
        plt.savefig(plotname + '.png', dpi=300)
        plt.close()

    def test_DC1D(self):
        """
        test if pyfastwell gives same results as dc1dwell for a simple vertical well with skin factor 0.
        It is the same case as in the dc1d input file example.xml (modified with skin factor 0)

        :return:
        """
        dc1dsettings = './input/dc1dwell_skin0.yml'
        settingfile = './input/npv_thermogis2024_138.yml'
        dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        dc1dwell.qvol =   -125/3600
        dc1dwell.dp = 30
        dc1dwell.calculateDP_qvol()
        dpres = dc1dwell.dpres

        economics = GeothermalEconomics.from_dc1d(settingfile, dc1dsettings, dc1dwell)
        # reference economics call (not using pyfastwell)

        fastmodel = WellFastModel(economics, k=500, khkv=1, segmentlength=50)
        # check if pyfastwell solution gives same result with flowrate set to same as in economics
        fastmodel.set_dP_BHP(dpres)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff  = fastmodel.compute_economics()
        # check for
        self.assertAlmostEqual(Tsurface ,56.5, delta=1e-2*56.5)
        self.assertAlmostEqual(DT_eff, 26.5, delta=1e-2*26.5)
        self.assertAlmostEqual( DP_eff, 30, delta=1e-2*30)
        self.assertAlmostEqual(fastmodel.flowrate, dc1dwell.qvol*3600, delta =1e-2*dc1dwell.qvol*3600)
        self.assertAlmostEqual(power, 4.03, delta=1e-2*4.03)





    def test_skins_varyingdip(self):
        """
        test varying skin factor calculations for deviated wells
        they are compared to the analytical solutions from:

        1. CH. Cinco, F.G. Miller, H.J. Ramey Jr., 1975. *Unsteady-state pressure distribution created by a directionally drilled well. J. Petrol. Technol., 27 (11) (1975), pp. 1-392*
        2. E.J. Rogers, M.J. Economides, 1996.  *The skin due to slant of deviated wells in permeability-anisotropic reservoirs International Conference on Horizontal Well Technology, Society of Petroleum Engineers*
        3. Besson H., 1990. *Performance of slanted and horizontal wells on an anisotropic medium. European Petroleum Conference, Society of Petroleum Engineers (1990)*


        """
        import numpy as np
        deviations = np.array([0, 10, 20, 30, 40, 50, 60, 70])

        khkvs = [1,3,5,10]


        xoff = 650
        L_dc1d = 2* xoff
        topdepth = 2300
        H = 100  # m
        botdepth = topdepth + H
        k =100 # mD
        rw = 0.1

        for i, khkv in enumerate( khkvs):
            skin_fastmodel = []
            skin_rogers = []
            skin_cinco = []
            skin_besson = []
            for deviation in deviations:
                output_dir = 'output/benchmark_aem/skin_deviation/deviation' + str(deviation)
                trajectoryfile = os.path.join(output_dir, 'inputXYZ.yml')
                if (i==0):
                    # create well trajectories with varying deviation
                    perforation_L = H/np.cos(np.radians(deviation))
                    perforation_X = perforation_L * np.sin(np.radians(deviation))
                    xi1,xi2 = xoff -0.5* perforation_X, xoff +0.5* perforation_X
                    xp1,xp2 = -xi1, -xi2
                    xyz_inj = np.array([
                        [0, 0, 0],
                        [0, 0, 0.5*topdepth],
                        [xi1, 0,topdepth],
                        [xi2, 0, botdepth],
                    ])
                    xyz_prd = np.array([
                        [0, 0, 0],
                        [0, 0, 0.5*topdepth],
                        [xp1, 0, topdepth],
                        [xp2, 0, botdepth],
                    ])
                    data = {
                        # "format": DoubleQuotedScalarString("XYZGENERIC"),
                        "format": 'XYZGENERIC',
                        "reservoir": {
                            "basic": {
                                "top_reservoir_depth_TVD": topdepth,
                                "bottom_reservoir_depth_TVD": botdepth
                            }
                        },
                        "well_trajectories": {
                            "INJ1": {
                                "main_wellbore": {
                                    "xyz": xyz_inj.tolist(),
                                    "radius": rw,
                                    "mindist": 25
                                }
                            },
                            "PRD1": {
                                "main_wellbore": {
                                    "xyz": xyz_prd.tolist(),
                                    "radius": rw,
                                    "mindist": 25
                                }
                            }
                        }
                    }
                    os.makedirs(output_dir, exist_ok=True)
                    with open(trajectoryfile, 'w') as f:
                        yaml.dump(data, f, default_flow_style=False)

                minimal_settings = './input/npv_thermogis2024.yml'
                economics = GeothermalEconomics.from_trajectory(minimal_settings, trajectoryfile=trajectoryfile)
                npv, lcoe_val, cashflow, *_ = economics.compute_economics()
                fastmodel = WellFastModel(economics,k=k, khkv=khkv, segmentlength=25)
                npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power,  Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
                fastmodel.plot_trajectories(outdir=output_dir)
                if deviation>0:
                   skin_cinco.append (fastmodel.skin_from_deviation(H, rw, khkv, deviation=deviation, type='cinco'))
                   skin_rogers.append(fastmodel.skin_from_deviation(H, rw, khkv, deviation=deviation, type='rogers'))
                   skin_besson.append(fastmodel.skin_from_deviation(H, rw, khkv, deviation=deviation, type='besson'))
                   skin_fastmodel.append(fastmodel.getSkinFactors_dc1d(L_dc1d,rw)[1])
                else:
                    skin_cinco.append(0.0)
                    skin_rogers.append(0.0)
                    skin_besson.append(0.0)
                    skin_fastmodel.append(fastmodel.getSkinFactors_dc1d(L_dc1d,rw)[1])
            plt.plot(deviations, skin_cinco, label='Cinco et al.,1975', marker='o',markerfacecolor='none')
            plt.plot(deviations, skin_rogers, label='Rogers and Economides, 1996', marker='s',markerfacecolor='none')
            plt.plot(deviations, skin_besson, label='Besson, 1990', marker='^',markerfacecolor='none')
            plt.plot(deviations, skin_fastmodel, label='pyfastwell-AEM')

            plt.xlabel('Deviation (degrees)')
            plt.ylabel('Skin Factor')
            plt.title('Skin Factor vs. Deviation for khkv=' + str(khkv))
            plt.legend()
            plt.grid(True)
            plt.savefig ('output/benchmark_aem/skin_deviation/skin_vs_deviation'+ str(khkv)+ '.png')
            plt.close()

if __name__ == '__main__':
    unittest.main()
import unittest

from pywellgeo.well_data.dc1dwell import *



class MyTestCase(unittest.TestCase):
    def test_dc1dwell(self):
        """
        Simple test case for Doubletcalc1d with two vertical wells
        The parameters are taken from the example.xml input file of Doubletcalc1D V1.5 (which is available in the docs/source)
        it checks if the calculated DP is close to the expected value of 30 bar

        :return:
        """
        k = 500 # mDarcy
        H = 70 #thickness m
        L = 1700 # distance between the wells m
        tvd = [1500, 1500] # true vertical depth of the wells m
        temp = [30, 10 + (tvd[0]+0.5*H)*0.031 ]
        salinity = [70000, 70000]
        skin = [0.5, 2.0]
        ahd = [ 2018, 2121]
        rwm  = 0.5 * 7 * 2.54/100 # inch to m
        rw = [ rwm, rwm]
        roughness = 1.38 # millinch
        dc1d = Dc1dwell(k, H, L, tvd, temp, salinity, skin, ahd, rw, roughness)
        qvol = 125/3600 # m3/s
        dc1d.calculateDP(qvol)
        print('DP of full loop [bar]', dc1d.dp)
        self.assertAlmostEqual(dc1d.dp, 30, delta=0.05)

    def test_dc1dwell2(self):
        """
        Simple test case for Doubletcalc1d with two vertical wells, same as test_dc1dwell but now reading the parameters from a config file

        :return:
        """
        configfile = "./input/dc1dwell.yml"
        dc1d = Dc1dwell.from_configfile(configfile)
        qvol = 125 / 3600
        dc1d.calculateDP(qvol)
        print('DP of full loop [bar]', dc1d.dp)
        self.assertAlmostEqual(dc1d.dp, 30, delta=0.05)

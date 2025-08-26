from datetime import datetime

import numpy as np

from pywellgeo.well_data.dc1dwell import Dc1dwell
from pythermonomics.geothermal_economics import GeothermalEconomics

from pyfastwell.stochastic.stochasticmodel import Stochasticmodel


ITIME = 0
IPOWER=1
ICOP =2
INPV =3
ILCOH=4


"""
this is a stochastic model wrapper around the dc1dwell model and the geothermal economics model
it can be used to generate ensembles of npv and cop values based on uncertain parameters

"""

class   Dc1d_lg_model(object):
    """
    StochasticWellModel class
    """

    def __init__(self, settingfile, dc1dsettings, locked_param=None, free_dist=None):
        """
        set the dc1dwell object
        :param dc1dsettings: file name yml with dc1dwell settings
        :param locked_param: dictionary with added or singular value parameters (defaulkt is none and considers the dedault from dc1dsettings
        :param free_dist: dictionary with the free parameters and their distributions
        """
        self.settingfile = settingfile
        self.dc1dsettings = dc1dsettings
        self.dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        self.ahd_min_tvd = np.asarray(self.dc1dwell.ahd) - np.asarray(self.dc1dwell.tvd)
        g = self.dc1dwell.get_params()
        if locked_param is not None:
            g.update(locked_param)  # add singular value parameters to the dictionary
        free_param = list(free_dist.keys())
        argdict = {Stochasticmodel.LOCKED: g, Stochasticmodel.FREEPARAM: free_param,
                        Stochasticmodel.DISTPARAM: free_dist}
        """
        Initialize the StochasticWellModel object
        :param argdict: dictionary with keys 'locked', 'free_param', and 'free_dist' containing the model parameters
        """

        self.stoch = Stochasticmodel(argdict=argdict)

    def forward(self, param):  # Forward model to which the parameters need to be fitted
        """
        :param  param: np.ndarray, including
        :param argdict: dictionary with the following keys in the dictionary:
            which can be used to update the parameters of the forward model in the param (if is included as distribition)
            k, H, topdepth, DP, return_temp, Q

        :return: forward model results in dimension of locked parameters 't'
        """

        t = self.stoch.get_from_dict('t', param)


        k = self.stoch.get_from_dict('k', param)
        H  = self.stoch.get_from_dict('H', param)
        topdepth = self.stoch.get_from_dict('topdepth', param)
        DP = self.stoch.get_from_dict('DP', param)
        return_temp = self.stoch.get_from_dict('return_temp', param)
        qvol = self.stoch.get_from_dict('qvol', param)
        tvd = np.array([topdepth, topdepth])
        ahd = tvd + self.ahd_min_tvd

        temp = np.array(self.dc1dwell.temp) * 1.0
        temp[0] = return_temp

        g = {'k': k, 'H': H, 'tvd': tvd, 'ahd': ahd, 'dp': DP, 'temp': temp, 'qvol': qvol}

        self.dc1dwell.update_params(**g)

        self.dc1dwell.calculateDP_qvol()

        dc1d = self.dc1dwell
        economics= GeothermalEconomics.from_dc1d(self.settingfile, self.dc1dsettings, dc1d)
        npv, lcoe, cashflow, *_ = economics.compute_economics()

        #economics= Geothermal_Lcoe.from_dc1d(self.settingfile, self.dc1dsettings, self.dc1dwell)
        #npv, lcoe, cashflow, simdata, wells = economics.compute_npv()

        cop = economics.cop
        power = economics.power
        #ahd = gt_lcoe.welltrajectory.tw['INJ1']['welltree'].cumulative_ahd()

        fwd = np.vstack((t, [power],[cop], [npv*1e-6], [lcoe]))
        return fwd

    def generate_ensemble(self, nsamples):
        """
        generate the ensemble of parameters
        :param nsamples: number of samples
        :return: ensemble of parameters
        """
        return self.stoch.generate_ensemble(nsamples)

    def run_ensemble(self, M):
        """
        generate the ensemble of results
        :param M: list of parameters
        :return: ensemble of parameters
        """
        return self.stoch.run_ensemble(M, self.forward)

def main():

    np.random.seed(12345)

    print("Testing Dc1dmodel ensemble running code")
    stime = datetime.now()
    print ('start at: ', stime.time())

    nsamples = 500  # number of values in a prior, for one parameter

    t = np.arange(1)  # basically one time only

    dc1dsettings = 'input/dc1dwell.yml'
    free_dist = {'k': {'dist': 'triangular', 'values': [350, 650, 500]}, \
                 'H': {'dist': 'triangular', 'values': [50,90,70]}
                 }
    #free_dist = {'k': {'dist': 'uniform', 'values': [490, 510]}, \
    #             'H': {'dist': 'uniform', 'values': [65,75]}
    #             }
    qvol = -125/3600
    locked_param = { 't': t, 'topdepth': 1500, 'DP': 30, 'return_temp': 30, 'qvol':qvol }

    settingfile = 'input/npv_thermogis2024_138.yml'
    dc1dmodel = Dc1d_lg_model(settingfile, dc1dsettings, locked_param=locked_param, free_dist=free_dist)

    # get median values for the parameters and test the forward function
    funcparam = dc1dmodel.stoch.getmedianparams()
    res = dc1dmodel.forward(funcparam)

    M = dc1dmodel.generate_ensemble(nsamples)
    fwi = np.array(dc1dmodel.run_ensemble(M))

    etime = datetime.now()
    print ('end at: ', etime.time())
    print ('duration: ', etime - stime)


    # first index of the output array is the index of the realizations  in the model prediction ensemble, typically all :
    # second index is the output index of the model predictions, typically INPV or ICOP
    # third index is the time index in the model predictions, ITIME
    pvals = [10,30,50,70,90]
    dc1dmodel.stoch.expectation_plot(fwi[:, IPOWER, ITIME], 'POWER [MW] ', pvals=pvals)
    dc1dmodel.stoch.expectation_plot(fwi[:, ICOP, ITIME], 'COP [-]', pvals=pvals)
    dc1dmodel.stoch.expectation_plot(fwi[:, INPV, ITIME], 'NPV [mln EUR]', pvals=pvals)
    dc1dmodel.stoch.expectation_plot(fwi[:, ILCOH, ITIME], 'LCOH [EUR/MWh]', pvals=pvals)


__main__ = main()

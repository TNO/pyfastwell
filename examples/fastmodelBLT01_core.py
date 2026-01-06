import os
from datetime import datetime

import numpy as np
import yaml

from pywellgeo.well_data.dc1dwell import Dc1dwell
from pythermonomics.geothermal_economics import GeothermalEconomics

from pyfastwell.stochastic.stochasticmodel import Stochasticmodel
import pyfastwell as pwf


ITIME = 0
IPOWER=1
ICOP =2
INPV =3
ILCOH=4
ITEMPPRD = 5
IFLOWRATE = 6




"""
this is a stochastic_dc1d model wrapper around the inputs_BLT01 tarjectory file, compose 
"""

def get_filename_noext(full_path):
    directory = os.path.dirname(full_path)
    base = os.path.basename(full_path)  # Gets "file.txt"
    root, ext = os.path.splitext(base)  # Splits into ("file", ".txt")
    res = os.path.join(directory,root)
    return res


def run_fastmodel(settingfile, trajectoryfile, k=100, khkv=3, dpBHP=50, segmentlength=50,cop_crit=1, verbose=False):
    """

    Parameters
    ----------
    settingfile: settings file for the techno-economic model
    trajectoryfile: trajectory file
    k:  permeability in x,y  [mDarcy]
    khkv: horizontal to vertical permeability ratio [-]
    dpBHP:  Bottom Hole Pressure difference between injector and producer [bar]
    segmentlength: segment length at reservoir depth for the perforated reservoir [m]
    cop_crit:  minimum COP.  If fricitional losses for given dpBHP is too high, it will be lowered accordingly
    verbose:  if true print results

    Returns
    -------

    """

    economics = GeothermalEconomics.from_trajectory(settingfile, trajectoryfile)
    npv, lcoe, cashflow, simdata, wells, well_results = economics.compute_economics()
    if verbose:
        print(f"Well cost {cashflow['wellcosts'].iloc[0] / 1e6} MEUR")

    fastmodel = pwf.WellFastModel(economics, k=k, khkv=khkv, segmentlength=segmentlength, verbose=verbose)

    fastmodel.set_dP_BHP(dpBHP)

    npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
    flowrate = fastmodel.flowrate
    if (cop<cop_crit):
        cop_old = cop
        flowrate_old = flowrate
        flowrate = flowrate * 1.1
        fastmodel.set_flowrate(flowrate)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
        iter = 0
        maxiter = 10000
        while ((abs(1/cop - 1/cop_crit) > 1e-4) and (iter < maxiter)):
            dfx = (1/cop - 1/cop_old) / (flowrate - flowrate_old)
            flowrate_old = flowrate
            cop_old = cop
            flowrate = flowrate -  (1/cop - 1/cop_crit) / dfx
            fastmodel.set_flowrate(flowrate)
            npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
            iter+=1

    if verbose:
        print(
            f"NPV {npv} EUR, LCOE {lcoe_val} EUR MWh-1, COP {cop}, Power {power} MW, Treservoir {fastmodel.production_temperature} C, \
            Tsurface {Tsurface} C, flowrate {flowrate} m3/h, DT_eff {DT_eff} C, DP_eff {DP_eff} bar")
        print(f"Well cost {cashflow['wellcosts'].iloc[0] / 1e6} MEUR")
    return npv, lcoe_val, cop, power, Tsurface, flowrate, fastmodel


def run_fastmodel_dploop(settingfile, trajectoryfile, k=100, khkv=3, dpLoop=50, segmentlength=50,  verbose=False):

    economics = GeothermalEconomics.from_trajectory(settingfile, trajectoryfile)
    npv, lcoe, cashflow, simdata, wells, well_results = economics.compute_economics()
    if verbose:
        print(f"Well cost {cashflow['wellcosts'].iloc[0] / 1e6} MEUR")

    fastmodel = pwf.WellFastModel(economics, k=k, khkv=khkv, segmentlength=segmentlength, verbose=verbose)

    fastmodel.set_dP_BHP(dpLoop)

    npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
    flowrate = fastmodel.flowrate
    if verbose:
        print(
            f"NPV {npv} EUR, LCOE {lcoe_val} EUR MWh-1, COP {cop}, Power {power} MW, Treservoir {fastmodel.production_temperature} C, \
            Tsurface {Tsurface} C, flowrate {flowrate} m3/h, DT_eff {DT_eff} C, DP_eff {DP_eff} bar")
        print(f"Well cost {cashflow['wellcosts'].iloc[0] / 1e6} MEUR")

    DP_eff_old = DP_eff
    flowrate_old = flowrate
    flowrate = flowrate * 1.1
    fastmodel.set_flowrate(flowrate)
    npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
    iter = 0
    maxiter = 10000
    while ( (abs(DP_eff-dpLoop)>1e-3) and (iter<maxiter) ):
        dfx = (DP_eff-DP_eff_old)/(flowrate-flowrate_old)
        flowrate_old = flowrate
        DP_eff_old = DP_eff
        flowrate = flowrate - (DP_eff-dpLoop)/dfx
        fastmodel.set_flowrate(flowrate)
        npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, Tsurface, DT_eff, DP_eff = fastmodel.compute_economics()
    return npv, lcoe_val, cop, power, Tsurface, flowrate, fastmodel

class   Fastmodel_stoch_BLT01(object):
    """
    StochasticFastModel class
    """

    def get_default_params(self):
        """
        get the parameters of the fastmodel
        :return: dictionary with the parameters
        """
        g = { 'top': 1925, 'H':  60, 'k': 50, 'khkv': 3,
              'DY': 800, 'Lh': 1000,  'DPBHP': 60, 'cop_crit': 10,
              'tinj': 36, 'tgrad': 0.03376 }
        return g

    def __init__(self, settingfile, trajectoryfile, locked_param=None, free_dist=None):
        """
        set the dc1dwell object
        :param dc1dsettings: file name yml with dc1dwell settings
        :param locked_param: dictionary with added or singular value parameters (defaulkt is none and considers the dedault from dc1dsettings
        :param free_dist: dictionary with the free parameters and their distributions
        """
        self.basesettingfile = settingfile
        self.settingfile = f'{get_filename_noext(settingfile)}_stoch.yml'
        self.basetrajectoryfile = trajectoryfile
        self.trajectoryfile = f'{get_filename_noext(trajectoryfile)}_stoch.yml'
        g = self.get_default_params()
        if locked_param is not None:
            g.update(locked_param)  # add singular value parameters to the dictionary
        free_param = list(free_dist.keys())
        argdict = {Stochasticmodel.LOCKED: g, Stochasticmodel.FREEPARAM: free_param,
                        Stochasticmodel.DISTPARAM: free_dist}
        """
        Initialize the StochasticFastmodel object
        :param argdict: dictionary with keys 'locked', 'free_param', and 'free_dist' containing the model parameters
        """

        self.stoch = Stochasticmodel(argdict=argdict)

    def forward(self, param, verbose=False):  # Forward model to which the parameters need to be fitted
        """
        :param  param: np.ndarray, including
        :param argdict: dictionary with the following keys in the dictionary:
            which can be used to update the parameters of the forward model in the param (if is included as distribition)
            k, H, topdepth, DP, return_temp, Q

        :return: forward model results in dimension of locked parameters 't'
        """

        t = self.stoch.get_from_dict('t', param)
        top = self.stoch.get_from_dict('top', param)
        H = self.stoch.get_from_dict('H', param)
        bottom = top + H
        k = self.stoch.get_from_dict('k', param)
        khkv = self.stoch.get_from_dict('khkv', param)
        DY  = self.stoch.get_from_dict('DY', param)
        Lh = self.stoch.get_from_dict('Lh', param)
        DPBHP = self.stoch.get_from_dict('DPBHP', param)
        cop_crit = self.stoch.get_from_dict('cop_crit', param)
        tinj = self.stoch.get_from_dict('tinj', param)
        tgrad = self.stoch.get_from_dict('tgrad', param)


        with open(self.basetrajectoryfile) as stream:
            try:
                dt = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        # make modifications to the  model yml to support sampled values
        dt['reservoir']['basic']['top_reservoir_depth_TVD'] = top
        dt['reservoir']['basic']['bottom_reservoir_depth_TVD'] = bottom
        dip = np.asin(H/Lh) * 180.0 / np.pi
        tomodify = dt['well_trajectories']['INJ1']['main_wellbore']['subs']
        tomodify['sub1']['dip'] = float(dip)
        tomodify['sub1']['y'] = float(-0.5 * DY)
        tomodify['sub1']['z'] = float(-top)
        tomodify['sub2']['L'] = float(Lh+0.1)
        tomodify = dt['well_trajectories']['PRD1']['main_wellbore']['subs']
        tomodify['sub1']['dip'] = float(dip)
        tomodify['sub1']['y'] = float(0.5 * DY)
        tomodify['sub1']['z'] = float(-top)
        tomodify['sub2']['L'] = float(Lh+0.1)
        with open(self.trajectoryfile, 'w') as f:
            yaml.dump(dt, f, default_flow_style=False)


        with open(self.basesettingfile) as stream:
            try:
                ds = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        # make modifications to the model yml to support sampled values
        ds['energy_loss_parameters']['tgrad'] = tgrad
        ds['reservoir_simulation_parameters']['injection_temperature'] = tinj
        with open(self.settingfile, 'w') as f:
            yaml.dump(ds, f, default_flow_style=False)

        npv, lcoe, cop, power, tempprd,flowrate, fastmodel = run_fastmodel(self.settingfile, self.trajectoryfile,  k=k, khkv=khkv, dpBHP=DPBHP, cop_crit=cop_crit, segmentlength=50, verbose=verbose)
        self.fastmodel = fastmodel
        fwd = np.vstack((t, [power],[cop], [npv*1e-6], [lcoe],  tempprd, flowrate))
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







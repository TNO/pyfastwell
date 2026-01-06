import os

from fastmodelBLT01_core import Fastmodel_stoch_BLT01,  IPOWER, ICOP, INPV, ILCOH, ITEMPPRD, IFLOWRATE, ITIME
from datetime import datetime

import numpy as np
import yaml

def main():

    np.random.seed(12345)

    print("Testing Dc1dmodel ensemble running code")
    stime = datetime.now()
    print ('start at: ', stime.time())

    nsamples = 500  # number of values in a prior, for one parameter

    t = np.arange(1)  # basically one time only

    trajectoryfile = 'inputBLT01/inputs_BLT01.yml'

    free_dist = {'k': {'dist': 'uniform', 'values': [240, 360]},
                 'Lh': {'dist': 'uniform', 'values': [80, 2000]} }


    locked_param =  { 't': t, 'top': 1925, 'H':  60, 'DY': 1200,   'DPBHP': 50,
              'tinj': 36, 'tgrad': 0.03376, 'cop_crit': 10 }

    settingfile = 'inputBLT01/npv_barros_138_BLT01.yml'

    model = Fastmodel_stoch_BLT01(settingfile, trajectoryfile, locked_param=locked_param, free_dist=free_dist)

    # get median values for the parameters and test the forward function
    funcparam = model.stoch.getmedianparams()
    res = model.forward(funcparam, verbose = True)
    param_units = ['[mDarcy]', '[m]']

    # for the median parameter values get the equivalent skin for Doubletcalc1D, plot the model
    L = locked_param['DY']
    dia_inch = 8
    rw = dia_inch * 0.0254 * 0.5
    skin_inj, skin_prd, ratio = model.fastmodel.getSkinFactors_dc1d(L, rw)
    print( f"for use in DC1D with L ={L} [m]  and diameter={dia_inch} [inch] Skin factors: inj {skin_inj}, prd {skin_prd}, ratio {ratio}")

    outdir = 'output/BLT01_300mDarcy'
    model.fastmodel.run_coldfront(simyears=30, rockdens=2700.0, rockcap=1000.0, porosity=0.21, gridrange=3000)
    model.fastmodel.plot_coldfront(outdir=outdir)
    model.fastmodel.plot_trajectories(outdir=outdir)


    M = model.generate_ensemble(nsamples)
    fwi = np.array(model.run_ensemble(M))

    etime = datetime.now()
    print ('end at: ', etime.time())
    print ('duration: ', etime - stime)


    # first index of the output array is the index of the realizations  in the model prediction ensemble, typically all :
    # second index is the output index of the model predictions, typically INPV or ICOP
    # third index is the time index in the model predictions, ITIME
    pvals = [10,30,50,70,90]
    resnames = ['POWER [MW]', 'COP [-]', 'NPV [mln EUR]', 'LCOE [EUR MWh-1]','TEMPPRD [C]', 'FLOWRATE [m3 h-1]']
    indices = [ IPOWER, ICOP, INPV, ILCOH, ITEMPPRD, IFLOWRATE]

    for ires, resname in enumerate(resnames):
        index = indices[ires]
        model.stoch.expectation_plot(fwi[:, index, ITIME], resname, pvals=pvals, filename_noext=os.path.join(outdir,'stoch'))
        for i,param in enumerate(model.stoch.argdict['free_param']):
            model.stoch.cross_plot(np.array(M)[:,i], f'{param} {param_units[i]}', fwi[:, index, ITIME], resname,filename_noext=os.path.join(outdir,'stoch'))



__main__ = main()

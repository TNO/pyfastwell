import copy

import numpy as np

from pyfastwell.plot.plotwells import plot_wells
from pyfastwell.wellflow.coldfront import Coldfront
from pythermonomics.geothermal_economics import GeothermalEconomics
import pywellgeo.well_data.water_properties as watprop
from pywellgeo.well_data.names_constants import Constants
from typing import Optional, Tuple

from pyfastwell.wellflow.wi_tno import WiTNO
import os


class WellFastModel:

    @staticmethod
    def skin_from_deviation(H: float, rw: float, khkv:float, deviation:float, type:str="rogers") -> float:
        """
        calculate the skin factor from the deviation of the well from vertical, it is valid for
        deviations from 10 up to 75 degrees




        Parameters
        ----------
        H:  vertical thickness of the reservoir in m

        rw: well radius in m

        khkv: ratio of kH to kV

        deviation: deviation from vertical in degrees

        type: type of correlation to use, options are "rogers", "abobaker", "cinco", "besson"

        Returns
        -------
        skin factor (dimensionless)

        Notes
        -----
        references:


        1. Abobaker, A., Al-Khulaifi, Y., Al-Menhali, A. S., & Al-Menhali, M. (2021). A New Correlation for Estimating the Productivity Index of Deviated Geothermal Wells. Energies, 14(17), 5585. https://doi.org/10.3390/en14175585
        2. CH. Cinco, F.G. Miller, H.J. Ramey Jr., 1975. *Unsteady-state pressure distribution created by a directionally drilled well. J. Petrol. Technol., 27 (11) (1975), pp. 1-392*
        3. E.J. Rogers, M.J. Economides, 1996.  *The skin due to slant of deviated wells in permeability-anisotropic reservoirs International Conference on Horizontal Well Technology, Society of Petroleum Engineers*
        4. Besson H., 1990. *Performance of slanted and horizontal wells on an anisotropic medium. European Petroleum Conference, Society of Petroleum Engineers (1990)*
        4. Hassan et al., 2021  https://www.sciencedirect.com/science/article/pii/S0920410520308317

        """

        h_d = H / rw
        theta_d = deviation/75.0
        if type=="rogers":
            iani = khkv**0.5
            if (iani>=1):
                skin = -2.48 * (np.sin(deviation * np.pi / 180))**5.87 * h_d**0.152 /  (iani**0.964)
            else:
                skin = -1.64 * (np.sin(deviation * np.pi / 180))**1.77 * h_d**0.184 /  (iani**0.821)
        elif type=="abobaker":
            skin = -(0.1 + 0.00025 * h_d + 1.852 * theta_d - 0.056 * khkv)**2
        elif type=="cinco":
            iani = khkv**0.5
            h_d *= iani
            theta_w = np.atan(np.tan(deviation * np.pi / 180)/iani)
            theta_w = theta_w * 180 / np.pi
            skin = - (theta_w/41)**2.06 - (theta_w/56)**1.865 * np.log10(h_d/100)
        elif type=="besson":
            L = H / np.cos(deviation * np.pi / 180)
            alfa = khkv**0.5
            gamma= ((1/alfa**2) + (H**2/L**2)* (1- 1/alfa**2))**0.5
            if (khkv==1):
               skin = np.log(4*rw/L) + (H/L) * np.log( (L*H)**0.5/(4*rw) )
            else:
               skin = np.log( (4 * rw / L) * (1/ (alfa*gamma))) + (H /(gamma*L)) * np.log( ((L * H) ** 0.5 / (4 * rw)) * (2 *alfa * gamma**0.5) /(1 + 1/gamma))
        return skin

    @staticmethod
    def pi_dc1d(L: float, rw: float, kH: float, mu: float, skin: float) -> float:
        """
        calculate the productivity index of a vertical well (flow rate [m3/s] for 1 bar pressure difference)
        in a reservoir with constant thickness (H) and permeability (k) and distance


        :param L:  distance between the wells in m
        :param rw: well radius in m
        :param kH: permeability times thickness in mD*m
        :param mu: viscosity in Pa*s
        :param skin: skin factor (dimensionless)
        :return: productivity index in m3/s/bar
        """

        scalefac = 2 * np.pi * kH * Constants.DARCY * 1e-3
        wf = mu * (np.log(L / rw) + skin)
        pi =  scalefac/(Constants.SI_BAR * wf)
        return pi

    def __init__(
        self,
        economics: GeothermalEconomics,
        k: float = 100.0,  # default permeability in mD
        khkv: float = 1.0,  # default ratio of kH to kV
        segmentlength: float = 50.0,  # default segment length in m
        verbose: Optional[bool] = False
    ) -> None:

        """

        Initialize the WellFastModel with the given economics and settings.
        This model computes the well flow and economics based on the provided parameters.
        It will always create a deep copy of the well trees to avoid modifying the original well trees in the economics instance.
        it will also modify the production temperature based on actual temperatures in the well tree if available.

        For the thickness of the reservoir, it will use the top and bottom reservoir depth from the INJ1 well.


        Parameters:
        -----------
        economics :  An instance of GeothermalEconomics containing the well trajectory and simulation results

        economics_settingfile :  settings file name for the economics model

        k : Horizontal permeability in mD

        khkv : Ratio of kH to kv

        segmentlength : Segment length in meters for the WiTNO model

        verbose :  If True, prints additional information during computation
        """


        self.economics = economics
        self.economics_config =  economics.economics_config
        self.verbose = verbose

        trj = self.economics.welltrajectory

        self.wells = [trj.tw['INJ1']['welltree'], trj.tw['PRD1']['welltree']]
        self.wells = [copy.deepcopy(w) for w in self.wells]  # deep copy to avoid modifying the original well trees
        #self.wells[0].init_ahd()
        #self.wells[1].init_ahd()
        self.topres = [trj.tw['INJ1']['main_wellbore']['top_reservoir_depth_TVD'],
                  trj.tw['PRD1']['main_wellbore']['top_reservoir_depth_TVD']]
        self.botres = [trj.tw['INJ1']['main_wellbore']['bottom_reservoir'],
                  trj.tw['PRD1']['main_wellbore']['bottom_reservoir']]


        self.production_temperature = trj.simresults.res_params.production_temperature
        try:
            # try to get the actual production temperature from the well tree
            tsoil= self.economics_config.energy_loss_parameters
            self.production_temperature = tsoil.tsurface + tsoil.tgrad * (self.topres[1]+self.botres[1]) *0.5
        except:
            pass
        self.injection_temperature = trj.simresults.res_params.injection_temperature
        self.k = k
        self.khkv = khkv
        self.segmentlength = segmentlength # default segment length in m
        salinity = trj.simresults.res_params.salinity
        self.flowrate = trj.simresults.res_params.flowrate
        if self.verbose:
            print(f"Flow rate set to: {self.flowrate} m3/s")
            print(f"AHD INJ1:  {self.wells[0].cumulative_ahd()} ")
            print(f"AHD PROD1: {self.wells[1].cumulative_ahd()} ")

        self.muprod = watprop.viscosity(self.production_temperature, salinity * 1e-6)
        self.muinj = watprop.viscosity(self.injection_temperature, salinity * 1e-6)
        self.wi_tno = WiTNO(self.wells[0], self.wells[1], -self.topres[0], -self.botres[0],
                    self.k, self.k, self.k/self.khkv, muinj=self.muinj, muprod=self.muprod,
                   segmentlength=self.segmentlength)
        self.H = self.botres[0] - self.topres[0]
        self.res, self.ii, self.pi = self.wi_tno.setupmatrix(n=0)
        if self.verbose:
            print(f"Productivity Index: {self.pi} m3/s/bar")
            print(f"Injectivity Index: {self.ii} m3/s/bar")



    def set_dP_BHP(self, dp_BHP: float):
        """
        Set the pressure difference (dP) for the well model, and updates the flow rate accordingly.
        The pressure difference is the difference between the bottom hole pressure (BHP) of the producer and the injector well.

        Parameters:
        -----------
        dp_BHP : Pressure difference in bar, which is the difference between the bottom hole pressure (BHP) of the producer and the injector well.

        Notes:
        ------
        raises ValueError if dp_BHP is not a positive value.

        """
        if dp_BHP <= 0:
            raise ValueError("Pressure difference must be a positive value.")
        w =  self.pi/self.ii
        flowrate = (dp_BHP / (w+1)) * self.pi * 3600  # convert m3/s to m3/h
        self.set_flowrate(flowrate)


    def get_BHPs(self):
        """
        get the injector and producer well BHP

        Returns
        -------
            BHPs [injector, producer] in bar
        """
        self.set_flowrate(self.flowrate)
        return np.asarray([ self.economics_config.reservoir_parameters.injection_BHP,  self.economics_config.reservoir_parameters.production_BHP])


    def set_flowrate(self, flowrate: float):
        """
        Set the flow rate for the well model, and updates the simulation results accordingly,
        including the actual productivity and injectivity indices

        Parameters
        ----------
        flowrate : Flow rate in m3/h

        Notes
        -----
        Raises ValueError if flowrate is not a positive value.
        This method also updates the economics configuration with the new flow rate and recalculates the simulation results.
        The bottom hole pressure (BHP) for both the injector and producer wells is calculated
        based on the flow rate and the injectivity/productivity indices.

        """

        if flowrate <= 0:
            raise ValueError("Flow rate must be a positive value.")
        self.flowrate = flowrate
        self.economics_config.reservoir_parameters.flowrate = flowrate
        refBHP = 0.5*( self.topres[0] + self.botres[0]) *9.81*1e-2  # in bar
        self.economics_config.reservoir_parameters.injection_BHP = refBHP +  (flowrate/3600)/ self.ii  # convert m3/h to m3/s and divide by injectivity index
        self.economics_config.reservoir_parameters.production_BHP = refBHP - (flowrate/3600)/ self.pi # convert m3/h to m3/s and divide by productivity index
        trj = self.economics.welltrajectory
        self.economics_config.reservoir_parameters.production_temperature = self.production_temperature
        self.economics.simresults = GeothermalEconomics._create_simresults(self.economics_config)
        self.economics.welltrajectory = trj


        if self.verbose:
            print(f"Flow rate set to: {self.flowrate} m3/h")


    def set_k_khkv(self, k: float, khkv: float, tryscale: bool = True):
        """
        Set the permeability and khkv for the well model. this will also update the WiTNO instance.
        It is couputationally expensive to set the permeability and khkv, so it should be done only when necessary, unless tryscale is True.
        If tryscale is True, it will scale the injectivity and productivity indices only if khkv is not altered.

        Parameters:
        ----------
        k:  Permeability in mD

        khkv: ratio of kH to kv

        tryscale: If True, it will simply scale the II and PI,  if khkv is not altered

        """

        if k <= 0 or khkv <= 0:
            raise ValueError("Permeability and khkv must be positive values.")

        if tryscale and self.khkv == khkv:
            # only scale the II and PI
            self.ii *= k / self.k
            self.pi *= k / self.k
            if self.verbose:
                print(f"Scaled II and PI by factor {k / self.k}")
            self.k = k
        else:
            # reset the WiTNO instance with new k and khkv
            self.k = k
            self.khkv = khkv
            self.wi_tno = WiTNO(self.wells[0], self.wells[1], -self.topres[0], -self.botres[0], self.k, self.k,
                                self.k / self.khkv, muinj=self.muinj, muprod=self.muprod,
                                segmentlength=self.segmentlength)
            self.res, self.ii, self.pi = self.wi_tno.setupmatrix(n=0)
            if self.verbose:
                print(f"Permeability set to: {self.k} mD*m")
                print(f"khkv set to: {self.khkv} mD*m")
        self.set_flowrate(self.flowrate)

    def getSkinFactors_dc1d(self, L:float, rw:float) -> Tuple[float, float]:
        """
        Get the skin factors for the injector and producer wells for use in DC1D based on the current well model.
        This method will use the WiTNO instance to calculate the skin factors for both wells. It will
        be consistent with the well distance (L) and well radius (rw), as to be used in the DC1D model.

        Parameters:
        ----------
        L : Distance between the wells in meters
        rw : Well radius in meters

        Returns:
        -------
        A tuple containing the skin factors for the injector and producer wells (skin_inj, skin_prod) and
        the ratio of actual II and II for a vertical injector well.

        Notes:
        -----
        The skin factor is a dimensionless number that represents the effect of near-wellbore conditions on the flow of fluid into or out of a well.
        A positive skin factor indicates a reduction in flow efficiency, while a negative skin factor indicates an enhancement in flow efficiency.
        The skin factors are calculated based on the current flow rate and the well model.

        """

        H = self.botres[0] - self.topres[0]
        ii_pi_vert = [ WellFastModel.pi_dc1d(L,rw, self.k*H, self.muinj, 0.0),
                       WellFastModel.pi_dc1d(L,rw, self.k*H, self.muprod, 0.0)]
        ratios = [self.ii / ii_pi_vert[0], self.pi / ii_pi_vert[1]]
        # ratios[i] =  (np.log(L / rw) + 0.0) / (np.log(L / rw) + skin[i])
        # skin[i] = np.log(L / rw) * (1 - ratios[i]) / ratios[i]
        skin_inj =  np.log(L / rw) * (1 - ratios[0]) / ratios[0]
        skin_prod = np.log(L / rw) * (1 - ratios[1]) / ratios[1]



        # check the skin factors by recalculating the II and PI
        ii_check = WellFastModel.pi_dc1d(L, rw, self.k *H, self.muinj, skin_inj)
        pi_check = WellFastModel.pi_dc1d(L, rw, self.k *H, self.muprod, skin_prod)
        if self.verbose:
            print(f"Calculated skin factors: skin_inj={skin_inj}, skin_prod={skin_prod}")
            print(f"Check II: {self.ii} vs {ii_check}, PI: {self.pi} vs {pi_check}")

        return skin_inj, skin_prod, ratios[0]



    def compute_economics(self):
        """
        Compute the economics of the well model using the configured parameters.
        This method will use the actual simulation results (which can be actualized by set_flowrate)
        and other parameters to calculate the NPV, LCOE, cash flow, simdata, wells, wellresults, COP, power, surface production temperature, DT_effective, and DP_loop.


        Returns:
        -------
        A tuple containing NPV [EUR], LCOE [EUR/kWh], cash flow [pd.dataframe], simdata, wells, wellresults, COP[-], power [MW],
        surface production temperature [C], DT_effective [C], DP_loop [bar]

        Notes:
        -----
        This method will print the NPV, LCOE, COP, and power output if verbose is set to True.
        The NPV is the net present value of the project, LCOE is the levelized cost of heat,
        COP is the coefficient of performance, and power is the power output in MW.
        The cash flow is also returned, which is the cash flow pandas dataframe of the project over the simulation period.
        :raises RuntimeError: If the economics model has not been set up correctly or if the simulation results are not available.

        """
        npv, lcoe_val, cashflow,  simdata, wells, well_results = self.economics.compute_economics()
        cop = self.economics.cop
        power = self.economics.power
        surface_production_temperature = self.injection_temperature+ well_results["dTemp_PRD1[C]"].values[-1]
        DT_effective = surface_production_temperature - self.injection_temperature
        DP_loop = well_results["dP_INJ1[bar]"].values[-1] + well_results["dP_PRD1[bar]"].values[-1]
        BHP = self.get_BHPs()
        if self.verbose:
            print(f"NPV: {npv}, LCOE: {lcoe_val}, COP: {cop}, Power: {power} MW, Surface Production Temperature: {surface_production_temperature} C, DT_effective: {DT_effective} C, DP_loop: {DP_loop} bar")
        return npv, lcoe_val, cashflow, simdata, wells, well_results, cop, power, surface_production_temperature, DT_effective, DP_loop

    def plot_trajectories(self, outdir: Optional[str] = None):
        """
        Plot the well trajectory and save it to a file if specified.
        This method will plot the well trajectories in 3D and 2D projections.

        Parameters:
        ----------
        outdir : If provided, the plot will be saved to the trajectory file base name in the specified outdir.
        Provide "" to save to the default directory.

        """
        trajectoryfile = self.economics.welltrajectory.trajectoryfile
        if outdir != None:
            filename_noext = os.path.join( outdir, str(os.path.basename(trajectoryfile)).split('.')[0])
        else:
            filename_noext = None
        plot_wells(self.wells, filename_noext=filename_noext)



    def run_coldfront(self, simyears: float, rockdens: float = 2700.0, rockcap: float = 1000.0, porosity: float = 0.21,
                      loadhours=8760,  maxlength: float  = 20, gridrange: float = 2000,  ncell:int = 200) -> Coldfront:
        """
        Run the cold front model for the well model.
        This method will create a Coldfront instance and run the cold front model for the specified number of years.
        The cold front model will use the well  injection temperature,
        salinity, rock density, rock heat capacity, porosity, and flow rate to calculate
        the cold front propagation in the reservoir. It does not take into thermal diffusivity in the reservoir and over and underburden layers.

        Parameters
        ----------
        simyears : Number of years to run the cold front model

        rockdens : Rock density in kg/m3

        rockcap : Rock heat capacity in J/(kg*K)

        porosity : Porosity of the rock

        loadhours : load hours

        maxlength : Maximum length of the cold front in meters

        gridrange : Range of the grid for the cold front model in meters

        ncell : Number of cells for the flow grid in the cold front model (in x and y direction)

        Returns:
        -------
        None

        """

        flowrate = self.flowrate

        self.coldfront = Coldfront (wellres=self.wi_tno, flowrate=flowrate,
                              tres=self.production_temperature, tinj=self.injection_temperature,
                              salinity=self.economics.simresults.res_params.salinity,
                              rockdens=rockdens, rockcap=rockcap, porosity=porosity, simyears=simyears, loadhours=loadhours, maxlength=maxlength, gridrange=gridrange, ncell=ncell)


    def plot_coldfront(self, outdir: Optional[str] = None, title: Optional[str]=None):
        """
        Plot the perforated section of the wells and cold front and save it to a file if specified.

        Parameters
        ----------
        outdir :  if provided, the plot will be saved to this directory.

        title :  if provided, the plot title, otherwise automatically generated

        Notes:
        -----
        raises RuntimeError: If the cold front model has not been run yet.

        """
        if not hasattr(self, 'coldfront'):
            raise RuntimeError("Cold front model has not been run yet. Please run run_coldfront() first.")

        trajectoryfile = self.economics.welltrajectory.trajectoryfile
        filename_noext = None
        if (outdir != None):
            filename_noext = os.path.join( outdir, str(os.path.basename(trajectoryfile)).split('.')[0])
        dps = self.get_BHPs()

        formatted = "{:.0f}".format(self.flowrate)
        dpformatted = "{:.2f}".format(dps[0]-dps[1])
        title = str(f"Thermal Front: {dpformatted} bar BHP, {formatted} m3/h, {self.coldfront.simyears} years ")
        self.coldfront.plot_coldfront(filename_noext=filename_noext, title=title )
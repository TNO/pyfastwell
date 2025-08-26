
import unittest
from scipy.interpolate import CubicSpline
from pythermonomics.geothermal_economics import GeothermalEconomics
from pywellgeo.well_data.dc1dwell import *
from pywellgeo.well_tree.well_tree_azim_dip import *
from pyfastwell.plot.plotwells import plot_wells



class MyTestCase(unittest.TestCase):
    def test_WellTreeTNO(self):
        wells = []
        radius =0.1
        w = WellTreeTNO.from_vertical(0, 0, -3000, radius)
        w.perforate(-2000, -2500)
        wells.append(w)
        w = WellTreeTNO.from_vertical(800, 0, -3000, radius)
        w.perforate(-2000, -2500)
        wells.append(w)

        xyz = np.asarray( [[200,0,0],[200,0,-800], [600,0,-2000],[1600,0,-2500]])
        xyz = xyz.transpose()
        w = WellTreeTNO.from_xyz(xyz[0], xyz[1], xyz[2], radius)
        w.perforate(-2000, -2500)
        wells.append(w)

        xyz = np.asarray( [[200,0,0],[200,0,-800], [600,500,-2000],[1600,500,-2500]])
        xyz = xyz.transpose()
        w.add_xyz(xyz[0], xyz[1], xyz[2], color='green', radius=radius)
        w.perforate(-2000, -2500)
        wells.append(w)

        plot_wells(wells, 'output/welltypes/basic/test_WellTree')





    def test_DC1D_trajectory(self):
        dc1dsettings = './input/dc1dwell.yml'
        settingfile = './input/npv_thermogis2024.yml'
        outdir = 'output/welltypes/dc1d'
        dc1dwell = Dc1dwell.from_configfile(dc1dsettings)
        dc1dwell.qvol =   125/3600
        dc1dwell.dp = 30
        dc1dwell.calculateDP_qvol()
        economics = GeothermalEconomics.from_dc1d(settingfile, dc1dsettings, dc1dwell)
        trj= economics.welltrajectory
        trj_inp = trj.trajectoryinput
        dc1d = trj_inp.dc1d
        tw = trj.tw
        ahd = tw['INJ1']['welltree'].cumulative_ahd()
        wells_and_states = trj_inp.getdefaultwellstates()
        for i,w in enumerate(wells_and_states.keys()):
            print(w, wells_and_states[w])
            if wells_and_states[w]=='prod':
                pass
        qvol = 125/3600
        templist = dc1d.temp
        salinity = dc1d.salinity[0]
        tubingdiameter_inch = dc1d.rw[0]*2 * Constants().SI_INCH
        roughness_minch = dc1d.roughness
        trj_inp.friction_all( qvol, templist, salinity, tubingdiameter_inch, roughness_minch, wells_and_states)
        self.assertAlmostEqual(trj_inp.dpSyphon, -2.001, delta=0.01)
        self.assertAlmostEqual(trj_inp.dp_frictioninj, 2.033, delta=0.01)
        self.assertAlmostEqual(trj_inp.dp_frictionprod,  1.931, delta=0.01)
        self.assertAlmostEqual(trj_inp.dpsum, 1.9633, delta=0.01)

        templosses = trj_inp.temploss_all(qvol, templist, salinity, wells_and_states)
        wells = [ tw['INJ1']['welltree'], tw['PRD1']['welltree']]
        filename_noext= os.path.join(outdir, str(os.path.basename(dc1dsettings)).split('.')[0])
        plot_wells(wells, filename_noext)

        dc1d.calculateDP(qvol)

        self.assertAlmostEqual(ahd, 2140, delta=1.0)
        self.assertAlmostEqual(dc1d.dp, 30, delta=0.05)
        self.assertAlmostEqual(dc1d.temp[1]-dc1d.tprod, templosses[1], delta=0.3)

    def test_WellTreeAzimDip(self):

        files = [
            './input/inputsDetailedTNOhor.yml',
            './input/inputsSimpleAzimDipEL_Tilburg2.yml'
        ]
        outdir = 'output/welltypes/azimdip'
        for trajectoryfile in files:
            settingfile = './input/npv_thermogis2024.yml'
            gt_economics = GeothermalEconomics.from_trajectory(
                settingfile,
                trajectoryfile=trajectoryfile,
            )
            wells = [ gt_economics.welltrajectory.tw['INJ1']['welltree'], gt_economics.welltrajectory.tw['PRD1']['welltree'] ]
            filename_noext= os.path.join(outdir, str(os.path.basename(trajectoryfile)).split('.')[0])
            plot_wells(wells, filename_noext )


    def test_xyzgeneric(self):
        trajectoryfile = './input/inputXYZ.yml'
        settingfile = './input/npv_thermogis2024.yml'
        outdir = 'output/welltypes/xyzgeneric'
        gt_economics = GeothermalEconomics.from_trajectory(
            settingfile,
            trajectoryfile=trajectoryfile,
        )
        wells = [ gt_economics.welltrajectory.tw['INJ1']['welltree'], gt_economics.welltrajectory.tw['PRD1']['welltree'] ]
        filename_noext= os.path.join(outdir, str(os.path.basename(trajectoryfile)).split('.')[0])
        #fname = os.path.basename(dc1dsettings)
        #outputnamebase = os.path.join(os.getcwd() + '\\', 'output', str(fname).split('.')[0])
        plot_wells(wells, filename_noext )

    def test_splinestest(self):
        outdir = 'output/welltypes/splinestest'
        xyz_tores   = np.asarray([[0, 0, 0], [0, 0, -500], [500, 500, -2300]])

        xyz_endres = np.asarray([1500, 500, 2400])
        xyz_prd = np.asarray([ [0,0,0],[-500,0,2300],[-500, 0, 2400]])

        # Define the t, x, y, and z coordinates of the points



        x = xyz_tores[:, 0]
        y = xyz_tores[:, 1]
        z = xyz_tores[:, 2]



        # Define the derivative at the first and last point for each dimension
        derivative_at_start = np.array([0, 0, -1])
        derivative_at_end = np.array([0.9, 0, -0.1])
        npoints =3
        diffs =  np.diff(xyz_tores, axis=0)
        distances = np.cumsum(np.linalg.norm(diffs, axis=1))
        t_ahd = np.insert(distances, 0,0 )
        # Create the cubic splines with the constrained derivative for each dimension
        cs_x = CubicSpline(t_ahd, x, bc_type=((1, derivative_at_start[0]), (1, derivative_at_end[0])))
        cs_y = CubicSpline(t_ahd, y, bc_type=((1, derivative_at_start[1]), (1, derivative_at_end[1])))
        cs_z = CubicSpline(t_ahd, z, bc_type=((1, derivative_at_start[2]), (1, derivative_at_end[2])))
        deltaseg = 100
        tpoints = np.arange(0, t_ahd[-1], deltaseg)
        x = cs_x(tpoints)
        y = cs_y(tpoints)
        z = cs_z(tpoints)
        xyz_inj = np.column_stack([(x, y, -z), xyz_endres]).T

        # Now you can use the cubic splines to estimate x, y, z at any point along t
        #
        data = {
            #"format": DoubleQuotedScalarString("XYZGENERIC"),
            "format": 'XYZGENERIC',
            "reservoir": {
                "basic": {
                    "top_reservoir_depth_TVD": 2300.0,
                    "bottom_reservoir_depth_TVD": 2400
                }
            },
            "well_trajectories": {
                "INJ1": {
                    "main_wellbore": {
                        "xyz": xyz_inj.tolist(),
                        "radius": 0.1,
                        "mindist": 25
                    }
                },
                "PRD1": {
                    "main_wellbore": {
                        "xyz":xyz_prd.tolist(),
                        "radius": 0.1,
                        "mindist": 25
                    }
                }
            }
        }
        output_dir = outdir
        os.makedirs(output_dir, exist_ok=True)
        trajectoryfile = os.path.join(output_dir, 'inputXYZspline.yml')
        with open(trajectoryfile, 'w') as f:
           yaml.dump(data, f, default_flow_style=False)


        minimal_settings = './input/npv_thermogis2024.yml'
        economics = GeothermalEconomics.from_trajectory(minimal_settings, trajectoryfile=trajectoryfile)

        npv, lcoe_val, cashflow, *_ = economics.compute_economics()

        wells = [economics.welltrajectory.tw['INJ1']['welltree'], economics.welltrajectory.tw['PRD1']['welltree']]
        filename_noext = os.path.join(output_dir, 'test_splinestest')
        plot_wells(wells, filename_noext)





if __name__ == '__main__':
    unittest.main()
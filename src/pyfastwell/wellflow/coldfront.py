
import numpy as np
import matplotlib.pyplot as plt

import pywellgeo.well_data.water_properties as watprop
import pyfastwell.wellflow.wi_tno as WiTNO
import xarray as xr
import time
import os


from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp

from shapely.ops import unary_union
from shapely.geometry import Polygon
from shapely.geometry import LineString

def _velocity_field(t, pos, vx_interp, vy_interp):
    # pos: [x, y]
    v = np.array([vx_interp(pos), vy_interp(pos)]).flatten()
    return v

def _rk45 (t_span: tuple[float, float], ini_pos: np.ndarray, vx_interp, vy_interp):
    """
    Adaptive RK4 integration for a velocity field defined by vx_interp and vy_interp.

    Parameters
    ----------
    t_span :  Time span \((0, sim\_time)\).

    ini_pos : Initial positions \([ [x, y], ... ]\).

    vx_interp : callable function of  Interpolated x-velocity field.

    vy_interp : callable function of Interpolated y-velocity field.

    Returns
    -------
    Updated positions after integration.
    """

    # RK45 integration
    results = []
    for pos in ini_pos:
        sol = solve_ivp(
            _velocity_field,
            t_span,
            pos,
            method='RK45',
            args=(vx_interp, vy_interp),
            rtol=1e-6,
            atol=1e-9)
        results.append (sol.y[:, -1])
    return np.array(results)



def _get_all_coords(geom):
    if geom.geom_type == 'Polygon':
        return [list(geom.exterior.coords)]
    elif geom.geom_type == 'MultiPolygon':
        return [list(p.exterior.coords) for p in geom.geoms]
    else:
        return []


def _subdivide_polygon_edges(polygon, max_length):
    coords = list(polygon.exterior.coords)
    new_coords = []
    for i in range(len(coords) - 1):
        p1 = np.array(coords[i])
        p2 = np.array(coords[i + 1])
        edge_len = np.linalg.norm(p2 - p1)
        n_segments = max(1, int(np.ceil(edge_len / max_length)))
        for j in range(n_segments):
            frac = j / n_segments
            new_point = p1 + frac * (p2 - p1)
            new_coords.append(tuple(new_point))
    new_coords.append(tuple(coords[-1]))  # close the polygon
    return Polygon(new_coords)

def _extend_line_start(x0, y0, x1, y1, length):
    dx = x1 - x0
    dy = y1 - y0
    norm = np.hypot(dx, dy)
    ux, uy = dx / norm, dy / norm
    new_x0 = x0 - ux * length
    new_y0 = y0 - uy * length
    new_x1 = x1
    new_y1 = y1
    return (new_x0, new_y0)

def _extend_line_end(x0, y0, x1, y1, length):
    dx = x1 - x0
    dy = y1 - y0
    norm = np.hypot(dx, dy)
    ux, uy = dx / norm, dy / norm
    new_x0 = x0
    new_y0 = y0
    new_x1 = x1 + ux * length
    new_y1 = y1 + uy * length
    return (new_x1, new_y1)

def _get_bounds_from_polygon(polygon):
    # Handles Polygon and MultiPolygon
    if polygon.geom_type == 'Polygon':
        coords = np.array(polygon.exterior.coords)
    elif polygon.geom_type == 'MultiPolygon':
        coords = np.vstack([np.array(p.exterior.coords) for p in polygon.geoms])
    else:
        raise ValueError("Input must be Polygon or MultiPolygon")
    xmin, ymin = coords.min(axis=0)
    xmax, ymax = coords.max(axis=0)
    return xmin, xmax, ymin, ymax


class Coldfront:

    def __init__(
        self,
        wellres: WiTNO,
        flowrate: float,
        tres: float,
        tinj: float,
        salinity: float,
        rockdens: float = 2700,
        rockcap: float = 1000,
        porosity: float = 0.21,
        simyears: int = 30,
        maxlength: float = 20,
        polwidth: float = 60,
        gridrange: float = 2000,
        ncell: int = 200,
        loadhours: int = 8760,
        verbose: bool = False
    ):
        """
        Simulates the cold front propagation in a reservoir flow grid.

        Parameters
        ----------

        wellres :   Well and reservoir model.

        flowrate : Flow rate in m³/h.

        tres : Reservoir temperature in K.

        tinj :  Injected water temperature in K.

        salinity : Salinity of injected water in ppm.

        rockdens :  Rock density in kg/m³.

        rockcap :   Rock heat capacity in J/(kg·K).

        porosity :    Rock porosity (fraction).

        simyears :    Number of simulation years.

        maxlength :    Max polygon edge length in m.

        polwidth :     Polygon width around well in m.

        gridrange :  Grid range around bounding box of injectors at reservoir level in m.

        ncell :  Number of grid cells in x and y for flow grid.

        loadhours :    Annual operating hours.

        verbose :  Print additional info.
        """

        self.H = wellres.ztop-wellres.zbottom
        self.wellres = wellres
        self.flowrate = flowrate  # m3/h
        self.tres = tres
        self.tinj= tinj
        brinedens = watprop.density(-wellres.ztop * 9.81 * 1060, tinj, salinity * 1e-6)
        brinecap = watprop.heatcapacity(tinj, salinity * 1e-6)  # J/kg/K
        bulkvolcap = porosity * brinedens * brinecap + (1 - porosity) * rockdens * rockcap # J/m3/K, assuming rock heat capacity of 2.5 MJ/m3/K
        self.lam = brinedens * brinecap   / bulkvolcap

        self.simyears = simyears
        if verbose:
            print ("lam ", self.lam, " brinedens ", brinedens, " brinecap ", brinecap, " bulvolkcap ", bulvolkcap)

        self.scaleflow = (self.flowrate/3600)/self.getwellQ()
        p_inj = self.scaleflow
        p_prod = self.scaleflow / self.wellres.pratio
        if verbose:
             print("pressure for injector, producer [bar] ", p_inj, p_prod)



        #inj_union = self.get_pol_from_pixels(max_length=maxlength)
        inj_union = self.get_inj_polygons(polwidth=polwidth, max_length=maxlength)

        xmin, xmax, ymin, ymax = _get_bounds_from_polygon(inj_union)
        xmin-= gridrange
        xmax+= gridrange
        ymin-= gridrange
        ymax+= gridrange

        # Create grid
        x = np.linspace(xmin, xmax, ncell)
        y = np.linspace(ymin, ymax, ncell)
        grid_x, grid_y = np.meshgrid(x, y, indexing='ij')
        da_x = xr.DataArray(grid_x, dims=('x', 'y'), coords={'x': x, 'y': y})
        da_y = xr.DataArray(grid_y, dims=('x', 'y'), coords={'x': x, 'y': y})

        self.da_x = da_x
        self.da_y = da_y
        inj_points = []
        pindex = [] # starting indices of the injector polygon tracers (front)
        for injector in _get_all_coords(inj_union):
             pindex.append(len(inj_points))
             for p in injector:
                 inj_points.append( (p[0], p[1]))
        pindex.append(len(inj_points))
        inj_points = np.array(inj_points)

        sources = self.setSources()



        # Initialize velocity fields
        vx = xr.zeros_like(da_x)
        vy = xr.zeros_like(da_y)

        for sx, sy, strength in sources:
            dx = da_x - sx
            dy = da_y - sy
            r2 = dx ** 2 + dy ** 2 + 1e-12  # avoid division by zero
            vx += strength * dx / r2
            vy += strength * dy / r2


        # Normalize or scale as needed (e.g., divide by 2*pi for 2D potential flow)
        vx /= (self.H * 2 * np.pi)
        vy /= (self.H * 2 * np.pi)

        self.vx = vx
        self.vy = vy

        # vx, vy: xarray.DataArray velocity fields
        vx_interp = RegularGridInterpolator((vx.x, vx.y), vx.values)
        vy_interp = RegularGridInterpolator((vy.x, vy.y), vy.values)

        tracer_positions = []

        tracer_positions = inj_points
        front_x = tracer_positions[:, 0]
        front_y = tracer_positions[:, 1]


        sim_time = simyears  # simulate for 30 years
        if verbose:
            print(f"Simulating tracer front for {sim_time} years")
        start = time.time()

        cooled_initial = inj_union.area* self.H
        brinevol_t0= cooled_initial/self.lam
        hours_t0 = brinevol_t0 / self.flowrate
        years_t0 = hours_t0/loadhours

        print("years_t0" , years_t0)

        self.years_t0 = years_t0

            # Use the RK4 adaptive method to compute the final positions
        result = _rk45(t_span=(0, max(0, (sim_time - years_t0)) * self.lam * loadhours * 3600),  # time span in seconds
                       ini_pos=tracer_positions,
                       vx_interp=vx_interp,
                       vy_interp=vy_interp
                       )

        end = time.time()

        if verbose:
            print(f"Simulation completed in {end - start:.2f} seconds")
            print(f"Number of tracers: {len(tracer_positions)}")


        self.pindex = pindex
        self.result = result
        self.front_x = front_x
        self.front_y = front_y
        self.inj_union = inj_union

    def plot_coldfront(self, filename_noext: str = None, title: str = None) -> None:
        """
        Plot the cold front of injected water in the flow grid.

        Parameters
        ----------

        filename_noext :  If provided, saves the plot to this filename (without extension), otherwise displays it.

        title : Title for the plot. If None, a default title is used.

        Returns
        -------
        None
        """

        fig, ax = plt.subplots()
        inj_union = self.inj_union
        result = self.result
        front_x = self.front_x
        front_y = self.front_y
        pindex = self.pindex
        # Initial positions

        for i1, i2 in zip(pindex[:-1], pindex[1:]):
            polygon = plt.Polygon(result[i1:i2], closed=True, facecolor='orange', edgecolor='orange')
            ax.add_patch(polygon)
        try:
            for pol in inj_union.geoms:
                patch = plt.Polygon(list(pol.exterior.coords), closed=True, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
                ax.add_patch(patch)
        except:
            patch = plt.Polygon(list(inj_union.exterior.coords), closed=True, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
            ax.add_patch(patch)
        injlines = self.get_perforation_lines(inj=True)
        prodlines = self.get_perforation_lines(inj=False)
        for idx, (x0,y0,x1,y1) in enumerate(injlines):
            label = 'Injector' if idx == 0 else None
            ax.plot([x0, x1], [y0, y1], color='blue', linewidth=1, label=label)
            ax.plot(x0, y0, marker='o', color='blue', markersize=2)
            ax.plot(x1, y1, marker='o', color='blue', markersize=2)
        for idx, (x0, y0, x1, y1) in enumerate(prodlines):
            label = 'Producer' if idx == 0 else None
            ax.plot([x0, x1], [y0, y1], color='red', linewidth=1, label=label)
            ax.plot(x0, y0, marker='o', color='red', markersize=2)
            ax.plot(x1, y1, marker='o', color='red', markersize=2)
        #ax.plot (front_x, front_y, color='orange', label='Initial Front')
        # Resulting positions (after advection)
        ax.plot(result[:, 0], result[:, 1], marker='o', color='black', markersize=1, label='Final Front')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.legend()
        if title == None:
            formatted = "{:.2f}".format(self.flowrate)
            plt.title(f"Thermal Front: {formatted} m3/h, {self.simyears} years ")
        else:
            plt.title(title)
        ax.set_aspect('equal')
        if filename_noext is not None:
            output_dir = os.path.dirname(filename_noext)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            plt.savefig(filename_noext + '_coldfront.png', dpi=300, bbox_inches='tight')
        else:
            plt.show()# Initial positions
        plt.close()




    def getwellQ(self) -> float:
        """
        Calculate the total flowrate of the wells in m³/s.

        Returns
        -------
        Total flowrate of the wells in m³/s.
        """

        wellres = self.wellres
        sumq = 0
        for i,xs in enumerate(wellres.wlist['xs']):
            xe = wellres.wlist['xe'][i]
            L =  wellres.wlist['L'][i]
            qL = wellres.wlist['res'][i]*L
            #print("branch xs, xe, L, qL ", xs, xe, L, qL)
            sumq += max(qL,0)
        # determine injection and production pressure achieve the flowrate

        return sumq


    def setSources(self) -> list:
        """
        Generate a list of source points for the reservoir flow model.

        Returns
        -------
        list of tuples, each tuple contains (x, y, strength) for a source point.
        """
        wellres = self.wellres
        sources = []
        for i, xs in enumerate(wellres.wlist['xs']):
            xe = wellres.wlist['xe'][i]
            L =  wellres.wlist['L'][i]
            qL = wellres.wlist['res'][i]*L* self.scaleflow
            #print("branch xs, xe, L, qL ", xs, xe, L, qL)
            xpos = 0.5* (xs + xe)
            sources.append((xpos[0], xpos[1], qL))

        return sources

    def get_perforation_lines(self, inj: bool = True) -> list:
        """
        Get perforation line segments for injectors or producers.

        Parameters
        ----------
        inj : If True, return injector lines; if False, return producer lines.

        Returns
        -------
        list of tuple, Each tuple contains (x0, y0, x1, y1) coordinates for a perforation line.
        """
        wellres = self.wellres
        lines = []
        for i, xs in enumerate(wellres.wlist['xs']):
            xe = wellres.wlist['xe'][i]
            L =  wellres.wlist['L'][i]
            qL = wellres.wlist['res'][i]
            if (inj and qL > 0) or (not inj and qL < 0):
                lines.append((xs[0], xs[1], xe[0], xe[1]))
        return lines

    def get_inj_polygons(self, polwidth=60, max_length=5):
        """
        Generate polygons around injector well segments.

        Parameters
        ----------
        polwidth : Width of the polygon around each injector line (default: 60).
        max_length : Maximum edge length for polygon subdivision (default: 5).

        Returns
        -------
        shapely.geometry.Polygon or MultiPolygon, Unified polygon(s) representing injector regions.
        """
        wellres = self.wellres
        lines =[]
        line= []
        for i, xs in enumerate(wellres.wlist['xs']):
            xe = wellres.wlist['xe'][i]
            L = wellres.wlist['L'][i]
            qL = wellres.wlist['res'][i]
            if qL > 0:
                if (len(line) == 0):
                    if ( (xs[0], xs[1]) != (xe[0], xe[1]) ):
                        line.append( (xs[0], xs[1]) )
                        line.append( (xe[0], xe[1]) )
                    else:
                        line.append( (xs[0]-0.1*polwidth, xs[1]) )
                        line.append(( xe[0]+0.1*polwidth, xe[1]))
                elif ( (xs[0],xs[1]) == line[-1] ):
                    line.append( (xe[0], xe[1]))
                else:
                    lines.append(line)
                    line = []
        if (len(line)>0):
            lines.append(line)

        # extend each of the poly lines with half the width
        for line in lines:
             x0,y0,x1,y1 = line[0][0], line[0][1], line[1][0], line[1][1]
             line[0] =  _extend_line_start(x0, y0, x1, y1, 0.5 * polwidth)
             x0, y0, x1, y1 = line[-2][0], line[-2][1], line[-1][0], line[-1][1]
             line[-1] = _extend_line_end(x0, y0, x1, y1, 0.5 * polwidth)

        pols = []
        for line in lines:
            linestring = LineString(line)
            pol = linestring.buffer(polwidth / 2, cap_style=2)  # cap_style=2 for flat ends
            pols.append(pol)
        inj_union = unary_union(pols)
        inj_union = inj_union.convex_hull
        #inj_union = inj_union.buffer(2)

        try:
            newpols = []
            for pol in inj_union.geoms:
                pol = _subdivide_polygon_edges(pol, max_length=max_length)  # Subdivide edges if needed
                newpols.append(pol)
            return unary_union(newpols)
        except:
            inj_union = _subdivide_polygon_edges(inj_union, max_length=max_length)
            return inj_union

















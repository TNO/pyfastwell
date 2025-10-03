from pywellgeo.well_tree.well_tree_tno import *
from pywellgeo.well_data.names_constants import Constants
import numpy as np
import pandas as pd

class WiTNO:
    '''
    well index calculation object
    '''



    def __init__ (self, inj: WellTreeTNO, prod: WellTreeTNO, ztop, zbottom,
                  kx, ky, kz, muinj=1e-3, muprod=1e-3, restrictbranch=None, segmentlength=100, verbose=False):
        """
        initialize the well index calculation object

        Parameters:
        -----------
        inj: WellTree model of the injection well

        prod: WellTree model of the production well

        ztop: perforated top of the reservoir [m] (minus sign for values below surface)

        zbottom: perforated bottom of the reservoir [m] (minus sign for values below surface)

        kx: horizontal permeability [mDarcy]

        ky: horizontal permeability [mDarcy]

        kz: vertical peremeability [mDarcy]

        muinj:  viscosity of brine  at injector [Pa s]

        muprod:  viscosity of brine at producer [Pa s]

        restrictbranch: restrict the active wells in the well Tree model to this name

        segmentlength: coarsening target for the perforated part of the reservoir (will not refine)

        verbose: if True, print debug information

        """
        self.injector = inj
        self.producer = prod

        self.ztop = ztop
        self.zbottom = zbottom

        mDarcy2SI = Constants.DARCY * 1e-3  # convert to mD*m
        self.kx = kx*mDarcy2SI
        self.ky = ky*mDarcy2SI
        self.kz = kz*mDarcy2SI
        self.k =  (self.kx*self.ky*self.kz)**(1.0/3.0)
        self.a = (self.ky * self.kz) ** 0.5 / self.k
        self.b = (self.kx * self.kz) ** 0.5 / self.k
        self.c = (self.kx * self.ky) ** 0.5 / self.k
        self.muinj = muinj
        self.muprod = muprod



        self.injector.perforate(self.ztop, self.zbottom)
        self.producer.perforate(self.ztop, self.zbottom)
        # get injector wells as list of segments
        self.injector.coarsen(segmentlength,  perforated=True)
        self.producer.coarsen(segmentlength,  perforated=True)
        self.s = np.asarray([(self.k/self.kx)**0.5, (self.k/self.ky)**0.5,(self.k/self.kz)**0.5])
        self._scale(self.s)

        wlistinj = self.injector.getbranch(name=restrictbranch, perforated=True)
        # set injector inj field to true
        wlistinj['inj'] = wlistinj['name']!=None
        wlistprod = self.producer.getbranch(name=restrictbranch,perforated=True)
        # set producer inj field to true
        wlistprod['inj'] = wlistprod['name'] == None
        self.wlist = pd.concat([wlistinj, wlistprod], ignore_index=True)

        self._calcRWeffall()
        self._calcL()


        if (verbose):
            print('xs', self.wlist['xs'] )
            print('xe', self.wlist['xe'])
            print('rw', self.wlist['rw'])
            print('rweff', self.wlist['rweff'])
            print('name', self.wlist['name'])
            print('inj', self.wlist['inj'])



    def _scale(self, s):
        self.producer.scale(s)
        self.injector.scale(s)
        self.ztop  *= s[2]
        self.zbottom *=  s[2]

    def _calcL(self):
        self.wlist['L'] = self.wlist['rw'] * 1.0
        for i,rw in enumerate(self.wlist['rw']):
            x1 = self.wlist.loc[i,'xs']
            x2 = self.wlist.loc[i,'xe']
            self.wlist.loc[i,'L'] = 0.5* np.linalg.norm(x1-x2)

    def _calcRWeffall(self):
        self.wlist['rweff'] = self.wlist['rw']*1.0
        for i,rw in enumerate(self.wlist['rw']):
            x1 = self.wlist.loc[i,'xs']
            x2 = self.wlist.loc[i,'xe']
            self.wlist.loc[i,'rweff'] = self._calcRweff(rw, x1, x2)



    def _calcRweff(self, rw, x1, x2):
        """
        calculate effective well radius for the segment x1,x2

        Notes:
        ------
        see Besson, J. (1997). "Well Index Calculation in a Boussinesq Medium." SPE Annual Technical Conference and Exhibition. doi:10.2118/38834-MS
        Note: this is a simplified version of the Besson equation, which assumes that the well is straight and the segment is horizontal.
        The Besson equation is given by:
        bb_radius_eq = 0.5* backbone.Radius/alpha^(1/3)*1/beta * sqrt((1+(beta^2/gamma))^2 + ((delta - 1/delta)* cos(theta)* cos(bb_azim)* sin(bb_azim)/gamma)^2);
        where:
        - bb_radius_eq is the effective well radius
        - backbone.Radius is the radius of the well
        - alpha is the ratio of the geometric mean of the horizontal permeabilities to the vertical permeability
        - beta is the ratio of the horizontal permeabilities
        - gamma is a function of the inclination and azimuth of the well
        - delta is the ratio of the horizontal permeabilities

        % bb_radius_eq
        a = (perm_y*perm_z)^0.5/perm;
        b = (perm_x*perm_z)^0.5/perm;
        c = (perm_x*perm_y)^0.5/perm;
        theta   = pi/2-bb_incl; % Note:Besson inclination definition is  90-our_inclination.
        alpha   = sqrt((perm_x*perm_y)^0.5/perm_z);
        beta    = sqrt(a/b * (cos(bb_azim))^2 + b/a * (sin(bb_azim))^2);
        gamma   = sqrt((cos(theta))^2 + (a^2/c^2 * (cos(bb_azim))^2 + b^2/c^2 * (sin(bb_azim))^2) * (sin(theta))^2);
        delta   = sqrt(perm_x/perm_y);
        bb_radius_eq = 0.5* backbone.Radius/alpha^(1/3)*1/beta * sqrt((1+(beta^2/gamma))^2 + ((delta - 1/delta)* cos(theta)* cos(bb_azim)* sin(bb_azim)/gamma)^2);


        :param rw: radius of the well
        :param x1: startung point of segment
        :param x2: end point of segment
        :return: rw effective
        """
        dvec = x2-x1
        dz = dvec[2]
        # dip/inclination, horizontal is 0.5 PI, vertical is 0
        theta = 0.5* np.pi - abs(np.arcsin(dvec[2]/np.dot(dvec,dvec)**0.5))
        # azimuth, x-axis is 0
        phi = np.arctan2(dvec[0], dvec[1])

        a = self.a
        b = self.b
        c = self.c
        alfa = ((self.kx * self.ky)**0.5 / self.kz)**0.5
        delta = (self.kx/self.ky)**0.5
        beta = ( (a/b) * np.cos(phi)**2 + (b/a) * np.sin(phi)**2 ) ** 0.5
        gamma = ( np.cos(theta)**2 +  ( (a**2/c**2) *np.cos(phi)**2 + (b**2/c**2) * np.sin(phi)**2 ) * np.sin(theta)**2 ) **0.5
        rweff = 0.5 * rw / alfa ** (1.0 / 3.0) * (1 / beta) * np.sqrt((1 + (beta ** 2 / gamma)) ** 2  + ((delta - 1 / delta) * np.cos(theta) * np.cos(phi) * np.sin(phi) / gamma) ** 2)
        return rweff


    def _get_Infast(self, n, a, b, c, s1, s2, I_0):
        """"

        :param n order of representation of pressure, 0 is constant, 1 linear variation, 2 second order polynomial
            normally mode 0 would be used
        :rtype: object


        """
        if n==0:
            return I_0
        elif n==1:
            return (s1-s2-0.5 *b *I_0)/c
        else:
            v1 = I_0;
            v2 = (s1 - s2 - 0.5 * b * v1) / c
            for i in  range(2,n+1):
                value = s1 + (-1) ** i * s2 - 0.5 * (2 * i - 1) * b * v2 - (i - 1) * a * v1
                value = value / (i * c)
                v1 = v2
                v2 = value
        return value





    def _setupcontrolpoints(self, n=0, dseg= 0.5):
        """
        on each segment it creates control points for the flow rate solution. Every segment contains n+1 control points
        (these are evenly distributed including the end points of the segment, and are moved dseg fraction to internal boundary)
        ( for one control point the location is chosen at the start of the segment!)

        it sets the folllowing attributes of self
        self.nSeg : number of segments
        self.ncontrol: n+1
        self.M :  number of equations, which is eaual to self.nSeg*self.ncontrol +1
        self.control:  the geometric location of the control points as list (size self.M-1) of one dimensional numpy arrays (size 3)
        :param n: number of control points on each segment
        :param dseg: fraction shrinkage for the segment ((default 4%)
        :return:
        """
        self.nSeg =  len(self.wlist)
        self.ncontrol = n+1
        self.M = self.nSeg*self.ncontrol + 1
        self.control = []

        for iseg in  range(0,self.nSeg):
            seg= []
            seg.append(self.wlist.loc[iseg, 'xs'])
            seg.append(self.wlist.loc[iseg, 'xe'])
            dvec = seg[1]-seg[0]
            # range starts with 0 to argument-1
            for i in range (self.ncontrol):
                if (n==0):
                    t = 0
                else:
                    t = i / (n*1.0)
                x = (1-t)* (seg[0] +dseg*dvec) + t*(seg[1]-dseg*dvec)
                # add the rw radius
                rweff =  self.wlist.loc[iseg,'rweff']
                # get the normal to the segment
                xr = dvec*1.0
                if (abs(dvec[0]+dvec[1])>1.0):
                    xr[1] = dvec[0]
                    xr[0] = -dvec[1];
                else:
                    xr[2] = dvec[1]
                    xr[1] = -dvec[2]
                # determine outer product of dvec*xr (which is in plane normal to dvec and xr)
                xr = np.cross(xr,dvec)
                xr = xr/np.linalg.norm(xr)
                # normalize
                x += xr*rweff
                self.control.append(x)

    def setupmatrix(self, n=0, refPprod=1e5, doprint=False):
        # setup system of equations, injection pressure given, producer pressure is dp lower
        self._setupcontrolpoints(n)
        self.bfac = []
        for iseg in range(self.nSeg):
            inj = self.wlist.loc[iseg, 'inj']
            if (inj):
                mu = self.muinj
            else:
                mu = self.muprod
            self.bfac.append(mu / (4 * np.pi * self.k))
            # self.bfac.append(np.euler_gamma* mu / (2 * np.pi * self.k))

        # setup the m matrix  and b vector
        m = np.zeros([self.M, self.M])
        m *= 0
        brhs = np.zeros([self.M])
        # last term in RHS is mass balance set to 0
        brhs[self.M - 1] = 0
        brhs *= 0

        for iseg in range(self.nSeg):
            seg = []
            seg.append(self.wlist.loc[iseg, 'xs'])
            seg.append(self.wlist.loc[iseg, 'xe'])
            # dvec = seg[1]-seg[0]
            # L = 0.5 * np.dot(dvec, dvec) ** 0.5
            # h = 0.5 * (dvec);
            L = self.wlist.loc[iseg, 'L']
            zmirror = np.asarray([self.ztop, self.zbottom])
            hres = self.ztop - self.zbottom
            nmirrors = int(3000 / hres)
            for im in range(-1, 2):
                # create three sources -1 is source itself, 0 and 1 the mirrors in ztop and zbottom
                segim = [seg[0] * 1.0, seg[1] * 1.0]
                for (imirror) in range(nmirrors):
                    if imirror == 0 or im >= 0:
                        if (im >= 0):
                            isign = 1
                            if (im > 0):
                                isign = -1
                            zm = zmirror[im] + isign * imirror * hres
                            for ip in range(2):
                                dzmirror = segim[ip][2] - zm
                                segim[ip][2] -= 2 * dzmirror

                        # midpoint of the (mirrored) source
                        mp = 0.5 * (segim[0] + segim[1]);
                        # L = 0.5 * np.dot(dvec, dvec) ** 0.5
                        if (n > 0):
                            dvec = segim[1] - segim[0]
                            h = 0.5 * (dvec);
                        for kseg in range(self.nSeg):
                            inj = self.wlist.loc[iseg, 'inj']
                            for icontrol in range(self.ncontrol):
                                # kseg, icontrol is the geometric control point (the row of the equations)
                                # kk is the row in the set of equations
                                kk = kseg * self.ncontrol + icontrol
                                x = self.control[kk]
                                if (n > 0):
                                    y = x - mp
                                    a = np.dot(y, y)
                                    b = -2 * np.dot(y, h)
                                    c = np.dot(h, h)
                                    s1 = np.sqrt(c + b + a);
                                    s2 = np.sqrt(c - b + a);

                                du = x - segim[0]
                                u = np.dot(du, du) ** 0.5
                                dv = x - segim[1]
                                v = np.dot(dv, dv) ** 0.5

                                I_0 = np.log((u + v + 2 * L) / (u + v - 2 * L)) / L
                                if (n > 0):
                                    I_0 = -np.log(max(1.0 - 4.0 * L / (u + v + 2.0 * L), 1e-12)) / L
                                # jcontrol is the integration control
                                for jcontrol in range(self.ncontrol):
                                    value = I_0
                                    if (n > 0):
                                        value = self._get_Infast(jcontrol, a, b, c, s1, s2, I_0)
                                    ii = iseg * self.ncontrol + jcontrol
                                    m[kk][ii] += self.bfac[kseg] * L * value

                                    m[self.M - 1][ii] = L
                                    # decide on RHS and inclusion of Pp
                                    if (inj):
                                        m[ii][self.M - 1] = 0
                                        brhs[ii] = refPprod
                                    else:
                                        m[ii][self.M - 1] = -1
                                        brhs[ii] = 0

        if (doprint):
            print(' dp ', brhs)
        res = np.linalg.solve(m, brhs)
        if (doprint):
            print('correct input ', np.dot(m, res))

        self.pratio = refPprod / res[-1]
        res[0:-1] *= 2
        injindex = prodindex = 0
        for iseg in range(self.nSeg):
            self.wlist.loc[iseg, 'res'] = res[iseg]
            inj = self.wlist.loc[iseg, 'inj']
            L = self.wlist.loc[iseg, 'L']
            if (inj):
                injindex += res[iseg] * L
            else:
                prodindex += res[iseg] * L * self.pratio

        self.m = m
        self.brhs = brhs

        sinv = 1.0 / self.s
        self._scale(sinv)

        return res, injindex, prodindex

    def setupmatrix_new(self, n:int =0, refPprod:float =1e5, verbose:bool =False) -> tuple:
        """
        setup the matrix for the well index calculation, using the control points

        Parameters:
        -----------
        n: order of the polynomial representation of the pressure, 0 is constant,
        1 is linear, 2 is quadratic

        refPprod: reference pressure [Pa] for injectivity and productivity calculation

        doprint: if True, print the matrix and the right hand side vector

        Returns
        -------
         tuple of:

         -  flowrates of the pressure solution in m3/s for each of the AE segments,
         positive values are injection, negative production

         - the injection index and the production index (the flow rate in m3/s at the injector and producer for the given reference pressure)

        """

        # setup system of equations, injection pressure given, producer pressure is dp lower
        self._setupcontrolpoints(n)
        self.bfac = []
        for iseg in range(self.nSeg):
            inj = self.wlist.loc[iseg, 'inj']
            if (inj):
                mu = self.muinj
            else:
                mu = self.muprod
            self.bfac.append(mu / (4 * np.pi * self.k))
            # self.bfac.append(np.euler_gamma* mu / (2 * np.pi * self.k))

        # setup the m matrix  and b vector
        m = np.zeros([self.M, self.M])
        m *= 0
        brhs = np.zeros([self.M])
        # last term in RHS is mass balance set to 0
        brhs[self.M - 1] = 0
        brhs *= 0

        for iseg in range(self.nSeg):
            seg = []
            seg.append(self.wlist.loc[iseg, 'xs'])
            seg.append(self.wlist.loc[iseg, 'xe'])
            # dvec = seg[1]-seg[0]
            # L = 0.5 * np.dot(dvec, dvec) ** 0.5
            # h = 0.5 * (dvec);
            L = self.wlist.loc[iseg, 'L']
            zmirror = np.asarray([self.ztop, self.zbottom])
            hres = self.ztop - self.zbottom
            nmirrors = int(3000 / hres)

            zmid = 0.5 * (self.ztop + self.zbottom)
            dzmid = [seg[0][2] - zmid, seg[1][2] - zmid]

            for imirror in range(-nmirrors, nmirrors + 1):
                # create all the mirrored and original sources in a single array
                segim = [seg[0] * 1.0, seg[1] * 1.0]
                sgn = 1 - 2 * (imirror % 2)  # sgn = +1 for even imirror; -1 for odd imirror
                for istartend in [0, 1]:
                    segim[istartend][2] = zmid + imirror * hres + sgn * dzmid[istartend]

                # midpoint of the (mirrored) source
                mp = 0.5 * (segim[0] + segim[1])
                # L = 0.5 * np.dot(dvec, dvec) ** 0.5
                if (n > 0):
                    dvec = segim[1] - segim[0]
                    h = 0.5 * (dvec);
                for kseg in range(self.nSeg):
                    inj = self.wlist.loc[iseg, 'inj']
                    for icontrol in range(self.ncontrol):
                        # kseg, icontrol is the geometric control point (the row of the equations)
                        # kk is the row in the set of equations
                        kk = kseg * self.ncontrol + icontrol
                        x = self.control[kk]
                        if (n > 0):
                            y = x - mp
                            a = np.dot(y, y)
                            b = -2 * np.dot(y, h)
                            c = np.dot(h, h)
                            s1 = np.sqrt(c + b + a);
                            s2 = np.sqrt(c - b + a);

                        du = x - segim[0]
                        u = np.dot(du, du) ** 0.5
                        dv = x - segim[1]
                        v = np.dot(dv, dv) ** 0.5

                        I_0 = np.log((u + v + 2 * L) / (u + v - 2 * L)) / L
                        if (n > 0):
                            I_0 = -np.log(max(1.0 - 4.0 * L / (u + v + 2.0 * L), 1e-12)) / L
                        # jcontrol is the integration control
                        for jcontrol in range(self.ncontrol):
                            value = I_0
                            if (n > 0):
                                value = self._get_Infast(jcontrol, a, b, c, s1, s2, I_0)
                            ii = iseg * self.ncontrol + jcontrol
                            m[kk][ii] += self.bfac[kseg] * L * value

                            m[self.M - 1][ii] = L
                            # decide on RHS and inclusion of Pp
                            if (inj):
                                m[ii][self.M - 1] = 0
                                brhs[ii] = refPprod
                            else:
                                m[ii][self.M - 1] = -1
                                brhs[ii] = 0

        if (verbose):
            print(' dp ', brhs)
        res = np.linalg.solve(m, brhs)
        if (verbose):
            print('correct input ', np.dot(m, res))

        self.pratio = refPprod / res[-1]
        res[0:-1] *= 2
        injindex = prodindex = 0
        for iseg in range(self.nSeg):
            self.wlist.loc[iseg, 'res'] = res[iseg]
            inj = self.wlist.loc[iseg, 'inj']
            L = self.wlist.loc[iseg, 'L']
            if (inj):
                injindex += res[iseg] * L
            else:
                prodindex += res[iseg] * L * self.pratio

        self.m = m
        self.brhs = brhs

        sinv = 1.0 / self.s
        self._scale(sinv)

        return res, injindex, prodindex



